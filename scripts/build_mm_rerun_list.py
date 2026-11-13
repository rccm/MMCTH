#!/usr/bin/env python3
"""Build a --modisid-file list from run manifests or an output/database audit."""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR))

from src.misc.run_status import expected_output_path, work_item_record  # noqa: E402


DEFAULT_DECODER = "modis_cloud_mask_byte0_bit0_gate_native_bits1_2_v2"
SUCCESS_STATUSES = {"SAVED", "SKIPPED_NO_OVERLAP"}


def read_jsonl_files(paths):
    for path in paths:
        with path.open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, start=1):
                line = line.strip()
                if not line:
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError as exc:
                    print(
                        f"WARNING: ignoring malformed {path}:{line_number}: {exc}",
                        file=sys.stderr,
                    )


def decoder_problem(path: str, required_decoder: str) -> str | None:
    if not required_decoder:
        return None
    try:
        from netCDF4 import Dataset

        with Dataset(path, "r") as dataset:
            variable = dataset.variables.get("MODIS_CloudMask")
            if variable is None:
                return "MODIS_CLOUD_MASK_MISSING"
            actual = getattr(variable, "cloud_mask_decoder_version", None)
            if actual != required_decoder:
                return f"DECODER_MISMATCH:{actual or 'missing'}"
    except Exception as exc:
        return f"UNREADABLE_OUTPUT:{type(exc).__name__}:{exc}"
    return None


def _decoder_problem_task(item):
    path, required_decoder = item
    return path, decoder_problem(path, required_decoder)


def audit_run(run_dir: Path, required_decoder: str):
    planned = {}
    for record in read_jsonl_files(sorted(run_dir.glob("planned_rank_*.jsonl"))):
        planned[str(record["gid"])] = record

    latest = {}
    for record in read_jsonl_files(sorted(run_dir.glob("events_rank_*.jsonl"))):
        latest[str(record["gid"])] = record

    details = []
    for gid, plan in planned.items():
        event = latest.get(gid)
        status = event.get("status") if event else "NOT_STARTED"
        output = plan.get("expected_output")

        if status not in SUCCESS_STATUSES:
            reason = status
        elif status == "SKIPPED_NO_OVERLAP":
            continue
        elif not output or not Path(output).is_file():
            reason = "SAVED_OUTPUT_MISSING"
        elif str(output).endswith("_MM.nc"):
            reason = decoder_problem(output, required_decoder)
            if reason is None:
                continue
        else:
            continue

        details.append(
            {
                **plan,
                "last_status": status,
                "reason": reason,
            }
        )
    return details


def database_groups(database: Path, start_date: str, end_date: str):
    start_ts = f"{start_date} 00:00:00"
    end_ts = f"{end_date} 23:59:59"
    query = """
        SELECT f.group_id, f.file_path, f.orbit_number
        FROM files AS f
        WHERE f.group_id IN (
            SELECT DISTINCT group_id
            FROM files
            WHERE date >= ? AND date <= ?
        )
        ORDER BY f.group_id, f.file_id
    """
    grouped = defaultdict(lambda: {"files": [], "orbit": None})
    with sqlite3.connect(database) as connection:
        for gid, file_path, orbit in connection.execute(query, (start_ts, end_ts)):
            grouped[gid]["files"].append(file_path)
            if grouped[gid]["orbit"] is None:
                grouped[gid]["orbit"] = orbit
    return grouped


def audit_database(
    database: Path,
    start_date: str,
    end_date: str,
    output_dir: Path,
    required_decoder: str,
    workers: int,
    trust_after_mtime: float | None,
):
    details = []
    existing = []
    for gid, group in database_groups(database, start_date, end_date).items():
        files = group["files"]
        orbit = group["orbit"]
        output = expected_output_path(
            files,
            orbit,
            str(output_dir),
            check_input_exists=False,
        )
        if not output or not output.endswith("_MM.nc"):
            continue

        output_path = Path(output)
        if not output_path.is_file():
            reason = "MISSING_OUTPUT"
        elif trust_after_mtime is not None and output_path.stat().st_mtime >= trust_after_mtime:
            continue
        else:
            existing.append(
                (
                    work_item_record(
                        gid,
                        files,
                        orbit,
                        str(output_dir),
                        check_input_exists=False,
                    ),
                    output,
                )
            )
            continue

        record = work_item_record(
            gid,
            files,
            orbit,
            str(output_dir),
            check_input_exists=False,
        )
        record.update({"last_status": "UNKNOWN_OLD_RUN", "reason": reason})
        details.append(record)

    task_args = ((output, required_decoder) for _record, output in existing)
    if workers <= 1:
        results = map(_decoder_problem_task, task_args)
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        results = executor.map(_decoder_problem_task, task_args, chunksize=64)

    try:
        for (record, output), (_checked_output, reason) in zip(existing, results):
            if reason is None:
                continue
            record.update({"last_status": "UNKNOWN_OLD_RUN", "reason": reason})
            details.append(record)
    finally:
        if workers > 1:
            executor.shutdown()
    return details


def write_outputs(details, output: Path, details_path: Path):
    modis_ids = sorted({d.get("modis_id", "") for d in details if d.get("modis_id")})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{modis_id}\n" for modis_id in modis_ids), encoding="utf-8")

    fields = [
        "gid",
        "modis_id",
        "modis_file",
        "misr_file",
        "orbit",
        "expected_output",
        "last_status",
        "reason",
    ]
    with details_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(details)

    counts = defaultdict(int)
    for detail in details:
        counts[detail["reason"].split(":", 1)[0]] += 1
    print(f"Wrote {len(modis_ids)} unique MODIS IDs to {output}")
    print(f"Wrote {len(details)} affected groups to {details_path}")
    for reason, count in sorted(counts.items()):
        print(f"  {reason}: {count}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--run-dir", type=Path, help="A ds_output/.../_run_logs/<run_id> directory")
    source.add_argument("--database", type=Path, help="SQLite inputfiles_all database for an output audit")
    parser.add_argument("--start-date", help="Required with --database: YYYY-MM-DD")
    parser.add_argument("--end-date", help="Required with --database: YYYY-MM-DD")
    parser.add_argument("--output-dir", type=Path, help="Required with --database: product directory")
    parser.add_argument("--output", type=Path, required=True, help="Output MODIS-ID list for main.py --modisid-file")
    parser.add_argument("--details", type=Path, help="Output TSV details (default: <output>.details.tsv)")
    parser.add_argument("--required-decoder", default=DEFAULT_DECODER)
    parser.add_argument("--workers", type=int, default=8, help="Parallel NetCDF header readers for --database (default: 8)")
    parser.add_argument(
        "--trust-outputs-newer-than",
        help=(
            "Optional ISO timestamp for --database. Outputs with an mtime at or after "
            "this known decoder-deployment time are trusted without opening their headers."
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    details_path = args.details or args.output.with_name(args.output.name + ".details.tsv")

    if args.run_dir:
        details = audit_run(args.run_dir, args.required_decoder)
    else:
        missing = [
            name
            for name, value in (
                ("--start-date", args.start_date),
                ("--end-date", args.end_date),
                ("--output-dir", args.output_dir),
            )
            if value is None
        ]
        if missing:
            raise SystemExit(f"Required with --database: {', '.join(missing)}")
        trust_after_mtime = None
        if args.trust_outputs_newer_than:
            trust_after_mtime = datetime.fromisoformat(
                args.trust_outputs_newer_than
            ).timestamp()
        details = audit_database(
            args.database,
            args.start_date,
            args.end_date,
            args.output_dir,
            args.required_decoder,
            max(1, args.workers),
            trust_after_mtime,
        )

    write_outputs(details, args.output, details_path)


if __name__ == "__main__":
    main()
