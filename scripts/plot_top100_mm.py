#!/usr/bin/env python3
"""Run makeplots/plotting.py for MM files listed in the top-100 CSV."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


REPO_DIR = Path(__file__).resolve().parents[1]
DEFAULT_CSV = REPO_DIR / "mm_flag_ratio_2015_top100.csv"
DEFAULT_PLOTTING = REPO_DIR / "makeplots" / "plotting.py"
K9_PYTHON = Path("/data/keeling/a/gzhao1/miniconda3/envs/k9/bin/python")
PLOT_SUFFIXES = ("_1_hgt.png", "_2_ctp.png", "_3_diff.png")


def parse_args() -> argparse.Namespace:
    default_python = str(K9_PYTHON if K9_PYTHON.exists() else sys.executable)

    parser = argparse.ArgumentParser(
        description="Make MM diagnostic plots for rows in mm_flag_ratio_2015_top100.csv."
    )
    parser.add_argument("--csv", default=str(DEFAULT_CSV), help="CSV with a path column.")
    parser.add_argument(
        "--plotting",
        default=str(DEFAULT_PLOTTING),
        help="Path to makeplots/plotting.py.",
    )
    parser.add_argument(
        "--python",
        default=default_python,
        help="Python executable used to run plotting.py.",
    )
    parser.add_argument(
        "--log-dir",
        default=str(REPO_DIR / "logs" / "plot_top100_mm"),
        help="Directory for one log file per plotted granule.",
    )
    parser.add_argument("--jobs", type=int, default=1, help="Parallel plotting jobs.")
    parser.add_argument("--limit", type=int, help="Only process the first N selected rows.")
    parser.add_argument("--start-rank", type=int, help="Skip rows with rank lower than this.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate plots even when all three expected PNGs already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running plotting.py.",
    )
    return parser.parse_args()


def mm_path_from_row(row: dict[str, str]) -> Path:
    mm_path = Path(row["path"])
    if not mm_path.is_absolute():
        mm_path = REPO_DIR / mm_path
    return mm_path


def expected_pngs(mm_path: Path) -> list[Path]:
    return [mm_path.parent / "plots" / f"{mm_path.stem}{suffix}" for suffix in PLOT_SUFFIXES]


def read_rows(csv_path: Path, start_rank: int | None, limit: int | None) -> list[dict[str, str]]:
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "path" not in (reader.fieldnames or []):
            raise ValueError(f"{csv_path} does not have a 'path' column")

        rows = []
        for row in reader:
            if start_rank is not None and int(row["rank"]) < start_rank:
                continue
            rows.append(row)
            if limit is not None and len(rows) >= limit:
                break
    return rows


def run_row(row: dict[str, str], args: argparse.Namespace) -> tuple[str, str, str]:
    rank = row.get("rank", "?")
    mm_path = mm_path_from_row(row)
    outputs = expected_pngs(mm_path)

    if not mm_path.exists():
        return "missing", rank, f"{mm_path} does not exist"

    if not args.force and all(path.exists() for path in outputs):
        return "skipped", rank, f"{mm_path.name}: all expected plots already exist"

    cmd = [args.python, args.plotting, "-f", str(mm_path), "--quiet"]
    command_text = shlex.join(cmd)

    if args.dry_run:
        return "dry-run", rank, command_text

    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"rank_{int(rank):03d}_{mm_path.stem}.log"

    start = time.monotonic()
    with log_file.open("w") as handle:
        handle.write(f"$ {command_text}\n\n")
        handle.flush()
        proc = subprocess.run(
            cmd,
            cwd=REPO_DIR,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    elapsed = time.monotonic() - start

    if proc.returncode == 0:
        return "done", rank, f"{mm_path.name} in {elapsed:.1f}s; log={log_file}"
    return "failed", rank, f"{mm_path.name} exit={proc.returncode}; log={log_file}"


def main() -> int:
    args = parse_args()
    rows = read_rows(Path(args.csv), args.start_rank, args.limit)
    if not rows:
        print("No rows selected.")
        return 0

    print(f"Selected {len(rows)} row(s) from {args.csv}")
    print(f"Using plotting script: {args.plotting}")
    print(f"Using Python: {args.python}")

    counts: dict[str, int] = {}
    if args.jobs == 1:
        for row in rows:
            status, rank, message = run_row(row, args)
            counts[status] = counts.get(status, 0) + 1
            print(f"[rank {rank}] {status}: {message}", flush=True)
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_row, row, args) for row in rows]
            for future in as_completed(futures):
                status, rank, message = future.result()
                counts[status] = counts.get(status, 0) + 1
                print(f"[rank {rank}] {status}: {message}", flush=True)

    summary = ", ".join(f"{key}={counts[key]}" for key in sorted(counts))
    print(f"Summary: {summary}")
    return 1 if counts.get("failed") or counts.get("missing") else 0


if __name__ == "__main__":
    raise SystemExit(main())
