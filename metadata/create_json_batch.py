#!/usr/bin/env python3
"""Create adjacent CMR JSON sidecars using creat_json.py in persistent workers."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import re
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


_generator = None
_options = None


def initialize_worker(script, options):
    global _generator, _options
    spec = importlib.util.spec_from_file_location("mmcth_metadata_generator", script)
    _generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(_generator)
    _options = options


def sidecar_is_current(sidecar, source, source_stat):
    try:
        if sidecar.stat().st_mtime_ns < max(source_stat.st_mtime_ns, _options["generator_mtime_ns"]):
            return False
        metadata = json.loads(sidecar.read_text())
        collection = metadata["CollectionReference"]
        archive = metadata["DataGranule"]["ArchiveAndDistributionInformation"][0]
        checksum = archive["Checksum"]
        return (
            metadata["GranuleUR"] == source.name
            and collection["ShortName"] == _options["collection_shortname"]
            and str(collection["Version"]) == _options["collection_version"]
            and archive["Name"] == source.name
            and archive["SizeUnit"] == "MB"
            and archive["Size"] == source_stat.st_size / (1024.0 * 1024.0)
            and checksum["Algorithm"] == "MD5"
            and re.fullmatch(r"[0-9a-fA-F]{32}", checksum["Value"]) is not None
            and metadata.get("Projects") == ([{"ShortName": _options["project"]}] if _options["project"] else None)
            and metadata.get("PGEVersionClass") == ({"PGEVersion": _options["pge_version"]} if _options["pge_version"] else None)
        )
    except (OSError, ValueError, KeyError, IndexError, TypeError, AttributeError):
        return False


def process_file(filename):
    source = Path(filename)
    sidecar = source.with_name(source.name + ".cmr.json")
    temporary = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    try:
        before = source.stat()
        if not _options["overwrite"] and sidecar_is_current(sidecar, source, before):
            return "skipped", filename, "", False
        metadata, mod03 = _generator.build_ummg(
            str(source), _options["mod03_root"],
            _options["collection_shortname"], _options["collection_version"],
            project=_options["project"], pge_version=_options["pge_version"],
        )
        after = source.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError("Source NetCDF changed while generating metadata; retry this file")
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(metadata, stream, indent=2)
            stream.write("\n")
        os.replace(temporary, sidecar)
        return "written", filename, "", mod03 is None
    except Exception as exc:
        temporary.unlink(missing_ok=True)
        return "failed", filename, f"{type(exc).__name__}: {exc}", False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--metadata-script", type=Path, default=Path(__file__).with_name("creat_json.py"))
    parser.add_argument("--mod03-root", default="/data/gdi/satellite/TerraDataArchive/MODIS/MOD03")
    parser.add_argument("--collection-shortname", default="MODMISR_L2_CP")
    parser.add_argument("--collection-version", default="V00")
    parser.add_argument("--project")
    parser.add_argument("--pge-version")
    parser.add_argument("--pattern", default="MODMISR_L2_CP*.nc")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--overwrite", action="store_true", help="Regenerate even current, valid JSON sidecars")
    parser.add_argument("--limit", type=int, help="Process only the first N files (for a small trial)")
    parser.add_argument("--log-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1 or (args.limit is not None and args.limit < 1):
        parser.error("--workers and --limit must be positive")
    if "/" in args.pattern:
        parser.error("--pattern must be a filename pattern, without directories")
    directory = args.input_dir.resolve()
    if not directory.is_dir():
        parser.error(f"Input directory does not exist: {directory}")
    script = args.metadata_script.resolve()
    if not script.is_file():
        parser.error(f"Metadata generator does not exist: {script}")
    files = sorted(entry.path for entry in os.scandir(directory)
                   if entry.is_file() and Path(entry.name).match(args.pattern))
    if args.limit:
        files = files[:args.limit]
    if not files:
        parser.error(f"No files match {args.pattern!r} directly inside {directory}")
    logs = args.log_dir.resolve()
    logs.mkdir(parents=True, exist_ok=True)
    options = {
        "mod03_root": args.mod03_root,
        "collection_shortname": args.collection_shortname,
        "collection_version": args.collection_version,
        "project": args.project, "pge_version": args.pge_version,
        "overwrite": args.overwrite, "generator_mtime_ns": script.stat().st_mtime_ns,
    }
    counts = Counter()
    fallback_count = 0
    start = time.monotonic()
    print(f"Input: {directory}; files: {len(files)}; workers: {args.workers}", flush=True)
    print("Sidecars: <source>.nc.cmr.json, beside each source NetCDF", flush=True)
    with (logs / "failures.tsv").open("w", newline="") as errors, (logs / "mod03_fallbacks.txt").open("w") as fallbacks:
        writer = csv.writer(errors, delimiter="\t")
        writer.writerow(["netcdf", "error"])
        with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize_worker,
                                 initargs=(str(script), options)) as pool:
            for number, (status, filename, message, fallback) in enumerate(pool.map(process_file, files, chunksize=8), 1):
                counts[status] += 1
                if status == "failed":
                    writer.writerow([filename, message])
                    errors.flush()
                    print(f"FAILED {filename}: {message}", flush=True)
                if fallback:
                    fallback_count += 1
                    fallbacks.write(filename + "\n")
                    fallbacks.flush()
                if number % 100 == 0 or number == len(files):
                    print(f"Processed {number}/{len(files)}; written={counts['written']}, "
                          f"skipped={counts['skipped']}, failed={counts['failed']}", flush=True)
    summary = {
        "input_directory": str(directory), "metadata_generator": str(script),
        "files_selected": len(files), "workers": args.workers,
        "written": counts["written"], "skipped_current": counts["skipped"],
        "failed": counts["failed"], "mod03_fallbacks": fallback_count,
        "elapsed_seconds": round(time.monotonic() - start, 1),
    }
    (logs / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)
    print(f"Logs: {logs}", flush=True)
    return 1 if counts["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
