#!/usr/bin/env python3
"""Count MM_Flag ratios and multilayer confidence totals for MM NetCDF granules.
python scripts/mm_flag_ratio.py ds_output/2015 \
>   --workers 8 \
>   --output-prefix ds_output/2015_mm_flag_ratio \
>   --progress-every 1000
For each granule:
  numerator   = pixels where MM_Flag is 8 or 9
  denominator = pixels where MM_Flag is 6, 7, 8, or 9
  ratio       = numerator / denominator

The script writes:
  * a per-granule CSV
  * a top-N granule CSV sorted by count(8 or 9)
  * a JSON summary with total counts, the overall ratio, and
    MM_MultilayerConfidenceFlag totals for values 0 through 4, including
    totals where MM_Flag is 8 or 9
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import ctypes.util
import heapq
import json
import math
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


NC_NOWRITE = 0
DEFAULT_PATTERN = "*V00*_MM.nc"
DEFAULT_VARIABLE = "MM_Flag"
DEFAULT_MULTILAYER_VARIABLE = "MM_MultilayerConfidenceFlag"
COUNT_FLAGS = (6, 7, 8, 9)
NUMERATOR_FLAGS = (8, 9)
MULTILAYER_CONFIDENCE_FLAGS = (0, 1, 2, 3, 4)
MULTILAYER_COUNT_FIELDS = tuple(
    f"count_multilayer_confidence_{flag}" for flag in MULTILAYER_CONFIDENCE_FLAGS
)
MM89_MULTILAYER_COUNT_FIELDS = tuple(
    f"count_mm_flag_89_multilayer_confidence_{flag}"
    for flag in MULTILAYER_CONFIDENCE_FLAGS
)
CSV_FIELDS = [
    "status",
    "granule",
    "path",
    "count_6",
    "count_7",
    "count_8",
    "count_9",
    "count_6789",
    "count_89",
    "ratio_89_over_6789",
    *MULTILAYER_COUNT_FIELDS,
    *MM89_MULTILAYER_COUNT_FIELDS,
    "multilayer_error",
    "error",
]

_NETCDF_LIB: ctypes.CDLL | None = None


def _load_netcdf_lib() -> ctypes.CDLL:
    global _NETCDF_LIB
    if _NETCDF_LIB is not None:
        return _NETCDF_LIB

    lib_name = ctypes.util.find_library("netcdf") or "libnetcdf.so"
    lib = ctypes.CDLL(lib_name)

    lib.nc_open.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
    lib.nc_open.restype = ctypes.c_int

    lib.nc_close.argtypes = [ctypes.c_int]
    lib.nc_close.restype = ctypes.c_int

    lib.nc_inq_varid.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
    lib.nc_inq_varid.restype = ctypes.c_int

    lib.nc_inq_varndims.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
    lib.nc_inq_varndims.restype = ctypes.c_int

    lib.nc_inq_vardimid.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
    lib.nc_inq_vardimid.restype = ctypes.c_int

    lib.nc_inq_dimlen.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.nc_inq_dimlen.restype = ctypes.c_int

    lib.nc_get_var_schar.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
    lib.nc_get_var_schar.restype = ctypes.c_int

    lib.nc_strerror.argtypes = [ctypes.c_int]
    lib.nc_strerror.restype = ctypes.c_char_p

    _NETCDF_LIB = lib
    return lib


def _check_nc(lib: ctypes.CDLL, rc: int, context: str) -> None:
    if rc == 0:
        return
    message = lib.nc_strerror(rc).decode("utf-8", errors="replace")
    raise RuntimeError(f"{context}: {message}")


def _read_int8_variable(path: str, variable: str) -> np.ndarray:
    lib = _load_netcdf_lib()

    ncid = ctypes.c_int()
    _check_nc(lib, lib.nc_open(os.fsencode(path), NC_NOWRITE, ctypes.byref(ncid)), "open")

    try:
        varid = ctypes.c_int()
        _check_nc(
            lib,
            lib.nc_inq_varid(ncid.value, variable.encode("utf-8"), ctypes.byref(varid)),
            f"find variable {variable}",
        )

        ndims = ctypes.c_int()
        _check_nc(
            lib,
            lib.nc_inq_varndims(ncid.value, varid.value, ctypes.byref(ndims)),
            f"read {variable} dimension count",
        )

        dimids = (ctypes.c_int * ndims.value)()
        _check_nc(
            lib,
            lib.nc_inq_vardimid(ncid.value, varid.value, dimids),
            f"read {variable} dimension ids",
        )

        shape: list[int] = []
        for dimid in dimids:
            dimlen = ctypes.c_size_t()
            _check_nc(
                lib,
                lib.nc_inq_dimlen(ncid.value, dimid, ctypes.byref(dimlen)),
                f"read {variable} dimension length",
            )
            shape.append(dimlen.value)

        data = np.empty(shape, dtype=np.int8)
        _check_nc(
            lib,
            lib.nc_get_var_schar(
                ncid.value,
                varid.value,
                data.ctypes.data_as(ctypes.c_void_p),
            ),
            f"read {variable}",
        )
        return data
    finally:
        lib.nc_close(ncid.value)


def count_fields_for_flags(
    data: np.ndarray,
    fields: tuple[str, ...],
    flags: tuple[int, ...],
) -> dict[str, int]:
    bins = np.bincount(data.ravel().view(np.uint8), minlength=256)
    return {field: int(bins[flag]) for field, flag in zip(fields, flags)}


def empty_multilayer_counts() -> dict[str, int]:
    return {field: 0 for field in MULTILAYER_COUNT_FIELDS}


def empty_mm89_multilayer_counts() -> dict[str, int]:
    return {field: 0 for field in MM89_MULTILAYER_COUNT_FIELDS}


def process_file(task: tuple[str, str, str]) -> dict[str, Any]:
    path, variable, multilayer_variable = task
    granule = Path(path).name
    try:
        data = _read_int8_variable(path, variable)
        bins = np.bincount(data.ravel().view(np.uint8), minlength=256)
        count_6 = int(bins[6])
        count_7 = int(bins[7])
        count_8 = int(bins[8])
        count_9 = int(bins[9])
        count_6789 = count_6 + count_7 + count_8 + count_9
        count_89 = count_8 + count_9
        ratio = count_89 / count_6789 if count_6789 else math.nan

        try:
            multilayer_data = _read_int8_variable(path, multilayer_variable)
            if multilayer_data.shape != data.shape:
                raise ValueError(
                    f"{multilayer_variable} shape {multilayer_data.shape} "
                    f"does not match {variable} shape {data.shape}"
                )
            multilayer_counts = count_fields_for_flags(
                multilayer_data,
                MULTILAYER_COUNT_FIELDS,
                MULTILAYER_CONFIDENCE_FLAGS,
            )
            mm89_multilayer_counts = count_fields_for_flags(
                multilayer_data[np.isin(data, NUMERATOR_FLAGS)],
                MM89_MULTILAYER_COUNT_FIELDS,
                MULTILAYER_CONFIDENCE_FLAGS,
            )
            multilayer_error = ""
        except Exception as exc:
            multilayer_counts = empty_multilayer_counts()
            mm89_multilayer_counts = empty_mm89_multilayer_counts()
            multilayer_error = f"{type(exc).__name__}: {exc}"

        return {
            "status": "ok",
            "granule": granule,
            "path": path,
            "count_6": count_6,
            "count_7": count_7,
            "count_8": count_8,
            "count_9": count_9,
            "count_6789": count_6789,
            "count_89": count_89,
            "ratio_89_over_6789": ratio,
            **multilayer_counts,
            **mm89_multilayer_counts,
            "multilayer_error": multilayer_error,
            "error": "",
        }
    except Exception as exc:
        return {
            "status": "error",
            "granule": granule,
            "path": path,
            "count_6": 0,
            "count_7": 0,
            "count_8": 0,
            "count_9": 0,
            "count_6789": 0,
            "count_89": 0,
            "ratio_89_over_6789": math.nan,
            **empty_multilayer_counts(),
            **empty_mm89_multilayer_counts(),
            "multilayer_error": "",
            "error": f"{type(exc).__name__}: {exc}",
        }


def ratio_to_csv(value: float) -> str:
    return "" if math.isnan(value) else f"{value:.12g}"


def csv_to_ratio(value: str) -> float:
    return float(value) if value else math.nan


def row_for_csv(row: dict[str, Any]) -> dict[str, Any]:
    out = dict(row)
    out["ratio_89_over_6789"] = ratio_to_csv(float(row["ratio_89_over_6789"]))
    return out


def row_from_csv(row: dict[str, str]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "status": row["status"],
        "granule": row["granule"],
        "path": row["path"],
        "multilayer_error": row.get("multilayer_error", ""),
        "error": row.get("error", ""),
        "ratio_89_over_6789": csv_to_ratio(row.get("ratio_89_over_6789", "")),
    }
    count_fields = (
        "count_6",
        "count_7",
        "count_8",
        "count_9",
        "count_6789",
        "count_89",
        *MULTILAYER_COUNT_FIELDS,
        *MM89_MULTILAYER_COUNT_FIELDS,
    )
    for field in count_fields:
        out[field] = int(row.get(field, 0) or 0)
    return out


def empty_totals() -> dict[str, int]:
    return {
        "count_6": 0,
        "count_7": 0,
        "count_8": 0,
        "count_9": 0,
        "count_6789": 0,
        "count_89": 0,
        **empty_multilayer_counts(),
        **empty_mm89_multilayer_counts(),
    }


def add_to_totals(totals: dict[str, int], row: dict[str, Any]) -> None:
    if row["status"] != "ok":
        return
    for field in totals:
        totals[field] += int(row[field])


def top_item(row: dict[str, Any]) -> tuple[int, int, int, int, float, str, str, int, int]:
    return (
        int(row["count_89"]),
        int(row["count_8"]),
        int(row["count_9"]),
        int(row["count_6789"]),
        float(row["ratio_89_over_6789"]),
        row["granule"],
        row["path"],
        int(row["count_6"]),
        int(row["count_7"]),
    )


def maybe_add_top(heap: list[tuple[Any, ...]], row: dict[str, Any], top_n: int) -> None:
    if row["status"] != "ok" or top_n <= 0 or int(row["count_89"]) <= 0:
        return
    item = top_item(row)
    if len(heap) < top_n:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)


def read_resume_rows(path: Path) -> tuple[set[str], dict[str, int], int, int, int]:
    completed: set[str] = set()
    totals = empty_totals()
    ok_count = 0
    error_count = 0
    multilayer_error_count = 0

    if not path.exists():
        return completed, totals, ok_count, error_count, multilayer_error_count

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            parsed = row_from_csv(row)
            completed.add(parsed["path"])
            if parsed["status"] == "ok":
                ok_count += 1
            else:
                error_count += 1
            if parsed["status"] == "ok" and parsed.get("multilayer_error"):
                multilayer_error_count += 1
            add_to_totals(totals, parsed)
    return completed, totals, ok_count, error_count, multilayer_error_count


def missing_csv_fields(path: Path, required_fields: list[str]) -> list[str]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
    return [field for field in required_fields if field not in fieldnames]


def write_top_csv(path: Path, heap: list[tuple[Any, ...]]) -> None:
    fields = [
        "rank",
        "granule",
        "path",
        "ratio_89_over_6789",
        "count_89",
        "count_6789",
        "count_6",
        "count_7",
        "count_8",
        "count_9",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for rank, item in enumerate(sorted(heap, reverse=True), start=1):
            count_89, c8, c9, count_6789, ratio, granule, granule_path, c6, c7 = item
            writer.writerow(
                {
                    "rank": rank,
                    "granule": granule,
                    "path": granule_path,
                    "ratio_89_over_6789": ratio_to_csv(float(ratio)),
                    "count_89": count_89,
                    "count_6789": count_6789,
                    "count_6": c6,
                    "count_7": c7,
                    "count_8": c8,
                    "count_9": c9,
                }
            )


def write_summary(path: Path, summary: dict[str, Any]) -> None:
    with path.open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")


def ensure_parent(path: Path) -> None:
    if path.parent != Path("."):
        path.parent.mkdir(parents=True, exist_ok=True)


def discover_files(input_dir: Path, pattern: str, limit: int | None) -> list[str]:
    files = sorted(str(path) for path in input_dir.glob(pattern) if path.is_file())
    if limit is not None:
        return files[:limit]
    return files


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count pixels where MM_Flag is 8/9 and 6/7/8/9, then write the "
            "granules with the largest count of 8/9 pixels."
        )
    )
    parser.add_argument(
        "input_dir",
        nargs="?",
        default="ds_output/2015",
        help="Directory containing MM NetCDF granules. Default: ds_output/2015",
    )
    parser.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help=f"Glob pattern under input_dir. Default: {DEFAULT_PATTERN}",
    )
    parser.add_argument(
        "--variable",
        default=DEFAULT_VARIABLE,
        help=f"Flag variable name. Default: {DEFAULT_VARIABLE}",
    )
    parser.add_argument(
        "--multilayer-variable",
        default=DEFAULT_MULTILAYER_VARIABLE,
        help=(
            "Multilayer confidence flag variable name. "
            f"Default: {DEFAULT_MULTILAYER_VARIABLE}"
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 8),
        help="Number of worker processes. Default: min(CPU count, 8)",
    )
    parser.add_argument("--top", type=int, default=100, help="Number of top granules to write.")
    parser.add_argument(
        "--output-prefix",
        default="mm_flag_ratio_2015",
        help="Prefix for output files. Default: mm_flag_ratio_2015",
    )
    parser.add_argument("--granule-csv", help="Per-granule CSV path.")
    parser.add_argument("--top-csv", help="Top granule CSV path.")
    parser.add_argument("--summary-json", help="Summary JSON path.")
    parser.add_argument("--limit", type=int, help="Only scan the first N sorted files.")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help="Print progress every N newly processed files. Default: 1000",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Append to an existing per-granule CSV and skip paths already present there.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=8,
        help="multiprocessing imap_unordered chunk size. Default: 8",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        print("ERROR: --workers must be at least 1", file=sys.stderr)
        return 2
    if args.chunk_size < 1:
        print("ERROR: --chunk-size must be at least 1", file=sys.stderr)
        return 2
    if args.top < 0:
        print("ERROR: --top must be nonnegative", file=sys.stderr)
        return 2

    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        print(f"ERROR: input directory not found: {input_dir}", file=sys.stderr)
        return 2

    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

    granule_csv = Path(args.granule_csv or f"{args.output_prefix}_granules.csv")
    top_csv = Path(args.top_csv or f"{args.output_prefix}_top{args.top}.csv")
    summary_json = Path(args.summary_json or f"{args.output_prefix}_summary.json")
    for output_path in (granule_csv, top_csv, summary_json):
        ensure_parent(output_path)

    files = discover_files(input_dir, args.pattern, args.limit)
    if not files:
        print(f"ERROR: no files matched {input_dir / args.pattern}", file=sys.stderr)
        return 2

    if args.resume:
        missing_fields = missing_csv_fields(granule_csv, CSV_FIELDS)
        if missing_fields:
            print(
                "ERROR: cannot --resume because the existing granule CSV is missing "
                f"required columns: {', '.join(missing_fields)}",
                file=sys.stderr,
            )
            print(
                "Run without --resume, or use a new --output-prefix/--granule-csv.",
                file=sys.stderr,
            )
            return 2

        completed, totals, ok_count, error_count, multilayer_error_count = read_resume_rows(
            granule_csv
        )
        heap = []
        if granule_csv.exists():
            with granule_csv.open(newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    maybe_add_top(heap, row_from_csv(row), args.top)
        files_to_scan = [path for path in files if path not in completed]
        csv_mode = "a"
        write_header = not granule_csv.exists() or granule_csv.stat().st_size == 0
    else:
        completed = set()
        totals = empty_totals()
        heap = []
        ok_count = 0
        error_count = 0
        multilayer_error_count = 0
        files_to_scan = files
        csv_mode = "w"
        write_header = True

    start = time.time()
    print(
        f"Matched {len(files)} files; scanning {len(files_to_scan)} "
        f"({len(completed)} already in CSV).",
        file=sys.stderr,
    )

    with granule_csv.open(csv_mode, newline="", buffering=1) as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        if write_header:
            writer.writeheader()

        if files_to_scan:
            pool = None
            if args.workers == 1:
                iterator = map(
                    process_file,
                    ((path, args.variable, args.multilayer_variable) for path in files_to_scan),
                )
            else:
                from multiprocessing import Pool

                pool = Pool(processes=args.workers)
                iterator = pool.imap_unordered(
                    process_file,
                    ((path, args.variable, args.multilayer_variable) for path in files_to_scan),
                    chunksize=args.chunk_size,
                )

            try:
                for scanned_new, row in enumerate(iterator, start=1):
                    writer.writerow(row_for_csv(row))
                    if row["status"] == "ok":
                        ok_count += 1
                    else:
                        error_count += 1
                    if row["status"] == "ok" and row.get("multilayer_error"):
                        multilayer_error_count += 1
                    add_to_totals(totals, row)
                    maybe_add_top(heap, row, args.top)

                    if args.progress_every > 0 and scanned_new % args.progress_every == 0:
                        elapsed = time.time() - start
                        rate = scanned_new / elapsed if elapsed else 0.0
                        print(
                            f"Processed {scanned_new}/{len(files_to_scan)} new files "
                            f"({rate:.2f} files/s); ok={ok_count}, errors={error_count}",
                            file=sys.stderr,
                        )
            except KeyboardInterrupt:
                if pool is not None:
                    pool.terminate()
                    pool.join()
                    pool = None
                raise
            finally:
                if pool is not None:
                    pool.close()
                    pool.join()

    elapsed = time.time() - start
    total_6789 = totals["count_6789"]
    total_89 = totals["count_89"]
    overall_ratio = total_89 / total_6789 if total_6789 else None
    multilayer_confidence_totals = {
        str(flag): totals[field]
        for flag, field in zip(MULTILAYER_CONFIDENCE_FLAGS, MULTILAYER_COUNT_FIELDS)
    }
    mm89_multilayer_confidence_totals = {
        str(flag): totals[field]
        for flag, field in zip(MULTILAYER_CONFIDENCE_FLAGS, MM89_MULTILAYER_COUNT_FIELDS)
    }

    write_top_csv(top_csv, heap)
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_dir": str(input_dir),
        "pattern": args.pattern,
        "variable": args.variable,
        "multilayer_variable": args.multilayer_variable,
        "files_matched": len(files),
        "files_already_completed": len(completed),
        "files_scanned_this_run": len(files_to_scan),
        "successful_files": ok_count,
        "failed_files": error_count,
        "multilayer_error_files": multilayer_error_count,
        "total_count_6": totals["count_6"],
        "total_count_7": totals["count_7"],
        "total_count_8": totals["count_8"],
        "total_count_9": totals["count_9"],
        "total_count_6789": total_6789,
        "total_count_89": total_89,
        "ratio_89_over_6789": overall_ratio,
        "total_mm_multilayer_confidence_flag_counts": multilayer_confidence_totals,
        "total_mm_flag_89_multilayer_confidence_flag_counts": (
            mm89_multilayer_confidence_totals
        ),
        "top_n": args.top,
        "top_sort": "count_89_desc",
        "elapsed_seconds": elapsed,
        "outputs": {
            "granule_csv": str(granule_csv),
            "top_csv": str(top_csv),
            "summary_json": str(summary_json),
        },
    }
    write_summary(summary_json, summary)

    print(f"Successful files: {ok_count}", file=sys.stderr)
    print(f"Failed files: {error_count}", file=sys.stderr)
    print(
        f"Files with MM_MultilayerConfidenceFlag read errors: {multilayer_error_count}",
        file=sys.stderr,
    )
    print(f"Total count MM_Flag in 6,7,8,9: {total_6789}", file=sys.stderr)
    print(f"Total count MM_Flag in 8,9: {total_89}", file=sys.stderr)
    for flag in MULTILAYER_CONFIDENCE_FLAGS:
        print(
            "Total count MM_MultilayerConfidenceFlag "
            f"== {flag}: {multilayer_confidence_totals[str(flag)]}",
            file=sys.stderr,
        )
    for flag in MULTILAYER_CONFIDENCE_FLAGS:
        print(
            "Total count where MM_Flag is 8 or 9 and "
            "MM_MultilayerConfidenceFlag "
            f"== {flag}: {mm89_multilayer_confidence_totals[str(flag)]}",
            file=sys.stderr,
        )
    ratio_text = "NA" if overall_ratio is None else f"{overall_ratio:.12g}"
    print(f"Overall ratio (8 or 9)/(6,7,8,9): {ratio_text}", file=sys.stderr)
    print(f"Wrote {granule_csv}", file=sys.stderr)
    print(f"Wrote {top_csv}", file=sys.stderr)
    print(f"Wrote {summary_json}", file=sys.stderr)
    return 1 if error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
