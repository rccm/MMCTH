#!/usr/bin/env python3
"""
Benchmark NetCDF I/O for reading ONLY time slice 0 (and/or valid_time slice 0)
for all 4D variables in a single file.

Measures per engine:
  - open_dataset() time (metadata only)
  - time to actually read slice 0 into memory (forces I/O)
  - RSS delta for the slice read
  - approximate GB read (from in-memory nbytes)

Example:
  python xr_io_time0.py /data/gdi/e/ERA5/multi/era5_profile_2011_03_13.nc --engines netcdf4 h5netcdf
  python xr_io_time0.py file.nc --nrepeat 3
  python xr_io_time0.py file.nc --no-decode-times
"""

from __future__ import annotations

import argparse
import gc
import os
import time
from typing import Dict, List

import xarray as xr

try:
    import psutil  # pip/conda install psutil
except Exception:
    psutil = None


def rss_mb() -> float:
    if psutil is None:
        return float("nan")
    return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)


def now() -> float:
    return time.perf_counter()


def select_4d_vars(ds: xr.Dataset) -> List[str]:
    return [k for k, da in ds.data_vars.items() if len(da.dims) == 4]


def slice0_dataset(ds: xr.Dataset, vars_4d: List[str]) -> xr.Dataset:
    """
    Build a dataset containing only time-like slice 0 for each 4D variable.
    For variables with dim 'time', selects time=0.
    For variables with dim 'valid_time', selects valid_time=0.
    If a var has both (unlikely), selects both.
    """
    out = {}
    for v in vars_4d:
        da = ds[v]
        isel_kwargs = {}
        if "time" in da.dims:
            isel_kwargs["time"] = 0
        if "valid_time" in da.dims:
            isel_kwargs["valid_time"] = 0
        out[v] = da.isel(**isel_kwargs) if isel_kwargs else da
    return xr.Dataset(out)


def force_io(ds_slice: xr.Dataset, vars_4d: List[str]) -> float:
    """
    Force actual reads. Using a reduction + .compute() ensures the backend
    must touch all values for each variable.
    Returns total bytes (nbytes) of loaded arrays.
    """
    # Trigger reads deterministically
    for v in vars_4d:
        _ = float(ds_slice[v].mean(skipna=False).compute())

    # Now ensure arrays are in memory
    ds_slice.load()

    total_nbytes = 0
    for v in vars_4d:
        data = ds_slice[v].data
        if hasattr(data, "nbytes"):
            total_nbytes += int(data.nbytes)
    return float(total_nbytes)


def bench_once(path: str, engine: str | None, decode_times: bool) -> Dict[str, float]:
    gc.collect()
    base_rss = rss_mb()

    # Open (metadata)
    t0 = now()
    ds = xr.open_dataset(path, engine=engine, decode_times=decode_times)
    t1 = now()
    rss_after_open = rss_mb()

    vars_4d = select_4d_vars(ds)

    # Slice 0 + force I/O
    gc.collect()
    rss_before = rss_mb()
    t2 = now()
    ds_slice = slice0_dataset(ds, vars_4d)
    nbytes = force_io(ds_slice, vars_4d)
    t3 = now()
    rss_after = rss_mb()

    # Cleanup
    ds.close()
    del ds_slice, ds
    gc.collect()

    return {
        "open_s": t1 - t0,
        "open_rss_mb": rss_after_open - base_rss,
        "n4d": float(len(vars_4d)),
        "slice0_s": t3 - t2,
        "slice0_rss_mb": rss_after - rss_before,
        "slice0_gb": nbytes / 1e9,
    }


def fmt(x: float) -> str:
    if x != x:  # NaN
        return "NA"
    return f"{x:.3f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument(
        "--engines",
        nargs="+",
        default=["netcdf4", "h5netcdf", None],
        help='Engines to try. Use None (or "none") for xarray autodetect.',
    )
    ap.add_argument("--nrepeat", type=int, default=1)
    ap.add_argument("--no-decode-times", action="store_true")
    args = ap.parse_args()

    decode_times = not args.no_decode_times

    print(f"File: {args.path}")
    print(f"decode_times={decode_times}")
    if psutil is None:
        print("Note: psutil not found; RSS memory numbers will be NA. Install: pip/conda install psutil")
    print()

    header = (
        "engine".ljust(12)
        + "rep".rjust(4)
        + "open_s".rjust(10)
        + "open_rssMB".rjust(12)
        + "n4d".rjust(6)
        + "slice0_s".rjust(12)
        + "slice0_rssMB".rjust(13)
        + "slice0_GB".rjust(11)
    )
    print(header)
    print("-" * len(header))

    for eng in args.engines:
        engine = None if eng in (None, "None", "none") else str(eng)
        for rep in range(1, args.nrepeat + 1):
            try:
                r = bench_once(args.path, engine=engine, decode_times=decode_times)
                line = (
                    (str(engine) if engine is not None else "default").ljust(12)
                    + f"{rep:>4d}"
                    + f"{fmt(r['open_s']):>10}"
                    + f"{fmt(r['open_rss_mb']):>12}"
                    + f"{int(r['n4d']):>6d}"
                    + f"{fmt(r['slice0_s']):>12}"
                    + f"{fmt(r['slice0_rss_mb']):>13}"
                    + f"{fmt(r['slice0_gb']):>11}"
                )
                print(line)
            except Exception as e:
                print(
                    (str(engine) if engine is not None else "default").ljust(12)
                    + f"{rep:>4d}"
                    + "  ERROR: " + repr(e)
                )


if __name__ == "__main__":
    main()