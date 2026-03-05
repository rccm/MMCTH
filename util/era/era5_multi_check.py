#!/usr/bin/env python3
import argparse
import csv
import os
import random
import sys
import traceback

import xarray as xr

# Variables you expect in multi files; adjust if yours differ
DEFAULT_VARS = ["temperature", "specific_humidity", "geopotential", "u", "v", "w"]

def safe_open(path: str):
    # Try engines in order; both can fail for corrupt filters.
    last = None
    for eng in ("h5netcdf", "netcdf4"):
        try:
            ds = xr.open_dataset(path, engine=eng, decode_times=False, mask_and_scale=False, cache=False)
            return ds, eng, None
        except Exception as e:
            last = (eng, e)
    return None, None, last

def pick_dim(da, candidates):
    for c in candidates:
        if c in da.dims:
            return c
    return None

def probe_reads(ds: xr.Dataset, var: str, probes: int):
    """Attempt many small reads at random indices; return first error (or None)."""
    if var not in ds.variables:
        return None  # missing var is not "corruption" here

    da = ds[var]

    tdim = pick_dim(da, ["time", "valid_time"])
    ldim = pick_dim(da, ["level", "pressure_level"])
    ydim = pick_dim(da, ["latitude", "lat", "y"])
    xdim = pick_dim(da, ["longitude", "lon", "x"])

    sizes = da.sizes

    # Build indexers for random point reads
    def rand_idx(d):
        n = sizes.get(d, 0)
        return 0 if n <= 1 else random.randrange(0, n)

    # If dims are not as expected, still try a generic small slice
    for p in range(probes):
        indexers = {}
        where = []

        for d in (tdim, ldim, ydim, xdim):
            if d is None:
                continue
            idx = rand_idx(d)
            indexers[d] = idx
            where.append(f"{d}={idx}")

        try:
            # This forces an actual HDF5 read
            _ = da.isel(**indexers).values
        except Exception as e:
            return ",".join(where) if where else "(unknown)", repr(e)

    # Also force a small contiguous slab (often triggers filter issues)
    try:
        slab = da
        if tdim and sizes[tdim] > 0:
            slab = slab.isel({tdim: 0})
        if ldim and sizes[ldim] > 0:
            slab = slab.isel({ldim: 0})
        if ydim and sizes[ydim] >= 4:
            slab = slab.isel({ydim: slice(0, 4)})
        if xdim and sizes[xdim] >= 4:
            slab = slab.isel({xdim: slice(0, 4)})
        _ = slab.values
    except Exception as e:
        return "slab(0/0/0:4/0:4)", repr(e)

    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--out", required=True)
    ap.add_argument("--append", action="store_true")
    ap.add_argument("--probes", type=int, default=300)
    ap.add_argument("--vars", nargs="*", default=DEFAULT_VARS)
    args = ap.parse_args()

    mode = "a" if args.append else "w"
    new_file = (mode == "w") or (not os.path.exists(args.out))

    with open(args.out, mode, newline="") as fcsv:
        w = csv.writer(fcsv)
        if new_file:
            w.writerow(["file", "var", "where", "error"])

        ds, eng, last = safe_open(args.path)
        if ds is None:
            w.writerow([args.path, "(open)", "", f"{last[0]}: {repr(last[1])}"])
            return 0

        try:
            # Probe expected vars + anything else that matches pattern (optional)
            vars_to_check = list(dict.fromkeys(args.vars))  # preserve order, unique
            for v in vars_to_check:
                err = probe_reads(ds, v, args.probes)
                if err is not None:
                    where, msg = err
                    w.writerow([args.path, v, where, msg])

        finally:
            try:
                ds.close()
            except Exception:
                pass

    return 0

if __name__ == "__main__":
    raise SystemExit(main())