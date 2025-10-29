#!/usr/bin/env python3
# mpi_merge_era5.py

import os,shutil
import sys
import time
from pathlib import Path
from datetime import datetime, timedelta
from mpi4py import MPI
import subprocess as sp

# Configuration
profile_dir = "/data/gdi/e/ERA5/multi"
vwind_dir   = "/data/gdi/e/ERA5/vwind"
output_dir  = "/data/gdi/e/ERA5/merged"
temp_dir    = "/data/gdi/e/ERA5/temp"

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def run(cmd, env=None):
    """Run a shell command robustly and return (stdout, stderr). Raise on error."""
    res = sp.run(cmd, shell=True, capture_output=True, text=True, env=env, check=True)
    return res.stdout, res.stderr

def ensure_exists(path, tries=2, delay=0.3):
    """Poll for file existence to ride out NFS metadata lag."""
    p = Path(path)
    for i in range(tries):
        if p.exists():
            return True
        time.sleep(delay)
    return p.exists()

def process_date(date_str, temp_dir_rank, base_env):
    dt = datetime.strptime(date_str, "%Y-%m-%d")
    year = dt.year
    month_str = f"{dt.month:02d}"
    day_str   = f"{dt.day:02d}"

    profile_file    = f"{profile_dir}/era5_profile_{year}_{month_str}_{day_str}.nc"
    vwind_file      = f"{vwind_dir}/era5_vwind_{year}_{month_str}_{day_str}.nc"
    temp_vwind_file = f"{temp_dir_rank}/temp_vwind_{year}_{month_str}_{day_str}.nc"
    output_file     = f"{output_dir}/era5_profile_{year}_{month_str}_{day_str}.nc"

    # Check inputs early
    if not os.path.exists(profile_file):
        return f"[{date_str}] SKIP: missing profile file: {profile_file}"
    if not os.path.exists(vwind_file):
        return f"[{date_str}] SKIP: missing vwind file: {vwind_file}"

    try:
        # Ensure dirs
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path(temp_dir_rank).mkdir(parents=True, exist_ok=True)

        # Make per-rank TMPDIR for NCO/tools (avoids internal temp collisions)
        env = dict(base_env)
        env["TMPDIR"] = temp_dir_rank

        # 1) Copy profile -> output
        run(f"cp -f {profile_file} {output_file}", env=env)


        # 2) Cast vwind coords to float, write to per-rank temp; force overwrite, netCDF4, compression
        ncap_cmd = (
            f"ncap2 -O -4 -L 1 "
            f"-s 'latitude=float(latitude);longitude=float(longitude);' "
            f"{vwind_file} {temp_vwind_file}"
        )
        try:
            run(ncap_cmd, env=env)
        except sp.CalledProcessError as e:
            return (f"[{date_str}] ERROR: ncap2 failed\n"
                    f"CMD: {ncap_cmd}\nSTDERR:\n{e.stderr}")

        # 2b) Wait for the temp file to materialize on NFS
        if not ensure_exists(temp_vwind_file, tries=4, delay=0.25):
            return (f"[{date_str}] ERROR: ncap2 output not found after creation: {temp_vwind_file}")

        # 3) Append temp to output
        ncks_cmd = f"ncks -A -4 -L 1 {temp_vwind_file} {output_file}"
        try:
            run(ncks_cmd, env=env)
        except sp.CalledProcessError as e:
            return (f"[{date_str}] ERROR: ncks failed\n"
                    f"CMD: {ncks_cmd}\nSTDERR:\n{e.stderr}")

        # 4) Cleanup
        try:
            os.remove(temp_vwind_file)
        except FileNotFoundError:
            pass

        return f"[{date_str}] OK"
    except Exception as e:
        return f"[{date_str}] ERROR: {e}"

def main():
    # Only rank 0 generates the date list
    if rank == 0:
        start_date = datetime(2019, 6, 3)
        end_date   = datetime(2019, 6, 3)

        date_list = []
        cur = start_date
        while cur <= end_date:
            date_list.append(cur.strftime("%Y-%m-%d"))
            cur += timedelta(days=1)

        # round-robin chunking
        chunks = [[] for _ in range(size)]
        for i, d in enumerate(date_list):
            chunks[i % size].append(d)
    else:
        chunks = None

    local_dates = comm.scatter(chunks, root=0)

    # Per-rank scratch + env
    temp_dir_rank = os.path.join(temp_dir, f"r{rank}")
    base_env = os.environ.copy()

    results = []
    for d in local_dates:
        results.append(process_date(d, temp_dir_rank, base_env))

    all_results = comm.gather(results, root=0)

    if rank == 0:
        for group in all_results:
            for line in group:
                print(line)

if __name__ == "__main__":
    main()
