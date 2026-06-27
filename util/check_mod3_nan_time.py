#!/usr/bin/env python3
"""
Scan one or more MOD03 files and report those with NaNs in the EV center time.

Usage examples:
    # Explicitly check two files
    python check_ev_time_nan.py /path/to/MOD03.A2015176.1650*.hdf /path/to/other.hdf

    # Recursively check all MOD03 files in a directory
    python check_ev_time_nan.py /data/gdi/satellite/TerraDataArchive/MODIS/MOD03
"""

import sys
import os
import glob
import numpy as np
from pyhdf.SD import SD, SDC

def file_contains_nan_ev_center_time(filepath: str) -> bool:
    """Return True if the EV center time SDS contains any NaN after fill replacement."""
    hdf = SD(filepath, SDC.READ)
    try:
        data = hdf.select("EV center time").get().astype(np.float64)
    finally:
        hdf.end()

    # Replace the MOD03 fill value (≤ -1e9) with NaN
    data[data <= -1.0e9] = np.nan
    return np.isnan(data).any()

def main(paths):
    # Build a flat list of .hdf files from given arguments.
    file_list = []
    for p in paths:
        if os.path.isdir(p):
            # Recursively search for MOD03 files under this directory
            file_list += glob.glob(os.path.join(p, '**', 'MOD03.*.hdf'), recursive=True)
        else:
            file_list.append(p)

    bad_files = []
    for f in file_list:
        try:
            if file_contains_nan_ev_center_time(f):
                bad_files.append(f)
        except Exception as err:
            print(f"Warning: could not inspect {f}: {err}", file=sys.stderr)

    if bad_files:
        print("Files containing NaN EV center time:")
        for bf in bad_files:
            print(bf)
    else:
        print("No files with NaN EV center time were found.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        # Default to scanning the entire MOD03 archive if no args are given
        main(['/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2015/'])