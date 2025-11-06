#!/usr/bin/env python3
"""
mmcth driver – distribute MODIS/MISR/ERA5 jobs with MPI
"""

import argparse
import logging
import os
import sqlite3
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from mpi4py import MPI

# Silence the cfgrib / ecCodes warning once, globally
warnings.filterwarnings(
    "ignore",
    message="Failed to load cfgrib - most likely there is a problem accessing the ecCodes library."
)

# ---------------------------------------------------------------------
#  Logging – one global config, no passing logger objects everywhere
# ---------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,                  # flip to DEBUG if needed
    format="%(asctime)s %(levelname)s %(name)s:%(lineno)d — %(message)s"
)
log = logging.getLogger("mmcth.driver")  # this file’s logger

# ---------------------------------------------------------------------
#  MPI
# ---------------------------------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ---------------------------------------------------------------------
#  Project imports (after we’ve amended sys.path)
# ---------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parent
sys.path.append(str(PROJECT_DIR / "src"))
from src.misc import fetch_methods                       
from src.data_processors.process_main import MainProcessor
# from src.data_processors.process_wMOD06_woTC import WMOD06WoTC  
# from src.data_processors.process_woMOD06_wTC import WoMOD06WTC   


# ---------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------
def create_connection(db_file: str) -> sqlite3.Connection | None:
    try:
        return sqlite3.connect(db_file)
    except sqlite3.Error as exc:
        log.error("SQLite connection error: %s", exc)
        return None


def determine_process_class(mod06_file: str, tccloud_file: str):
    # has_mod06 = os.path.basename(mod06_file) != "X"
    # has_tc    = os.path.basename(tccloud_file) != "X"
    return MainProcessor

def fetch_and_process_group(inputfile_list: list[str], orbit_number: int):
    """Run one processor on one group of files with orbit number."""
    mod06_file, tccloud_file = inputfile_list[1], inputfile_list[3]
    proc_cls = determine_process_class(mod06_file, tccloud_file)
    if proc_cls is None:
        log.info("[rank %d] Skipping group (both MOD06 & TC-Cloud exist): %s",
                rank, os.path.basename(inputfile_list[0]))
        return
    proc_log = logging.LoggerAdapter(
        logging.getLogger(proc_cls.__name__), {"rank": rank}
    )
    start = time.time()
    
    # Pass orbit number to the processor if it accepts it
    try:
        # Try creating processor with orbit number if the constructor supports it
        proc = proc_cls(inputfile_list, orbit=str(orbit_number), logger=None)
    except TypeError:
        # Fall back to original constructor if orbit_number is not supported
        proc = proc_cls(inputfile_list, logger=None)
    proc.run_process()
    proc_log.info("Finished job for orbit %s in %.1f s", orbit_number, time.time() - start)


# ---------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------
def main() -> None:
    db_path = f"/data/keeling/a/gzhao1/f/Database/inputfiles_all.sqlite"

    p = argparse.ArgumentParser(
        description="Process MODIS/MISR/ERA5 granule groups."
    )
    p.add_argument("-y", "--year",   type=int)
    p.add_argument("-m", "--month",  type=int)
    p.add_argument("-d", "--date",   type=str, help="MM-DD")
    p.add_argument("-o", "--orbit",  type=str)
    p.add_argument("-i", "--modisid", type=str)
    p.add_argument("--debug", action="store_true", help="enable debug logging")
    p.add_argument("--quiet", action="store_true", help="only warnings and errors")
    p.add_argument('-s', '--start-date', type=str, help='Start date YYYY-MM-DD (range mode)')
    p.add_argument('-e', '--end-date', type=str, help='End date YYYY-MM-DD (inclusive, range mode)')
    p.add_argument('-n', '--months', type=int, help='If set with --start-date, span this many months (inclusive)')


    args = p.parse_args()

    conn = create_connection(db_path)
    if conn is None:
        sys.exit(1)
    conn = create_connection(db_path)
    if conn is None:
        sys.exit(1)

    # ---------------- fetch file groups ----------------
    if args.start_date and args.end_date:
        groups, orbit_info = fetch_methods.fetch_files_by_date_range(conn, args.start_date, args.end_date)
    elif args.start_date and args.months:
        groups, orbit_info = fetch_methods.fetch_files_next_n_months(conn, args.start_date, args.months)
    elif args.date:
        month, day = map(int, args.date.split("-"))
        stamp      = f"{args.year}-{month:02d}-{day:02d}"
        groups, orbit_info     = fetch_methods.fetch_files_by_date(conn, stamp)
    elif args.month:
        groups, orbit_info = fetch_methods.fetch_files_by_month(conn, args.year, args.month)
    elif args.orbit:
        groups, orbit_info = fetch_methods.fetch_files_by_orbit(conn, args.orbit, False)
    elif args.modisid:
        groups,orbit_info = fetch_methods.fetch_files_by_modisid(conn, args.modisid) #.A2003001.2015.
    else:
        groups,orbit_info = fetch_methods.fetch_files_by_year(conn, args.year)
    group_ids = list(groups.keys())
    
    level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s  %(levelname)-8s  %(name)s: %(message)s"
    )
    log.info("Found %d groups", len(group_ids))

    # ---------------- scatter work with MPI ------------
    chunks = np.array_split(group_ids, size) if rank == 0 else None
    my_ids = comm.scatter(chunks, root=0)

    for gid in my_ids:
        # Get orbit number for this group
        orbit_number = None
        if isinstance(orbit_info, dict) and gid in orbit_info:
            orbit_number = orbit_info[gid]
        elif isinstance(orbit_info, (int, str)):
            orbit_number = orbit_info
        elif isinstance(orbit_info, list) and len(orbit_info) > 0:
            # If we have a list of orbits, use the first one or find a matching one
            orbit_number = orbit_info[0]
        
        log.debug("Rank %s processing group %s (orbit: %s)", rank, gid, orbit_number)
        print(f"Rank {rank} is processing group {gid} with orbit {orbit_number}")
        print(f"Files: {groups[gid]}")
        fetch_and_process_group(groups[gid], orbit_number)
    
    conn.close()

if __name__ == "__main__":
    main()
