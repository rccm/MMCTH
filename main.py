#!/usr/bin/env python3
"""
mmcth driver – distribute MODIS/MISR/ERA5 jobs with MPI

Key improvement:
- For --modisid-file, MODIS ID lookup is parallelized across ranks.
  Each rank:
    1) receives a subset of MODIS IDs
    2) opens its own SQLite connection
    3) fetches its own file groups
    4) processes its own work immediately

This avoids the long serial rank-0 prefetch bottleneck.
"""

import argparse
import logging
import os
import sqlite3
import sys
import time
import warnings
from pathlib import Path
import gc
from mpi4py import MPI

# Silence the cfgrib / ecCodes warning once, globally
warnings.filterwarnings(
    "ignore",
    message="Failed to load cfgrib - most likely there is a problem accessing the ecCodes library."
)

# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s:%(lineno)d — %(message)s"
)
log = logging.getLogger("mmcth.driver")

# ---------------------------------------------------------------------
# MPI
# ---------------------------------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ---------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parent
sys.path.append(str(PROJECT_DIR / "src"))

from src.misc import fetch_methods
from src.data_processors.process_main import MainProcessor

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def create_connection(db_file: str):
    """
    Open a read-mostly SQLite connection and apply a few pragmas that can
    help performance for heavy lookup workloads.
    """
    try:
        conn = sqlite3.connect(db_file)
        # Read-only-ish performance pragmas
        conn.execute("PRAGMA journal_mode=OFF;")
        conn.execute("PRAGMA synchronous=OFF;")
        conn.execute("PRAGMA temp_store=MEMORY;")
        conn.execute("PRAGMA cache_size=-200000;")  # ~200 MB cache
        return conn
    except sqlite3.Error as exc:
        log.error("SQLite connection error: %s", exc)
        return None


def determine_process_class(mod06_file: str, tccloud_file: str):
    return MainProcessor


def read_modisid_file(path: str) -> list[str]:
    """Read one MODIS ID per line; ignore blank lines and comments."""
    modisids = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            modisids.append(s)
    return modisids


def split_list(items, n_ranks):
    """Split a list into n_ranks striped chunks."""
    return [items[i::n_ranks] for i in range(n_ranks)]


def build_work_items(groups, orbit_info):
    """
    Convert fetched groups/orbit info into a flat list of work items.

    Returns
    -------
    work_items : list[tuple]
        [(gid, files, orbit_number), ...]
    """
    work_items = []

    for gid, files in groups.items():
        orbit_number = None
        if isinstance(orbit_info, dict):
            orbit_number = orbit_info.get(gid)
        elif isinstance(orbit_info, (int, str)):
            orbit_number = orbit_info
        elif isinstance(orbit_info, list) and len(orbit_info) > 0:
            orbit_number = orbit_info[0]

        work_items.append((gid, files, orbit_number))

    return work_items


def split_work_items(work_items, n_ranks):
    """Split work items into n_ranks striped chunks."""
    return [work_items[i::n_ranks] for i in range(n_ranks)]


def fetch_groups_for_modisids(conn, modisids: list[str]):
    """
    Fetch groups for a subset of MODIS IDs on one MPI rank.
    """
    all_groups = {}
    all_orbit_info = {}

    n_ids = len(modisids)
    if n_ids == 0:
        return all_groups, all_orbit_info

    for i, modisid in enumerate(modisids, start=1):
        if i == 1 or i % 25 == 0 or i == n_ids:
            log.info("[rank %d] Fetching MODIS ID %d/%d: %s", rank, i, n_ids, modisid)

        groups, orbit_info = fetch_methods.fetch_files_by_modisid(conn, modisid)

        if not groups:
            log.warning("[rank %d] No groups found for MODIS ID: %s", rank, modisid)
            continue

        for gid, files in groups.items():
            if gid in all_groups:
                log.warning(
                    "[rank %d] Duplicate gid=%s encountered for MODIS ID %s; overwriting",
                    rank, gid, modisid
                )
            all_groups[gid] = files

        if isinstance(orbit_info, dict):
            all_orbit_info.update(orbit_info)
        else:
            for gid in groups.keys():
                all_orbit_info[gid] = orbit_info

    log.info("[rank %d] Finished fetching %d MODIS IDs -> %d groups",
             rank, n_ids, len(all_groups))
    return all_groups, all_orbit_info


def fetch_groups_by_modisid_file_serial(conn, modisid_file: str):
    """
    Original serial implementation. Kept for reference and for possible
    non-MPI or debugging use.
    """
    all_groups = {}
    all_orbit_info = {}

    modisids = read_modisid_file(modisid_file)
    log.info("Read %d MODIS IDs from %s", len(modisids), modisid_file)

    for i, modisid in enumerate(modisids, start=1):
        if i == 1 or i % 25 == 0 or i == len(modisids):
            log.info("Fetching MODIS ID %d/%d: %s", i, len(modisids), modisid)

        groups, orbit_info = fetch_methods.fetch_files_by_modisid(conn, modisid)

        if not groups:
            log.warning("No groups found for MODIS ID: %s", modisid)
            continue

        for gid, files in groups.items():
            if gid in all_groups:
                log.warning(
                    "Duplicate gid=%s encountered for MODIS ID %s; overwriting",
                    gid, modisid
                )
            all_groups[gid] = files

        if isinstance(orbit_info, dict):
            all_orbit_info.update(orbit_info)
        else:
            for gid in groups.keys():
                all_orbit_info[gid] = orbit_info

    log.info("Finished fetching all MODIS IDs")
    return all_groups, all_orbit_info


def fetch_and_process_group(
    inputfile_list: list[str],
    orbit_number,
    output_dir: str,
    destripe_flag: bool = True,
):
    mod06_file = inputfile_list[1]
    tccloud_file = inputfile_list[3]
    proc_cls = determine_process_class(mod06_file, tccloud_file)

    if proc_cls is None:
        log.info(
            "[rank %d] Skipping group: %s",
            rank, os.path.basename(inputfile_list[0])
        )
        return

    start = time.time()
    proc = None

    try:
        try:
            proc = proc_cls(
                inputfile_list,
                orbit=str(orbit_number),
                logger=None,
                output_dir=output_dir,
                destripe_flag=destripe_flag,
            )
        except TypeError:
            proc = proc_cls(inputfile_list, logger=None, output_dir=output_dir)

        proc.run_process(save_flag="not_debug")

        log.info(
            "[rank %d] Finished job for orbit %s in %.1f s",
            rank, orbit_number, time.time() - start
        )

    finally:
        if proc is not None and hasattr(proc, "close"):
            try:
                proc.close()
            except Exception:
                log.exception("[rank %d] Error while closing processor", rank)

        del proc
        gc.collect()

def process_work_items(work_items, output_dir: str, destripe_flag: bool = True):
    """
    Process a flat list of (gid, files, orbit_number) on the current rank.
    """
    log.info("Rank %d received %d work items", rank, len(work_items))

    for gid, files, orbit_number in work_items:
        log.debug("Rank %d processing gid=%s orbit=%s", rank, gid, orbit_number)
        print(f"Rank {rank} is processing group {gid} with orbit {orbit_number}", flush=True)
        print(f"Files: {files}", flush=True)

        try:
            fetch_and_process_group(files, orbit_number, output_dir, destripe_flag)
        except ValueError as e:
            log.warning(
                "Non-fatal: gid=%s skipped (%s) orbit=%s files=%s",
                gid, e, orbit_number, files
            )
            continue
        except Exception:
            mod021km = os.path.basename(files[0]) if len(files) > 0 else "NA"
            misr_tc = os.path.basename(files[3]) if len(files) > 3 else "NA"
            log.exception(
                "FAILED gid=%s orbit=%s modis=%s misr_tc=%s files=%s",
                gid, orbit_number, mod021km, misr_tc, files
            )
            continue


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    db_path = "/data/keeling/a/gzhao1/f/Database/inputfiles_all.sqlite"

    p = argparse.ArgumentParser(
        description="Process MODIS/MISR/ERA5 granule groups."
    )
    p.add_argument("-y", "--year", type=int)
    p.add_argument("-m", "--month", type=int)
    p.add_argument("-d", "--date", type=str, help="MM-DD")
    p.add_argument("-o", "--orbit", type=str)
    p.add_argument("-i", "--modisid", type=str)
    p.add_argument("--modisid-file", type=str, help="Text file with one MODIS ID per line")
    p.add_argument("--debug", action="store_true", help="enable debug logging")
    p.add_argument("--quiet", action="store_true", help="only warnings and errors")
    destripe_group = p.add_mutually_exclusive_group()
    destripe_group.add_argument(
        "--destripe",
        dest="destripe_flag",
        action="store_true",
        default=True,
        help="enable MODIS EDF destriping for supported emissive bands (default)",
    )
    destripe_group.add_argument(
        "--no-destripe",
        dest="destripe_flag",
        action="store_false",
        help="disable MODIS EDF destriping",
    )
    p.add_argument("-s", "--start-date", type=str, help="Start date YYYY-MM-DD (range mode)")
    p.add_argument("-e", "--end-date", type=str, help="End date YYYY-MM-DD (inclusive, range mode)")
    p.add_argument("-n", "--months", type=int, help="If set with --start-date, span this many months (inclusive)")
    p.add_argument(
        "-O",
        "--output-dir",
        default="ds_output",
        help="Directory for output NetCDF files (default: ds_output)",
    )

    args = p.parse_args()

    level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    logging.getLogger().setLevel(level)

    # -----------------------------------------------------------------
    # Fast parallel path for --modisid-file
    # -----------------------------------------------------------------
    if args.modisid_file:
        if rank == 0:
            modisids = read_modisid_file(args.modisid_file)
            if not modisids:
                raise SystemExit(f"No MODIS IDs found in {args.modisid_file}")
            log.info("Read %d MODIS IDs from %s", len(modisids), args.modisid_file)
            modisid_chunks = split_list(modisids, size)
            for i, chunk in enumerate(modisid_chunks):
                log.info("Assigned %d MODIS IDs to rank %d", len(chunk), i)
        else:
            modisid_chunks = None

        my_modisids = comm.scatter(modisid_chunks, root=0)
        log.info("Rank %d received %d MODIS IDs", rank, len(my_modisids))

        conn = create_connection(db_path)
        if conn is None:
            log.error("Rank %d could not open SQLite database: %s", rank, db_path)
            comm.Abort(1)

        try:
            my_groups, my_orbit_info = fetch_groups_for_modisids(conn, my_modisids)
        finally:
            conn.close()

        my_work_items = build_work_items(my_groups, my_orbit_info)
        log.info("Rank %d built %d work items from its MODIS IDs", rank, len(my_work_items))

        process_work_items(my_work_items, args.output_dir, args.destripe_flag)
        log.info("[rank %d] Finished all assigned work; entering final barrier", rank)
        comm.Barrier()
        log.info("[rank %d] All ranks reached final barrier; exiting cleanly", rank)
        return

    # -----------------------------------------------------------------
    # Original rank-0 setup path for other query modes
    # -----------------------------------------------------------------
    setup_ok = True
    setup_error = ""
    work_chunks = None

    if rank == 0:
        conn = None
        try:
            conn = create_connection(db_path)
            if conn is None:
                raise RuntimeError(f"Could not open SQLite database: {db_path}")

            if args.start_date and args.end_date:
                groups, orbit_info = fetch_methods.fetch_files_by_date_range(
                    conn, args.start_date, args.end_date
                )
            elif args.start_date and args.months:
                groups, orbit_info = fetch_methods.fetch_files_next_n_months(
                    conn, args.start_date, args.months
                )
            elif args.date:
                if args.year is None:
                    raise ValueError("--year is required when using --date MM-DD")
                month, day = map(int, args.date.split("-"))
                stamp = f"{args.year}-{month:02d}-{day:02d}"
                groups, orbit_info = fetch_methods.fetch_files_by_date(conn, stamp)
            elif args.month is not None:
                if args.year is None:
                    raise ValueError("--year is required when using --month")
                groups, orbit_info = fetch_methods.fetch_files_by_month(
                    conn, args.year, args.month
                )
            elif args.orbit:
                groups, orbit_info = fetch_methods.fetch_files_by_orbit(
                    conn, args.orbit, False
                )
            elif args.modisid:
                groups, orbit_info = fetch_methods.fetch_files_by_modisid(
                    conn, args.modisid
                )
            else:
                if args.year is None:
                    raise ValueError(
                        "Provide one of --year, --month, --date, --orbit, "
                        "--modisid, --modisid-file, or date-range options."
                    )
                groups, orbit_info = fetch_methods.fetch_files_by_year(
                    conn, args.year
                )

            work_items = build_work_items(groups, orbit_info)
            log.info("Found %d groups", len(work_items))

            work_chunks = split_work_items(work_items, size)
            for i, chunk in enumerate(work_chunks):
                log.info("Assigned %d groups to rank %d", len(chunk), i)

        except Exception as exc:
            setup_ok = False
            setup_error = str(exc)
            log.exception("Rank 0 setup failed")
        finally:
            if conn is not None:
                conn.close()

    setup_ok = comm.bcast(setup_ok, root=0)
    setup_error = comm.bcast(setup_error, root=0)

    if not setup_ok:
        if rank != 0:
            log.error("Aborting because rank 0 setup failed: %s", setup_error)
        sys.exit(1)

    my_chunk = comm.scatter(work_chunks, root=0)
    process_work_items(my_chunk, args.output_dir, args.destripe_flag)
    log.info("[rank %d] Finished all assigned work; entering final barrier", rank)
    comm.Barrier()
    log.info("[rank %d] All ranks reached final barrier; exiting cleanly", rank)

if __name__ == "__main__":
    main()
