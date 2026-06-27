#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mmcth lookup tool

Non-MPI version of the MMCTH driver.

Purpose:
    Given the same query options as the processing driver, fetch the matching
    MODIS/MISR/ERA5 input-file groups from the SQLite database and print them.

This script does NOT process anything:
    - no mpi4py
    - no MainProcessor
    - no proc.run_process()

Examples:
    python lookup_mmcth_groups.py --modisid MOD021KM.A2021061.0240.061.2021061204412.hdf

    python lookup_mmcth_groups.py --modisid-file modisids.txt

    python lookup_mmcth_groups.py --orbit 112783

    python lookup_mmcth_groups.py --year 2021 --date 03-02

    python lookup_mmcth_groups.py --start-date 2021-03-01 --end-date 2021-03-03
"""

import argparse
import logging
import os
import sqlite3
import sys
import warnings
from pathlib import Path


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
log = logging.getLogger("mmcth.lookup")


# ---------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parent
sys.path.append(str(PROJECT_DIR / "src"))

from src.misc import fetch_methods


# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------
DEFAULT_DB_PATH = "/data/keeling/a/gzhao1/f/Database/inputfiles_all.sqlite"


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def print_inputfiles_block(gid, files, orbit_number, group_index=None):

    """

    Print one group as a Python-ready inputfiles list.

    Example output:

        inputfiles = ['file1', \

                      'file2']

    """

    if not isinstance(files, (list, tuple)):

        return

    print("=" * 100)

    if group_index is not None:

        print(f"Python inputfiles block for group {group_index} | gid={gid} | orbit={orbit_number}")

    else:

        print(f"Python inputfiles block | gid={gid} | orbit={orbit_number}")

    print("-" * 100)

    print("inputfiles = [", end="")

    for i, path in enumerate(files):

        path = str(path)

        if i == 0:

            prefix = ""

        else:

            prefix = " " * 13

        if i < len(files) - 1:

            print(f"{prefix}{path!r}, \\")

        else:

            print(f"{prefix}{path!r}]")

    print()
def create_connection(db_file: str):
    """
    Open a read-mostly SQLite connection and apply a few pragmas that can
    help performance for lookup workloads.
    """
    try:
        conn = sqlite3.connect(db_file)
        conn.execute("PRAGMA journal_mode=OFF;")
        conn.execute("PRAGMA synchronous=OFF;")
        conn.execute("PRAGMA temp_store=MEMORY;")
        conn.execute("PRAGMA cache_size=-200000;")
        return conn
    except sqlite3.Error as exc:
        log.error("SQLite connection error: %s", exc)
        return None


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


def fetch_groups_for_modisids(conn, modisids: list[str]):
    """
    Fetch groups for a list of MODIS IDs serially.
    """
    all_groups = {}
    all_orbit_info = {}

    n_ids = len(modisids)
    if n_ids == 0:
        return all_groups, all_orbit_info

    for i, modisid in enumerate(modisids, start=1):
        if i == 1 or i % 25 == 0 or i == n_ids:
            log.info("Fetching MODIS ID %d/%d: %s", i, n_ids, modisid)

        groups, orbit_info = fetch_methods.fetch_files_by_modisid(conn, modisid)

        if not groups:
            log.warning("No groups found for MODIS ID: %s", modisid)
            continue

        for gid, files in groups.items():
            if gid in all_groups:
                log.warning(
                    "Duplicate gid=%s encountered for MODIS ID %s; overwriting",
                    gid,
                    modisid,
                )
            all_groups[gid] = files

        if isinstance(orbit_info, dict):
            all_orbit_info.update(orbit_info)
        else:
            for gid in groups.keys():
                all_orbit_info[gid] = orbit_info

    log.info("Finished fetching %d MODIS IDs -> %d groups", n_ids, len(all_groups))
    return all_groups, all_orbit_info


def classify_file(path: str) -> str:
    """
    Simple label for common input files.
    """
    base = os.path.basename(str(path))

    if "MOD021KM" in base:
        return "MOD021KM"
    if "MOD03" in base:
        return "MOD03"
    if "MOD06" in base:
        return "MOD06_L2"
    if "TC_CLOUD" in base or "TC_CLOUD" in str(path) or "TCCloud" in str(path):
        return "MISR_TC_CLOUD"
    if "GP_GMP" in base or "GMP" in base:
        return "MISR_GMP"
    if "AGP" in base:
        return "MISR_AGP"
    if "ERA5" in base or "era5" in str(path).lower():
        return "ERA5"
    return "OTHER"


def print_group(gid, files, orbit_number, group_index=None):
    """
    Print one group cleanly.
    """
    header = f"Group {group_index}: {gid}" if group_index is not None else f"Group: {gid}"

    print("=" * 100)
    print(header)
    print(f"Orbit: {orbit_number}")
    print("-" * 100)

    if files is None:
        print("files = None")
        return

    if isinstance(files, dict):
        for key, value in files.items():
            print(f"{key}: {value}")
        return

    for i, path in enumerate(files):
        label = classify_file(path)
        print(f"[{i}] {label}")
        print(f"    {path}")

    mod021_paths = [
        str(p) for p in files
        if "MOD021KM" in os.path.basename(str(p))
    ]

    if mod021_paths:
        print("-" * 100)
        print("MOD021KM path(s):")
        for p in mod021_paths:
            print(f"    {p}")

    print()


def print_all_groups(groups, orbit_info):
    """
    Print all groups and a compact summary.
    """
    work_items = build_work_items(groups, orbit_info)

    print()
    print("Lookup summary")
    print("--------------")
    print(f"Number of groups: {len(work_items)}")
    print()

    for idx, (gid, files, orbit_number) in enumerate(work_items, start=1):
        print_group(gid, files, orbit_number, group_index=idx)

    print("=" * 100)
    print("Compact MOD021KM summary")
    print("------------------------")

    found_any_mod021 = False

    for idx, (gid, files, orbit_number) in enumerate(work_items, start=1):
        if not isinstance(files, (list, tuple)):
            continue

        for path in files:
            if "MOD021KM" in os.path.basename(str(path)):
                found_any_mod021 = True
                print(f"group {idx} | gid={gid} | orbit={orbit_number}")
                print(f"    {path}")

    if not found_any_mod021:
        print("No MOD021KM paths found in returned groups.")
    print()
    print("=" * 100)
    print("Python inputfiles block(s)")
    print("-------------------------")
    for idx, (gid, files, orbit_number) in enumerate(work_items, start=1):

        print_inputfiles_block(gid, files, orbit_number, group_index=idx)


def fetch_by_args(conn, args):
    """
    Dispatch to the same fetch_methods calls used by the original driver.
    """
    if args.modisid_file:
        modisids = read_modisid_file(args.modisid_file)
        if not modisids:
            raise SystemExit(f"No MODIS IDs found in {args.modisid_file}")
        log.info("Read %d MODIS IDs from %s", len(modisids), args.modisid_file)
        return fetch_groups_for_modisids(conn, modisids)

    if args.start_date and args.end_date:
        return fetch_methods.fetch_files_by_date_range(
            conn,
            args.start_date,
            args.end_date,
        )

    if args.start_date and args.months:
        return fetch_methods.fetch_files_next_n_months(
            conn,
            args.start_date,
            args.months,
        )

    if args.date:
        if args.year is None:
            raise ValueError("--year is required when using --date MM-DD")
        month, day = map(int, args.date.split("-"))
        stamp = f"{args.year}-{month:02d}-{day:02d}"
        return fetch_methods.fetch_files_by_date(conn, stamp)

    if args.month is not None:
        if args.year is None:
            raise ValueError("--year is required when using --month")
        return fetch_methods.fetch_files_by_month(
            conn,
            args.year,
            args.month,
        )

    if args.orbit:
        return fetch_methods.fetch_files_by_orbit(
            conn,
            args.orbit,
            False,
        )

    if args.modisid:
        return fetch_methods.fetch_files_by_modisid(
            conn,
            args.modisid,
        )

    if args.year is not None:
        return fetch_methods.fetch_files_by_year(
            conn,
            args.year,
        )

    raise ValueError(
        "Provide one of --year, --month, --date, --orbit, "
        "--modisid, --modisid-file, or date-range options."
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Lookup MODIS/MISR/ERA5 input-file groups from the MMCTH SQLite "
            "database. This script does not process files."
        )
    )

    p.add_argument(
        "--db",
        default=DEFAULT_DB_PATH,
        help=f"SQLite database path. Default: {DEFAULT_DB_PATH}",
    )

    p.add_argument("-y", "--year", type=int)
    p.add_argument("-m", "--month", type=int)
    p.add_argument("-d", "--date", type=str, help="MM-DD")
    p.add_argument("-o", "--orbit", type=str)
    p.add_argument("-i", "--modisid", type=str)
    p.add_argument("--modisid-file", type=str, help="Text file with one MODIS ID per line")

    p.add_argument("-s", "--start-date", type=str, help="Start date YYYY-MM-DD")
    p.add_argument("-e", "--end-date", type=str, help="End date YYYY-MM-DD inclusive")
    p.add_argument("-n", "--months", type=int, help="If set with --start-date, span this many months")

    p.add_argument("--debug", action="store_true", help="enable debug logging")
    p.add_argument("--quiet", action="store_true", help="only warnings and errors")

    args = p.parse_args()

    level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    logging.getLogger().setLevel(level)

    print("Lookup settings")
    print("---------------")
    print(f"Database: {args.db}")

    conn = create_connection(args.db)
    if conn is None:
        raise SystemExit(f"Could not open SQLite database: {args.db}")

    try:
        groups, orbit_info = fetch_by_args(conn, args)
    finally:
        conn.close()

    if not groups:
        print()
        print("No groups found.")
        return

    print_all_groups(groups, orbit_info)


if __name__ == "__main__":
    main()