#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot MMCTH MODIS/MISR/MM diagnostic panels.

Main change:
    Display MM_FailureReason instead of MM_Flag in the multipanel plots.

MM_FailureReason metadata:

    _FillValue = -128b

    flag_values =
        0, 10, 20, 21, 22, 23, 24,
        30,
        40, 41, 42,
        50, 51, 52,
        60, 61, 62,
        70, 71

    flag_meanings =
        ok
        not_attempted
        invalid_misr_ctp
        invalid_surface_pressure
        invalid_modis_radiance
        invalid_atmospheric_profile
        pressure_grid_mapping_failed
        no_vertical_search_domain
        signal_below_noise
        no_valid_ratio
        no_ratio_crossing
        candidate_not_above_misr
        candidate_outside_band_range
        no_selectable_band_pair
        candidate_not_higher_than_one_layer
        reference_one_layer_failed
        candidate_has_no_valid_ctp_digit
        ctp_ok_emissivity_denominator_bad
        ctp_ok_emissivity_ratio_bad
"""

import argparse
import logging
import os
import sqlite3
import sys
import time
import warnings
from pathlib import Path
import re
import glob
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from mpi4py import MPI

import matplotlib
matplotlib.use("Agg")   # non-interactive backend for HPC
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from get_data import MODISMISRProcessor
from gen_plot import plot_mm_multipanel


# ---------------------------------------------------------------------
#  Paths
# ---------------------------------------------------------------------
MM_DIR = "/data/gdi/f/gzhao1/mmcth/ds_output/"
MM_DIR = "/data/gdi/f/gzhao1/mmcth/test_output/"


# ---------------------------------------------------------------------
#  Optional overlay line
# ---------------------------------------------------------------------
x_line = []
y_line = []


# ---------------------------------------------------------------------
#  Silence warnings
# ---------------------------------------------------------------------
warnings.filterwarnings(
    "ignore",
    message="Failed to load cfgrib - most likely there is a problem accessing the ecCodes library."
)


# ---------------------------------------------------------------------
#  Logging
# ---------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s:%(lineno)d — %(message)s"
)
log = logging.getLogger("mmcth.driver")


# ---------------------------------------------------------------------
#  MPI
# ---------------------------------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# ---------------------------------------------------------------------
#  Project imports
# ---------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]   # mmcth/
sys.path.insert(0, str(PROJECT_DIR))

from src.misc import fetch_methods


# ---------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------
BAD_DIMS_ERA5 = {
    "era5_cell_along_swath_10km",
    "era5_cell_cross_swath_10km",
}

MM_FAILURE_REASON_CODES = (
    0, 10,
    20, 21, 22, 23, 24,
    30,
    40, 41, 42,
    50, 51, 52,
    60, 61, 62,
    70, 71,
)
MM_FAILURE_REASON_LABELS = (
    "0", "10",
    "20", "21", "22", "23", "24",
    "30",
    "40", "41", "42",
    "50", "51", "52",
    "60", "61", "62",
    "70", "71",
)
MM_FAILURE_REASON_PLOT_VAR = "MM_FailureReason_plot"
MM_MULTILAYER_CONFIDENCE_CATEGORIES = (0, 1, 2, 3, 4)
MM_MULTILAYER_CONFIDENCE_LABELS = ("0", "1", "2", "3", "4")
MM_MULTILAYER_CONFIDENCE_COLORS = (
    "white",
    "blue",
    "yellow",
    "orange",
    "red",
)


# ---------------------------------------------------------------------
#  Helper functions
# ---------------------------------------------------------------------
def parse_modis_id(s: str) -> str:
    """
    Convert YYYYDDD.HHMM, optionally prefixed with A, to a timestamp.

    Example
    -------
    2015062.1320 -> 2015-03-03 (062) 13:20
    """
    match = re.fullmatch(
        r"A?(?P<year>\d{4})(?P<doy>\d{3})\.(?P<hhmm>\d{4})",
        s,
    )
    if match is None:
        raise ValueError(f"Invalid MODIS ID: {s!r}")

    year = int(match.group("year"))
    doy = int(match.group("doy"))
    hhmm = match.group("hhmm")
    hour = int(hhmm[:2])
    minute = int(hhmm[2:])

    dt = datetime(year, 1, 1) + timedelta(days=doy - 1)
    dt = dt.replace(hour=hour, minute=minute)

    return dt.strftime(f"%Y-%m-%d ({doy:03d}) %H:%M")


def create_connection(db_file: str) -> sqlite3.Connection | None:
    try:
        return sqlite3.connect(db_file)
    except sqlite3.Error as exc:
        log.error("SQLite connection error: %s", exc)
        return None


def get_mm(mm_file: str):
    ds = xr.open_dataset(mm_file)
    return ds


def merge_rad_into_mm_1km(
    mm_file: str,
    rad: dict[str, np.ndarray],
    output_flag: bool = False,
) -> xr.Dataset:
    """
    Read an existing MM NetCDF file, keep only variables that do NOT use the
    ERA5 10 km dimensions, add 1 km RGB radiance fields from `rad`, and
    optionally write a new NetCDF.

    Parameters
    ----------
    mm_file:
        Existing MM NetCDF file.

    rad:
        Dictionary of 2D 1 km radiance arrays.

    output_flag:
        If True, write a new file with `_RGB` added before `V00.nc`.

    Returns
    -------
    xr.Dataset
        Dataset containing 1 km MM variables plus added RGB/radiance fields.
    """

    ds_mm = xr.open_dataset(mm_file)

    # Start output with 1 km coordinates.
    ds_out = xr.Dataset(
        coords={
            "cell_along_swath_1km": ds_mm["cell_along_swath_1km"],
            "cell_cross_swath_1km": ds_mm["cell_cross_swath_1km"],
        }
    )

    # Keep 1 km variables that do not use ERA5 10 km dimensions.
    for name, da in ds_mm.data_vars.items():
        dims = set(da.dims)

        if BAD_DIMS_ERA5.isdisjoint(dims):
            if (
                "cell_along_swath_1km" in dims
                and "cell_cross_swath_1km" in dims
            ):
                ds_out[name] = da

    # Attach radiance fields.
    ny = ds_out.sizes["cell_along_swath_1km"]
    nx = ds_out.sizes["cell_cross_swath_1km"]

    for name, arr in rad.items():
        arr = np.asarray(arr)

        if arr.shape != (ny, nx):
            raise ValueError(
                f"{name}: shape {arr.shape} != "
                f"(cell_along_swath_1km={ny}, cell_cross_swath_1km={nx})"
            )

        ds_out[name] = xr.DataArray(
            arr,
            dims=("cell_along_swath_1km", "cell_cross_swath_1km"),
            attrs={
                "long_name": name,
                "units": "W m-2 sr-1 µm-1",
                "coordinates": "Latitude Longitude",
                "grid_mapping": "crs",
            },
        )

    # Copy global attributes.
    ds_out.attrs = ds_mm.attrs.copy()

    # Encoding/compression.
    encoding = {}
    for vname, da in ds_out.data_vars.items():
        enc = {
            "zlib": True,
            "complevel": 4,
        }

        fv = da.attrs.get("_FillValue")
        if fv is not None:
            enc["_FillValue"] = fv

        encoding[vname] = enc

    if output_flag:
        out_file = mm_file.replace("_MM_V00.nc", "_MM_RGB_V00.nc")
        ds_out.to_netcdf(
            out_file,
            format="NETCDF4_CLASSIC",
            encoding=encoding,
        )

    ds_mm.close()
    return ds_out


def build_mm_failure_reason_cmap():
    """Build an evenly spaced categorical scale for MM_FailureReason."""

    failure_categories = list(range(len(MM_FAILURE_REASON_CODES)))
    failure_labels = list(MM_FAILURE_REASON_LABELS)

    # tab20 is designed for categorical data. Reserve green for success and
    # remove its green pair so failure categories cannot look like code 0.
    tab20_colors = [
        color
        for index, color in enumerate(plt.get_cmap("tab20").colors)
        if index not in (4, 5)
    ]
    failure_colors = [
        (0.0, 0, 0.0, 1.0),  # 0: ok
        (1.0, 1.0, 1.0, 1.0),   # 10: not attempted
        *[(*color, 1.0) for color in tab20_colors[:17]],
    ]

    failure_cmap = ListedColormap(
        failure_colors,
        name="mm_failure_reason",
    )
    # Keep fill/unknown pixels distinct from white "not attempted" pixels.
    failure_cmap.set_bad(color="#bdbdbd")

    return failure_categories, failure_labels, failure_cmap


def mask_mm_failure_reason_fillvalue(ds: xr.Dataset) -> xr.Dataset:
    """Mask the MM_FailureReason fill value before category conversion."""

    if "MM_FailureReason" in ds:
        source = ds["MM_FailureReason"]
        fill_value = source.encoding.get(
            "_FillValue",
            source.attrs.get("_FillValue", -128),
        )
        ds["MM_FailureReason"] = ds["MM_FailureReason"].where(
            ds["MM_FailureReason"] != fill_value
        )
    else:
        log.warning("MM_FailureReason not found in dataset.")

    return ds


def add_mm_failure_reason_plot_variable(ds: xr.Dataset) -> xr.Dataset:
    """Map sparse failure codes to evenly spaced plotting categories."""

    if "MM_FailureReason" not in ds:
        return ds

    source = ds["MM_FailureReason"]
    source_values = np.asarray(source.values)
    category_indices = np.full(source_values.shape, np.nan, dtype=np.float32)

    for category_index, code in enumerate(MM_FAILURE_REASON_CODES):
        category_indices[source_values == code] = category_index

    known_codes = np.isin(source_values, MM_FAILURE_REASON_CODES)
    unknown_codes = np.unique(
        source_values[np.isfinite(source_values) & ~known_codes]
    )
    if unknown_codes.size:
        log.warning(
            "Unknown MM_FailureReason codes will be shown as missing: %s",
            unknown_codes.tolist(),
        )

    ds = ds.copy()
    ds[MM_FAILURE_REASON_PLOT_VAR] = xr.DataArray(
        category_indices,
        coords=source.coords,
        dims=source.dims,
        attrs={
            "long_name": "MM two-layer retrieval failure reason",
            "source_variable": "MM_FailureReason",
        },
    )
    return ds


def fetch_and_process_group(inputfile_list: list[str], orbit_number: int):
    """
    Fetch/process one MODIS/MISR/ERA5 group, merge RGB/radiance data into
    MM output, and generate diagnostic plots.
    """

    # -----------------------------------------------------------------
    #  1. Run processor to get 1 km RGB/radiance bands
    # -----------------------------------------------------------------
    mod = MODISMISRProcessor(inputfile_list)
    rad = mod.mm_process()

    print(rad.keys())

    # -----------------------------------------------------------------
    #  2. Extract MODIS ID from MOD021KM filename
    # -----------------------------------------------------------------
    mod_id_match = re.search(
        r"MOD021KM\.(A\d{7}\.\d{4})",
        inputfile_list[0],
    )

    if mod_id_match is None:
        log.error("Could not extract MODIS ID from %s", inputfile_list[0])
        return None

    mod_id = mod_id_match.group(1)
    print(mod_id)

    # -----------------------------------------------------------------
    #  3. Find existing MM file
    # -----------------------------------------------------------------
    matches = glob.glob(f"{MM_DIR}MODMISR*{mod_id}*MM*V00.nc")

    if not matches:
        log.error("MM file does not exist for MODIS ID %s", mod_id)
        return None

    mm_file = matches[0]
    print(mm_file)

    mm_path = Path(mm_file)

    # -----------------------------------------------------------------
    #  4. Merge MM 1 km fields with RGB/radiance bands
    # -----------------------------------------------------------------
    ds = merge_rad_into_mm_1km(str(mm_file), rad["bands"])

    # Mask fill value for MM_FailureReason before plotting.
    ds = mask_mm_failure_reason_fillvalue(ds)
    ds = add_mm_failure_reason_plot_variable(ds)

    # -----------------------------------------------------------------
    #  5. Parse MODIS ID and orbit for plot title
    # -----------------------------------------------------------------
    m = re.search(
        r"_A(\d{7}\.\d{4})(?:_\d{4})?_(O\d{6})(?:_|\.|$)",
        mm_path.name,
    )

    if m:
        mod_id = m.group(1)
        orbit = m.group(2)
    else:
        log.warning(
            "Could not parse title metadata from %s; using input MODIS ID "
            "and orbit %s",
            mm_path.name,
            orbit_number,
        )
        orbit = str(orbit_number)

    mod_id = parse_modis_id(mod_id)

    plots_dir = mm_path.parent / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------
    #  6. Build categorical colormaps
    # -----------------------------------------------------------------
    failure_categories, failure_labels, failure_cmap = build_mm_failure_reason_cmap()

    # -----------------------------------------------------------------
    #  7. Plot 1: height / failure reason / multilayer confidence
    # -----------------------------------------------------------------
    png_path_1 = plots_dir / (mm_path.stem + "_1_hgt.png")

    fig = plot_mm_multipanel(
        ds,
        varnames=[
            "bt_31",
            "MM_CloudTopHeight",
            "MODIS_CloudTopHeight",
            "MISR_CloudTopHeight",
            MM_FAILURE_REASON_PLOT_VAR,
            "MM_MultilayerConfidenceFlag",
        ],
        nrows=1,
        ncols=7,   # 1 RGB + 6 scalar panels
        add_latlon_lines=True,
        show_lat_labels=True,
        show_lon_labels=True,
        lat_step=5.0,
        lon_step=10.0,
        discrete={
            MM_FAILURE_REASON_PLOT_VAR: {
                "categories": failure_categories,
                "cmap": failure_cmap,
                "labels": failure_labels,
            },
            "MM_MultilayerConfidenceFlag": {
                "categories": MM_MULTILAYER_CONFIDENCE_CATEGORIES,
                "colors": MM_MULTILAYER_CONFIDENCE_COLORS,
                "labels": MM_MULTILAYER_CONFIDENCE_LABELS,
            },
        },
        vlims={
            "bt_31": (200, 300),
            "MM_CloudTopHeight": (0, 16000),
            "MISR_CloudTopHeight": (0, 16000),
            "MODIS_CloudTopHeight": (0.01, 16000),
            "MM_CloudTopPressure": (100, 1000),
            "MODIS_CloudOpticalDepth": (0, 7),
        },
        flip_cbar_vars={
            "bt_31",
        },
        cmaps={
            "bt_31": "turbo_r",
        },
        titles={
            "bt_31": "BT 11µm (K)",
            "MM_CloudTopHeight": "MM Cloud Top Height (m)",
            "MM_CloudTopPressure": "MM Cloud Top Pressure (hPa)",
            "MISR_CloudTopHeight": "MISR Cloud Top Height (m)",
            "MODIS_CloudTopHeight": "MODIS Cloud Top Height (m)",
            "MODIS_CloudOpticalDepth": "MODIS Cloud Optical Depth",
            MM_FAILURE_REASON_PLOT_VAR: "MM Failure Reason",
            "MM_MultilayerConfidenceFlag": "MM Multilayer Confidence",
        },
        rgb_title=f"{mod_id} -- {orbit}",
        overlay_xy=(x_line, y_line),
        overlay_kwargs={
            "color": "fuchsia",
            "linewidth": 2,
            "marker": None,
        },
        output=png_path_1,
    )

    plt.close(fig)
    log.info("Done writing %s", png_path_1)

    # -----------------------------------------------------------------
    #  8. Plot 2: pressure / emissivity / optical depth / failure reason
    # -----------------------------------------------------------------
    png_path_2 = plots_dir / (mm_path.stem + "_2_ctp.png")

    fig = plot_mm_multipanel(
        ds,
        varnames=[
            "MM_CloudTopPressure",
            "MODIS_CloudTopPressure",
            "MM_CloudEffectiveEmissivity",
            "MODIS_CloudEffectiveEmissivity",
            "MM_CloudOpticalDepth",
            "MODIS_CloudOpticalDepth",
            MM_FAILURE_REASON_PLOT_VAR,
        ],
        cmaps={
            "MM_CloudTopPressure": "turbo_r",
            "MODIS_CloudTopPressure": "turbo_r",
        },
        flip_cbar_vars={
            "MM_CloudTopPressure",
            "MODIS_CloudTopPressure",
        },
        nrows=1,
        ncols=8,   # 1 RGB + 7 scalar panels
        add_latlon_lines=True,
        show_lat_labels=True,
        show_lon_labels=True,
        lat_step=5.0,
        lon_step=10.0,
        discrete={
            MM_FAILURE_REASON_PLOT_VAR: {
                "categories": failure_categories,
                "cmap": failure_cmap,
                "labels": failure_labels,
            },
        },
        vlims={
            "MM_CloudEffectiveEmissivity": (0, 1),
            "MODIS_CloudEffectiveEmissivity": (0, 1),
            "MM_CloudTopPressure": (100, 1000),
            "MODIS_CloudTopPressure": (100, 1000),
            "MM_CloudOpticalDepth": (0, 7),
            "MODIS_CloudOpticalDepth": (0, 7),
        },
        titles={
            "bt_31": "BT 11µm (K)",
            "MM_CloudTopPressure": "MM Cloud Top Pressure (hPa)",
            "MM_CloudEffectiveEmissivity": "MM Cloud Effective Emissivity",
            "MODIS_CloudEffectiveEmissivity": "MODIS Cloud Effective Emissivity",
            "MODIS_CloudTopPressure": "MODIS Cloud Top Pressure (hPa)",
            "MODIS_CloudOpticalDepth": "MODIS Cloud Optical Depth",
            "MM_CloudOpticalDepth": "MM Cloud Optical Depth",
            MM_FAILURE_REASON_PLOT_VAR: "MM Failure Reason",
        },
        rgb_title=f"{mod_id} -- {orbit}",
        overlay_xy=(x_line, y_line),
        overlay_kwargs={
            "color": "fuchsia",
            "linewidth": 2,
            "marker": None,
        },
        output=png_path_2,
    )

    plt.close(fig)
    log.info("Done writing %s", png_path_2)

    # -----------------------------------------------------------------
    #  9. Plot 3: MM - MODIS differences / failure reason
    # -----------------------------------------------------------------
    png_path_3 = plots_dir / (mm_path.stem + "_3_diff.png")

    ds = ds.copy()

    ds["MM_minus_MODIS_CloudTopPressure"] = (
        ds["MM_CloudTopPressure"] - ds["MODIS_CloudTopPressure"]
    )
    ds["MM_minus_MODIS_CloudTopHeight"] = (
        ds["MM_CloudTopHeight"] - ds["MODIS_CloudTopHeight"]
    )
    ds["MM_minus_MODIS_CloudEffectiveEmissivity"] = (
        ds["MM_CloudEffectiveEmissivity"] - ds["MODIS_CloudEffectiveEmissivity"]
    )
    ds["MM_minus_MODIS_CloudOpticalDepth"] = (
        ds["MM_CloudOpticalDepth"] - ds["MODIS_CloudOpticalDepth"]
    )

    # Keep this as MM_Flag unless you want to redefine the difference-mask logic.
    # This only controls which pixels are shown in the difference panels.
    if "MM_Flag" in ds:
        keep = ds["MM_Flag"] > 5
    else:
        log.warning("MM_Flag not found. Difference panels will not be MM_Flag-masked.")
        keep = xr.ones_like(ds["MM_minus_MODIS_CloudTopPressure"], dtype=bool)

    for v in [
        "MM_minus_MODIS_CloudTopPressure",
        "MM_minus_MODIS_CloudTopHeight",
        "MM_minus_MODIS_CloudEffectiveEmissivity",
        "MM_minus_MODIS_CloudOpticalDepth",
    ]:
        ds[v] = ds[v].where(keep, other=np.nan)

    fig = plot_mm_multipanel(
        ds,
        varnames=[
            "MM_minus_MODIS_CloudTopHeight",
            "MM_minus_MODIS_CloudTopPressure",
            "MM_minus_MODIS_CloudEffectiveEmissivity",
            "MM_minus_MODIS_CloudOpticalDepth",
            MM_FAILURE_REASON_PLOT_VAR,
        ],
        cmap="PRGn",
        flip_cbar_vars={
            "MM_CloudTopPressure",
            "MODIS_CloudTopPressure",
        },
        nrows=1,
        ncols=6,   # 1 RGB + 5 scalar panels
        add_latlon_lines=True,
        show_lat_labels=True,
        show_lon_labels=True,
        lat_step=5.0,
        lon_step=10.0,
        discrete={
            MM_FAILURE_REASON_PLOT_VAR: {
                "categories": failure_categories,
                "cmap": failure_cmap,
                "labels": failure_labels,
            },
        },
        vlims={
            "MM_minus_MODIS_CloudTopHeight": (-2000, 2000),
            "MM_minus_MODIS_CloudTopPressure": (-100, 100),
            "MM_minus_MODIS_CloudEffectiveEmissivity": (-1, 1),
            "MM_minus_MODIS_CloudOpticalDepth": (-2, 2),
        },
        titles={
            "MM_minus_MODIS_CloudTopHeight": "MM - MODIS CTH",
            "MM_minus_MODIS_CloudTopPressure": "MM - MODIS CTP",
            "MM_minus_MODIS_CloudEffectiveEmissivity": "MM - MODIS CEE",
            "MM_minus_MODIS_CloudOpticalDepth": "MM - MODIS OPT",
            MM_FAILURE_REASON_PLOT_VAR: "MM Failure Reason",
        },
        rgb_title=f"{mod_id} -- {orbit}",
        overlay_xy=(x_line, y_line),
        overlay_kwargs={
            "color": "fuchsia",
            "linewidth": 2,
            "marker": None,
        },
        output=png_path_3,
    )

    plt.close(fig)
    log.info("Done writing %s", png_path_3)

    return None


# ---------------------------------------------------------------------
#  Main
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
    p.add_argument("--debug", action="store_true", help="enable debug logging")
    p.add_argument("--quiet", action="store_true", help="only warnings and errors")
    p.add_argument("-s", "--start-date", type=str, help="Start date YYYY-MM-DD")
    p.add_argument("-e", "--end-date", type=str, help="End date YYYY-MM-DD inclusive")
    p.add_argument(
        "-n",
        "--months",
        type=int,
        help="If set with --start-date, span this many months inclusive",
    )
    p.add_argument("-f", "--mmfile", type=str, help="MM file")

    args = p.parse_args()

    # Logging level from args.
    level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING

    logging.getLogger().setLevel(level)

    conn = create_connection(db_path)
    if conn is None:
        sys.exit(1)

    # -----------------------------------------------------------------
    #  Fetch file groups
    # -----------------------------------------------------------------
    if args.start_date and args.end_date:
        groups, orbit_info = fetch_methods.fetch_files_by_date_range(
            conn,
            args.start_date,
            args.end_date,
        )

    elif args.start_date and args.months:
        groups, orbit_info = fetch_methods.fetch_files_next_n_months(
            conn,
            args.start_date,
            args.months,
        )

    elif args.date:
        if args.year is None:
            log.error("--year is required when using --date MM-DD")
            sys.exit(1)

        month, day = map(int, args.date.split("-"))
        stamp = f"{args.year}-{month:02d}-{day:02d}"

        groups, orbit_info = fetch_methods.fetch_files_by_date(conn, stamp)

    elif args.month:
        if args.year is None:
            log.error("--year is required when using --month")
            sys.exit(1)

        groups, orbit_info = fetch_methods.fetch_files_by_month(
            conn,
            args.year,
            args.month,
        )

    elif args.orbit:
        groups, orbit_info = fetch_methods.fetch_files_by_orbit(
            conn,
            args.orbit,
            False,
        )

    elif args.modisid:
        groups, orbit_info = fetch_methods.fetch_files_by_modisid(
            conn,
            args.modisid,
        )

    elif args.mmfile:
        m = re.search(r"_A(\d{7}\.\d{4})_(O\d{6})_", args.mmfile)

        if m:
            modis_id = f".A{m.group(1)}."
        else:
            log.error("Could not parse MODIS ID from mmfile: %s", args.mmfile)
            sys.exit(1)

        groups, orbit_info = fetch_methods.fetch_files_by_modisid(
            conn,
            modis_id,
        )

    else:
        if args.year is None:
            log.error(
                "No valid selection supplied. Use --year, --date, --month, "
                "--orbit, --modisid, --start-date/--end-date, or --mmfile."
            )
            sys.exit(1)

        groups, orbit_info = fetch_methods.fetch_files_by_year(conn, args.year)

    group_ids = list(groups.keys())

    log.info("Found %d groups", len(group_ids))

    # -----------------------------------------------------------------
    #  Scatter work with MPI
    # -----------------------------------------------------------------
    chunks = np.array_split(group_ids, size) if rank == 0 else None
    my_ids = comm.scatter(chunks, root=0)

    for gid in my_ids:
        orbit_number = None

        if isinstance(orbit_info, dict) and gid in orbit_info:
            orbit_number = orbit_info[gid]
        elif isinstance(orbit_info, (int, str)):
            orbit_number = orbit_info
        elif isinstance(orbit_info, list) and len(orbit_info) > 0:
            orbit_number = orbit_info[0]

        log.debug(
            "Rank %s processing group %s orbit=%s",
            rank,
            gid,
            orbit_number,
        )

        print(f"Rank {rank} is processing group {gid} with orbit {orbit_number}")
        print(f"Files: {groups[gid]}")

        try:
            fetch_and_process_group(groups[gid], orbit_number)

        except ValueError as e:
            log.warning(
                "Non-fatal: gid=%s skipped: %s orbit=%s",
                gid,
                e,
                orbit_number,
            )
            continue

        except Exception:
            log.exception("gid=%s failed; continuing", gid)
            continue

    conn.close()


if __name__ == "__main__":
    main()
