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
import re
import numpy as np
from mpi4py import MPI
from get_data import MODISMISRProcessor
import glob
import xarray as xr
from datetime import datetime, timedelta
from gen_plot import plot_mm_multipanel
import matplotlib
from matplotlib.colors import BoundaryNorm, ListedColormap
matplotlib.use("Agg")   # non-interactive backend for HPC
import matplotlib.pyplot as plt
MM_DIR = '/data/gdi/f/gzhao1/mmcth/ds_output/'
x_line=[22, 24, 26, 28, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 52, 54, 56, 58, 61, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 95, 96, 99, 100, 102, 105, 107, 109, 110, 113, 115, 117, 119, 121, 123, 125, 127, 130, 131, 134, 135, 138, 140, 142, 144, 146, 148, 150, 152, 155, 157, 159, 161, 163, 165, 167, 170, 171, 174, 175, 178, 180, 182, 184, 186, 189, 191, 193, 195, 197, 199, 201, 204, 206, 208, 210, 212, 214, 217, 219, 221, 223, 225, 227, 230, 232, 234, 236, 238, 241, 243, 245, 247, 249, 251, 253, 256, 258, 260, 262, 264, 267, 269, 271, 273, 275, 277, 280, 282, 284, 286, 288, 290, 293, 295, 297, 299, 301, 304, 306, 308, 310, 312, 314, 316, 319, 321, 323, 325, 327, 329, 332, 334, 336, 338, 340, 342, 344, 347, 349, 351, 353, 355, 358, 359, 362, 364, 366, 368, 370, 372, 374, 377, 378, 381, 383, 385, 387, 389, 391, 393, 395, 398, 399, 402, 403, 406, 408, 410, 412, 414, 416, 418, 420, 422, 424, 426, 428, 430, 432, 434, 437, 438, 440, 442, 444, 447, 448, 451, 452, 454, 456, 458, 461, 462, 464, 466, 468, 470, 472, 474, 476, 478, 480, 482, 484, 486, 488, 489, 492, 494, 495, 497, 499, 501, 503, 505, 507, 508, 511, 512, 514, 516, 518, 520, 521, 523, 525, 527, 529, 531, 533, 534, 536, 538, 540, 542, 543]
y_line=[1084, 1080, 1076, 1072, 1067, 1063, 1059, 1055, 1051, 1047, 1043, 1038, 1034, 1030, 1026, 1022, 1018, 1014, 1009, 1005, 1001, 997, 993, 989, 985, 981, 976, 972, 968, 964, 960, 956, 952, 947, 943, 939, 935, 931, 927, 923, 918, 914, 910, 906, 902, 897, 893, 889, 885, 881, 877, 873, 868, 864, 860, 856, 852, 848, 844, 840, 835, 831, 827, 823, 819, 815, 811, 806, 802, 798, 794, 790, 786, 781, 777, 773, 769, 765, 761, 756, 752, 748, 744, 740, 736, 732, 727, 723, 719, 715, 711, 707, 703, 698, 694, 690, 686, 682, 678, 673, 669, 665, 661, 657, 653, 649, 644, 640, 636, 632, 628, 624, 619, 615, 611, 607, 603, 599, 594, 590, 586, 582, 578, 574, 570, 565, 561, 557, 553, 549, 545, 540, 536, 532, 528, 524, 520, 515, 511, 507, 503, 499, 495, 491, 486, 482, 478, 474, 470, 466, 461, 457, 453, 449, 445, 441, 436, 432, 428, 424, 420, 416, 411, 407, 403, 399, 395, 391, 386, 382, 378, 374, 370, 366, 361, 357, 353, 349, 345, 341, 336, 332, 328, 324, 320, 315, 311, 307, 303, 299, 295, 291, 286, 282, 278, 274, 269, 265, 261, 257, 253, 249, 245, 241, 236, 232, 228, 224, 219, 215, 211, 207, 203, 198, 194, 190, 186, 182, 178, 174, 169, 165, 161, 157, 153, 148, 144, 140, 136, 132, 127, 124, 119, 115, 111, 107, 103, 98, 94, 90, 86, 82, 77, 73, 69, 65, 61, 56, 52, 48, 44, 40, 35, 32, 27]


# Silence the cfgrib / ecCodes warning once, globally
warnings.filterwarnings(
    "ignore",
    message="Failed to load cfgrib - most likely there is a problem accessing the ecCodes library."
)
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
PROJECT_DIR = Path(__file__).resolve().parents[1]   # mmcth/
sys.path.insert(0, str(PROJECT_DIR))
from src.misc import fetch_methods

# ---------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------
def parse_modis_id(s: str) -> str:
    year = int(s[:4])
    doy  = int(s[4:7])
    hhmm = s.split(".")[1]

    hour   = int(hhmm[:2])
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
    return 

import xarray as xr
import numpy as np
from pathlib import Path

BAD_DIMS_ERA5 = {"era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"}

def merge_rad_into_mm_1km(
    mm_file: str,
    rad: dict[str, np.ndarray],
    output_flag=False,
) -> None:
    """
    Read an existing MM NetCDF (mm_file), keep ONLY variables that do NOT use
    the ERA5 10 km dims, add 1 km RGB radiance fields from `rad`, and write
    a new NetCDF.
    """

    # Open original MM dataset
    ds_mm = xr.open_dataset(mm_file)

    # --- 1) Start output with only 1 km index coordinates ---
    ds_out = xr.Dataset(
        coords={
            "cell_along_swath_1km": ds_mm["cell_along_swath_1km"],
            "cell_cross_swath_1km": ds_mm["cell_cross_swath_1km"],
        }
    )

    # You may also want to keep Lat/Lon/VZA etc. that are 1 km:
    for name, da in ds_mm.data_vars.items():
        dims = set(da.dims)
        # keep if it does NOT use ERA5 dims
        if BAD_DIMS_ERA5.isdisjoint(dims):
            # and it’s on the 1 km grid (you can relax this if needed)
            if "cell_along_swath_1km" in dims and "cell_cross_swath_1km" in dims:
                ds_out[name] = da

    # --- 2) Attach your radiance fields (1 km grid) ---
    ny = ds_out.dims["cell_along_swath_1km"]
    nx = ds_out.dims["cell_cross_swath_1km"]

    for name, arr in rad.items():
        arr = np.asarray(arr)
        if arr.shape != (ny, nx):
            raise ValueError(
                f"{name}: shape {arr.shape} != (cell_along_swath_1km={ny}, cell_cross_swath_1km={nx})"
            )

        ds_out[name] = xr.DataArray(
            arr,
            dims=("cell_along_swath_1km", "cell_cross_swath_1km"),
            attrs={
                "long_name": name,  # adjust to something nicer
                "units": "W m-2 sr-1 µm-1",
                "coordinates": "Latitude Longitude",
                "grid_mapping": "crs",
            },
        )

    # --- 3) Copy global attrs from original MM file (optional) ---
    ds_out.attrs = ds_mm.attrs.copy()

    # --- 4) Simple encoding/compression ---
    encoding = {}
    for vname, da in ds_out.data_vars.items():
        enc = {"zlib": True, "complevel": 4}
        fv = da.attrs.get("_FillValue")
        if fv is not None:
            enc["_FillValue"] = fv
        encoding[vname] = enc

    if output_flag:
        out_file = mm_file.replace("_MM_V00.nc", "_MM_RGB_V00.nc")
        ds_out.to_netcdf(out_file, format="NETCDF4_CLASSIC", encoding=encoding)
    ds_mm.close()
    return ds_out

def fetch_and_process_group(inputfile_list: list[str], orbit_number: int):
     
    # 1) Run your processor to get 1 km RGB radiances
    mod = MODISMISRProcessor(inputfile_list)
    rad = mod.mm_process()       # dict of 2D arrays on 1 km grid
    print(rad.keys())

    # 2) Extract MODIS ID from the MOD021KM filename
    mod_id = re.search(r'MOD021KM\.(A\d{7}\.\d{4})', inputfile_list[0]).group(1)
    print(mod_id)

    # 3) Find the existing MM file that has all the 1 km fields
    matches = glob.glob(f'{MM_DIR}MODMISR*{mod_id}*MM*V00.nc')
    if not matches:
        log.error("MM File does not exist")
        return None
    
    mm_file = matches[0]
    print(mm_file)
    # 4) Decide name for the output file (you can overwrite or make a new one)
    #    Example: append `_RGB` before V00.nc
    mm_path = Path(mm_file)
    

    # 5) Merge all 1 km vars from mm_file with your RGB radiances and dump out
    ds = merge_rad_into_mm_1km(str(mm_file), rad['bands'])
    m = re.search(r'_A(\d{7}\.\d{4})_(\d{4})_(O\d{6})', mm_file)
    if m:
        mod_id = m.group(1)   # '2015062.1320'
        orbit  = m.group(2)   # 'O08
    mod_id = parse_modis_id(mod_id)
    plots_dir = mm_path.parent / "plots"   # /data/gdi/f/gzhao1/mmcth/plots
    # plots_dir =  Path('/data/gdi/f/gzhao1/mmcth/ds_output/plots/')
    # plots_dir.mkdir(parents=True, exist_ok=True)  # make sure it exists
    png_path_1 = plots_dir / (mm_path.stem + "_1_hgt.png")
    cmap8 = plt.get_cmap("Set1", 8)

    flag_cmap = ListedColormap([
    "white",          # 0
    "yellow",   # 1
    "deepskyblue", # 2
    "goldenrod",      # 3
    "steelblue",      # 4
    "lime",           # 5
    "deeppink",      # 6
    "red",        # 7
    ])
    
    fig = plot_mm_multipanel(
        ds,
        varnames=[
            "bt_31",
            "MM_CloudTopHeight",
            "MODIS_CloudTopHeight",
            "MISR_CloudTopHeight",
            # "MODIS_CloudOpticalDepth",
            "MM_Flag",
        ],
        nrows=1,
        ncols=6,   # 1 RGB + 4 scalar panels = 5 total
        add_latlon_lines=True,
        show_lat_labels=True,
        show_lon_labels=True,
        lat_step=5.0,
        lon_step=10.0,
        discrete={
            "MM_Flag": {
                "categories": list(range(0, 8)),
                "cmap": flag_cmap,  # or provide "colors": [...]
                "labels": [str(i) for i in range(0, 8)],  # optional
            }
        },
        vlims={
            "bt_31": (200, 300),
            "MM_CloudTopHeight": (0, 16000),
            "MISR_CloudTopHeight": (0, 16000),
            "MODIS_CloudTopHeight": (0, 16000),
            "MM_CloudTopPressure": (100, 1000),
            "MODIS_CloudOpticalDepth": (0, 7),
            "MM_Flag": (0, 7),

        },
        flip_cbar_vars={"bt_31"},
        cmaps={"bt_31":"turbo_r"},
        titles={
            "bt_31": "BT 11µm (K)",
            "MM_CloudTopHeight": "MM Cloud Top Height (m)",
            "MM_CloudTopPressure": "MM Cloud Top Pressure (m)",
            "MISR_CloudTopHeight": "MISR Cloud Top Height (m)",
            "MODIS_CloudTopHeight":"MODIS Cloud Top Height (m)",
            "MODIS_CloudOpticalDepth": "MODIS Cloud Optical Depth",
            "MM_Flag": "MM Flag",
        },
        rgb_title=f'{mod_id} -- {orbit}',
        overlay_xy=(x_line, y_line),   # turn ON
    # overlay_xy=None,             # turn OFF
        overlay_kwargs={"color": "fuchsia", "linewidth": 2, "marker": None},
        output=png_path_1,
    
    )
    plt.close(fig)
    log.info(f"Done")
    
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
            "MM_Flag",
        ],
        cmaps={"MM_CloudTopPressure":"turbo_r","MODIS_CloudTopPressure":"turbo_r"},
        flip_cbar_vars={"MM_CloudTopPressure","MODIS_CloudTopPressure"},
        nrows=1,
        ncols=8,   # 1 RGB + 4 scalar panels = 5 total
        add_latlon_lines=True,
        show_lat_labels=True,
        show_lon_labels=True,
        lat_step=5.0,
        lon_step=10.0,
        discrete={
            "MM_Flag": {
                "categories": list(range(0, 8)),
                "cmap": flag_cmap,  # or provide "colors": [...]
                "labels": [str(i) for i in range(0, 8)],  # optional
            }
        },
        vlims={
            "MM_CloudEffectiveEmissivity": (0, 1),
            "MODIS_CloudEffectiveEmissivity": (0, 1),
            "MM_CloudTopPressure": (100, 1000),
            "MODIS_CloudTopPressure": (100, 1000),
            "MM_CloudOpticalDepth": (0, 7),
            "MODIS_CloudOpticalDepth": (0, 7),
            "MM_Flag": (0, 8),

        },
        titles={
            "bt_31": "BT 11µm (K)",
            "MM_CloudTopPressure": "MM Cloud Top Pressure (hpa)",
            "MM_CloudEffectiveEmissivity": "MM Cloud Effective Emissivity",
            "MODIS_CloudEffectiveEmissivity": "MODIS Cloud Effective Emissivity",
            "MODIS_CloudTopPressure":"MODIS Cloud Top Pressure (hpa)",
            "MODIS_CloudOpticalDepth": "MODIS Cloud Optical Depth",
            "MM_CloudOpticalDepth": "MM Cloud Optical Depth",
            "MODIS_CloudOpticalDepth": "MODIS Cloud Optical Depth",
            "MM_Flag": "MM Flag",
        },
        rgb_title=f'{mod_id} -- {orbit}',
        overlay_xy=(x_line, y_line),   # turn ON
    # overlay_xy=None,             # turn OFF
        overlay_kwargs={"color": "fuchsia", "linewidth": 2, "marker": None},
        output=png_path_2,
    )
    plt.close(fig)
    log.info(f"Done")



    png_path_3 = plots_dir / (mm_path.stem + "_3_diff.png")
    ds = ds.copy()
    ds["MM_minus_MODIS_CloudTopPressure"] = ds["MM_CloudTopPressure"] - ds["MODIS_CloudTopPressure"]
    ds["MM_minus_MODIS_CloudTopHeight"] = ds["MM_CloudTopHeight"] - ds["MODIS_CloudTopHeight"]
    ds["MM_minus_MODIS_CloudEffectiveEmissivity"] = ds["MM_CloudEffectiveEmissivity"] - ds["MODIS_CloudEffectiveEmissivity"]
    ds["MM_minus_MODIS_CloudOpticalDepth"] = ds["MM_CloudOpticalDepth"] - ds["MODIS_CloudOpticalDepth"]
    keep = ds["MM_Flag"] > 5
    for v in [
        "MM_minus_MODIS_CloudTopPressure",
        "MM_minus_MODIS_CloudTopHeight",
        "MM_minus_MODIS_CloudEffectiveEmissivity",
        "MM_minus_MODIS_CloudOpticalDepth",
    ]:
        ds[v] = ds[v].where(keep,other=np.nan)
    fig = plot_mm_multipanel(
        ds,
        varnames=[
            "MM_minus_MODIS_CloudTopHeight",
            "MM_minus_MODIS_CloudTopPressure",
            "MM_minus_MODIS_CloudEffectiveEmissivity",
            "MM_minus_MODIS_CloudOpticalDepth",
            "MM_Flag",
        ],
        cmap = "PRGn",
        #cmaps={ "MM_minus_MODIS_CloudTopPressure": "RdBu_r", },
        flip_cbar_vars={"MM_CloudTopPressure","MODIS_CloudTopPressure"},
        nrows=1,
        ncols=6,   # 1 RGB + 4 scalar panels = 5 total
        add_latlon_lines=True,
        show_lat_labels=True,
        show_lon_labels=True,
        lat_step=5.0,
        lon_step=10.0,
        discrete={
            "MM_Flag": {
                "categories": list(range(0, 8)),
                "cmap": flag_cmap,  # or provide "colors": [...]
                "labels": [str(i) for i in range(0, 8)],  # optional
            }
        },
        vlims={
            "MM_minus_MODIS_CloudTopHeight": (-2000,2000),
            "MM_minus_MODIS_CloudTopPressure": (-100,100),
            "MM_minus_MODIS_CloudEffectiveEmissivity":(-1,1),
            "MM_minus_MODIS_CloudOpticalDepth":(-2,2), 
            "MM_Flag": (0, 8),

        },
        titles={
            "MM_minus_MODIS_CloudTopHeight": "MM-MODIS CTH",
            "MM_minus_MODIS_CloudTopPressure": "MM-MODIS CTP",
            "MM_minus_MODIS_CloudEffectiveEmissivity": "MM-MODIS CEE",
            "MM_minus_MODIS_CloudOpticalDepth": "MM-MODIS OPT",
            "MM_Flag": "MM Flag",
        },
        rgb_title=f'{mod_id} -- {orbit}',
        overlay_xy=(x_line, y_line),   # turn ON
    # overlay_xy=None,             # turn OFF
        overlay_kwargs={"color": "fuchsia", "linewidth": 2, "marker": None},
        output=png_path_3,
    )
    plt.close(fig)
    log.info(f"Done")     

    return  

 
# ---------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------
def main() -> None:
    import re
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
    p.add_argument('-f', '--mmfile', type=str, help='MM Files')
 

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
    elif args.mmfile:
        m = re.search(r'_A(\d{7}\.\d{4})_(O\d{6})_', args.mmfile)
        if m: 
            modis_id=f'.A{m.group(1)}.'
        else:
            sys.exit(1)
        groups,orbit_info = fetch_methods.fetch_files_by_modisid(conn, modis_id) #.A2003001.2015.
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
        try:
            fetch_and_process_group(groups[gid], orbit_number)
        except ValueError as e:
            log.warning("Non-fatal: gid=%s skipped (%s) orbit=%s", groups[gid], e, orbit_number)
            continue
        except Exception:
            log.exception("gid=%s failed; continuing", gid)
            continue
    
    conn.close()

if __name__ == "__main__":
    main()
