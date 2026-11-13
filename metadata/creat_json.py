#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a minimal-plus UMM-G v1.6.5 CMR JSON sidecar for a MODMISR_L2_CP_V00
NetCDF granule.

Example usage:

    python creat_json.py \
      --mmcth ./temp/MODMISR_L2_CP_V00_061_A2015112.1010_0422_O081609_MM.nc

Expected filename convention:

    MODMISR_L2_CP_V00_061_Ayyyyddd.hhmm_MMDD_Ooooooo_MM.nc

Expected JSON sidecar filename:

    MODMISR_L2_CP_V00_061_Ayyyyddd.hhmm_MMDD_Ooooooo_MM.nc.cmr.json

This script generates UMM-G granule metadata, not UMM-C collection metadata.

Key choices:
- TemporalExtent is the actual sensing beginning/end time of this individual granule.
- DataGranule.ProductionDateTime describes the derived MODMISR NetCDF product,
  using the NetCDF global attribute date_created if available.
- DayNightFlag and sensing time are read from the corresponding local MOD03 file.
- SpatialExtent is read from the MODMISR NetCDF global bounding attributes.

Assumptions:
- The MODMISR filename contains _Ayyyyddd.hhmm_.
- The MODMISR granule corresponds 1:1 to a MOD03 granule.
- MOD03 root is:
  /data/gdi/satellite/TerraDataArchive/MODIS/MOD03/<YYYY>/<DDD>/MOD03.Ayyyyddd.hhmm*.hdf
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone, timedelta
from pathlib import Path


# ---------- tiny utils ----------
def iso_z(dtobj: datetime) -> str:
    """
    Convert a datetime object to UTC ISO-8601 with trailing Z and no microseconds.
    """
    if dtobj.tzinfo is None:
        dtobj = dtobj.replace(tzinfo=timezone.utc)
    return (
        dtobj.astimezone(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def md5sum(path, buf=1 << 20):
    """
    Compute MD5 checksum for a file.
    """
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()


def provider_dates_from_mtime(path):
    """
    Use filesystem modification time for CMR ProviderDates Create/Update.

    Note:
    ProviderDates are archive/provider metadata timestamps, not sensing time.
    """
    t = datetime.fromtimestamp(os.path.getmtime(path), tz=timezone.utc)
    ts = iso_z(t)
    return [
        {"Date": ts, "Type": "Create"},
        {"Date": ts, "Type": "Update"},
    ]


def parse_mmcth_stamp(basename):
    """
    Parse _Ayyyyddd.hhmm_ from the MODMISR filename.

    Example:
        MODMISR_L2_CP_V00_061_A2015176.1650_0625_O082545_MM.nc
        -> 2015, doy 176, 16:50 UTC
    """
    m = re.search(r"_A(\d{4})(\d{3})\.(\d{2})(\d{2})_", basename)
    if not m:
        raise ValueError(f"Cannot find _Ayyyyddd.hhmm_ in: {basename}")

    year, doy, hh, minute = map(int, m.groups())
    t0 = (
        datetime(year, 1, 1, tzinfo=timezone.utc)
        + timedelta(days=doy - 1, hours=hh, minutes=minute)
    )
    return year, doy, f"{hh:02d}{minute:02d}", t0


# ---------- very fast MOD03 lookup ----------
def _find_by_prefix(dirpath, prefix):
    """
    Find first regular .hdf file in dirpath starting with prefix.
    """
    try:
        with os.scandir(dirpath) as it:
            for e in it:
                if not e.is_file():
                    continue
                name = e.name
                if name.startswith(prefix) and name.endswith(".hdf"):
                    return os.path.join(dirpath, name)
    except FileNotFoundError:
        return None
    return None


def find_mod03(mod03_root, t0):
    """
    Locate corresponding MOD03 granule using exact time first, then small offsets.

    This handles rare cases where the filename time and MOD03 directory minute
    may differ slightly.
    """
    for off in (0, 5, -5, 10, -10):
        tt = t0 + timedelta(minutes=off)
        y = f"{tt.year:04d}"
        d = f"{tt.timetuple().tm_yday:03d}"
        hhmm = f"{tt.hour:02d}{tt.minute:02d}"

        dpath = os.path.join(mod03_root, y, d)
        prefix = f"MOD03.A{y}{d}.{hhmm}"

        p = _find_by_prefix(dpath, prefix)
        if p:
            return p

    return None


# ---------- read ECS PVL bits from MOD03 ----------
def _read_core_archive_text(mod03_path):
    """
    Read CoreMetadata.0 and ArchiveMetadata.0 from MOD03.

    First tries pyhdf. If pyhdf is unavailable, falls back to ncdump -h.
    """
    if not mod03_path:
        return ""

    try:
        from pyhdf.SD import SD, SDC

        sd = SD(mod03_path, SDC.READ)
        attrs = sd.attributes()
        txt = f"{attrs.get('CoreMetadata.0', '')}\n{attrs.get('ArchiveMetadata.0', '')}"
        sd.end()

        if txt.strip():
            return txt
    except Exception:
        pass

    try:
        return subprocess.check_output(
            ["ncdump", "-h", mod03_path],
            text=True,
            errors="ignore",
        )
    except Exception:
        return ""


_PVL_FLAGS = re.IGNORECASE | re.DOTALL


def pvl_obj_num(txt: str, key: str):
    """
    Find numeric VALUE inside an ECS PVL OBJECT block:

        OBJECT = <KEY>
          ...
          VALUE = <num>
          ...
        END_OBJECT = <KEY>

    Falls back to KEY = <num> if present.
    """
    if not txt:
        return None

    pat = (
        rf"OBJECT\s*=\s*{re.escape(key)}\s.*?"
        rf"VALUE\s*=\s*([\-0-9\.Ee\+]+)\s.*?"
        rf"END_OBJECT\s*=\s*{re.escape(key)}"
    )
    m = re.search(pat, txt, _PVL_FLAGS)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            return None

    m = re.search(rf"{re.escape(key)}\s*=\s*([\-0-9\.Ee\+]+)", txt, re.IGNORECASE)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            return None

    return None


def pvl_obj_str(txt: str, key: str):
    """
    Find quoted string VALUE inside an ECS PVL OBJECT block.

    Falls back to KEY = "str" if present.
    """
    if not txt:
        return None

    pat = (
        rf"OBJECT\s*=\s*{re.escape(key)}\s.*?"
        rf"VALUE\s*=\s*\"([^\"]+)\"\s.*?"
        rf"END_OBJECT\s*=\s*{re.escape(key)}"
    )
    m = re.search(pat, txt, _PVL_FLAGS)
    if m:
        return m.group(1)

    m = re.search(rf"{re.escape(key)}\s*=\s*\"([^\"]+)\"", txt, re.IGNORECASE)
    return m.group(1) if m else None


def extract_temporal_spatial_daynight(mod03_path, fallback_center_dt):
    """
    Extract sensing time, spatial bounds, day/night flag, and MOD03 production time
    from MOD03 metadata.

    Returns:
        ((beg_iso, end_iso), (west, south, east, north), daynight, mod03_prod_iso)

    For this UMM-G file:
    - beg/end are used for TemporalExtent.
    - daynight is used for DataGranule.DayNightFlag.
    - mod03_prod_iso is returned for reference but not used as the derived
      MODMISR ProductionDateTime.
    """
    txt = _read_core_archive_text(mod03_path) if mod03_path else ""

    # Times from CoreMetadata.0
    bd = pvl_obj_str(txt, "RANGEBEGINNINGDATE")
    bt = pvl_obj_str(txt, "RANGEBEGINNINGTIME")
    ed = pvl_obj_str(txt, "RANGEENDINGDATE")
    et = pvl_obj_str(txt, "RANGEENDINGTIME")

    if bd and bt and ed and et:
        try:
            b = datetime.fromisoformat(f"{bd}T{bt}").replace(tzinfo=timezone.utc)
            e = datetime.fromisoformat(f"{ed}T{et}").replace(tzinfo=timezone.utc)
            beg_iso, end_iso = iso_z(b), iso_z(e)
        except Exception:
            beg_iso = iso_z(fallback_center_dt - timedelta(minutes=2, seconds=30))
            end_iso = iso_z(fallback_center_dt + timedelta(minutes=2, seconds=30))
    else:
        beg_iso = iso_z(fallback_center_dt - timedelta(minutes=2, seconds=30))
        end_iso = iso_z(fallback_center_dt + timedelta(minutes=2, seconds=30))

    # Bounding box from ArchiveMetadata.0
    west = pvl_obj_num(txt, "WESTBOUNDINGCOORDINATE")
    east = pvl_obj_num(txt, "EASTBOUNDINGCOORDINATE")
    south = pvl_obj_num(txt, "SOUTHBOUNDINGCOORDINATE")
    north = pvl_obj_num(txt, "NORTHBOUNDINGCOORDINATE")

    if None in (west, east, south, north):
        west, south, east, north = -180.0, -90.0, 180.0, 90.0

    daynight = pvl_obj_str(txt, "DAYNIGHTFLAG") or "Unspecified"

    mod03_prod_iso = None
    prod = pvl_obj_str(txt, "PRODUCTIONDATETIME")
    if prod:
        try:
            mod03_prod_iso = iso_z(datetime.fromisoformat(prod.replace("Z", "+00:00")))
        except Exception:
            mod03_prod_iso = None

    return (beg_iso, end_iso), (west, south, east, north), daynight, mod03_prod_iso


# ---------- read MODMISR NetCDF metadata ----------
def exact_geo_range(path):
    """
    Read exact spatial bounds from MODMISR NetCDF global attributes.
    """
    import xarray as xr

    with xr.open_dataset(path, decode_times=False) as ds:
        north = float(ds.attrs["NORTHBOUNDINGCOORDINATE"])
        south = float(ds.attrs["SOUTHBOUNDINGCOORDINATE"])
        west = float(ds.attrs["WESTBOUNDINGCOORDINATE"])
        east = float(ds.attrs["EASTBOUNDINGCOORDINATE"])

    return north, south, west, east


def netcdf_date_created(path):
    """
    Return global attribute date_created from the MODMISR NetCDF file, if present.

    This is preferred for UMM-G DataGranule.ProductionDateTime because the JSON
    describes the derived MODMISR product, not the source MOD03 granule.

    Falls back to None if unavailable or unparsable.
    """
    import xarray as xr

    try:
        with xr.open_dataset(path, decode_times=False) as ds:
            val = ds.attrs.get("date_created")

        if not val:
            return None

        return iso_z(datetime.fromisoformat(str(val).replace("Z", "+00:00")))

    except Exception:
        return None


def derived_product_production_time(path):
    """
    ProductionDateTime for the derived MODMISR granule.

    Priority:
    1. NetCDF global date_created
    2. File modification time
    """
    return netcdf_date_created(path) or iso_z(
        datetime.fromtimestamp(os.path.getmtime(path), tz=timezone.utc)
    )


# ---------- UMM-G builder ----------
def build_ummg(
    mmcth_path,
    mod03_root,
    collection_shortname,
    collection_version,
    project=None,
    pge_version=None,
):
    base = os.path.basename(mmcth_path)
    year, doy, hhmm, t0 = parse_mmcth_stamp(base)

    mod03_path = find_mod03(mod03_root, t0)

    (beg, end), (_mod03_west, _mod03_south, _mod03_east, _mod03_north), daynight, _mod03_prod = (
        extract_temporal_spatial_daynight(mod03_path, t0)
    )

    # Use exact bounds from the derived MODMISR NetCDF, not from MOD03.
    north, south, west, east = exact_geo_range(mmcth_path)

    # ProductionDateTime should describe the derived MODMISR NetCDF granule.
    prod_iso = derived_product_production_time(mmcth_path)

    file_size_mb = os.path.getsize(mmcth_path) / (1024.0 * 1024.0)

    ummg = {
        "AccessConstraints": {
            "Description": "WORLD",
            "Value": 4,
        },
        "GranuleUR": base,
        "ProviderDates": provider_dates_from_mtime(mmcth_path),
        "CollectionReference": {
            "ShortName": collection_shortname,
            "Version": str(collection_version),
        },
        "MetadataSpecification": {
            "Name": "UMM-G",
            "URL": "https://cdn.earthdata.nasa.gov/umm/granule/v1.6.5",
            "Version": "1.6.5",
        },
        "DataGranule": {
            "ProductionDateTime": prod_iso,
            "DayNightFlag": daynight,
            "Identifiers": [
                {
                    "Identifier": base,
                    "IdentifierType": "ProducerGranuleId",
                }
            ],
            "ArchiveAndDistributionInformation": [
                {
                    "Name": base,
                    "Size": file_size_mb,
                    "SizeUnit": "MB",
                    "Checksum": {
                        "Algorithm": "MD5",
                        "Value": md5sum(mmcth_path),
                    },
                }
            ],
        },
        "TemporalExtent": {
            "RangeDateTime": {
                "BeginningDateTime": beg,
                "EndingDateTime": end,
            }
        },
        "SpatialExtent": {
            "HorizontalSpatialDomain": {
                "Geometry": {
                    "BoundingRectangles": [
                        {
                            "WestBoundingCoordinate": west,
                            "SouthBoundingCoordinate": south,
                            "EastBoundingCoordinate": east,
                            "NorthBoundingCoordinate": north,
                        }
                    ]
                }
            }
        },
    }

    if project:
        ummg["Projects"] = [{"ShortName": project}]

    if pge_version:
        ummg["PGEVersionClass"] = {"PGEVersion": str(pge_version)}

    return ummg, mod03_path


# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(
        description="Minimal-plus UMM-G generator for MODMISR_L2_CP_V00 granules."
    )

    ap.add_argument(
        "--mmcth",
        required=True,
        help="Path to MODMISR/MMCTH .nc file. Filename must contain _Ayyyyddd.hhmm_.",
    )

    ap.add_argument(
        "--mod03-root",
        default="/data/gdi/satellite/TerraDataArchive/MODIS/MOD03",
        help=(
            "Root of MOD03 tree: "
            "<root>/<YYYY>/<DDD>/MOD03.Ayyyyddd.hhmm*.hdf"
        ),
    )

    ap.add_argument(
        "--collection-shortname",
        default="MODMISR_L2_CP",
        help="CollectionReference.ShortName.",
    )

    ap.add_argument(
        "--collection-version",
        default="V00",
        help="CollectionReference.Version.",
    )

    ap.add_argument(
        "--project",
        help="Projects[0].ShortName, optional.",
    )

    ap.add_argument(
        "--pge-version",
        help="PGEVersionClass.PGEVersion, optional.",
    )

    ap.add_argument(
        "-o",
        "--out",
        default=None,
        help="Output JSON filename. Must end with .cmr.json.",
    )

    ap.add_argument(
        "-od",
        "--outdir",
        default=None,
        help="Output directory. Default: same directory as --mmcth.",
    )

    args = ap.parse_args()

    mmcth_path = Path(args.mmcth).resolve()
    if not mmcth_path.exists():
        raise FileNotFoundError(f"Input NetCDF file not found: {mmcth_path}")

    out_dir = Path(args.outdir).resolve() if args.outdir else mmcth_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.out:
        out_name = args.out
    else:
        if mmcth_path.suffix == ".nc":
            out_name = mmcth_path.name.replace(".nc", ".nc.cmr.json")
        else:
            print(f"Warning: input file does not end with .nc: {mmcth_path.name}")
            out_name = f"{mmcth_path.name}.cmr.json"

    if not out_name.endswith(".cmr.json"):
        raise ValueError(
            f"Output metadata filename must end with '.cmr.json', got: {out_name}"
        )

    out_path = out_dir / out_name

    ummg, mod03_path = build_ummg(
        str(mmcth_path),
        args.mod03_root,
        args.collection_shortname,
        args.collection_version,
        project=args.project,
        pge_version=args.pge_version,
    )

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(ummg, f, indent=2, sort_keys=False)

    print(f"Wrote {out_path}")
    print(f"MOD03 matched: {mod03_path if mod03_path else 'NONE (used fallback temporal/daynight)'}")
    print(f"GranuleUR: {ummg['GranuleUR']}")
    print(f"CollectionReference: {ummg['CollectionReference']['ShortName']} {ummg['CollectionReference']['Version']}")
    print(f"ProductionDateTime: {ummg['DataGranule']['ProductionDateTime']}")
    print(
        "TemporalExtent: "
        f"{ummg['TemporalExtent']['RangeDateTime']['BeginningDateTime']} to "
        f"{ummg['TemporalExtent']['RangeDateTime']['EndingDateTime']}"
    )


if __name__ == "__main__":
    main()