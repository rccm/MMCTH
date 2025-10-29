#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimal-plus UMM-G v1.6.5 generator for MMCTH using LOCAL MOD03.

Required (UMM-G): GranuleUR, ProviderDates, CollectionReference, MetadataSpecification
Minimal-plus extras: DataGranule.ProductionDateTime, Archive+MD5, TemporalExtent, BoundingRectangle, DayNightFlag

Assumptions:
- MMCTH filename contains _Ayyyyddd.hhmm_ and corresponds 1:1 to a MOD03 granule.
- MOD03 root is /data/gdi/satellite/TerraDataArchive/MODIS/MOD03/<YYYY>/<DDD>/MOD03.Ayyyyddd.hhmm*.hdf
"""

import argparse, os, re, json, hashlib, subprocess
from datetime import datetime, timezone, timedelta
from pathlib import Path

# ---------- tiny utils ----------
def iso_z(dtobj: datetime) -> str:
    if dtobj.tzinfo is None:
        dtobj = dtobj.replace(tzinfo=timezone.utc)
    return dtobj.astimezone(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00","Z")

def md5sum(path, buf=1<<20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()

def provider_dates_from_mtime(path):
    t = datetime.fromtimestamp(os.path.getmtime(path), tz=timezone.utc)
    ts = iso_z(t)
    return [{"Date": ts, "Type": "Create"}, {"Date": ts, "Type": "Update"}]

def parse_mmcth_stamp(basename):
    m = re.search(r"_A(\d{4})(\d{3})\.(\d{2})(\d{2})_", basename)
    if not m:
        raise ValueError(f"Cannot find _Ayyyyddd.hhmm_ in: {basename}")
    year, doy, HH, MM = map(int, m.groups())
    t0 = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy-1, hours=HH, minutes=MM)
    return year, doy, f"{HH:02d}{MM:02d}", t0

# ---------- very fast MOD03 lookup (direct folder, scandir) ----------
def _find_by_prefix(dirpath, prefix):
    try:
        with os.scandir(dirpath) as it:
            for e in it:
                if not e.is_file(): continue
                name = e.name
                if name.startswith(prefix) and name.endswith(".hdf"):
                    return os.path.join(dirpath, name)
    except FileNotFoundError:
        return None
    return None

def find_mod03(mod03_root, t0):
    # exact folder + minute, then small fallbacks that can cross day/year
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

# ---------- read ECS PVL bits from MOD03 (CoreMetadata.0 / ArchiveMetadata.0) ----------
def _read_core_archive_text(mod03_path):
    try:
        from pyhdf.SD import SD, SDC  # faster path if available
        sd = SD(mod03_path, SDC.READ)
        attrs = sd.attributes()
        txt = f"{attrs.get('CoreMetadata.0','')}\n{attrs.get('ArchiveMetadata.0','')}"
        sd.end()
        if txt.strip():
            return txt
    except Exception:
        pass
    # fallback: ncdump -h
    try:
        return subprocess.check_output(["ncdump", "-h", mod03_path], text=True, errors="ignore")
    except Exception:
        return ""

_PVL_FLAGS = re.IGNORECASE | re.DOTALL

def pvl_obj_num(txt: str, key: str):
    """
    Find numeric VALUE inside an OBJECT block:
      OBJECT = <KEY> ... VALUE = <num> ... END_OBJECT = <KEY>
    Falls back to 'KEY = <num>' if present.
    """
    if not txt:
        return None
    pat = rf"OBJECT\s*=\s*{re.escape(key)}\s.*?VALUE\s*=\s*([\-0-9\.Ee\+]+)\s.*?END_OBJECT\s*=\s*{re.escape(key)}"
    m = re.search(pat, txt, _PVL_FLAGS)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            return None
    # fallback: KEY = <num>
    m = re.search(rf"{re.escape(key)}\s*=\s*([\-0-9\.Ee\+]+)", txt, re.IGNORECASE)
    return float(m.group(1)) if m else None

def pvl_obj_str(txt: str, key: str):
    """
    Find quoted string VALUE inside an OBJECT block, else fallback to KEY="str".
    """
    if not txt:
        return None
    pat = rf"OBJECT\s*=\s*{re.escape(key)}\s.*?VALUE\s*=\s*\"([^\"]+)\"\s.*?END_OBJECT\s*=\s*{re.escape(key)}"
    m = re.search(pat, txt, _PVL_FLAGS)
    if m:
        return m.group(1)
    # fallback: KEY = "str"
    m = re.search(rf"{re.escape(key)}\s*=\s*\"([^\"]+)\"", txt, re.IGNORECASE)
    return m.group(1) if m else None

# UPDATE your metadata extraction to use the PVL-aware helpers:
def extract_temporal_spatial_daynight(mod03_path, fallback_center_dt):
    txt = _read_core_archive_text(mod03_path) if mod03_path else ""

    # Times from CoreMetadata.0 (PVL)
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

    # Bounding box from ArchiveMetadata.0 (PVL)
    west  = pvl_obj_num(txt, "WESTBOUNDINGCOORDINATE")
    east  = pvl_obj_num(txt, "EASTBOUNDINGCOORDINATE")
    south = pvl_obj_num(txt, "SOUTHBOUNDINGCOORDINATE")
    north = pvl_obj_num(txt, "NORTHBOUNDINGCOORDINATE")

    if None in (west, east, south, north):
        # safe fallback if a rare granule is missing coords
        west, south, east, north = -180.0, -90.0, 180.0, 90.0

    # Day/Night + MOD03 production time (optional)
    daynight = pvl_obj_str(txt, "DAYNIGHTFLAG") or "Unspecified"
    prod = pvl_obj_str(txt, "PRODUCTIONDATETIME")
    prod_iso = None
    if prod:
        try:
            prod_iso = iso_z(datetime.fromisoformat(prod.replace("Z","+00:00")))
        except Exception:
            prod_iso = None

    return (beg_iso, end_iso), (west, south, east, north), daynight, prod_iso

def exact_geo_range(path):
    import xarray as xr
    ds = xr.open_dataset(path)
    north = ds.attrs["NORTHBOUNDINGCOORDINATE"]
    south = ds.attrs["SOUTHBOUNDINGCOORDINATE"]
    west  = ds.attrs["WESTBOUNDINGCOORDINATE"]
    east  = ds.attrs["EASTBOUNDINGCOORDINATE"]
    return north, south, west, east


# ---------- UMM-G builder ----------
def build_ummg(mmcth_path, mod03_root, collection_shortname, collection_version,
               project=None, pge_version=None):
    base = os.path.basename(mmcth_path)
    year, doy, hhmm, t0 = parse_mmcth_stamp(base)

    mod03_path = find_mod03(mod03_root, t0)
    (beg, end), (west, south, east, north), daynight, mod03_prod = extract_temporal_spatial_daynight(mod03_path, t0)
    north, south, west, east = exact_geo_range(mmcth_path)
    # print(north, south, west, east)
    # ProductionDateTime: prefer MOD03 prod time; else this file's mtime
    prod_iso = mod03_prod or iso_z(datetime.fromtimestamp(os.path.getmtime(mmcth_path), tz=timezone.utc))

    ummg = {
        "AccessConstraints": {
            "Description": "WORLD",
            "Value": 4
        },
        "GranuleUR": base,
        "ProviderDates": provider_dates_from_mtime(mmcth_path),
        "CollectionReference": {
            "ShortName": collection_shortname,
            "Version": str(collection_version)
        },
        "MetadataSpecification": {
            "Name": "UMM-G",
            "URL": "https://cdn.earthdata.nasa.gov/umm/granule/v1.6.5",
            "Version": "1.6.5"
        },
        "DataGranule": {
            # "ProductionDateTime": prod_iso,
            "ProductionDateTime": provider_dates_from_mtime(mmcth_path),
            "DayNightFlag": daynight,
            "Identifiers": [
                {"Identifier": base, "IdentifierType": "ProducerGranuleId"}
            ],
            "ArchiveAndDistributionInformation": [
                {
                    "Name": base,
                    "Size": os.path.getsize(mmcth_path) / (1024.0*1024.0),
                    "SizeUnit": "MB",
                    "Checksum": {"Algorithm": "MD5", "Value": md5sum(mmcth_path)}
                }
            ]
        },
        "TemporalExtent": {
            "RangeDateTime": {"BeginningDateTime": beg, "EndingDateTime": end}
        },
        "SpatialExtent": {
            "HorizontalSpatialDomain": {
                "Geometry": {
                    "BoundingRectangles": [{
                        "WestBoundingCoordinate": west,
                        "SouthBoundingCoordinate": south,
                        "EastBoundingCoordinate": east,
                        "NorthBoundingCoordinate": north
                    }]
                }
            }
        }
    }
    if project:
        ummg["Projects"] = [{"ShortName": project}]
    if pge_version:
        ummg["PGEVersionClass"] = {"PGEVersion": str(pge_version)}
    return ummg, mod03_path

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Minimal-plus UMM-G generator (fast local MOD03 lookup; no DB).")
    ap.add_argument("--mmcth", required=True, help="Path to MMCTH .nc (must contain _Ayyyyddd.hhmm_)")
    ap.add_argument("--mod03-root", default="/data/gdi/satellite/TerraDataArchive/MODIS/MOD03",
                    help="Root of MOD03 tree: <root>/<YYYY>/<DDD>/MOD03.Ayyyyddd.hhmm*.hdf")
    ap.add_argument("--collection-shortname", default="MMCTH", help="CollectionReference.ShortName (e.g., MMCTH_L2_CP)")
    ap.add_argument("--collection-version",  default="V01", help="CollectionReference.Version (e.g., V01)")
    ap.add_argument("--project", help="Projects[0].ShortName (optional)")
    ap.add_argument("--pge-version", help="PGEVersionClass.PGEVersion (optional)")
    ap.add_argument("-o",  "--out",    default=None, help="Output file name (default: <input>.nc.cmr.json)")
    ap.add_argument("-od", "--outdir", default=None, help="Output directory (default: same dir as --mmcth)")

    args = ap.parse_args()

    # Resolve paths
    mmcth_path = Path(args.mmcth).resolve()
    out_dir = Path(args.outdir).resolve() if args.outdir else mmcth_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Choose output filename
    if args.out:
        out_name = args.out
    else:
        # Keep your original behavior: file.nc -> file.nc.cmr.json; otherwise fall back to <name>.cmr.json
        if mmcth_path.suffix == ".nc":
            out_name = mmcth_path.name.replace(".nc", ".nc.cmr.json")
        else:
            print(f"Warning: input file does not end with .nc : {mmcth_path.name}")
            out_name = f"{mmcth_path.name}.cmr.json"

    out_path = out_dir / out_name

    # Build UMM-G
    ummg, mod03_path = build_ummg(
        str(mmcth_path),
        args.mod03_root,
        args.collection_shortname,
        args.collection_version,
        project=args.project,
        pge_version=args.pge_version
    )

    # Write JSON
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(ummg, f, indent=2, sort_keys=False)

    print(f"Wrote {out_path}")
    print(f"MOD03 matched: {mod03_path if mod03_path else 'NONE (used fallback temporal/bbox)'}")
if __name__ == "__main__":
    main()