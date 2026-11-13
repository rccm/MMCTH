# test_netcdf_saver.py
import numpy as np
import xarray as xr
import logging
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from src.data_processors.xarray_saver import XarraySaver
import numpy as np
import xarray as xr
from pathlib import Path
import logging

# ---------- sample data ----------
def sample_data():
    rng = np.random.default_rng(42)

    # 2D scene size for MM/MODIS/MISR
    ny, nx = 10, 10

    # ERA5 grid (coarser)
    eyy, exx = 20, 20   # any sizes are fine; just consistent

    data = {
        "mm_variables": {
            "cth":        rng.integers(0, 20000, size=(ny, nx), dtype=np.int32),
            "ctp":        rng.uniform(100, 1100, size=(ny, nx)).astype(np.float32),
            "emissivity": rng.uniform(0, 1, size=(ny, nx)).astype(np.float32),
            "opt":        rng.uniform(0, 100, size=(ny, nx)).astype(np.float32),
            "cflag":      rng.integers(0, 4, size=(ny, nx), dtype=np.int32),
            "llflag":      rng.integers(0, 4, size=(ny, nx), dtype=np.int32),
        },
        "mm_geo": {
            "lat": np.linspace(-40, 40, ny).reshape(ny, 1).repeat(nx, axis=1).astype(np.float32),
            "lon": np.linspace(-180, 180, nx).reshape(1, nx).repeat(ny, axis=0).astype(np.float32),
            "vza": rng.uniform(0, 90, size=(ny, nx)).astype(np.float32),
            "sza": rng.uniform(0, 180, size=(ny, nx)).astype(np.float32),
            "time2d": rng.uniform(7.0e8, 7.0e8 + 300, size=(ny, nx)).astype(np.float64),
        },
        # NOTE: The class expects MISR key 'misrcth'
        "misr_variables": {
            "misrcth": rng.integers(0, 20000, size=(ny, nx), dtype=np.int32)
        },
        "modis_variables": {
            "cth":        rng.integers(0, 20000, size=(ny, nx), dtype=np.int32),
            # "ctp":        rng.uniform(100, 1100, size=(ny, nx)).astype(np.float32),
            "emissivity": rng.uniform(0, 1, size=(ny, nx)).astype(np.float32),
            "opt":        rng.uniform(0, 100, size=(ny, nx)).astype(np.float32),
            "cphase":     rng.integers(0, 4, size=(ny, nx), dtype=np.int32),
        },
        # ERA5: mix 3D (z, y, x) with z=37 and 2D (y, x)
        # You can use either mapped keys (e.g., 'temperature_org') or final names (e.g., 'ERA5Temperature').
        "era5_variables": {
            "temperature_org":     rng.uniform(150, 350, size=(37, eyy, exx)).astype(np.float32),  # 3D
            "surface_pressure":    rng.uniform(500, 1100, size=(eyy, exx)).astype(np.float32),     # 2D
            "geopotential_org":    rng.uniform(-5e4, 5e4, size=(37, eyy, exx)).astype(np.float32), # 3D
            "ERA5Humidity":        rng.uniform(0, 1, size=(eyy, exx)).astype(np.float32),          # 2D
        },
        "input_files": [
            "/in/MOD06.hdf",
            "/in/MISR_L2.nc",
            "/in/ERA5_T.nc",
        ],
        "global_attrs": {
            "scene_utc": "2023-05-14T10:24:00Z",
            "pipeline": "mm-v1",
        }
    }
    return ny, nx, data

# ---------- quick test ----------
def test_save_mm(NetCDFSaverClass):
    ny, nx, data = sample_data()
    filename = Path("mm_test.nc")

    logger = logging.getLogger("NetCDFSaverTest")
    logger.setLevel(logging.INFO)

    saver = NetCDFSaverClass(str(filename), complevel=5, logger=logger)
    saver.save_mm(
        mm_variables=data["mm_variables"],
        mm_geo=data["mm_geo"],
        misr_variables=data["misr_variables"],
        modis_variables=data["modis_variables"],
        era5_variables=data["era5_variables"],
        input_files=data["input_files"],
        global_attrs=data["global_attrs"],
        # var_attrs={"MMFlag": {"long_name": "MM cloud processing flag"}}
    )

    # Read back and sanity-check
    ds = xr.open_dataset(filename, decode_coords=False)

    # existence
    assert "MM_CloudTopHeight" in ds, "MM_CloudTopHeight not saved"
    assert "Latitude" in ds, "Latitude not saved"
    assert "Longitude" in ds, "Longitude not saved"
    assert "SZA" in ds, "SZA not saved"

    # shapes for MM plane (must match lat/lon)
    assert ds["MM_CloudTopHeight"].shape == (ny, nx), "Incorrect MM_CloudTopHeight shape"
    assert ds["Latitude"].shape == (ny, nx), "Incorrect Latitude shape"
    assert ds["Longitude"].shape == (ny, nx), "Incorrect Longitude shape"
    assert ds["SZA"].shape == (ny, nx), "Incorrect SZA shape"
    assert ds["SZA"].attrs["long_name"] == "solar viewing zenith angle", "Incorrect SZA long_name"
    assert ds["SZA"].attrs["units"] == "degrees", "Incorrect SZA units"
    assert ds["SZA"].attrs["coordinates"] == "Latitude Longitude pixel_time", "Incorrect SZA coordinates"

    # ERA5 3D: check z length and first/last values
    assert "ERA5_Temperature" in ds, "ERA5_Temperature not saved"
    zcoord = ds.coords.get("z", None)
    assert zcoord is not None, "z coord missing for ERA5 3D variables"
    assert zcoord.size == 37, "ERA5 z must have 37 pressure levels"
    assert int(zcoord.values[0]) == 1000 and int(zcoord.values[-1]) == 1, "ERA5 z levels unexpected"

    ds.close()
    print("✅ basic round-trip checks passed:", filename)

# ---------- run locally ----------
if __name__ == "__main__":
    # from your_module import NetCDFSaver  # if defined elsewhere
    # For inline testing, just pass the class object you defined earlier:
    test_save_mm(XarraySaver)
