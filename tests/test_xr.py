import xarray as xr

fn = "/data/gdi/e/ERA5///single/era5_single_2015_09_23.nc"

for eng in ["h5netcdf", "netcdf4"]:
    try:
        ds = xr.open_dataset(fn, engine=eng)
        ds["sst"].isel(time=0).load()
        print("OK", eng)
    except Exception as e:
        print("FAIL", eng, repr(e))