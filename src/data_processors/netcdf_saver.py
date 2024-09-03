
import xarray as xr
import numpy as np
class NetCDFSaver:
    def __init__(self, filename,logger):
        self.filename = filename
        self.logger = logger

    def save(self, bands_BT,misr_cth, mod_geo, mod06,era5_variables_misrswath):
        try:
            # Assuming mod_geo contains 'lat' and 'lon'
            mod_lat = np.array(mod_geo['lat'])
            mod_lon = np.array(mod_geo['lon'])
            
            lat_dim = mod_lat.shape[0]
            lon_dim = mod_lon.shape[1]
            mod_lat_flat = mod_lat[:, 0]  # Assuming lat is the same for each column
            mod_lon_flat = mod_lon[0, :]  # Assuming lon is the same for each row
            
            # Convert the dictionaries to numpy arrays if necessary

            bands_BT_arrays = np.array([v for v in bands_BT.values()])
            bands_BT_keys = list(bands_BT.keys())
           

            mod06 = {k: np.array(v) for k, v in mod06.items()} if isinstance(mod06, dict) else {}
            era5_variables_misrswath_2d = {k: np.array(v) for k, v in era5_variables_misrswath.items() if v.ndim == 2}
            era5_variables_misrswath_3d = {k: np.array(v) for k, v in era5_variables_misrswath.items() if v.ndim == 3}

            # Create a dataset
            ds = xr.Dataset(
                {
                    "misr_cth": (["lat", "lon"], misr_cth),
                    **{k: (["lat", "lon"], v) for k, v in era5_variables_misrswath_2d.items()},
                    **{k: (["level", "lat", "lon"], v) for k, v in era5_variables_misrswath_3d.items()},
                    "bands_BT": (["band", "lat", "lon"], bands_BT_arrays),
                    "mod06_product": (["product", "lat", "lon"], np.array(list(mod06.values()))),
                    "landsea": (["lat", "lon"], mod_geo['landsea']),
                    "vza": (["lat", "lon"], mod_geo['vza'])
                },
                coords={
                    "latitude": (["lat"], mod_lat_flat),
                    "longitude": (["lon"], mod_lon_flat),
                    "band": list(bands_BT.keys()),
                    "product": list(mod06.keys()),
                    "variable_2d": list(era5_variables_misrswath_2d.keys()),
                    "variable_3d": list(era5_variables_misrswath_3d.keys()),
                    "level": np.arange(1, 102)  # Assuming pressure levels are from 1 to 101
                },
            )

            # Add metadata
            ds.attrs["description"] = "Processed MISR, MODIS, and ERA5 data"
            ds.attrs["history"] = "Created " + str(np.datetime64('now', 's'))
            ds.attrs["source"] = "MISR, MODIS, and ERA5 data"

            ds["latitude"].attrs["units"] = "degrees_north"
            ds["longitude"].attrs["units"] = "degrees_east"
            ds["misr_cth"].attrs["units"] = "meters"
            for var in era5_variables_misrswath_2d:
                ds[var].attrs["units"] = "Unknown"
            for var in era5_variables_misrswath_3d:
                ds[var].attrs["units"] = "Unknown"
            ds["bands_BT"].attrs["units"] = "Kelvin"
            ds["mod06_product"].attrs["units"] = "Unknown"
            ds["landsea"].attrs["units"] = "Unknown"
            ds["vza"].attrs["units"] = "degrees"

            # Save to NetCDF file
            comp = dict(zlib=True, complevel=8)
            encoding = {var: comp for var in ds.data_vars}

            # Save to NetCDF file with compression
            ds.to_netcdf(self.filename, encoding=encoding)
            self.logger.info(f"Successfully saved data to {self.filename}")

        except Exception as e:
            self.logger.error(f"Failed to save data to {self.filename}: {e}", exc_info=True)
