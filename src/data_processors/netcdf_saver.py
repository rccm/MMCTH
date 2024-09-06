import xarray as xr
import numpy as np

class NetCDFSaver:
    def __init__(self, filename, logger):
        self.filename = filename
        self.logger = logger

    def save(self, bands_BT=None, misr_cth=None, mod_geo=None, mod06=None, era5_variables_misrswath=None):
        try:
            # Ensure mod_geo exists
            if mod_geo is None:
                raise ValueError("mod_geo must be provided.")
            
            # Extract latitude and longitude
            mod_lat = np.array(mod_geo['lat'])
            mod_lon = np.array(mod_geo['lon'])
            
            lat_dim = mod_lat.shape[0]
            lon_dim = mod_lon.shape[1]
            mod_lat_flat = mod_lat[:, 0]  # Assuming lat is the same for each column
            mod_lon_flat = mod_lon[0, :]  # Assuming lon is the same for each row

            # Prepare bands_BT arrays
            if bands_BT:
                bands_BT_arrays = np.array([v for v in bands_BT.values()])
                bands_BT_keys = list(bands_BT.keys())
            else:
                bands_BT_arrays = None

            # Prepare mod06 data
            if mod06:
                mod06 = {k: np.array(v) for k, v in mod06.items()} if isinstance(mod06, dict) else {}
            
            # Prepare ERA5 data
            era5_variables_misrswath_2d = {k: np.array(v) for k, v in era5_variables_misrswath.items() if v.ndim == 2} if era5_variables_misrswath else {}
            era5_variables_misrswath_3d = {k: np.array(v) for k, v in era5_variables_misrswath.items() if v.ndim == 3} if era5_variables_misrswath else {}

            # Create dataset with optional variables
            data_vars = {}
            
            if misr_cth:
                data_vars["misr_cth"] = (["lat", "lon"], misr_cth["misrcth"])
                data_vars["misr_cth_qa"] = (["lat", "lon"], misr_cth["misrcth_qa"])

            if bands_BT_arrays is not None:
                data_vars["bands_BT"] = (["band", "lat", "lon"], bands_BT_arrays)

            if mod06:
                data_vars["mod06_product"] = (["product", "lat", "lon"], np.array(list(mod06.values())))

            # Add ERA5 variables
            data_vars.update({k: (["lat", "lon"], v) for k, v in era5_variables_misrswath_2d.items()})
            data_vars.update({k: (["level", "lat", "lon"], v) for k, v in era5_variables_misrswath_3d.items()})
            
            # Add geographic data (always present)
            data_vars["landsea"] = (["lat", "lon"], mod_geo['landsea'])
            data_vars["vza"] = (["lat", "lon"], mod_geo['vza'])

            # Create a dataset
            ds = xr.Dataset(
                data_vars,
                coords={
                    "latitude": (["lat"], mod_lat_flat),
                    "longitude": (["lon"], mod_lon_flat),
                    "band": bands_BT_keys if bands_BT else [],
                    "product": list(mod06.keys()) if mod06 else [],
                    "variable_2d": list(era5_variables_misrswath_2d.keys()) if era5_variables_misrswath_2d else [],
                    "variable_3d": list(era5_variables_misrswath_3d.keys()) if era5_variables_misrswath_3d else [],
                    "level": np.arange(1, 102)  # Assuming pressure levels are from 1 to 101
                }
            )

            # Add metadata
            ds.attrs["description"] = "Processed MISR, MODIS, and ERA5 data"
            ds.attrs["history"] = "Created " + str(np.datetime64('now', 's'))
            ds.attrs["source"] = "MISR, MODIS, and ERA5 data"

            ds["latitude"].attrs["units"] = "degrees_north"
            ds["longitude"].attrs["units"] = "degrees_east"
            if "misr_cth" in ds:
                ds["misr_cth"].attrs["units"] = "meters"
            for var in era5_variables_misrswath_2d:
                ds[var].attrs["units"] = "Unknown"
            for var in era5_variables_misrswath_3d:
                ds[var].attrs["units"] = "Unknown"
            if "bands_BT" in ds:
                ds["bands_BT"].attrs["units"] = "Kelvin"
            if "mod06_product" in ds:
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
