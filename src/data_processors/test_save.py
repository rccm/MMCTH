import xarray as xr
import numpy as np
import time
from dask.diagnostics import ProgressBar

def test_netcdf_write():
    # Define dimensions based on your data
    idx_size = 3976111  # Size of your main dimension
    level_size = 37     # Size of the 'level' dimension

    # Create coordinate arrays
    idx = np.arange(idx_size)
    level = np.arange(level_size)

    # Create synthetic data arrays
    data_vars = {}
    for var_name in ['z', 'u', 'v']:
        data = np.random.rand(idx_size, level_size).astype('float32')
        data_array = xr.DataArray(
            data,
            dims=('idx', 'level'),
            coords={'idx': idx, 'level': level},
            name=var_name
        )
        data_vars[var_name] = data_array

    # Create a dataset
    dataset = xr.Dataset(data_vars)

    # Optionally, add attributes (or leave them empty)
    dataset.attrs = {}
    for var in dataset.data_vars:
        dataset[var].attrs = {}

    # Optimize NetCDF writing
    comp = dict(zlib=False)  # Start with no compression for testing
    encoding = {var: comp for var in dataset.data_vars}

    # Define the output file
    outputfile = 'test_output.nc'
    print(f"Writing dataset to {outputfile}...")

    # Measure the time taken to write the file
    start_time = time.time()

    # Use ProgressBar to monitor the write process
    with ProgressBar():
        dataset.to_netcdf(
            path=outputfile,
            mode="w",
            engine="netcdf4",
            encoding=encoding,
            compute=True
        )

    end_time = time.time()
    print(f"Finished writing in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    test_netcdf_write()
