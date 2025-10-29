from mpi4py import MPI
import xarray as xr
import numpy as np
import os
import netCDF4 as nc

def calculate_scale_offset(data, int_range=(-32766, 32767)):
    x_min = np.nanmin(data)
    x_max = np.nanmax(data)
    scale_factor = (x_max - x_min) / (int_range[1] - int_range[0])
    add_offset = x_min - int_range[0] * scale_factor
    return scale_factor, add_offset

def pack_data(data, scale_factor, add_offset, fill_value=-32767):
    packed = np.round((data - add_offset) / scale_factor).astype(np.int16)
    packed[np.isnan(data)] = fill_value
    return packed

def write_uv_fast(output_file, ds_existing, u_data, v_data):
    with nc.Dataset(output_file, 'w', format='NETCDF4') as dst:
        # Create dimensions from the existing dataset
        for dim_name, dim_len in ds_existing.dims.items():
            dst.createDimension(dim_name, dim_len)

        # Copy coordinate variables: longitude, latitude, level, time
        for var_name in ['longitude', 'latitude', 'level', 'time']:
             # pub
            src_var = ds_existing[var_name]
            dst_var = dst.createVariable(var_name, 'f4', src_var.dims)
            dst_var[:] = src_var.values
            dst_var.setncatts(src_var.attrs)

        # Copy other existing variables (t, q, z)
        for var_name in ['t', 'q', 'z']:
            src_var = ds_existing[var_name]
            dst_var = dst.createVariable(
                var_name, src_var.dtype, src_var.dims,
                zlib=True, complevel=4, shuffle=True
            )
            dst_var[:] = src_var.values
            dst_var.setncatts(src_var.attrs)

        # Create and add u and v variables
        for var_name, data in zip(['u', 'v'], [u_data, v_data]):
            dst_var = dst.createVariable(
                var_name, 'f4', ('time', 'level', 'latitude', 'longitude'),
                zlib=True, complevel=4, shuffle=True, 
            )
            print(data.shape)
            dst_var[:] = data
            dst_var.units = 'm s**-1'
            dst_var.long_name = f'{var_name.upper()} component of wind'

def extract_daily_uv(monthly_file, daily_files_dir, output_dir, write_flag_fast=True):
    try:
        ds_monthly = xr.open_dataset(monthly_file, engine='cfgrib')
    except Exception as e:
        print(f"Error opening {monthly_file}: {e}")
        return

    # Identify unique days in the monthly data
    unique_days = np.unique(ds_monthly['time'].dt.floor('D'))
    for day in unique_days:
        day_str = np.datetime_as_string(day, unit='D')
        # if day_str.split('-')[2][0:2] != '29':
        #     continue
        print(f"Processing day: {day_str} from file: {monthly_file}")
        ds_daily = ds_monthly.sel(time=slice(str(day_str), str(day_str)))

        u_daily = ds_daily['u']
        v_daily = ds_daily['v']

        file_day_str = day_str.replace('-', '_')
        daily_file_name = f"era5_profile_{file_day_str}.nc"
        daily_file_path = os.path.join(daily_files_dir, daily_file_name)
        output_daily_file = os.path.join(output_dir, daily_file_name)

        if os.path.exists(daily_file_path):
            print(f"Found daily file: {daily_file_path}")
            ds_daily_existing = xr.open_dataset(daily_file_path)
            # Rename dimension if necessary
            if 'valid_time' in ds_daily_existing.dims:
                ds_daily_existing = ds_daily_existing.rename({'valid_time': 'time'})
                print("Renamed 'valid_time' dimension to 'time'")
            if 'pressure_level' in ds_daily_existing.dims:
                ds_daily_existing = ds_daily_existing.rename({'pressure_level': 'level'})
                print("Renamed 'pressure_level' dimension to 'level'")
            u_daily_expanded = u_daily.rename({'isobaricInhPa': 'level'})
            v_daily_expanded = v_daily.rename({'isobaricInhPa': 'level'})

            if write_flag_fast:
                try: 
                    write_uv_fast(output_daily_file, ds_daily_existing,
                              u_daily_expanded.values, v_daily_expanded.values)
                except Exception as e:
                    print(f"Error opening {monthly_file} {daily_files_dir}: {e}")
                    return

            else:
                ds_daily_existing['u'] = u_daily_expanded
                ds_daily_existing['v'] = v_daily_expanded
                encoding = {
                    'u': {'zlib': True, 'complevel': 4, 'dtype': 'float32'},
                    'v': {'zlib': True, 'complevel': 4, 'dtype': 'float32'}
                }
                ds_daily_existing.to_netcdf(output_daily_file, mode='w', encoding=encoding)
            print(f"Updated daily file: {output_daily_file}")
        else:
            print(f"Daily file not found: {daily_file_path}")
        # Process all days in the file (remove break if you want to loop through all days)
        # break

def generate_monthly_files(start_year, start_month, end_year, end_month, base_dir='/data/gdi/f/gzhao1/UV/daily'):
    """
    Generates a list of monthly file paths assuming file naming convention:
    'Era5_<year>_<month>_v.grib'
    """
    monthly_files = []
    year = start_year
    month = start_month
    while (year < end_year) or (year == end_year and month <= end_month):
        file_name = f"Era5_{year}_{month}_v.grib"
        file_path = os.path.join(base_dir, file_name)
        monthly_files.append(file_path)
        # Increment month and adjust year if necessary
        month += 1
        if month > 12:
            month = 1
            year += 1
    return monthly_files

if __name__ == "__main__":
    # Initialize MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Generate the list of monthly files between 2000/03 and 2022/02

    monthly_files = generate_monthly_files(2016, 1, 2021, 12)
    files_for_this_rank = [f for i, f in enumerate(monthly_files) if i % size == rank]


    daily_files_dir = '/data/gdi/e/ERA5/multi/'
    output_dir = '/data/gdi/f/gzhao1/UV/test'
    print(f"MPI Rank {rank}/{size} processing {len(files_for_this_rank)} files.")
    for monthly_file in files_for_this_rank:
        if os.path.exists(monthly_file):
            print(f"Rank {rank} processing file: {monthly_file}")
            print(monthly_file)
            extract_daily_uv(monthly_file, daily_files_dir, output_dir, write_flag_fast=True)
        else:
            print(f"Rank {rank} file not found: {monthly_file}")

    #  syn or not?
    # comm.Barrier()
