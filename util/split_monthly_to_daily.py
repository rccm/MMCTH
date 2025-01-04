import xarray as xr
import pandas as pd
import os
import argparse
import traceback

def split_monthly_to_daily(input_file, output_dir, compress=True):
    """
    Splits a monthly NetCDF file into daily NetCDF files.

    Parameters:
    - input_file (str): Path to the input monthly NetCDF file.
    - output_dir (str): Directory where daily NetCDF files will be saved.
    - compress (bool): Whether to apply compression to output files.
    """
    print(f"Processing file: {input_file}")

    # Open the dataset with chunking based on 'valid_time'
    ds = xr.open_dataset(input_file)  # Adjust chunks based on data resolution

    # Convert 'valid_time' to datetime objects
    ds['valid_time'] = pd.to_datetime(ds['valid_time'].values, unit='s', origin='unix')
    print(f"valid_time dtype: {ds['valid_time'].dtype}")

    # Ensure 'valid_time' is recognized as a coordinate
    if 'valid_time' not in ds.coords:
        ds = ds.set_coords('valid_time')

    # Create a new coordinate with date strings (e.g., '2021_01_01')
    dates_str = pd.to_datetime(ds['valid_time'].values).strftime('%Y_%m_%d')
    ds = ds.assign_coords(date_str=('valid_time', dates_str))

    # Group by the new 'date_str' coordinate
    grouped = ds.groupby('date_str')

    # Extract base name for output files (e.g., 'eradaily_2021_02')
    base_name = os.path.basename(input_file).replace('_v.nc', '')

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define encoding for compression if needed
    if compress:
        encoding = {var: {'zlib': True, 'complevel': 5} for var in ds.data_vars}
    else:
        encoding = None

    # Function to save each group
    def save_group(date, group_ds):
        
        output_filename = f"{base_name}_{date}.nc"
        output_path = os.path.join(output_dir, output_filename)
        group_ds['valid_time'] = (
            (group_ds['valid_time'].values.astype('int64') // 10**9)
        ).astype('int64')
        # Update attributes
        group_ds['valid_time'].attrs['units'] = "seconds since 1970-01-01 00:00:00"
        group_ds['valid_time'].attrs['long_name'] = "time"
        group_ds['valid_time'].attrs['standard_name'] = "time"
        group_ds['valid_time'].attrs['calendar'] = "proleptic_gregorian"

        if compress:
            group_ds.to_netcdf(output_path, encoding=encoding)
        else:
            group_ds.to_netcdf(output_path)
        print(f"Created {output_path}")

    # Iterate over each group and save to a file
    for date, group_ds in grouped:
        save_group(date, group_ds)

    # Close the dataset
    ds.close()

def main():
    parser = argparse.ArgumentParser(description="Split monthly NetCDF files into daily files.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing monthly NetCDF files.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save daily NetCDF files.")
    parser.add_argument('--pattern', type=str, default="eradaily_*_v.nc", help="Glob pattern to match monthly files.")
    parser.add_argument('--compress', action='store_true', help="Enable compression for output files.")

    args = parser.parse_args()

    # Find all matching files based on the pattern
    monthly_files = [
        os.path.join(args.input_dir, f) 
        for f in os.listdir(args.input_dir) 
        if f.startswith("eradaily_") and f.endswith("_v.nc")
    ]

    if not monthly_files:
        print(f"No files matched the pattern 'eradaily_*_v.nc' in directory {args.input_dir}. Exiting.")
        return

    # Process each file sequentially
    for monthly_file in monthly_files:
        try:
            split_monthly_to_daily(monthly_file, args.output_dir, compress=args.compress)
        except Exception as e:
            print(f"Error processing {monthly_file}: {e}")
            traceback.print_exc()

    print("All monthly files have been processed.")

if __name__ == "__main__":
    main()
