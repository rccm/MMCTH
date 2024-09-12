# merge_netcdf.py
import xarray as xr
import argparse
import sys

def merge_netcdf(file1, file2, output_file):
    try:
        # Open the two NetCDF datasets
        ds1 = xr.open_dataset(file1)
        ds2 = xr.open_dataset(file2)

 
        merged_ds = xr.merge([ds1, ds2])

        merged_ds.to_netcdf(output_file)

        print(f"Files successfully merged into {output_file}")

    except Exception as e:
        print(f"Error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Create an argument parser for command-line arguments
    parser = argparse.ArgumentParser(description="Merge two NetCDF files.")
    parser.add_argument("file1", help="Path to the first NetCDF file")
    parser.add_argument("file2", help="Path to the second NetCDF file")
    parser.add_argument("output_file", help="Path for the output merged NetCDF file")

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the merge function with provided arguments
    merge_netcdf(args.file1, args.file2, args.output_file)
