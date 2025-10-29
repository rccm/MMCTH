#!/usr/bin/env python3
import sys, glob, os
from netCDF4 import Dataset
from datetime import date, timedelta

vari = 'single'

if vari == 'profile':
    req_vars = {'q', 'z', 'w', 'u', 'v', 't'}
elif vari == 'single':
    req_vars = {'d2m', 't2m', 'sst', 'skt', 'sp', 'msl'}
else:
    print('Invalid Variable Name ')
    exit

if len(sys.argv) < 3:
    print("Usage: python check_vars.py <directory> <years_or_ranges> [output_file]")
    print("Examples:")
    print("  python check_vars.py /data 2007 2014 2022 missing.txt")
    print("  python check_vars.py /data 2002-2004 2007 missing.txt")
    sys.exit(1)

dir_path = sys.argv[1]
output_file = sys.argv[-1] if sys.argv[-1].endswith('.txt') else "missing_files.txt"
year_args = sys.argv[2:-1] if sys.argv[-1].endswith('.txt') else sys.argv[2:]

# Parse year ranges
years = set()
for arg in year_args:
    if '-' in arg:
        start, end = map(int, arg.split('-'))
        years.update(range(start, end + 1))
    else:
        years.add(int(arg))


missing = []

print(f"Checking years: {sorted(years)}")

for year in sorted(years):
    print(f"Checking year {year}...")
    
    # Get all files for this year
    # existing_files = glob.glob(f"{dir_path}/era5_profile_{year}_*.nc")
    existing_files = glob.glob(f"{dir_path}/era5_{vari}_{year}_*.nc")
    existing_dates = set()
    
    # Extract dates from existing files
    for file in existing_files:
        try:
            # Extract date from filename: era5_profile_YYYY_MM_DD.nc
            basename = os.path.basename(file)
            parts = basename.split('_')
            if len(parts) >= 5:
                file_year = int(parts[2])
                file_month = int(parts[3])
                file_day = int(parts[4].split('.')[0])
                existing_dates.add((file_year, file_month, file_day))
        except (ValueError, IndexError):
            # Skip files with unexpected naming pattern
            continue
    
    # Check for all days in the year
    start_date = date(year, 1, 1)
    end_date = date(year, 12, 31)
    current_date = start_date
    
    while current_date <= end_date:
        # expected_file = f"{dir_path}/era5_profile_{current_date.year}_{current_date.month:02d}_{current_date.day:02d}.nc"
        expected_file = f"{dir_path}/era5_{vari}_{current_date.year}_{current_date.month:02d}_{current_date.day:02d}.nc"
        
        # Check if file exists
        if (current_date.year, current_date.month, current_date.day) not in existing_dates:
            missing.append(f"{expected_file} (File missing)")
            print(f"Missing: {expected_file}")
        else:
            # File exists, check variables
            try:
                with Dataset(expected_file, 'r') as nc:
                    if not req_vars.issubset(nc.variables.keys()):
                        missing_vars = req_vars - set(nc.variables.keys())
                        missing.append(f"{expected_file} (Missing variables: {', '.join(sorted(missing_vars))})")
                        print(f"{expected_file} (Missing variables: {', '.join(sorted(missing_vars))})")
            except Exception as e:
                missing.append(f"{expected_file} (Error: {str(e)})")
                print(f"{expected_file} (Error: {str(e)})")
        
        current_date += timedelta(days=1)

# Write results
with open(output_file, 'w') as f:
    for item in missing:
        f.write(item + '\n')

print(f"\nFound {len(missing)} files with issues")
print(f"Output saved to: {output_file}")