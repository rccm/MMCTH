#!/bin/bash

output_file="/data/gdi/f/gzhao1/Database/MOD021KM_completelist.txt"

source_dir='/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM'

# Clear the output file
> "$output_file"

# Function to list files within a specified range
list_files() {
  local year=$1
  local start_day=$2
  local end_day=$3
  for day in $(seq -w $start_day $end_day); do
    if [ -d "$source_dir/$year/$day" ]; then
      find "$source_dir/$year/$day" -type f  -name "M*.hdf" >> "$output_file"
    fi
  done
}

# List files for 2000 from March 1 (day 061) to December 31 (day 366)
list_files 2000 061 366

# List files for 2001 to 2021 for all days in the year
for year in {2001..2021}; do
  list_files "$year" 001 366
done

# List files for 2022 up to April 1 (day 091)
list_files 2022 001 091
