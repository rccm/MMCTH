#!/bin/bash

#------------------------------------------------------------------------------
# Get the correct surface albedo date based upon the input data time
#
#------------------------------------------------------------------------------

# Check number of arguments
if [ $# != 1 ]; then
  echo "Usage: get_surface_albedo_date.sh date"
  echo "where"
  echo "  file_date is the required date (yyyyddd)"
  exit 1
fi


# Set date ranges (file every 32 days)

date_range=(353 321 289 257 225 \
                  193 161 129 97 65 33 1)

# Get arguments
file_date=$1
if [[ $file_date < 2000001 || $file_date > 2100365 ]]; then
  echo "Invalid date string:" $file_date
  exit 1
fi

tmp_doy=`echo $file_date | cut -c5-7`

#Remove leading zeros from julian date if they exist
doy=`expr $tmp_doy + 0`

for delta in "${date_range[@]}"
do
  if [ $doy -ge  $delta ]; then
     if [ $delta -lt 10 ]; then
       JULDAY=00${delta}
       echo $JULDAY
       exit 0 
     fi
     if [ $delta -lt 100 ]; then
       JULDAY=0${delta}
       echo $JULDAY
       exit 0 
     else 
       JULDAY=${delta}
       echo $JULDAY
       exit 0 
     fi
  fi
done

exit 1 
