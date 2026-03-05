#!/bin/bash

#------------------------------------------------------------------------------
# Get the Collection 6 MODIS 8 Day Clear Sky ancillary file for Aqua or Terra 
# closest to an imput MODIS date (yyyyddd)
#
#  Updated March 2011 to use "wget" instead of "ncftpget".  
#  Kathleen Strabala  UW-Madison, SSEC, kathy.strabala@ssec.wisc.edu
#
# Note that the following environment variables must be defined:
# LOCAL_ANC_DIR is the local ancillary data directory
#   e.g. $HOME/ancillary
# REMOTE_ANC_DIR is the remote FTP ancillary data directory (URL format)
#   e.g. ftp://ftp.ssec.wisc.edu/pub/eosdb/ancillary
#
# Uses dateplus.exe command found in /imapp_modisl2/bin directory. The
# source code can be found at: http://www.orlandokuntao.com/dateplus_c.html
#------------------------------------------------------------------------------

# Check number of arguments
if [ $# != 2 ]; then
  echo "Usage: get_anc_clearsky.sh satellite date"
  echo "where"
  echo "  satellite is Aqua or Terra"
  echo "  date is the required date (yyyyddd)"
  exit 1
fi

# Check if should download ancillary data form external site if not found locally
download=$DOWNLOAD_ANC

if [[ $download == true || $download == TRUE || $download == True ]]; then

  # Before we begin, make sure that wget is available
  which wget > /dev/null
  if [ $? != 0 ]; then
      echo " "  > /dev/stderr
      echo "**** ERROR:  Could not find the wget utility. The ancillary" > /dev/stderr
      echo "**** data fetchers will not work without it. Please install" > /dev/stderr
      echo "**** and try again." > /dev/stderr
      echo " " > /dev/stderr
      exit 1
  fi
fi

# Extract file names
SAT=$1
date=$2

if [[ "$SAT" == "Aqua" || "$SAT" == "aqua" ]]; then
  PREFIX="MYD"
elif [[ "$SAT" == "Terra" || "$SAT" == "terra" ]]; then
  PREFIX="MOD"
else
  echo "Satellite name incorrect" $SAT
  exit 1
fi

#--- Dataset specific information ---
# Set NISE filename head and tail
file_head="${PREFIX}CSR_B.A"
file_tail="hdf"

# Set Clear Sky File date range 
# Because the dates on the files represent the first day of the 
#  compositing, you really want files that have date stamps 
#  8 days previous to the given date

# Check local directory
if [ ! -d $LOCAL_ANC_DIR ]; then
  echo "Local ancillary data directory does not exist:" $LOCAL_ANC_DIR
  exit 1
fi

if [[ $date < 2000001 || $date > 2100365 ]]; then
  echo "Invalid date string:" $date
  exit 1
fi

# Compute Gregorian date (yyyymmdd)
greg_date=`dateplus.exe -J $date`

# Search date range for an acceptable file
for delta in -8 -9 -7 -10 -6 -11 -5 -12 -4 -13 -3 -14 -2 -15 -1 -16 0 -17 -18 -19 -20
do
  
  # Compute Gregorian date adjusted by delta (yymmdd)
  tmp_file_date=`dateplus.exe $delta $greg_date | cut -c 1-8`
  year=`echo $tmp_file_date | cut -c1-4`
  month=`echo $tmp_file_date | cut -c5-6`
  dd=`echo $tmp_file_date | cut -c7-8`

  # Now convert back to Julian Date
  file_date=`date +"%Y%j" --date="$year/$month/$dd"`
  
  # Set file name
  file_name=${file_head}${file_date}".006."
  echo "Searching for file with prefix "$file_name > /dev/stderr
    
  # Get year and day of year (yyyy, ddd)
  day=`echo $file_date | cut -c5-8`
  
  # Set local file path and name
  DAY_DIR=${year}_${month}_${dd}_${day}
  local_file=$LOCAL_ANC_DIR/$DAY_DIR

  # Check for file in local directory
  success=`find $LOCAL_ANC_DIR/$DAY_DIR -name ${file_name}\*`
  if [ "$success" != "" ]; then

    # File was found locally
    echo "File was found on local disk" > /dev/stderr
    echo ${local_file}/${file_name}*
    exit 0

  else

    if [[ $download == true || $download == TRUE || $download == True ]]; then

      # File was not found locally, so try to download it
      echo "Trying to download file from "$REMOTE_ANC_DIR/$DAY_DIR/${file_name}\* > /dev/stderr
      #Replaced ncftpget with wget  --  more flexible use for proxies
      wget -q -N -t 30 -nd --no-parent -r -l1 "$REMOTE_ANC_DIR/$DAY_DIR/${file_name}*.hdf"
      success1=`find . -name ${file_name}\*`
      if [ "$success1" != "" ]; then
        if [ ! -d $LOCAL_ANC_DIR/$DAY_DIR ]; then mkdir $LOCAL_ANC_DIR/$DAY_DIR; fi
        mv ${file_name}* ${local_file}/.
        if [ $? == 0 ]; then
          echo "File was downloaded successfully " > /dev/stderr
          echo ${local_file}/${file_name}*
          exit 0
        fi
      fi

    fi

  fi

done

# Ancillary file was not found
echo "ERROR: Ancillary Clear Sky Bias file could not be found" > /dev/stderr
exit 1
