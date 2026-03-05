#!/bin/bash

# Get the NCEP Numerical Weather Prediction GRIB File that matches the 
#  input date and time.  Real time data sets will use the global model
#  forecast fields (gfs).  Other data sets will be processed using the global
#  analysis files (gdas1).  The analysis fields should be used whenever
#  possible, but these files can be as much as 7-9 hours behind real-time.
#
# This script searches for two envirnmental variables: 
#
#  $LOCAL_ANC_DIR - The ancillary directory on your local machine.  This
#                        This is where the script will look for a matching
#                        gdas or gfs file, and also where it will download
#                        it to if fetching from a remote ftp site. The 
#                        default directory is the current directory.The 
#
#  $REMOTE_ANC_DIR - The remote ancillary data location to search for a
#                         matching file.  The default site is the 
#                         UW SSEC ancillary data archive at:
#                         ftp://ftp.ssec.wisc.edu/pub/eosdb/ancillary
#
# Note: the script will search for files in a directory structure like this:
#        Main_directory/Date_directory (YYYY_MM_DD_DDD)/ for gdas 
#        Main_directory/Date_directory (YYYY_MM_DD_DDD)/forecast for gfs 
#        files.  It will also create directories on your local machine
#        if the files are downloaded.
#
#  Updated July 2011 to use "wget" instead of "ncftpget".  
#  Translated from csh to bash Jan 2014 by Rebecca Cintineo.
#
# Uses dateplus utility created by Bob Orlando.   For more information,
#    please see his web site: http://www.orlandokuntao.com/dateplus_c.html
#
# 
#  Kathleen Strabala  UW-Madison, SSEC, kathy.strabala@ssec.wisc.edu
#
# --------------------------------------------------------------------------
# Check number of arguments
if [ $# != 3 ]; then
  echo "Usage: get_anc_gdas_gfs12hr.csh type date time"
  echo "where"
  echo "  type is gfs or gdas"
  echo "  date is the required date (yyyyddd)"
  echo "  date is the required time (hhmm)"
  exit 1
fi

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

# If the local ancillary directory variable is not set, then default to 
#  current directory

if [ ! $LOCAL_ANC_DIR ]; then
   export LOCAL_ANC_DIR=${PWD}
fi   

# If the remote ancillary directory variable is not set, then default to 
#  the 
if [ ! $REMOTE_ANC_DIR ]; then
  export REMOTE_ANC_DIR=ftp://ftp.ssec.wisc.edu/pub/eosdb/ancillary
fi   

# Get arguments (note: leading zero is removed from time argument)
type=$1
date=$2
time=`expr $3 + 0`

# check type
if ! [[  $type == "gdas" || $type == "gfs" ]]; then
  echo "Invalid type: $type"
  exit 1
fi   

#Fix problem with time = 0000 
if [[ $time -gt -1 && $time -lt 1 ]]; then
    time=1
fi   

if [[ $date -lt 2000001 || $date -gt 2100365 ]]; then
  echo "Invalid date string:" $date
  exit 1
fi   

# Check local directory - make it if it doesn't exist
mkdir -p $LOCAL_ANC_DIR
if [ ! -d $LOCAL_ANC_DIR ]; then
  echo "Local ancillary data directory does not exist:" $LOCAL_ANC_DIR
  exit 1
fi   

gfs_delta=0


#Now see if any of the files in the listing matches our data times
# Set the day and time of the closest GDAS1 file
localcount=0
found=( )

for relative in before after
do
  if [ $relative == "before" ]; then
    if [ $time -ge 1 ] && [ $time -le  559 ]; then
      gdas_time="00"
      gdas_delta=0
      gfs_time="12"
      gfs_delta=-1
    elif [ $time -ge  600 ] && [ $time -le 1159 ]; then
      gdas_time="06"
      gdas_delta=0
      gfs_time="18"
      gfs_delta=-1
    elif [ $time -ge 1200 ] && [ $time -le 1759 ]; then
      gdas_time="12"
      gdas_delta=0
      gfs_time="00"
      gfs_delta=0
    elif [ $time -ge 1800 ] && [ $time -le 2359 ]; then
      gdas_time="18"
      gdas_delta=0
      gfs_time="06"
      gfs_delta=0
    fi   
  
  elif [ $relative == "after" ]; then
    if [ $time -ge 1 ] && [ $time -le  559 ]; then
      gdas_time="06"
      gdas_delta=0
      gfs_time="18"
      gfs_delta=-1
    elif [ $time -ge  600 ] && [ $time -le 1159 ]; then
      gdas_time="12"
      gdas_delta=0
      gfs_time="00"
      gfs_delta=0
    elif [ $time -ge 1200 ] && [ $time -le 1759 ]; then
      gdas_time="18"
      gdas_delta=0
      gfs_time="06"
      gfs_delta=0
    elif [ $time -ge 1800 ] && [ $time -le 2359 ]; then
      gdas_time="00"
      gdas_delta=1
      gfs_time="12"
      gfs_delta=0
    fi   
  fi   
 
  jday1=`echo $date | cut -c5-7`
  jday=`expr $jday1 - 1`
  data_year=`echo $date | cut -c1-4`
  greg_date=`date --date="$data_year-01-01 + $jday days" "+%y%m%d"`
  greg_month=`date --date="$data_year-01-01 + $jday days" "+%m"`
  greg_day=`date --date="$data_year-01-01 + $jday days" "+%d"`

  # Compute Gregorian date adjusted by delta (yymmdd)
  file_date=`date --date="$data_year-$greg_month-$greg_day $gdas_delta days" "+%y%m%d"`
  year=`echo "20"$file_date | cut -c1-4`
  month=`echo $file_date | cut -c3-4`
  dd=`echo $file_date | cut -c5-6`
  day=`date --date="$year-$greg_month-$greg_day $gdas_delta days" "+%j"`

  # Find the file we want
  if [ $type == "gdas" ]; then
    subdir=""
    file_name="gdas1.PGrbF00."${file_date}.${gdas_time}"z"
  else
    subdir="forecast/"
    gfs_date=`date --date="$data_year-$greg_month-$greg_day $gfs_delta days" "+%y%m%d"`
    year=`echo "20"$gfs_date | cut -c1-4`
    month=`echo $gfs_date | cut -c3-4`
    dd=`echo $gfs_date | cut -c5-6`
    day=`date --date="$year-$month-$dd" "+%j"`
    time_step=12
    file_name="gfs.t"${gfs_time}"."${gfs_date}".pgrbf"${time_step}
  fi   

  #set file_name="gdas1.pgrb00.1p0deg."${file_date}"_"${gdas_time}"_000.grib2"
  echo "Searching for file: "$file_name > /dev/stderr
  
  DAY_DIR=${year}_${month}_${dd}_${day}
  local_file=$LOCAL_ANC_DIR/$DAY_DIR/${subdir}$file_name

  # Check for file in local directory
  if [ -e $local_file ]; then
 
    echo -n "$local_file " 
    localcount=`expr $localcount + 1`
    continue  

  else

    wget --spider -q $REMOTE_ANC_DIR/$DAY_DIR/${subdir}${file_name} > /dev/stderr
    if [ $? == 0 ]; then
      found=( $found "$DAY_DIR/${subdir}${file_name}" )
      #found=("$DAY_DIR/${subdir}${file_name}" ${found[@]})
      continue
    fi   
  fi   
done

total=`expr ${#found[@]} + $localcount`

if [ $localcount == 2 ]; then
  echo "Files were found on local disk" > /dev/stderr
  echo $local_file
  exit 0
fi 

if [ $total == 2 ]; then
  echo "downloading: "$found > /dev/stderr
  for url in ${found[@]}
  do
    local_file=$LOCAL_ANC_DIR/$url
    DAY_DIR=`dirname $url`
    file_name=`basename $url`
    echo "Trying to download file from "$REMOTE_ANC_DIR/$url > /dev/stderr
    wget -q -N -t 30 -O ${file_name} $REMOTE_ANC_DIR/$url > /dev/stderr
    if [ $? == 0 ]; then
      if [ ! -d $LOCAL_ANC_DIR/$DAY_DIR ]; then mkdir -p $LOCAL_ANC_DIR/$DAY_DIR ; fi
      mv $file_name $local_file
      if [ $? == 0 ]; then
        echo -n "$local_file "    # printing file locations
        localcount=`expr $localcount + 1`
      else
        break
      fi   
    fi   
  done
fi   

if [ $localcount == 2 ]; then
  exit 0
fi

echo ""

# If we get here, no matching grib files were found
echo "ERROR: Could not find any gfs or gdas model grib files matching data time" > /dev/stderr

exit 1
