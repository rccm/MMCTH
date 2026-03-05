#!/bin/bash

# This script gets the ancillary data and sets up the input files and environment for running
# the cloudmask algorithm (MOD35).  Additional scripts are called to create the skeleton 
# cloudmask file (run_mod35.sh) and run the LEOCAT algorthim for creating the cloudmask (run_leocat_mod35.sh).
# An executable is called to convert the cloud mask file that is produced to a cloud mask file
# that has the cloudmask data in a more readable format (*.maskbyte1.hdf).  

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "MOD35 cloud mask processing started at "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Set top level variables
ROOT=$MODIS_L2_HOME

# Check arguments
if [ $# != 12 ]; then
  echo "Usage: run_modis_cloudmask.sh SAT FIL1KM FILQKM FILGEO GDAS1 GDAS2 RSST NISE ICEC REAGG"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  echo "FILQKM is the Level 1B  250 meter radiance HDF file (optional: enter MISSING if not used)"
  echo "FILGEO is the Level 1B 1000 meter geolocation HDF file"
  echo "GDAS1 is the first GDAS GRIB file"
  echo "GDAS2 is the second GDAS GRIB file"
  echo "RSST is the Reynolds SST binary file"
  echo "NISE is the snow/ice HDFEOS file"
  echo "ICEC is the sea ice GRIB file"
  echo "OUTDIR is the product output directory"
  echo "REAGG is set to 1 to create re-aggregated products, 0 otherwise"
  echo "COMP is set to 1 to compress output product, 0 otherwise"
  exit 1
fi

# Extract arguments
SAT=$1
FIL1KM=$2
FILQKM=$3
FILGEO=$4
FILGDAS1=$5
FILGDAS2=$6
DIRGDAS=`dirname $5`
FILRSST=$7
FILNISE=$8
FILICEC=$9
OUTDIR=${10}
REAGG=${11}
COMP=${12}

# Get platform header (MOD or MYD)
if [ $SAT == "terra" ] || [ $SAT == "Terra" ] || [ $SAT == "TERRA" ]; then
  HEADER="MOD"
  SAT="terra"
elif [ $SAT == "aqua" ] || [ $SAT == "Aqua" ] || [ $SAT == "AQUA" ]; then
  HEADER="MYD"
  SAT="aqua"
else
  echo "Satellite name not recognized: "$SAT
  exit 1
fi

BASE_FIL1KM=`basename $FIL1KM`
# Get date and time string from L1B 1KM meter file
DATE_TIME=`echo $BASE_FIL1KM:t | cut -d. -f2,3 | cut -dA -f2`
DATE=`echo $DATE_TIME | cut -d. -f1 | cut -dA -f2`
TEMP_TIME=`echo $DATE_TIME | cut -d. -f2`
TIME=`echo $TEMP_TIME | awk '{print $1 + 0}'`

DATE_ORIG=$DATE
if [ ${#DATE} -eq 5 ]; then
  DATE=20$DATE
fi

#-------------------------------------------------------------------------------
# CHECK FOR LEOCAT COEFFICIENT DIRECTORY
#-------------------------------------------------------------------------------

if [ ! -e $MODIS_L2_HOME/src/leocat/coeff ]; then
  echo
  echo "WARNING: No LEOCAT coefficient directory found. MOD35 Cloud Mask"
  echo " processing cannot be run."
  exit 1
fi
 
#-------------------------------------------------------------------------------
# GET ANCILLARY DATA 
#-------------------------------------------------------------------------------
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Getting ancillary data for MODIS cloud mask processing"
echo

# Update leapsec.dat and utcpole.dat files if more than 7 days old
# ----------------------------------------------------------------
get_anc_leapsec_utcpole.bash 

# GDAS1 numerical weather prediction model grid field analysis
# ----------------------------------------------------------------
echo "Date and time from L1B 1KM file: "$DATE $TIME
GDASnearest=`get_anc_gdas_gfs.sh $DATE $TIME`

# Get GDAS or GFS files before and after the given time - RMC 2 Jan 2014
# Updated RMC 16 Sept 2015
# ----------------------------------------------------------------
GDASname=`basename $GDASnearest` 
rm $GDASnearest

GDAStype=`echo $GDASname | cut -c1-3`
if [ "$GDAStype" == "gda" ]; then
  GDAS_before_after=`get_anc_gdas12hrgfs.sh gdas $DATE $TIME`
elif [ "$GDAStype" == "gfs" ]; then
  GDAS_before_after=`get_anc_gdas12hrgfs.sh gfs $DATE $TIME`
else
    echo "ERROR: Global forecast file type not recognized: must be GDAS or GFS"
    exit 1
fi

GDASbefore=`echo $GDAS_before_after | cut -d ' ' -f 1`
GDASafter=`echo $GDAS_before_after | cut -d ' ' -f 2`

# If you can't find a file, set it equal to missing
status=$?
if [ $status != 0 ]; then GDASbefore='MISSING' GDASafter='MISSING'; fi

echo "Nearest GDAS file before is " $GDASbefore
echo "Nearest GDAS file after is " $GDASafter

# GDAS1 numerical weather prediction model grid field analysis
# ----------------------------------------------------------------
echo "Date and time from L1B 1KM file: "$DATE $TIME
GDASnearest=`get_anc_gdas_gfs.sh $DATE $TIME`
status=$?
# If you can't find a file, then try an older gdas file.
if [ "$status" != "0" ]; then
  GDASnearest=`$ROOT/scripts/get_nearest_gdas.sh $DATE $TIME`
  status=$?
  if [ "$status" != "0" ]; then
    GDASnearest="MISSING"
  fi
fi

echo "Nearest GDAS file is " $GDASnearest

# Reynolds Blended SST
# ----------------------------------------------------------------
# Run script which checks local directory first, and then ftp site
RSST=`get_anc_rsst_wday.sh $DATE`
# If you can't find a file, set it equal to missing
status=$?
if [ $status != 0 ]; then RSST='MISSING'; fi
# Convert backslashes so that they will be properly inserted using sed
ROISST=`echo ${RSST} | sed s/"\/"/"\\\/"/g`

# NISE snow and ice files
# ----------------------------------------------------------------
NISE=`get_anc_nise.sh $DATE`
# If you can't find a file, set it equal to missing
status=$?
if [ $status != 0 ]; then NISE="MISSING"; fi

# NCEP ice concentration file
# ----------------------------------------------------------------
ICEC=`get_anc_icec.sh $DATE`
# If you can't find a file, set it equal to missing
status=$?
if [ $status != 0 ]; then ICEC="MISSING"; fi


#-------------------------------------------------------------------------------
# RUN THE ALGORITHMS
#-------------------------------------------------------------------------------

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running MOD_PR35 skeleton file code..."
echo run_mod35.sh $SAT $FIL1KM $FILQKM $FILGEO $GDASnearest $GDASbefore $GDASafter $RSST $NISE $ICEC $REAGG
run_mod35.sh $SAT $FIL1KM $FILQKM $FILGEO $GDASnearest $GDASbefore $GDASafter $RSST $NISE $ICEC $REAGG

# Set name of PCF file
FILPCF=${HEADER}_PR35.$DATE_TIME.pcf
# Set names of L2 output file
FILOUT_L2=${HEADER}35_L2.$DATE_TIME.006.hdf
# output LEOCAT file name looks like MOD_LEO29_L2_2006240_2120.hdf
FILOUT_ALG29_L2=./${HEADER}_LEO29_L2_${DATE_ORIG}_${TEMP_TIME}.hdf

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running LEOCAT for cloud mask..."
run_leocat_mod35.sh $SAT $FILOUT_L2 $FILPCF $FILOUT_ALG29_L2 $REAGG

PGE_STATUS=$?

if [ "$PGE_STATUS" == "0" ]; then
  rm -f $FILPCF
  rm -f $FILOUT_ALG29_L2
#else
#  echo "Cloudmask did not finish processing"
  #exit 1
fi

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Renaming product files to direct broadcast convention..."

if [ $HEADER == "MOD" ]; then
  DB_HEADER="t1"
else
  DB_HEADER="a1"
fi

if [ `echo $DATE_TIME | cut -c1-2` == "XX" ]; then
  DB_DATETIME=`echo $DATE_TIME | cut -c3-12`
else
  DB_DATETIME=$DATE_TIME
fi

cp $FILOUT_L2 ${DB_HEADER}.${DB_DATETIME}.mod35.hdf

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running convert_cloudmask..."
$ROOT/bin/convert_cloudmask.exe ${DB_HEADER}.${DB_DATETIME}.mod35.hdf ${DB_HEADER}.${DB_DATETIME}.mask_byte1.hdf
echo

if [ "$COMP" == "1" ]; then
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Compressing output files..."
  hrepack.sh ${DB_HEADER}.${DB_DATETIME}.mod35.hdf ${DB_HEADER}.${DB_DATETIME}.mask_byte1.hdf
fi

# Move output files to output directory
/bin/cp ${DB_HEADER}.${DB_DATETIME}.mod35.hdf $OUTDIR
/bin/mv ${DB_HEADER}.${DB_DATETIME}.mask_byte1.hdf $OUTDIR

# Print finish message
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Cloud mask processing ended at  "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo
exit 0
