#!/bin/bash
#
# Run the IMAPP MODIS cloud top properties software
#

echo
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "MOD06 cloud top properties processing started at  "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Set top level variables
ROOT=${MODIS_L2_HOME}
DATEPLUS=$ROOT/scripts/dateplus.exe

# Check arguments
if [ $# != 13 ]; then
  echo "Usage: run_modis_cloudtop.sh SAT FIL1KM FILGEO FILMOD35 FILCSKB GDAS1 GDAS2 RSST NISE ICEC REAGG"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  echo "FILGEO is the Level 1B 1000 meter geolocation HDF file"
  echo "FILMOD35 is the C6 MOD35 cloud mask product file"
  echo "FILCSRB is the MODIS clear sky radiance bias HDF file"
  echo "GDAS1 is the first GDAS GRIB file"
  echo "GDAS2 is the second GDAS GRIB file"
  echo "RSST is the Reynolds SST binary file"
  echo "NISE is the snow/ice HDFEOS file"
  echo "ICEC is the sea ice GRIB file"
  echo "OUTDIR is the output product directory"
  echo "REAGG is set to 1 to create re-aggregated products, 0 otherwise"
  echo "COMP is set to 1 to compress output product, 0 otherwise"
  exit 1
fi

# Extract arguments
SAT=$1
FIL1KM=$2
FILGEO=$3
FILMOD35=$4
FILCSRB=$5
FILGDAS1=$6
FILGDAS2=$7
DIRGDAS=`dirname $6`
FILRSST=$8
FILNISE=$9
FILICEC=${10}
OUTDIR=${11}
REAGG=${12}
COMP=${13}

# Get platform header (MOD or MYD)
if [ $SAT == "terra" ] || [ $SAT == "Terra" ] || [ $SAT == "TERRA" ]; then
  HEADER="MOD"
  FILESAT="Terra"
  SAT="terra"
elif [ $SAT == "aqua" ] || [ $SAT == "Aqua" ] || [ $SAT == "AQUA" ]; then
  HEADER="MYD"
  FILESAT="Aqua"
  SAT="aqua"
else
  echo "Satellite name not recognized: "$SAT
  exit 1
fi

BASE_FIL1KM=`basename $FIL1KM`

# Get date and TIME string from L1B 1KM meter file
DATE_TIME=`echo $BASE_FIL1KM:t | cut -d. -f2,3 | cut -dA -f2`
DATE=`echo $DATE_TIME | cut -d. -f1 | cut -dA -f2`
TEMP_TIME=`echo $DATE_TIME | cut -d. -f2`
TIME=`echo $TEMP_TIME | awk '{print $1 + 0}'`

DATE_ORIG=$DATE
if [ ${#DATE} -eq 5 ]; then
  DATE=20$DATE
fi

#-------------------------------------------------------------------------------
# SET FILE NAMES
#-------------------------------------------------------------------------------

# Output LEOCAT file name looks like MOD_LEO17_L2_2006240_2120.hdf
FILOUT_ALG17_L2=./${HEADER}_LEO17_L2_${DATE}_${TEMP_TIME}.hdf

# Set name of PCF file
FIL_PR06_PCF=${HEADER}_PR06CT.$DATE_TIME.pcf

# Set names of L2 output file
FIL_PR06_OUT=${HEADER}06CT_L2.$DATE_TIME.006.hdf

#-------------------------------------------------------------------------------
# CHECK FOR LEOCAT COEFFICIENT DIRECTORY
#-------------------------------------------------------------------------------

if [ ! -e $MODIS_L2_HOME/src/leocat/coeff ]; then
  echo
  echo "WARNING: No LEOCAT coefficient directory found. MOD06 Cloud Top"
  echo " processing cannot be run."
  exit 1
fi

#-------------------------------------------------------------------------------
# GET ANCILLARY DATA 
#-------------------------------------------------------------------------------

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Getting ancillary data for MODIS cloud top processing"
echo

# Update leapsec.dat and utcpole.dat files if more than 7 days old
# ----------------------------------------------------------------
get_anc_leapsec_utcpole.bash

# Clear sky radiance files
# ----------------------------------------------------------------
FILCSRB=`get_anc_clearsky.sh $SAT $DATE`
status=$?
if [ $status != 0 ]; then
    FILCSRB="MISSING"
    echo "Warning: clear sky radiance file not found"
    exit 1
fi

echo "Clear sky radiance file is " $FILCSRB

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

# Get GDAS or GFS files before and after the given time - RMC 2 Jan 2014
# ----------------------------------------------------------------
GDASname=`basename $GDASnearest`
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
if [ $status != 0 ]; then ICEC="MISSING"; fi


#-------------------------------------------------------------------------------
# RUN THE ALGORITHMS
#-------------------------------------------------------------------------------
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running MOD_PR06 file creation code..."
$ROOT/scripts/run_mod06cr.sh $SAT $FIL1KM $FIL_PR06_OUT

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running MOD_PR06 5-km code..."
$ROOT/scripts/run_mod06ct.sh $SAT $FIL1KM $FILGEO $FILMOD35 $FIL_PR06_OUT $FILCSRB $GDASnearest $GDASbefore $GDASafter $RSST $NISE $ICEC $REAGG

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running MOD_PR06 cirrus detection code..."
$ROOT/scripts/run_mod06cd.sh $SAT $FIL1KM $FILGEO $FILMOD35 $FIL_PR06_OUT

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running LEOCAT for cloud top properties..."
$ROOT/scripts/run_leocat_mod06.sh $SAT $FIL_PR06_OUT $FIL_PR06_PCF $FILOUT_ALG17_L2

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

/bin/mv $FIL_PR06_OUT ${DB_HEADER}.${DB_DATETIME}.mod06ct.hdf

if [ "$COMP" == "1" ]; then
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Compressing output file..."
  hrepack.sh ${DB_HEADER}.${DB_DATETIME}.mod06ct.hdf
fi

/bin/cp ${DB_HEADER}.${DB_DATETIME}.mod06ct.hdf $OUTDIR

# Print finish message
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Cloud top properties processing ended at  "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo
exit 0
