#!/bin/bash
#
# Run the IMAPP MODIS regression retrieval software, created
# by Eva Borabas, CIMSS.
#
# Set top level variables
ROOT=$MODIS_L2_HOME
PROD=MOD_PR07
MDIR=$ROOT/src/profiles

echo
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "MOD07 profiles processing started at "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Check arguments
if [ $# != 8 ]; then
  echo "Usage: run_profiles.sh SAT FIL1KM FILGEO FILCM GDAS"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  echo "FILGEO is the Level 1B 1000 meter geolocation HDF file"
  echo "FILCM  is the Level 2 cloud mask HDF file"
  echo "GDAS is the GDAS1 GRIB file"
  echo "OUTDIR is the output products directory"
  echo "REAGG is set to 1 to create re-aggregated products, 0 otherwise"
  echo "COMP is set to 1 to compress output product, 0 otherwise"
  exit 1
fi

# Extract arguments
SAT=$1
FIL1KM=`basename $2`
FILGEO=`basename $3`
FILCM=`basename $4`
FILGDAS=`basename $5`
DIR1KM=`dirname $2`
DIRGEO=`dirname $3`
DIRCM=`dirname $4`
DIRGDAS=`dirname $5`
OUTDIR=$6
REAGG=$7
COMP=$8

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

#-------------------------------------------------------------------------------
# SET UP INPUT FILES
#-------------------------------------------------------------------------------

# Update leapsec.dat and utcpole.dat files if more than 7 days old
# ----------------------------------------------------------------
get_anc_leapsec_utcpole.bash

# Create links to static input files
ln -f -s $MDIR/src/${HEADER}*.MCF .
ln -f -s $MDIR/coeff/*${SAT}* .
ln -f -s $MDIR/coeff/*${SAT}*.dat.* .
ln -f -s $MDIR/coeff/MODIS_senzen.bin.le .

#-------------------------------------------------------------------------------
# Copy leapsec.dat and utcpole.dat files to run directory
#-------------------------------------------------------------------------------

cp $LOCAL_ANC_DIR/util/leapsec.dat .
cp $LOCAL_ANC_DIR/util/utcpole.dat .

#-------------------------------------------------------------------------------
# GET ANCILLARY DATA 
#-------------------------------------------------------------------------------

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

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Getting ancillary data for MODIS profiles processing"
echo

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

FILGDAS=`basename $GDASnearest`
DIRGDAS=`dirname $GDASnearest`

# NISE snow and ice files
# ----------------------------------------------------------------
NISE=`get_anc_nise.sh $DATE`
# If you can't find a file, set it equal to missing
status=$?
if [ $status != 0 ]; then NISE="MISSING"; fi
FILNISE=`basename $NISE`
DIRNISE=`dirname $NISE`

# NCEP ice concentration file
# ----------------------------------------------------------------
ICEC=`get_anc_icec.sh $DATE`
# If you can't find a file, set it equal to missing
if [ $status != 0 ]; then ICEC="MISSING"; fi
FILICEC=`basename $ICEC`
DIRICEC=`dirname $ICEC`

#-------------------------------------------------------------------------------
# SET UP PCF FILE
#-------------------------------------------------------------------------------

# Set current date/time
CURRENT=`date -u +%Y%j%H%M%S`

# Set name of template file
TEMPLATE=$MDIR/template/${HEADER}_PR07.pcf.template

# Set name of new PCF file
FILPCF=${HEADER}_PR07.$DATE_TIME.pcf

# Set names of output files
FILOUT_L2=${HEADER}07_L2.$DATE_TIME.hdf
FILOUT_QC=${HEADER}07_L2.$DATE_TIME.QC

# Set dummy start and stop times for PCF
START_TIME='2000-01-01T00:00:00'
STOP_TIME=`date -u +%Y-%m-%dT%H:%M:%S`

# Create new PCF file from the template
sed \
  -e "s?FIL1KM?${FIL1KM:t}?g" \
  -e "s?FILGEO?${FILGEO:t}?g" \
  -e "s?FILCM?${FILCM:t}?g" \
  -e "s?FILGDAS?${FILGDAS:t}?g" \
  -e "s?FILNISE?${FILNISE}?g" \
  -e "s?FILICEC?${FILICEC}?g" \
  -e "s?DIR1KM?${DIR1KM}?g" \
  -e "s?DIRGEO?${DIRGEO}?g" \
  -e "s?DIRCM?${DIRCM}?g" \
  -e "s?DIRGDAS?${DIRGDAS}?g" \
  -e "s?DIRNISE?${DIRNISE}?g" \
  -e "s?DIRICEC?${DIRICEC}?g" \
  -e "s?FILOUT_L2?${FILOUT_L2}?g" \
  -e "s?FILOUT_QC?${FILOUT_QC}?g" \
  -e "s?START_TIME?${START_TIME}?g" \
  -e "s?STOP_TIME?${STOP_TIME}?g" \
  $TEMPLATE > $FILPCF

#-------------------------------------------------------------------------------
# RUN THE ALGORITHM
#-------------------------------------------------------------------------------

# Set name of PCF file
export PGS_PC_INFO_FILE=./$FILPCF

# Print start message
echo $PROD "processing started at "`date`

# Run the algorithm
if [ $REAGG == 1 ]; then
  time $ROOT/bin/profiles.exe $REAGG
else
  time $ROOT/bin/profiles.exe
fi

if [ $? != 0 ]; then
  echo "********** Error encountered **********"
  exit 1
else
  echo "********** Successful completion **********"
fi

PGE_STATUS=$?

# Clean up
if [ $PGE_STATUS == 0 ]; then
  rm -f $FILOUT_L2.met $FILOUT_QC $FILOUT_QC.met
  rm -f $FILPCF ShmMem LogUser LogReport GetAttr.temp
  rm -f ${HEADER}*.MCF
  rm -f pc*
  rm -f MODIS_senzen.bin.le
  rm -f MODIS_REGCOEF_FACTORS.${SAT}.*
  rm -f MODIS_REGCOEF_col5.${SAT}
  rm -f ${SAT}_det.dat.v2
  rm -f ${SAT}_bias.dat.v2
  rm -f ${SAT}_bias.dat.v2.test
  rm -f ${SAT}_bias.dat.v2_old
  rm -f ${SAT}_bias.dat.v2.old
  rm -f ${SAT}_bias.dat.v6
  rm -f ${SAT}_zero_bias.dat.v2
  rm -f ${SAT}_bias_zeros.dat
  rm -f leapsec.dat utcpole.dat
fi

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Renaming product files to direct broadcast convention..."
echo

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

/bin/mv $FILOUT_L2 ${DB_HEADER}.${DB_DATETIME}.mod07.hdf

if [ "$COMP" == "1" ]; then
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Compressing output file..."
  hrepack.sh ${DB_HEADER}.${DB_DATETIME}.mod07.hdf
fi

/bin/cp ${DB_HEADER}.${DB_DATETIME}.mod07.hdf $OUTDIR

# Print finish message
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Profiles processing ended at "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo
exit 0
