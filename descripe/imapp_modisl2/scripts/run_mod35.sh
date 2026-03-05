#!/bin/bash

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Set top level variables
ROOT=$MODIS_L2_HOME
PROD=MOD_PR35
MDIR=$ROOT/src/cloudmask

# Check arguments
if [ $# != 11 ]; then
  echo "Usage: run_mod35.sh SAT FIL1KM FILQKM FILGEO FILGDAS GDAS1 GDAS2 RSST NISE ICEC REAGG"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  echo "FILQKM is the Level 1B  250 meter radiance HDF file (optional: enter MISSING if not used)"
  echo "FILGEO is the Level 1B 1000 meter geolocation HDF file"
  echo "FILGDAS is the nearest GDAS GRIB file"
  echo "GDAS1 is the first GDAS GRIB file used by LEOCAT"
  echo "GDAS2 is the second GDAS GRIB file used by LEOCAT"
  echo "RSST is the Reynolds SST binary file"
  echo "NISE is the snow/ice HDFEOS file"
  echo "ICEC is the sea ice GRIB file"
  echo "REAGG is set to 1 to create re-aggregated products, 0 otherwise"
  exit 1
fi

# Extract arguments
SAT=$1
FIL1KM=`basename $2`
FILQKM=`basename $3`
FILGEO=`basename $4`
FILGDAS0=`basename $5`
FILGDAS1=`basename $6`
FILGDAS2=`basename $7`
FILRSST=`basename $8`
FILNISE=`basename $9`
FILICEC=`basename ${10}`
DIR1KM=`dirname  $2`
DIRQKM=`dirname $3`
DIRGEO=`dirname  $4`
DIRGDAS0=`dirname  $5`
DIRGDAS1=`dirname $6`
DIRGDAS2=`dirname $7`
DIRRSST=`dirname  $8`
DIRNISE=`dirname  $9`
DIRICEC=`dirname ${10}`
REAGG=${11}

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
 
# Get date and time string from L1B 1KM meter file
DATE_TIME=`echo $FIL1KM:t | cut -d. -f2,3 | cut -dA -f2`

# Set current date/time
CURRENT=`date -u +%Y%j%H%M%S`

# Set name of template file
TEMPLATE=$MDIR/template/${HEADER}_PR35.pcf.template

# Set name of new PCF file
FILPCF=${HEADER}_PR35.$DATE_TIME.pcf

# Set names of output files
FILOUT_L2=${HEADER}35_L2.$DATE_TIME.006.hdf
FILOUT_QC=${HEADER}35_L2.$DATE_TIME.QC
FILOUT_CS=${HEADER}CSR_G.$DATE_TIME.hdf

# Set dummy start and stop times for PCF
START_TIME='2000-01-01T00:00:00'
STOP_TIME=`date -u +%Y-%m-%dT%H:%M:%S`

# Set platform name 
PLATFORM_NAME=`uname -a`

# Create new PCF file from the template
sed \
  -e "s?FIL1KM?${FIL1KM}?g" \
  -e "s?FILQKM?${FILQKM}?g" \
  -e "s?FILGEO?${FILGEO}?g" \
  -e "s?FILGDAS0?${FILGDAS0:t}?g" \
  -e "s?FILGDAS1?${FILGDAS1:t}?g" \
  -e "s?FILGDAS2?${FILGDAS2:t}?g" \
  -e "s?FILRSST?${FILRSST}?g" \
  -e "s?FILNISE?${FILNISE}?g" \
  -e "s?FILICEC?${FILICEC}?g" \
  -e "s?DIR1KM?${DIR1KM}?g" \
  -e "s?DIRQKM?${DIRQKM}?g" \
  -e "s?DIRGEO?${DIRGEO}?g" \
  -e "s?DIRGDAS0?${DIRGDAS0:t}?g" \
  -e "s?DIRGDAS1?${DIRGDAS1:t}?g" \
  -e "s?DIRGDAS2?${DIRGDAS2:t}?g" \
  -e "s?DIRRSST?${DIRRSST}?g" \
  -e "s?DIRNISE?${DIRNISE}?g" \
  -e "s?DIRICEC?${DIRICEC}?g" \
  -e "s?FILOUT_L2?${FILOUT_L2}?g" \
  -e "s?FILOUT_QC?${FILOUT_QC}?g" \
  -e "s?FILOUT_CS?${FILOUT_CS}?g" \
  -e "s?START_TIME?${START_TIME}?g" \
  -e "s?STOP_TIME?${STOP_TIME}?g" \
  -e "s?PLATFORM_NAME?${PLATFORM_NAME}?g" \
  $TEMPLATE > $FILPCF

#-------------------------------------------------------------------------------
# Create links to static input files
#-------------------------------------------------------------------------------

ln -f -s $MDIR/src/${HEADER}*.MCF .
ln -f -s $MDIR/coeff/*${SAT}* .
ln -f -s $MDIR/coeff/*img.v1 .

#-------------------------------------------------------------------------------
# Copy leapsec.dat and utcpole.dat files to run directory
#-------------------------------------------------------------------------------

cp $LOCAL_ANC_DIR/util/leapsec.dat .
cp $LOCAL_ANC_DIR/util/utcpole.dat .

#-------------------------------------------------------------------------------
# RUN THE ALGORITHM
#-------------------------------------------------------------------------------

# Set toolkit environment variables
export PGS_PC_INFO_FILE=./$FILPCF
export PGSMSG=$PGSHOME/message

# Run the algorithm
if [ $REAGG == 1 ]; then
  if [ $FILQKM == "MISSING" ]; then
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "MOD_PRAGG re-aggregation code not run since no 250m L1B file available"
   
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo $PROD "processing started at "`date` 
    $ROOT/bin/cloudmask.exe
  else 
    echo 
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Running MOD_PRAGG re-aggregation code..."
    export PGS_PC_INFO_FILE=$FILPCF
    $ROOT/bin/reaggregate.exe $REAGG
    
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo $PROD "processing started at "`date`
    $ROOT/bin/cloudmask.exe 1 
  fi
else
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo $PROD "processing started at "`date`
  $ROOT/bin/cloudmask.exe
fi

PGE_STATUS=$?

# Clean up
if [ $PGE_STATUS == 0 ]; then
  rm -f $FILOUT_L2.met $FILOUT_QC $FILOUT_QC.met
  rm -f ShmMem LogUser LogReport GetAttr.temp
  rm -f thresholds.dat.$SAT.*
  rm -f *img.v1
  rm -f leapsec.dat utcpole.dat
  rm -f pc*
  rm -f oisst*
fi

# Print failure or success message
if [ $PGE_STATUS != 0 ]; then
  echo "********** Error encountered creating cloud mask skeleton file **********"
  exit 1
else
  echo "********** Successful creation of cloud mask skeleton file **********"
fi

# Print finish message
echo $PROD "processing ended at  "`date`
exit 0
