#!/bin/bash

# Set top level variables
ROOT=$MODIS_L2_HOME
PROD=MOD_PR06CR
MDIR=$ROOT/src/cloudtop/$PROD

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Check arguments
if [ $# != 3 ]; then
  echo "Usage: run_mod06cr.sh SAT FIL1KM"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  echo "FILOUT is the Level 2 cloud product HDF file"
  exit 1
fi

# Extract arguments
SAT=$1
FIL1KM=`basename $2`
DIR1KM=`dirname $2`
FILOUT_L2=$3

# Get platform header (MOD or MYD)
if [ $SAT == "terra" ]; then
  HEADER="MOD"
elif [ $SAT == "aqua" ]; then
  HEADER="MYD"
else
  echo "Satellite name not recognized: "$SAT
  exit 1
fi

#-------------------------------------------------------------------------------
# SET UP INPUT FILES
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SET UP PCF FILE
#-------------------------------------------------------------------------------

# Get date and time string from L1B 1KM meter file
DATE_TIME=`echo $FIL1KM:t | cut -d. -f2,3 | cut -dA -f2`

# Set current date/time
CURRENT=`date -u +%Y%j%H%M%S`

# Set name of template file
TEMPLATE=$MDIR/template/${HEADER}_PR06CR.pcf.template

# Set name of new PCF file
FILPCF=${HEADER}_PR06CR.$DATE_TIME.pcf
echo $TEMPLATE $FILPCF $FIL1KM $DIR1KM $2

# Set names of output files
FILOUT_QC=${HEADER}06_L2.$DATE_TIME.QC

# Set dummy start and stop times for PCF
START_TIME='2000-01-01T00:00:00'
STOP_TIME=`date -u +%Y-%m-%dT%H:%M:%S`

# Create new PCF file from the template
sed \
  -e "s?FIL1KM?${FIL1KM}?g" \
  -e "s?DIR1KM?${DIR1KM}?g" \
  -e "s?FILOUT_L2?${FILOUT_L2}?g" \
  -e "s?START_TIME?${START_TIME}?g" \
  -e "s?STOP_TIME?${STOP_TIME}?g" \
  $TEMPLATE > $FILPCF

#-------------------------------------------------------------------------------
# RUN THE ALGORITHM
#-------------------------------------------------------------------------------

# Set name of PCF file
export PGS_PC_INFO_FILE=./$FILPCF
export PGSMSG=$PGSHOME/message

# Print start message
echo $PROD "processing started at "`date`

# Run the algorithm
$ROOT/bin/createcloudfile.exe
PGE_STATUS=$?

# Clean up
if [ $PGE_STATUS == 0 ]; then
  rm -f $FILOUT_L2.met $FILOUT_QC $FILOUT_QC.met
  rm -f $FILPCF ShmMem GetAttr.temp
fi

# Print failure or success message
if [ $PGE_STATUS != 0 ]; then
  echo "********** Error encountered **********"
  exit 1
else
  echo "********** Successful completion **********"
fi

# Print finish message
echo $PROD "processing ended at  "`date`
exit 0
