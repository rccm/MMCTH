#!/bin/bash

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Set top level variables
ROOT=$MODIS_L2_HOME
PROD=leocat/src/MOD_PRAlg29St
MDIR=$ROOT/src/$PROD

# Check arguments
if [ $# != 5 ]; then
  echo "Usage: run_alg29stats.sh SAT FILALG29 FILC5 FILC5PCF"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FILALG29 is the LEOCAT Algorithm 29 Cloud Mask Output hdf file"
  echo "FILC5 is the skeleton MOD35 skeleton file from the Collection 5 run"
  echo "FILC5PCF is the Process Control File (PCF) from the Collection 5 run"
  echo "REAGG is set to 1 to create re-aggregated products, 0 otherwise"
  exit 1
fi

# Extract arguments
SAT=$1
FILOUT_LEOCAT=`basename $2`
DIROUT_LEOCAT=`dirname $2`
FILC5=`basename $3`
DIRC5=`dirname $3`
FILC5PCF=$4
REAGG=$5

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
# CREATE LINKS TO STATIC INPUT FILES
#-------------------------------------------------------------------------------

ln -f -s $MDIR/${HEADER}*.mcf .
ln -f -s $PGSHOME/database/common/TD/leapsec.dat .
ln -f -s $PGSHOME/database/common/CSC/utcpole.dat .

#-------------------------------------------------------------------------------
# CREATE PCF FILE
#-------------------------------------------------------------------------------

# Get date time string from Collection 5 output file
DATE_TIME=`echo $FILC5:t | cut -d. -f2,3`

# Set current date/time
CURRENT=`date -u +%Y%j%H%M%S`

# Get name of PCF file
FILPCF=Alg29Stats.$DATE_TIME.pcf

# Set dummy start and stop times for PCF
START_TIME='2000-01-01T00:00:00'
STOP_TIME=`date -u +%Y-%m-%dT%H:%M:%S`

FILOUT_L2=${HEADER}35_L2.$DATE_TIME.hdf
FILOUT_QC=${HEADER}35_L2.$DATE_TIME.QC

# Set platform name
PLATFORM_NAME=`uname -a`

# Create new PCF file from the template
# Replace the first occurence of $FILOUT_L2 with FILC5|DIRC5
sed \
  -e "s?FILOUT_LEOCAT?${FILOUT_LEOCAT}?g" \
  -e "s?DIROUT_LEOCAT?${DIROUT_LEOCAT}?g" \
  -e "s?${FILOUT_L2}|?${FILC5}|${DIRC5}?" \
  -e "s?FILOUT_QC?${FILOUT_QC}?g" \
  -e "s?START_TIME?${START_TIME}?g" \
  -e "s?STOP_TIME?${STOP_TIME}?g" \
  -e "s?PLATFORM_NAME?${PLATFORM_NAME}?g" \
  $FILC5PCF > $FILPCF

#-------------------------------------------------------------------------------
# RUN THE ALGORITHM
#-------------------------------------------------------------------------------

# Set toolkit environment variables
export PGS_PC_INFO_FILE=./$FILPCF
export PGSMSG=$PGSHOME/message

# Print start message
echo
echo "alg29stats processing started at "`date`

if [ -e $DIROUT_LEOCAT/$FILOUT_LEOCAT ]; then
  echo $DIROUT_LEOCAT/$FILOUT_LEOCAT
else
  echo $DIROUT_LEOCAT/$FILOUT_LEOCAT "Does not exist"
fi

# Run the algorithm
echo "Running " MOD_PRAlg29St.exe
if [ $REAGG == 1 ]; then
  $ROOT/bin/leocat_cloudmaskstats.exe $REAGG
else
  $ROOT/bin/leocat_cloudmaskstats.exe
fi
PGE_STATUS=$?

# Clean up
if [ $PGE_STATUS == 0 ]; then
  rm -f $DIRC5/$FILC5.met $FILOUT_QC $FILOUT_QC.met
  rm -f $FILC5PCF $FILPCF ShmMem LogUser LogReport GetAttr.temp
  rm -f ${HEADER}*.mcf
fi

# Print failure or success message
if [ $PGE_STATUS != 0 ]; then
  echo "********** Error encountered **********"
  echo
  exit 1
else
  echo "********** Successful completion **********"
  echo
fi

# Print finish message
echo "alg29stats processing ended at  "`date`
exit 0
