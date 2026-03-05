#!/bin/bash

# Set top level variables
ROOT=$MODIS_L2_HOME
PROD=MOD_PR06CD
MDIR=$ROOT/src/cloudtop/$PROD

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Check arguments
if [ $# != 5 ]; then
  echo "Usage: run_mod06ct.sh SAT FIL1KM FILGEO FILCM FILOUT FILCSRB FILGDAS GDAS1 GDAS2 RSST NISE ICEC"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  echo "FILGEO is the Level 1B 1000 meter geolocation HDF file"
  echo "FILCM  is the Level 2 cloud mask HDF file"
  echo "FILOUT is the Level 2 cloud product HDF file"
  exit 1
fi

# Extract arguments
SAT=$1
FIL1KM=`basename $2`
FILGEO=`basename $3`
FILCM=`basename $4`
FILOUT=$5
DIR1KM=`dirname $2`
DIRGEO=`dirname $3`
DIRCM=`dirname $4`

# Get platform header (MOD or MYD)
if [ $SAT == "terra" ] || [ $SAT == "Terra" ] || [ $SAT == "TERRA" ]; then
  HEADER="MOD"
elif [ $SAT == "aqua" ] || [ $SAT == "Aqua" ] || [ $SAT == "AQUA" ]; then
  HEADER="MYD"
else
  echo "Satellite name not recognized: "$SAT
  exit 1
fi

#-------------------------------------------------------------------------------
# SET UP INPUT FILES
#-------------------------------------------------------------------------------

# Create links to static input files
ln -f -s $MDIR/src/${HEADER}*.mcf .
ln -f -s $MDIR/coeff/* .

#-------------------------------------------------------------------------------
# Copy leapsec.dat and utcpole.dat files to run directory
#-------------------------------------------------------------------------------

cp $LOCAL_ANC_DIR/util/leapsec.dat .
cp $LOCAL_ANC_DIR/util/utcpole.dat .

#-------------------------------------------------------------------------------
# SET UP PCF FILE
#-------------------------------------------------------------------------------

# Get date and time string from L1B 1KM meter file
DATE_TIME=`echo $FIL1KM:t | cut -d. -f2,3 | cut -dA -f2`

# Set current date/time
CURRENT=`date -u +%Y%j%H%M%S`

# Set name of template file
TEMPLATE=$MDIR/template/${HEADER}_PR06CD.pcf.template

# Set name of new PCF file
FILPCF=${HEADER}_PR06CD.$DATE_TIME.pcf

# Set names of output files
FILOUT_L2=$FILOUT
FILOUT_QC=${HEADER}06_L2.$DATE_TIME.QC

# Set dummy start and stop times for PCF
START_TIME='2000-01-01T00:00:00'
STOP_TIME=`date -u +%Y-%m-%dT%H:%M:%S`

# Create new PCF file from the template
echo $DIRGDAS
sed \
  -e "s?FIL1KM?${FIL1KM:t}?g" \
  -e "s?FILGEO?${FILGEO:t}?g" \
  -e "s?FILCM?${FILCM:t}?g" \
  -e "s?DIR1KM?${DIR1KM}?g" \
  -e "s?DIRGEO?${DIRGEO}?g" \
  -e "s?DIRCM?${DIRCM}?g" \
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
export PGSMSG=$PGSHOME/message

# Print start message
echo $PROD "processing started at "`date`

# Run the algorithm
$ROOT/bin/cirruscloud.exe
PGE_STATUS=$?

# Clean up
if [ $PGE_STATUS == 0 ]; then
  rm -f $FILOUT_L2.met $FILOUT_QC $FILOUT_QC.met
  rm -f $FILPCF ShmMem LogUser LogReport GetAttr.temp
  rm -f ${HEADER}*.mcf *MCF
  rm -f TRANS* 
  rm -f pc*
  rm -f leapsec.dat utcpole.dat
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
