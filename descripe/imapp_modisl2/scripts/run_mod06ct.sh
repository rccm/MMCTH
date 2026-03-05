#!/bin/bash

# Set top level variables
ROOT=$MODIS_L2_HOME
PROD=MOD_PR06CT
MDIR=$ROOT/src/cloudtop/$PROD

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Check arguments
if [ $# != 13 ]; then
  echo "Usage: run_mod06ct.sh SAT FIL1KM FILGEO FILCM FILOUT FILCSRB FILGDAS GDAS1 GDAS2 RSST NISE ICEC REAGG"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  echo "FILGEO is the Level 1B 1000 meter geolocation HDF file"
  echo "FILCM  is the Level 2 cloud mask HDF file"
  echo "FILOUT is the Level 2 cloud product HDF file"
  echo "FILCSRB is the MODIS clear sky radiance bias HDF file"
  echo "FILGDAS is the nearest GDAS GRIB file used by the 5-km code"
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
FILGEO=`basename $3`
FILCM=`basename $4`
FILOUT=$5
FILCSRB=`basename $6`
FILGDAS0=`basename $7`
FILGDAS1=`basename $8`
FILGDAS2=`basename $9`
FILRSST=`basename ${10}`
FILNISE=`basename ${11}`
FILICEC=`basename ${12}`
DIR1KM=`dirname $2`
DIRGEO=`dirname $3`
DIRCM=`dirname $4`
DIRCSRB=`dirname $6`
DIRGDAS=`dirname $7`
DIRRSST=`dirname ${10}`
DIRNISE=`dirname ${11}`
DIRICEC=`dirname ${12}`
REAGG=${13}

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
ln -f -s $ROOT/src/cloudtop/SHR_MCF/${HEADER}*.MCF .
ln -f -s $MDIR/coeff/lit/modisdet.*.v2 .
ln -f -s $MDIR/coeff/lit/global_emiss* .
ln -f -s $MDIR/coeff/goge1_2_img.v1 .

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
TEMPLATE=$MDIR/template/${HEADER}_PR06CT.pcf.template

# Set name of new PCF file
FILPCF=${HEADER}_PR06CT.$DATE_TIME.pcf

# Set names of output files
FILOUT_L2=$FILOUT
FILOUT_QC=${HEADER}06_L2.$DATE_TIME.QC

# Set dummy start and stop times for PCF
START_TIME='2000-01-01T00:00:00'
STOP_TIME=`date -u +%Y-%m-%dT%H:%M:%S`

# Create new PCF file from the template
sed \
  -e "s?FIL1KM?${FIL1KM:t}?g" \
  -e "s?FILGEO?${FILGEO:t}?g" \
  -e "s?FILCM?${FILCM:t}?g" \
  -e "s?FILCSRB?${FILCSRB:t}?g" \
  -e "s?FILGDAS0?${FILGDAS0:t}?g" \
  -e "s?FILGDAS1?${FILGDAS1:t}?g" \
  -e "s?FILGDAS2?${FILGDAS2:t}?g" \
  -e "s?FILRSST?${FILRSST:t}?g" \
  -e "s?FILNISE?${FILNISE:t}?g" \
  -e "s?FILICEC?${FILICEC:t}?g" \
  -e "s?DIR1KM?${DIR1KM}?g" \
  -e "s?DIRGEO?${DIRGEO}?g" \
  -e "s?DIRCM?${DIRCM}?g" \
  -e "s?DIRCSRB?${DIRCSRB}?g" \
  -e "s?DIRGDAS?${DIRGDAS}?g" \
  -e "s?DIRRSST?${DIRRSST}?g" \
  -e "s?DIRNISE?${DIRNISE}?g" \
  -e "s?DIRICEC?${DIRICEC}?g" \
  -e "s?FILOUT_L2?${FILOUT_L2}?g" \
  -e "s?FILOUT_QC?${FILOUT_QC}?g" \
  -e "s?START_TIME?${START_TIME}?g" \
  -e "s?STOP_TIME?${STOP_TIME}?g" \
  -e "s?ACT_LIB_PATH?$MDIR/coeff/lit?g" \
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
$ROOT/bin/cloudtop.exe $REAGG
PGE_STATUS=$?

# Clean up
if [ $PGE_STATUS == 0 ]; then
  rm -f ShmMem LogUser LogReport GetAttr.temp temp_flist
  rm -f modisdet* goge* global_emiss* oisst*
  rm -f pc*
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
