#!/bin/bash

echo
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "MOD04 aerosol processing started at "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

# Set top level variables
ROOT=$MODIS_L2_HOME
MDIR=$ROOT/src/aerosol

# Check arguments
if [ $# != 10 ]; then
  echo "Usage: run_modis_aerosol.sh SAT FIL1KM FILQKM FILGEO GDAS1 GDAS2 RSST NISE ICEC REAGG COMP"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  echo "FILHKM is the Level 1B  500 meter radiance HDF file (optional: enter MISSING if not used)"
  echo "FILQKM is the Level 1B  250 meter radiance HDF file (optional: enter MISSING if not used)"
  echo "FILGEO is the Level 1B 1000 meter geolocation HDF file"
  echo "FILPROF is the Level 2 retrieved profiles HDF file"
  echo "FILCM is the Level 2 cloud mask HDF file"
  echo "OUTDIR is the output products directory"
  echo "REAGG is set to 1 to create re-aggregated products, 0 otherwise"
  echo "COMP is set to 1 to compress output product, 0 otherwise"
  exit 1
fi

# Extract arguments
SAT=$1
FIL1KM=`basename $2`
FILHKM=`basename $3`
FILQKM=`basename $4`
FILGEO=`basename $5`
FILPROF=`basename $6`
FILCM=`basename $7`

# directory names
DIR1KM=`dirname $2`
DIRHKM=`dirname $3`
DIRQKM=`dirname $4`
DIRGEO=`dirname $5`
DIRPROF=`dirname $6`
DIRCM=`dirname $7`

OUTDIR=$8
REAGG=$9
COMP=${10}

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
# CHECK FOR AEROSOL COEFFICIENT DIRECTORY
#-------------------------------------------------------------------------------

if [ ! -e $MODIS_L2_HOME/src/aerosol/coeff ]; then
  echo
  echo "WARNING: No MOD04 aerosol coefficient directory found. MOD04 Aerosol"   
  echo " processing cannot be run."
  exit 1
fi

#-------------------------------------------------------------------------------
# SET UP INPUT FILES
#-------------------------------------------------------------------------------

# Update leapsec.dat and utcpole.dat files if more than 7 days old
# ----------------------------------------------------------------
get_anc_leapsec_utcpole.bash

# Create links to static input files
ln -f -s $MDIR/SHR_MCF/${HEADER}* .

# links for MOD_PR04_05
ln -f -s $MDIR/MOD_PR04_05/src/${HEADER}*.mcf .

# links for MOD04_3K
ln -f -s $MDIR/MOD04_3K/src/${HEADER}*.mcf

#-------------------------------------------------------------------------------
# Copy leapsec.dat and utcpole.dat files to run directory
#-------------------------------------------------------------------------------

cp $LOCAL_ANC_DIR/util/leapsec.dat .
cp $LOCAL_ANC_DIR/util/utcpole.dat .

#-------------------------------------------------------------------------------
# GET ANCILLARY DATA 
#-------------------------------------------------------------------------------
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Getting ancillary data for MODIS aerosol processing"
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

echo "GDAS = " $GDASnearest

FILGDAS=`basename $GDASnearest`
DIRGDAS=`dirname $GDASnearest`

# TOAST file
# ----------------------------------------------------------------
# GDAS1 numerical weather prediction model grid field analysis
# ----------------------------------------------------------------
echo "Date and time from L1B 1KM file: "$DATE $TIME
TOASTfile=`get_anc_toast_ozone.sh $DATE`
status=$?

FILTOAST=`basename $TOASTfile`
DIRTOAST=`dirname $TOASTfile`

#-------------------------------------------------------------------------------
# SET UP PCF FILE
#-------------------------------------------------------------------------------

# Set current date/time
CURRENT=`date -u +%Y%j%H%M%S`

# Set names of output files
FILOUT_MOD04_L2=${HEADER}04_L2.$DATE_TIME.hdf
FILOUT_MOD04_QC=${HEADER}04_QC.$DATE_TIME.txt
FILOUT_MOD05_L2=${HEADER}05_L2.$DATE_TIME.hdf
FILOUT_MOD05_QC=${HEADER}05_QC.$DATE_TIME.txt
FILOUT_MOD05C_QC=${HEADER}05C_QC.$DATE_TIME.txt
FILOUT_MOD04_3K=${HEADER}04_3K.$DATE_TIME.hdf

# Set dummy start and stop times for PCF
START_TIME='2000-01-01T00:00:00'
STOP_TIME=`date -u +%Y-%m-%dT%H:%M:%S`

#-------------------------------------------------------------------------------
# SET UP MOD04_3K PCF FILE
#-------------------------------------------------------------------------------

# Set name of template file
TEMPLATE=$MDIR/MOD04_3K/template/${HEADER}04_3K.pcf.template

# Set name of new PCF file
FILPCF_MOD04_3K=${HEADER}_PR04_3K.$DATE_TIME.pcf

# Set coefficient directory path
COEFDIR=${MDIR}/coeff

# Create new PCF file from the template
sed \
  -e "s?COEFDIR?${COEFDIR}?g" \
  -e "s?FIL1KM?${FIL1KM:t}?g" \
  -e "s?DIR1KM?${DIR1KM}?g" \
  -e "s?FILHKM?${FILHKM}?g" \
  -e "s?DIRHKM?${DIRHKM}?g" \
  -e "s?FILQKM?${FILQKM}?g" \
  -e "s?DIRQKM?${DIRQKM}?g" \
  -e "s?FILGEO?${FILGEO}?g" \
  -e "s?DIRGEO?${DIRGEO}?g" \
  -e "s?FILPROF?${FILPROF}?g" \
  -e "s?DIRPROF?${DIRPROF}?g" \
  -e "s?FILCM?${FILCM}?g" \
  -e "s?DIRCM?${DIRCM}?g" \
  -e "s?FILGDAS?${FILGDAS}?g" \
  -e "s?DIRGDAS?${DIRGDAS}?g" \
  -e "s?FILTOAST?${FILTOAST}?g" \
  -e "s?DIRTOAST?${DIRTOAST}?g" \
  -e "s?FILOUT_MOD04_L2?${FILOUT_MOD04_L2}?g" \
  -e "s?FILOUT_MOD05_L2?${FILOUT_MOD05_L2}?g" \
  -e "s?FILOUT_MOD04_3K?${FILOUT_MOD04_3K}?g" \
  -e "s?FILOUT_MOD04_QC?${FILOUT_MOD04_QC}?g" \
  -e "s?FILOUT_MOD05_QC?${FILOUT_MOD05_QC}?g" \
  -e "s?FILOUT_MOD05C_QC?${FILOUT_MOD05C_QC}?g" \
  -e "s?START_TIME?${START_TIME}?g" \
  -e "s?STOP_TIME?${STOP_TIME}?g" \
   $TEMPLATE > $FILPCF_MOD04_3K

#-------------------------------------------------------------------------------
# SET UP MOD_PR04DB PCF FILE
#-------------------------------------------------------------------------------

# Set name of template file
TEMPLATE=$MDIR/MOD_PR04DB/template/${HEADER}_PR04DB.pcf.template

# Set name of new PCF file
FILPCF_MOD04_DB=${HEADER}_PR04DB.$DATE_TIME.pcf

# Create new PCF file from the template
sed \
  -e "s?COEFDIR?${COEFDIR}?g" \
  -e "s?FIL1KM?${FIL1KM:t}?g" \
  -e "s?DIR1KM?${DIR1KM}?g" \
  -e "s?FILHKM?${FILHKM}?g" \
  -e "s?DIRHKM?${DIRHKM}?g" \
  -e "s?FILQKM?${FILQKM}?g" \
  -e "s?DIRQKM?${DIRQKM}?g" \
  -e "s?FILGEO?${FILGEO}?g" \
  -e "s?DIRGEO?${DIRGEO}?g" \
  -e "s?FILPROF?${FILPROF}?g" \
  -e "s?DIRPROF?${DIRPROF}?g" \
  -e "s?FILCM?${FILCM}?g" \
  -e "s?DIRCM?${DIRCM}?g" \
  -e "s?FILGDAS?${FILGDAS}?g" \
  -e "s?DIRGDAS?${DIRGDAS}?g" \
  -e "s?FILTOAST?${FILTOAST}?g" \
  -e "s?DIRTOAST?${DIRTOAST}?g" \
  -e "s?FILOUT_MOD04_L2?${FILOUT_MOD04_L2}?g" \
  -e "s?FILOUT_MOD05_L2?${FILOUT_MOD05_L2}?g" \
  -e "s?FILOUT_MOD04_3K?${FILOUT_MOD04_3K}?g" \
  -e "s?FILOUT_MOD04_QC?${FILOUT_MOD04_QC}?g" \
  -e "s?FILOUT_MOD05_QC?${FILOUT_MOD05_QC}?g" \
  -e "s?FILOUT_MOD05C_QC?${FILOUT_MOD05C_QC}?g" \
  -e "s?START_TIME?${START_TIME}?g" \
  -e "s?STOP_TIME?${STOP_TIME}?g" \
   $TEMPLATE > $FILPCF_MOD04_DB
  
#-------------------------------------------------------------------------------
# RUN THE ALGORITHM
#-------------------------------------------------------------------------------

# Set toolkit environment variables
export PGS_PC_INFO_FILE=./$FILPCF_MOD04_3K

# Print start message
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo $PROD "Aerosol processing started at "`date`

# Run MOD_PR04CR1 - create the HDF skeleton file 
if [ $REAGG == 1 ]; then
  $ROOT/bin/createaerosolfile1.exe $REAGG
else
  $ROOT/bin/createaerosolfile1.exe
fi

if [ $? != 0 ]; then
  echo "********** Error encountered running MOD_PR04CR1 **********"
  exit 1
else
  echo "********** Successful completion of MOD_PR04CR1 **********"
fi 
  
# Run MOD_PR04_05 - procuce MODIS L2 Aerosol (MOD04 L2) and NIR Precipitable Water (MOD05 L2) products 
if [ $REAGG == 1 ]; then
  $ROOT/bin/aerosol.exe $REAGG
else
  $ROOT/bin/aerosol.exe
fi

if [ $? != 0 ]; then
  echo "********** Error encountered running MOD_PR04_05 **********"
  exit 1
else
  echo "********** Successful completion of MOD_PR04_05 **********"
fi

# Run MOD_PR04CR2 - create HDF skeleton file for MOD04 3KM
if [ $REAGG == 1 ]; then
  $ROOT/bin/createaerosolfile2.exe $REAGG
else
  $ROOT/bin/createaerosolfile2.exe
fi

if [ $? != 0 ]; then
  echo "********** Error encountered running MOD_PR04CR2 **********"
  exit 1
else
  echo "********** Successful completion of MOD_PR04CR2 **********"
fi

# Run MOD04_3K - produces MODIS L2 Aerosol at higher 3km resolution (no MOD05 product)
if [ $REAGG == 1 ]; then
  $ROOT/bin/aerosol3km.exe $REAGG
else
  $ROOT/bin/aerosol3km.exe
fi

if [ $? != 0 ]; then
  echo "********** Error encountered running MOD04_3K **********"
  exit 1
else
  echo "********** Successful completion of MOD04_3K **********"
fi

# Run MOD_PR04DB (DeepBlue) - produce Deep Blue update to MODIS L2 aerosol product
export PGS_PC_INFO_FILE=./$FILPCF_MOD04_DB

#$ROOT/src/aerosol/MOD_PR04DB/src/MOD_PR04DB.exe
if [ $REAGG == 1 ]; then
  $ROOT/bin/aerosoldeepblue.exe $REAGG
else
  $ROOT/bin/aerosoldeepblue.exe
fi

if [ $? != 0 ]; then
  echo "********** Error encountered running MOD_PR04DB **********"
  exit 1
else
  echo "********** Successful completion of MOD_PR04DB **********"
fi

PGE_STATUS=$?

# Clean up
if [ $PGE_STATUS == 0 ]; then
  rm -f $FILPCF_MOD04_3K $FILPCF_MOD04_DB
  rm -f $FILOUT_MOD04_QC $FILOUT_MOD05_QC $FILOUT_MOD05C_QC
  rm -f $FILOUT_MOD04_L2.met $FILOUT_MOD04_3K.met $FILOUT_MOD04_QC.met $FILOUT_MOD05_L2.met $FILOUT_MOD05_QC.met $FILOUT_MOD05C_QC.met
  rm -f MOD04_DB_temp.bin
  rm -f ShmMem LogUser LogReport GetAttr.temp temp_flist #LogStatus
  rm -f ${HEADER}*.mcf ${HEADER}*.MCF
  rm -f pc*
  rm -f fort.36
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

/bin/mv $FILOUT_MOD04_L2 $OUTDIR/${DB_HEADER}.${DB_DATETIME}.mod04.hdf
/bin/mv $FILOUT_MOD04_3K $OUTDIR/${DB_HEADER}.${DB_DATETIME}.mod04_3k.hdf 

if [ "$COMP" == "1" ]; then
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Compressing output file..."
  hrepack.sh ${DB_HEADER}.${DB_DATETIME}.mod07.hdf
fi

# Print finish message
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "MOD04 aerosol processing ended at  "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo
exit 0
