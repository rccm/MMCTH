#!/bin/bash

#
# Wrapper script for creating quicklook images for MODIS Level 2 products
# using python
#
# April 2014 Rebecca Cintineo
#

echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Starting quicklook creation script"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo

# Set top level variables
ROOT=$MODIS_L2_HOME

# Check arguments
if [ $# != 6 ]; then
  echo "Usage: run_quicklook.sh SAT FIL1KM FILQKM FILGEO GDAS1 GDAS2 RSST NISE ICEC REAGG"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "DATE_TIME is the date and time of the granule being plotted"
  echo "FILGEO is the Level 1B 1000 meter geolocation HDF file"
  echo "FIL1KM is the Level 1B 1000 meter HDF file"
  echo "OUTDIR is the product output directory"
  echo "POLAR is TRUE if the granule is a polar granule and FALSE otherwise"
  exit 1
fi

SAT=$1
DATE_TIME=$2
FILGEO=$3
FIL1KM=$4
OUTDIR=$5
POLAR=$6

if [ "$SAT" == "TERRA" ] || [ "$SAT" == "terra" ] || [ "$SAT" == "Terra" ]; then
  DB_PREFIX="t1" 
elif [ "$SAT" == "AQUA" ] || [ "$SAT" == "aqua" ] || [ "$SAT" == "Aqua" ]; then
  DB_PREFIX="a1"
fi
  
if [ `echo $DATE_TIME | cut -c1-2` == "XX" ]; then
   DB_DATETIME=`echo $DATE_TIME | cut -c3-12`
else 
  DB_DATETIME=$DATE_TIME
fi
  
FILROOT=$DB_PREFIX.$DB_DATETIME

if [ $POLAR == "TRUE" ]; then

  echo ${MODIS_L2_HOME}/ShellB3/bin/python -m IMAPPquicklooks -root $FILROOT -outdir $OUTDIR -geo $FILGEO -l1b1km $FIL1KM -imagedir $MODIS_L2_IMAGES -res $MODIS_L2_RES -polar

  ${MODIS_L2_HOME}/ShellB3/bin/python -m IMAPPquicklooks -root $FILROOT -outdir $OUTDIR -geo $FILGEO -l1b1km $FIL1KM -imagedir $MODIS_L2_IMAGES -res $MODIS_L2_RES -polar

else

  echo ${MODIS_L2_HOME}/ShellB3/bin/python -m IMAPPquicklooks -root $FILROOT -outdir $OUTDIR -geo $FILGEO -l1b1km $FIL1KM -imagedir $MODIS_L2_IMAGES -res $MODIS_L2_RES

  ${MODIS_L2_HOME}/ShellB3/bin/python -m IMAPPquicklooks -root $FILROOT -outdir $OUTDIR -geo $FILGEO -l1b1km $FIL1KM -imagedir $MODIS_L2_IMAGES -res $MODIS_L2_RES

fi
