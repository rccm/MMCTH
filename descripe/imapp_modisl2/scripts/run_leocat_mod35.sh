#!/bin/bash

ROOT=$MODIS_L2_HOME

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------
if [ $# != 5 ]; then
  echo "Usage: run_leocat_mod35.sh SAT FILC5 FILPCF FILLEO"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FILC5 is the full path of the skeleton MOD35 C5 file with fields to be filled in for C6"
  echo "FILPCF is the full path of the MOD35 C5 Process Control File (PCF)"
  echo "FILLEO is the full path of the Intermediate Ouput File resulting from a LEOCAT run"
  echo "REAGG is set to 1 to create re-aggregated products, 0 otherwise"
  exit 1
fi

# Extract arguments
SAT=$1
FILC5=`basename $2`
DIRC5=`dirname  $2`
FILPCF=`basename $3`
DIRPCF=`dirname  $3`
FILLEO=$4
REAGG=$5

# Assume FILPCF is of the form  MYD_PR35.AYYYYDDD.HHMM.pcf
DATE=`echo $FILPCF | cut -d. -f2 | cut -dA -f2`
TIME=`echo $FILPCF | cut -d. -f3`

# Get the following information from the PCF file
GEOFILE=`awk '/600000/' $DIRPCF/$FILPCF | cut -d'|' -f2`
GEODIR=`awk '/600000/' $DIRPCF/$FILPCF | cut -d'|' -f3`

# GDASFILE is the nearest file used by MOD_PR35 skeleton
#  GDAS1FILE and GDAS2FILE are what's used by LEOCAT
GDASFILE=`awk '/900000/' $DIRPCF/$FILPCF | cut -d'|' -f2`
GDAS1FILE=`awk '/900001/' $DIRPCF/$FILPCF | cut -d'|' -f2`
GDAS2FILE=`awk '/900002/' $DIRPCF/$FILPCF | cut -d'|' -f2`
GDASDIR=`awk '/900000/' $DIRPCF/$FILPCF | cut -d'|' -f3`
GDASDIR1=`awk '/900001/' $DIRPCF/$FILPCF | cut -d'|' -f3`
GDASDIR2=`awk '/900002/' $DIRPCF/$FILPCF | cut -d'|' -f3`

NISEFILE=`awk '/900100/' $DIRPCF/$FILPCF | cut -d'|' -f2`
NISEDIR=`awk '/900100/' $DIRPCF/$FILPCF | cut -d'|' -f3`

RSSTFILE=`awk '/900030/' $DIRPCF/$FILPCF | cut -d'|' -f2`
RSSTDIR=`awk '/900030/' $DIRPCF/$FILPCF | cut -d'|' -f3`

FIL1KMDS=`awk '/430000/' $DIRPCF/$FILPCF | cut -d'|' -f2`
DIR1KMDS=`awk '/430000/' $DIRPCF/$FILPCF | cut -d'|' -f3`

# FILQKM could be missing
FILQKM=`awk '/700000/' $DIRPCF/$FILPCF | cut -d'|' -f2`
DIRQKM=`awk '/700000/' $DIRPCF/$FILPCF | cut -d'|' -f3`

echo "FILQKM is $FILQKM"
if [ $FILQKM == "" ]; then
   QKMARG="MISSING"
else
   QKMARG=${DIRQKM}/${FILQKM}
fi
echo "QKMARG is $QKMARG"

SD="$ROOT/src/leocat/coeff/MODIS_NPOESS_STATIC/collection6"
WD=$SD

#-------------------------------------------------------------------------------
# RUN THE ALGORITHMS
#-------------------------------------------------------------------------------
$ROOT/scripts/run_leocat_noftp_mod35.sh $WD $SD $GDASDIR $GEODIR/$GEOFILE $DIR1KMDS/$FIL1KMDS $QKMARG $FILLEO $SAT $DATE $TIME $GDASDIR1/$GDAS1FILE $GDASDIR2/$GDAS2FILE $NISEDIR/$NISEFILE $RSSTDIR/$RSSTFILE

# Copy the C6 output to the C5 skeleton file
$ROOT/bin/leocat_cloudmask.exe $FILLEO $DIRC5/$FILC5

if [ $? == 1 ]; then
  echo
  echo "********** WARNING: Error copying LEOCAT output to MOD35 skeleton file **********"
  echo
  exit 1
else
  echo
  echo "********** Success copying LEOCAT output to MOD35 skeleton file **********"
  echo
fi

# Write the statistics to the C5 skeleton file
$ROOT/scripts/run_alg29stats.sh $SAT $FILLEO $DIRC5/$FILC5 $DIRPCF/$FILPCF $REAGG

# Print finish message
echo "Script complete"
exit 0
