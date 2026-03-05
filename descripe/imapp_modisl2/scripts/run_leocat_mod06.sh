#!/bin/bash

ROOT=$MODIS_L2_HOME

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------
if [ $# != 4 ]; then
  echo "Usage: run_leocat_mod06.sh SAT FILC5 FILPCF FILLEO"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FILC5 is the full path of the MOD06CT C5 file with fields to be filled in for C6"
  echo "FILPCF is the full path of the MOD06CT C5 Process Control File (PCF)"
  echo "FILLEO is the full path of the Intermediate Ouput File resulting from a LEOCAT run"
  exit 1
fi

# Extract arguments
SAT=$1
FILC5=`basename $2`
DIRC5=`dirname  $2`
FILPCF=`basename $3`
DIRPCF=`dirname  $3`
FILLEO=$4

# Assume FILPCF is of the form  MYD_PR06CT.AYYYYDDD.HHMM.pcf
DATE=`echo $FILPCF | cut -d. -f2 | cut -dA -f2`
TIME=`echo $FILPCF | cut -d. -f3`

# Get the following information from the PCF file
GEOFILE=`awk '/600000/' $DIRPCF/$FILPCF | cut -d'|' -f2`
GEODIR=`awk '/600000/' $DIRPCF/$FILPCF | cut -d'|' -f3`

# GDASFILE is the nearest file used for the 5-km
# GDAS1FILE and GDAS2FILE are what's used by LEOCAT
GDASFILE=`awk '/900000/' $DIRPCF/$FILPCF | cut -d'|' -f2`
GDAS1FILE=`awk '/900001/' $DIRPCF/$FILPCF | cut -d'|' -f2`
GDAS2FILE=`awk '/900002/' $DIRPCF/$FILPCF | cut -d'|' -f2`
GDASDIR=`awk '/900000/' $DIRPCF/$FILPCF | cut -d'|' -f3`

NISEFILE=`awk '/900100/' $DIRPCF/$FILPCF | cut -d'|' -f2`
NISEDIR=`awk '/900100/' $DIRPCF/$FILPCF | cut -d'|' -f3`

RSSTFILE=`awk '/900030/' $DIRPCF/$FILPCF | cut -d'|' -f2`
RSSTDIR=`awk '/900030/' $DIRPCF/$FILPCF | cut -d'|' -f3`

FIL1KMDS=`awk '/430000/' $DIRPCF/$FILPCF | cut -d'|' -f2`
DIR1KMDS=`awk '/430000/' $DIRPCF/$FILPCF | cut -d'|' -f3`

FILCSRB=`awk '/477500/' $DIRPCF/$FILPCF | cut -d'|' -f2`
DIRCSRB=`awk '/477500/' $DIRPCF/$FILPCF | cut -d'|' -f3`

FILCM=`awk '/422500/' $DIRPCF/$FILPCF | cut -d'|' -f2`
DIRCM=`awk '/422500/' $DIRPCF/$FILPCF | cut -d'|' -f3`

SD="$ROOT/src/leocat/coeff/MODIS_NPOESS_STATIC/collection6"
WD=$SD

# Run the algorithm
$ROOT/scripts/run_leocat_noftp_mod06.sh $WD $SD $GDASDIR $GEODIR/$GEOFILE $DIR1KMDS/$FIL1KMDS $DIRCSRB/$FILCSRB $FILLEO $SAT $DATE $TIME $GDASDIR/$GDAS1FILE $GDASDIR/$GDAS2FILE $NISEDIR/$NISEFILE $RSSTDIR/$RSSTFILE $DIRCM/$FILCM

# Copy the C6 output to the C5 skeleton file
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Copying LEOCAT 1-km results to C6 MOD_PR06 output file..."
$ROOT/bin/leocat_cloudtop.exe $FILLEO $DIRC5/$FILC5

PGE_STATUS=$?

if [ $PGE_STATUS == 1 ]; then
  echo
  echo "********** WARNING: Error copying LEOCAT output to MOD06 skeleton file **********"
  echo
  exit 1
else
  echo
  echo "********** Success copying LEOCAT output to MOD06 skeleton file **********"
  echo
fi

# Clean up
if [ $PGE_STATUS == 0 ]; then
  rm -f $FILPCF $FILLEO
fi

# Print finish message
echo "run_leocat_mod06.sh script complete"
exit 0
