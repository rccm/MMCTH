#!/bin/bash

#---------------  MODIS DESTRIPING RUN SCRIPT --------------------------------
#                      21 September  2017
#
#  Significant revisions based upon a re-evaluation of the behavior
#  of the Terra instrument over its lifetime.  No changes have been
#  made to the destriping source code, but there are now 24 
#  Terra coefficient files.
#
#---------------  MODIS DESTRIPING RUN SCRIPT --------------------------------
#                      19 September  2016
#
#  Updated version 4 and verion 5 of Terra coefficient files.  Version 5 added
#  as a consequence of the Terra Safe Hold event that occurred in February 2016.
#
#

#-------------------------------------------------------------------------------
# SETUP AND CHECK ARGUMENTS
#-------------------------------------------------------------------------------

echo
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running destriping code..."

# Set top level variables
ROOT=$MODIS_L2_HOME
PROD=destripe
MDIR=$ROOT/src/$PROD

# Check arguments
if [ $# != 2 ]; then
  echo "Usage: run_modis_destripe.sh SAT FIL1KM"
  echo "where"
  echo "SAT    is the satellite name (terra or aqua)"
  echo "FIL1KM is the Level 1B 1000 meter radiance HDF file"
  exit -1
fi

# Extract arguments
SAT=$1
FIL1KM=`basename $2`
DIR1KM=`dirname  $2`

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
# CREATE LINKS TO STATIC INPUT FILES
#-------------------------------------------------------------------------------

# Create links to static input files
ln -f -s $MDIR/coeff/destripe_config_${SAT}.dat.* .

#-------------------------------------------------------------------------------
# SET UP PCF FILE
#-------------------------------------------------------------------------------

# Get date and time string from L1B 1KM meter file
DATE_TIME=`echo $FIL1KM:t | cut -d. -f2,3`

DATA_DATE=`echo $FIL1KM:t | cut -d. -f2`

FILETYPE=`echo $FIL1KM | cut -c1-3`
if [[ $FILETYPE == "MOD" ]] || [[ $FILETYPE == "MYD" ]]; then
   DATA_YEAR=`echo $DATA_DATE | cut -c2-5`
   JDAY1=`echo $DATA_DATE | cut -c6-8`
elif [[ $FILETYPE == "t1." ]] || [[ $FILETYPE == "a1." ]]; then
   DATA_YEAR=20`echo $DATA_DATE | cut -c1-2`
   JDAY1=`echo $DATA_DATE | cut -c3-5`
fi

DATA_JDAY=`expr $JDAY1  -  1`
DATA_MONTH=`date --date="$DATA_YEAR-01-01 + $DATA_JDAY days" "+%m"`
DATA_DAY=`date --date="$DATA_YEAR-01-01 + $DATA_JDAY days" "+%d"`
EPOCH_TIME=`date --date="$DATA_YEAR-$DATA_MONTH-$DATA_DAY" "+%s"` 
EPOCH_DATE=`expr $EPOCH_TIME/86400 | bc`


# Logistics to figure out which coefficient file to use
if [ $SAT == "aqua" ]; then
   TIME_V2=`date --date="2005-01-10" "+%s"`
   DATA_V2=`expr $TIME_V2/86400 | bc`

   if [ $EPOCH_DATE -lt $DATA_V2 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v1"
   else
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v2"
   fi

elif [ $SAT == "terra" ]; then

   TIME_V1=`date --date="2000-02-24" "+%s"`
   DATA_V1=`expr $TIME_V1/86400 | bc`

   TIME_V2=`date --date="2001-01-19" "+%s"`
   DATA_V2=`expr $TIME_V2/86400 | bc`

   TIME_V3=`date --date="2002-03-19" "+%s"`
   DATA_V3=`expr $TIME_V3/86400 | bc`

   TIME_V4=`date --date="2003-12-16" "+%s"`
   DATA_V4=`expr $TIME_V4/86400 | bc`

   TIME_V5=`date --date="2004-04-27" "+%s"`
   DATA_V5=`expr $TIME_V5/86400 | bc`

   TIME_V6=`date --date="2004-06-23" "+%s"`
   DATA_V6=`expr $TIME_V6/86400 | bc`

   TIME_V7=`date --date="2005-10-01" "+%s"`
   DATA_V7=`expr $TIME_V7/86400 | bc`

   TIME_V8=`date --date="2005-11-05" "+%s"`
   DATA_V8=`expr $TIME_V8/86400 | bc`

   TIME_V9=`date --date="2006-02-19" "+%s"`
   DATA_V9=`expr $TIME_V9/86400 | bc`

   TIME_V10=`date --date="2006-06-04" "+%s"`
   DATA_V10=`expr $TIME_V10/86400 | bc`

   TIME_V11=`date --date="2006-08-29" "+%s"`
   DATA_V11=`expr $TIME_V11/86400 | bc`

   TIME_V12=`date --date="2007-07-12" "+%s"`
   DATA_V12=`expr $TIME_V12/86400 | bc`

   TIME_V13=`date --date="2008-03-01" "+%s"`
   DATA_V13=`expr $TIME_V13/86400 | bc`

   TIME_V14=`date --date="2008-11-03" "+%s"`
   DATA_V14=`expr $TIME_V14/86400 | bc`

   TIME_V15=`date --date="2010-12-01" "+%s"`
   DATA_V15=`expr $TIME_V15/86400 | bc`

   TIME_V16=`date --date="2011-07-01" "+%s"`
   DATA_V16=`expr $TIME_V16/86400 | bc`

   TIME_V17=`date --date="2012-05-16" "+%s"`
   DATA_V17=`expr $TIME_V17/86400 | bc`

   TIME_V18=`date --date="2012-07-01" "+%s"`
   DATA_V18=`expr $TIME_V18/86400 | bc`

   TIME_V19=`date --date="2013-05-05" "+%s"`
   DATA_V19=`expr $TIME_V19/86400 | bc`

   TIME_V20=`date --date="2013-09-08" "+%s"`
   DATA_V20=`expr $TIME_V20/86400 | bc`

   TIME_V21=`date --date="2014-04-15" "+%s"`
   DATA_V21=`expr $TIME_V21/86400 | bc`

   TIME_V22=`date --date="2014-08-10" "+%s"`
   DATA_V22=`expr $TIME_V22/86400 | bc`

   TIME_V23=`date --date="2015-04-02" "+%s"`
   DATA_V23=`expr $TIME_V23/86400 | bc`

   TIME_V24=`date --date="2016-02-25" "+%s"`
   DATA_V24=`expr $TIME_V24/86400 | bc`

   if [ $EPOCH_DATE -ge $DATA_V24 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v24r2"
   elif [ $EPOCH_DATE -lt $DATA_V24 ] && [ $EPOCH_DATE -ge $DATA_V23 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v23r2"
   elif [ $EPOCH_DATE -lt $DATA_V23 ] && [ $EPOCH_DATE -ge $DATA_V22 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v22r2"
   elif [ $EPOCH_DATE -lt $DATA_V22 ] && [ $EPOCH_DATE -ge $DATA_V21 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v21r2"
   elif [ $EPOCH_DATE -lt $DATA_V21 ] && [ $EPOCH_DATE -ge $DATA_V20 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v20r2"
   elif [ $EPOCH_DATE -lt $DATA_V20 ] && [ $EPOCH_DATE -ge $DATA_V19 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v19r2"
   elif [ $EPOCH_DATE -lt $DATA_V19 ] && [ $EPOCH_DATE -ge $DATA_V18 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v18r2"
   elif [ $EPOCH_DATE -lt $DATA_V18 ] && [ $EPOCH_DATE -ge $DATA_V17 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v17r2"
   elif [ $EPOCH_DATE -lt $DATA_V17 ] && [ $EPOCH_DATE -ge $DATA_V16 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v16r2"
   elif [ $EPOCH_DATE -lt $DATA_V16 ] && [ $EPOCH_DATE -ge $DATA_V15 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v15r2"
   elif [ $EPOCH_DATE -lt $DATA_V15 ] && [ $EPOCH_DATE -ge $DATA_V14 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v14r2"
   elif [ $EPOCH_DATE -lt $DATA_V14 ] && [ $EPOCH_DATE -ge $DATA_V13 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v13r2"
   elif [ $EPOCH_DATE -lt $DATA_V13 ] && [ $EPOCH_DATE -ge $DATA_V12 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v12r2"
   elif [ $EPOCH_DATE -lt $DATA_V12 ] && [ $EPOCH_DATE -ge $DATA_V11 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v11r2"
   elif [ $EPOCH_DATE -lt $DATA_V11 ] && [ $EPOCH_DATE -ge $DATA_V10 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v10r2"
   elif [ $EPOCH_DATE -lt $DATA_V10 ] && [ $EPOCH_DATE -ge $DATA_V9 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v9r2"
   elif [ $EPOCH_DATE -lt $DATA_V9 ] && [ $EPOCH_DATE -ge $DATA_V8 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v8r2"
   elif [ $EPOCH_DATE -lt $DATA_V8 ] && [ $EPOCH_DATE -ge $DATA_V7 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v7r2"
   elif [ $EPOCH_DATE -lt $DATA_V7 ] && [ $EPOCH_DATE -ge $DATA_V6 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v6r2"
   elif [ $EPOCH_DATE -lt $DATA_V6 ] && [ $EPOCH_DATE -ge $DATA_V5 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v5r2"
   elif [ $EPOCH_DATE -lt $DATA_V5 ] && [ $EPOCH_DATE -ge $DATA_V4 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v4r2"
   elif [ $EPOCH_DATE -lt $DATA_V4 ] && [ $EPOCH_DATE -ge $DATA_V3 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v3r2"
   elif [ $EPOCH_DATE -lt $DATA_V3 ] && [ $EPOCH_DATE -ge $DATA_V2 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v2r2"
   elif [ $EPOCH_DATE -lt $DATA_V2 ]; then
      DESTRIPE_COEFF="destripe_config_${SAT}.dat.v1r2"
   fi

else
     echo "Satellite does not match aqua or terra " $SAT
     exit 1
fi

echo "Using coefficient file: " $DESTRIPE_COEFF

# Set current date/time
CURRENT=`date -u +%Y%j%H%M%S`

# Set name of template file
TEMPLATE=$MDIR/template/${HEADER}_PRDS.pcf.template

# Set name of new PCF file
FILPCF=${HEADER}_PRDS.$DATE_TIME.pcf

# Set dummy start and stop times for PCF
START_TIME='2000-01-01T00:00:00'
STOP_TIME=`date -u +%Y-%m-%dT%H:%M:%S`

# Create new PCF file from the template
sed \
  -e "s?FIL1KM?${FIL1KM}?g" \
  -e "s?DIR1KM?${DIR1KM}?g" \
  -e "s?DESTRIPE_CONFIG?${DESTRIPE_COEFF}?g" \
  -e "s?START_TIME?${START_TIME}?g" \
  -e "s?STOP_TIME?${STOP_TIME}?g" \
  $TEMPLATE > $FILPCF

#-------------------------------------------------------------------------------
# RUN THE ALGORITHM
#-------------------------------------------------------------------------------

# Set toolkit environment variables
export PGS_PC_INFO_FILE=./$FILPCF

# Print start message
echo $PROD "processing started at "`date`

# Run the algorithm
$ROOT/bin/destripe.exe
PGE_STATUS=$?

# Clean up and print failure or success message
if [ $PGE_STATUS != 0 ]; then
  echo "Error encountered"
  rm -f $FILPCF ShmMem
  rm -f destripe_config_${SAT}.dat.*
  exit 1
  echo "********** Failure running destriping code **********"
else
  rm -f destripe_config_${SAT}.dat.*
  echo "********** Successful completion **********"
fi

# Print finish message
echo $PROD "processing ended at  "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo
exit 0
