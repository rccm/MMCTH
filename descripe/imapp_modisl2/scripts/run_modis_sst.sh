#!/bin/bash
#
# Run the IMAPP MODIS sea surface temperature software.
#
# You must first set up the environmental variables to match your system and
# directory structure prior to execution.  These variables can be 
# automatically set by editing the one of the scripts found in the 
# imapp_modisl2/env directory.  Set the top level
# directory for the modis level 2 processing path and then type:
#  source imapp_modisl2.bash_env for bash shells or
#  source imapp_modisl2.csh_env for csh shells.
#
# The variables that need to be set to run this code are:
#
# MODIS_L2_HOME - home directory for MODIS Level 2 processing
# LOCAL_ANC_DIR - Local directory for ancillary data.  The default is:
#                 ${MODIS_L2_HOME}/ancillary
# REMOTE_ANC_DIR - Remote directory from which to fetch ancillary data if
#                 it is not found locally.  The default is: 
#                 ftp://aqua.ssec.wisc.edu/pub/terra/ancillary
# MODIS_L2_CFG   - Directory where the configuration files can be found. The 
#                 default is:  ${MODIS_L2_HOME}/config
# MODIS_L2_COEFF - Directory where the software coefficients can be found. The
#                 default is: ${MODIS_L2_HOME}/coeff
# MODIS_L2_BIN   - Directory where the statically linked pre-compiled binaries
#                 can be found.  The default is: ${MODIS_L2_HOME}/bin
#  In addtion, sourcing the bash or csh environmental scripts will set the path
#  so that the binaries can be executed from any directory.
# PATH=.:${MODIS_L2_BIN}:${MODIS_L2_HOME}/scripts:${PATH} 
#
# Also new with this release is the option to create binary, hdf or both types
# of output files.  With this inclusion, the IDL software to convert the binary
# to hdf is no longer needed.  
#
#------------------------------------------------------------------------------

echo
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running MODIS SST" 
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

# Check arguments
if [ $# != 5 ]; then
  echo "Usage: run_modis_sst.sh SAT OUTPUT_TYPE FIL1KM OUTDIR COMP"
  echo "where"
  echo "SAT    is the satellite name (Aqua or Terra)"
  echo "OUTPUT_TYPE is the format of the output product file"
  echo "   1 = binary only, 2 = hdf only, 3 = binary and hdf"
  echo "FIL1KM is the MODIS L1B 1000 meter resolution radiance flat file"
  echo "OUTDIR is the directory where the products will be stored"
  echo "COMP is set to 1 to compress output product, 0 otherwise"
  exit -1
fi
#------------------------------------------------------------------------------

# Extract file names
SAT=$1
OUTPUT_TYPE=$2
FIL1KM=$3
OUTDIR=$4
COMP=$5

# Set root of output file name (e.g. 't1.02001.1815')
BASE_FIL1KM=`basename $FIL1KM`
ROOT=`echo $BASE_FIL1KM:t | cut -d. -f1,2,3`
#ROOT=$FIL1KM:r
#ROOT=$ROOT:r
#ROOT=$ROOT:t

# Get year and date for ancillary data extraction
BASE_FIL1KM=`basename $FIL1KM`
# Get date and time string from L1B 1KM meter file
DATE_TIME=`echo $BASE_FIL1KM:t | cut -d. -f2,3 | cut -dA -f2`
DATE=`echo $DATE_TIME | cut -d. -f1 | cut -dA -f2`

DATE_ORIG=$DATE
if [ ${#DATE} -eq 5 ]; then
  DATE=20$DATE
fi

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Getting ancillary data for MODIS sst processing"

# get required ancillary data - Reynolds Blended SST
# ----------------------------------------------------------------
# Run script which checks local directory first, and then ftp site
RSST=`get_anc_rsst_wday.sh $DATE`
# If you can't find a file, set it equal to missing
if [ $? != 0 ]; then RSST="MISSING" ; fi

# Convert backslashes so that they will be properly inserted using sed
ROISST=`echo ${RSST} | sed s/"\/"/"\\\/"/g`
ROISST=`basename $ROISST`

# Link ancillary file to current directory
ln -s $RSST .

# Print start message for processing log
echo ""
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Started MODIS sea surface temperature processing at "`date`

# Set sea surface temperature output file names
OUTSSTFIL=$ROOT.mod28.img
OUTSSTHDR=$ROOT.mod28.hdr
OUTSSTHDF=$ROOT.mod28.hdf
SSTDEBUG=$ROOT.mod28debug.txt

# Create a configuration file from the template
TEMPLATE=${MODIS_L2_CFG}/sst.cfg
/bin/sed \
  -e "s/FIL1KM_IMG/${ROOT}.1000m.img/g" \
  -e "s/FIL1KM_HDR/${ROOT}.1000m.hdr/g" \
  -e "s/FILGEO_IMG/${ROOT}.geo.img/g" \
  -e "s/FILGEO_HDR/${ROOT}.geo.hdr/g" \
  -e "s/MOD28.img/${OUTSSTFIL}/g" \
  -e "s/MOD28.hdr/${OUTSSTHDR}/g" \
  -e "s/MOD28.hdf/${OUTSSTHDF}/g" \
  -e "s/sstdebug.dat/${SSTDEBUG}/g" \
  -e "s/oisst_weekly/${ROISST}/g" \
  $TEMPLATE > sst.cfg

# Run the sea surface temperature algorithm
sst.exe sst.cfg $SAT $OUTPUT_TYPE

PGE_STATUS=$?
# Clean up
if [ $PGE_STATUS == 0 ]; then
#  rm -f $FILPCF ShmMem LogUser LogReport LogStatus GetAttr.temp
  rm -f ShmMem LogUser LogReport GetAttr.temp
  rm -f $ROISST
fi

# Rename the configuration file
/bin/mv sst.cfg $ROOT.mod28.cfg

# Set name of SST HDF file
OUTSSTHDF=$ROOT.mod28.hdf

# Move the output files to the product directory
/bin/mv $SSTDEBUG $ROOT.mod28.cfg $OUTDIR

if [ $OUTPUT_TYPE == 1 ] || [ $OUTPUT_TYPE == 3 ]; then
    /bin/mv $OUTSSTFIL $OUTSSTHDR $OUTDIR
fi

if [ $OUTPUT_TYPE == 2 ] || [ $OUTPUT_TYPE == 3 ]; then
    echo "Renaming product files to direct broadcast convention..."

    if [ $SAT == "terra" ] || [ $SAT == "Terra" ] || [ $SAT == "TERRA" ]; then
      DB_HEADER="t1"
    elif [ $SAT == "aqua" ] || [ $SAT == "Aqua" ] || [ $SAT == "AQUA" ]; then
      DB_HEADER="a1"
    else
      echo "Satellite name not recognized: "$SAT
      exit 1
    fi

    if [ `echo $DATE_TIME | cut -c1-2` == "XX" ]; then
       DB_DATETIME=`echo $DATE_TIME | cut -c3-12`
    else
      DB_DATETIME=$DATE_TIME
    fi
    DB_SSTFIL=${DB_HEADER}.${DB_DATETIME}.mod28.hdf

    /bin/mv $OUTSSTHDF $DB_SSTFIL

    if [ "$COMP" == "1" ]; then
      echo
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo "Compressing output files..."
      hrepack.sh $DB_SSTFIL
    fi

    /bin/mv $DB_SSTFIL $OUTDIR
fi

# Print end message for processing log
echo ""
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Finished MODIS sea surface temperature processing at "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo

# Exit gracefully
exit 0
