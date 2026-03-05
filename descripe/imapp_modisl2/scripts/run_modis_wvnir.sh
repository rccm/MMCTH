#!/bin/bash
#
# Adapted for NIR water vapor by Peter Albert September 2004
# peter.albert@gmx.dw
#
# Run the IMAPP MODIS algorithm for retrieving total water vapor
# from visilble channels.
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

echo
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Running MODIS WVNIR" 
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

#-------------------------------------------------------------------------------
# SETUP
#-------------------------------------------------------------------------------

# Check arguments
if [ $# != 5 ]; then
  echo "Usage: modis_wvnir.sh SAT OUTPUT_TYPE FIL1KM OUTDIR COMP"
  echo "where"
  echo "SAT is the MODIS platform name (Terra or Aqua)"
  echo "OUTPUT_TYPE is the format of the output product file"
  echo "  1 = binary only, 2 = hdf only, 3 = binary and hdf"
  echo "FIL1KM is the MODIS L1B 1000 meter resolution radiance flat file"
  echo "OUTDIR is the directory where the products will be stored"
  echo "COMP is set to 1 to compress output product, 0 otherwise"
  exit -1
fi

# Extract file names
SAT=$1
OUTPUT_TYPE=$2
FIL1KM=$3
OUTDIR=$4
COMP=$5

# Set root of output file name (e.g. 't1.02001.1815')
BASE_FIL1KM=`basename $FIL1KM`
ROOT=`echo $BASE_FIL1KM:t | cut -d. -f1,2,3`
#ROOT=`echo $ROOT:r`
#ROOT=`echo $ROOT:t`

#------------------------
# RUN THE WVNIR ALGORITHM
#------------------------

# Print start message for processing log
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Started IMAPP MODIS NIR water vapor processing at "`date`

# Set WVNIR output file names
OUTNIRWVFIL=$ROOT.wvnir.img
OUTNIRWVHDR=$ROOT.wvnir.hdr
OUTWVNIRHDF=$ROOT.wvnir.hdf

# Create a configuration file from the template
TEMPLATE=${MODIS_L2_CFG}/wvnir.cfg
/bin/sed \
  -e "s/FIL1KM_IMG/${ROOT}.1000m.img/g" \
  -e "s/FIL1KM_HDR/${ROOT}.1000m.hdr/g" \
  -e "s/FILGEO_IMG/${ROOT}.geo.img/g" \
  -e "s/FILGEO_HDR/${ROOT}.geo.hdr/g" \
  -e "s/FILWVNIR_IMG/${OUTNIRWVFIL}/g" \
  -e "s/FILWVNIR_HDR/${OUTNIRWVHDR}/g" \
  -e "s/FILWVNIR_HDF/${OUTWVNIRHDF}/g" \
  $TEMPLATE > wvnir.cfg

# Run the wvnir executable
wvnir.exe wvnir.cfg $OUTPUT_TYPE

# Move the output files to the product directory
/bin/mv wvnir.cfg $OUTDIR/${ROOT}.wvnir.cfg

if [ $OUTPUT_TYPE == 1 ] || [ $OUTPUT_TYPE == 3 ]; then
    /bin/mv $OUTNIRWVFIL $OUTNIRWVHDR $OUTDIR
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

    # Get date and time string from L1B 1KM meter file
    BASE_FIL1KM=`basename $FIL1KM`
    DATE_TIME=`echo $BASE_FIL1KM:t | cut -d. -f2,3 | cut -dA -f2`
    if [ `echo $DATE_TIME | cut -c1-2` == "XX" ]; then
      DB_DATETIME=`echo $DATE_TIME | cut -c3-12`
    else
      DB_DATETIME=$DATE_TIME
    fi
    DB_OUTFIL=${DB_HEADER}.${DB_DATETIME}.wvnir.hdf

    /bin/mv $OUTWVNIRHDF $DB_OUTFIL

    if [ "$COMP" == "1" ]; then
      echo
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo "Compressing output files..."
      hrepack.sh $DB_OUTFIL 
    fi

    /bin/mv $DB_OUTFIL $OUTDIR
fi

#-------------------------------------------------------------------------------
# CLEANUP AND EXIT
#-------------------------------------------------------------------------------

# Print end message for processing log
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Finished IMAPP MODIS NIR water vapor processing at "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo

# Exit gracefully
exit 0
