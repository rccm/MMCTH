#!/bin/bash

#-----------------------------------------------------------------------
# Script for MODIS Level-2 processing using direct broadcast data
#  Default directories are set within the env files located within the
#  $MODIS_L2_HOME/env sub-directory
#
#   Environmental variables that need to be set prior to execution of this
#    script:
#
#    MODIS_L2_HOME
#    LOCAL_ANC_DIR
#    REMOTE_ANC_DIR
#    MODIS_L2_BIN
#    MODIS_L2_CFG
#    PGSHOME
#    PGSMSG
#    PGSINC
#    PGSLIB
#    API_HOME
#    API_DIR
#    API_INC
#    API_SRC
#    HDF5LIB
#
# In addition, the PATH environmental variable must be set to include the 
# ${MODIS_L2_HOME}/ShellB3/bin, ${MODIS_L2_BIN}, ${MODIS_L2_HOME}/scripts, ${PGSHOME}, 
# and ${PGSMSG} directories.  The ${HDF5LIB} must be added to the LD_LIBRARY_PATH 
# environmental variable.  Again, these can be set by sourcing the env file located
# within the $MODIS_L2_HOME/env sub-directory.
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# GET ARGUMENTS
#-----------------------------------------------------------------------

run_cloudmask=1
run_sst=1
run_cloudtop=1
run_cloudopt=1
run_profiles=1
run_physretrieval=0 
run_aerosol=1
run_wvnir=1
run_ndsi=1
run_icec=1
run_ist=1
run_inv=1
run_quicklook=1
do_deltmp=0
do_reagg=1
do_comp=0

allargs="$*"

while [ 1 ]; do
  echo checking arg $1
  case $1 in
    "-nocloudmask")
      shift
      run_cloudmask=0
      ;;
    "-nosst")
      shift
      run_sst=0
      ;;
    "-nocloudtop")
      shift
      run_cloudtop=0
      ;;
    "-nocloudopt")
      shift
      run_cloudopt=0
      ;;
    "-noprofiles")
      shift
      run_profiles=0
      ;;
    "-nophysretrieval")
      shift
      run_physretrieval=0
      ;;
    "-noaerosol")
      shift
      run_aerosol=0
      ;;
    "-nowvnir")
      shift
      run_wvnir=0
      ;;
    "-nondsi")
      shift
      run_ndsi=0
      ;;
    "-noicec")
      shift
      run_icec=0
      ;;
    "-noist")
      shift
      run_ist=0
      ;;
    "-noinv")
      shift
      run_inv=0
      ;;
    "-noquicklook")
      shift
      run_quicklook=0
      ;;
    "-deltmp")
      shift
      do_deltmp=1
      ;;
    "-noreagg")
      shift
      do_reagg=0
      ;;
    "-compress")
      shift
      do_comp=1
      ;;
    *)
      break
  esac
done

echo $*

# Check arguments 
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
  echo "--------------------------------------------------------------------"
  echo "Usage: modis_level2.sh [Options] L1BFILE OUTDIR [TMPDIR]"
  echo ""
  echo "Where:"
  echo ""
  echo " [Options] - These flags are optional. Default is to run all IMAPP"
  echo "             product executables, and for the temporary run directory"
  echo "             to remain upon completion."
  echo "   -nocloudmask     - do not run cloudmask"
  echo "   -nosst           - do not run sst"
  echo "   -nocloudtop      - do not run cloudtop properties"
  echo "   -nocloudopt      - do not run cloud optical properties"
  echo "   -noprofiles      - do not run profiles"
  echo "   -noaerosol       - do not run aerosols"
  echo "   -nowvnir         - do not run wvnir"
  echo "   -noist           - do not run ice surface temperature"
  echo "   -nondsi          - do not run snow mask"
  echo "   -noicec          - do not run ice concentration"
  echo "   -noinv           - do not run temperature inversion software"
  echo "   -noreagg         - do not re-aggregate L1B data before creating products"
  echo "   -noquicklook     - do not run quicklook plotting script"
  echo "   -compress        - Internally compress output product files (uses hrepack)"
  echo "   -deltmp          - delete temporary directory when complete"
  echo ""
  echo " L1BFILE - Level-1B 1KM file"
  echo " OUTDIR  - output file directory"
  echo " TMPDIR  - [Optional]  Name of temporary directory to create"
  echo "             or re-use for flat files. Default is date-based name."
  echo ""
  exit 1
  echo "--------------------------------------------------------------------"
fi
  echo ""

# Get arguments
FIL1KM=$1
OUTDIR=$2
if [ $# == 3 ]; then
  TMPDIR=$3
fi

# Set the Latitude Thresholds for Polar Products to be created
# If the center latitude of the pass is North or South of 
# these thresholds, the polar products will be created if the
# user choose to produce those products.
NH_THRESHOLD=40 #60
SH_THRESHOLD=-40 #-60

# Check for Level-1B 1KM file
if [ ! -e $FIL1KM ]; then
  echo "Input Level-1B 1KM file does not exist: "$FIL1KM
  exit 1
fi

# Get satellite name given input BASE_FIL1KM form imapp naming convention (a1 and t1)
# or NASA naming convention (MYD and MOD)
BASE_FIL1KM=`basename $FIL1KM`
F1=`echo $BASE_FIL1KM:t | cut -c1`
F3=`echo $BASE_FIL1KM:t | cut -c1-3`
if [ $F1 == "t" ] || [ $F3 == "MOD" ]; then
  SAT="Terra"
  SATMC="TERRA"
  LSAT="terra"
  PREFIX="t1"
elif [ "$F1" == "a" ] || [ "$F3" == "MYD" ]; then
  SAT="Aqua"
  SATMC="AQUA"
  LSAT="aqua"
  PREFIX="a1"
else
  echo "Unable to determine satellite name from file name" $BASE_FIL1KM
  exit 1
fi

# Get directory paths
DIR_NAME=`dirname $FIL1KM`
STARTDIR=$PWD

if [ "$OUTDIR" == "." ] || [ "$OUTDIR" == "./" ]; then
  OUTDIR=$PWD
fi
if [ "$DIR_NAME" == "." ]; then
  DIR_NAME=$PWD
  FIL1KM=$DIR_NAME/$FIL1KM
fi

# Print start message
echo "#################################################################"
echo "MODIS Level-2 processing started at "`date -u`
echo "#################################################################"
echo " "
echo " Input Arguments: "
echo "    Input L1B Filename: "$FIL1KM
echo "    Output File Directory Name: "$OUTDIR
echo " "

# Get date and time string from L1B 1KM meter file
DATE_TIME=`echo $BASE_FIL1KM:t | cut -d. -f2,3 | cut -dA -f2`
DATE=`echo $DATE_TIME | cut -d. -f1 | cut -dA -f2`
TEMP_TIME=`echo $DATE_TIME | cut -d. -f2`
TIME=`echo $TEMP_TIME | awk '{print $1 + 0}'`

DATE_ORIG=$DATE 
if [ ${#DATE} -eq 5 ]; then
  DATE=20$DATE
fi

# Determine if Day or Night data set defined as Day scans > 0
DAY_SCANS=`ncdump -h ${FIL1KM} | grep "Number_of_Day_mode_scans" | awk  '{ print $3 }'`
if [ "${DAY_SCANS}" == "" ]; then
  DAY_SCANS=`ncdump -h ${FIL1KM} | grep "Number of Day mode scans" | awk  '{ print $7 }'`
fi
if [ ${DAY_SCANS} -gt 100 ]; then
  LIGHT="DAY"
else
  LIGHT="NIGHT"
fi
echo " "
echo " MODIS granule is a $LIGHT Granule"
echo " "

# Get number of 1 km lines of data
NSCANS=`ncdump -h ${FIL1KM} | grep "Number_of_Scans" | awk  '{ print $3 }'`
if ["${NSCANS}" == ""]; then
  NSCANS=`ncdump -h ${FIL1KM} | grep "Number of Scans" | awk  '{ print $5 }'`
fi
echo " Number of Scans in File " $NSCANS
NLINS=`expr 10 "*" $NSCANS`
echo " Number of Lines in File " $NLINS
echo " "

# Get the center Latitude of the data segment for determining the
#  projection of the HDF Quicklooks

NLAT=`ncdump -h ${FIL1KM} | awk '/NORTHBOUNDINGCOORDINATE/,/END_OBJECT/ {print $4}' | grep "\." |cut -d. -f1`
echo "Extracted Northernmost Latitude: " $NLAT
echo " "

SLAT=`ncdump -h ${FIL1KM} | awk '/SOUTHBOUNDINGCOORDINATE/,/END_OBJECT/ {print $4}' | grep "\." |cut -d. -f1`
echo "Extracted Southernmost Latitude: " $SLAT
echo " "

if [ $NLAT -gt 0 ]; then
   HEMI=1
elif [ $NLAT -lt 0 ]  &&  [ $NLAT -gt -90 ]; then
   HEMI=0
elif [ $NLAT -eq 0 ] && [ $SLAT -lt 0 ]; then
   HEMI=0
else
   echo "ERROR:  Could not determine data hemisphere " $NLAT
   exit 1
fi

# Get the center Latitude of the data segment for determining the
#  projection of the HDF Quicklooks
CENTER_LAT=`expr $NLAT + $SLAT`
CENTER_LAT=`expr $CENTER_LAT / 2`

echo "Center Latitude of Data Segment :" $CENTER_LAT
echo " "

POLAR="FALSE"
if [ $CENTER_LAT -gt $NH_THRESHOLD ] || [ $CENTER_LAT -lt $SH_THRESHOLD ];  then
   POLAR="TRUE"
   suffix="_polar"
else
   suffix=""
fi

#-----------------------------------------------------------------------
# CREATE TEMPORARY DIRECTORY
#-----------------------------------------------------------------------

# default is to use a date stamped directory name (YYYY_MM_DD_sat_HHMM_hhmmss)
if [ ! $TMPDIR ]; then
  TMPDIR=${DATE}_${SAT}_${TIME}_`date -u +%H%M%S`
fi

# Check for a temporary directory, otherwise create one
if [ -d "$TMPDIR" ]; then
  run_extract=0
else
  mkdir $TMPDIR
  if [ $? -ne 0 ]; then
    echo "Could not create temporary directory: "$TMPDIR
    exit 1
  fi
  run_extract=1
fi

# Enter temporary directory
cd $TMPDIR

#-----------------------------------------------------------------------
# GET MODIS LEVEL-1 GEO AND QKM FILENAMES
#-----------------------------------------------------------------------

# Set file identifier string (e.g., A2004089.1633)
FILEID=`echo $BASE_FIL1KM | cut -d. -f2,3`
FILEID2=`echo $FILEID | cut -dA -f2`

if [ $F1 == "t" ]; then
  PREFIX="t1"
elif [ $F1 == "a" ]; then
  PREFIX="a1"
elif [ $F3 == "MOD" ] || [ $F3 == "MYD" ] ; then
  PREFIX=$F3
fi

# Set geolocation filename
if [ $PREFIX == "MOD" ] || [ $PREFIX == "MYD" ]; then
  FILGEO=`find $DIR_NAME -mindepth 1 -maxdepth 1 -name "${PREFIX}03.${FILEID}*.hdf"`
elif [ $PREFIX == "t1" ] || [ $PREFIX == "a1" ]; then
  FILGEO=`find $DIR_NAME -mindepth 1 -maxdepth 1 -name "${PREFIX}.${FILEID}.geo.hdf"`
fi

if [ "$FILGEO" == '' ]; then
  echo "Input geolocation file does not exist"
  exit 1
fi

# Set QKM filename
if [ $PREFIX == "MOD" ] || [ $PREFIX == "MYD" ]; then
  FILQKM=`find $DIR_NAME  -name "${PREFIX}02QKM.${FILEID}*.hdf" | tail -1` 
elif [ $PREFIX == "t1" ] || [ $PREFIX == "a1" ]; then
  FILQKM=`find $DIR_NAME  -name "${PREFIX}.${FILEID}.250m.hdf" | tail -1` 
fi

if [ "$FILQKM" == '' ]; then 
  FILQKM='MISSING'
fi

# Set HKM filename
if [ $PREFIX == "MOD" ] || [ $PREFIX == "MYD" ]; then
  FILHKM=`find $DIR_NAME  -name "${PREFIX}02HKM.${FILEID}*.hdf" | tail -1`
elif [ $PREFIX == "t1" ] || [ $PREFIX == "a1" ]; then
  FILHKM=`find $DIR_NAME  -name "${PREFIX}.${FILEID}.500m.hdf" | tail -1`
fi

if [ "$FILHKM" == '' ]; then
  FILHKM='MISSING'
fi

# Print MODIS Level-1 filenames for log
echo "   INPUT FILES  "
echo "MOD021KM file: "$FIL1KM
echo "MOD02HKM file: "$FILHKM
echo "MOD02QKM file: "$FILQKM
echo "MOD03 file: "$FILGEO
echo "Temporary Directory: "$TMPDIR

#-----------------------------------------------------------------------
# DESTRIPE DATA IF NEEDED
#-----------------------------------------------------------------------
# Uncompress the data
unhrepack.sh $FIL1KM
# Check if file has already been destriped before running destriping code
destripe=`ncdump -h $FIL1KM | grep -i destripe`
if [ $? == 1 ]; then
  run_modis_destripe.sh $SAT $FIL1KM
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "L1B file already destriped - destriping code not run"
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo
fi

#-----------------------------------------------------------------------
# IF USING REAGGREGATED FILES, COPY 1KM L1B FILE TO RUN DIRECTORY PRIOR
# TO REAGGREGATION SO ORIGINAL DATA IN 1KM L1B FILE NOT OVERWRITTEN
#-----------------------------------------------------------------------

if [ $SAT == "Terra" ]; then
  do_reagg=0
fi

if [ $do_reagg == 1 ]; then
  cp $FIL1KM .
  FIL1KM=`pwd`/`basename $FIL1KM`
fi

#-----------------------------------------------------------------------
# FIRST RUN THE FLATFILE EXTRACTORS 
#-----------------------------------------------------------------------

if [ $run_sst == 1 ] || [ $run_wvnir == 1 ]; then
  if [ $run_extract == 1 ]; then

    if [ -e ${MODIS_L2_BIN}/modis_extract_1km.exe ]; then
       if [ $do_reagg == 1 ];then
         run_modis_flatfile.sh $SAT `basename $FIL1KM` $FILQKM $FILGEO
       else
         run_modis_flatfile.sh $SAT $FIL1KM $FILQKM $FILGEO
       fi
    else
      # You need the flat binary extracted files in order to the run the 
      # rest of the software.
      echo "Flatfile extractors did not execute correctly."
      exit -1
    fi

  else
    echo "Using previously created flat files from $TMPDIR"
  fi
fi

#-----------------------------------------------------------------------
# RUN THE CLOUD MASK (MOD35)
#-----------------------------------------------------------------------

if [ $run_cloudmask == 1 ]; then
  # Run the cloudmask if the executable exists
  if [ -e ${MODIS_L2_BIN}/cloudmask.exe ]; then
    run_modis_cloudmask.sh $SAT $FIL1KM $FILQKM $FILGEO gdas1 gdas2 rsst nise icec $OUTDIR $do_reagg $do_comp
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Cloudmask disabled - no executable"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Cloudmask disabled - by option"    
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
fi

#-----------------------------------------------------------------------
# RUN THE IMAPP SST ALGLORITHM
#-----------------------------------------------------------------------

if [ $run_sst == 1 ]; then
  # Run the sst software if the executable exists
  if [ -e ${MODIS_L2_BIN}/sst.exe ]; then
     run_modis_sst.sh $SAT $OUTPUT_TYPE $FIL1KM $OUTDIR $do_comp
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "SST disabled - no executable"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "SST disabled - by option"    
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
fi

#-----------------------------------------------------------------------
# RUN THE IMAPP CLOUD TOP PROPERTIES AND CLOUD PHASE ALGLORITHM (MOD06CT)
#-----------------------------------------------------------------------

# Get cloudmask name - how should I deal with the cloudmask file if it's created by the MOD_PR35 script or
# if it's given by the user

FILCM=`find . -name "${PREFIX}35_L2.*.hdf"`
if [ ! $FILCM ]; then
  FILCM=`find . -name "${PREFIX}.*.mod35.hdf"`
fi

if [ ! $FILCM ]; then
  echo "WARNING: Cloudmask file does not exist"
  exit 1
fi


if [ $run_cloudtop == 1 ]; then
  # Run the ctp software if the executables exist
  if [ -e ${MODIS_L2_BIN}/createcloudfile.exe ] && [ -e ${MODIS_L2_BIN}/cloudtop.exe ] && [ -e ${MODIS_L2_BIN}/cirruscloud.exe ]; then
    run_modis_cloudtop.sh $SAT $FIL1KM $FILGEO $FILCM filcsrb gdas1 gdas2 rsst nise icec $OUTDIR $do_reagg $do_comp
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Cloud top properties disabled - no executable"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Cloud top properties disabled - by option"    
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo
fi

#-----------------------------------------------------------------------
# RUN THE IMAPP ATMOSPHERIC PROFILES STATISTICAL RETRIEVAL ALGORITHM (MOD07)
#-----------------------------------------------------------------------

if [ $run_profiles == 1 ]; then
  # Run the profiles algorithm if the executable exists
  if [ -e ${MODIS_L2_BIN}/profiles.exe ]; then
    run_modis_profiles.sh $SAT $FIL1KM $FILGEO $FILCM gdas $OUTDIR $do_reagg $do_comp
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Profiles processing disabled - no executable"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Profiles processing disabled - by option"    
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo 
fi

#-----------------------------------------------------------------------
# RUN THE IMAPP ICE SURFACE TEMPERATURE ALGLORITHM - Jeff Key, NOAA
#-----------------------------------------------------------------------

if [ $POLAR == "TRUE" ]; then
  if [ $run_ist == 1 ]; then
    # Run the ice surface temperature software if the executable exists
    if [ -e ${MODIS_L2_BIN}/ice_surface_temperature.exe ]; then
      run_modis_ist.bash $HEMI $FIL1KM $FILCM $OUTDIR $OUTPUT_TYPE $do_comp
    else
      echo
      echo
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo "Ice surface temperature processing disabled - no executable"
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo
      echo
    fi
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Ice surface temperature processing disabled - by option"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Ice surface temperature not run - granule did not meet polar"
  echo "requirements or ice surface temperature software disabled." 
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo
fi

#-----------------------------------------------------------------------
# RUN THE IMAPP ICE MASK AND ICE CONCENTRATION ALGLORITHM - Yinghui Liu,
#  UW-Madison, CIMSS/SSEC and Jeffry Key, NOAA/NESDIS/ASPB
#-----------------------------------------------------------------------

if [ $POLAR == "TRUE" ]; then
  if [ $run_icec == 1 ]; then
    # Run the ice mask and ice concentration algorithm if it exists
    if [ -e ${MODIS_L2_BIN}/ice_concentration.exe ]; then
      run_modis_ice_concentration.bash $FIL1KM $FILCM $FILGEO $OUTDIR $OUTPUT_TYPE $do_comp
    else
      echo
      echo
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo "Ice concentration processing disabled - no executable"
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo
      echo
    fi
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Ice concentration processing disabled - by option"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Ice concentration not run - granule did not meet polar"
  echo "requirements or ice concentration software disabled." 
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo
fi

#-----------------------------------------------------------------------
# RUN THE MODIS IMAPP TEMPERATURE INVERSION RETRIEVAL ALGLORITHM - 
#  Yinghui Liu, William Straka III UW-Madison, CIMSS/SSEC and 
#  Jeffry Key, NOAA/NESDIS/ASPB
#-----------------------------------------------------------------------

if [ $POLAR == "TRUE" ]; then
  if [ $run_inv == 1 ]; then
    # Run the temperature inversion algorithm if it exists
    if [ -e ${MODIS_L2_BIN}/inversions.exe ]; then
      run_modis_inversions.bash $HEMI $FIL1KM $FILCM $OUTDIR $OUTPUT_TYPE $do_comp
    else
      echo
      echo
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo "Temperature inversions processing disabled - no executable"
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo
      echo
    fi
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Temperature inversions processing disabled - by option"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Polar inversion software not run - granule did not meet polar"
  echo "requirements or polar inversion software disabled." 
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo
fi

#-----------------------------------------------------------------------
# RUN THE IMAPP DAYTIME ONLY PROCESSES
#-----------------------------------------------------------------------
    
if [ ${LIGHT} == "DAY" ]; then

#-----------------------------------------------------------------------
# RUN THE MODIS CLOUD OPTICAL PROPERTIES ALGLORITHM (MOD06OD) 
#-----------------------------------------------------------------------
 
if [ $run_cloudopt == 1 ]; then
 
  FILMOD06CT=`find . -name "${PREFIX}06CT_L2.*.hdf"`
  if [ ! $FILMOD06CT ]; then
    if [ $PREFIX == "MYD" ] || [ $PREFIX == "a1" ]; then
      FILMOD06CT=`find . -name "a1.*.mod06ct.hdf"`
    else
      FILMOD06CT=`find . -name "t1.*.mod06ct.hdf"`
    fi
    if [ ! $FILMOD06CT ]; then
      FILMOD06CT='MISSING'
    fi
  fi

#   if [ $run_cloudopt == 1 ]; then
   # Run the cloud optical properties software if the executable exists
      if [ -e ${MODIS_L2_BIN}/cloudoptical.exe ]; then
        # You must have unpacked the MOD06OD coefficients file
        if [ -e ${MODIS_L2_HOME}/src/cloudoptical/coeff/staticdata/Forward_Libraries/MODIS_Ice_library.hdf ]; then
          if [ $FILMOD06CT != ''MISSING'' ]; then
            # Copy the MOD06CT File to MOD06 for use by the cloud optical properties
            #  software
            FILMOD06=${PREFIX}.${FILEID}.mod06.hdf
            cp $FILMOD06CT $FILMOD06
            run_modis_cloudoptical.sh $SAT $FIL1KM $FILCM $FILGEO $FILMOD06 $OUTDIR
          else
            echo "Cloud optical properties not run - no cloud top products file found"
          fi
        else
          echo "Cloud optical properties coefficient files not found"
        fi
      else
        echo "Cloud optical properties disabled - no executable"
      fi
else
  echo "Cloud optical properties disabled - by option"
fi

#-----------------------------------------------------------------------
# RUN THE IMAPP AEROSOL ALGORITHM (MOD04)
#-----------------------------------------------------------------------

if [ $run_aerosol == 1 ]; then
  FILPROF=`find . -name "${PREFIX}07_L2.*.hdf"`
  if [ ! $FILPROF ]; then
    if [ $PREFIX == "MYD" ] || [ $PREFIX == "a1" ]; then
      FILPROF=`find . -name "a1.*.mod07.hdf"`
    else
      FILPROF=`find . -name "t1.*.mod07.hdf"`
    fi
  fi

  # Run the aerosol algorithm if the executable exists
  if [ -e ${MODIS_L2_BIN}/aerosol.exe ]; then
    if [ -e $FILHKM ] && [ -e $FILQKM ]; then
      if [ "$FILPROF" != "" ] && [ -e $FILPROF ]; then
        run_modis_aerosol.sh $SAT $FIL1KM $FILHKM $FILQKM $FILGEO $FILPROF $FILCM $OUTDIR $do_reagg $do_comp
      else
        echo
        echo
        echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        echo "Aerosol processing disabled - MOD07 profiles file does not exist"
        echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        echo
        echo
      fi
    else
      echo
      echo
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo "Aerosol processing disabled - 250m and/or 500m L1B file"
      echo " not available"
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo
      echo
    fi
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Aerosol processing disabled - no executable"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Aerosol processing disabled - by option"
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo
fi

#-----------------------------------------------------------------------
#   RUN THE PETER ALBERT WATER VAPOR RETRIEVAL ALGORITHM
#-----------------------------------------------------------------------
        
if [ $run_wvnir == 1 ]; then
  if [ -e ${MODIS_L2_BIN}/wvnir.exe ]; then
    run_modis_wvnir.sh $SAT $OUTPUT_TYPE $FIL1KM $OUTDIR $do_comp
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "WVNIR processing disabled - no executable"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "WVNIR processing disabled - by option"
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo
fi

#-----------------------------------------------------------------------
# RUN THE IMAPP NDSI SNOW MASK ALGORITHM - William Straka III, CIMSS
#-----------------------------------------------------------------------

if [ $POLAR == "TRUE" ]; then
  if [ $run_ndsi == 1 ]; then
    # Run the snow mask software if the executable exists
    if [ -e ${MODIS_L2_BIN}/snow_mask.exe ]; then
      run_modis_snowmask.bash $SAT $FIL1KM $FILCM $OUTDIR $OUTPUT_TYPE $do_comp
    else
      echo
      echo
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo "Snow mask processing disabled - no executable"
      echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      echo
      echo
    fi
  else
    echo
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "Snow mask processing disabled - by option"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    echo
  fi
else
  echo
  echo
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo "Snow mask software not run - granule did not meet polar"
  echo "requirements or snow mask software disabled." 
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo
  echo
fi

#-----------------------------------------------------------------------
# END THE DAYTIME ONLY PROCESSING
#-----------------------------------------------------------------------

fi

#-----------------------------------------------------------------------
# CREATE QUICKLOOK IMAGES
#-----------------------------------------------------------------------

if [ $run_quicklook == 1 ]; then

  run_quicklook.sh $SAT $DATE_TIME $FILGEO $FIL1KM $OUTDIR $POLAR 

fi

#-----------------------------------------------------------------------
# CLEANUP AND EXIT
#-----------------------------------------------------------------------

# Remove pc files from tookit runtime directory
rm -rf $PGSHOME/runtime/pc*

# Remove cloud mask file from run directory
#rm -rf $TMPDIR/$FILCM
#rm -rf $TMPDIR/temp_flist
#rm -rf $TMPDIR/*.hdf

# Remove the temporary directory if option is set
if [ $do_deltmp == 1 ]; then
  cd $STARTDIR
  rm -rf $TMPDIR
fi   
      
# Print finish message
echo
echo "#################################################################"        
echo "MODIS Level-2 processing finished at "`date -u`
echo "#################################################################"
    
exit 0
