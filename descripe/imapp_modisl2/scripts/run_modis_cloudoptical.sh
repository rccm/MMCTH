#!/bin/bash -x

# Run the IMAPP MODIS atmopheric profiles (vertical profiles of
# temperature and moisture plus total precipitable water vapor
# and stability indices.
#
# You must first set up the environmental variables to match your system and
# directory structure prior to execution.  These variables can be 
# automatically set by editing the script found in the 
# imapp_modisl2/env directory.  Set the top level
# directory for the modis level 2 processing path and then type:
#  source imapp_modisl2.bash_env for bash shells
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
#  In addtion, sourcing the bash environmental script will set the path
#  so that the binaries can be executed from any directory.
# PATH=.:${MODIS_L2_BIN}:${MODIS_L2_HOME}/scripts:${PATH} 
#
# The option to save both the binary and/or hdf files in now included
# with this release. 
#
#-------------------------------------------------------------------------------
# SETUP
#-------------------------------------------------------------------------------

# Set coefficient directories
COEFFDIR=${MODIS_L2_HOME}/src/cloudoptical/coeff/staticdata
CTCOEFFDIR=${MODIS_L2_HOME}/src/cloudtop/MOD_PR06CT/coeff/lit

# Check arguments
if [ $# != 6 ]; then
  echo "Usage: run_cloud_optical.sh SAT FIL1KM FILMASK FILGEO FILEMOD06 OUTDIR"
  echo "where"
  echo "SAT is the satellite platform (aqua or terra)"
  echo "FIL1KM is MODIS L1B 1000 meter resolution HDF image file"
  echo "FILMASK is the MODIS L1B 1000 meter Cloud Mask HDF file"
  echo "FILGEO is MODIS L1B 1000 meter resolution geolocation HDF file"
  echo "OUTDIR is the directory where products will be stored"
  exit 1 
fi

# Extract file names
SAT=$1
FIL1KM=$2
FILMASK=$3
FILGEO=$4
FILMOD06=$5
OUTDIR=$6

# Get platform header (MOD or MYD)
if [ $SAT == "terra" ] || [ $SAT == "Terra" ] || [ $SAT == "TERRA" ]; then
  PREFIX="MOD"
elif [ $SAT == "aqua" ] || [ $SAT == "Aqua" ] || [ $SAT == "AQUA" ]; then
  PREFIX="MYD"
else
  echo "Satellite name not recognized: "$SAT
  exit 1
fi

# Set root of output file name (e.g. 't1.02001.1815')
ROOT=`basename $FIL1KM`
ROOT=`echo $ROOT | cut -d. -f1-3`

# Get date and time string from L1B 1KM meter file
DATE_TIME=`echo $ROOT | cut -d. -f2,3 | cut -dA -f2`
DATE=`echo $DATE_TIME | cut -d. -f1 | cut -dA -f2`
if [ ${#DATE} -eq 5 ]; then
  DATE=20$DATE
fi
YEAR=`echo $DATE | cut -c 1-4`
DAY=`echo $DATE | cut -c 5-7`

TIME=`echo $DATE_TIME | cut -d. -f2`
HOUR=`echo $TIME | cut -c 1-2`
MINUTE=`echo $TIME | cut -c 3-4`

INPUTS="${YEAR} ${DAY} ${TIME}"

#-------------------------------------------------------------------------------
# CHECK FOR CLOUD OPTICAL PROPERTIES COEFFICIENT DIRECTORY
#-------------------------------------------------------------------------------

if [ ! -e $MODIS_L2_HOME/src/cloudoptical/coeff ]; then
  echo
  echo "WARNING: No MOD06 cloud optical properties coefficient directory found." 
  echo " MOD06 cloud optical properties processing cannot be run."
  exit 1
fi

#-------------------------------------------------------------------------------
# SET UP INPUT FILES
#-------------------------------------------------------------------------------

# Get new filnames
DAAC_FIL1KM=${PREFIX}021KM.A${YEAR}${DAY}.${HOUR}${MINUTE}.hdf
DAAC_FILMASK=${PREFIX}35_L2.A${YEAR}${DAY}.${HOUR}${MINUTE}.hdf
DAAC_FILGEO=${PREFIX}03.A${YEAR}${DAY}.${HOUR}${MINUTE}.hdf
DAAC_FILMOD06=${PREFIX}06_L2.A${YEAR}${DAY}.${HOUR}${MINUTE}.hdf

# Create output directiry within temporary run directory
OUTPUTDIRNAME=./outputs-${YEAR}${DAY}.${TIME}
mkdir $OUTPUTDIRNAME
mkdir $OUTPUTDIRNAME/OD

# Create input directory within temporary run directory
INPUTDIRNAME=./inputs-${YEAR}${DAY}.${TIME}
mkdir $INPUTDIRNAME
cd $INPUTDIRNAME

#Link DB file naming convention to MYD and MOD standard DAAC naming convention
ln -s $FIL1KM $DAAC_FIL1KM
ln -s ../$FILMASK $DAAC_FILMASK
ln -s $FILGEO $DAAC_FILGEO
#ln -s ../$FILMOD06 ../$OUTPUTDIRNAME/OD/$DAAC_FILMOD06
cp ../$FILMOD06 ../$OUTPUTDIRNAME/OD/$DAAC_FILMOD06

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Getting ancillary data for MODIS cloud properties processing"

# Clear sky reflectance file
# ----------------------------------------------------------------
CSR=`get_anc_clearsky.sh $SAT $DATE`
echo "CSR = "$CSR
# If you can't find a file, set it equal to missing
status=$?
if [ $status != 0 ]; then
  CSR="none"
fi
echo "CSR = "$CSR

# NCEP ice concentration file
# ----------------------------------------------------------------
ICEC=`get_anc_icec.sh $DATE`
# If you can't find a file, set it equal to missing
status=$?
if [ $status != 0 ]; then 
  ICEC="MISSING"
fi
# Link file to this directory
GRIBICE=`basename $ICEC`
ln -s ${ICEC} .
# Convert backslashes so that they will be properly inserted using sed
ICEC=`echo ${PWD}/${GRIBICE} | sed s/"\/"/"\\\/"/g`
ICEC=`basename $ICEC`

# NISE snow and ice files
# ----------------------------------------------------------------
NISE=`get_anc_nise.sh $DATE`
# If you can't find a file, set it equal to missing
status=$?
if [ $status != 0 ]; then
  NISE="MISSING"
fi
ln -s ${NISE} .
# Convert backslashes so that they will be properly inserted using sed
HDFNISE=`echo ${NISE} | sed s/"\/"/"\\\/"/g`
HDFNISE=`basename $HDFNISE`

# GDAS1 numerical weather prediction model grid field analysis nearest before
# ---------------------------------------------------------------------------
gdas_datetime_before=( $(get_anc_gdas_gfs_before_and_after_time.sh $DATE $TIME BACKWARD) )
DATE_BEFORE=`echo $gdas_datetime_before | cut -c1-7`
TIME_BEFORE=`echo $gdas_datetime_before | cut -c8-11`

GDASbefore=`get_anc_gdas_gfs.sh $DATE_BEFORE $TIME_BEFORE`
status=$?
# If you can't find a file, then try an older gdas file.
if [ "$status" != "0" ]; then
  GDASbefore=`get_nearest_gdas.sh $DATE_BEFORE $TIME_BEFORE`
  status=$?
  if [ "$status" != "0" ]; then
    GDASbefore="MISSING"
  fi
fi  

ln -s ${GDASbefore} .
GDASbefore=`basename ${GDASbefore}`
    
echo "Nearest GDAS file before is " $GDASbefore
#
## GDAS1 numerical weather prediction model grid field analysis nearest after
## --------------------------------------------------------------------------
gdas_datetime_after=( $(get_anc_gdas_gfs_before_and_after_time.sh $DATE $TIME FORWARD) )
DATE_AFTER=`echo $gdas_datetime_after | cut -c1-7`
TIME_AFTER=`echo $gdas_datetime_after | cut -c8-11`

GDASafter=`get_anc_gdas_gfs.sh $DATE_AFTER $TIME_AFTER`
status=$?
# If you can't find a file, then try an older gdas file.
if [ "$status" != "0" ]; then
  GDASafter=`get_nearest_gdas.sh $DATE_AFTER $TIME_AFTER`
  status=$?
  if [ "$status" != "0" ]; then
    GDASafter="MISSING" 
  fi
fi

ln -s ${GDASafter} .
GDASafter=`basename ${GDASafter}`

echo "Nearest GDAS file after is " $GDASafter



# Prepare the coeffcient file names for insertion into config file
# ----------------------------------------------------------------
#IGBP base map
E1=`find ${COEFFDIR}/Surface_Albedo/IGBPmap -follow -name "IGBP.EcoMap.NtoS.2004.149.v004.hdf" \
  -print` 
status=$?
if [ $status == 0 ]; then
  ECOM=`echo ${E1} | sed s/"\/"/"\\\/"/g`
fi
# Snow Albedo Map
S1=`find ${COEFFDIR}/Surface_Albedo/SnowAlbedoStats -follow -name "AlbSnwSts.ByNISE.W90.D90.WS.Hemi.2000-2004.YrAvg.hdf" -print` 
status=$?
if [ $status == 0 ]; then
  ALBSL=`echo ${S1} | sed s/"\/"/"\\\/"/g`
fi
# Transmittance Map
T1=`find ${COEFFDIR}/Water_Vapor_Transmittance -follow -name "Transmittance.hdf" -print` 
status=$?
if [ $status == 0 ]; then
  TRANW=`echo ${T1} | sed s/"\/"/"\\\/"/g`
fi
# Forward Model Library for Ice
I1=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Ice_library.hdf" -print` 
status=$?
if [ $status == 0 ]; then
  ICEL=`echo ${I1} | sed s/"\/"/"\\\/"/g`
fi

I1=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Ice_library_ws3.hdf" -print`
status=$?
if [ $status == 0 ]; then
  ICELWS3=`echo ${I1} | sed s/"\/"/"\\\/"/g`
fi

I1=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Ice_library_ws3sd.hdf" -print`
status=$?
if [ $status == 0 ]; then
  ICELWS3SD=`echo ${I1} | sed s/"\/"/"\\\/"/g`
fi

I1=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Ice_library_ws7.hdf" -print`
status=$?
if [ $status == 0 ]; then
  ICELWS7=`echo ${I1} | sed s/"\/"/"\\\/"/g`
fi

I1=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Ice_library_ws7sd.hdf" -print`
status=$?
if [ $status == 0 ]; then
  ICELWS7SD=`echo ${I1} | sed s/"\/"/"\\\/"/g`
fi

I1=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Ice_library_ws15.hdf" -print`
status=$?
if [ $status == 0 ]; then
  ICELWS15=`echo ${I1} | sed s/"\/"/"\\\/"/g`
fi

I1=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Ice_library_ws15sd.hdf" -print`
status=$?
if [ $status == 0 ]; then
  ICELWS15SD=`echo ${I1} | sed s/"\/"/"\\\/"/g`
fi

# Forward Model Library for Water
W2=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Water_library.hdf" -print`
status=$?
if [ $status == 0 ]; then
  WATERL=`echo ${W2} | sed s/"\/"/"\\\/"/g`
fi

W2=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Water_library_ws3.hdf" -print`
status=$?
if [ $status == 0 ]; then
  WATERLWS3=`echo ${W2} | sed s/"\/"/"\\\/"/g`
fi

W2=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Water_library_ws3sd.hdf" -print`
status=$?
if [ $status == 0 ]; then
  WATERLWS3SD=`echo ${W2} | sed s/"\/"/"\\\/"/g`
fi

W2=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Water_library_ws7.hdf" -print`
status=$?
if [ $status == 0 ]; then
  WATERLWS7=`echo ${W2} | sed s/"\/"/"\\\/"/g`
fi

W2=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Water_library_ws7sd.hdf" -print`
status=$?
if [ $status == 0 ]; then
  WATERLWS7SD=`echo ${W2} | sed s/"\/"/"\\\/"/g`
fi

W2=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Water_library_ws15.hdf" -print`
status=$?
if [ $status == 0 ]; then
  WATERLWS15=`echo ${W2} | sed s/"\/"/"\\\/"/g`
fi

W2=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Water_library_ws15sd.hdf" -print`
status=$?
if [ $status == 0 ]; then
  WATERLWS15SD=`echo ${W2} | sed s/"\/"/"\\\/"/g`
fi

# Forward Model Library for Ice Water Phase Function
W2=`find ${COEFFDIR}/Forward_Libraries -follow -name "MODIS_Ice_WaterPhaseFunc.hdf" -print`
status=$?
if [ $status == 0 ]; then
  ICEWATERPHASE=`echo ${W2} | sed s/"\/"/"\\\/"/g`
fi


#Locate Surface Albedo files for 5 wavelengths
# ----------------------------------------------------------------

ALBDATE=`get_surface_albedo_date_monthly.sh $DATE`
if [ $ALBDATE == "" ]; then
  $ALBDATE = "MISSING"
fi

A1=`find ${COEFFDIR}/Surface_Albedo/AlbedoMaps/2013/monthlyavg/00-05.${ALBDATE} -follow -name "MCD43GF_wsa_Band1_${ALBDATE}_2013_avg.hdf" -print`
echo "A1 = "$A1
ATEMP1=$A1
status=$?
if [ $status == 0 ]; then
  ALB659=`echo $ATEMP1 | sed s/"\/"/"\\\/"/g`
fi
A2=`find ${COEFFDIR}/Surface_Albedo/AlbedoMaps/2013/monthlyavg/00-05.${ALBDATE} -follow -name "MCD43GF_wsa_Band2_${ALBDATE}_2013_avg.hdf" -print`
echo "A2 = "$A2
ATEMP2=$A2
status=$?
if [ $status == 0 ]; then
  ALB858=`echo $ATEMP2 | sed s/"\/"/"\\\/"/g`
fi
A3=`find ${COEFFDIR}/Surface_Albedo/AlbedoMaps/2013/monthlyavg/00-05.${ALBDATE} -follow -name "MCD43GF_wsa_Band5_${ALBDATE}_2013_avg.hdf" -print`
echo "A3 = "$A3
ATEMP3=$A3
status=$?
if [ $status == 0 ]; then
  ALB124=`echo $ATEMP3 | sed s/"\/"/"\\\/"/g`
fi
A4=`find ${COEFFDIR}/Surface_Albedo/AlbedoMaps/2013/monthlyavg/00-05.${ALBDATE} -follow -name "MCD43GF_wsa_Band6_${ALBDATE}_2013_avg.hdf" -print`
echo "A4 = "$A4
ATEMP4=$A4 
status=$?
if [ $status == 0 ]; then
  ALB164=`echo $ATEMP4 | sed s/"\/"/"\\\/"/g`
fi
A5=`find ${COEFFDIR}/Surface_Albedo/AlbedoMaps/2013/monthlyavg/00-05.${ALBDATE} -follow -name "MCD43GF_wsa_Band7_${ALBDATE}_2013_avg.hdf" -print`
echo "A5 = "$A5
ATEMP5=$A5 
status=$?
if [ $status == 0 ]; then
  ALB213=`echo $ATEMP5 | sed s/"\/"/"\\\/"/g`
fi

cd ..

# Set name of emissivity file
EMISS=${CTCOEFFDIR}/global_emiss_intABI_2005001.hdf

#-------------------------------------------------------------------------------
# SET UP PCF FILE
#-------------------------------------------------------------------------------

# Create a configuration file from the template
TEMPLATE=${MODIS_L2_CFG}/cloud_optical.cfg
sed \
  -e "s?TRACK_DATE?${INPUTS}?g" \
  -e "s?FIL1KM_HDF?${DAAC_FIL1KM}?g" \
  -e "s?FILMASK_HDF?${DAAC_FILMASK}?g" \
  -e "s?FILGEO_HDF?${DAAC_FILGEO}?g" \
  -e "s?FILMOD06_HDF?${DAAC_FILMOD06}?g" \
  -e "s?GDASBEFORE_NAME?${GDASbefore}?g" \
  -e "s?GDASAFTER_NAME?${GDASafter}?g" \
  -e "s?ANC_ICEC_NAME?${ICEC}?g" \
  -e "s?ANC_NISE_NAME?${HDFNISE}?g" \
  -e "s?CSR_BIAS_FILE_NAME?${CSR}?g" \
  -e "s?FORICE_NAME?${ICEL}?g" \
  -e "s?FORWAT_NAME?${WATERL}?g" \
  -e "s?FORICE_WS3_NAME?${ICELWS3}?g" \
  -e "s?FORICE_WS3SD_NAME?${ICELWS3SD}?g" \
  -e "s?FORICE_WS7_NAME?${ICELWS7}?g" \
  -e "s?FORICE_WS7SD_NAME?${ICELWS7SD}?g" \
  -e "s?FORICE_WS15_NAME?${ICELWS15}?g" \
  -e "s?FORICE_WS15SD_NAME?${ICELWS15SD}?g" \
  -e "s?FORWAT_WS3_NAME?${WATERLWS3}?g" \
  -e "s?FORWAT_WS3SD_NAME?${WATERLWS3SD}?g" \
  -e "s?FORWAT_WS7_NAME?${WATERLWS7}?g" \
  -e "s?FORWAT_WS7SD_NAME?${WATERLWS7SD}?g" \
  -e "s?FORWAT_WS15_NAME?${WATERLWS15}?g" \
  -e "s?FORWAT_WS15SD_NAME?${WATERLWS15SD}?g" \
  -e "s?FORICE_WATERPHASEFUNC_NAME?${ICEWATERPHASE}?g" \
  -e "s?TRANSLIB_NAME?${TRANW}?g" \
  -e "s?IGBPMAP_NAME?${ECOM}?g" \
  -e "s?ASNOW_NAME?${ALBSL}?g" \
  -e "s?SALB659_NAME?${ALB659}?g" \
  -e "s?SALB858_NAME?${ALB858}?g" \
  -e "s?SALB124_NAME?${ALB124}?g" \
  -e "s?SALB164_NAME?${ALB164}?g" \
  -e "s?SALB213_NAME?${ALB213}?g" \
  -e "s?CT_LIB_PATH_NAME?${CTCOEFFDIR}?g" \
  -e "s?EMISS_NAME?${EMISS}?g" \
  $TEMPLATE > cloud_optical.cfg

#-------------------------------------------------------------------------------
# RUN THE ALGORITHM
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Print start message for processing log
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Started MODIS optical properties processing at "`date`

# Run the cloud optical properties software
time ${MODIS_L2_BIN}/cloudoptical.exe cloud_optical.cfg

STATUS=$?

# Print failure or success message
if [ $STATUS != 0 ]; then
  echo "********** Error encountered creating cloud optical depth product **********"
  exit 1
else 
  echo "********** Successful creation of cloud optical depth product **********"
fi


echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Renaming product files to direct broadcast convention..."

if [ $PREFIX == "MOD" ]; then
  DB_HEADER="t1"
else
  DB_HEADER="a1"
fi

if [ `echo $DATE_TIME | cut -c1-2` == "XX" ]; then
  DB_DATETIME=`echo $DATE_TIME | cut -c3-12`
else
  DB_DATETIME=$DATE_TIME
fi

/bin/mv $OUTPUTDIRNAME/OD/$DAAC_FILMOD06 ${DB_HEADER}.${DB_DATETIME}.mod06.hdf

# Rename the configuration file
/bin/mv cloud_optical.cfg ${DB_HEADER}.${DB_DATETIME}.mod06od.cfg

# Move the output files into the output directory
#/bin/mv $FILMOD06 $ROOT.mod06od.cfg $OUTDIR
/bin/mv ${DB_HEADER}.${DB_DATETIME}.mod06.hdf ${DB_HEADER}.${DB_DATETIME}.mod06od.cfg $OUTDIR

# Clean up
if [ $STATUS == 0 ]; then
  rm -rf inputs-${YEAR}${DAY}.${TIME} outputs-${YEAR}${DAY}.${TIME} 
  rm -f ${FILEMOD06}
fi

# Print finish message for processing log
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Cloud optical properties processing ended at "`date`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo
#-------------------------------------------------------------------------------

# Exit gracefully
exit 0
