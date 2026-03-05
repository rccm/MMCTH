#!/bin/bash

ROOT=$MODIS_L2_HOME

echo $'\n'"#############################################################"
echo      "###### LEOCAT Wrapper script for qualified algorithms #######"
echo      "#############################################################"$'\n'

echo "Number of arguments:" $#

args=("$@")
i=0
while [ $i -lt $# ];
do
	echo arg#$i = ${args[$i]}
	let i=i+1
done


if [ $# -ne 15 ]
then
	echo "wrong number of args" 
	echo $'\t\t'"1) Location of required binaries" 
	echo $'\t\t'"2) Location of static data"
	echo $'\t\t'"3) Location of Dynamic ancilliary data"

	echo $'\t\t'"1) Full path of input geolocation file"
	echo $'\t\t'"2) Full path of level 1b input LAC"
	echo $'\t\t'"3) Full path of level 1b input Clear-sky radiance bias file"
	echo $'\t\t'"4) Full path of intermediate output level 2 HDF file"

	echo $'\t\t'"5) Satellite Name (aqua/terra)"
	echo $'\t\t'"6) Date (YYYYDDD)"
	echo $'\t\t'"7) Time (HHMM)"
	
        echo $'\t\t'"8) Full Path of first GDAS file"
        echo $'\t\t'"9) Full Path of second GDAS file"
	echo $'\t\t'"10) Full Path of NISE file"
	echo $'\t\t'"11) Full Path of Reynolds SST file"
	echo $'\t\t'"12) Full Path to C6 Cloud Mask file"

	exit 1
fi

export BINARY_DIR=$1
export STATIC_DATA_DIR=$2
export DYNAMIC_DATA_DIR=$3

export GEO_INPUT_FULLPATH=$4
export LEVEL1_INPUT_FULLPATH=$5
export CSRB_INPUT_FULLPATH=$6
export LEVEL2_OUTPUT_FULLPATH=$7

export SAT_NAME=$8
export SAT_DATE=$9
export SAT_TIME=${10}

export GDAS_FILE1=${11}
export GDAS_FILE2=${12}
export NISE_FILE=${13}
export RSST_FILE=${14}
export MSK_FILE=${15}

# Make sure all NWP files are in the same directory
cp $GDAS_FILE1 $GDAS_FILE2 $DYNAMIC_DATA_DIR

# Algorithms to run (in order given)

ALG_NUMBER="17"
#ALG_NUMBER="29"
#ALG_NUMBER="20"
#ALG_NUMBER="20 24"
#ALG_NUMBER="20 24 22"

echo $'\n'"ALG_NUMBER:" $ALG_NUMBER

echo $'\n'"#############################################################"
echo      "######                Binary Paths                    #######"
echo      "#############################################################"$'\n\n'

# Set some binary paths which persist only for the duration of this script
export PATH=$BINARY_DIR/modis-atmos/scripts:$PATH
export PATH=$BINARY_DIR/scripts:$PATH

# Directory containing leocat executable
LEOCAT_BIN_DIR=$ROOT/src/leocat/src

echo "Location of binaries      : " $BINARY_DIR
echo "Location of LEOCAT binary : " $LEOCAT_BIN_DIR
echo "PATH : "$'\n'$PATH

echo $'\n'"#############################################################"
echo      "######    Input and output granule files and paths    #######"
echo      "#############################################################"$'\n\n'

# Full path of input geolocation file
relPath=$4;
GEO_INPUT_FILE=$(basename $relPath);
echo "Geolocation file =      "$GEO_INPUT_FILE
GEO_INPUT_FULLPATH=$(cd ${relPath%$GEO_INPUT_FILE};pwd)"/"$GEO_INPUT_FILE;
echo "Geolocation Full path = "$GEO_INPUT_FULLPATH;
GEO_INPUT_DIR=$(cd ${GEO_INPUT_FULLPATH%$GEO_INPUT_FILE};pwd);
echo "Geolocation path =      "$GEO_INPUT_DIR;

# Full path of input level 1b input HDF files
relPath=$5;
LEVEL1B_INPUT_FILE=$(basename $relPath);
echo $'\n'"Level 1b input file ="$'\t'$LEVEL1B_INPUT_FILE
LEVEL1_INPUT_FULLPATH=$(cd ${relPath%$LEVEL1B_INPUT_FILE};pwd)"/"$LEVEL1B_INPUT_FILE;
echo "Level 1b Full path ="$'\t'$LEVEL1_INPUT_FULLPATH;
LEVEL1_INPUT_DIR=$(cd ${LEVEL1_INPUT_FULLPATH%$LEVEL1B_INPUT_FILE};pwd);
echo "Level 1b path ="$'\t\t'$LEVEL1_INPUT_DIR;

# Full path of input level 1b input clear-sky radiacne bias files
relPath=$6;
if [ "$relPath" != "MISSING" ]; then
  CSRB_INPUT_FILE=$(basename $relPath);
  echo $'\n'"CSRB input file ="$'\t'$CSRB_INPUT_FILE
  CSRB_INPUT_FULLPATH=$(cd ${relPath%$CSRB_INPUT_FILE};pwd)"/"$CSRB_INPUT_FILE;
  echo "CSRB Full path ="$'\t'$CSRB_INPUT_FULLPATH;
  CSRB_INPUT_DIR=$(cd ${CSRB_INPUT_FULLPATH%$CSRB_INPUT_FILE};pwd);
  echo "CSRB path ="$'\t\t'$CSRB_INPUT_DIR;
fi

# Full path of output level 2 output HDF files
relPath=$7;
LEVEL2_OUTPUT_FILE=$(basename $relPath);
echo $'\n'"Level 2 output file ="$'\t'$LEVEL2_OUTPUT_FILE
LEVEL2_OUTPUT_FULLPATH=$(cd ${relPath%$LEVEL2_OUTPUT_FILE};pwd)"/"$LEVEL2_OUTPUT_FILE;
echo "Level 2 Full path ="$'\t'$LEVEL2_OUTPUT_FULLPATH;
LEVEL2_OUTPUT_DIR=$(cd ${LEVEL2_OUTPUT_FULLPATH%$LEVEL2_OUTPUT_FILE};pwd);
echo "Level 2 path ="$'\t\t'$LEVEL2_OUTPUT_DIR;

echo $'\n'"#############################################################"
echo      "######        Get GDAS and NISE ancillary data        #######"
echo      "#############################################################"$'\n\n'

# Location of static surface data
SFC_DATA_DIR=$STATIC_DATA_DIR/sfc_data

# Where to find algorithm specific data
ALGDATA_DIR=$STATIC_DATA_DIR/sfc_data

## Location of NWP dynamic ancilliary data
NWP_DIR=$DYNAMIC_DATA_DIR

# ssec modis software
echo "Location of static ancillary data         : " $STATIC_DATA_DIR
echo "Location of static surface ancillary data : " $SFC_DATA_DIR
echo "Location of dynamic ancillary data        : " $DYNAMIC_DATA_DIR
echo "Location of NWP dynamic ancillary data    : " $NWP_DIR

echo $'\n'"#############################################################"
echo      "######        Get GDAS and NISE ancillary data        #######"
echo      "#############################################################"$'\n'

export time=${wholeTime:0:4}
echo "sat_name           :" $SAT_NAME
echo "date               :" $SAT_DATE
echo "time               :" $SAT_TIME
echo ""

case $SAT_NAME in
	"terra") 
		DATA_SRC="modis";
		;;
	"aqua") 
		DATA_SRC="modis";
		;;
	"npp") 
		DATA_SRC="viirs";
		;;
	*       ) 
		;;
esac

echo "Data source: "$DATA_SRC;
	
# Check for required ancillary data (GDAS and NISE)
if [ ! -r $GDAS_FILE1 ]; then
  echo "First GDAS file not found: "$GDAS_FILE1
  exit 1
fi
if [ ! -r $GDAS_FILE2 ]; then
  echo "Second GDAS file not found: "$GDAS_FILE2
  exit 1
fi
if [ ! -r $NISE_FILE ]; then
  echo "NISE file not found: "$NISE_FILE
  exit 1
fi

if [ ! -r $RSST_FILE ]; then
  echo "RSST file not found: "$RSST_FILE
  exit 1
fi

NWP_TYPE="gdas"

echo $'\n'"#############################################################"
echo      "######           Parameters to pass to LEOCAT         #######"
echo      "#############################################################"$'\n'

#	-satType		Satellite name/type
#	-satYYYDDD		Year and Day of granule
#	-satHHMM		Hour and minute of granule
#	-verbose :		Spit out more information
#	-wd :			Where to run the executable in
#	-nol3 :			No Level 3 data to be created
#	-a :			Algorithm number
#	-geo_file :		Full path of geolocation file
#	-gdas_file :	Ancillary input gdas file
#	-nise_file :	Ancillary input nise file
#	-l2_dir :		Location of level 2 files
#	-l2_name :		Name of level 2 output file
#	-nwp :			Type of nwp
#	-nwp_dir :		Location of nwp file
#	-l1b_dir :		Location of level 1b files
#	-f :			Level 1b input file
#   -algData_dir    Location of algorithm specific data
#   -staticData_dir Location of static data
#	-sfc_dir :		Location of surface data within static data
echo "Satellite name          (-satType) :"	$SAT_NAME
echo "Satellite date        (-satYYYDDD) :"	$SAT_DATE
echo "Satellite time          (-satHHMM) :"	$SAT_TIME
echo "Verbosity               (-verbose)  " 
echo "Working directory            (-wd) :" $LEOCAT_BIN_DIR
echo "No Level 3 output:         (-nol3)  "		   
echo "Algorithm number(s)           (-a) :" $ALG_NUMBER
echo "Algorithm data dir  (-algData_dir) :" $ALGDATA_DIR
echo "Static data dir  (-staticData_dir) :" $STATIC_DATA_DIR
echo "Surface data dir        (-sfc_dir) :" $SFC_DATA_DIR
echo "NWP type                    (-nwp) :"	$NWP_TYPE
echo "NWP directory           (-nwp_dir) :"	$NWP_DIR
echo "First GDAS file      (-gdas_file1) :" $GDAS_FILE1
echo "Second GDAS file     (-gdas_file2) :" $GDAS_FILE2
echo "NISE file             (-nise_file) :" $NISE_FILE
echo "Reynolds SST file      (-sst_file) :" $RSST_FILE
echo "Geolocation file       (-geo_file) :" $GEO_INPUT_FULLPATH
echo "Level 2 output directory (-l2_dir) :" $LEVEL2_OUTPUT_DIR
echo "Level 2 output file     (-l2_name) :" $LEVEL2_OUTPUT_FILE
echo "Level 1 input directory (-l1b_dir) :"	$LEVEL1_INPUT_DIR
echo "Level 1b input file           (-f) :"	$LEVEL1B_INPUT_FILE
echo "CSRB directory       (-csrb_dir):"	$CSRB_INPUT_DIR
echo "CSRB input file         (-csrb_file) :"	$CSRB_INPUT_FILE

BASE_GDAS_FILE1=/$(basename $GDAS_FILE1);
BASE_GDAS_FILE2=/$(basename $GDAS_FILE2);
echo $'\n'"#############################################################"
echo      "######             Command line to execute...         #######"
echo      "#############################################################"$'\n'

if [ "$SAT_NAME" != "aqua" ]; then
  echo "$ROOT/bin/leocat.exe -verbose -satType $SAT_NAME -satYYYDDD $SAT_DATE -satHHMM $SAT_TIME  -wd $LEOCAT_BIN_DIR -nol3 -a $ALG_NUMBER -geo_file $GEO_INPUT_FULLPATH -nise_file $NISE_FILE -sst_file $RSST_FILE -l2_dir $LEVEL2_OUTPUT_DIR -l2_name $LEVEL2_OUTPUT_FILE -nwp $NWP_TYPE -nwp_dir $NWP_DIR -l1b_dir $LEVEL1_INPUT_DIR -f $LEVEL1B_INPUT_FILE -cmask_file $MSK_FILE -algData_dir $ALGDATA_DIR -sfc_dir $SFC_DATA_DIR -rtm plod -csrb_dir $CSRB_INPUT_DIR -csrb_file $CSRB_INPUT_FILE -gdas_file1 $BASE_GDAS_FILE1 -gdas_file2 $BASE_GDAS_FILE2 -terra_plod_spectral_shift"$'\n\n'

  $ROOT/bin/leocat.exe -verbose -satType $SAT_NAME -satYYYDDD $SAT_DATE -satHHMM $SAT_TIME  -wd $LEOCAT_BIN_DIR -nol3 -a $ALG_NUMBER -geo_file $GEO_INPUT_FULLPATH -nise_file $NISE_FILE -sst_file $RSST_FILE -l2_dir $LEVEL2_OUTPUT_DIR -l2_name $LEVEL2_OUTPUT_FILE -nwp $NWP_TYPE -nwp_dir $NWP_DIR -l1b_dir $LEVEL1_INPUT_DIR -f $LEVEL1B_INPUT_FILE -cmask_file $MSK_FILE -algData_dir $ALGDATA_DIR -sfc_dir $SFC_DATA_DIR -rtm plod -csrb_dir $CSRB_INPUT_DIR -csrb_file $CSRB_INPUT_FILE -gdas_file1 $BASE_GDAS_FILE1 -gdas_file2 $BASE_GDAS_FILE2 -terra_plod_spectral_shift

else
  echo "$ROOT/bin/leocat.exe -verbose -satType $SAT_NAME -satYYYDDD $SAT_DATE -satHHMM $SAT_TIME  -wd $LEOCAT_BIN_DIR -nol3 -a $ALG_NUMBER -geo_file $GEO_INPUT_FULLPATH -nise_file $NISE_FILE -sst_file $RSST_FILE -l2_dir $LEVEL2_OUTPUT_DIR -l2_name $LEVEL2_OUTPUT_FILE -nwp $NWP_TYPE -nwp_dir $NWP_DIR -l1b_dir $LEVEL1_INPUT_DIR -f $LEVEL1B_INPUT_FILE -cmask_file $MSK_FILE -algData_dir $ALGDATA_DIR -sfc_dir $SFC_DATA_DIR -rtm plod -csrb_dir $CSRB_INPUT_DIR -csrb_file $CSRB_INPUT_FILE -gdas_file1 $BASE_GDAS_FILE1 -gdas_file2 $BASE_GDAS_FILE2 -aqua_plod_spectral_shift"$'\n\n'

  $ROOT/bin/leocat.exe -verbose -satType $SAT_NAME -satYYYDDD $SAT_DATE -satHHMM $SAT_TIME  -wd $LEOCAT_BIN_DIR -nol3 -a $ALG_NUMBER -geo_file $GEO_INPUT_FULLPATH -nise_file $NISE_FILE -sst_file $RSST_FILE -l2_dir $LEVEL2_OUTPUT_DIR -l2_name $LEVEL2_OUTPUT_FILE -nwp $NWP_TYPE -nwp_dir $NWP_DIR -l1b_dir $LEVEL1_INPUT_DIR -f $LEVEL1B_INPUT_FILE -cmask_file $MSK_FILE -algData_dir $ALGDATA_DIR -sfc_dir $SFC_DATA_DIR -rtm plod -csrb_dir $CSRB_INPUT_DIR -csrb_file $CSRB_INPUT_FILE -gdas_file1 $BASE_GDAS_FILE1 -gdas_file2 $BASE_GDAS_FILE2 -aqua_plod_spectral_shift

fi

if [ $? == 1 ]; then
  echo
  echo "********** WARNING: LEOCAT failed for cloud top properties **********"
  echo
  exit 1
else
  echo
  echo "********** Successful completion of LEOCAT for cloud top properties **********"
  echo
fi

if [ -e $LEVEL2_OUTPUT_FILE ]; then
	mv -v $LEVEL2_OUTPUT_FILE $LEVEL2_OUTPUT_DIR;
fi

exit 0
