#!/bin/bash

#########################################################################
#									#
# This script sets the environmental variables needed for compiling	#
# products in the IMAPP Level 2 software with gfortran 4.4.7		#
# R. Cintineo 2013							#
#                                                                       #
#########################################################################

#--- Set home directory of MODIS Level 2 processing ---#
export MODIS_L2_HOME=/home/imapp_modisl2
export MODIS_STORE=$MODIS_L2_HOME

#--- Set environmental variables for compiling MOD_PR07 --#
export TOOLKIT_HOME=$MODIS_L2_HOME/toolkits/TOOLKIT_gfortran44_f90
export PGSINC=$TOOLKIT_HOME/include
export PGSLIB=$TOOLKIT_HOME/lib/linux
export PGSBIN=$TOOLKIT_HOME/bin/linux

#--- Run toolkit setup script ---#
source $TOOLKIT_HOME/bin/linux/pgs-dev-env.ksh

#--- Set MAPI toolkit environment variables ---#
export API_HOME=$MODIS_L2_HOME/toolkits/mapi6.0.1
export API_INC=$API_HOME/h
export API_LIB=$API_HOME/lib/linux

#--- Set HDF environment variables ---#
export HDFEOS_INC=$TOOLKIT_HOME/hdfeos/include
export HDFEOS_LIB=$TOOLKIT_HOME/hdfeos/lib/linux
export HDFEOS5_LIB=$TOOLKIT_HOME/hdfeos5/lib/linux
export HDFINC=$TOOLKIT_HOME/hdf/linux64/hdf-4.2.6/include
export HDFLIB=$TOOLKIT_HOME/hdf/linux64/hdf-4.2.6/lib
export HDF5LIB=$TOOLKIT_HOME/hdf5/linux/hdf5-1.8.8/lib
export HDF5INC=$TOOLKIT_HOME/hdf5/linux/hdf5-1.8.8/include

if [ $LD_LIBRARY_PATH ]; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$HDF5LIB # For hdf5 Shared libs
else
    export LD_LIBRARY_PATH=$HDF5LIB
fi

#--- Make sure MODIS_39500 and _____ files are created ---#
$PGSBIN/smfcompile -f $MODIS_L2_HOME/src/cloudmask/shared_src/src_L2/MODIS_39500.t -all
$PGSBIN/smfcompile -f $TOOLKIT_HOME/message/MODIS_39500.t -all -i

export BRAND="linux"

export CFLAGS="-DLINUX -Df2cFortran -DSYSV"

# Set to compile with gfortran
#export F77='/usr/bin/gfortran -v'
#export F90='/usr/bin/gfortran -v'
export F77='/usr/bin/gfortran44 -v'
export F90='/usr/bin/gfortran44 -v'

export F77FLAGS="-DLINUX -ffixed-line-length-none"
export F90FLAGS="-DLINUX -ffixed-line-length-none"
