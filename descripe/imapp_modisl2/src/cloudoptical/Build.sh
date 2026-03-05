#!/bin/bash

#---------------------------------------------------------------------
#
# Build script for Collection 6 MODIS-atmos Level2 cloud optical 
# properties software 
#
# R. Cintineo 2013
#
# NOTE: All products are compiled with gfortran gcc version 4.4.7.
# The toolkit must be compiled with the same compiler.  You can contact 
# Rebecca Cintineo at Rebecca.Cintineo@ssec.wisc.edu if you want help to 
# recompile the cloud optical properties software with a different compiler.
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

source $ROOT/src/Setup_gfortran.sh

# Compile cloud optical properties product
cd $ROOT/src/cloudoptical/src/MOD_PR06OD; make -f MOD_PR06OD.mk.linux

cd $ROOT
