#!/bin/bash

#---------------------------------------------------------------------
#
# Build script for Collection 6 MODIS-atmos Level2 aerosol software 
# Run Clean.sh before running this script.
#
# R. Cintineo 2013
#
# NOTE: All products are compiled with gfortran gcc version 4.4.7 
#
#---------------------------------------------------------------------

ROOT=${MODIS_L2_HOME}

source $ROOT/Setup_gfortran.sh

# Compile aerosol product
cd $ROOT/src/aerosol/MOD_PR04CR/src; make -f MOD_PR04CR1.mk
cd $ROOT/src/aerosol/MOD_PR04CR/src; make -f MOD_PR04CR2.mk
cd $ROOT/src/aerosol/MOD_PR04_05/src; make -f MOD_PR04_05.mk
cd $ROOT/src/aerosol/MOD04_3K/src; make -f MOD04_3K.mk
cd $ROOT/src/aerosol/MOD_PR04DB/src; make -f MOD_PR04DB.mk

cd $ROOT
