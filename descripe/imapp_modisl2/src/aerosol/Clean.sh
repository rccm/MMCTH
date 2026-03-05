#!/bin/bash

#---------------------------------------------------------------------
#
# Clean script for Collection 6 MODIS-atmos Level2 aerosol software.
# Run before running Build.sh to build the executables. 
#
# R. Cintineo 2013
#
#---------------------------------------------------------------------

ROOT=${MODIS_L2_HOME}

cd $ROOT/src/aerosol/MOD_PR04CR/src; make clean -f MOD_PR04CR1.mk
cd $ROOT/src/aerosol/MOD_PR04CR/src; make clean -f MOD_PR04CR2.mk
cd $ROOT/src/aerosol/MOD_PR04_05/src; make clean -f MOD_PR04_05.mk
cd $ROOT/src/aerosol/MOD04_3K/src; make clean -f MOD04_3K.mk
cd $ROOT/src/aerosol/MOD_PR04DB/src; make clean -f MOD_PR04DB.mk

cd $ROOT
