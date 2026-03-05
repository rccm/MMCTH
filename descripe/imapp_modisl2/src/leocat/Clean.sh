#!/bin/bash

#---------------------------------------------------------------------
#
# Clean script for Collection 6 MODIS-atmos Level2 LEOCAT software. 
# Run before running Build.sh script to build executables.
# R. Cintineo 2013
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

cd $ROOT/src/leocat/src/MOD_PRAlg17; make clean -f MOD_PRAlg17.mk
cd $ROOT/src/leocat/src/MOD_PRAlg29; make clean -f MOD_PRAlg29.mk
cd $ROOT/src/leocat/src/MOD_PRAlg29St; make clean -f MOD_PRAlg29St.mk
cd $ROOT/src/leocat/src/MOD_PRL2M35; make clean -f MOD_PRL2M35.mk
cd $ROOT/src/leocat/src/MOD_PRL2M06; make clean -f MOD_PRL2M06.mk
cd $ROOT/src/leocat/src; make clean -f MOD_PRLCAT.mk

cd $ROOT
