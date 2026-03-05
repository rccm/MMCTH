#!/bin/bash

#---------------------------------------------------------------------
#
# Clean script for Collection 6 MODIS-atmos Level2 cloudmask software.
# Run before running Build.sh script to build executable. 
#
# R. Cintineo 2013
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

cd $ROOT/src/cloudmask/src; make clean -f MOD_PR35.mk

cd $ROOT
