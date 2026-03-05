#!/bin/bash

#---------------------------------------------------------------------
#
# Build script for Collection 6 MOD07 code in the 
# MODIS-atmos Level2 software package 
#
# R. Cintineo 2013
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

source $ROOT/Setup_gfortran.sh

# Compile profile retrievals product
cd $ROOT/src/profiles/src; make -f MOD_PR07.mk

cd $ROOT
