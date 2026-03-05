#!/bin/bash

#---------------------------------------------------------------------
#
# Build script for Collection 6 MODIS-atmos Level2 software package 
#
# R. Cintineo 2013
#
# NOTE: All products have been compiled with gfortran gcc version 4.4.7 
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

source $ROOT/Setup_gfortran.sh

# Compile dateplus
#cd $ROOT/src/dateplus/src; gcc -static -O -o dateplus.exe dateplus.c
cd $ROOT/src/dateplus/src; make clean -f makefile.static

cd $ROOT
