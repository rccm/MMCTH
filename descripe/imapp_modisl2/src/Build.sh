#!/bin/bash

#---------------------------------------------------------------------
#
# Build script for Collection 6 MODIS-atmos Level2 software package 
# 
# This script only needs to be run if you would like to recompile
# the software.  Newly compiled executables will be created in the
# /src directory for each product. Pre-compiled, statically-linked 
# executables are available in the imapp_modisl2/bin directory.
#
# R. Cintineo 2013
#
# NOTE: All products have been compiled with gfortran gcc version 4.4.7.
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

source $ROOT/Setup_gfortran.sh

# Compile cloud mask product
cd $ROOT/src/cloudmask/src; make -f MOD_PR35.mk

# Compile profile retrievals product
cd $ROOT/src/profiles/src; make -f MOD_PR07.mk

# Compile cloud top properties products
cd $ROOT/src/cloudtop/MOD_PR06CD/src; make -f MOD_PR06CD.mk
cd $ROOT/src/cloudtop/MOD_PR06CR/src; make -f MOD_PR06CR.mk
cd $ROOT/src/cloudtop/MOD_PR06CT/src; make -f MOD_PR06CT.mk 

# Compile cloud optical propertires product
cd $ROOT/src/cloudoptical/MOD_PR06OD/src; make -f MOD_PR06OD.mk.linux 

# Compile cloud top aerosol product
cd $ROOT/src/aerosol/MOD_PR04CR/src; make -f MOD_PR04CR1.mk
cd $ROOT/src/aerosol/MOD_PR04CR/src; make -f MOD_PR04CR2.mk
cd $ROOT/src/aerosol/MOD_PR04_05/src; make -f MOD_PR04_05.mk
cd $ROOT/src/aerosol/MOD04_3K/src; make -f MOD04_3K.mk
cd $ROOT/src/aerosol/MOD_PR04DB/src; make -f MOD_PR04DB.mk

# Compile LEOCAT algorithms
cd $ROOT/src/leocat/src/MOD_PRAlg17; make -f MOD_PRAlg17.mk
cd $ROOT/src/leocat/src/MOD_PRAlg29; make -f MOD_PRAlg29.mk
cd $ROOT/src/leocat/src/MOD_PRAlg29St; make -f MOD_PRAlg29St.mk
cd $ROOT/src/leocat/src/MOD_PRL2M35; make -f MOD_PRL2M35.mk
cd $ROOT/src/leocat/src/MOD_PRL2M06; make -f MOD_PRL2M06.mk
cd $ROOT/src/leocat/src; make -f MOD_PRLCAT.mk

# Compile dateplus
#cd $ROOT/src/dateplus/src; gcc -static -O -o dateplus.exe dateplus.c
cd $ROOT/src/dateplus/src; make -f makefile.static

cd $ROOT
