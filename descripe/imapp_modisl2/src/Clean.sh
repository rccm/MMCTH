#!/bin/bash

#---------------------------------------------------------------------
#
# Clean script for Collection 6 MODIS-atmos Level2 software package 
#
# Run this script to clean up the /src directories for each product
# before recompiling the code.
#
# R. Cintineo 2013
#
# NOTE: All products are compiled with gfortran gcc version 4.4.7.
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

# Clean cloud mask product directory
cd $ROOT/src/cloudmask/src; make clean -f MOD_PR35.mk

# Clean profile retrievals product directory
cd $ROOT/src/profiles/src; make clean -f MOD_PR07.mk

# Clean cloud top properties products directory
cd $ROOT/src/cloudtop/MOD_PR06CD/src; make clean -f MOD_PR06CD.mk
cd $ROOT/src/cloudtop/MOD_PR06CR/src; make clean -f MOD_PR06CR.mk
cd $ROOT/src/cloudtop/MOD_PR06CT/src; make clean -f MOD_PR06CT.mk 

# Clean cloud optical propertires product directory
cd $ROOT/src/cloudoptical/src/MOD_PR06OD; make clean -f MOD_PR06OD.mk.linux

# Clean cloud top aerosol product directory
cd $ROOT/src/aerosol/MOD_PR04CR/src; make clean -f MOD_PR04CR1.mk
cd $ROOT/src/aerosol/MOD_PR04CR/src; make clean -f MOD_PR04CR2.mk
cd $ROOT/src/aerosol/MOD_PR04_05/src; make clean -f MOD_PR04_05.mk
cd $ROOT/src/aerosol/MOD04_3K/src; make clean -f MOD04_3K.mk
cd $ROOT/src/aerosol/MOD_PR04DB/src; make clean -f MOD_PR04DB.mk

# Clean LEOCAT algorithms directory
cd $ROOT/src/leocat/src/MOD_PRAlg17; make clean -f MOD_PRAlg17.mk
cd $ROOT/src/leocat/src/MOD_PRAlg29; make clean -f MOD_PRAlg29.mk
cd $ROOT/src/leocat/src/MOD_PRAlg29St; make clean -f MOD_PRAlg29St.mk
cd $ROOT/src/leocat/src/MOD_PRL2M35; make clean -f MOD_PRL2M35.mk
cd $ROOT/src/leocat/src/MOD_PRL2M06; make clean -f MOD_PRL2M06.mk
cd $ROOT/src/leocat/src; make clean -f MOD_PRLCAT.mk

# Clean dateplus directory directory
cd $ROOT/src/dateplus/src; make clean -f makefile.static

cd $ROOT
