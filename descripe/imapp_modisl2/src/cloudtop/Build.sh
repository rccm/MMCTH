#!/bin/bash

#---------------------------------------------------------------------
#
# Build script for Collection 6 MODIS-atmos Level2 cloud top software 
#
# R. Cintineo 2013
#
# NOTE: All products are compiled with gfortran gcc version 4.4.7.
# The toolkit must be compiled with the same compiler.  You can contact 
# Rebecca Cintineo at Rebecca.Cintineo@ssec.wisc.edu if you want help to 
# recompile the cloud top software with a different compiler.
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

source $ROOT/Setup_gfortran.sh

# Compile cloud top properties products
cd $ROOT/src/cloudtop/MOD_PR06CD/src; make -f MOD_PR06CD.mk
cd $ROOT/src/cloudtop/MOD_PR06CR/src; make -f MOD_PR06CR.mk
cd $ROOT/src/cloudtop/MOD_PR06CT/src; make -f MOD_PR06CT.mk

cd $ROOT
