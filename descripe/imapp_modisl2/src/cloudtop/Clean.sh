#!/bin/bash

#---------------------------------------------------------------------
#
# Build script for Collection 6 MODIS-atmos Level2 software package 
#
# R. Cintineo 2013
#
# NOTE: All products are compiled with gfortran gcc version 4.1.2 
# except for MOD_PR06CT which is compiled with gfortran gcc version 4.4.7.  
# The toolkit must be compiled with the same compiler.  You can contact 
# Rebecca Cintineo at Rebecca.Cintineo@ssec.wisc.edu if you want help to 
# recompile MOD_PR06CT with the gfortran 4.4.7 compiler.
#
#---------------------------------------------------------------------

ROOT=$MODIS_L2_HOME

cd $ROOT/src/cloudtop/MOD_PR06CD/src; make clean -f MOD_PR06CD.mk
cd $ROOT/src/cloudtop/MOD_PR06CR/src; make clean -f MOD_PR06CR.mk
cd $ROOT/src/cloudtop/MOD_PR06CT/src; make clean -f MOD_PR06CT.mk

cd $ROOT
