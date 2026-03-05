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

source $ROOT/Setup_gfortran.sh

# Compile LEOCAT algorithms
#cd $ROOT/src/MOD_PRLCAT/MOD_PRAlg17; make -f MOD_PRAlg17.mk
#cd $ROOT/src/MOD_PRLCAT/MOD_PRAlg29; make -f MOD_PRAlg29.mk
#cd $ROOT/src/MOD_PRLCAT/MOD_PRAlg29St; make -f MOD_PRAlg29St.mk
#cd $ROOT/src/MOD_PRLCAT/MOD_PRL2M35; make -f MOD_PRL2M35.mk
#cd $ROOT/src/MOD_PRLCAT/MOD_PRL2M06; make -f MOD_PRL2M06.mk
#cd $ROOT/src/MOD_PRLCAT; make -f MOD_PRLCAT.mk

cd $ROOT/src/leocat/src/MOD_PRAlg17; make -f MOD_PRAlg17.mk
cd $ROOT/src/leocat/src/MOD_PRAlg29; make -f MOD_PRAlg29.mk
cd $ROOT/src/leocat/src/MOD_PRAlg29St; make -f MOD_PRAlg29St.mk
cd $ROOT/src/leocat/src/MOD_PRL2M35; make -f MOD_PRL2M35.mk
cd $ROOT/src/leocat/src/MOD_PRL2M06; make -f MOD_PRL2M06.mk
cd $ROOT/src/leocat/src; make -f MOD_PRLCAT.mk

cd $ROOT
