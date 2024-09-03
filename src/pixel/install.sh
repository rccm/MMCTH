!/bin/bash

gfortran -c modis_planck.f90 modis_surfemis.f90 modis_modrad_emis.f90 
gfortran -c modis_transmission_module.f90
f2py -c  -m modis *.f90
rm *.o *.mod
