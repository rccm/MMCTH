#gfortran -c modis_planck.f90 modis_surfemis.f90 modis_modrad_emis.f90 
#gfortran -c modis_transmission_module.f90
#f2py -c  -m modis  *.f90
#rm *.o *.mod
rm -f *.o *.mod *.so
f2py -c --f90flags="-g -O0  -fallow-argument-mismatch  -fcheck=all -fbacktrace -fno-inline" -m modis modis_planck.f90 modis_surfemis.f90 modis_modrad_emis.f90 modis_transmission_module.f90 modis_sub4.f90  modis_co2_slice.f90 -DF2PY_REPORT_ON_ARRAY_COPY=1  
