rm -f modis_edf_destripe*.so
FC=gfortran F77=gfortran F90=gfortran \
python -m numpy.f2py -c -m modis_edf_destripe modis_edf_destripe.f90 interp.f --opt="-O3"
