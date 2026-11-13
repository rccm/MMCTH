/data/keeling/a/gzhao1/miniconda3/envs/k9/bin/python -m numpy.f2py -c -m interpolate_to_pressure_levels era5_to_modis_pressurelevels.f90 -DF2PY_REPORT_ON_ARRAY_COPY=1
/data/keeling/a/gzhao1/miniconda3/envs/k9/bin/python -m numpy.f2py -c -m misr_cth_to_pressure misr_cth_to_pressure.f90 -DF2PY_REPORT_ON_ARRAY_COPY=1
