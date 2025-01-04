import sys
import os
import xarray as xr
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from src.data_processors import process_main
import numpy as np
import pytest 
# from memory_profiler import profile
import logging
from src.data_processors.process_main import MainProcessor
 
@pytest.fixture
def input_data():
    file =  '/data/keeling/a/gzhao1/f/mmcth/output/MODMISR_L2_CP_A2016074.0330_O086367_v01.nc'
    ds  = xr.open_dataset(file) 
    print(ds.variables)
    lat  =  ds['latitude'].values
    lon = ds['longitude'].values
    lats,lons = np.meshgrid(lon,lat)
    era_variables = {
     'sst': ds['sst'].values,
     't': ds['temperature_exp'].values,
     'lat': lats,
     'lon': lons,
    }
    print(era_variables['lat'].shape)
    return era_variables

def test_block(input_data):
    inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2016/074/MOD021KM.A2016074.0330.061.2017325033257.hdf', \
                  '/data/gdi/d/MOD06/2016/074/MOD06_L2.A2016074.0330.061.2017325190341.hdf', \
                    '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2016/074/MOD03.A2016074.0330.061.2017325025726.hdf', \
                        '/data/gdi/e/MISR/TCCloud/2016.03.14/MISR_AM1_TC_CLOUD_P126_O086367_F01_0001.hdf', \
                              '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P126_F01_24.hdf', \
                                '/data/gdi/e/ERA5//single/era5_single_2016_03_14.nc', \
                                    '/data/gdi/e/ERA5//multi/era5_profile_2016_03_14.nc', \
                                        '/data/gdi/e/ERA5/single/era5_single_2016_03_15.nc', \
                                              '/data/gdi/e/ERA5/multi/era5_profile_2016_03_15.nc']
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    # era = ERA5Processor(inputfiles,logger)
    # results = era.era5_process()
    mm = MainProcessor(inputfiles,logger)
    a = mm.block_average_dict(input_data)
    print(a)