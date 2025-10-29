
import logging
import sys
import os

# Adjust the path as needed
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from src.data_processors.process_main import ERA5Processor
from src.data_processors.process_main import MODISMISRProcessor 
from src.data_processors.process_main import MainProcessor

if __name__ == "__main__":
    import logging
    inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2016/338/MOD021KM.A2016338.1615.061.2017328222457.hdf', \
                  '/data/gdi/d/MOD06/2016/338/MOD06_L2.A2016338.1615.061.2017329080242.hdf', \
                    '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2016/338/MOD03.A2016338.1615.061.2017328214312.hdf', \
                        '/data/gdi/e/MISR/TCCloud/2016.12.03/MISR_AM1_TC_CLOUD_P013_O090219_F01_0001.hdf', \
                            '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P013_F01_24.hdf', \
                                '/data/gdi/e/ERA5//single/era5_single_2016_12_03.nc', \
                                    '/data/gdi/e/ERA5//multi/era5_profile_2016_12_03.nc', \
                                        '/data/gdi/e/ERA5/single/era5_single_2016_12_04.nc', \
                                            '/data/gdi/e/ERA5/multi/era5_profile_2016_12_04.nc']
    inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2016/074/MOD021KM.A2016074.0330.061.2017325033257.hdf', \
                  '/data/gdi/d/MOD06/2016/074/MOD06_L2.A2016074.0330.061.2017325190341.hdf', \
                    '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2016/074/MOD03.A2016074.0330.061.2017325025726.hdf', \
                        '/data/gdi/e/MISR/TCCloud/2016.03.14/MISR_AM1_TC_CLOUD_P126_O086367_F01_0001.hdf', \
                              '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P126_F01_24.hdf', \
                                '/data/gdi/e/ERA5//single/era5_single_2016_03_14.nc', \
                                    '/data/gdi/e/ERA5//multi/era5_profile_2016_03_14.nc', \
                                        '/data/gdi/e/ERA5/single/era5_single_2016_03_15.nc', \
                                              '/data/gdi/e/ERA5/multi/era5_profile_2016_03_15.nc']
    
    # 2040 size not 2030 figure out why
    inputfiles =  ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2007/003/MOD021KM.A2007003.1345.061.2017245234514.hdf',\
                   '/data/gdi/d/MOD06/2007/003/MOD06_L2.A2007003.1345.061.2017279021020.hdf', \
                   '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2007/003/MOD03.A2007003.1345.061.2017245074830.hdf', \
                   '/data/gdi/e/MISR/TCCloud/2007.01.03/MISR_AM1_TC_CLOUD_P220_O037472_F01_0001.hdf', \
                   '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P220_F01_24.hdf', '/data/gdi/e/ERA5//single/era5_single_2007_01_03.nc', \
                   '/data/gdi/e/ERA5//multi/era5_profile_2007_01_03.nc', '/data/gdi/e/ERA5/single/era5_single_2007_01_04.nc', \
                   '/data/gdi/e/ERA5/multi/era5_profile_2007_01_04.nc']
    
    inputfiles =  ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2007/003/MOD021KM.A2007003.1345.061.2017245234514.hdf',\
                   '/data/gdi/d/MOD06/2007/003/MOD06_L2.A2007003.1345.061.2017279021020.hdf', \
                   '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2007/003/MOD03.A2007003.1345.061.2017245074830.hdf', \
                   'X', \
                   '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P220_F01_24.hdf', '/data/gdi/e/ERA5//single/era5_single_2007_01_03.nc', \
                   '/data/gdi/e/ERA5//multi/era5_profile_2007_01_03.nc', '/data/gdi/e/ERA5/single/era5_single_2007_01_04.nc', \
                   '/data/gdi/e/ERA5/multi/era5_profile_2007_01_04.nc']
    
 
    
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    # era = ERA5Processor(inputfiles,logger = logger)
    # results = era.era5_process()
    mm = MainProcessor(inputfiles,logger = logger)
    mm.run_process(save_flag = 'normal')