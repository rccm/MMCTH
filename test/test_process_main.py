
import logging
import sys
import os

# Adjust the path as needed
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

# from src.data_processors.process_main import ERA5Processor
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
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    # era = ERA5Processor(inputfiles,logger)
    # results = era.era5_process()
    mm = MainProcessor(inputfiles,logger)
    mm.run_process()