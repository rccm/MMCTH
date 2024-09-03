
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
    inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2002/337/MOD021KM.A2002337.0215.061.2017185135616.hdf', \
                 '/data/gdi/d/MOD06/2002/337/MOD06_L2.A2002337.0215.061.2017280150853.hdf',\
                 '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2002/337/MOD03.A2002337.0215.061.2017185133149.hdf',\
                 '/data/gdi/e/MISR/TCCloud/2002.12.03/MISR_AM1_TC_CLOUD_P112_O015738_F01_0001.hdf',\
                 '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P112_F01_24.hdf',\
                 '/data/gdi/e/ERA5/single/era5_single_2002_12_03.nc',\
                 '/data/gdi/e/ERA5/multi/era5_profile_2002_12_03.nc',\
                 'X',\
                 'X',
            ]
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    # era = ERA5Processor(inputfiles,logger)
    # results = era.era5_process()
    mm = MainProcessor(inputfiles,logger)
    mm.run_process()