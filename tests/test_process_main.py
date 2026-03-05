
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
    inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2015/176/MOD021KM.A2015176.1650.061.2017321013930.hdf', \
                   '/data/gdi/d/MOD06/2015/176/MOD06_L2.A2015176.1650.061.2017321130959.hdf', \
                    '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2015/176/MOD03.A2015176.1650.061.2017321001513.hdf', \
                        '/data/gdi/e/MISR/TCCloud/2015.06.25/MISR_AM1_TC_CLOUD_P020_O082545_F01_0001.hdf', \
                            '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P020_F01_24.hdf', \
                                '/data/gdi/e/ERA5//single/era5_single_2015_06_25.nc', \
                                    '/data/gdi/e/ERA5//multi/era5_profile_2015_06_25.nc', \
                                        '/data/gdi/e/ERA5/single/era5_single_2015_06_26.nc', \
                                            '/data/gdi/e/ERA5/multi/era5_profile_2015_06_26.nc']
    # inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2015/144/MOD021KM.A2015144.1140.061.2017320181302.hdf', '/data/gdi/d/MOD06/2015/144/MOD06_L2.A2015144.1140.061.2017321025424.hdf', '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2015/144/MOD03.A2015144.1140.061.2017320173030.hdf', '/data/gdi/e/MISR/TCCloud/2015.05.24/MISR_AM1_TC_CLOUD_P205_O082076_F01_0001.hdf', '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P205_F01_24.hdf', '/data/gdi/e/ERA5//single/era5_single_2015_05_24.nc', '/data/gdi/e/ERA5//multi/era5_profile_2015_05_24.nc', '/data/gdi/e/ERA5/single/era5_single_2015_05_25.nc', '/data/gdi/e/ERA5/multi/era5_profile_2015_05_25.nc']
    # inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2015/347/MOD021KM.A2015347.2310.061.2017323105434.hdf', '/data/gdi/d/MOD06/2015/347/MOD06_L2.A2015347.2310.061.2017324012417.hdf', '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2015/347/MOD03.A2015347.2310.061.2017323063532.hdf', '/data/gdi/e/MISR/TCCloud/2015.12.13/MISR_AM1_TC_CLOUD_P081_O085039_F01_0001.hdf', '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P081_F01_24.hdf', \
    #                '/data/gdi/e/ERA5//single/era5_single_2015_12_13.nc', '/data/gdi/e/ERA5//multi/era5_profile_2015_12_13.nc',\
    #                   '/data/gdi/e/ERA5/single/era5_single_2015_12_14.nc', '/data/gdi/e/ERA5/multi/era5_profile_2015_12_14.nc']   
    
    # inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2015/267/MOD021KM.A2015267.2315.061.2017322052133.hdf', '/data/gdi/d/MOD06/2015/267/MOD06_L2.A2015267.2315.061.2017322235315.hdf', '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2015/267/MOD03.A2015267.2315.061.2017322044744.hdf', '/data/gdi/e/MISR/TCCloud/2015.09.24/MISR_AM1_TC_CLOUD_P081_O083874_F01_0001.hdf', '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P081_F01_24.hdf', '/data/gdi/e/ERA5//single/era5_single_2015_09_24.nc', '/data/gdi/e/ERA5//multi/era5_profile_2015_09_24.nc', '/data/gdi/e/ERA5/single/era5_single_2015_09_25.nc', '/data/gdi/e/ERA5/multi/era5_profile_2015_09_25.nc']
    # inputfiles =  ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2016/184/MOD021KM.A2016184.0000.061.2017326171357.hdf',  \
    #                '/data/gdi/d/MOD06/2016/184/MOD06_L2.A2016184.0000.061.2017327045045.hdf', \
    #                 '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2016/184/MOD03.A2016184.0000.061.2017326114536.hdf', \
    #                       '/data/gdi/e/MISR/TCCloud/2016.07.02/MISR_AM1_TC_CLOUD_P096_O087967_F01_0001.hdf', \
    #                         '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P096_F01_24.hdf', \
    #                             '/data/gdi/e/ERA5//single/era5_single_2016_07_02.nc', \
    #                                   '/data/gdi/e/ERA5//multi/era5_profile_2016_07_02.nc', '/data/gdi/e/ERA5/single/era5_single_2016_07_03.nc', '/data/gdi/e/ERA5/multi/era5_profile_2016_07_03.nc']
    # # 2040 size not 2030 figure out why
    # inputfiles =  ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2007/003/MOD021KM.A2007003.1345.061.2017245234514.hdf',\
    #                '/data/gdi/d/MOD06/2007/003/MOD06_L2.A2007003.1345.061.2017279021020.hdf', \
    #                '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2007/003/MOD03.A2007003.1345.061.2017245074830.hdf', \
    #                '/data/gdi/e/MISR/TCCloud/2007.01.03/MISR_AM1_TC_CLOUD_P220_O037472_F01_0001.hdf', \
    #                '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P220_F01_24.hdf', '/data/gdi/e/ERA5//single/era5_single_2007_01_03.nc', \
    #                '/data/gdi/e/ERA5//multi/era5_profile_2007_01_03.nc', '/data/gdi/e/ERA5/single/era5_single_2007_01_04.nc', \
    #                '/data/gdi/e/ERA5/multi/era5_profile_2007_01_04.nc']
    
    # inputfiles =  ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2007/003/MOD021KM.A2007003.1345.061.2017245234514.hdf',\
    #                '/data/gdi/d/MOD06/2007/003/MOD06_L2.A2007003.1345.061.2017279021020.hdf', \
    #                '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2007/003/MOD03.A2007003.1345.061.2017245074830.hdf', \
    #                'X', \
    #                '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P220_F01_24.hdf', '/data/gdi/e/ERA5//single/era5_single_2007_01_03.nc', \
    #                '/data/gdi/e/ERA5//multi/era5_profile_2007_01_03.nc', '/data/gdi/e/ERA5/single/era5_single_2007_01_04.nc', \
    #                '/data/gdi/e/ERA5/multi/era5_profile_2007_01_04.nc']
    # inputfiles = ['/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2015/168/MOD021KM.A2015168.1835.061.2017321012713.hdf', \
    #               '/data/gdi/d/MOD06/2015/168/MOD06_L2.A2015168.1835.061.2017321090544.hdf', \
    #                 '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2015/168/MOD03.A2015168.1835.061.2017321000710.hdf', \
    #                       '/data/gdi/e/MISR/TCCloud/2015.06.17/MISR_AM1_TC_CLOUD_P044_O082430_F01_0001.hdf', \
    #                         '/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P044_F01_24.hdf', \
    #                             '/data/gdi/e/ERA5//single/era5_single_2015_06_17.nc', \
    #                                   '/data/gdi/e/ERA5//multi/era5_profile_2015_06_17.nc', \
    #                                     '/data/gdi/e/ERA5/single/era5_single_2015_06_18.nc', \
    #                                         '/data/gdi/e/ERA5/multi/era5_profile_2015_06_18.nc'] 
    
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    # era = ERA5Processor(inputfiles,logger = logger)
    # results = era.era5_process()
    mm = MainProcessor(inputfiles,logger = logger, save_flag='non_debug', output_dir="/data/keeling/a/gzhao1/f/mmcth/test_output/")
    mm.run_process(save_flag = 'non_debug')