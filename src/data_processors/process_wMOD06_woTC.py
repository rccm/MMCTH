from .process_main import MODISMISRProcessor,ERA5Processor
import os
from .netcdf_saver import NetCDFSaver
import numpy as np
import re
class WMOD06WoTC:
    '''Process MODIS '''
    def __init__(self, inputfile_list, logger):
        if len(inputfile_list) != 9:
            logger.error("Error! The input file list is WRONG")
            os.exit()
        self.logger = logger
        self.timestamp = re.search(r'\.A\d{7}\.(\d{4})\.', inputfile_list[0]).group(1)
        self.mm_processor = MODISMISRProcessor(inputfile_list, logger)
        self.era5_processor = ERA5Processor(inputfile_list, logger)
        output_dir: str = "/data/keeling/a/gzhao1/f/mmcth/output/",
        self.output_dir = output_dir
    def run_process(self):
        try:
            self.logger.debug("Running the MISR/MODIS data processing pipeline")
            res = self.mm_processor.mm_process(process_misr_cth=False)
            if res is None:
                self.logger.error("Processing MISR/MODIS failed")
                return
            #!  Finish the part results only return two datasets 
            bands_BT = res['bands_BT']   # may be None
            misr_cth = res['misr_cth']   # may be None
            mod_geo  = res['mod_geo']    # always present
            mod06    = res['mod06']      # may 
    
            era5_lats, era5_lons = self.era5_processor.era5_lat_lon()
            era5_variables = self.era5_processor.era5_process()
            self.logger.debug("ERA5 processing completed successfully")
            era5_variables_misrswath = self.era5_processor.map_era5_to_modis(mod_geo['lat'], mod_geo['lon'], era5_lats, era5_lons, era5_variables)
            output_dir = '/data/keeling/a/gzhao1/f/mmcth/output/'
            outputfile_name = f'{output_dir}MODMISR_L2_CP_{self.mm_processor.id}_O{self.mm_processor.orbit}_v01.nc'

            mm_variables = np.copy(mod06)
            mm_variables['cflag'] = np.zeros_like(mod06['cth']) + self.cflag_modis
            self.save_pixels(mm_variables,mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
            mmcth = NetCDFSaver(outputfile_name, self.logger)
            mmcth.save(None,None,mod_geo,mod06,era5_variables_misrswath)
            # self.logger.debug("Output file saved successfully")           
        except Exception as e:
            self.logger.error(f" Main Processing Problems: {e}")
            return None


 