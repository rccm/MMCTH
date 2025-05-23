from .process_main import MODISMISRProcessor,ERA5Processor
import os
from .netcdf_saver import NetCDFSaver

class WoMOD06WTC:
    '''Process MODIS '''
    def __init__(self, inputfile_list, logger):
        if len(inputfile_list) != 9:
            logger.error("Error! The input file list is WRONG")
            os.exit()
        self.logger = logger
        self.timestamp = re.search(r'\.A\d{7}\.(\d{4})\.', inputfile_list[0]).group(1)
        self.mm_processor = MODISMISRProcessor(inputfile_list, logger)
        self.era5_processor = ERA5Processor(inputfile_list, logger)
    def run_process(self):
        try:
            self.logger.debug("Running the MISR/MODIS data processing pipeline")
            misr_cth, mod_geo  = self.mm_processor.mm_process(process_mod06_cth=False)
            era5_lats, era5_lons = self.era5_processor.era5_lat_lon()
            era5_variables = self.era5_processor.era5_process()
            self.logger.debug("ERA5 processing completed successfully")
            era5_variables_misrswath = self.era5_processor.map_era5_to_modis(mod_geo['lat'], mod_geo['lon'], era5_lats, era5_lons, era5_variables)
            output_dir = '/data/keeling/a/gzhao1/f/mmcth/output/'
            outputfile_name = f'{output_dir}MODMISR_L2_CP_{self.mm_processor.id}_O{self.mm_processor.orbit}_v01.nc'
            mmcth = NetCDFSaver(outputfile_name, self.logger)
            mmcth.save(None,misr_cth,mod_geo,None,era5_variables_misrswath)

            # self.logger.debug("Output file saved successfully")           
        except Exception as e:
            self.logger.error(f" Main Processing Problems: {e}")
            return None
