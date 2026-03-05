import os,sys
import importlib
import re
import logging 
import numpy as np
import scipy
from scipy.interpolate import interp1d,PchipInterpolator
import pandas as pd
import xarray as xr
from pyresample import kd_tree, geometry
from pyresample.geometry import SwathDefinition
from typing import Optional  
from pyproj import Transformer
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]   # mmcth/
sys.path.insert(0, str(PROJECT_DIR))



class MODISMISRProcessor:
    '''  Process MODIS data  and MISR TC_Cloud data if avaible    '''
    def __init__(self, input_files, logger: Optional[logging.Logger] = None):
        if len(input_files) != 9:
            logger.error("Error! The input file list is WRONG")
            sys.exit()    
        self.mod21_file = input_files[0]
        self.mod06_file = input_files[1]
        self.mod03_file = input_files[2]
        self.tccloud_file = input_files[3]
        self.agp_file = input_files[4]
        self.log  = logger or logging.getLogger(self.__class__.__name__)
        self.mod_read = importlib.import_module('src.data_readers.mod_read')
        self.misr_read = importlib.import_module('src.data_readers.misr_read')
        self.id = re.search(r'A\d{7}\.\d{4}', self.mod21_file).group() if re.search(r'A\d{7}\.\d{4}', self.mod21_file) else None
        self.orbit = re.search(r'O\d{6}', self.tccloud_file).group() if re.search(r'O\d{6}', self.tccloud_file) else None
    
    def process_mod21(self, save_format='normal', bands = [1,4, 3]):
        try:
            self.log.debug("Processing MOD21 data")
            mod21 = self.mod_read.MODL1Granule(self.mod21_file)
            bands_BT = {f'bt_{band}': mod21.get_BT(str(band)) for band in [31,32]}
            if (save_format == 'org'):
                print('write reflecance as well...')
                bands_rad = {f'rad_{band}': mod21.get_band_data(str(band), 'radiance') for band in [2,36, 35, 34, 33, 31]}
            else:
                bands_rgb = {f'rad_{band}': mod21.get_band_data(str(band), 'radiance') for band in bands}
            self.log.debug("MOD21 data processed successfully")
            combined_bands = {**bands_BT, **bands_rgb}
            return combined_bands
        except Exception as e:
            self.log.error(f"Error processing MOD21: {self.mod21_file} {e}")
            return None

    def process_mod06(self):
        try:
            self.log.debug("Processing MOD06 data")
            mod06 = self.mod_read.MOD06Granule(self.mod06_file)
            result = {
                'ctp': mod06.get_ctp(),
                'opt': mod06.get_opt(),
                'cphase': mod06.get_cphase(),
                'cth': mod06.get_cth(),
                'emissivity': mod06.get_emissivity(),
                'ctm': mod06.get_ctm(),
            }
            self.log.debug("MOD06 data processed successfully")
            return result
        except Exception as e:
            self.log.error(f"Error processing MOD06: {self.mod06_file} {e}")
            return None
        
    def process_mod03(self):
        try:
            self.log.debug("Processing MOD03 data")
            mod03 = self.mod_read.MOD03Granule(self.mod03_file)
            result = mod03.get_lat(), mod03.get_lon(), mod03.get_landsea_mask(), mod03.get_vza()
            self.log.debug("MOD03 data processed successfully")
            return result
        except Exception as e:
            self.log.error(f"Error processing MOD03: {self.mod03_file} {e}")
            return None
        
    def get_misr_geo(self):
        try:
            agp = self.misr_read.MISRGranule(self.agp_file)
            misr_lat, misr_lon = agp.get_geo()
            return misr_lat, misr_lon
        except Exception as e:
            self.log.error(f"Error getting MISR geo: {e}")
            return None   
    
    def get_misr_cth(self):
            try:
                tc_cloud = self.misr_read.MISRGranule(self.tccloud_file)
                cth,cth_qa = tc_cloud.get_cth()
                return  cth,cth_qa
            except Exception as e:
                self.log.error(f"Error getting MISR CTH: {e}")
                return None
            
    def misr_to_modis(self, source_lat, source_lon, target_lat, target_lon, proj_data):
        try:
            source_def = SwathDefinition(lons=source_lon, lats=source_lat)
            target_def = SwathDefinition(lons=target_lon, lats=target_lat)
            reproj_data = kd_tree.resample_nearest(
                source_def,
                proj_data,
                target_def,
                radius_of_influence=1100,
                fill_value=-999
            )
            return reproj_data
        except Exception as e:
            self.log.error(f"Error in MISR to MODIS reprojection: {e}")
            return None

    def cutoff_misr_swath(self, misr_swath):
        try:
            misr_swath[misr_swath < -990] = np.nan
            mask = ~np.isnan(misr_swath)
            valid_rows = np.any(mask, axis=1)
            valid_cols = np.any(mask, axis=0)
            return valid_rows, valid_cols
        except Exception as e:
            self.log.error(f"Error in cutoff MISR swath: {e}")
            return None
        
    def apply_valid_indices(self, array, valid_rows, valid_cols):
        try:
            return array[valid_rows][:, valid_cols]
        except Exception as e:
            self.log.error(f"Error applying valid indices: {e}")
            return None
    
    def _to_qi_int8(self,a, *, missing_sentinel=-999.0, fill=-128):
        """Map NaN/±Inf/-999 → -128, round valid values, clip to [0,100], cast to int8."""
        # handle masked arrays from upstream (just in case)
        if np.ma.isMaskedArray(a):
            a = a.filled(np.nan)

        a = np.asarray(a, dtype=np.float32)
        bad = (~np.isfinite(a)) | (a == missing_sentinel)    # NaN/Inf or -999
        a[bad] = float(fill)

        keep = a != float(fill)
        # round and clip valid QI values to [0, 100]
        a[keep] = np.clip(np.rint(a[keep]), 0, 100)

        # guard against overflow before casting
        np.clip(a, -128, 127, out=a)
        return a.astype(np.int8, copy=False)
    
    def mm_process(self, *, process_misr_cth: bool = True, process_mod06_cth: bool = True):
        """
        Re‑project MISR → MODIS, subset to the overlapping swath, and (optionally)
        assemble MISR‑CTH and MOD06 fields.
        Returns
        -------
        dict
            {
                'bands_BT' : dict | None,   # 5 IR BT + radiance bands on MODIS grid
                'misr_cth' : dict | None,   # {'misrcth', 'misrcth_qa'} on MODIS grid
                'mod_geo'  : dict,          # lat/lon/landsea/vza  (always present)
                'mod06'    : dict | None    # MOD06 variables on MODIS grid
            }
            Any unavailable entry is set to None.
        """
        # ---------- initialise the return container ----------
        out = {'bands': None, 'mod_geo': None, 'mod06': None}

        # ---------- 1. Read MOD03 geo ----------
        mod_lat, mod_lon, mod_landsea, mod_vza = self.process_mod03()
        if any(v is None for v in (mod_lat, mod_lon, mod_landsea)):
            self.log.error("Error retrieving MODIS lat/lon/land‑sea mask")
            return None
        if np.isnan(mod_lat).any() or np.isnan(mod_lon).any() or (mod_lat<-444.0).any() or (mod_lon <-444.0).any() :
            self.log.error("NaN detected in MODIS lat/lon; returning None.")
            return None         
        
        # ---------- 2. Read MISR geo ----------
        misr_lat, misr_lon = self.get_misr_geo()
        # print(misr_lat.shape, mod_lat.shape)
        if any(v is None for v in (misr_lat, misr_lon)):
            self.log.error("Error retrieving MISR lat/lon")
            return None

        # ---------- 3. Overlap mask ----------
        misr_mod_swath = self.misr_to_modis(misr_lat, misr_lon, mod_lat, mod_lon, misr_lat)
        if misr_mod_swath is None:
            self.log.error("Error finding MISR swath on MODIS grid")
            return None    
        valid_rows, valid_cols = self.cutoff_misr_swath(misr_mod_swath)
        if any(v is None for v in (valid_rows, valid_cols)):
            self.log.error("Error retrieving MISR lat/lon")
            return None 
        
        # Helper to slice MODIS arrays to the MISR footprint
        s = lambda arr: self.apply_valid_indices(arr, valid_rows, valid_cols)

        # ---------- 4. Always build MODIS geo subset ----------
        out['mod_geo'] = {
            'lat'     : s(mod_lat),
            'lon'     : s(mod_lon),
            'landsea' : s(mod_landsea),
            'vza'     : s(mod_vza),
        }

        # ---------- 5. Bands BT/radiance (needed only if MOD06 is used) ----------
        if process_mod06_cth:
            bands_rad = self.process_mod21()
            if bands_rad is None:
                self.log.error("Error reading MOD21")
                return None
            out['bands'] = {k: s(v) for k, v in bands_rad.items()}
            mod06 = self.process_mod06()
            if mod06 is None:
                self.log.error("Error reading MOD06")
                return None
            out['mod06'] = {k: s(v) for k, v in mod06.items()}
        else:
            out['bands_rad'] = None
            out['mod06'] = None
        return out
            