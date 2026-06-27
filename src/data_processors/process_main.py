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
from .netcdf_saver import NetCDFSaver
from .xarray_saver import XarraySaver 
from src.pixel import modis
from pyproj import Transformer

from .misr_cth_to_pressure import height_to_log_pressure
from .interpolate_to_pressure_levels import interpolate_to_pressure_levels

class MODISMISRProcessor:
    '''  Process MODIS data  and MISR TC_Cloud data if avaible    '''
    def __init__(
        self,
        input_files,
        logger: Optional[logging.Logger] = None,
        destripe_flag: bool = True,
    ):
        if len(input_files) != 9:
            logger.error("Error! The input file list is WRONG")
            sys.exit()    
        self.mod21_file = input_files[0]
        self.mod06_file = input_files[1]
        self.mod03_file = input_files[2]
        self.tccloud_file = input_files[3]
        self.agp_file = input_files[4]
        self.log  = logger or logging.getLogger(self.__class__.__name__)
        self.destripe_flag = destripe_flag
        self.mod_read = importlib.import_module('src.data_readers.mod_read')
        self.misr_read = importlib.import_module('src.data_readers.misr_read')
        self.id = re.search(r'A\d{7}\.\d{4}', self.mod21_file).group() if re.search(r'A\d{7}\.\d{4}', self.mod21_file) else None
        self.orbit = re.search(r'O\d{6}', self.tccloud_file).group() if re.search(r'O\d{6}', self.tccloud_file) else None
    
    def process_mod21(self, save_format='normal'):
        try:
            self.log.debug("Processing MOD21 data")
            mod21 = self.mod_read.MODL1Granule(
                self.mod21_file,
                destripe_flag=self.destripe_flag,
            )
            bands_BT = {f'bt_{band}': mod21.get_BT(str(band)) for band in [36, 35, 34, 33, 31]}
            if (save_format == 'org'):
                print('write reflecance as well...')
                bands_radiance = {f'rad_{band}': mod21.get_band_data(str(band), 'radiance') for band in [2,36, 35, 34, 33, 31]}
                combined_bands = {**bands_BT, **bands_radiance}
            elif (save_format == 'debug'):
                bands_radiance = {f'rad_{band}': mod21.get_band_data(str(band), 'radiance') for band in [2,36, 35, 34, 33, 31]}
                bands_dn = {f'dn_{band}': mod21.get_band_data(str(band), 'radiance',scale_flag=False) for band in [2,36, 35, 34, 33, 31]} 
                combined_bands = {**bands_BT, **bands_radiance,**bands_dn}
            else:
                bands_radiance = {f'rad_{band}': mod21.get_band_data(str(band), 'radiance') for band in [36, 35, 34, 33, 31]}
                combined_bands = {**bands_BT, **bands_radiance}

            self.log.debug("MOD21 data processed successfully")
           
            # combined_bands = {**bands_BT}
            return combined_bands
        except Exception as e:
            self.log.error(f"Error processing MOD21: {self.mod21_file} {e}")
            return None

    def process_mod06(self):
        try:
            self.log.debug("Processing MOD06 data")
            mod06 = self.mod_read.MOD06Granule(self.mod06_file)
            qa = mod06.get_qa()   # shape (ny, nx, 9), int8 or byte
            qa_u8 = np.asarray(qa, dtype=np.int8).view(np.uint8)
            byte4 = qa_u8[:, :, 4]
            multilayer_flag = ((byte4 >> 3) & 0b111).astype(np.int8)
            # byte5 = qa_u8[:, :, 5]
            multilayer_flag = ((byte4 >> 3) & 0b111).astype(np.int8)
            # ml_test_phase_diff = ((byte5 >> 0) & 0b1).astype(np.int8)
            # ml_test_delta_pwv = ((byte5 >> 1) & 0b1).astype(np.int8)
            # ml_test_delta_pwv_900 = ((byte5 >> 2) & 0b1).astype(np.int8)
            # ml_test_tau_diff = ((byte5 >> 3) & 0b1).astype(np.int8)
            # ml_test_pavolonis = ((byte5 >> 4) & 0b1).astype(np.int8)

            result = {
                'ctp': mod06.get_ctp(),
                'opt': mod06.get_opt(),
                'cphase': mod06.get_cphase(),
                'cth': mod06.get_cth(),
                'emissivity': mod06.get_emissivity(),
                'ctm': mod06.get_ctm(),
                # 'cmulti': mod06.get_multi(),
                # 'cqa': multilayer_flag,
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
            result = (
                mod03.get_lat(),
                mod03.get_lon(),
                mod03.get_landsea_mask(),
                mod03.get_vza(),
                mod03.get_time_2d_since_2000(use_center_time=True, fast=True),
                mod03.get_time_1d_since_2000(use_center_time=True),
            )
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
                sdcm = tc_cloud.get_sdcm()
                agp = self.misr_read.MISRGranule(self.agp_file)
                landwater = agp.get_landmask()
                land_mask = np.isin(landwater,[1,2,3,4])
                sdcm_mask = np.isin(sdcm, [3,4])
                cth[land_mask & sdcm_mask] = -9999.0
                return  cth,cth_qa
            except Exception as e:
                self.log.error(f"Error getting MISR CTH: {e}")
                return None
            
    def misr_to_modis(self, source_lat, source_lon, target_lat, target_lon, proj_data):
        try:
            source_lat = self._as_float_array(source_lat)
            source_lon = self._as_float_array(source_lon)
            target_lat = self._as_float_array(target_lat)
            target_lon = self._as_float_array(target_lon)
            source_valid = self._valid_latlon_mask(source_lat, source_lon)
            target_valid = self._valid_latlon_mask(target_lat, target_lon)

            if not np.any(source_valid):
                self.log.error("No valid source geolocation for MISR to MODIS reprojection")
                return None
            if not np.any(target_valid):
                self.log.error("No valid target MODIS geolocation for MISR to MODIS reprojection")
                return None

            proj_data = np.asarray(proj_data)
            if proj_data.shape != source_lat.shape:
                self.log.error(
                    "Projection data shape %s does not match source geolocation shape %s",
                    proj_data.shape,
                    source_lat.shape,
                )
                return None

            target_lat_safe = np.where(target_valid, target_lat, 0.0)
            target_lon_safe = np.where(target_valid, target_lon, 0.0)
            source_values = proj_data[source_valid]
            if np.issubdtype(source_values.dtype, np.integer):
                dtype_info = np.iinfo(source_values.dtype)
                if -999 < dtype_info.min or -999 > dtype_info.max:
                    source_values = source_values.astype(np.int16, copy=False)

            source_def = SwathDefinition(lons=source_lon[source_valid], lats=source_lat[source_valid])
            target_def = SwathDefinition(lons=target_lon_safe, lats=target_lat_safe)
            reproj_data = kd_tree.resample_nearest(
                source_def,
                source_values,
                target_def,
                radius_of_influence=1100,
                fill_value=-999
            )
            if np.ma.isMaskedArray(reproj_data):
                reproj_data = reproj_data.filled(-999)
            reproj_data = np.asarray(reproj_data)
            reproj_data[~target_valid] = -999
            return reproj_data
        except Exception as e:
            self.log.error(f"Error in MISR to MODIS reprojection: {e}")
            return None

    @staticmethod
    def _as_float_array(arr):
        if np.ma.isMaskedArray(arr):
            arr = arr.filled(np.nan)
        return np.asarray(arr, dtype=np.float64)

    def _valid_latlon_mask(self, lat, lon):
        lat = self._as_float_array(lat)
        lon = self._as_float_array(lon)
        return (
            np.isfinite(lat) &
            np.isfinite(lon) &
            (lat >= -90.0) &
            (lat <= 90.0) &
            (lon >= -180.0) &
            (lon <= 180.0)
        )

    def _prune_invalid_geo_from_rect(self, valid_rows, valid_cols, geo_valid):
        rows = np.asarray(valid_rows, dtype=bool).copy()
        cols = np.asarray(valid_cols, dtype=bool).copy()

        while np.any(rows) and np.any(cols):
            selected_rows = np.flatnonzero(rows)
            selected_cols = np.flatnonzero(cols)
            invalid = ~geo_valid[np.ix_(selected_rows, selected_cols)]
            if not np.any(invalid):
                return rows, cols

            bad_rows = np.any(invalid, axis=1)
            bad_cols = np.any(invalid, axis=0)
            if np.count_nonzero(bad_rows) <= np.count_nonzero(bad_cols):
                rows[selected_rows[bad_rows]] = False
            else:
                cols[selected_cols[bad_cols]] = False

        return rows, cols

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
    
    def mm_process(self, *, process_misr_cth: bool = True, process_mod06_cth: bool = True, scale_flag: str = 'not_debug'):
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
        out = {'bands_BT': None, 'misr_cth': None, 'mod_geo': None, 'mod06': None}

        # ---------- 1. Read MOD03 geo ----------
        mod_lat, mod_lon, mod_landsea, mod_vza,time2d,time1d = self.process_mod03()
        if any(v is None for v in (mod_lat, mod_lon, mod_landsea)):
            self.log.error("Error retrieving MODIS lat/lon/land‑sea mask")
            return None
        mod_geo_valid = self._valid_latlon_mask(mod_lat, mod_lon)
        if not np.any(mod_geo_valid):
            self.log.error("No valid MODIS lat/lon values; returning None.")
            return None
        invalid_geo_count = int(mod_geo_valid.size - np.count_nonzero(mod_geo_valid))
        if invalid_geo_count:
            self.log.warning(
                "Ignoring %d invalid MODIS geolocation pixels before reprojection.",
                invalid_geo_count,
            )
        
        # ---------- 2. Read MISR geo ----------
        misr_lat, misr_lon = self.get_misr_geo()
        if any(v is None for v in (misr_lat, misr_lon)):
            self.log.error("Error retrieving MISR lat/lon")
            return None

        # ---------- 3. Overlap mask ----------
        misr_mod_swath = self.misr_to_modis(misr_lat, misr_lon, mod_lat, mod_lon, misr_lat)
        if misr_mod_swath is None or np.all(misr_mod_swath <= -999):
            print('return')
            self.log.error("Error finding MISR swath on MODIS grid")
            return None
        valid_rows, valid_cols = self.cutoff_misr_swath(misr_mod_swath)
        if any(v is None for v in (valid_rows, valid_cols)):
            self.log.error("Error retrieving MISR lat/lon")
            return None 
        rows_before = np.count_nonzero(valid_rows)
        cols_before = np.count_nonzero(valid_cols)
        valid_rows, valid_cols = self._prune_invalid_geo_from_rect(
            valid_rows,
            valid_cols,
            mod_geo_valid,
        )
        if not np.any(valid_rows) or not np.any(valid_cols):
            self.log.error("No valid MODIS/MISR overlap remains after removing invalid MODIS geolocation.")
            return None
        rows_after = np.count_nonzero(valid_rows)
        cols_after = np.count_nonzero(valid_cols)
        if rows_after != rows_before or cols_after != cols_before:
            self.log.warning(
                "Removed %d rows and %d columns from the MODIS/MISR overlap because of invalid MODIS geolocation.",
                rows_before - rows_after,
                cols_before - cols_after,
            )
        # Helper to slice MODIS arrays to the MISR footprintprint(valid_rows)
        s = lambda arr: self.apply_valid_indices(arr, valid_rows, valid_cols)

        # ---------- 4. Always build MODIS geo subset ----------
        out['mod_geo'] = {
            'lat'     : s(mod_lat),
            'lon'     : s(mod_lon),
            'landsea' : s(mod_landsea),
            'vza'     : s(mod_vza),
            'time2d'  : s(time2d),
            'time1d'  : time1d[valid_rows],
            'orig_row_index' : np.flatnonzero(valid_rows).astype(np.int32),
            'orig_col_index' : np.flatnonzero(valid_cols).astype(np.int32),
            'orig_nrows'     : np.int32(mod_lat.shape[0]),
            'orig_ncols'     : np.int32(mod_lat.shape[1]),
        }
        # ---------- 5. Bands BT/radiance (needed only if MOD06 is used) ----------
        if process_mod06_cth:
            bands_BT = self.process_mod21(save_format=scale_flag)
            if bands_BT is None:
                self.log.error("Error reading MOD21")
                return None
            out['bands_BT'] = {k: s(v) for k, v in bands_BT.items()}
            mod06 = self.process_mod06()
            if mod06 is None:
                self.log.error("Error reading MOD06")
                return None
            out['mod06'] = {k: s(v) for k, v in mod06.items()}
        else:
            out['bands_BT'] = None
            out['mod06'] = None

        # ---------- 7. MISR CTH subset ----------'
        if process_misr_cth:
            misr_cth, misr_cth_qa = self.get_misr_cth()
            if misr_cth is None:
                self.log.error("Error reading MISR CTH")
                return None
            
        # Converted Ellipsoid Heights to Sea Surface Height
            
            # Build a pipeline: ellipsoidal (WGS84) -> EGM2008 orthometric height via 1′ grid
   
            MISR_CTH_FILL = -9999.0

            # EGM2008 1′ grid with 2.5′ fallback (make sure the files are in PROJ's data dir)
            pipeline = (
                "+proj=pipeline "
                "+step +proj=longlat +datum=WGS84 "
                "+step +proj=vgridshift +grids=us_nga_egm2008_1.tif,egm08_25.gtx"
            )
            tr = Transformer.from_pipeline(pipeline)
            # Flatten to consistent 1-D views
            lonf = np.asarray(misr_lon, dtype=np.float64).ravel(order="C")
            latf = np.asarray(misr_lat, dtype=np.float64).ravel(order="C")
            hf   = np.asarray(misr_cth, dtype=np.float64).ravel(order="C")
            
            # Valid mask: real coords and not fill
            valid = (
                np.isfinite(hf) & (hf != MISR_CTH_FILL) &
                np.isfinite(latf) & (latf >= -90.0) & (latf <= 90.0) &
                np.isfinite(lonf) & (lonf >= -180.0) & (lonf <= 180.0)
            )
            # Output initialized with FILL
            h_msl_f = np.full(hf.shape, MISR_CTH_FILL, dtype=np.float64)
            # Transform only valid points (H = h - N)
            if valid.any():
                _, _, h_valid = tr.transform(lonf[valid], latf[valid], hf[valid])
                h_msl_f[valid] = h_valid

            misr_cth = h_msl_f.reshape(misr_lon.shape, order="C")
            qi_data = self.misr_to_modis(misr_lat, misr_lon, out['mod_geo']['lat'], out['mod_geo']['lon'], misr_cth_qa) 
            misrcth_data = self.misr_to_modis(misr_lat, misr_lon, out['mod_geo']['lat'], out['mod_geo']['lon'], misr_cth)
            if qi_data is None or misrcth_data is None:
                self.log.error("Error resampling MISR CTH/QA to MODIS grid")
                return None
            qi_data_int8 = self._to_qi_int8(qi_data)  
            out['misr_cth'] = {
                'misrcth'    : misrcth_data,
                'misrcth_qa' : qi_data_int8
            }
            del lonf, latf, hf,h_msl_f
        else:
            out['misr_cth'] = {
                'misrcth'    : None,
                'misrcth_qa' : None,
            }

        return out
            

class ERA5Processor:
    def __init__(self, input_files, logger: Optional[logging.Logger] = None, latlon_idx = None):
        self.log  = logger or logging.getLogger(self.__class__.__name__)
        self.era5_single_file = input_files[5]
        self.era5_multi_file = input_files[6]
        self.era5_single_file_next_day = input_files[7]
        self.era5_multi_file_next_day = input_files[8]
        self.mod21_file  = input_files[0]
        self.era5_read = importlib.import_module('src.data_readers.era5_read')
        self.month = int(re.search(r"_\d{4}_(\d{2})_\d{2}\.nc", self.era5_multi_file).group(1))
        self.MODIS = importlib.import_module('src.pixel.modis')
        self.P_levels = self.MODIS.transmission.pstd
        self.P_levels_gt_1mb = self.P_levels[self.P_levels>=1]
        self.era5_P_levels = np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000], dtype=int)
        self.latlon_idx = latlon_idx
        self.timestamp = re.search(r"\.(\d{4})\.",self.mod21_file).group(1) 
        self.profile_file = '/data/gdi/f/gzhao1/mmcth/externaldata/std_atmos_profile.nc'

    def era5_timeinterpolate_bk(self, era5_vars_dict, era5_vars_dict_next_day=None):
        try:
            self.log.debug("Interpolating ERA5 data to target time")
            interpolated_data = {}
            target_hour = float(self.timestamp[0:2])
            target_minute = float(self.timestamp[2:])
            target_time = target_hour + target_minute / 60.0
            
            # Determine if next day's data is needed
            next_day = target_time > 23.0
            
            for key, era5_var in era5_vars_dict.items():
                era5_var = np.array(era5_var)
                hr_interval = 24 / era5_var.shape[0]
                before_idx = int(target_time // hr_interval)
                after_idx = (before_idx + 1) % era5_var.shape[0]
                weight_after = (target_time - before_idx * hr_interval) / hr_interval
                weight_before = 1 - weight_after
                
                # If next day's data is required, use it for interpolation
                if next_day and era5_vars_dict_next_day:
                    next_day_var = np.array(era5_vars_dict_next_day[key])
                    interpolated = era5_var[-1] * weight_before + next_day_var[0] * weight_after
                else:
                    interpolated = era5_var[before_idx] * weight_before + era5_var[after_idx] * weight_after
                interpolated_data[key] = interpolated
            self.log.debug("ERA5 data interpolated to target time successfully")
            return interpolated_data
        except Exception as e:
            self.log.error(f"Error interpolating ERA5 data: {e}")
            return None
    
    def era5_lat_lon(self):
        try:
            self.log.debug("Retrieving ERA5 lat/lon grid")
            era5_single = self.era5_read.ERA5Single(self.era5_single_file)
            result = era5_single.get_latitude(), era5_single.get_longitude()
            self.log.debug("ERA5 lat/lon grid retrieved successfully")
            return result
        except Exception as e:
            self.log.error(f"Error retrieving ERA5 lat/lon grid: {e}")
            return None

    def ear5_read_single(self,era5_single_file):
        ''' Read ERA5 surface/near-surface variables'''
        try:
            self.log.debug("Processing ERA5 data (surface layer)")
            era5_single = self.era5_read.ERA5Single(era5_single_file)

            dew2m = era5_single.get_d2m()
            sp = era5_single.get_sp()
            # convert dew2m to mixing ratio
            dew2m_mixingratio = 622.0 *6.113e2 * np.exp(5423.0 * (dew2m - 273.15)/(dew2m * 273.15))/(1000.0 * sp)
            result = {
                'swmr': dew2m_mixingratio,
                'sst': era5_single.get_sst(),
                'temp2m': era5_single.get_t2m(),
                'surface_pressure': era5_single.get_sp(),
                'skint': era5_single.get_skt(),  
                'msp' : era5_single.get_msl(),      
            }
            self.log.debug("ERA5 data (surface layer) read successfully")
            return result
        except Exception as e:
            self.log.error(f"Error processing ERA5 data (surface layer): {e}")
            return None

    def ear5_read_multi(self,era5_multi_file):
        ''' Read ERA5 pressure-slevel variables'''
        try:
            self.log.debug("Processing ERA5 data (pressure level)")
            era5_multi = self.era5_read.ERA5Multi(era5_multi_file)
            result = {
                'temperature': era5_multi.get_temperature(),
                'specific_humidity': era5_multi.get_specific_humidity(),  # convert the unit of kg/kg to g/kg
                'geopotential': era5_multi.get_geopotential(),
                'u': era5_multi.get_u(),
                'v': era5_multi.get_v(),
                'w': era5_multi.get_w(),
            }
            self.log.debug("ERA5 data (pressure level) read successfully")
            return result
        except Exception as e:
            self.log.error(f"Error processing ERA5 data (pressure level): {e}")
            return None

    def save_to_netcdf(self, filename, **variables):
        """
        Save variables to a NetCDF file using xarray.

        Parameters:
            filename (str): The name of the NetCDF file.
            **variables: The variables to save, provided as keyword arguments.
        """
        latitude, longitude = self.era5_lat_lon()

        data_vars = {}
        for name, data in variables.items():
            if len(data.shape) == 3:
                num_levels = data.shape[0]
                coords = {'latitude': latitude, 'longitude': longitude, f'level_{name}': np.arange(num_levels)}
                data_vars[name] = (['level_' + name, 'latitude', 'longitude'], data)
            elif len(data.shape) == 2:
                coords = {'latitude': latitude, 'longitude': longitude}
                data_vars[name] = (['latitude', 'longitude'], data)
        ds = xr.Dataset(data_vars, coords=coords)
        ds.to_netcdf(filename)
        self.log.info(f"Saved data to {filename}")
        return
    
    def convert_lon(self,lon):
        """Convert longitude from -180-180 to 0-360 range if necessary"""
        return np.where(lon < 0, lon + 360, lon)

    def find_nearest_indices(self,modis_lats, modis_lons, era5_lats, era5_lons):
        """ 
        Find the nearest ERA5 grid indices for each MODIS pixel.
        
        Parameters:
        modis_lats (np.ndarray): MODIS latitude array
        modis_lons (np.ndarray): MODIS longitude array
        era5_lats (np.ndarray): ERA5 latitude grid
        era5_lons (np.ndarray): ERA5 longitude grid
        
        Returns: 
        np.ndarray, np.ndarray: Indices of the nearest ERA5 grid points
        """
        # Convert MODIS longitudes to 0-360 range
        modis_lons = self.convert_lon(modis_lons)
        
        # Round MODIS latitudes and longitudes to the nearest 0.25 degrees
        rounded_modis_lats = np.round(modis_lats / 0.25) * 0.25
        rounded_modis_lons = np.round(modis_lons / 0.25) * 0.25
        
        # Find the indices of the rounded values in the ERA5 latitude and longitude arrays
        # lat_indices = np.searchsorted(era5_lats, rounded_modis_lats) - 1
        era5_lats_sorted = era5_lats[::-1]
        lat_indices = len(era5_lats) - np.searchsorted(era5_lats_sorted, rounded_modis_lats, side='right')
        lon_indices = np.searchsorted(era5_lons, rounded_modis_lons) - 1
        
        return lat_indices, lon_indices
    
    def map_era5_to_modis(self, modis_lats, modis_lons, era5_lats, era5_lons, era5_vars_dict):
        try:
            self.log.debug("Mapping ERA5 variables to MODIS grid points")
            lat_indices, lon_indices = self.find_nearest_indices(modis_lats, modis_lons, era5_lats, era5_lons)

            modis_vars = {}
            for key, era5_var in era5_vars_dict.items():
                if len(era5_var.shape) == 3:
                    interpolated_var = np.empty((era5_var.shape[0], modis_lats.shape[0], modis_lats.shape[1]))
                    for level in range(era5_var.shape[0]):
                        interpolated_var[level, :, :] = era5_var[level, lat_indices, lon_indices]
                    modis_vars[key] = interpolated_var
                else:
                    modis_var = era5_var[lat_indices, lon_indices]
                    modis_vars[key] = modis_var
            self.log.debug("ERA5 variables mapped to MODIS grid points successfully")
            return modis_vars
        except Exception as e:
            self.log.error(f"Error mapping ERA5 to MODIS: {e}")
            return None

    def load_std_atmos_profile(self):
        netcdf_file = self.profile_file
        month = self.month
        ds = xr.open_dataset(netcdf_file)
        if month in (4, 5, 6, 7, 8, 9):
            dataset = ds['summer_spring']
        else:
            dataset = ds['fall_winter']
        
        return dataset


    def era5_process(self, no_interploate = False):
        """
        Interpolates multi-level ERA5 data to 101 pressure .
        
        Parameters:
            profiles (optional): Profile data for determining profiles based on latitude.
        
        Return:

            turple: EAR5 surface/near-surface data along with profile data that have been interpolated into 101 pressure levels 
        """
        

        latitudes, _ = self.era5_lat_lon()
        latitudes = latitudes.values
        
        def _normalize_time_dim(d: dict, *, ntime=24):
            out = {}
            for k, v in d.items():
                if not hasattr(v, "dims"):
                    out[k] = v
                    continue

                rename = {}

                # normalize time dim
                for cand in ("time", "valid_time"):
                    if cand in v.dims and v.sizes.get(cand, None) == ntime:
                        if cand != "time":
                            rename[cand] = "time"
                        break

                # normalize vertical dim (critical for next-day mixing)
                if "pressure_level" in v.dims and "level" not in v.dims:
                    rename["pressure_level"] = "level"
                # some cfgrib files use this name
                if "isobaricInhPa" in v.dims and "level" not in v.dims:
                    rename["isobaricInhPa"] = "level"

                if rename:
                    v = v.rename(rename)

                out[k] = v
            return out
        # def _try_load(d, label, tindex):
        #     for k, v in d.items():
        #         try:
        #             src = getattr(v, "encoding", {}).get("source", "NA")
        #             print(f"[{label}] trying {k} time={tindex} source={src}")
        #             v.isel(time=tindex).load()   # forces read
        #             print(f"[{label}] OK {k}")
        #         except Exception as e:
        #             print(f"[{label}] FAIL {k} time={tindex} source={src}\n  -> {repr(e)}")
        #             raise
        var_single = _normalize_time_dim(self.ear5_read_single(self.era5_single_file), ntime=24)
        var_multi  = _normalize_time_dim(self.ear5_read_multi(self.era5_multi_file),  ntime=24)
        assert var_multi["temperature"].sizes["time"] == 24  
        # --- target time in hours (UTC) ---
        hh = int(self.timestamp[:2])
        mm = int(self.timestamp[2:])
        w  = mm / 60.0          # weight on next hour
        wb = 1.0 - w            # weight on current hour
        before = hh
        after  = (hh + 1) % 24
        need_next_day = (hh == 23 and mm > 0)
        if need_next_day:
            var_single_next_day = _normalize_time_dim(self.ear5_read_single(self.era5_single_file_next_day), ntime=24)
            var_multi_next_day  = _normalize_time_dim(self.ear5_read_multi(self.era5_multi_file_next_day),  ntime=24)
            var_single = {k: (v.isel(time=before) * wb + var_single_next_day[k].isel(time=0) * w) for k, v in var_single.items()}
            var_multi  = {k: (v.isel(time=before) * wb + var_multi_next_day[k].isel(time=0)  * w) for k, v in var_multi.items()}
        else:
            var_single = {k: (v.isel(time=before) * wb + v.isel(time=after) * w) for k, v in var_single.items()}
            var_multi  = {k: (v.isel(time=before) * wb + v.isel(time=after) * w) for k, v in var_multi.items()}
        # unit conversion after interpolation (kg/kg -> g/kg)
        if no_interploate: 
            result = var_single.copy()
            result.update(var_multi)
            return result
        var_multi["specific_humidity"] = var_multi["specific_humidity"] * 1e3
        T_multi = var_multi['temperature'] 
        W_multi = var_multi['specific_humidity']/1000/(1-var_multi['specific_humidity']/1000) # check again .....
        Z_multi = var_multi['geopotential'] / 9.80665e3
        self.log.debug("ERA5 Temporal Interoplation is done successfully") 
         
        # Define dimensions explicitly
        num_era5_levels = int(T_multi.shape[0])  # 37
        num_modis_levels = int(self.P_levels_gt_1mb.shape[0])  # 91
        lat_size = int(T_multi.shape[1])  # 721
        lon_size = int(T_multi.shape[2])  # 1440
        
        
        era5_levels_1d = np.asfortranarray(np.asarray(self.era5_P_levels, dtype=np.float32).reshape(-1))
        modis_levels_1d = np.asfortranarray(np.asarray(self.P_levels_gt_1mb, dtype=np.float32).reshape(-1))

        # Prepare 3D arrays and ensure they are Fortran contiguous
        T_multi = np.asfortranarray(T_multi.astype(np.float32))
        W_multi = np.asfortranarray(W_multi.astype(np.float32))
        Z_multi = np.asfortranarray(Z_multi.astype(np.float32))
 
        # Call the Fortran subroutine with explicitly shaped arrays

        surface_pressures = np.asfortranarray((var_single['surface_pressure'].astype(np.float32)) / 100.0)  # Convert to hPa
        W_surface = np.asfortranarray(var_single['swmr'].astype(np.float32))
        T_surface_air = np.asfortranarray(var_single['temp2m'].astype(np.float32))
        T_sst = np.asfortranarray(var_single['sst'].astype(np.float32))

    
        num_modis_levels = np.int32(modis_levels_1d.shape[0])
        lat_size = np.int32(T_multi.shape[1])
        lon_size = np.int32(T_multi.shape[2])
         
        T_modis = np.empty((num_modis_levels, lat_size, lon_size), dtype=np.float32, order='F')
        W_modis = np.empty((num_modis_levels, lat_size, lon_size), dtype=np.float32, order='F')
        Z_modis = np.empty((num_modis_levels, lat_size, lon_size), dtype=np.float32, order='F')

        interpolate_to_pressure_levels(
            era5_levels_1d, modis_levels_1d, T_multi, W_multi, Z_multi,
            surface_pressures, W_surface, T_surface_air, T_sst, T_modis, W_modis, Z_modis,
        )
        T_interpolated = np.asfortranarray(T_modis)
        W_interpolated = np.asfortranarray(W_modis)
        Z_interpolated = np.asfortranarray(Z_modis)
            
        # Combine the standard atmospheric profile data for pressure levels smaller than 1mb


        profile_data = self.load_std_atmos_profile()

        temperature_profile = profile_data.sel(variable="T(K)").values
        # Not shifting 
        # offset_T = T_interpolated[0, :, :] - temperature_profile[-1, :, :]
        # temperature_profile += offset_T[np.newaxis, :, :]

        T_expand = np.concatenate((temperature_profile,T_interpolated),axis=0)
        
        # need to verify each variable
        h20_profile = profile_data.sel(variable="h2o(cm-3)").values #unit g/kg water mixing ratio not water mass any more 
        w_profile = h20_profile/1e3   

        # offset_Q = W_interpolated[0, :, :] - sq_profile[-1, :, :]        
        # sq_profile += offset_Q[np.newaxis, :, :] 

        W_expand = np.concatenate((w_profile,W_interpolated),axis=0)
        W_expand = W_expand*1e3

        Z_profile = profile_data.sel(variable="z(km)").values
        offset_Z =  Z_interpolated[0, :, :] - Z_profile[-1, :, :]
        Z_profile += offset_Z[np.newaxis, :, :]
        Z_expand = np.concatenate((Z_profile,Z_interpolated),axis=0)
          
        # Set data below SST

        result = var_single.copy()
        result.update({
            'temperature_exp': T_expand,
            'mixingratio_exp': W_expand,
            'geopotential_exp': Z_expand,
            'temperature_org': T_multi,
            'mixingratio_org': W_multi*1e3,
            'geopotential_org': Z_multi,
            'u': var_multi['u'],
            'v': var_multi['v'],
            'w': var_multi['w'],
        })
        for k, v in result.items():
            if hasattr(v, "values"):
                result[k] = v.values
                #! For Debugging  
        # self.log.info("ERA5 result contains %d variables", len(result))
        # for k, v in result.items():
        #     if hasattr(v, "shape"):
        #         self.log.info("ERA5 %-22s shape=%s dtype=%s", k, v.shape, v.dtype)
        #     else:
        #         self.log.info("ERA5 %-22s type=%s", k, type(v).__name__)

        # internal_variable = {
        #     't_org': T_multi,
        #     'z_org': Z_multi,
        #     'q_org': W_multi*1000.0,
        # }
        # expand_variable = {
        #     'temperature': T_expand,
        #     'mixingratio': W_expand,
        #     'geopotential': Z_expand,
            
        # }
        # self.save_to_netcdf(f'../output/check_org.nc',**internal_variable)
        # self.save_to_netcdf(f'../output/check_int.nc',**expand_variable)
        # self.save_to_netcdf(f'../output/check_single.nc',**var_single)
 
        return result 


class MainProcessor: 

    def __init__(
        self,
        inputfile_list: list[str],
        orbit: str = None,
        logger: Optional[logging.Logger] = None,
        output_dir: str = "ds_output",
        save_flag: str = 'not_debug',
        destripe_flag: bool = True,
    ):
        # ------------------------ 2.  local logger -------------------------
        # if caller passes None we fall back to a class‑scoped one
        self.log = logger or logging.getLogger(self.__class__.__name__)

        if len(inputfile_list) != 9:
            self.log.error("Error! The input file list is WRONG")
            sys.exit()

        self.input_files = inputfile_list
        self.destripe_flag = destripe_flag
        output_dir = os.path.abspath(os.path.expanduser(output_dir))
        os.makedirs(output_dir, exist_ok=True)
        self.output_dir = output_dir + os.sep
        self.mod21_file = self.input_files[0]
        self.mod06_file = self.input_files[1]
        self.mod03_file = self.input_files[2]
        self.tccloud_file = self.input_files[3]
        self.agp_file = self.input_files[4]
        self.has_mod06 = (os.path.basename(self.mod06_file) != "X") and os.path.exists(self.mod06_file)
        self.has_tc    = (os.path.basename(self.tccloud_file) != "X") and os.path.exists(self.tccloud_file)

        # ------------------------------------------------------------------
        #  downstream processors (they, too, now accept logger=None)
        # ------------------------------------------------------------------
        self.mm_processor  = MODISMISRProcessor(
            inputfile_list,
            logger=self.log,
            destripe_flag=self.destripe_flag,
        )
        self.era5_processor = ERA5Processor(inputfile_list, logger=self.log)

        self.timestamp = re.search(r"\.A\d{7}\.(\d{4})\.", inputfile_list[0]).group(1)
        self.doy       = re.search(r"\.A\d{4}(\d{3})\.", inputfile_list[0]).group(1)
        self.met_date  = np.array(
            [*map(int, re.search(r"era5_single_(\d{4})_(\d{2})_(\d{2})", inputfile_list[5]).groups()), 0],
            dtype="int32",
        ).reshape(4)
        self.mmdd = "".join(re.search(r"era5_single_(\d{4})_(\d{2})_(\d{2})", inputfile_list[5]).groups()[1:])  
        self.cth_diff = 1000 # m threshold for MISR CTH - MODID CTH
        self.misr_cth_invalid = -500
        self.modis_cth_invalid = 0
        self.cflag_nan = 0 
        self.cflag_mod_only = 1
        self.cflag_misr_only = 2
        self.cflag_misr_mod = 3
        self.cflag_misr_mod_noco2 = 4
        self.cflag_mm_valid = 5
        self.cflag_mm_invalid = 6 
        self.orbit = 'O' + str(orbit).zfill(6)

    def extract_modis_metadata(self):
        # Extract the global attributes from the MODIS 06 files... (Furture work....)
        pass
        modis_granule = self.input_files[1]
        modis_metadata = {
            'modis_granule_id': modis_granule.id,
            'modis_start_time': modis_granule.start_time,
            'modis_end_time': modis_granule.end_time,
            'satellite': modis_granule.satellite,
            'instrument': modis_granule.instrument
        }
        return modis_metadata
    def block_average_dict(self, variable_dict, block_size=10, fill_nan_when_too_small=True):
        """
        Robust block-averaging over last two axes (y, x).
        - Never raises on small arrays; returns (0,0) or full-NaN grids instead.
        - Avoids 'Mean of empty slice' by using count-aware means.
        - Casts to float32 so NaNs are always supported.
        """
        averaged_dict = {}

        def _block_reduce_last2(a, B):
            a = np.asarray(a)
            if a.ndim not in (2, 3):
                raise ValueError("array must be 2D (y,x) or 3D (z,y,x)")

            # Ensure float for NaN-safe ops (ERA5 fields are continuous anyway)
            if not np.issubdtype(a.dtype, np.floating):
                a = a.astype(np.float32, copy=False)

            if a.ndim == 2:
                H, W = a.shape
                r, c = (H // B) * B, (W // B) * B
                if r == 0 or c == 0:
                    # No full blocks available
                    return (np.empty((0, 0), dtype=np.float32)
                            if not fill_nan_when_too_small else
                            np.full((H // B, W // B), np.nan, np.float32))
                v = a[:r, :c].reshape(r // B, B, c // B, B)  # (nby,B,nbx,B)
                with np.errstate(invalid="ignore"):
                    cnt = np.sum(~np.isnan(v), axis=(1, 3))       # (nby, nbx)
                    s   = np.nansum(v,           axis=(1, 3))
                    out = np.divide(s, cnt, out=np.full_like(s, np.nan, dtype=np.float32), where=cnt > 0)
                return out.astype(np.float32, copy=False)

            else:  # 3D (z,y,x)
                Z, H, W = a.shape
                r, c = (H // B) * B, (W // B) * B
                if r == 0 or c == 0:
                    # No full blocks available
                    return (np.empty((a.shape[0], 0, 0), dtype=np.float32)
                            if not fill_nan_when_too_small else
                            np.full((a.shape[0], H // B, W // B), np.nan, np.float32))
                v = a[:, :r, :c].reshape(Z, r // B, B, c // B, B)  # (Z,nby,B,nbx,B)
                with np.errstate(invalid="ignore"):
                    cnt = np.sum(~np.isnan(v), axis=(2, 4))         # (Z, nby, nbx)
                    s   = np.nansum(v,           axis=(2, 4))
                    out = np.divide(s, cnt, out=np.full_like(s, np.nan, dtype=np.float32), where=cnt > 0)
                return out.astype(np.float32, copy=False)

        for key, variable in variable_dict.items():
            if key.endswith('exp'):
                continue
            try:
                averaged_dict[key] = _block_reduce_last2(variable, block_size)
            except Exception as e:
                # Never let one field kill the rank
                self.log.warning("block_average_dict: key=%s shape=%s fallback to None (%s)",
                                    key, getattr(variable, "shape", None), e)
                averaged_dict[key] = None
        return averaged_dict

    def convert_output_pressures_to_hpa(self, variable_dict):
        """
        ERA5 single-level pressure fields are read and used internally in Pa,
        while the saved product schema declares these fields as hPa.
        """
        for key in ("surface_pressure", "msp"):
            value = variable_dict.get(key)
            if value is None:
                continue

            arr = np.asarray(value, dtype=np.float32)
            finite = np.isfinite(arr)
            if np.any(finite) and np.nanmedian(arr[finite]) > 2000.0:
                variable_dict[key] = (arr / 100.0).astype(np.float32, copy=False)
            else:
                variable_dict[key] = arr
        return variable_dict


    def calculate_block_lat_lon(self, era5_lats, era5_lons, block_size=10, center_choice="upper"):
        """
        Center-pixel lat/lon for each full block, matching block-average shapes.
        center_choice: "upper" -> index = B//2 (e.g., 5 for B=10)
                    "lower" -> index = B//2 - 1 (e.g., 4 for B=10)
        """
        assert era5_lats.shape == era5_lons.shape, "lat/lon shapes must match"
        nrows, ncols = era5_lats.shape
        rblocks = nrows // block_size
        cblocks = ncols // block_size
        if rblocks == 0 or cblocks == 0:
            return (np.empty((0, 0), dtype=era5_lats.dtype),
                    np.empty((0, 0), dtype=era5_lons.dtype))

        if center_choice == "lower":
            offset = block_size // 2 - 1  # 4 when B=10
            if offset < 0:  # for B=1 this would be -1; clamp
                offset = 0
        else:  # "upper" (your original behavior)
            offset = block_size // 2      # 5 when B=10

        row_idx = offset + np.arange(rblocks) * block_size
        col_idx = offset + np.arange(cblocks) * block_size

        block_lat = era5_lats[np.ix_(row_idx, col_idx)]
        block_lon = era5_lons[np.ix_(row_idx, col_idx)]
        return block_lat, block_lon
    
    def save_pixels(self, mm_variables, mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name):
        # 1) Aggregate ERA5 to 10 km
        from os.path import basename
        era5_variables_10km = self.block_average_dict(era5_variables_misrswath)
        era5_variables_10km = self.convert_output_pressures_to_hpa(era5_variables_10km)
        era5_lat_10km, era5_lon_10km = self.calculate_block_lat_lon(mod_geo['lat'], mod_geo['lon'], block_size=10)
        era5_variables_10km['ERA5Latitude']  = era5_lat_10km.astype(np.float32)
        era5_variables_10km['ERA5Longitude'] = era5_lon_10km.astype(np.float32)
        # assert some2d is None or era5_lat_10km.shape == some2d.shape

        # 2) Ensure MISR input has the expected dict shape
        misr_vars = {'misrcth': misr_cth} if isinstance(misr_cth, np.ndarray) else misr_cth

        # 3) Ensure input_files is a flat list[str]
        files = self.input_files if isinstance(self.input_files, list) else [self.input_files]
        base_files = [basename(str(p)) for p in files]

        # 4) Write
        if self.mm_processor.orbit:
            orbit_num = self.mm_processor.orbit.replace('O','')
        else:
            orbit_num = self.orbit

        if self.has_mod06 and self.has_tc:
            missing_flag = 'None'
        elif not self.has_mod06:
            missing_flag = 'MODIS'
        elif not self.has_tc:
            missing_flag = 'MISR'
        else:
            missing_flag = 'No Way'
        mmcth = XarraySaver(outputfile_name, logger=self.log)
        ga = {
            "orbit_number": orbit_num,
            "missing_flag": missing_flag,
        }
        if "orig_nrows" in mod_geo:
            ga["original_modis_nrows"] = int(mod_geo["orig_nrows"])
        if "orig_ncols" in mod_geo:
            ga["original_modis_ncols"] = int(mod_geo["orig_ncols"])
        mmcth.save_mm(
            mm_variables=mm_variables,
            mm_geo=mod_geo,
            misr_variables=misr_vars,
            modis_variables=mod06,
            era5_variables=era5_variables_10km,
            input_files=base_files,
            global_attrs=ga,
            var_attrs={
                "ERA5Latitude":  {"standard_name": "latitude",  "units": "degrees_north"},
                "ERA5Longitude": {"standard_name": "longitude", "units": "degrees_east"},
            }
        )
        return
   
    def convert_misr_cth_ctp(self,era5_variables_misrswath,misr_cth):
        z_prof = np.asarray(era5_variables_misrswath['geopotential_exp'],
                            dtype=np.float32)      # (level, x, y)
        level, nx, ny = z_prof.shape
        npix_full = nx * ny
        misr_cth_2d = np.asarray(misr_cth['misrcth'], dtype=np.float32)  # (x, y), meters
        misr_cth_km = misr_cth_2d / 1e3                                  # km

        # ------------------------------------------------------------------
        # 2) Full-scene MISR CTP via Fortran routine (log-pressure interpolation)
        # ------------------------------------------------------------------
        misr_ctp_full = height_to_log_pressure(
            z_prof,          # (level, x, y)
            misr_cth_km,     # (x, y) in km
            np.int32(nx),
            np.int32(ny)
        ).astype(np.float32)  # (x, y) in hPa
        bad = (
            (~np.isfinite(misr_cth_2d)) |
            (misr_cth_2d <= self.misr_cth_invalid) |
            (~np.isfinite(misr_ctp_full)) |
            (misr_ctp_full <= 0.0)
        )
        misr_ctp_full[bad] = np.nan
        return misr_ctp_full
    
    def process_pixel_level(self,
                        bands_BT,
                        mod06,
                        mod_geo,
                        misr_cth,
                        era5_variables_misrswath,
                        met_date,
                        diagnostic_path: Optional[str] = None,
                        max_diagnostic_pixels: int = 2000):
        """
        Scene-wide call to modis.process_selected_pixels (Fortran/f2py).
        This version is careful about dtype + Fortran order to avoid f2py copies.
        """

        # ----------------------------
        # 1) Dimensions + MISR CTH/CTP
        # ----------------------------
        z_prof = np.asfortranarray(
            np.asarray(era5_variables_misrswath["geopotential_exp"], dtype=np.float32)
        )  # (level, nx, ny), F-order
        level, nx, ny = z_prof.shape
        npix_full = nx * ny

        misr_cth_2d = np.asfortranarray(
            np.asarray(misr_cth["misrcth"], dtype=np.float32)
        )  # (nx, ny), meters
        misr_cth_km = misr_cth_2d / 1e3

        misr_ctp_full = height_to_log_pressure(
            z_prof,                  # (level, nx, ny)
            misr_cth_km,             # (nx, ny)
            np.int32(nx),
            np.int32(ny),
        ).astype(np.float32, copy=False)

        # ----------------------------
        # 2) Flatten scene arrays (F)
        # ----------------------------
        w3 = np.asfortranarray(np.asarray(era5_variables_misrswath["mixingratio_exp"], dtype=np.float32))
        t3 = np.asfortranarray(np.asarray(era5_variables_misrswath["temperature_exp"], dtype=np.float32))
        h3 = z_prof  # already F float32

        wprof = w3.reshape(level, npix_full, order="F")
        tprof = t3.reshape(level, npix_full, order="F")
        hprof = h3.reshape(level, npix_full, order="F")

        psfc = (np.asfortranarray(np.asarray(era5_variables_misrswath["surface_pressure"], dtype=np.float32)) / 1e2) \
            .reshape(npix_full, order="F")
        pmsl = (np.asfortranarray(np.asarray(era5_variables_misrswath["msp"], dtype=np.float32)) / 1e2) \
            .reshape(npix_full, order="F")

        sst   = np.asfortranarray(np.asarray(era5_variables_misrswath["sst"],   dtype=np.float32)).reshape(npix_full, order="F")
        skint = np.asfortranarray(np.asarray(era5_variables_misrswath["skint"], dtype=np.float32)).reshape(npix_full, order="F")

        landsea_2d = np.asfortranarray(np.asarray(mod_geo["landsea"], dtype=np.int32))
        landsea = landsea_2d.reshape(npix_full, order="F")

        view = np.asfortranarray(np.asarray(mod_geo["vza"], dtype=np.float32)).reshape(npix_full, order="F")
        rlat = np.asfortranarray(np.asarray(mod_geo["lat"], dtype=np.float32)).reshape(npix_full, order="F")
        rlon = np.asfortranarray(np.asarray(mod_geo["lon"], dtype=np.float32)).reshape(npix_full, order="F")

        surftmp = sst.copy()  # 1D float32
        surftmp[landsea == 1] = skint[landsea == 1]

        misr_ctp_flat = np.asfortranarray(misr_ctp_full).reshape(npix_full, order="F")
        misr_ctp_raw_flat = misr_ctp_flat.copy()
        bad = (misr_ctp_flat > 1050.0)
        misr_ctp_flat[bad] = psfc[bad]

        misr_cth_flat = np.asfortranarray(misr_cth_2d).reshape(npix_full, order="F")

        # ----------------------------
        # 3) MOD06 fields (flatten F)
        # ----------------------------
        mod_cth_flat = np.asfortranarray(np.asarray(mod06["cth"], dtype=np.float32)).reshape(npix_full, order="F")
        mod_ctp_2d   = np.asfortranarray(np.asarray(mod06["ctp"], dtype=np.float32))
        mod_ctp_flat = mod_ctp_2d.reshape(npix_full, order="F")

        mod_opt_flat = np.asfortranarray(np.asarray(mod06["opt"], dtype=np.float32)).reshape(npix_full, order="F")
        mod_emi_flat = np.asfortranarray(np.asarray(mod06["emissivity"], dtype=np.float32)).reshape(npix_full, order="F")

        # mod_method: make sure it is int32, no NaNs
        ctm = np.asarray(mod06["ctm"])
        if np.issubdtype(ctm.dtype, np.floating):
            ctm = np.nan_to_num(ctm, nan=69.0, posinf=69.0, neginf=69.0)
        ctm = np.asfortranarray(ctm)  # keep Fortran layout if 2D
        mod_method_flat = np.asarray(ctm, dtype=np.int32, order="F").reshape(npix_full, order="F")

        # ----------------------------
        # 4) Build trad_scene directly: (5, npix_full) F
        # ----------------------------
        trad_scene = np.empty((5, npix_full), dtype=np.float32, order="F")
        # Here we fill row-by-row to avoid building a big (nbands,nx,ny) cube.
        band_items = list(bands_BT.values())
        for b in range(5):
            band2d = np.asfortranarray(np.asarray(band_items[5 + b], dtype=np.float32))
            trad_scene[b, :] = band2d.reshape(npix_full, order="F")

        # ----------------------------
        # 5) Clear-sky index list and clr_obs_scene
        # ----------------------------
        mod_ctp_nan = mod_ctp_2d.copy()
        mod_ctp_nan[mod_ctp_nan < 0] = np.nan

        clr_mask = np.isnan(mod_ctp_nan)
        clr_x, clr_y = np.where(clr_mask)      # 0-based indices (x,y)
        clr_lin0 = clr_x + nx * clr_y          # 0-based Fortran-linear index
        cs_idx = (clr_lin0 + 1).astype(np.int32, copy=False)  # 1-based for Fortran

        clr_obs_scene = np.full((5, npix_full), np.nan, dtype=np.float32, order="F")
        # reuse the already-built trad_scene (same 5 channels)
        clr_obs_scene[:, clr_lin0] = trad_scene[:, clr_lin0]

        # ----------------------------
        # 6) Call Fortran
        # ----------------------------
        (
            ctp_out,
            cth_out,
            emis_out,
            od_out,
            qf_out,
            ml_z_out,
            ml_flag_out,
            failure_out,
            bias_out,
        ) = modis.process_selected_pixels(
            wprof=wprof,
            tprof=tprof,
            hprof=hprof,
            psfc=psfc,
            pmsl=pmsl,
            surftmp=surftmp,
            view=view,
            trad=trad_scene,
            rlat=rlat,
            rlon=rlon,
            landsea=landsea,
            misr_ctp=misr_ctp_flat,
            misr_cth=misr_cth_flat,
            mod_cth=mod_cth_flat,
            mod_ctp=mod_ctp_flat,
            mod_method=mod_method_flat,
            mod_opt=mod_opt_flat,
            mod_emi=mod_emi_flat,
            met_date=np.asarray(met_date, dtype=np.int32),
            npix=np.int32(npix_full),
            clr_obs=clr_obs_scene,
            cs_idx=cs_idx,
        )

        # ----------------------------
        # 7) Reshape back to maps (F)
        # ----------------------------
        ctp_map  = np.asarray(ctp_out).reshape(nx, ny, order="F")
        cth_map  = np.asarray(cth_out).reshape(nx, ny, order="F")
        emis_map = np.asarray(emis_out).reshape(nx, ny, order="F")
        od_map   = np.asarray(od_out).reshape(nx, ny, order="F")
        qf_map   = np.asarray(qf_out).reshape(nx, ny, order="F")
        ml_z_map = np.asarray(ml_z_out, dtype=np.float32).reshape(nx, ny, order="F")
        ml_flag_map = np.asarray(ml_flag_out, dtype=np.int8).reshape(nx, ny, order="F")
        failure_map = np.asarray(failure_out, dtype=np.int8).reshape(nx, ny, order="F")

        if diagnostic_path is not None:
            self._save_co2_failure_diagnostics(
                diagnostic_path=diagnostic_path,
                max_pixels=max_diagnostic_pixels,
                nx=nx,
                ny=ny,
                wprof=wprof,
                tprof=tprof,
                hprof=hprof,
                psfc=psfc,
                pmsl=pmsl,
                surftmp=surftmp,
                view=view,
                trad_scene=trad_scene,
                clear_sky_bias=np.asarray(bias_out, dtype=np.float32),
                rlat=rlat,
                rlon=rlon,
                landsea=landsea,
                misr_ctp=misr_ctp_flat,
                misr_ctp_raw=misr_ctp_raw_flat,
                misr_cth=misr_cth_flat,
                mod_cth=mod_cth_flat,
                mod_ctp=mod_ctp_flat,
                mod_method=mod_method_flat,
                mod_opt=mod_opt_flat,
                mod_emi=mod_emi_flat,
                met_date=np.asarray(met_date, dtype=np.int32),
                qf_out=np.asarray(qf_out, dtype=np.int16),
                failure_out=np.asarray(failure_out, dtype=np.int16),
                ml_z_out=np.asarray(ml_z_out, dtype=np.float32),
                ml_flag_out=np.asarray(ml_flag_out, dtype=np.int16),
                mod_geo=mod_geo,
            )

        return ctp_map, cth_map, emis_map, od_map, qf_map, ml_z_map, ml_flag_map, failure_map
    
    def _save_co2_failure_diagnostics(
        self,
        *,
        diagnostic_path: str,
        max_pixels: int,
        nx: int,
        ny: int,
        wprof: np.ndarray,
        tprof: np.ndarray,
        hprof: np.ndarray,
        psfc: np.ndarray,
        pmsl: np.ndarray,
        surftmp: np.ndarray,
        view: np.ndarray,
        trad_scene: np.ndarray,
        clear_sky_bias: np.ndarray,
        rlat: np.ndarray,
        rlon: np.ndarray,
        landsea: np.ndarray,
        misr_ctp: np.ndarray,
        misr_ctp_raw: np.ndarray,
        misr_cth: np.ndarray,
        mod_cth: np.ndarray,
        mod_ctp: np.ndarray,
        mod_method: np.ndarray,
        mod_opt: np.ndarray,
        mod_emi: np.ndarray,
        met_date: np.ndarray,
        qf_out: np.ndarray,
        failure_out: np.ndarray,
        ml_z_out: np.ndarray,
        ml_flag_out: np.ndarray,
        mod_geo: dict,
    ) -> None:
        """Write profile-rich diagnostics for nonzero MM diagnostic reasons."""

        os.makedirs(os.path.dirname(diagnostic_path), exist_ok=True)

        valid_misr = np.isfinite(misr_cth) & (misr_cth > -500.0) & (misr_cth <= 20000.0)
        valid_mod = np.isfinite(mod_cth) & (mod_cth > 0.0) & (mod_cth <= 20000.0)
        eligible = valid_misr & valid_mod & (mod_cth > misr_cth)
        attempted_target = (failure_out > 0) & (failure_out != 10) & eligible
        skipped_target = (failure_out == 10) & eligible
        max_pixels = max(0, int(max_pixels))
        selected = np.flatnonzero(attempted_target)[:max_pixels]
        if selected.size < max_pixels:
            skipped = np.flatnonzero(skipped_target)[: max_pixels - selected.size]
            selected = np.concatenate([selected, skipped])

        pressure = np.asarray(self.era5_processor.P_levels, dtype=np.float32)
        bands = np.array([36, 35, 34, 33, 31], dtype=np.int16)
        pair = np.array(["36_35", "35_34", "35_33", "34_33"], dtype=object)
        n = selected.size
        nlev = pressure.size
        nb = bands.size
        npair = pair.size

        def f1(fill=np.nan, dtype=np.float32):
            a = np.empty(n, dtype=dtype)
            a[...] = fill
            return a

        def f2(shape, fill=np.nan, dtype=np.float32):
            a = np.empty(shape, dtype=dtype)
            a[...] = fill
            return a

        row = (selected % nx).astype(np.int32)
        col = (selected // nx).astype(np.int32)
        orig_row = row.copy()
        orig_col = col.copy()
        if "orig_row_index" in mod_geo:
            orig_row = np.asarray(mod_geo["orig_row_index"], dtype=np.int32)[row]
        if "orig_col_index" in mod_geo:
            orig_col = np.asarray(mod_geo["orig_col_index"], dtype=np.int32)[col]

        one_ctp = f1()
        one_cth = f1()
        one_emi = f1()
        one_od = f1()
        one_mask = f1(fill=0, dtype=np.int16)
        one_reason = f1(fill=0, dtype=np.int16)

        two_ctp = f1()
        two_cth = f1()
        two_emi = f1()
        two_od = f1()
        two_mask = f1(fill=0, dtype=np.int16)
        two_reason = f1(fill=0, dtype=np.int16)
        misr_ctp_used = f1()

        diag_isp = f1(fill=0, dtype=np.int16)
        diag_imisr = f1(fill=0, dtype=np.int16)
        diag_ltrp = f1(fill=0, dtype=np.int16)
        diag_lco2 = f1(fill=0, dtype=np.int16)
        diag_ipco2 = f1(fill=0, dtype=np.int16)
        pair_status = f2((n, npair), fill=0, dtype=np.int16)
        krto = f2((n, npair), fill=0, dtype=np.int16)
        lev = f2((n, npair))
        ctp_pres = f2((n, npair))
        tmisr = f1()
        emis_num = f1()
        emis_den = f1()
        emis_ratio = f1()
        rwcld = f1()
        tcold = f2((n, nb))
        robs = f2((n, nb))
        rclr = f2((n, nb))
        delr = f2((n, nb))
        tpad = f2((n, nlev))
        wpad = f2((n, nlev))
        ozpad = f2((n, nlev))
        taup = f2((n, nlev, nb))
        ra = f2((n, nlev, nb))
        radiance_input = f2((n, nb))
        radiance_adjusted = f2((n, nb))
        wprof_out = f2((n, nlev))
        tprof_out = f2((n, nlev))
        hprof_out = f2((n, nlev))

        for out_i, pix in enumerate(selected):
            pix = int(pix)
            trad_raw = np.asarray(trad_scene[:, pix], dtype=np.float32)
            trad_adj = np.asarray(trad_raw + clear_sky_bias, dtype=np.float32)
            radiance_input[out_i, :] = trad_raw
            radiance_adjusted[out_i, :] = trad_adj
            wprof_out[out_i, :] = wprof[:, pix]
            tprof_out[out_i, :] = tprof[:, pix]
            hprof_out[out_i, :] = hprof[:, pix]

            ctp_l = misr_ctp[pix]
            if (not np.isfinite(ctp_l)) or ctp_l <= 0.0 or ctp_l < pressure[0] or ctp_l > psfc[pix]:
                ctp_l = psfc[pix]
            misr_ctp_used[out_i] = ctp_l

            one = modis.co2cld_onepixel_misr_debug(
                wprof[:, pix], tprof[:, pix], hprof[:, pix],
                psfc[pix], pmsl[pix], surftmp[pix], view[pix],
                trad_adj, met_date, rlat[pix], rlon[pix], landsea[pix], psfc[pix],
            )
            two = modis.co2cld_onepixel_misr_debug(
                wprof[:, pix], tprof[:, pix], hprof[:, pix],
                psfc[pix], pmsl[pix], surftmp[pix], view[pix],
                trad_adj, met_date, rlat[pix], rlon[pix], landsea[pix], ctp_l,
            )

            one_ctp[out_i], one_cth[out_i], one_emi[out_i], one_od[out_i] = one[:4]
            one_mask[out_i] = one[4]
            one_reason[out_i] = one[5]

            two_ctp[out_i], two_cth[out_i], two_emi[out_i], two_od[out_i] = two[:4]
            two_mask[out_i] = two[4]
            two_reason[out_i] = two[5]
            diag_isp[out_i], diag_imisr[out_i], diag_ltrp[out_i] = two[6:9]
            diag_lco2[out_i], diag_ipco2[out_i] = two[9:11]
            pair_status[out_i, :] = two[11]
            krto[out_i, :] = two[12]
            lev[out_i, :] = two[13]
            ctp_pres[out_i, :] = two[14]
            tmisr[out_i], emis_num[out_i], emis_den[out_i] = two[15:18]
            emis_ratio[out_i], rwcld[out_i] = two[18:20]
            tcold[out_i, :] = two[20]
            robs[out_i, :] = two[21]
            rclr[out_i, :] = two[22]
            delr[out_i, :] = two[23]
            tpad[out_i, :] = two[24]
            wpad[out_i, :] = two[25]
            ozpad[out_i, :] = two[26]
            taup[out_i, :, :] = two[27]
            ra[out_i, :, :] = two[28]

        coords = {
            "diagnostic_pixel": np.arange(n, dtype=np.int32),
            "level": pressure,
            "band": bands,
            "pair": pair,
        }
        ds = xr.Dataset(coords=coords)
        dims_pix = ("diagnostic_pixel",)
        dims_prof = ("diagnostic_pixel", "level")
        dims_band = ("diagnostic_pixel", "band")
        dims_pair = ("diagnostic_pixel", "pair")
        dims_prof_band = ("diagnostic_pixel", "level", "band")

        for name, data in {
            "flat_index_fortran0": selected.astype(np.int32),
            "swath_row": row,
            "swath_col": col,
            "original_modis_row": orig_row,
            "original_modis_col": orig_col,
            "Latitude": rlat[selected].astype(np.float32),
            "Longitude": rlon[selected].astype(np.float32),
            "MM_Flag": qf_out[selected].astype(np.int16),
            "MM_FailureReason": failure_out[selected].astype(np.int16),
            "MM_MultilayerZScore": ml_z_out[selected].astype(np.float32),
            "MM_MultilayerConfidenceFlag": ml_flag_out[selected].astype(np.int16),
            "MISR_CloudTopHeight": misr_cth[selected].astype(np.float32),
            "MISR_CloudTopPressure_Input": misr_ctp_raw[selected].astype(np.float32),
            "MISR_CloudTopPressure_Used": misr_ctp_used.astype(np.float32),
            "MODIS_CloudTopHeight": mod_cth[selected].astype(np.float32),
            "MODIS_CloudTopPressure": mod_ctp[selected].astype(np.float32),
            "MODIS_CloudTopMethod": mod_method[selected].astype(np.int16),
            "MODIS_CloudOpticalDepth": mod_opt[selected].astype(np.float32),
            "MODIS_CloudEffectiveEmissivity": mod_emi[selected].astype(np.float32),
            "SurfacePressure": psfc[selected].astype(np.float32),
            "MeanSeaLevelPressure": pmsl[selected].astype(np.float32),
            "SurfaceTemperature": surftmp[selected].astype(np.float32),
            "SensorZenithAngle": view[selected].astype(np.float32),
            "LandSeaFlag": landsea[selected].astype(np.int16),
            "OneLayer_CTP": one_ctp,
            "OneLayer_CTH": one_cth,
            "OneLayer_Emissivity": one_emi,
            "OneLayer_OpticalDepth": one_od,
            "OneLayer_InternalProcessingMask": one_mask,
            "OneLayer_FailureReason": one_reason,
            "TwoLayer_CandidateCTP": two_ctp,
            "TwoLayer_CandidateCTH": two_cth,
            "TwoLayer_CandidateEmissivity": two_emi,
            "TwoLayer_CandidateOpticalDepth": two_od,
            "TwoLayer_InternalProcessingMask": two_mask,
            "TwoLayer_FailureReason": two_reason,
            "diag_isp": diag_isp,
            "diag_imisr": diag_imisr,
            "diag_ltrp": diag_ltrp,
            "diag_lco2": diag_lco2,
            "diag_ipco2": diag_ipco2,
            "diag_tmisr": tmisr,
            "emissivity_numerator": emis_num,
            "emissivity_denominator": emis_den,
            "emissivity_ratio": emis_ratio,
            "window_cloud_radiance": rwcld,
        }.items():
            ds[name] = xr.DataArray(data, dims=dims_pix)

        ds["Radiance_Input"] = xr.DataArray(radiance_input, dims=dims_band)
        ds["Radiance_Adjusted"] = xr.DataArray(radiance_adjusted, dims=dims_band)
        ds["BrightnessTemperature_FromRadiance"] = xr.DataArray(tcold, dims=dims_band)
        ds["ObservedCloudyRadiance"] = xr.DataArray(robs, dims=dims_band)
        ds["ModeledLowerBoundaryRadiance"] = xr.DataArray(rclr, dims=dims_band)
        ds["ObservedMinusModeledRadiance"] = xr.DataArray(delr, dims=dims_band)
        ds["WaterVaporProfile"] = xr.DataArray(wprof_out, dims=dims_prof)
        ds["TemperatureProfile"] = xr.DataArray(tprof_out, dims=dims_prof)
        ds["HeightProfile"] = xr.DataArray(hprof_out, dims=dims_prof)
        ds["TemperatureProfile_Padded"] = xr.DataArray(tpad, dims=dims_prof)
        ds["WaterVaporProfile_Padded"] = xr.DataArray(wpad, dims=dims_prof)
        ds["OzoneProfile_Padded"] = xr.DataArray(ozpad, dims=dims_prof)
        ds["TransmittanceProfile"] = xr.DataArray(taup, dims=dims_prof_band)
        ds["RadianceIntegral_RA"] = xr.DataArray(ra, dims=dims_prof_band)
        ds["PairStatus"] = xr.DataArray(pair_status, dims=dims_pair)
        ds["PairCrossingFlag"] = xr.DataArray(krto, dims=dims_pair)
        ds["PairCrossingLevel"] = xr.DataArray(lev, dims=dims_pair)
        ds["PairCloudTopPressure"] = xr.DataArray(ctp_pres, dims=dims_pair)
        ds["SceneClearSkyBias"] = xr.DataArray(clear_sky_bias.astype(np.float32), dims=("band",))

        units = {
            "level": "hPa",
            "Latitude": "degrees_north",
            "Longitude": "degrees_east",
            "MISR_CloudTopHeight": "m",
            "MISR_CloudTopPressure_Input": "hPa",
            "MISR_CloudTopPressure_Used": "hPa",
            "MODIS_CloudTopHeight": "m",
            "MODIS_CloudTopPressure": "hPa",
            "MODIS_CloudOpticalDepth": "1",
            "MODIS_CloudEffectiveEmissivity": "1",
            "SurfacePressure": "hPa",
            "MeanSeaLevelPressure": "hPa",
            "SurfaceTemperature": "K",
            "SensorZenithAngle": "degrees",
            "OneLayer_CTP": "hPa",
            "OneLayer_CTH": "km",
            "OneLayer_Emissivity": "1",
            "OneLayer_OpticalDepth": "1",
            "TwoLayer_CandidateCTP": "hPa",
            "TwoLayer_CandidateCTH": "km",
            "TwoLayer_CandidateEmissivity": "1",
            "TwoLayer_CandidateOpticalDepth": "1",
            "diag_tmisr": "K",
            "WaterVaporProfile": "g kg-1",
            "TemperatureProfile": "K",
            "HeightProfile": "km",
            "TemperatureProfile_Padded": "K",
            "WaterVaporProfile_Padded": "g kg-1",
            "OzoneProfile_Padded": "ppmv",
            "TransmittanceProfile": "1",
            "BrightnessTemperature_FromRadiance": "K",
            "emissivity_ratio": "1",
            "PairCloudTopPressure": "hPa",
        }
        for name, unit in units.items():
            if name in ds:
                ds[name].attrs["units"] = unit

        for name in [
            "Radiance_Input", "Radiance_Adjusted", "ObservedCloudyRadiance",
            "ModeledLowerBoundaryRadiance", "ObservedMinusModeledRadiance",
            "RadianceIntegral_RA", "SceneClearSkyBias", "window_cloud_radiance",
        ]:
            if name in ds:
                ds[name].attrs["units"] = "radiance units used internally by modis_co2_slice.f90"

        failure_flag_values = np.array(
            [0, 10, 20, 21, 22, 23, 24, 30, 40, 41, 42, 50, 51, 52, 60, 61, 62, 70, 71],
            dtype=np.int16,
        )
        failure_flag_meanings = (
            "ok "
            "not_attempted "
            "invalid_misr_ctp "
            "invalid_surface_pressure "
            "invalid_modis_radiance "
            "invalid_atmospheric_profile "
            "pressure_grid_mapping_failed "
            "no_vertical_search_domain "
            "signal_below_noise "
            "no_valid_ratio "
            "no_ratio_crossing "
            "candidate_not_above_misr "
            "candidate_outside_band_range "
            "no_selectable_band_pair "
            "candidate_not_higher_than_one_layer "
            "reference_one_layer_failed "
            "candidate_has_no_valid_ctp_digit "
            "ctp_ok_emissivity_denominator_bad "
            "ctp_ok_emissivity_ratio_bad"
        )
        for name in ["MM_FailureReason", "OneLayer_FailureReason", "TwoLayer_FailureReason"]:
            ds[name].attrs.update({
                "flag_values": failure_flag_values,
                "flag_meanings": failure_flag_meanings,
            })
        ds["PairStatus"].attrs.update({
            "flag_values": np.array([0, 1, 2, 3, 4, 5], dtype=np.int16),
            "flag_meanings": (
                "not_evaluated skipped_radiance_signal no_valid_ratio_levels "
                "no_zero_crossing zero_crossing_found selected_solution"
            ),
        })
        ds.attrs.update({
            "title": "CO2-slicing MM failure diagnostics",
            "description": (
                "Profile-rich diagnostics for pixels with a nonzero MM "
                "two-layer diagnostic reason."
            ),
            "source": "src/pixel/modis_co2_slice.f90 co2cld_onepixel_misr_debug",
            "diagnostic_pixel_selection": (
                "MM_FailureReason > 0, valid MISR/MODIS CTH, and MODIS CTH > MISR CTH; "
                "non-not_attempted reasons are selected first"
            ),
            "max_diagnostic_pixels": int(max_pixels),
            "selected_diagnostic_pixels": int(n),
            "algorithm_reference": "Mitra et al. 2023, JGR Atmospheres, doi:10.1029/2022JD038135",
        })
        ds.to_netcdf(diagnostic_path)
        self.log.info("Saved CO2 failure diagnostics to %s with %d pixels", diagnostic_path, n)

    
    def _nan2d_like(self, geo2d, dtype=np.float32):
        a = np.empty_like(geo2d, dtype=dtype)
        a[:] = np.nan
        return a

    def _zeros2d_like(self, geo2d, dtype=np.int8):
        return np.zeros_like(geo2d, dtype=dtype)
    
    def run_process(self,save_flag: str = "not_debug"):
        """
        Run the MISR/MODIS/ERA5 data processing pipeline and save outputs.

        Parameters:
            save_flag (str): Indicates the type of output saving strategy. Default is 'not_debug'.
        """
        if self.has_mod06 and self.has_tc:
            file_flag = 'MM'
        elif not self.has_mod06:
            file_flag = 'MI'
        elif not self.has_tc:
            file_flag = 'MO'
        else:
            file_flag = 'NW'

        if self.mm_processor.orbit is not None:
            outputfile_name =  f'{self.output_dir}MODMISR_L2_CP_061_{self.mm_processor.id}_{self.mmdd}_{self.mm_processor.orbit}_{file_flag}_V00.nc'
        elif self.orbit is not None:
            outputfile_name =  f'{self.output_dir}MODMISR_L2_CP_061_{self.mm_processor.id}_{self.mmdd}_{self.orbit}_{file_flag}_V00.nc'
        else:
            return
          
        misr_cth_invalid = -600
        self.log.debug("Running the MISR/MODIS data processing pipeline")
        co2_diagnostic = save_flag in {"co2_debug", "diagnostic", "mm_debug"}
        
        #!  Finish the part results only return two datasets 
        res = self.mm_processor.mm_process(
            process_misr_cth=self.has_tc,
            process_mod06_cth=self.has_mod06,
            scale_flag = save_flag, 
        )
        if res is None:
            self.log.error("Processing MISR/MODIS failed")
            return

        bands_BT = res['bands_BT']      # may be None
        misr_cth = res['misr_cth']      # may be None (or dict)
        mod_geo  = res['mod_geo']       # always present (dict with lat/lon/landsea/vza)
        mod06    = res['mod06']    
            
        # # Starting Processing ERA5 data
        # print(mod_geo['lat'].shape,mod_geo['lon'].shape)
        era5_lats_1d, era5_lons_1d = self.era5_processor.era5_lat_lon()
        era5_variables = self.era5_processor.era5_process()
        
        era5_variables_misrswath = self.era5_processor.map_era5_to_modis(
            mod_geo['lat'], mod_geo['lon'],
            era5_lats_1d.values, era5_lons_1d.values,
            era5_variables
        )

        if era5_variables_misrswath is None:
            self.log.error(
                "ERA5->MODIS mapping failed; skipping. "
                "mod21km =%s mod06=%s mod06=%s misr=%s",
                self.mod21_file,
                self.mod06_file,
                self.mod03_file,
                self.tccloud_file, 
            )
            return None
        del era5_variables #release MEM

        self.log.debug("ERA5 processing completed successfully")

        if isinstance(misr_cth, dict) and misr_cth.get('misrcth') is not None:
            misr_arr = np.asarray(misr_cth['misrcth'])
            if np.any(misr_arr > misr_cth_invalid):
                misr_cth['misrctp'] = self.convert_misr_cth_ctp(era5_variables_misrswath, misr_cth)

        # save_flag is for debugging 
        if save_flag == 'debug':
            outputfile_name = f'{self.output_dir}/debug/MODMISR_L2_CP_{self.mm_processor.id}_{self.mm_processor.orbit}_debug.nc'
            print(outputfile_name)
            mmcth = NetCDFSaver(outputfile_name, self.log)
            mmcth.save_org(bands_BT,misr_cth, mod_geo, mod06,era5_variables_misrswath,outputfile = outputfile_name)
            self.log.debug("Output file saved successfully")                  
            return 
        H, W = mod_geo['lat'].shape

        # Initializing the MM Variables 
        if mod06 is None:
            mm_variables = {
                            'ctp':        self._nan2d_like(mod_geo['lat']),
                            'emissivity': self._nan2d_like(mod_geo['lat']),
                            'opt':        self._nan2d_like(mod_geo['lat']),
                            'cflag':      self._zeros2d_like(mod_geo['lat'], dtype=np.int8),
                            'cth':        self._nan2d_like(mod_geo['lat']),
                        }
            misr_arr = misr_cth['misrcth'] if (isinstance(misr_cth, dict) and 'misrcth' in misr_cth) else None
            if (misr_arr is None) or (np.all(misr_arr) <= misr_cth_invalid) :
                self.log.error("MISR CTH dict missing 'misrcth'")
                return
            else:
                valid = (misr_arr > self.misr_cth_invalid)
                mm_variables['cth'][valid]   = misr_arr[valid]
                misr_ctp = misr_cth['misrctp']
                mm_variables['ctp'][valid]   = misr_ctp[valid]
                mm_variables['cflag'][valid] = self.cflag_misr_only
                self.save_pixels(mm_variables, None, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
                return
        mod_ctp = np.asarray(mod06['ctp'])
        mm_variables = {
            'ctp':        mod_ctp,          # view
            'emissivity': np.asarray(mod06['emissivity']),
            'opt':        np.asarray(mod06['opt']),
            'cflag':      np.zeros_like(mod06['ctp'], dtype=np.int8),
            'cth':        np.asarray(mod06['cth']),                
        }        
        mm_variables['cflag'][mod_ctp>0] = 1   
        
        # === Case: MOD06 only (no MISR CTH) ===
        if (misr_cth['misrcth'] is None) or (np.all(misr_cth['misrcth']) <= misr_cth_invalid):
            self.log.info("MISR CTH missing; writing MODIS-only product.")
            misr_cth = {}
            misr_cth['misrcth'] = self._zeros2d_like(mod_geo['lat']) - np.nan# Need to be modifed 
            misr_cth['misrctp'] = self._nan2d_like(mod_geo['lat'])
            misr_cth['misrcth_qa'] = self._zeros2d_like(mod_geo['lat'], dtype=np.int8)
            self.save_pixels(mm_variables, mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
            return
        # === Case: MISR and MODIS Co-Ex ===
         
        # Proceed with pixel-level processing  
        
        self.log.debug("Proceeding with pixel-level processing...")     
        diagnostic_path = None
        if co2_diagnostic:
            diagnostic_path = (
                f"{self.output_dir}debug/"
                f"MM_CO2FailureDiagnostics_{self.mm_processor.id}_{self.mm_processor.orbit}.nc"
            )

        ctp_map, cth_map, emis_map, od_map, qf_map, ml_z_map, ml_flag_map, failure_map = self.process_pixel_level(
            bands_BT,
            mod06,
            mod_geo,
            misr_cth,
            era5_variables_misrswath,
            self.met_date,
            diagnostic_path=diagnostic_path,
        )
        mm_variables['ctp']  =  ctp_map
        mm_variables['cth']  =  cth_map
        mm_variables['emissivity'] = emis_map
        mm_variables['opt'] = od_map
        mm_variables['cflag']= qf_map 
        mm_variables['ml_zscore'] = ml_z_map
        mm_variables['ml_flag'] = ml_flag_map
        mm_variables['failure_reason'] = failure_map
        self.save_pixels(mm_variables, mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
        return

   
