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
    
    def process_mod21(self, save_format='normal'):
        try:
            self.log.debug("Processing MOD21 data")
            mod21 = self.mod_read.MODL1Granule(self.mod21_file)
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
        mod_lat, mod_lon, mod_landsea, mod_vza = self.process_mod03()
        if any(v is None for v in (mod_lat, mod_lon, mod_landsea)):
            self.log.error("Error retrieving MODIS lat/lon/land‑sea mask")
            return None
        if np.isnan(mod_lat).any() or np.isnan(mod_lon).any() or (mod_lat<-444.0).any() or (mod_lon <-444.0).any() :
            self.log.error("NaN detected in MODIS lat/lon; returning None.")
            return None         
        
        # ---------- 2. Read MISR geo ----------
        misr_lat, misr_lon = self.get_misr_geo()
        if any(v is None for v in (misr_lat, misr_lon)):
            self.log.error("Error retrieving MISR lat/lon")
            return None

        # ---------- 3. Overlap mask ----------
        misr_mod_swath = self.misr_to_modis(misr_lat, misr_lon, mod_lat, mod_lon, misr_lat)
        if np.all(misr_mod_swath<=-999) or misr_mod_swath is None :
            print('return')
            self.log.error("Error finding MISR swath on MODIS grid")
            return None
        valid_rows, valid_cols = self.cutoff_misr_swath(misr_mod_swath)
        if any(v is None for v in (valid_rows, valid_cols)):
            self.log.error("Error retrieving MISR lat/lon")
            return None 
        
        # Helper to slice MODIS arrays to the MISR footprintprint(valid_rows)
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
            qi_data_int8 = self._to_qi_int8(qi_data)  
            out['misr_cth'] = {
                'misrcth'    : self.misr_to_modis(misr_lat, misr_lon, out['mod_geo']['lat'], out['mod_geo']['lon'], misr_cth),
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
        T_skin = np.asfortranarray(var_single['skint'].astype(np.float32))
        T_sst = np.asfortranarray(var_single['sst'].astype(np.float32))

    
        num_modis_levels = np.int32(modis_levels_1d.shape[0])
        lat_size = np.int32(T_multi.shape[1])
        lon_size = np.int32(T_multi.shape[2])
         
        T_modis = np.empty((num_modis_levels, lat_size, lon_size), dtype=np.float32, order='F')
        W_modis = np.empty((num_modis_levels, lat_size, lon_size), dtype=np.float32, order='F')
        Z_modis = np.empty((num_modis_levels, lat_size, lon_size), dtype=np.float32, order='F')

        interpolate_to_pressure_levels(
            era5_levels_1d, modis_levels_1d, T_multi, W_multi, Z_multi,
            surface_pressures, W_surface, T_skin, T_sst, T_modis, W_modis, Z_modis,
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
        output_dir: str = "/data/keeling/a/gzhao1/f/mmcth/ds_output/",
        save_flag: str = 'not_debug',
    ):
        # ------------------------ 2.  local logger -------------------------
        # if caller passes None we fall back to a class‑scoped one
        self.log = logger or logging.getLogger(self.__class__.__name__)

        if len(inputfile_list) != 9:
            self.log.error("Error! The input file list is WRONG")
            sys.exit()

        self.input_files = inputfile_list
        self.output_dir = output_dir
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
        self.mm_processor  = MODISMISRProcessor(inputfile_list, logger=self.log)
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
        self.output_dir = output_dir

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
                self.logger.warning("block_average_dict: key=%s shape=%s fallback to None (%s)",
                                    key, getattr(variable, "shape", None), e)
                averaged_dict[key] = None
        return averaged_dict


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
        mmcth.save_mm(
            mm_variables=mm_variables,
            mm_geo=mod_geo,
            misr_variables=misr_vars,
            modis_variables=mod06,
            era5_variables=era5_variables_10km,
            input_files=base_files,
            global_attrs={
                # CF/global attrs usually avoid spaces in keys; underscores are safer:
                "orbit_number": orbit_num,
                "Missing Flag": missing_flag,
            },
            var_attrs={
                "ERA5Latitude":  {"standard_name": "latitude",  "units": "degrees_north"},
                "ERA5Longitude": {"standard_name": "longitude", "units": "degrees_east"},
            }
            # var_attrs=...  # optional per-variable attr overrides
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
        return misr_ctp_full
    
    def process_pixel_level(self,
                        bands_BT,
                        mod06,
                        mod_geo,
                        misr_cth,
                        era5_variables_misrswath,
                        met_date):
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
        ctp_out, cth_out, emis_out, od_out, qf_out = modis.process_selected_pixels(
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

        return ctp_map, cth_map, emis_map, od_map, qf_map
    
    
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
                misr_ctp = self.convert_misr_cth_ctp(era5_variables_misrswath,misr_cth)
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
            misr_cth['misrcth_qa'] = self._zeros2d_like(mod_geo['lat'], dtype=np.int8)
            self.save_pixels(mm_variables, mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
            return
        # === Case: MISR and MODIS Co-Ex ===
         
        # Proceed with pixel-level processing  
        
        self.log.debug("Proceeding with pixel-level processing...")     
        ctp_map,cth_map,emis_map,od_map,qf_map = self.process_pixel_level(bands_BT, mod06,mod_geo, misr_cth, era5_variables_misrswath,self.met_date)
        mm_variables['ctp']  =  ctp_map
        mm_variables['cth']  =  cth_map
        mm_variables['emissivity'] = emis_map
        mm_variables['opt'] = od_map
        mm_variables['cflag']= qf_map 
        self.save_pixels(mm_variables, mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
        return

   


