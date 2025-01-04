import os,sys
import importlib
import re

import numpy as np
import scipy
from scipy.interpolate import interp1d,PchipInterpolator
import pandas as pd
import xarray as xr
from pyresample import kd_tree, geometry
from pyresample.geometry import SwathDefinition

from .netcdf_saver import NetCDFSaver
from src.pixel import modis
from .misr_cth_to_pressure import height_to_log_pressure
from .interpolate_to_pressure_levels import interpolate_to_pressure_levels

class MODISMISRProcessor:
    '''  Process MODIS data  and MISR TC_Cloud data if avaible    '''
    def __init__(self, input_files, logger):
        if len(input_files) != 9:
            logger.error("Error! The input file list is WRONG")
            sys.exit()    
        self.mod21_file = input_files[0]
        self.mod06_file = input_files[1]
        self.mod03_file = input_files[2]
        self.tccloud_file = input_files[3]
        self.agp_file = input_files[4]
        self.logger = logger
        self.mod_read = importlib.import_module('src.data_readers.mod_read')
        self.misr_read = importlib.import_module('src.data_readers.misr_read')
        self.id = re.search(r'A\d{7}\.\d{4}', self.mod21_file).group() if re.search(r'A\d{7}\.\d{4}', self.mod21_file) else None
        self.orbit = re.search(r'O\d{6}', self.tccloud_file).group() if re.search(r'O\d{6}', self.tccloud_file) else None
    
    def process_mod21(self):
        try:
            self.logger.debug("Processing MOD21 data")
            mod21 = self.mod_read.MODL1Granule(self.mod21_file)
            bands_BT = {f'band_{band}': mod21.get_BT(str(band)) for band in [36, 35, 34, 33, 31]}
            # bands_radiance = {f'band_{band}': mod21.get_band_data(str(band), 'radiance') for band in [29, 31, 32, 4]}
            self.logger.debug("MOD21 data processed successfully")
            # combined_bands = {**bands_BT, **bands_radiance}
            combined_bands = {**bands_BT}
            return combined_bands
        except Exception as e:
            self.logger.error(f"Error processing MOD21: {self.mod21_file} {e}")
            return None

    def process_mod06(self):
        try:
            self.logger.debug("Processing MOD06 data")
            mod06 = self.mod_read.MOD06Granule(self.mod06_file)
            result = {
                'ctp': mod06.get_ctp(),
                'opt': mod06.get_opt(),
                'cphase': mod06.get_cphase(),
                'cth': mod06.get_cth(),
                'emissivity': mod06.get_emissivity(),
                'ctm': mod06.get_ctm(),
            }
            self.logger.debug("MOD06 data processed successfully")
            return result
        except Exception as e:
            self.logger.error(f"Error processing MOD06: {self.mod06_file} {e}")
            return None
        
    def process_mod03(self):
        try:
            self.logger.debug("Processing MOD03 data")
            mod03 = self.mod_read.MOD03Granule(self.mod03_file)
            result = mod03.get_lat(), mod03.get_lon(), mod03.get_landsea_mask(), mod03.get_vza()
            self.logger.debug("MOD03 data processed successfully")
            return result
        except Exception as e:
            self.logger.error(f"Error processing MOD03: {self.mod03_file} {e}")
            return None
        
    def get_misr_geo(self):
        try:
            agp = self.misr_read.MISRGranule(self.agp_file)
            misr_lat, misr_lon = agp.get_geo()
            return misr_lat, misr_lon
        except Exception as e:
            self.logger.error(f"Error getting MISR geo: {e}")
            return None   
        self.logger.error(f"Error getting MISR geo: {e}")
        return None
    
    def get_misr_cth(self):
            try:
                tc_cloud = self.misr_read.MISRGranule(self.tccloud_file)
                cth,cth_qa = tc_cloud.get_cth()
                return  cth,cth_qa
            except Exception as e:
                self.logger.error(f"Error getting MISR CTH: {e}")
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
            self.logger.error(f"Error in MISR to MODIS reprojection: {e}")
            return None

    def cutoff_misr_swath(self, misr_swath):
        try:
            misr_swath[misr_swath < -990] = np.nan
            mask = ~np.isnan(misr_swath)
            valid_rows = np.any(mask, axis=1)
            valid_cols = np.any(mask, axis=0)
            return valid_rows, valid_cols
        except Exception as e:
            self.logger.error(f"Error in cutoff MISR swath: {e}")
            return None
        
    def apply_valid_indices(self, array, valid_rows, valid_cols):
        try:
            return array[valid_rows][:, valid_cols]
        except Exception as e:
            self.logger.error(f"Error applying valid indices: {e}")
            return None

    def mm_process_modis_only(self):
        mod_lat, mod_lon, mod_landsea, mod_vza = self.process_mod03()
        if mod_lat is None or mod_lon is None or mod_landsea is None:
            self.logger.error("Error retrieving MODIS Lat, Lon, or Landsea mask")
            return None

        misr_lat, misr_lon = self.get_misr_geo()
        if misr_lat is None or misr_lon is None:
            self.logger.error("Error retrieving MISR Lat, Lon")
            return None

        misr_mod_swath = self.misr_to_modis(misr_lat, misr_lon, mod_lat, mod_lon, misr_lat)
        if misr_mod_swath is None:
            self.logger.error("Error finding MISR Swath over MODIS grid")
            return None
        
        valid_rows, valid_cols = self.cutoff_misr_swath(misr_mod_swath)
        if valid_rows is None or valid_cols is None:
            self.logger.error("Error cutting off invalid parts of MISR swath")
            return None
 
        mod_lat_misrswath = self.apply_valid_indices(mod_lat, valid_rows, valid_cols)
        mod_lon_misrswath = self.apply_valid_indices(mod_lon, valid_rows, valid_cols)
        if mod_lat_misrswath is None or mod_lon_misrswath is None:
            self.logger.error("Error applying valid indices to MODIS swath")
            return None

        mod06_product = self.process_mod06()
        if mod06_product:
            mod06_misrswath = [self.apply_valid_indices(mod06_product[key], valid_rows, valid_cols) for key in mod06_product]
            mod06_misrswath_dict = {key: value for key, value in zip(mod06_product.keys(), mod06_misrswath)}

        else:
            self.logger.error("Error retrieving MOD06 data")
            return None

        bands_BT = self.process_mod21()
        bands_BT_misrswath = [self.apply_valid_indices(bands_BT[key], valid_rows, valid_cols) for key in bands_BT]
        bands_BT_misrswath_dict = {key: value for key, value in zip(bands_BT.keys(), bands_BT_misrswath)}

        mod_geo_misrswath = {
            'lat': self.apply_valid_indices(mod_lat, valid_rows, valid_cols),
            'lon': self.apply_valid_indices(mod_lon, valid_rows, valid_cols),
            'landsea': self.apply_valid_indices(mod_landsea, valid_rows, valid_cols),
            'vza': self.apply_valid_indices(mod_vza, valid_rows, valid_cols)
        }

        return bands_BT_misrswath_dict, mod_geo_misrswath, mod06_misrswath_dict
                
    def mm_process(self, process_misr_cth=True, process_mod06_cth=True ):

        mod_lat, mod_lon, mod_landsea, mod_vza = self.process_mod03()
        if mod_lat is None or mod_lon is None or mod_landsea is None:
            self.logger.error("Error retrieving MODIS Lat, Lon, or Landsea mask")
            return None

        misr_lat, misr_lon = self.get_misr_geo()
    
        if misr_lat is None or misr_lon is None:
            self.logger.error("Error retrieving MISR Lat, Lon")
            return None
        # The reprojection code requires the padding values np.nan instead of any acually numbers like -999
        misr_mod_swath = self.misr_to_modis(misr_lat, misr_lon, mod_lat, mod_lon, misr_lat)
        # print(misr_mod_swath[misr_mod_swath>-990])
        if misr_mod_swath is None:
            self.logger.error("Error finding MISR Swath over MODIS grid")
            return None
        
        valid_rows, valid_cols = self.cutoff_misr_swath(misr_mod_swath)
        if valid_rows is None or valid_cols is None:
            self.logger.error("Error cutting off invalid parts of MISR swath")
            return None
 
        mod_lat_misrswath = self.apply_valid_indices(mod_lat, valid_rows, valid_cols)
        mod_lon_misrswath = self.apply_valid_indices(mod_lon, valid_rows, valid_cols)
        if mod_lat_misrswath is None or mod_lon_misrswath is None:
            self.logger.error("Error applying valid indices to MODIS swath")
            return None
        
        mod_geo_misrswath = {
            'lat': self.apply_valid_indices(mod_lat, valid_rows, valid_cols),
            'lon': self.apply_valid_indices(mod_lon, valid_rows, valid_cols),
            'landsea': self.apply_valid_indices(mod_landsea, valid_rows, valid_cols),
            'vza': self.apply_valid_indices(mod_vza, valid_rows, valid_cols)
        }
        
        misr_cth_dict = None
        mod06_misrswath_dict = None
        bands_BT_misrswath_dict = None

        if process_mod06_cth and process_misr_cth:
            bands_BT = self.process_mod21()
            bands_BT_misrswath = [self.apply_valid_indices(bands_BT[key], valid_rows, valid_cols) for key in bands_BT]
            bands_BT_misrswath_dict = {key: value for key, value in zip(bands_BT.keys(), bands_BT_misrswath)}
            mod06_product = self.process_mod06()
            if mod06_product:
                mod06_misrswath = [self.apply_valid_indices(mod06_product[key], valid_rows, valid_cols) for key in mod06_product]
                mod06_misrswath_dict = {key: value for key, value in zip(mod06_product.keys(), mod06_misrswath)}
            else:
                self.logger.error("Error retrieving MOD06 data")
                return None        
            misr_cth,misr_cth_qa = self.get_misr_cth()
            if misr_cth is None:
                self.logger.error("Error retrieving MISR cloud top height data")
                return None
            misr_cth_modisswath = self.misr_to_modis(misr_lat, misr_lon, mod_lat_misrswath, mod_lon_misrswath, misr_cth)
            misr_cth_modisswath_qa = self.misr_to_modis(misr_lat, misr_lon, mod_lat_misrswath, mod_lon_misrswath, misr_cth_qa)
            misr_cth_dict = {
                 'misrcth': misr_cth_modisswath,
                 'misrcth_qa': misr_cth_modisswath_qa,
            }
            return bands_BT_misrswath_dict, misr_cth_dict, mod_geo_misrswath, mod06_misrswath_dict
        elif not process_mod06_cth and process_misr_cth: 
            misr_cth,misr_cth_qa = self.get_misr_cth()
            if misr_cth is None:
                self.logger.error("Error retrieving MISR cloud top height data")
                return None
            misr_cth_modisswath = self.misr_to_modis(misr_lat, misr_lon, mod_lat_misrswath, mod_lon_misrswath, misr_cth)
            misr_cth_modisswath_qa = self.misr_to_modis(misr_lat, misr_lon, mod_lat_misrswath, mod_lon_misrswath, misr_cth_qa)
            misr_cth_dict = {
                 'misr_cth': misr_cth_modisswath,
                 'misr_cth_qa': misr_cth_modisswath_qa,
            }
            return misr_cth_modisswath, mod_geo_misrswath 
        elif not process_misr_cth and process_mod06_cth : 
            mod06_product = self.process_mod06()
            if mod06_product:
                mod06_misrswath = [self.apply_valid_indices(mod06_product[key], valid_rows, valid_cols) for key in mod06_product]
                mod06_misrswath_dict = {key: value for key, value in zip(mod06_product.keys(), mod06_misrswath)}
            else:
                self.logger.error("Error retrieving MOD06 data")
                return None                 
            return  mod06_misrswath_dict, mod_geo_misrswath
                

class ERA5Processor:
    def __init__(self, input_files, logger, latlon_idx = None):
        self.logger = logger
        self.era5_single_file = input_files[5]
        self.era5_multi_file = input_files[6]
        self.era5_single_file_next_day = input_files[7]
        self.era5_multi_file_next_day = input_files[8]
        self.mod21_file  = input_files[0]
        self.era5_read = importlib.import_module('src.data_readers.era5_read')
        self.month =re.search(r"_\d{4}_(\d{2})_\d{2}\.nc", self.era5_multi_file).group(1) 
        self.MODIS = importlib.import_module('src.pixel.modis')
        self.P_levels = self.MODIS.transmission.pstd
        self.P_levels_gt_1mb = self.P_levels[self.P_levels>=1]
        self.era5_P_levels = np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000], dtype=int)
        self.latlon_idx = latlon_idx
        self.timestamp = re.search(r"\.(\d{4})\.",self.mod21_file).group(1) 
        self.profile_file = '/data/gdi/f/gzhao1/mmcth/externaldata/std_atmos_profile.nc'

    def era5_timeinterpolate(self, era5_vars_dict, era5_vars_dict_next_day=None):
        try:
            self.logger.debug("Interpolating ERA5 data to target time")
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
            self.logger.debug("ERA5 data interpolated to target time successfully")
            return interpolated_data
        except Exception as e:
            self.logger.error(f"Error interpolating ERA5 data: {e}")
            return None
    
    def era5_lat_lon(self):
        try:
            self.logger.debug("Retrieving ERA5 lat/lon grid")
            era5_single = self.era5_read.ERA5Single(self.era5_single_file)
            result = era5_single.get_latitude(), era5_single.get_longitude()
            self.logger.debug("ERA5 lat/lon grid retrieved successfully")
            return result
        except Exception as e:
            self.logger.error(f"Error retrieving ERA5 lat/lon grid: {e}")
            return None

    def ear5_read_single(self,era5_single_file):
        ''' Read ERA5 surface/near-surface variables'''
        try:
            self.logger.debug("Processing ERA5 data (surface layer)")
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
            self.logger.debug("ERA5 data (surface layer) processed successfully")
            return result
        except Exception as e:
            self.logger.error(f"Error processing ERA5 data (surface layer): {e}")
            return None

    def ear5_read_multi(self,era5_multi_file):
        ''' Read ERA5 pressure-slevel variables'''
        try:
            self.logger.debug("Processing ERA5 data (pressure level)")
            era5_multi = self.era5_read.ERA5Multi(era5_multi_file)
            result = {
                'temperature': era5_multi.get_temperature(),
                'specific_humidity': era5_multi.get_specific_humidity()*1e3,  # convert the unit of kg/kg to g/kg
                'geopotential': era5_multi.get_geopotential()
            }
            self.logger.debug("ERA5 data (pressure level) processed successfully")
            return result
        except Exception as e:
            self.logger.error(f"Error processing ERA5 data (pressure level): {e}")
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
        self.logger.info(f"Saved data to {filename}")
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
            self.logger.debug("Mapping ERA5 variables to MODIS grid points")
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
            self.logger.debug("ERA5 variables mapped to MODIS grid points successfully")
            return modis_vars
        except Exception as e:
            self.logger.error(f"Error mapping ERA5 to MODIS: {e}")
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
        # print(f'W profile: {modis.transmission.wstd}')

        # print(f'T profile: {modis.transmission.tstd}')

        # print(f'psfc profile: {modis.transmission.pstd}')



        latitudes, _ = self.era5_lat_lon()
        latitudes = latitudes.values
        
        var_single = self.ear5_read_single(self.era5_single_file)
        var_multi = self.ear5_read_multi(self.era5_multi_file)
        
        var_single_next_day = None
        var_multi_next_day = None
        
        target_hour = float(self.timestamp[0:2])
        if target_hour >= 23.0:
            var_multi_next_day = self.ear5_read_multi(self.era5_multi_file_next_day) 
            var_single_next_day = self.ear5_read_single(self.era5_single_file_next_day) 
        
        var_single = self.era5_timeinterpolate(var_single, var_single_next_day)
        var_multi = self.era5_timeinterpolate(var_multi, var_multi_next_day)
    
        if no_interploate: 
            result = var_single.copy()
            result.update(var_multi)
            return result



        T_multi = var_multi['temperature'] 
        W_multi = var_multi['specific_humidity']/1000/(1-var_multi['specific_humidity']/1000) # check again .....
        Z_multi = var_multi['geopotential'] / 9.80665e3
         
        # Interpolate the data into pressure_levels larger than 1mb
        # 'dew2m': dew2m_mixingratio,
        # 'sst': era5_single.get_sst(),
        # 'temp2m': era5_single.get_t2m(),
        # 'surface_pressure': era5_single.get_sp(),
        # 'skint': era5_single.get_skt(),  
        # 'msp' : era5_single.get_msl(),      
 

        # for key in ['surface_pressure', 'dew2m', 'skint', 'sst']:
        #     if key in var_single:
        #         value = var_single[key][200, 400]
        #         print(f"{key} at [200, 400]: {value.values}")
        #     else:
        #         print(f"{key} is not found in var_single")

         
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
        offset_T = T_interpolated[0, :, :] - temperature_profile[-1, :, :]
  
        temperature_profile += offset_T[np.newaxis, :, :]
        T_expand = np.concatenate((temperature_profile,T_interpolated),axis=0)
        
        # need to verify each variable
        h20_profile = profile_data.sel(variable="h2o(cm-3)").values #unit g/kg
        air_profile = profile_data.sel(variable="air(cm-3)").values #unit g/kg
        sq_profile = h20_profile/air_profile*18.016/28.966*1e3 
        w_profile = sq_profile/1000.0/(1-sq_profile/1000.0)

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
      
        })
         
        internal_variable = {
            't_org': T_multi,
            'z_org': Z_multi,
            'q_org': W_multi*1000.0,
        }
        expand_variable = {
            'temperature': T_expand,
            'mixingratio': W_expand,
            'geopotential': Z_expand,
            
        }
        self.save_to_netcdf(f'../output/check_org.nc',**internal_variable)
        self.save_to_netcdf(f'../output/check_int.nc',**expand_variable)
        self.save_to_netcdf(f'../output/check_single.nc',**var_single)
 
        return result 


class MainProcessor:

    def __init__(self, inputfile_list, logger, output_dir = '/data/keeling/a/gzhao1/f/mmcth/output/'):
        if len(inputfile_list) != 9:
            logger.error("Error! The input file list is WRONG")
            sys.exit()
        self.inputfiles = inputfile_list
        self.logger = logger
        self.timestamp = re.search(r'\.A\d{7}\.(\d{4})\.', inputfile_list[0]).group(1)
        self.mm_processor = MODISMISRProcessor(inputfile_list, logger)
        self.era5_processor = ERA5Processor(inputfile_list, logger)
        self.doy = re.search(r'\.A\d{4}(\d{3})\.', inputfile_list[0]).group(1)
        
        self.cth_diff = 1000
        self.misr_cth_invalid = -500
        self.modis_cth_invalid = 0
        
        self.cflag_corrected = 1
        self.cflag_misr = 2
        self.cflag_modis = 3
        self.cflag_invalid = 0 

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

    def block_average_dict(self,variable_dict, block_size=10):
        averaged_dict = {}
        for key, variable in variable_dict.items():
            if key.endswith('exp'):
                continue    
            if variable.ndim == 2:  # Surface variables (x, y)
                x, y = variable.shape
                y_truncated = (y // block_size) * block_size  # Truncate y to be divisible by block_size
                
                # Reshape into blocks and average, ignoring NaNs
                variable_truncated = variable[:, :y_truncated]
                variable_blocked = variable_truncated.reshape(x // block_size, block_size, 
                                                            y_truncated // block_size, block_size)
                averaged_variable = np.nanmean(variable_blocked, axis=(1, 3))  # Average, ignoring NaNs
                
            elif variable.ndim == 3:  # Profile data (37, x, y)
                profiles, x, y = variable.shape
                y_truncated = (y // block_size) * block_size  # Truncate y to be divisible by block_size
                
                # Reshape into blocks and average, ignoring NaNs
                variable_truncated = variable[:, :, :y_truncated]
                variable_blocked = variable_truncated.reshape(profiles, x // block_size, block_size,
                                                            y_truncated // block_size, block_size)
                averaged_variable = np.nanmean(variable_blocked, axis=(2, 4))  # Average, ignoring NaNs
                
            else:
                raise ValueError(f"Variable '{key}' must be 2D or 3D")
            averaged_dict[key] = averaged_variable

        return averaged_dict


    def calculate_block_lat_lon(self,era5_lats, era5_lons, block_size=10):
        """
        Calculate the latitude and longitude for the center pixel of each ERA5 block.
        
        Parameters:
            era5_lats (np.ndarray): Original 1 km latitude array for ERA5 data.
            era5_lons (np.ndarray): Original 1 km longitude array for ERA5 data.
            block_size (int): Size of the block (e.g., 10x10). Default is 10.
            
        Returns:
            block_lat (np.ndarray): Latitude array for the center of each block.
            block_lon (np.ndarray): Longitude array for the center of each block.
        """
       
        assert era5_lats.shape == era5_lons.shape, "Latitude and longitude arrays must have the same shape."

        # Get the dimensions of the original lat/lon arrays
        lat_dim, lon_dim = era5_lats.shape

        # Calculate the indices for the center pixels of each block
        block_center_lat_idx = np.arange(block_size // 2, lat_dim, block_size)
        block_center_lon_idx = np.arange(block_size // 2, lon_dim, block_size)

        # Extract the center pixel lat/lon for each block
        block_lat = era5_lats[block_center_lat_idx, :][:, block_center_lon_idx]
        block_lon = era5_lons[block_center_lat_idx, :][:, block_center_lon_idx]

        return block_lat, block_lon

    def valid_pixel_condition(self, bands_BT, mod06, misr_cth):
        mod_cth = np.array(mod06['cth'])
        # print(f'mod_cth: {np.shape(mod_cth)}')
        # print(f'misr_cth: {np.shape(mod_cth)}')
        mod_method = mod06['ctm']
        bands_stack = np.stack(list(bands_BT.values()), axis=0)  # Shape: (num_bands, height, width)
        mod_bt_all_bands = np.all(bands_stack > 0, axis=0) 
        height_diff = mod_cth - np.array(misr_cth)

        condition = (height_diff > self.cth_diff) & (mod_method < 5) & mod_bt_all_bands
        return condition

    def save_pixels(self, mm_variables,mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name):
        # Perform block averaging for ERA5 variables and calculate block lat/lon
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        era5_variables_10km = self.block_average_dict(era5_variables_misrswath)
        # era5_lats_10km, era5_lons_10km = self.calculate_block_lat_lon(mod_geo['lat'], mod_geo['lon'])
        # era5_geo_10km = {'era5_lat': era5_lats_10km, 'era5_lon': era5_lons_10km}
        # print(era5_lats_10km.shape,era5_lons_10km.shape)
        # # Save the processed data using NetCDFSaver
        mmcth = NetCDFSaver(outputfile_name, self.logger)
        print('**********************************')
        mmcth.save_mm(mm_variables, mod_geo, misr_cth, mod06, era5_variables_10km)
        return      

    def process_pixel_level(self,selected_pixels,bands_BT, mod_geo, misr_cth, era5_variables_misrswath):
        z_prof = era5_variables_misrswath['geopotential'].values
        selected_z_prof = z_prof[selected_pixels[0],selected_pixels[1]]     
        
        # Initialize the output array
        misr_ctp = np.empty_like(misr_cth)
        misr_ctp = height_to_log_pressure(z_prof, misr_cth/1000., z_prof.shape[1], z_prof.shape[2])
        selected_sst = era5_variables_misrswath['sst'].values[selected_pixels[0],selected_pixels[1]] 
        selected_dew2m  = era5_variables_misrswath['dew2m'].values[selected_pixels[0],selected_pixels[1]] 
        selected_landsea = mod_geo['landsea'].values[selected_pixels[0],selected_pixels[1]]  # landsea = 11
        selected_skintmp = era5_variables_misrswath['skint'].values[selected_pixels[0],selected_pixels[1]]
        selected_sp_prof = era5_variables_misrswath['mixingratio'].values[::-1, selected_pixels[0],selected_pixels[1]]/100 #Wprof - 1
        selected_t_prof = era5_variables_misrswath['temperature'].values[::-1, selected_pixels[0],selected_pixels[1]] #tprof - 2
        selected_surfp  = era5_variables_misrswath['surface_pressure'].values[selected_pixels[0],selected_pixels[1]] /100 # psfc - 3
        selected_seasurfp  = era5_variables_misrswath['msp'].values[selected_pixels[0],selected_pixels[1]] /100 # pmsl - 4
        selected_surftemp = np.copy(selected_sst)
        selected_surftemp[selected_landsea == 1] =  selected_skintmp[selected_landsea == 1] #surftmp - 5
        selected_vza  = mod_geo['vza'].values[selected_pixels[0],selected_pixels[1]] #view -6
        selected_bt = bands_BT[:,selected_pixels[0],selected_pixels[1]]   #tcold -7
        # selected_misrcth = misr_cth[selected_pixels[0],selected_pixels[1]]   #MISR_CTH 
        selected_misrctp = misr_ctp[selected_pixels[0],selected_pixels[1]]  #MISR_CTP - 13
        # selected_modcth = mod_cth[selected_pixels[0],selected_pixels[1]]   #MISR_CTH 
        # selected_modctp = mod_ctp[selected_pixels[0],selected_pixels[1]]   #MISR_CTP
        selected_lat = mod_geo['latitude'].values[selected_pixels[0],selected_pixels[1]] #latitude -6
        selected_lon = mod_geo['longitude'].values[selected_pixels[0],selected_pixels[1]] #latitude -6
    
        npix = selected_sst.shape[0]  # npix should be 27944
        ctp,_ = modis.process_selected_pixels(selected_sp_prof/100, selected_t_prof, selected_surfp, selected_seasurfp,
                                  selected_surftemp, selected_vza, selected_bt,selected_lat,
                                   selected_lon, selected_landsea,selected_misrctp,self.doy,npix)
        
        return ctp
    
    def run_process(self,save_flag = 'org'):
        """
        Run the MISR/MODIS/ERA5 data processing pipeline and save outputs.

        Parameters:
            save_flag (str): Indicates the type of output saving strategy. Default is 'org'.
        """
        outputfile_name =  f'{self.output_dir}MODMISR_L2_CP_{self.mm_processor.id}_{self.mm_processor.orbit}_v01.nc'
        misr_cth_invalid = -500

        self.logger.debug("Running the MISR/MODIS data processing pipeline")
        results = self.mm_processor.mm_process()

        if results is None:
            self.logger.error("Processing MISR/MODIS failed")
            return

        bands_BT, misr_cth, mod_geo, mod06 = results

        # # Starting Processing ERA5 data

        era5_lats_1d, era5_lons_1d = self.era5_processor.era5_lat_lon()
        era5_variables = self.era5_processor.era5_process()
        era5_variables_misrswath = self.era5_processor.map_era5_to_modis(
            mod_geo['lat'], mod_geo['lon'], era5_lats_1d.values, era5_lons_1d.values, era5_variables
        )
          
        
        self.logger.debug("ERA5 processing completed successfully")
   
        if save_flag == 'org':
            output_dir = '/data/keeling/a/gzhao1/f/mmcth/output/'
            outputfile_name = f'{output_dir}MODMISR_L2_CP_{self.mm_processor.id}_{self.mm_processor.orbit}_debug.nc'
            print(outputfile_name)
            mmcth = NetCDFSaver(outputfile_name, self.logger)
            mmcth.save_org(bands_BT,misr_cth, mod_geo, mod06,era5_variables_misrswath,outputfile = outputfile_name)
            self.logger.debug("Output file saved successfully")                  
            return 
        
        if np.all(misr_cth['misrcth'] <= misr_cth_invalid):
            self.logger.debug("No valid MISR CTH or pixel condition; Skip pixel-level processing")
            mm_variables = np.copy(mod06)
            mm_variables['cflag'] = np.zeros_like(mod06['cth']) + self.cflag_modis
            mm_variables['cflag'][mod06['cth']>self.mod_cth_invalid] =  self.cflag_modis
            self.save_pixels(mm_variables,mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
            return
        
        #  Check if no valid pixels for processing
        valid_pixels = self.valid_pixel_condition(bands_BT,mod06, misr_cth['misrcth'])  
        if not np.any(valid_pixels):
            self.logger.debug("No valid pixel for pixel processing condition; Skip pixel-level processing")
            mm_variables = np.copy(mod06)
            mm_variables['cflag'] = np.zeros_like(mod06['cth']) + self.cflag_modis
            mm_variables['cth'][misr_cth['cth'] > self.misr_cth_invalid] = misr_cth['cth'][misr_cth['cth'] > self.misr_cth_invalid]
            mm_variables['cflag'][misr_cth['cth'] > self.misr_cth_invalid] = self.cflag_misr 
            self.save_pixels(mm_variables,mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
            return
    
        # Proceed with pixel-level processing  
        self.logger.debug("Proceeding with pixel-level processing...")
        cflag = np.zeros_like(mod06['cth']) + self.cflag_modis
        cflag[misr_cth['misrcth'] > self.misr_cth_invalid] = self.cflag_misr
        selected_pixels = np.where(valid_pixels)  
         
        # print(f'selected_pixel: {np.shape(selected_pixels)}') 
        # corrected_ctp = self.process_pixel_level(selected_pixels,bands_BT, mod_geo, misr_cth, era5_variables_misrswath,self.doy)
        corrected_ctp = np.full_like(selected_pixels[0].shape,99)

        cflag = np.zeros_like(mod06['cth'])
        # print(f'cflag: {np.shape(cflag[selected_pixels])}') 
        cflag[selected_pixels] = np.where(corrected_ctp > 0, self.cflag_corrected, cflag[selected_pixels])
         
        mm_variables = mod06.copy()
        misr_cth['misrcth'] = np.array(misr_cth['misrcth'])

        # print("Shape of mm_variables['cth']:", mm_variables['cth'].shape)
        # print("Shape of misr_cth['misrcth']:", misr_cth['misrcth'].shape)
        mm_variables['cth'][misr_cth['misrcth'] > self.misr_cth_invalid] = misr_cth['misrcth'][misr_cth['misrcth'] > self.misr_cth_invalid]
        mm_variables['cth'][selected_pixels] = np.where(corrected_ctp > 0, corrected_ctp, mm_variables['cth'][selected_pixels])
        mm_variables['cflag'] = cflag
        self.save_pixels(mm_variables,mod06, mod_geo, misr_cth, era5_variables_misrswath, outputfile_name)
        return
         


