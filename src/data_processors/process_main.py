import os,sys
import importlib
import numpy as np
from pyresample import kd_tree, geometry
from pyresample.geometry import SwathDefinition
import re
from .netcdf_saver import NetCDFSaver
from scipy.interpolate import interp1d,PchipInterpolator
import scipy
import pandas as pd
import xarray as xr

class MODISMISRProcessor:
    '''    '''
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
            
            bands_BT = {f'band_{band}': mod21.get_BT(str(band)) for band in [31, 32, 33, 34, 35, 36]}
            bands_radiance = {f'band_{band}': mod21.get_band_data(str(band), 'radiance') for band in [29, 31, 32, 4]}
            
            self.logger.debug("MOD21 data processed successfully")
            
            combined_bands = {**bands_BT, **bands_radiance}
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
                'emissivity': mod06.get_emissivity()
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
            'lat': self.apply_valid_indices(mod_vza, valid_rows, valid_cols),
            'lon': self.apply_valid_indices(mod_vza, valid_rows, valid_cols),
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
        
        mod_geo_misrswath = {
            'lat': self.apply_valid_indices(mod_vza, valid_rows, valid_cols),
            'lon': self.apply_valid_indices(mod_vza, valid_rows, valid_cols),
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
            next_day = target_time >= 23.0
            
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
        lat_indices = np.searchsorted(era5_lats, rounded_modis_lats) - 1
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
            result = {
                'dew2m': era5_single.get_d2m(),
                'sst': era5_single.get_sst(),
                'temp2m': era5_single.get_t2m(),
                'surface_pressure': era5_single.get_sp(),
                'skint': era5_single.get_skt(),        
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
                'specific_humidity': era5_multi.get_specific_humidity(),
                'geopotential': era5_multi.get_geopotential()
            }
            self.logger.debug("ERA5 data (pressure level) processed successfully")
            return result
        except Exception as e:
            self.logger.error(f"Error processing ERA5 data (pressure level): {e}")
            return None

    def remove_below_surface_levels(self,sp_single,P_multi, T_multi, Z_multi, Q_multi):
        """
        Removes data levels below the surface pressure.

        Parameters:
            P_multi (ndarray): Pressure data.
            sp_single (ndarray): Surface pressure data.
            *args (ndarray): Variable number of arrays containing data at different pressure levels.
        Returns:
            tuple: Modified P_multi, followed by modified arrays in args and kwargs.
        """
        P_multi = np.array(P_multi, dtype=np.float64)
        Z_multi = np.array(Z_multi, dtype=np.float64)
        sp_single = np.array(sp_single, dtype=np.float64)

        valid_mask =((P_multi <= sp_single[None, :, :]) & (Z_multi > 0))

        # Process P_multi
        P_multi[~valid_mask] = np.nan
        # Process args
        T_multi[~valid_mask] = np.nan
        Z_multi[~valid_mask] = np.nan
        Q_multi[~valid_mask] = np.nan
        return P_multi, T_multi, Z_multi, Q_multi

    def avoid_duplicate_z_levels(self, Z_multi, T_multi, Q_multi):
        """
        Removes duplicate Z levels due to rounding.
        
        Parameters:
            Z_multi (ndarray): Geopotential height data.
            T_multi (ndarray): Temperature data.
            Q_multi (ndarray): Specific humidity data.
            
        Returns:
            tuple: Modified T_multi, Z_multi, and Q_multi arrays.
        """
        roundzero_mask = np.round(Z_multi[-1, :, :], 3) == 0
        Z_multi[-1, roundzero_mask] = np.nan
        Q_multi[-1, roundzero_mask] = np.nan
        T_multi[-1, roundzero_mask] = np.nan
        return T_multi, Z_multi, Q_multi

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

    def interpolate_to_pressure_levels_direct(self, variable):
        pchip_interpolator = PchipInterpolator((np.log(self.era5_P_levels)), variable,axis=0,extrapolate=True)
        return pchip_interpolator(np.log(self.P_levels_gt_1mb))
         
        
    def load_std_atmos_profile(self):
        netcdf_file = self.profile_file
        month = self.month
        ds = xr.open_dataset(netcdf_file)
        if month in (4, 5, 6, 7, 8, 9):
            dataset = ds['summer_spring']
        else:
            dataset = ds['fall_winter']
        
        return dataset

    def era5_process(self):
        """
        Interpolates multi-level ERA5 data to 101 pressure .
        
        Parameters:
            profiles (optional): Profile data for determining profiles based on latitude.
        
        Return:

            turple: EAR5 surface/near-surface data along with profile data that have been interpolated into 101 pressure levels 
        """
        latitudes, _ = self.era5_lat_lon()
        latitudes = latitudes.values
        
        var_multi = self.ear5_read_multi(self.era5_multi_file)
        var_single = self.ear5_read_single(self.era5_single_file)

        var_multi_next_day = None
        var_single_next_day = None
        target_hour = float(self.timestamp[0:2])
        if target_hour >= 23.0:
            var_multi_next_day = self.ear5_read_multi(self.era5_multi_file_next_day) 
            var_single_next_day = self.ear5_read_multi(self.era5_single_file_next_day) 
        
        var_single = self.era5_timeinterpolate(var_single, var_single_next_day)
        var_multi = self.era5_timeinterpolate(var_multi, var_multi_next_day)
    
        T_multi = var_multi['temperature']
        Q_multi = var_multi['specific_humidity']
        Z_multi = var_multi['geopotential'] / 9.80665e3

       
        # Avoid duplicated Z levels due to rounding.
        # T_multi, Z_multi, Q_multi = self.avoid_duplicate_z_levels(Z_multi, T_multi, Q_multi)

         
        # Interpolate the data into pressure_levels larger than 1mb
        
        T_interpolated = self.interpolate_to_pressure_levels_direct(T_multi)
        Q_interpolated = self.interpolate_to_pressure_levels_direct(Q_multi)
        Z_interpolated = self.interpolate_to_pressure_levels_direct(Z_multi)
        
        # Combine the standard atmospheric profile data for pressure levels smaller than 1mb
        profile_data = self.load_std_atmos_profile()

        temperature_profile = profile_data.sel(variable="T(K)").values
        print(temperature_profile.shape,T_interpolated.shape)
        offset_T = T_interpolated[0, :, :] - temperature_profile[-1, :, :]
  
        temperature_profile += offset_T[np.newaxis, :, :]
        T_expand = np.concatenate((temperature_profile,T_interpolated),axis=0)
        
        sq_profile = profile_data.sel(variable="h2o(cm-3)").values
        offset_Q = Q_interpolated[0, :, :] - sq_profile[-1, :, :]
        sq_profile += offset_Q[np.newaxis, :, :] 
        Q_expand = np.concatenate((sq_profile,Q_interpolated),axis=0)

        Z_profile = profile_data.sel(variable="z(km)").values
        offset_Z =  Z_interpolated[0, :, :] - Z_profile[-1, :, :]
        Z_profile += offset_Z[np.newaxis, :, :]
        Z_expand = np.concatenate((Z_profile,Z_interpolated),axis=0)
         
        # Set data below SST

        result = var_single.copy()
        result.update({
            'temperature': T_expand,
            'specific_humidity': Q_expand,
            'geopotential': Z_expand,
        })

        # internal_variable = {
        #     't_org': T_multi,
        #     'z_org': Z_multi,
        #     'q_org': Q_multi,
        # }
        # expand_variable = {
        #     'temperature': T_expand,
        #     'specific_humidity': Q_expand,
        #     'geopotential': Z_expand,
        # }

        # self.save_to_netcdf('../output/era5_org_variables.nc',**internal_variable)
        # self.save_to_netcdf('../output/era5_int_variables.nc',**expand_variable)
 
        return result

class MainProcessor:
    def __init__(self, inputfile_list, logger):
        print(inputfile_list)
        print(f'the input filelist lengith is {len(inputfile_list)}')
        if len(inputfile_list) != 9:
            logger.error("Error! The input file list is WRONG")
            sys.exit()
        self.logger = logger
        self.timestamp = re.search(r'\.A\d{7}\.(\d{4})\.', inputfile_list[0]).group(1)
        self.mm_processor = MODISMISRProcessor(inputfile_list, logger)
        self.era5_processor = ERA5Processor(inputfile_list, logger)

    def run_process(self):
        # try:
            self.logger.debug("Running the MISR/MODIS data processing pipeline")
            bands_BT, misr_cth, mod_geo, mod06 = self.mm_processor.mm_process()
            if  bands_BT is None:
                self.logger.error(f" Processing MISR/MODIS Problems: {e}")
            era5_lats, era5_lons = self.era5_processor.era5_lat_lon()
            era5_variables = self.era5_processor.era5_process()
            self.logger.debug("ERA5 processing completed successfully")
            era5_variables_misrswath = self.era5_processor.map_era5_to_modis(mod_geo['lat'], mod_geo['lon'], era5_lats, era5_lons, era5_variables)
            
             
            output_dir = '/data/keeling/a/gzhao1/f/mmcth/output/'
            outputfile_name = f'{output_dir}MODMISR_L2_CP_{self.mm_processor.id}_O{self.mm_processor.orbit}_v01.nc'
            print(outputfile_name)
            mmcth = NetCDFSaver(outputfile_name, self.logger)
            mmcth.save(bands_BT,misr_cth, mod_geo, mod06,era5_variables_misrswath)
            self.logger.debug("Output file saved successfully")           

        # except Exception as e:
        #     self.logger.error(f" Main Processing Problems: {e}")
        #     return None


