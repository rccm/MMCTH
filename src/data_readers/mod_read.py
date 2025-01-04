
from pyhdf.SD import SD, SDC
from pyhdf.HDF import *
from pyhdf.VS import *
import numpy as np
import re
import os
from datetime import datetime
class MODL1Granule():
    def __init__(self, mod_l1b_fullpath):
        self.fullpath = mod_l1b_fullpath  # Use the passed full path
        if os.path.isfile(self.fullpath):
            self.hdf_file = SD(self.fullpath)
        else:
            print(f'Corrupted file: {self.fullpath}')    
            raise FileNotFoundError("File not found: " + self.fullpath)
    def get_id(self):
        match = re.search(r"\.A(\d{7})\.(\d{4})\.", self.fullpath)  # Use self.fullpath to extract date/time
        if match:
            return match.group(0)
    def get_datetime(self):
        match = re.search(r"\.A(\d{7})\.(\d{4})\.", self.fullpath)  # Use self.fullpath to extract date/time
        if match:
            julian_date = match.group(1)
            time_str = match.group(2)
            date = datetime.strptime(julian_date, "%Y%j")
            time = datetime.strptime(time_str, "%H%M").time()
            return julian_date, date, time
        else:
            print("No match found in the filename.")
            return None, None, None
    def get_yr(self):
        julian_date, date, time = self.get_datetime()
        return julian_date[0:4]
    def get_doy(self):
        julian_date, date, time = self.get_datetime()
        return julian_date[4:]
    def get_month(self):
        julian_date, date, time = self.get_datetime()
        return date[4:6]
    def get_time(self):
        julian_date, date, time = self.get_datetime()
        return time
    def get_coremetadata(self):
        hdf = HDF(self.fullpath, HC.READ)  # Open the HDF file
        try:
            vs = hdf.vstart()  # Initialize the VS (Vdata) interface
            vdata = vs.attach('CoreMetadata.0')  # Access the CoreMetadata.0 Vdata by name
            metadata_records = vdata.read()  # Read the metadata fields
            metadata = ''.join(sum(metadata_records, []))  # Process the metadata
        finally:
            vdata.detach()  # Detach the Vdata
            vs.end()  # End the VS interface
            hdf.close()  # Close the HDF file
        return metadata
    def get_misrtime(self):
        # Parse the input datetime string
        _,dt,_ = self.get_datetime()
        if dt:
            formatted_str = dt.strftime("%Y.%m.%d")
            return formatted_str
        else:
            return None

    def parse_metadata_for_day_night_flag(self, metadata):
        match = re.search(r'OBJECT\s*=\s*DAYNIGHTFLAG\s*NUM_VAL\s*=\s*1\s*VALUE\s*=\s*"([^"]+)"', metadata)
        return match.group(1) if match else "DayNightFlag not found"

    def parse_metadata_for_orbit_number(self, metadata):
        match = re.search(r'OBJECT\s*=\s*ORBITNUMBER\s*CLASS\s*=\s*"1"\s*NUM_VAL\s*=\s*1\s*VALUE\s*=\s*(\d+)', metadata)
        return int(match.group(1)) if match else "Orbit Number not found"

    def get_orbit(self):
        metadata = self.get_coremetadata()
        return self.parse_metadata_for_orbit_number(metadata)

    def get_daynight(self):
        metadata = self.get_coremetadata()
        return self.parse_metadata_for_day_night_flag(metadata)
    
    
    def get_earth_sun_dist(self):
        file_ = SD(self.fullpath)
        earthsundist = getattr(file_, 'Earth-Sun Distance')
        file_.end()

        return earthsundist

    def get_data(self, fieldname, SD_field_rawData):
        hdf_file = SD(self.fullpath)
        if SD_field_rawData==0:
            data = hdf_file #science data set
        elif SD_field_rawData==1:
            data = hdf_file.select(fieldname) #field
        else:
            data = hdf_file.select(fieldname).get() #raw data
        return data, hdf_file
    
    def get_scale_and_offset(self,data_field, rad_or_ref=True):
        '''
        INPUT
            data:       numpy float array - get_data(filename,fieldname,SD_or_rawData=1)
            rad_or_ref: boolean           - True if radaince or False if reflectance offsets/scale desired
        RETURN
            2 numpy float arrays, scale factor & offset of size=number of bands
            for the chosen field
        '''
        if rad_or_ref:
            offset_name = 'radiance_offsets'
            scale_name  = 'radiance_scales'
        else:
            offset_name = 'reflectance_offsets'
            scale_name  = 'reflectance_scales'

        #grab scale and offset
     
        scale_factor = data_field.attributes()[scale_name]
        offset = data_field.attributes()[offset_name]
        return scale_factor, offset

    def get_radiance_or_reflectance(self,data_raw, data_field, rad_or_ref = True, scale_factor=True):
        '''
        INPUT
            data_raw:   get_data(filename, fieldname, SD_field_rawData=2)
            data_field: get_data(filename, fieldname, SD_field_rawData=1)
            rad_or_ref: boolean - True if radiance, False if reflectance
            scale_factor: boolean - return this as well if True
        RETURN
            radiance: numpy float array - shape=(number of bands, horizontal, vertical)
        '''
        #get dimensions of raw data
        num_bands = np.ma.size(data_raw, axis=0)
        num_horizontal = np.ma.size(data_raw, axis=1)
        num_vertical = np.ma.size(data_raw, axis=2)

        #reshape raw data to perform scale and offset correction
        data_raw_temp = np.reshape(data_raw,(num_bands, num_horizontal * num_vertical))
        scale_factor, offset = self.get_scale_and_offset(data_field, rad_or_ref=rad_or_ref)

        #fill values found in data
        detector_saturated = 65533
        detector_dead      = 65531
        missing_data       = 65535
        max_DN             = 32767
        min_DN             = 0

        #to replace DN outside of range and not any of the other fill vals
        fill_val_bad_data  = -999

        #save indices of where bad values occured occured
        detector_saturated_idx = np.where(data_raw_temp == detector_saturated)
        detector_dead_idx      = np.where(data_raw_temp == detector_dead)
        missing_data_idx       = np.where(data_raw_temp == missing_data)
        over_DN_max_idx        = np.where(data_raw_temp >  max_DN)
        below_min_DN_idx       = np.where(data_raw_temp <  min_DN)

        #mark all invalid data as nan
        data_raw_temp = data_raw_temp.astype(np.float)
        data_raw_temp[over_DN_max_idx]  = np.nan
        data_raw_temp[below_min_DN_idx] = np.nan
        #correct raw data to get radiance/reflectance values
        #correct first band manually
        data_corrected_total = (data_raw_temp[0,:] - offset[0]) * scale_factor[0]
        #for loop to put all the bands together
        for i in range(1,num_bands):
            #corrected band
            data_corrected = (data_raw_temp[i,:] - offset[i]) * scale_factor[i]

            #aggregate bands
            data_corrected_total = np.vstack((data_corrected_total, data_corrected))

        #add fill values back in
        data_corrected_total[over_DN_max_idx]        = fill_val_bad_data
        data_corrected_total[below_min_DN_idx]       = fill_val_bad_data
        data_corrected_total[detector_saturated_idx] = detector_saturated
        data_corrected_total[detector_dead_idx]      = detector_dead
        data_corrected_total[missing_data_idx]       = missing_data

        #get original shape and return radiance/reflectance
        if not scale_factor:
            return data_corrected_total.reshape((num_bands, num_horizontal, num_vertical))
        else:
            return data_corrected_total.reshape((num_bands, num_horizontal, num_vertical)) 
        
    def get_band_data(self, band_name, data_type):
        """
        Retrieve radiance or reflectance data for a specified band.

        Parameters:
        - band_name (str): Name of the band to retrieve.
        - data_type (str): Type of data to retrieve ('radiance' or 'reflectance').

        Returns:
        - numpy.ndarray: Scaled data for the specified band.
        
        Raises:
        - ValueError: If the band name is not found or if an invalid data type is specified.
        """
        hdf_file = SD(self.fullpath, SDC.READ)
        
        emissive_dataset = hdf_file.select('EV_1KM_Emissive')
        refsb_dataset = hdf_file.select('EV_1KM_RefSB')
        refsb_dataset_250m = hdf_file.select('EV_250_Aggr1km_RefSB')
        refsb_dataset_500m = hdf_file.select('EV_500_Aggr1km_RefSB')

        emissive_band_names = emissive_dataset.attributes()['band_names'].split(',')
        refsb_band_names = refsb_dataset.attributes()['band_names'].split(',')
        refsb_band_250m_names = refsb_dataset_250m.attributes()['band_names'].split(',')
        refsb_band_500m_names = refsb_dataset_500m.attributes()['band_names'].split(',')

        if band_name in emissive_band_names:
            dataset = emissive_dataset
            band_names = emissive_band_names
            scales = dataset.attributes()['radiance_scales']
            offsets = dataset.attributes()['radiance_offsets']
        elif band_name in refsb_band_names:
            dataset = refsb_dataset
            band_names = refsb_band_names
            if data_type == 'radiance':
                scales = dataset.attributes()['radiance_scales']
                offsets = dataset.attributes()['radiance_offsets']
            elif data_type == 'reflectance':
                scales = dataset.attributes()['reflectance_scales']
                offsets = dataset.attributes()['reflectance_offsets']
            else:
                raise ValueError(f"Invalid data_type '{data_type}' specified. Use 'radiance' or 'reflectance'.")
        elif  band_name in  refsb_band_250m_names: 
            dataset = refsb_dataset_250m
            band_names =refsb_band_250m_names
            if data_type == 'radiance':
                scales = dataset.attributes()['radiance_scales']
                offsets = dataset.attributes()['radiance_offsets']
            elif data_type == 'reflectance':
                scales = dataset.attributes()['reflectance_scales']
                offsets = dataset.attributes()['reflectance_offsets']
            else:
                raise ValueError(f"Invalid data_type '{data_type}' specified. Use 'radiance' or 'reflectance'.")
        elif  band_name in  refsb_band_500m_names: 
            dataset = refsb_dataset_500m
            band_names =refsb_band_500m_names
            if data_type == 'radiance':
                scales = dataset.attributes()['radiance_scales']
                offsets = dataset.attributes()['radiance_offsets']
            elif data_type == 'reflectance':
                scales = dataset.attributes()['reflectance_scales']
                offsets = dataset.attributes()['reflectance_offsets']
            else:
                raise ValueError(f"Invalid data_type '{data_type}' specified. Use 'radiance' or 'reflectance'.")
 
        else:
            raise ValueError(f"Band name '{band_name}' not found in either dataset")

        # Find the index of the requested band name
        band_index = band_names.index(band_name)

        # Extract the data for the specific band
        band_data = dataset[band_index, :, :]

        # Set values outside the valid range to np.nan
        valid_range = dataset.attributes()['valid_range']
        band_data = np.where((band_data >= valid_range[0]) & (band_data <= valid_range[1]), band_data, np.nan)

        # Scale the data using the appropriate scales and offsets
        band_data = (band_data - offsets[band_index]) * scales[band_index]
        return band_data
    
    def get_BT(self, band_name, units=1, satellite='TERRA'):
        """
        Compute brightness temperature for a MODIS infrared band.
        
        Parameters:
        rad       : 2D array of Planck radiance (units determined by UNITS)
        band      : MODIS IR band number (20-25, 27-36)
        units     : Flag defining radiance units
                    0 => milliWatts per square meter per steradian per inverse centimeter
                    1 => Watts per square meter per steradian per micron
        satellite : "TERRA" or "AQUA"
        
        Returns:
        Brightness temperature (Kelvin) as a 2D array
        """
        rad = self.get_band_data(str(band_name),'radiance')
        band = int(band_name)
        if band < 20 or band > 36 or band == 26:
            raise ValueError("Argument BAND must be in the range [20-25, 27-36]")
        
        if satellite.upper() == 'TERRA':
            # Effective central wavenumber (inverse centimeters)
            cwn = [
                2641.775, 2505.277, 2518.028, 2465.428, 2235.815, 2200.346, 0.0, 1477.967,
                1362.737, 1173.190, 1027.715, 908.0884, 831.5399, 748.3394, 730.8963, 718.8681, 704.5367
            ]
            # Temperature correction slope (no units)
            tcs = [
                0.9993411, 0.9998646, 0.9998584, 0.9998682, 0.9998819, 0.9998845, 0.0, 0.9994877,
                0.9994918, 0.9995495, 0.9997398, 0.9995608, 0.9997256, 0.9999160, 0.9999167, 0.9999191, 0.9999281
            ]
            # Temperature correction intercept (Kelvin)
            tci = [
                0.4770532, 0.09262664, 0.09757996, 0.08929242, 0.07310901, 0.07060415, 0.0, 0.2204921,
                0.2046087, 0.1599191, 0.08253401, 0.1302699, 0.07181833, 0.01972608, 0.01913568, 0.01817817, 0.01583042
            ]
        elif satellite.upper() == 'AQUA':
            # Effective central wavenumbers (inverse centimeters)
            cwn = [
                2647.409, 2511.760, 2517.908, 2462.442, 2248.296, 2209.547, 0.0, 1474.262,
                1361.626, 1169.626, 1028.740, 907.6813, 830.8411, 748.2978, 730.7766, 718.2094, 703.5007
            ]
            # Temperature correction slopes (no units)
            tcs = [
                0.9993363, 0.9998626, 0.9998627, 0.9998707, 0.9998737, 0.9998770, 0.0, 0.9995694,
                0.9994867, 0.9995270, 0.9997382, 0.9995270, 0.9997271, 0.9999173, 0.9999070, 0.9999198, 0.9999233
            ]
            # Temperature correction intercepts (Kelvin)
            tci = [
                0.4818401, 0.09426663, 0.09458604, 0.08736613, 0.07873285, 0.07550804, 0.0, 0.1848769,
                0.2064384, 0.1674982, 0.08304364, 0.1343433, 0.07135051, 0.01948513, 0.02131043, 0.01804156, 0.01683156
            ]
        else:
            raise ValueError('You must select either "TERRA" or "AQUA" as your satellite')

        def bright_m(w, r):
            """
            Compute brightness temperature given monochromatic Planck radiance.
            """
            # Constants
            h = 6.6260755e-34  # Planck constant (Joule second)
            c = 2.9979246e+8  # Speed of light in vacuum (meters per second)
            k = 1.380658e-23  # Boltzmann constant (Joules per Kelvin)

            c1 = 2.0 * h * c * c
            c2 = (h * c) / k

            # Convert wavelength to meters
            ws = 1.0E-6 * w

            # Compute brightness temperature
            return c2 / (ws * np.log(c1 / (1.0e+6 * r * ws**5) + 1.0))

        def brite_m(v, r):
            """
            Compute brightness temperature given monochromatic Planck radiance.
            """
            # Constants
            h = 6.6260755e-34  # Planck constant (Joule second)
            c = 2.9979246e+8  # Speed of light in vacuum (meters per second)
            k = 1.380658e-23  # Boltzmann constant (Joules per Kelvin)

            c1 = 2.0 * h * c * c
            c2 = (h * c) / k

            # Convert wavenumber to inverse meters
            vs = 1.0e+2 * v

            # Compute brightness temperature
            return c2 * vs / np.log(c1 * vs**3 / (1.0e-5 * r) + 1.0)

        # Compute brightness temperature
        if units == 1:
            # Radiance units are Watts per square meter per steradian per micron
            BT = (bright_m(1.0e+4 / cwn[band - 20], rad) - tci[band - 20]) / tcs[band - 20]
        else:
            # Radiance units are milliWatts per square meter per steradian per wavenumber
            BT = (brite_m(cwn[band - 20], rad) - tci[band - 20]) / tcs[band - 20]

        return BT


class MOD06Granule():
    def __init__(self, mod_l1b_fullpath):
        self.fullpath = mod_l1b_fullpath  # Use the passed full path
    
    def get_data(self, fieldname,upscale_1km = False, scale_factor_flag = True):

        hdf_file = SD(self.fullpath)
        data = hdf_file.select(fieldname)
        if scale_factor_flag:
            scale_factor = data.attributes()['scale_factor']
            offset = data.attributes()['add_offset']
            data_value = data.get()
            valid_range = data.attributes()['valid_range']
            data_value = np.where((data_value >= valid_range[0]) & (data_value <= valid_range[1]), data_value, np.nan)
            data_value = (data_value - offset) *scale_factor 
        else:
            data_value = data.get()

        along_sampling = data.attributes()['Cell_Along_Swath_Sampling']
        if upscale_1km and along_sampling[-1] != 1:
            # Upscale the data from 5 km to 1 km by copying
            repeat_factor = 5  # Copy each 5 km pixel 5 times to reach 1 km resolution
            target_shape = (2030, 1354)
            # Resample the data by repeating each value in both dimensions
            data_value = np.repeat(np.repeat(data_value, repeat_factor, axis=0), repeat_factor, axis=1)
            current_shape = data_value.shape

            if current_shape[1] < target_shape[1]:
                # Pad the across-swath dimension to reach 1354 columns if needed
                padding_across = target_shape[1] - current_shape[1]
                data_value = np.pad(data_value, ((0, 0), (0, padding_across)), mode='edge')  # Pads the end by repeating values
            data_value = data_value[:target_shape[0], :target_shape[1]]
        return data_value
        
    def get_ctp(self):
        return self.get_data('cloud_top_pressure_1km')    

    def get_cttemp(self):
        return self.get_data('cloud_top_temperature_1km') 
    
    def get_opt(self):
        return self.get_data('Cloud_Optical_Thickness')
   
    def get_cer(self):
        return self.get_data('Cloud_Effective_Radius') 
    
    def get_ctm(self):
        return self.get_data('cloud_top_method_1km') 
                             
    def get_cphase(self):
        '''
        0: Cloud Mask unavaible, missing data
        1: Cloud Mask Clear or Probably Clear, or Pixel Restored to clear sky
        2: Liquid Water
        3: Ice
        4: Undertermined
        '''
        return self.get_data('Cloud_Phase_Optical_Properties')    
    
    def get_multilayer_withPH(self):
        '''
        0 indicates no cloud was detected, 
        1 indicates a monolayer cloud and values
        2–10 indicate the cumulative weight of the positive multilayer tests. 
        '''

        return self.get_data('Cloud_Multi_Layer_Flag')
    
    def get_multilayer_withoutPH(self, sum_flag = False):

        '''
        Use quality-assurance at 1km  to turn off the Pavolonis–Heidinger (PH) test
        0 indicates no cloud was detected, 
        1 indicates a monolayer cloud and values
        2 indicates multi-layer clouds without the PH test

        '''

        cloud_multi_layer_flag = self.get_data('Cloud_Multi_Layer_Flag')
        quality_assurance_1km  = self.get_data('Quality_Assurance_1km')
        byte5 = quality_assurance_1km[:, :, 5]


        # Create a bitmask to isolate Bits 0-3 (other multilayer tests)
         
        other_tests_mask = 0b00001111

        other_tests_result = byte5 & other_tests_mask
 
        
        multilayer_detected = other_tests_result != 0

        new_cloud_multi_layer_flag = np.copy(cloud_multi_layer_flag)

        # Set pixels with multilayer detection (excluding PH test) to 2 (or another appropriate value)
        if  sum_flag:
            positive_tests_count = (
                (other_tests_result & 0b0001) +
                ((other_tests_result & 0b0010) >> 1) +
                ((other_tests_result & 0b0100) >> 2) +
                ((other_tests_result & 0b1000) >> 3)
            )
            new_cloud_multi_layer_flag = positive_tests_count  # Multilayer cloud detected
        else:
            new_cloud_multi_layer_flag[~multilayer_detected] = 2  # Single-layer cloud
        # Set pixels without multilayer detection to 1 (single-layer cloud)
        new_cloud_multi_layer_flag[~multilayer_detected] = 1  # Single-layer cloud
        # Handle pixels where no cloud is detected using Cloud_Multi_Layer_Flag
        no_cloud = cloud_multi_layer_flag == 0
        new_cloud_multi_layer_flag[no_cloud] = 0  # No cloud detected

        return new_cloud_multi_layer_flag

    def get_cth(self):
        return self.get_data('cloud_top_height_1km',scale_factor_flag=False) 
    
    def get_emissivity(self,upscale_1km = True):
        return self.get_data('Cloud_Effective_Emissivity',upscale_1km = upscale_1km) 
          
class MOD03Granule():
    def __init__(self, mod_l1b_fullpath):
        self.fullpath = mod_l1b_fullpath  # Use the passed full path
    def get_data(self, fieldname,scale_factor_flag = True):
        hdf_file = SD(self.fullpath)
        data = hdf_file.select(fieldname)
        if scale_factor_flag:
            scale_factor = data.attributes()['scale_factor']
            data_value = data.get()
            return data_value*scale_factor
        else:
            return data[()]
    def get_sza(self):
        return self.get_data('SolarZenith')  
    def get_saa(self):
        return self.get_data('SolarAzimuth')
    def get_vza(self):
        return self.get_data('SensorZenith')
    def get_vaa(self):
        return self.get_data('SensorAzimuth')
    def get_lat(self):
        return self.get_data('Latitude',scale_factor_flag=False)
    def get_lon(self):
        return self.get_data('Longitude',scale_factor_flag=False)
    def get_landsea_mask(self,binary_output=True):
        '''
        0: Shallow ocean
        1: Land
        2: Coastline
        3: Shallow inland water
        4: Ephemeral water
        5: Deep inland water
        6: Moderate or continental ocean
        7: Deep ocean
        '''
        landsea_mask = self.get_data('Land/SeaMask', scale_factor_flag=False)
        if binary_output:
            # Initialize a binary mask with default value of 0 (non-land)
            binary_mask = np.zeros_like(landsea_mask, dtype=int)
            
            # Set binary values: 1 for land (value 1), 2 for coastline (value 2)
            binary_mask[landsea_mask == 1] = 1
            binary_mask[landsea_mask == 2] = 2
            
            return binary_mask
        
        # Return the original Land/SeaMask values if binary_output is False
        return landsea_mask


if __name__ == "__main__":
    a = MOD06Granule('/data/gdi/d/MOD06/2002/203/MOD06_L2.A2002203.2315.061.2017280021147.hdf')
    b = a.get_ctp(upscale_factor = 5)
    print(b.shape)