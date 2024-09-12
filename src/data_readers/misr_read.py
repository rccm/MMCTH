#!/usr/bin/env python3
import os
import pyhdf.SD as SD
from pyhdf.HDF import *
from pyhdf.V import *
import numpy as np

class MISRGranule:
    def __init__(self, input_path):
        if input_path.startswith('/'):
            # Assume input_path is the full file path
            self.file_path = input_path
        else:
            # Assume input_path is the path number
            path_no = int(input_path)
            self.file_path = f'/data/gdi/e/MISR/AGP/MISR_AM1_AGP_P{str(path_no).zfill(3)}_F01_24.hdf'
        self.data = None
        self.block_shifts = self.generate_block_shifts()
        if os.path.isfile(self.file_path):
            self.hdf_file = SD.SD(self.file_path)
        else:
            raise FileNotFoundError("File not found: " + self.file_path)
        
    def close_file(self):
        if hasattr(self, 'hdf_file'):
            self.hdf_file.end()

    def get_dataset(self, dataset_name):
        hdf_file = SD.SD(self.file_path)
        if hasattr(self, 'hdf_file'):
            return self.hdf_file.select(dataset_name)
        else:
            raise Exception("HDF file not open. Call open_file() first.")
        
    def read_data(self, dataset_name, block_range=[1,180]):
        hdf_file = SD.SD(self.file_path)
        dataset = hdf_file.select(dataset_name)
        data = dataset.get()   
        block_start = np.array(block_range[0]) - 1
        block_end = np.array(block_range[1])
        data = data[block_start:block_end,:,:]
        return data
    
    def get_metadata(self, dataset_name):
       
        dataset = self.get_dataset(dataset_name)
        metadata = dataset.attributes(full=1)
        return metadata
    
    def generate_block_shifts(self):
        # Generate the block shift dictionary based on block numbers
        # Hard coped from the meta data
        block_shifts = np.zeros((179),dtype='int')
        block_shifts[[  1,   3,   7,  12, 171, 174, 177]] = 16
        block_shifts[[ 26,  30,  33,  36,  38,  40,  42,  43,  45,  47,  48,  50,  51,
                52,  54,  55,  56,  57,  59,  60,  61,  62,  63,  64,  65,  66,
                67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
                80,  81,  82,  83,  85,  86,  87,  88,  89,  90,  91,  92,  93,
                94,  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107,
                108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 120, 121,
                122, 123, 124, 126, 127, 128, 130, 131, 133, 135, 136, 138, 140,
                143, 145, 148, 153]] = -16
        block_shifts[[84,95]] = -32
        return block_shifts
    
    def align_blks(self, full_data, block_range, f_11km=1):
        block_start = np.array(block_range[0])
        block_end = np.array(block_range[1])
        nblocks = block_end - block_start
        z, x, y = full_data.shape
        
        if block_end <= block_start + 1:
            return full_data
        
        cumblock_shift = np.cumsum(self.block_shifts[block_start-1:block_end] * f_11km)
        
        left_bound = np.min(np.append(0, cumblock_shift))
        right_bound = np.max(np.append(0, cumblock_shift))
        
        width = int(y + abs(left_bound) + abs(right_bound))
        height = int(x * nblocks)
        
        combined_data = np.full((height, width), np.nan)       
        for i in range(nblocks):
            start_x = i * x
            end_x = (i + 1) * x
            shift = cumblock_shift[i-1] if i > 0 else 0
            left_y = abs(left_bound) + shift
            right_y = left_y + y
            combined_data[start_x:end_x, left_y:right_y] = full_data[i, :, :]
        return combined_data
    
    def get_scale_factor(self, group_name, attribute_name):
      try:
          if hasattr(self, 'hdf_file'):
              hdf = self.hdf_file.file  # Access the HDF file object
  
              # Check if the group exists
              if group_name in hdf:
                  group = hdf[group_name]
                  print(group.attribute_name)
                  # Check if the attribute exists in the group
                  if attribute_name in group.attrs:
                      scale_factor = group.attrs[attribute_name]
                      return scale_factor
                  else:
                      print(f"Attribute '{attribute_name}' not found in group '{group_name}'.")
              else:
                  print(f"Group '{group_name}' not found in the HDF file.")
      except Exception as e:
          print(f"Error: {e}")
  
      return None
    
    def get_rad(self, band_name,block_range,f_11km=1,align = False):
        band_name = band_name.lower()
        band_to_dataset = {'red': ('Red Radiance/RDQI',0.03853199630975723), \
                           'nir': ('NIR Radiance/RDQI', 0.024670248851180077), \
                           'green': ('Green Radiance/RDQI',0.04653150960803032),\
                           'blue': ('Blue Radiance/RDQI', 0.047203224152326584)   }  
        band_info = band_to_dataset.get(band_name, None)
        dataset_name, scale_factor = band_info
        if dataset_name is None:
            raise ValueError('Wrong dataset name')
        sdn_data = self.read_data(dataset_name,block_range)
        sdn_data = sdn_data >> 2
        sdn_data = sdn_data.astype('float')
        sdn_data[sdn_data>= 16377] = np.nan               
        if align:
            sdn = self.align_blks(sdn_data,block_range,f_11km=f_11km)
            sdn = sdn.astype('float')
            rad_aligned = sdn*scale_factor
            rad_raw = sdn_data*scale_factor
            return rad_aligned, rad_raw
        return sdn_data*scale_factor
        
        
    
    def get_brf(self, band_name,block_range,f_11km=1, align = False):
        
        band_name = band_name.lower()
        band_to_dataset = {'red': 'RedConversionFactor', \
                           'nir': 'NIRConversionFactor',  \
                           'green':'GreenConversionFactor',\
                           'blue': 'BlueConversionFactor'} 
        brffactor_name = band_to_dataset.get(band_name, None)  

        rad_raw = self.get_rad(band_name,block_range,f_11km=f_11km) 
        brf_factor_raw = self.read_data(brffactor_name,block_range = block_range)
        ratio = int(rad_raw.shape[1] / brf_factor_raw.shape[1])
        brf_factor = np.repeat(brf_factor_raw, ratio, axis=1)
        brf_factor = np.repeat(brf_factor, ratio, axis=2)
        brf_raw = rad_raw * brf_factor
        if align:
            return  self.align_blks(brf_raw, block_range,f_11km=f_11km)
        else:
            return brf_raw
        
    def get_cth(self, quality_flag = True, upscale_factor = 1):
        cth = self.read_data('CloudTopHeight')
        start_blk, end_blk = self.get_block_range()
        if start_blk > 1:
            cth[0:start_blk+1,:,:] = -777  #including the starting block
        if end_blk < 180:
            cth[end_blk:,:,:] = -777 #including the ending block
        if  upscale_factor >1: 
            cth = np.repeat(np.repeat(cth, upscale_factor, axis=1), upscale_factor, axis=2)
        if quality_flag:
            stereoquality = self.read_data('StereoQualityIndicator')
            return self.align_blks(cth,[1,179]), self.align_blks(stereoquality,[1,179])
        else:
            return self.align_blks(cth,[1,179])
        
    def get_geo(self, upscale_factor = 1):
        latitude = self.read_data('GeoLatitude')
        longitude = self.read_data('GeoLongitude')
        return self.align_blks(latitude,[1,179]), self.align_blks(longitude,[1,179])
    
    def get_block_range(self):
        hdf = SD.SD(self.file_path)
        global_attrs = hdf.attributes()
       # Get all global attributes
        start_block = global_attrs.get('Start_block')
        end_block = global_attrs.get('End block')  # Adjusted based on your output
        return start_block,end_block

 
    

    

if __name__ == "__main__":
    
    misr_l1file = "/Users/gzhao1/Downloads/MISR_AM1_GRP_ELLIPSOID_GM_P002_O072787_AN_F03_0024.hdf"
    
    misr_l1b = MISRGranule('003')
    lat,lon = misr_l1b.get_geo()
    
    # misr_reader.open_file()
    # dataset_name = "NIR Radiance/RDQI"  # Replace with the actual dataset name
    # Now you can work with the data and metadata as needed
    # Don't forget to close the file when you're done
    #data=misr_l1b.align_blks(dataset_name,[63,65],f_11km=4)
    #data,_= misr_l1b.get_rad('blue',[66,75],f_11km=4)
    # brf = misr_l1b.get_brf('blue',[40,85],f_11km=4)
    # scale = misr_l1b.get_scale_factor("NIRBand","Scale_factor") 
    # bitmask = (1 << 2) - 1# Extract bits 2-15
    # radiance = (data >> 2) & bitmask
    # radiance[data >= 17000] = -9999.
     