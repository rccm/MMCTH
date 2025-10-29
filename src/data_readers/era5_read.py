import xarray as xr
import numpy as np

class ERA5Single:
    def __init__(self, file_path):
        self.file_path = file_path
        self.dataset = self.load_data()

    def load_data(self):
        try:
            dataset = xr.open_dataset(self.file_path)
            return dataset
        except Exception as e:
            print(f"Failed to load data: {e}")
            return None

    def get_variable(self, var_name):
        if self.dataset is not None:
            if var_name in self.dataset.variables:
                return self.dataset[var_name]
            else:
                print(f"Variable {var_name} not found in the dataset.")
        return None

    def get_metadata(self):
        if self.dataset is not None:
            return self.dataset.attrs
        return None

    def get_coordinates(self):
        if self.dataset is not None:
            return {coord: self.dataset[coord] for coord in self.dataset.coords}
        return None

    def get_time(self):
        return self.get_variable('time')

    def get_longitude(self):
        return self.get_variable('longitude')

    def get_latitude(self):
        return self.get_variable('latitude')

    def get_d2m(self):
        return self.get_variable('d2m')

    def get_t2m(self):
        return self.get_variable('t2m')

    def get_sst(self):
        return self.get_variable('sst')

    def get_skt(self):
        '''Skin Temperature'''
        return self.get_variable('skt')
    
    def get_sp(self):
        '''Surface Pressure'''
        return self.get_variable('sp')
    
    def get_msl(self):
        '''Mean Sea-level Pressure'''
        return self.get_variable('msl')
    
    
import xarray as xr

class ERA5Multi:
    def __init__(self, file_path):
        self.file_path = file_path
        self.dataset = self.load_data()

    def load_data(self):
        try:
            dataset = xr.open_dataset(self.file_path)
            return dataset
        except Exception as e:
            print(f"Failed to load data: {e}")
            return None

    def get_variable(self, var_name):
        if self.dataset is not None:
            if var_name in self.dataset.variables:
                return self.dataset[var_name]
            else:
                print(f"Variable {var_name} not found in the dataset.")
        return None

    def get_metadata(self):
        if self.dataset is not None:
            return self.dataset.attrs
        return None

    def get_coordinates(self):
        if self.dataset is not None:
            return {coord: self.dataset[coord] for coord in self.dataset.coords}
        return None

    def get_time(self):
        return self.get_variable('time')

    def get_longitude(self):
        return self.get_variable('longitude')

    def get_latitude(self):
        return self.get_variable('latitude')

    def get_level(self):
        return self.get_variable('level')

    def get_temperature(self):
        return self.get_variable('t')

    def get_specific_humidity(self):
        return self.get_variable('q')

    def get_geopotential(self):
        return self.get_variable('z')
    
    def get_u(self):
        return self.get_variable('u')
    
    def get_v(self):
        return self.get_variable('v')
    
    def get_w(self):
        return self.get_variable('w')
