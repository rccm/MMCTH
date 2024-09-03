import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import xarray as xr
import os

# Define new pressure levels and latitude range
new_pressure_levels = np.array([0.005, 0.0161, 0.0384, 0.0769, 0.137, 0.2244, 0.3454, 0.5064, 0.714, 0.9753])
latitudes = np.arange(-90, 90.25, 0.25)

def interpolate_file(file_path):
    # Read the dataset
    data = pd.read_csv(file_path, delim_whitespace=True)
    # Extract columns
    columns = ['z(km)', 'p(mb)', 'T(K)', 'air(cm-3)', 'o3(cm-3)', 'o2(cm-3)', 'h2o(cm-3)', 'co2(cm-3)', 'no2(cm-3)']
  
  # Read the dataset with proper column alignment
    data = pd.read_csv(file_path, delim_whitespace=True, names=columns, skiprows=1)
  
  # Extract columns
    line_pressure = data['p(mb)']
    z_km = data['z(km)']
    h2o_ori = data['h2o(cm-3)']
    
    air = data['air(cm-3)']
    t_k = data['T(K)']
    
    h2o = h2o_ori/air * 18.016 / 28.966 * 1e3
    print(h2o_ori)
    print(h2o)
    pressure = np.log(line_pressure)

    # Interpolation functions
    print(file_path)
    print(new_pressure_levels)
   
    
    interp_z_km = interp1d(pressure, z_km, kind='linear', bounds_error=False, fill_value=np.nan)
    interp_h2o = interp1d(pressure, h2o, kind='linear', bounds_error=False, fill_value=np.nan)
    interp_air = interp1d(pressure, air, kind='linear', bounds_error=False, fill_value=np.nan)
    interp_t_k = interp1d(pressure, t_k, kind='linear', bounds_error=False, fill_value=np.nan)
    
    # Interpolate data
    
    z_km_new = interp_z_km(np.log(new_pressure_levels))
    h2o_new = interp_h2o(np.log(new_pressure_levels))
    air_new = interp_air(np.log(new_pressure_levels))
    t_k_new = interp_t_k(np.log(new_pressure_levels))
    
    # Combine data into one array
    combined_data = np.vstack((new_pressure_levels, z_km_new, h2o_new, air_new, t_k_new)).T
    
    return combined_data

def process_files(file_paths):
    interpolated_data = {}
    for file_path in file_paths:
        filename = os.path.basename(file_path)
        # Label the datasets based on file names (e.g., afglms.dat -> afglms)
        label = filename.replace('.dat', '')
        interpolated_data[label] = interpolate_file(file_path)
    
    return interpolated_data

def save_to_hdf5(interpolated_data, output_file):
    with pd.HDFStore(output_file, mode='w') as store:
        for key, value in interpolated_data.items():
            df = pd.DataFrame(value, columns=['p(mb)', 'z(km)', 'h2o(cm-3)', 'air(cm-3)', 'T(K)'])
            store.put(key, df)

def load_from_hdf5(output_file):
    with pd.HDFStore(output_file, mode='r') as store:
        interpolated_data = {}
        for key in store.keys():
            label = key.strip('/')
            df = store.get(key)
            interpolated_data[label] = df.to_numpy()
    
    return interpolated_data

def get_profile_by_latitude_and_month(profiles, latitude, month):
    abs_lat = np.abs(latitude)
    if abs_lat < 30.0:
        return profiles['afglt']
    elif ((latitude >= 30.0 and month in (4, 5, 6, 7, 8, 9)) or 
          (latitude <= -30.0 and month in (1, 2, 3, 10, 11, 12))):
        if abs_lat >= 60.0:
            return profiles['afglss']
        else:
            return profiles['afglms']
    else:
        if abs_lat >= 60.0:
            return profiles['afglsw']
        else:
            return profiles['afglmw']

def create_seasonal_datasets(profiles):
    summer_spring = np.zeros((len(latitudes), len(new_pressure_levels), 5))
    fall_winter = np.zeros((len(latitudes), len(new_pressure_levels), 5))
    
    for i, lat in enumerate(latitudes):
        summer_spring[i, :, :] = get_profile_by_latitude_and_month(profiles, lat, 7)  # Example month in summer
        fall_winter[i, :, :] = get_profile_by_latitude_and_month(profiles, lat, 1)   # Example month in winter
    
    return summer_spring, fall_winter

# List of file paths
file_paths = [
    'afglms.dat',  # Replace with actual file paths
    'afglmw.dat',
    'afglss.dat',
    'afglsw.dat',
    'afglt.dat'
]
# Process files and get interpolated data
interpolated_data = process_files(file_paths)

# Save the interpolated data to a single HDF5 file
output_file = 'interpolated_profiles.h5'
save_to_hdf5(interpolated_data, output_file)

# Load the interpolated data from the HDF5 file
loaded_data = load_from_hdf5(output_file)

# Create seasonal datasets
summer_spring, fall_winter = create_seasonal_datasets(loaded_data)

# Save the seasonal datasets to a NetCDF file
ds = xr.Dataset(
    {
        "summer_spring": (["latitude", "pressure_level", "variable"], summer_spring),
        "fall_winter": (["latitude", "pressure_level", "variable"], fall_winter),
    },
    coords={
        "latitude": latitudes,
        "pressure_level": new_pressure_levels,
        "variable": ["p(mb)", "z(km)", "h2o(cm-3)", "air(cm-3)", "T(K)"]
    },
)

# ds.to_netcdf("seasonal_profiles.nc")

# # Example: Read the NetCDF file using xarray
# ds_loaded = xr.open_dataset("seasonal_profiles.nc")

