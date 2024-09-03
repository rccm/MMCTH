import numpy as np 
import h5py

text_file_path = "/data/gdi/f/gzhao1/Database/tc_cloud_h5.list"

def extract_ref_cameras_from_file(hdf_file_path):
    try:
        with h5py.File(hdf_file_path, 'r') as hdf_file:
            # Navigate to the MISRCameras attribute
            attr_values = hdf_file['/Stereo_1.1_km/Grid Attributes/MISRReferenceCamera'][0][0]
            # Convert bytes to string and split into camera names
            cameras_string = attr_values.decode('utf-8')
            cameras = [cameras_string[i:i+2] for i in range(0, len(cameras_string), 2)]
            return cameras
    except Exception as e:
        print(f"Error processing file {hdf_file_path}: {e}")
        return None
    
def extract_cameras_from_file(hdf_file_path):
    try:
        with h5py.File(hdf_file_path, 'r') as hdf_file:
            # Navigate to the MISRCameras attribute
            attr_values = hdf_file['/Stereo_1.1_km/Grid Attributes/MISRCameras'][0][0]
            # Convert bytes to string and split into camera names
            cameras_string = attr_values.decode('utf-8')
            cameras = [cameras_string[i:i+2] for i in range(0, len(cameras_string), 2)]
            return cameras
    except Exception as e:
        print(f"Error processing file {hdf_file_path}: {e}")
        return None
    
with open(text_file_path, 'r') as file_list:
    hdf5_files = file_list.readlines()


with open('none_nadir_ref.list', 'w') as output_file:
    for hdf5_file in hdf5_files:
        hdf5_file = hdf5_file.strip()  # Remove any extra whitespace
        if hdf5_file == '/data/gdi/e/MISR/TCCloud/2022.07.01/MISR_AM1_TC_CLOUD_P178_O119864_F01_0001.h5':
            continue
        cameras = extract_cameras_from_file(hdf5_file)
        other_elements = [element for element in cameras if element not in ['Aa', 'An', 'Af']]
        if other_elements:
            output_file.write(f'{hdf5_file} {cameras[0]}\n')
            print(f'{hdf5_file} {cameras[0]}')  # Optional: you can remove this line if not needed

# Print out the cameras for each file
# for file, cameras in all_cameras.items():
#     print(f"File: {file}")
#     print(f"Cameras: {', '.join(cameras)}\n")
