# test_netcdf_saver.py
import pytest
import numpy as np
import xarray as xr
import logging
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from src.data_processors.netcdf_saver import NetCDFSaver

@pytest.fixture
def logger():
    logging.basicConfig(level=logging.INFO)
    return logging.getLogger("NetCDF Test Logger")

@pytest.fixture
def sample_data():
    return {
        'mm_variables': {
            'cth': np.random.randint(0, 20000, size=(10, 10)),
            'ctp': np.random.uniform(100, 1100, size=(10, 10)),
            'emissivity': np.random.uniform(0, 1, size=(10, 10)),
            'opt': np.random.uniform(0, 100, size=(10, 10)),
            'cflag': np.random.randint(0, 4, size=(10, 10))
        },
        'mm_geo': {
            'lat': np.linspace(-90, 90, 10).reshape(10, 1).repeat(10, axis=1),
            'lon': np.linspace(-180, 180, 10).reshape(1, 10).repeat(10, axis=0),
            'vza': np.random.uniform(0, 90, size=(10, 10))
        },
        'misr_variables': {
            'cphase': np.random.randint(0, 4, size=(10, 10))
        },
        'modis_variables': {
            'cth': np.random.randint(0, 20000, size=(10, 10)),
            'ctp': np.random.uniform(100, 1100, size=(10, 10)),
            'emissivity': np.random.uniform(0, 1, size=(10, 10)),
            'opt': np.random.uniform(0, 100, size=(10, 10)),
            'cphase': np.random.randint(0, 4, size=(10, 10))
        },
        'era5_variables': {
            'ERA5Geopotential': np.random.uniform(-50000, 50000, size=(37，100, 100)),
            'ERA5Temperature': np.random.uniform(150, 350, size=(37，100, 100)),
            'ERA5SurfacePressure': np.random.uniform(500, 1100, size=(37，100, 100)),
            'ERA5Humidity': np.random.uniform(100, 100, size=(100, 100))
        },
        'era5_geo': {
            'era5_lat': np.linspace(-90, 90, 5).reshape(5, 1).repeat(5, axis=1),
            'era5_lon': np.linspace(-180, 180, 5).reshape(1, 5).repeat(5, axis=0)
        }
    }

def test_save_mm(logger, sample_data):
    # Define test filename in temporary directory
    filename =  "test.nc"
    saver = NetCDFSaver(filename=filename, logger=logger)

    # Run the save_mm method
    saver.save_mm_group(
        mm_variables=sample_data['mm_variables'],
        mm_geo=sample_data['mm_geo'],
        misr_variables=sample_data['misr_variables'],
        modis_variables=sample_data['modis_variables'],
        era5_variables=sample_data['era5_variables'],
    )

    # Load and check if file was created and data is correct
    dataset = xr.open_dataset(filename,group='mm_variables')
    assert 'MMCloudTopHeight' in dataset, "MMCloudTopHeight not saved"
    assert 'Latitude' in dataset, "Latitude not saved"
    assert 'Longitude' in dataset, "Longitude not saved"
    assert dataset['MMCloudTopHeight'].shape == (10, 10), "Incorrect shape for MMCloudTopHeight"
