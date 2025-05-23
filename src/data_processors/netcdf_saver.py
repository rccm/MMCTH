import xarray as xr
import numpy as np
import json
from copy import deepcopy
   
class NetCDFSaver:
    def __init__(self, filename, logger):
        self.filename = filename
        self.logger = logger
        self.era5_P_levels = np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000], dtype=int)

    def save_mm_group(self, mm_variables, mm_geo, misr_variables, modis_variables, era5_variables, era5_geo, modis_metadata=None):
        """
        Save the provided variables into a NetCDF file with associated metadata and attributes,
        grouping them into the specified group (e.g., 'mm_variables').
        """
        # Create data dictionaries for each group
        data_dicts = self.create_group_data_dict(mm_variables, mm_geo, misr_variables, modis_variables, era5_variables, era5_geo)
        
        # Get attributes and metadata for all variables
        attributes, metadata, common_metadata = self.define_metadata_attributes()
        
        # If additional MODIS metadata is provided, add it to the global attributes
        if modis_metadata:
            attributes.update(modis_metadata)
        
        # Define encoding for all variables
        encoding = {
            'MMCloudTopHeight': {'dtype': 'int16', 'zlib': True, 'complevel': 5},
            'MMCloudTopPressure': {'dtype': 'float32', 'zlib': True, 'complevel': 5},
            'MMCloudEffectiveEmissivity': {'dtype': 'float32', 'zlib': True, 'complevel': 5},
            'MMCloudOpticalDepth': {'dtype': 'float32', 'zlib': True, 'complevel': 5},
            'MMFlag': {'dtype': 'int8', 'zlib': True, 'complevel': 5},
            # Add encoding for other variables if needed
        }

        try:
            first_group = True  # Use 'w' mode for the first group to create the file, 'a' for the rest
            for group_name, group_data_dict in data_dicts.items():
                # Create an empty xarray Dataset for the group
                dataset = xr.Dataset()
                
                # Loop over each variable in the group and add it to the dataset
                for field_name, field_info in group_data_dict.items():
                    data, dims, coords = field_info
        
                    # Start with common metadata and update with field-specific metadata if available
                    field_metadata = deepcopy(common_metadata) if common_metadata else {}
                    field_metadata.update(metadata.get(field_name, {}))
        
                    # Create an xarray DataArray with dimensions, coordinates, and metadata
                    data_array = xr.DataArray(
                        data, 
                        dims=dims, 
                        coords=coords, 
                        attrs=field_metadata
                    )
        
                    # Add the DataArray to the dataset
                    dataset[field_name] = data_array
        
                # Add global attributes to the dataset
                dataset.attrs = attributes
                
                # Define encoding for variables in this group
                group_encoding = {var_name: encoding.get(var_name, {'dtype': 'float32', 'zlib': True, 'complevel': 5}) for var_name in dataset.data_vars}
                    
                # Determine the file mode
                mode = 'w' if first_group else 'a'
                first_group = False  # Subsequent groups will use 'a'
                
                # Save the dataset to the NetCDF file under the specified group
                dataset.to_netcdf(
                    self.filename, 
                    mode=mode,
                    group=group_name,
                    encoding=group_encoding
                )
            
            # Log the success
            self.logger.info(f'Successfully saved data to {self.filename}')
        except Exception as e:
            # Log any errors that occur during the saving process
            self.logger.error(f'Error saving NetCDF file: {str(e)}')

    def create_group_data_dict(self, mm_variables, mm_geo, misr_variables, modis_variables, era5_variables, era5_geo):
        """
        Creates a grouped data dictionary for saving variables into groups in NetCDF4.
        """
        # Define Latitude and Longitude as 2D arrays for MODIS/MISR and ERA5 data
        lat = mm_geo['lat']
        lon = mm_geo['lon']
        era5_lat = era5_geo['era5_lat']
        era5_lon = era5_geo['era5_lon']
        x_cell_along_swath_1km = np.arange(lat[:, 0].shape[0])
        y_cell_cross_swath_1km = np.arange(lon[0, :].shape[0])
    
        # MM Cloud Variables Group
        mm_data_dict = {
            'MMCloudTopHeight': (
                mm_variables['cth'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MMCloudTopPressure': (
                mm_variables['ctp'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MMCloudEffectiveEmissivity': (
                mm_variables['emissivity'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MMCloudOpticalDepth': (
                mm_variables['opt'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MMFlag': (
                mm_variables['cflag'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
             'Latitude': (
                lat,
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'Longitude': (
                lon,
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
        }

        # MM Cloud Variables Group
        modis_data_dict = {
            'MODISCloudTopHeight': (
                modis_variables['cth'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MODISCloudTopPressure': (
                modis_variables['ctp'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MODISCloudEmissivity': (
                modis_variables['emissivity'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MODISCloudOpticalDepth': (
                modis_variables['opt'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MODISCloudPhase': (
                modis_variables['cphase'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MISRCloudTopHeights': (
                misr_variables['cphase'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),

        }
        misr_data_dict = {
            'MISRCloudTopHeights': (
                misr_variables['cth'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            # 'MISRCloudTopHeights_QA': (
            #     misr_variables['cth'],
            #     ['cell_along_swath_1km', 'cell_cross_swath_1km'],
            #     {'cell_along_swath_1km':  x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            # ),
        # Add other MISR variables if needed
        }
 
        era5_data_dict = {}
        for var_name, var_data in era5_variables.items():
            dims = var_data.shape
            if len(dims) == 3:
                era5_data_dict[var_name] = (
                    var_data,
                    ['z', 'era5_cell_along_swath_10km', 'era5_cell_cross_swath_10km'],
                    {'z': self.era5_P_levels,'era5_cell_along_swath_1km': np.arange(dims[1]), 'era5_cell_cross_swath_1km': np.arange(dims[2])},  # 1D coordinates
                )
            elif len(dims) == 2:
                era5_data_dict[var_name] = (
                    var_data,
                    ['era5_cell_along_swath_10km', 'era5_cell_cross_swath_10km'],
                    {'era5_cell_along_swath_10km': np.arange(dims[0]), 'era5_cell_cross_swath_1km': np.arange(dims[1])},
                )

       
        return {
            'MM_Variables': mm_data_dict,
            'MODIS_Variables': modis_data_dict,
            'MISR_Variables': misr_data_dict,
            'ERA5_Variables': era5_data_dict,
            # Add other groups like 'modis_variables', 'misr_variables', etc. if needed
        }

    def save_mm(
        self, 
        mm_variables, 
        mm_geo, 
        misr_variables=None, 
        modis_variables=None, 
        era5_variables=None, 
        modis_metadata=None
    ):
        """
        Save the provided variables into a NetCDF file with associated metadata and attributes.
        
        :param mm_variables: Dictionary of MM variables.
        :param mm_geo: Dictionary of geographical information.
        :param misr_variables: Dictionary of MISR variables (optional).
        :param modis_variables: Dictionary of MODIS variables (optional).
        :param era5_variables: Dictionary of ERA5 reanalysis variables (optional).
        :param modis_metadata: Dictionary of metadata extracted from MODIS granules (optional).
        """
        # Validate required inputs
        if not mm_variables or not mm_geo:
            raise ValueError("mm_variables and mm_geo are required parameters.")
        
        # Create data dictionary
        data_dict = self.create_data_dict(
            mm_variables=mm_variables, 
            mm_geo=mm_geo, 
            misr_variables=misr_variables, 
            modis_variables=modis_variables, 
            era5_variables=era5_variables
        )
        
        # Get attributes and metadata for all variables
        attributes, metadata, common_metadata = self.define_metadata_attributes()
        
        # If additional MODIS metadata is provided, add it to the global attributes
        if modis_metadata:
            attributes.update(modis_metadata)
        
        try:
            # Create an empty xarray Dataset to hold the variables
            dataset = xr.Dataset()

            # Loop over each variable in data_dict and add it to the dataset
            for field_name, field_info in data_dict.items():
                data, dims, coords = field_info

                # Start with common metadata and update with field-specific metadata if available
                field_metadata = deepcopy(common_metadata) if common_metadata else {}
                field_metadata.update(metadata.get(field_name, {}))

                # Create an xarray DataArray with dimensions, coordinates, and metadata
                data_array = xr.DataArray(
                    data, 
                    dims=dims, 
                    coords=coords, 
                    attrs=field_metadata
                )
                # Add the DataArray to the dataset
                dataset[field_name] = data_array

            # Add global attributes to the dataset
            dataset.attrs = attributes

            # Define base encoding
            base_encoding = {
                'zlib': True,
                'complevel': 5
            }

            # Define specific encodings for known variables
            specific_encodings = {
                'MMCloudTopHeight': {'dtype': 'int16'},
                'MMCloudTopPressure': {'dtype': 'float32'},
                'MMCloudEffectiveEmissivity': {'dtype': 'float32'},
                'MMCloudOpticalDepth': {'dtype': 'float32'},
                'MMFlag': {'dtype': 'int8'},
                'MODISCloudTopHeight': {'dtype': 'int16'},
                'MODISCloudTopPressure': {'dtype': 'float32'},
                'MODISCloudEmissivity': {'dtype': 'float32'},
                'MODISCloudOpticalDepth': {'dtype': 'float32'},
                'MODISCloudPhase': {'dtype': 'int8'},
                'MISRCloudTopHeights': {'dtype': 'int16'},
                'Latitude': {'dtype': 'float32'},
                'Longitude': {'dtype': 'float32'},
                'VZA': {'dtype': 'float32'},
                'ERA5Geopotential': {'dtype': 'float32'},
                'ERA5Temperature': {'dtype': 'float32'},
                'ERA5SpecificHumidity': {'dtype': 'float32'},
                'ERA5DewPoint2m': {'dtype': 'float32'},
                'ERA5SeaSurfaceTemperature': {'dtype': 'float32'},
                'ERA5Temperature2m': {'dtype': 'float32'},
                'ERA5SurfacePressure': {'dtype': 'float32'},
                'ERA5SkinTemperature': {'dtype': 'float32'},
                'ERA5MeanSeaLevelPressure': {'dtype': 'float32'}
            }

            # Dynamically build the encoding dictionary based on available data_dict keys
            encoding = {}
            for var in data_dict.keys():
                if var in specific_encodings:
                    encoding[var] = {**base_encoding, **specific_encodings[var]}
                else:
                    # Provide a default encoding if specific encoding is not defined
                    encoding[var] = base_encoding.copy()

            # Save the dataset to a NetCDF file with compression and data types
            print(self.filename)
            dataset.to_netcdf(self.filename, format='NETCDF4', encoding=encoding)
            # Log the success
            self.logger.info(f'Successfully saved data to {self.filename}')

        except Exception as e:
            # Log any errors that occur during the saving process
            self.logger.error(f'Error saving NetCDF file: {str(e)}')

    def create_data_dict(
        self, 
        mm_variables, 
        mm_geo, 
        misr_variables=None, 
        modis_variables=None, 
        era5_variables=None
    ):
        """
        Creates the data dictionary for all variables, including handling both 2D and 3D arrays for ERA5 variables,
        with separate 2D lat/lon arrays for MODIS/MISR and ERA5 data.
        
        :param mm_variables: Dictionary of MM variables.
        :param mm_geo: Dictionary of geographical information.
        :param misr_variables: Dictionary of MISR variables (optional).
        :param modis_variables: Dictionary of MODIS variables (optional).
        :param era5_variables: Dictionary of ERA5 reanalysis variables (optional).
        :return: Dictionary suitable for creating an xarray Dataset.
        """
        # Validate required inputs
        if not mm_variables or not mm_geo:
            raise ValueError("mm_variables and mm_geo are required parameters.")

        required_mm_keys = ['cth', 'ctp', 'emissivity', 'opt', 'cflag']
        for key in required_mm_keys:
            if key not in mm_variables:
                raise KeyError(f"Missing required MM variable '{key}' in mm_variables.")
        
        required_mm_geo_keys = ['lat', 'lon', 'vza']
        for key in required_mm_geo_keys:
            if key not in mm_geo:
                raise KeyError(f"Missing required geographical key '{key}' in mm_geo.")

        # Define Latitude and Longitude as 2D arrays for MODIS/MISR and ERA5 data
        lat = mm_geo['lat']
        lon = mm_geo['lon']

        x_cell_along_swath_1km = np.arange(lat.shape[0])
        y_cell_cross_swath_1km = np.arange(lon.shape[1])

        # Initialize data_dict with MM variables and geolocation
        data_dict = {
            'MMCloudTopHeight': (
                mm_variables['cth'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MMCloudTopPressure': (
                mm_variables['ctp'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MMCloudEffectiveEmissivity': (
                mm_variables['emissivity'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MMCloudOpticalDepth': (
                mm_variables['opt'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'MMFlag': (
                mm_variables['cflag'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'Latitude': (
                lat,
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'Longitude': (
                lon,
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
            'VZA': (
                mm_geo['vza'],
                ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
            ),
        }

        # Add MODIS variables if available
        if modis_variables:
            modis_mapping = {
                'cth': 'MODISCloudTopHeight',
                'ctp': 'MODISCloudTopPressure',
                'emissivity': 'MODISCloudEmissivity',
                'opt': 'MODISCloudOpticalDepth',
                'cphase': 'MODISCloudPhase'
            }

            for key, new_var_name in modis_mapping.items():
                if key in modis_variables:
                    data_dict[new_var_name] = (
                        modis_variables[key],
                        ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                        {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
                    )
                    self.logger.debug(f"Added MODIS variable '{new_var_name}'")
                else:
                    self.logger.warning(f"MODIS variable '{key}' is missing and will be skipped.")

        # Add MISR variables if available
        if misr_variables:
            misr_mapping = {
                'misrcth': 'MISRCloudTopHeights'
            }

            for key, new_var_name in misr_mapping.items():
                if key in misr_variables:
                    data_dict[new_var_name] = (
                        misr_variables[key],
                        ['cell_along_swath_1km', 'cell_cross_swath_1km'],
                        {'cell_along_swath_1km': x_cell_along_swath_1km, 'cell_cross_swath_1km': y_cell_cross_swath_1km}
                    )
                    self.logger.debug(f"Added MISR variable '{new_var_name}'")
                else:
                    self.logger.warning(f"MISR variable '{key}' is missing and will be skipped.")

        # Automatically handle ERA5 variables (both 2D and 3D) with separate lat/lon for ERA5 (2D arrays)
        if era5_variables:
            variable_mapping = {
                'geopotential_org': 'ERA5Geopotential',
                'temperature_org': 'ERA5Temperature',
                'specific_humidity_org': 'ERA5SpecificHumidity',
                'dew2m': 'ERA5DewPoint2m',
                'sst': 'ERA5SeaSurfaceTemperature',
                'temp2m': 'ERA5Temperature2m',
                'surface_pressure': 'ERA5SurfacePressure',
                'skint': 'ERA5SkinTemperature',
                'msp': 'ERA5MeanSeaLevelPressure'
            }

            for var_name, var_data in era5_variables.items():
                new_var_name = variable_mapping.get(var_name, var_name)
                dims = var_data.shape

                if len(dims) == 3:
                    # Assuming the first dimension is 'z'
                    data_dict[new_var_name] = (
                        var_data,
                        ['z', 'era5_cell_along_swath_10km', 'era5_cell_cross_swath_10km'],
                        {
                            'era5_cell_along_swath_10km': np.arange(dims[1]), 
                            'era5_cell_cross_swath_10km': np.arange(dims[2])
                        },
                    )
                    self.logger.debug(f"Added 3D ERA5 variable '{new_var_name}' with dimensions {dims}")
                elif len(dims) == 2:
                    data_dict[new_var_name] = (
                        var_data,
                        ['era5_cell_along_swath_10km', 'era5_cell_cross_swath_10km'],
                        {
                            'era5_cell_along_swath_10km': np.arange(dims[0]), 
                            'era5_cell_cross_swath_10km': np.arange(dims[1])
                        },
                    )
                    self.logger.debug(f"Added 2D ERA5 variable '{new_var_name}' with dimensions {dims}")
                else:
                    self.logger.warning(f"Skipping ERA5 variable '{var_name}' with unsupported dimensions: {dims}")

        return data_dict        

    

    def define_metadata_attributes(self):
        """
        Defines the metadata attributes for the NetCDF file, including handling both 2D and 3D ERA5 variables
        and separate lat/lon for ERA5 data.
        
        :return: attributes, metadata, and common_metadata dictionaries.
        """
        
        attributes = {
            'title': 'MODIS/MISR and ERA5 Cloud Top Height, Pressure, and Emissivity Data',
            'institution': 'University of Illinois at Urbana-Champaign',
            'source': 'MISR, MODIS, ERA5 Data',
            'description': 'This file contains cloud top height, pressure, emissivity, and phase data from MODIS and MISR, along with ERA5 reanalysis data for the same scene.',
            'history': 'Generated by custom pipeline for combining MODIS, MISR, and ERA5 data',
            'references': 'Data generated using NASA MISR and MODIS instruments, and ECMWF ERA5 data'
        }

        common_metadata = {
        }

        metadata = {
            # MMCloud fields
            'MMCloudTopHeight': {
                'units': 'm',
                'valid_range': [0, 20000],
                '_FillValue': -9999
            },
            'MMCloudTopPressure': {
                'units': 'hPa',
                'valid_range': [100, 1100],
                '_FillValue': -9999
            },
            'MMCloudEffectiveEmissivity': {
                'units': 'n/a',
                'valid_range': [0, 1],
                '_FillValue': -9999
            },
            'MMCloudOpticalDepth': {
                'units': 'n/a',
                'valid_range': [0, 100],
                '_FillValue': -9999
            },
            'MMFlag': {
                'units': 'n/a',
                'flag_values': [0, 1, 2, 3],
                '_FillValue': -9999
            },

            # MODISCloud fields
            'MODISCloudTopHeight': {
                'units': 'm',
                'valid_range': [0, 20000],
                '_FillValue': -9999
            },
            'MODISCloudTopPressure': {
                'units': 'hPa',
                'valid_range': [100, 1100],
                '_FillValue': -9999
            },
            'MODISCloudEmissivity': {
                'units': 'n/a',
                'valid_range': [0, 1],
                '_FillValue': -9999
            },
            'MODISCloudOpticalDepth': {
                'units': 'n/a',
                'valid_range': [0, 100],
                '_FillValue': -9999
            },
            'MODISCloudPhase': {
                'units': 'n/a',
                'flag_values': [0, 1, 2, 3],
                '_FillValue': 0
            },

            # MISRCloud fields
            'MISRCloudTopHeights': {
                'units': 'm',
                'valid_range': [0, 20000],
                '_FillValue': -9999
            },

            # Geolocation fields for MODIS/MISR (now using 2D lat/lon)
            'Latitude': {
                'units': 'degrees_north',
                'valid_range': [-90, 90],
                '_FillValue': -9999
            },
            'Longitude': {
                'units': 'degrees_east',
                'valid_range': [-180, 180],
                '_FillValue': -9999
            },
            'VZA': {
                'units': 'degrees',
                'valid_range': [0, 90],
                '_FillValue': -9999
            },

            # ERA5 fields with 2D lat/lon arrays
            'ERA5Geopotential': {
                'units': 'm2/s2',
                'valid_range': [-50000, 50000],
                '_FillValue': -9999
            },
            'ERA5Temperature': {
                'units': 'K',
                'valid_range': [150, 350],
                '_FillValue': -9999
            },
            'ERA5SurfacePressure': {
                'units': 'hPa',
                'valid_range': [500, 1100],
                '_FillValue': -9999
            },
            'ERA5WindSpeed': {
                'units': 'm/s',
                'valid_range': [0, 100],
                '_FillValue': -9999
            },
            'ERA5Humidity': {
                'units': '%',
                'valid_range': [0, 100],
                '_FillValue': -9999
            }
        }

        return attributes, metadata, common_metadata

    def save_metadata_to_json(self, metadata, filename='metadata.json'):
        with open(filename, 'w') as json_file:
            json.dump(metadata, json_file, indent=4)

    def load_metadata_from_json(self, filename='metadata.json'):
        with open(filename, 'r') as json_file:
            return json.load(json_file)

    def save_inter(self, bands_BT=None, misr_cth=None, mod_geo=None, mod06=None, era5_variables_misrswath=None, era5_lat_lon=None):
        try:
            # Ensure mod_geo exists
            if mod_geo is None:
                raise ValueError("mod_geo must be provided.")
            
            # Extract MODIS latitude and longitude (original dimensions)
            mod_lat = np.array(mod_geo['lat'])
            mod_lon = np.array(mod_geo['lon'])

            cell_along_swath_1km = mod_lat.shape[0]
            cell_cross_swath_1km = mod_lon.shape[1]

            # Prepare bands_BT arrays
            if bands_BT:
                bands_BT_arrays = np.array([v for v in bands_BT.values()])
                bands_BT_keys = list(bands_BT.keys())
            else:
                bands_BT_arrays = None

            # Prepare mod06 data (original dimensions)
            if mod06:
                mod06 = {k: np.array(v) for k, v in mod06.items()} if isinstance(mod06, dict) else {}

            # Prepare ERA5 data (block-averaged)
            era5_variables_misrswath_2d = {k: np.array(v) for k, v in era5_variables_misrswath.items() if v.ndim == 2} if era5_variables_misrswath else {}
            era5_variables_misrswath_3d = {k: np.array(v) for k, v in era5_variables_misrswath.items() if v.ndim == 3} if era5_variables_misrswath else {}

            # Create dataset with optional variables
            data_vars = {}
            
            if misr_cth:
                data_vars["misr_cth"] = (["latitude", "longitude"], misr_cth["misrcth"])
                data_vars["misr_cth_qa"] = (["latitude", "longitude"], misr_cth["misrcth_qa"])

            if bands_BT_arrays is not None:
                data_vars["bands_BT"] = (["band", "latitude", "longitude"], bands_BT_arrays)

            if mod06:
                data_vars["mod06_product"] = (["product", "latitude", "longitude"], np.array(list(mod06.values())))

            # Add ERA5 variables (block-averaged)
            data_vars.update({k: (["latitude_block", "longitude_block"], v) for k, v in era5_variables_misrswath_2d.items()})
            data_vars.update({k: (["level", "latitude_block", "longitude_block"], v) for k, v in era5_variables_misrswath_3d.items()})

            # Add geographic data for ERA5 (center pixel of each block)
            if era5_lat_lon:
                data_vars["era5_latitude"] = (["latitude_block", "longitude_block"], era5_lat_lon['lat'])
                data_vars["era5_longitude"] = (["latitude_block", "longitude_block"], era5_lat_lon['lon'])

            # Add MODIS geographic data (original dimensions)
            data_vars["landsea"] = (["latitude", "longitude"], mod_geo['landsea'])
            data_vars["latitude"] = (["latitude"], cell_along_swath_1km)  # Lat for each row (assuming uniform for each row)
            data_vars["longitude"] = (["longitude"], cell_cross_swath_1km)  # Lon for each column (assuming uniform for each column)
            data_vars["vza"] = (["latitude", "longitude"], mod_geo['vza'])

            # Create a dataset
            ds = xr.Dataset(
                data_vars,
                coords={
                    "latitude": (["latitude"], cell_along_swath_1km),
                    "longitude": (["longitude"],cell_cross_swath_1km),
                    "latitude_block": (["latitude_block"], era5_lat_lon['lat'][:, 0]),
                    "longitude_block": (["longitude_block"], era5_lat_lon['lon'][0, :]),
                    "band": bands_BT_keys if bands_BT else [],
                    "product": list(mod06.keys()) if mod06 else [],
                    "variable_2d": list(era5_variables_misrswath_2d.keys()) if era5_variables_misrswath_2d else [],
                    "variable_3d": list(era5_variables_misrswath_3d.keys()) if era5_variables_misrswath_3d else [],
                    "level": np.arange(1, 102)  # Assuming pressure levels are from 1 to 101
                }
            )

            # Add x and y as attributes
            ds.attrs["x"] = np.arange(cell_along_swath_1km).tolist()
            ds.attrs["y"] = np.arange(cell_cross_swath_1km).tolist()
            ds.attrs["x_block"] = np.arange(era5_lat_lon['lat'].shape[0]).tolist()
            ds.attrs["y_block"] = np.arange(era5_lat_lon['lon'].shape[1]).tolist()

            # Add metadata
            ds.attrs["description"] = "Processed MISR, MODIS (original dimensions), and ERA5 (block-averaged) data"
            ds.attrs["history"] = "Created " + str(np.datetime64('now', 's'))
            ds.attrs["source"] = "MISR, MODIS, and ERA5 data"

            ds["latitude"].attrs["units"] = "degrees_north"
            ds["longitude"].attrs["units"] = "degrees_east"
            ds["era5_latitude"].attrs["units"] = "degrees_north"
            ds["era5_longitude"].attrs["units"] = "degrees_east"
            if "misr_cth" in ds:
                ds["misr_cth"].attrs["units"] = "meters"
            for var in era5_variables_misrswath_2d:
                ds[var].attrs["units"] = "Unknown"
            for var in era5_variables_misrswath_3d:
                ds[var].attrs["units"] = "Unknown"
            if "bands_BT" in ds:
                ds["bands_BT"].attrs["units"] = "Kelvin"
            if "mod06_product" in ds:
                ds["mod06_product"].attrs["units"] = "Unknown"
            ds["landsea"].attrs["units"] = "Unknown"
            ds["vza"].attrs["units"] = "degrees"

            # Save to NetCDF file
            comp = dict(zlib=True, complevel=8)
            encoding = {var: comp for var in ds.data_vars}

            # Save to NetCDF file with compression
            ds.to_netcdf(self.filename, encoding=encoding)
            self.logger.info(f"Successfully saved data to {self.filename}")
        except Exception as e:
            self.logger.error(f"Failed to save data to {self.filename}: {e}", exc_info=True)
    
    def save_org(self, bands_BT=None, misr_cth=None, mod_geo=None, mod06=None, era5_variables_misrswath=None, outputfile = 'test'):
        try:
            # Ensure mod_geo exists
            if mod_geo is None:
                raise ValueError("mod_geo must be provided.")
            
            # Extract latitude and longitude
            mod_lat = np.array(mod_geo['lat'])
            mod_lon = np.array(mod_geo['lon'])
            
            cell_along_swath_1km = mod_lat.shape[0]
            cell_cross_swath_1km = mod_lon.shape[1]
            mod_lat_flat = mod_geo['lat'][:,0] # Assuming lat is the same for each column
            mod_lon_flat = mod_geo['lon'][0,:]  # Assuming lon is the same for each row

            # Prepare bands_BT arrays
            if bands_BT:
                bands_BT_arrays = np.array([v for v in bands_BT.values()])
                bands_BT_keys = list(bands_BT.keys())
            else:
                bands_BT_arrays = None

            # Prepare mod06 data
            if mod06:
                mod06 = {k: np.array(v) for k, v in mod06.items()} if isinstance(mod06, dict) else {}
            
            # Prepare ERA5 data
            era5_variables_misrswath_2d = {k: np.array(v) for k, v in era5_variables_misrswath.items() if v.ndim == 2} if era5_variables_misrswath else {}
            era5_variables_misrswath_3d = {k: np.array(v) for k, v in era5_variables_misrswath.items() if k.endswith("exp") and v.ndim == 3} if era5_variables_misrswath else {}
            # Create dataset with optional variables
            data_vars = {}
            
            if misr_cth:
                data_vars["misr_cth"] = (["lat", "lon"], misr_cth["misrcth"])
                data_vars["misr_cth_qa"] = (["lat", "lon"], misr_cth["misrcth_qa"])

            if bands_BT_arrays is not None:
                data_vars["bands_BT"] = (["band", "lat", "lon"], bands_BT_arrays)

            if mod06:
                data_vars["mod06_product"] = (["product", "lat", "lon"], np.array(list(mod06.values())))

            # Add ERA5 variables
            data_vars.update({k: (["lat", "lon"], v) for k, v in era5_variables_misrswath_2d.items()})
            data_vars.update({k: (["level", "lat", "lon"], v) for k, v in era5_variables_misrswath_3d.items()})
            
            # Add geographic data (always present)
            data_vars["landsea"] = (["lat", "lon"], mod_geo['landsea'])
            data_vars["vza"] = (["lat", "lon"], mod_geo['vza'])

            # Create a dataset
            ds = xr.Dataset(
                data_vars,
                coords={
                    "latitude": (["lat"], mod_lat_flat),
                    "longitude": (["lon"], mod_lon_flat),
                    "band": bands_BT_keys if bands_BT else [],
                    "product": list(mod06.keys()) if mod06 else [],
                    "variable_2d": list(era5_variables_misrswath_2d.keys()) if era5_variables_misrswath_2d else [],
                    "variable_3d": list(era5_variables_misrswath_3d.keys()) if era5_variables_misrswath_3d else [],
                    "level": np.arange(1, 102)  # Assuming pressure levels are from 1 to 101
                }
            )

            # Add metadata
            ds.attrs["description"] = "Processed MISR, MODIS, and ERA5 data"
            ds.attrs["history"] = "Created " + str(np.datetime64('now', 's'))
            ds.attrs["source"] = "MISR, MODIS, and ERA5 data"

            ds["latitude"].attrs["units"] = "degrees_north"
            ds["longitude"].attrs["units"] = "degrees_east"
            if "misr_cth" in ds:
                ds["misr_cth"].attrs["units"] = "meters"
            for var in era5_variables_misrswath_2d:
                ds[var].attrs["units"] = "Unknown"
            for var in era5_variables_misrswath_3d:
                ds[var].attrs["units"] = "Unknown"
            if "bands_BT" in ds:
                ds["bands_BT"].attrs["units"] = "Kelvin"
            if "mod06_product" in ds:
                ds["mod06_product"].attrs["units"] = "Unknown"
            ds["landsea"].attrs["units"] = "Unknown"
            ds["vza"].attrs["units"] = "degrees"

            # Save to NetCDF file
            comp = dict(zlib=True, complevel=8)
            encoding = {var: comp for var in ds.data_vars}

            # Save to NetCDF file with compression

            ds.to_netcdf(outputfile, encoding=encoding)
            self.logger.info(f"Successfully saved data to {self.filename}")
        except Exception as e:
            self.logger.error(f"Failed to save data to {self.filename}: {e}", exc_info=True)