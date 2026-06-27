       SUBROUTINE GET_SPEC_INFO (PROGRAMNAME, INPUTNAME,
     $        LONG_NAME, DIMLIST, UNITS_STR, VALID_RANGE, FILL_VALUE,
     $        NUMBERTYPE, SCALE_FACTOR, ADD_OFFSET, PARAM_TYPE,
     $        CELL_ACROSS_SWATH_SAMPLING, CELL_ALONG_SWATH_SAMPLING,
     $        DESCRIPTION, DESCCOUNT, TITLE, TITLECOUNT,
     $        HISTORY, HISTORYCOUNT, GEOLOCATION_POINTER)

       implicit none
       save

c---------------------------------------------------------------------
C!F77
c
c!DESCRIPTION:
c   Subroutine which extracts the correct file spec information
c   and creates an HDF output file for MOD35, MOD07, or MOD06.
c
c!INPUT PARAMETERS:
c   ProgramName   Ouput hdf file name
c   InputName     Name of Input SDS
c
c!OUTPUT PARAMETERS:
c   long_name                   HDF long name of SDS
c   dimlist                     Dimension list for SDS
c   units_str                   Unit of current SDS
c   valid_range                 Valid range of current SDS
c   Fill_Value                  Missing data value for given SDS
c   NumberType                  Type of data for current SDS
c   scale_factor                Scale factor for current SDS
c   add_offset                  Add/offset values for current SDS
c   Param_Type                  Parameter type of current SDS
c   Cell_Across_Swath_Sampling  Sampling for current SDS
c   Cell_Along_Swath_Sampling   Sampling for current SDS
c   description                 SDS description
c   DescCount                   Number of lines to be included in description
c   title                       SDS title
c   TitleCount                  Number of title lines
c   history                     SDS history summary information
c   HistoryCount                Number of history lines
c   Geolocation_Pointer         SDS geolocation pointer
c
c!REVISION HISTORY:
c   This FORTRAN source file is created by spec2fort.c
c   spec2fort.c was designed by Walter.Wolf@ssec.wisc.edu
c
c!TEAM-UNIQUE HEADER:
c
c!REFERENCES AND CREDITS:
c
c!END
c---------------------------------------------------------------------

       include 'mapi.inc'

       include 'hdf.inc'

       include 'dffunc.inc'

       character*(*) ProgramName, InputName, Geolocation_Pointer
       character*(*) long_name, units_str, Param_Type, dimlist
       double precision valid_range(*), Fill_Value
       double precision scale_factor, add_offset
       integer Cell_Across_Swath_Sampling(*)
       integer Cell_Along_Swath_Sampling(*)
       character*(*) title, history
       character*(*) description(*)
       integer TitleCount, HistoryCount, DescCount, NumberType

       integer unlimited

       unlimited = 0
       TitleCount = 0
       HistoryCount = 0
       if (ProgramName .EQ. "MOD04_L2") then
          
          
          if (InputName .EQ. "MODIS_Band_Land") then
             NumberType = DFNT_INT32
             DimList = 
     $ "MODIS_Band_Land"
             long_name = 
     $ "Center Wavelengths of MODIS Bands Used in Land Retrieval Algorithms" 
             units_str = "Nanometers" 
          End If

          if (InputName .EQ. "MODIS_Band_Ocean") then
             NumberType = DFNT_INT32
             DimList = 
     $ "MODIS_Band_Ocean"
             long_name = 
     $ "Center Wavelengths of MODIS Bands Used in Ocean Retrieval Algorithms" 
             units_str = "Nanometers" 
          End If

          if (InputName .EQ. "Solution_1_Land") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Solution_1_Land"
             long_name = 
     $ "Central Wavelength of MODIS Bands Used in Continental Model Retrieval" 
             units_str = "Nanometers" 
          End If

          if (InputName .EQ. "Solution_2_Land") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Solution_2_Land"
             long_name = 
     $ "Central Wavelength of MODIS Bands Used in Corrected Model Retrieval" 
             units_str = "Nanometers" 
          End If

          if (InputName .EQ. "Solution_3_Land") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Solution_3_Land"
             long_name = 
     $ "Central Wavelength of MODIS Bands Used in Scatter Plot Solution" 
             units_str = "Nanometers" 
          End If

          if (InputName .EQ. "Solution_Ocean") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Solution_Ocean"
             long_name = 
     $ "Index of Ocean Solution Types  1 - Best 2 - Average" 
             units_str = "None" 
          End If

          if (InputName .EQ. "Solution_Index") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Solution_Index"
             long_name = 
     $ "Solution Index for Ocean Best Small (1-4) and Large (5-9) Modes" 
             units_str = "None" 
          End If

          if (InputName .EQ. "Longitude") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Geodetic Longitude" 
             units_str = "Degrees_east" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Geolocation data not applicable" 
             Fill_Value = -999. 
             valid_range(1) = -180. 
             valid_range(2) = 180. 
          End If

          if (InputName .EQ. "Latitude") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Geodetic Latitude" 
             units_str = "Degrees_north" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Geolocation data not applicable" 
             Fill_Value = -999. 
             valid_range(1) = -90. 
             valid_range(2) = 90. 
          End If

          if (InputName .EQ. "Scan_Start_Time") then
             NumberType = DFNT_FLOAT64
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "TAI Time at Start of Scan replicated across the swath" 
             units_str = "Seconds since 1993-1-1 00:00:00.0 0" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -999 
             valid_range(1) = 0 
             valid_range(2) = 3155800000 
          End If

          if (InputName .EQ. "Solar_Zenith") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Solar Zenith Angle, Cell to Sun" 
             units_str = "Degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 18000 
          End If

          if (InputName .EQ. "Solar_Azimuth") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Solar_Azimuth Angle, Cell to Sun" 
             units_str = "Degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
          End If

          if (InputName .EQ. "Sensor_Zenith") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Sensor_Zenith Angle, Cell to Sensor" 
             units_str = "Degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 18000 
          End If

          if (InputName .EQ. "Sensor_Azimuth") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Sensor_Azimuth Angle, Cell to Sensor" 
             units_str = "Degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
          End If

          if (InputName .EQ. "Scattering_Angle") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Scattering Angle" 
             units_str = "Degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 18000 
          End If
          
           if (InputName .EQ. "Glint_Angle") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Glint Angle" 
             units_str = "Degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 18000 
          End If

          if (InputName .EQ. "Land_Ocean_Quality_Flag") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             valid_range(1) = 0 
             valid_range(2) = 3 
             Fill_Value = -9999 
             long_name = 
     $ "Quality Flag for Land and ocean Aerosol retreivals 0= bad  1 " //
     $ " = Marginal 2= Good 3=Very Good)"
             units_str = "None" 
             scale_factor = 1 
             add_offset = 0 
             Param_Type = "Output" 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Land_sea_Flag") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             valid_range(1) = 0 
             valid_range(2) = 1 
             Fill_Value = -9999 
             long_name = 
     $ "Land_sea_Flag(based on MOD03 Landsea mask 0 = Ocean, 1 = Land and Ephemeral water 2 =Coastal)" 
             units_str = "None" 
             scale_factor = 1 
             add_offset = 0 
             Param_Type = "Output" 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If 

          if (InputName .EQ. "Wind_Speed_Ncep_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             valid_range(1) = 0 
             valid_range(2) = 8000
             Fill_Value = -9999 
             long_name = 
     $ "Wind Speed based on NCEP reanalysis for Ocean " 
             units_str = "Meters" 
             scale_factor = 0.01 
             add_offset = 0 
             Param_Type = "Output" 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Optical_Depth_Land_And_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "AOT at 0.55 micron for both ocean (Average) (Quality flag=1,2,3) and land (corrected) " //
     $ "(Quality flag=3)"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Image_Optical_Depth_Land_And_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "AOT at 0.55 micron for both ocean (Average) and land (corrected) " //
     $ " with all quality data (Quality flag=0,1,2,3)"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Aerosol_Type_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Aerosol Type: 1 = Continental, 2 = Moderate Absorption Fine, 3 = Strong Absorption Fine,"//
     $  "4 = Weak Absorption Fine, 5 = Dust Coarse"
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 4 
          End If

          if (InputName .EQ. "Fitting_Error_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Spectral Fitting error for inversion over land" 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 1000 
          End If

          if (InputName .EQ. "Surface_Reflectance_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_2_Land"
             long_name = 
     $ "Estimated Surface Reflectance at 0.47,0.66 and 2.13micron" 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Corrected_Optical_Depth_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_3_Land"
             long_name = 
     $ "Retrieved  AOT at 0.47, 0.55,0.66   micron" 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Corrected_Optical_Depth_Land_wav2p1") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Retrieved  AOT at 2.13   micron" 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Optical_Depth_Ratio_Small_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Fraction of AOT contributed by fine dominated model" 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 1000 
          End If

          if (InputName .EQ. "Number_Pixels_Used_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_1_Land"
             long_name = 
     $ "Number of pixels used for land retrieval at 0.47 & 0.66 microns"
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 1 
             valid_range(2) = 400 
          End If

          if (InputName .EQ. "Mean_Reflectance_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Land"
             long_name = 
     $ "Mean reflectance of pixels used for land retrieval at 0.47,0.55,0.65,0.86,1.24,1.63,2.11 microns "
             units_str = "None" 
             scale_factor = 0.0001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 10000 
          End If

          if (InputName .EQ. "STD_Reflectance_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Land"
             long_name = 
     $ "Standard deviation of reflectance of pixels used for land retrieval at 0.47,0.55,0.65,0.86,1.24,1.63,2.11 microns "
             units_str = "None" 
             scale_factor = 0.0001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 20000 
          End If

          if (InputName .EQ. "Mass_Concentration_Land") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Estimated Column Mass(per area) using assumed mass extinction efficiency " 
             units_str = "1.0e-6g/cm^2" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -999. 
             valid_range(1) = 0. 
             valid_range(2) = 1000. 
          End If

          if (InputName .EQ. "Aerosol_Cloud_Fraction_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Cloud fraction from Land aerosol cloud mask from retrieved and " //
     $ "overcast pixels not including cirrus mask"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 1000 
          End If

          if (InputName .EQ. "Quality_Assurance_Land") then
             NumberType = DFNT_INT8
             DimList = 
     $ "QA_Byte_Land,Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Runtime QA flags" 
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             description(1) = "see MODIS atmosphere QA plan for details"
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = 0 
             valid_range(1) = 0 
             valid_range(2) = 255 
             DescCount = 1 
          End If

          if (InputName .EQ. "Solution_Index_Ocean_Small") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
             long_name = 
     $ "index identifying fine mode from Look Up Table for 'best' solution" 
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 1 
             valid_range(2) = 4 
          End If

          if (InputName .EQ. "Solution_Index_Ocean_Large") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
             long_name = 
     $ "index identifying coarse mode from Look Up Table for 'best' solution " 
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 5 
             valid_range(2) = 9 
          End If

          if (InputName .EQ. "Effective_Optical_Depth_Best_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Retrieved AOT for 'best' solution at 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um " 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Effective_Optical_Depth_Average_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Retrieved AOT for 'average' solution at 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Optical_Depth_Small_Best_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Retreived optical thickness for fine mode (best solution) for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Optical_Depth_Small_Average_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Retreived optical thickness for fine mode (Average solution) for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um" 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Optical_Depth_Large_Best_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Retrieved AOT of large mode for 'best' solution at 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um "
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Optical_Depth_Large_Average_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Retrieved AOT of large mode for 'average' solution at 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Mass_Concentration_Ocean") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
             long_name = 
     $ "Estimated Column Mass (per area) using assumed mass extinction coefficients for 'best' (1) and 'average' (2) solutions" 
             units_str = "1.0e-6g/cm^2" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -999. 
             valid_range(1) = 0. 
             valid_range(2) = 1000. 
          End If

          if (InputName .EQ. "Aerosol_Cloud_Fraction_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Cloud fraction from Land aerosol cloud mask from retrieved and " //
     $ " overcast pixels not including cirrus mask"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 1000 
          End If

          if (InputName .EQ. "Effective_Radius_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
             long_name = 
     $ "Effective_Radius at 0.55 micron for 'best' (1) and 'average' (2) solutions" 
             units_str = "micron" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "PSML003_Ocean") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
             long_name = 
     $ "Inferred column number concentration (number per area) of particles larger than 0.03 micron for"//
     $ "'best' (1) and 'average' (2) solutions"
             units_str = "Particles/cm^2"  
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -999. 
             valid_range(1) = 0. 
             valid_range(2) = 9.9999998e+10 
          End If

          if (InputName .EQ. "Asymmetry_Factor_Best_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Inferred Asymmetry_Factor for 'best' solution at 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um " 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 3000 
          End If

          if (InputName .EQ. "Asymmetry_Factor_Average_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Inferred Asymmetry_Factor for 'average' solution at 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 3000 
          End If

          if (InputName .EQ. "Backscattering_Ratio_Best_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Inferred Backscattering_Ratio for 'best' solution at 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 3000 
          End If

          if (InputName .EQ. "Backscattering_Ratio_Average_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $"Inferred Backscattering_Ratio for 'average' solution at 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 3000 
          End If

      if (InputName .EQ. "Angstrom_Exponent_1_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Calculated Angstrom Exponent for 0.55 vs 0.86 micron  for Average Solution " 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 10 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 10 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -1000 
             valid_range(2) = 5000 
          End If
 
          if (InputName .EQ. "Angstrom_Exponent_2_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Calculated Angstrom Exponent for 0.86 vs 2.13 micron for Average Solution "
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 10 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 10 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -1000 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Least_Squares_Error_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
             long_name = 
     $ "Residual of least squares fitting for inversion over land for best (1) and average (2) solutions" 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 1000 
          End If

          if (InputName .EQ. "Optical_Depth_Ratio_Small_Ocean_0.55micron") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Fraction of AOT (at 0.55 micron) contributed by fine mode for average solution"
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 10 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 10 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 1000 
          End If

          if (InputName .EQ. "Optical_Depth_by_models_ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Solution_Index"
             long_name = 
     $ "Retrieved AOT (at 0.55 micron) partioned by mode index (for Average solution)" 
             units_str = "None" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = -100 
             valid_range(2) = 5000 
          End If

          if (InputName .EQ. "Number_Pixels_Used_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Number of pixels used for ocean retrieval at 865 nm" 
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 1 
             valid_range(2) = 400 
          End If

          if (InputName .EQ. "Mean_Reflectance_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Mean reflectance of pixels used for ocean retrieval at 0.47,0.55,0.65,0.86,1.24,1.63,2.11 microns"
             units_str = "None" 
             scale_factor = 0.0001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 10000 
          End If

          if (InputName .EQ. "STD_Reflectance_Ocean") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
             long_name = 
     $ "Standard deviation of reflectance of pixels used for ocean retrieval at 0.47,0.55,0.65,0.86,1.24,1.63,2.11 microns"
             units_str = "None" 
             scale_factor = 0.0001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = -9999 
             valid_range(1) = 0 
             valid_range(2) = 20000 
          End If

          if (InputName .EQ. "Quality_Assurance_Ocean") then
             NumberType = DFNT_INT8
             DimList = 
     $ "QA_Byte_Ocean,Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Run time QA flags" 
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             description(1) = "(see MODIS atmosphere QA plan for details)"
             Geolocation_Pointer = "Internal geolocation arrays" 
             Fill_Value = 0 
             valid_range(1) = 0 
             valid_range(2) = 255 
             DescCount = 1 
          End If 
            if (InputName .EQ. "Topographic_Altitude_Land") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             valid_range(1) = 0 
             valid_range(2) = 1400
             Fill_Value = -9999 
             long_name = 
     $ "Averaged topographic altitude (in km) for Land" 
             units_str = "KM" 
             scale_factor = 0.01 
             add_offset = 0 
             Param_Type = "Output" 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2060 
             Cell_Along_Swath_Sampling(3) = 3 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 3 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If
       End If

       return
       end



       SUBROUTINE GET_SPEC_DATA(PROGRAMNAME, FIELDNAMES, DIMNAMES,
     $        DIMNUMBERS, DIMCOUNT, DATANAMES, DATACOUNT,
     $        DATAVALUES, NUMOFDATA, LIMIT, FIELDCOUNT)

       implicit none
       save

c---------------------------------------------------------------------
C!F77
c
c!DESCRIPTION:
c  Subroutine which extracts the correct file spec information
c  for the given product (MOD35, MOD07, or MOD06)
c
c!INPUT PARAMETERS:
c ProgramName   Ouput hdf file name
c Limit         Maximum size of Arrays
c
c!OUTPUT PARAMETERS:
c FieldNames    HDF Field Names for SDS
c FieldCount    Number of Dimension Values
c DimNames      Dimension Names for SDS
c DimNumbers    Dimension Values
c DimCount      Number of Dimension Values
c DataNames     Data Names for SDS
c DataNumbers   Data Values
c NumOfData     Number of Data Values
c
c!REVISION HISTORY:
c   This FORTRAN source file is created by spec2fort.c
c   spec2fort.c was designed by Walter.Wolf@ssec.wisc.edu
c
c!TEAM-UNIQUE HEADER:
c
c!REFERENCES AND CREDITS:
c
c!END
c---------------------------------------------------------------------

       character*(*) ProgramName
       character*(*) FieldNames(*), DimNames(*), DataNames(*)
       integer Limit, FieldCount, NumOfData
       integer DimNumbers(*), DataCount(*), DimCount

       real DataValues(Limit, Limit)
       integer unlimited

       NumOfData = 0
       DimCount = 0
       FieldCount = 0

       unlimited = 0

       if (ProgramName .EQ. "MOD04_L2") then
          DimNames(1) = "Cell_Along_Swath"
          DimNumbers(1) = 686
          DimNames(2) = "Cell_Across_Swath"
          DimNumbers(2) = 500
          DimNames(3) = "Solution_3_Land"
          DimNumbers(3) = 3
          DimNames(4) = "Solution_1_Land"
          DimNumbers(4) = 2
          DimNames(5) = "Solution_2_Land"
          DimNumbers(5) = 3
          DimNames(6) = "Solution_4_Land"
          DimNumbers(6) = 4
          DimNames(7) = "MODIS_Band_Land"
          DimNumbers(7) = 7
          DimNames(8) = "QA_Byte_Land"
          DimNumbers(8) = 5
          DimNames(9) = "Num_By_Products"
          DimNumbers(9) = 7
          DimNames(10) = "Solution_Ocean"
          DimNumbers(10) = 2
          DimNames(11) = "MODIS_Band_Ocean"
          DimNumbers(11) = 7
          DimNames(12) = "Solution_Index"
          DimNumbers(12) = 9
          DimNames(13) = "QA_Byte_Ocean"
          DimNumbers(13) = 5
          DimCount = 13 
          DataNames(1) = "MODIS_Band_Land" 
          DataValues(1,1) = 47
          DataValues(1,2) = 55
          DataValues(1,3) = 65
          DataValues(1,4) = 86
          DataValues(1,5) = 1240
          DataValues(1,6) = 1630
          DataValues(1,7) = 2110
          DataCount(1) = 7
          DataNames(2) = "MODIS_Band_Ocean" 
          DataValues(2,1) = 47
          DataValues(2,2) = 55
          DataValues(2,3) = 65
          DataValues(2,4) = 86
          DataValues(2,5) = 1240
          DataValues(2,6) = 1630
          DataValues(2,7) = 2110
          DataCount(2) = 7
          DataNames(3) = "Solution_1_Land" 
          DataValues(3,1) = 47
          DataValues(3,2) = 65
          DataCount(3) = 2
          DataNames(4) = "Solution_2_Land" 
          DataValues(4,1) = 47
          DataValues(4,2) = 55
          DataValues(4,3) = 65
          DataCount(4) = 3
          DataNames(5) = "Solution_3_Land" 
          DataValues(5,1) = 47
          DataValues(5,2) = 55
          DataValues(5,3) = 65
          DataCount(5) = 3
          DataNames(6) = "Solution_Ocean" 
          DataValues(6,1) = 1
          DataValues(6,2) = 2
          DataCount(6) = 2
          DataNames(7) = "Solution_Index" 
          DataValues(7,1) = 1
          DataValues(7,2) = 2
          DataValues(7,3) = 3
          DataValues(7,4) = 4
          DataValues(7,5) = 5
          DataValues(7,6) = 6
          DataValues(7,7) = 7
          DataValues(7,8) = 8
          DataValues(7,9) = 9
          DataCount(7) = 9
          NumOfData = 7 
          FieldNames(1) = "MODIS_Band_Land"
          FieldNames(2) = "MODIS_Band_Ocean"
          FieldNames(3) = "Solution_1_Land"
          FieldNames(4) = "Solution_2_Land"
          FieldNames(5) = "Solution_3_Land"
          FieldNames(6) = "Solution_Ocean"
          FieldNames(7) = "Solution_Index"
          FieldNames(8) = "Longitude"
          FieldNames(9) = "Latitude"
          FieldNames(10) = "Scan_Start_Time"
          FieldNames(11) = "Solar_Zenith"
          FieldNames(12) = "Solar_Azimuth"
          FieldNames(13) = "Sensor_Zenith"
          FieldNames(14) = "Sensor_Azimuth"
          FieldNames(15) = "Scattering_Angle"
          FieldNames(16) = "Glint_Angle"
          FieldNames(17) = "Land_Ocean_Quality_Flag"
          FieldNames(18) = "Land_sea_Flag"
          FieldNames(19) = "Wind_Speed_Ncep_Ocean"
          FieldNames(20) = "Optical_Depth_Land_And_Ocean"
          FieldNames(21) = "Image_Optical_Depth_Land_And_Ocean"
          FieldNames(22) = "Aerosol_Type_Land"
          FieldNames(23) = "Fitting_Error_Land"
          FieldNames(24) = "Surface_Reflectance_Land"
          FieldNames(25) = "Corrected_Optical_Depth_Land"
          FieldNames(26) = "Corrected_Optical_Depth_Land_wav2p1"
          FieldNames(27) = "Optical_Depth_Ratio_Small_Land"
          FieldNames(28) = "Number_Pixels_Used_Land"
          FieldNames(29) = "Mean_Reflectance_Land"
          FieldNames(30) = "STD_Reflectance_Land"
          FieldNames(31) = "Mass_Concentration_Land"
          FieldNames(32) = "Aerosol_Cloud_Fraction_Land"
          FieldNames(33) = "Quality_Assurance_Land"
          FieldNames(34) = "Solution_Index_Ocean_Small"
          FieldNames(35) = "Solution_Index_Ocean_Large"
          FieldNames(36) = "Effective_Optical_Depth_Best_Ocean"
          FieldNames(37) = "Effective_Optical_Depth_Average_Ocean"
          FieldNames(38) = "Optical_Depth_Small_Best_Ocean"
          FieldNames(39) = "Optical_Depth_Small_Average_Ocean"
          FieldNames(40) = "Optical_Depth_Large_Best_Ocean"
          FieldNames(41) = "Optical_Depth_Large_Average_Ocean"
          FieldNames(42) = "Mass_Concentration_Ocean"
          FieldNames(43) = "Aerosol_Cloud_Fraction_Ocean"
          FieldNames(44) = "Effective_Radius_Ocean"
          FieldNames(45) = "PSML003_Ocean"
          FieldNames(46) = "Asymmetry_Factor_Best_Ocean"
          FieldNames(47) = "Asymmetry_Factor_Average_Ocean"
          FieldNames(48) = "Backscattering_Ratio_Best_Ocean"
          FieldNames(49) = "Backscattering_Ratio_Average_Ocean"
          FieldNames(50) = "Angstrom_Exponent_1_Ocean"
          FieldNames(51) = "Angstrom_Exponent_2_Ocean"
          FieldNames(52) = "Least_Squares_Error_Ocean"
          FieldNames(53) = "Optical_Depth_Ratio_Small_Ocean_0.55micron"
          FieldNames(54) = "Optical_Depth_by_models_ocean"
          FieldNames(55) = "Number_Pixels_Used_Ocean"
          FieldNames(56) = "Mean_Reflectance_Ocean"
          FieldNames(57) = "STD_Reflectance_Ocean"
          FieldNames(58) = "Quality_Assurance_Ocean"
          FieldNames(59) = "Topographic_Altitude_Land"

          FieldCount = 59

       End If

       return
       end
