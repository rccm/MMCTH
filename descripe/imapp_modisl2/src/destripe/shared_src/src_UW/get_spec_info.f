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
       integer Number_of_Instrument_Scans
       integer Maximum_Number_of_1km_Frames
       character*48 ScaleFactor_AddOffset_Application
       character*46 Band_Number
       character*96 Pressure_Levels
       character*5 Satellite_Name
       character*5 Instrument_Name
       character*17 Grid_Type
       character*5 Grid_Resolution
       character*9 Start_of_Compositing_Period
       character*9 End_of_Compositing_Period
       integer Files_Used_In_Compositing
       character*9 File_Creation_Date
       character*18 Bands
       character*17 Lat_Resolution
       integer Granule_Start_TAI
       integer Granule_End_TAI

       unlimited = 0
       TitleCount = 0
       HistoryCount = 0
       if (ProgramName .EQ. "MOD06") then
          
          Number_of_Instrument_Scans = 203
          Maximum_Number_of_1km_Frames = 1354
          TitleCount = 1
          title = "MODIS Level 2 Cloud Properties"
          HistoryCount = 1
          history = "$Id: MOD06_L2.CDL.fs,v 1.3 2012/02/07 20:59:34 wind Exp $"
          
          if (InputName .EQ. "Band_Number") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Band_Number"
             long_name = 
     $ "MODIS Band Number" 
             units_str = "none" 
             Geolocation_Pointer = "Geolocation data not applicable" 
          End If

          if (InputName .EQ. "Statistics_1km") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistic_Parameter_1km"
             long_name = 
     $ "Granule Statistics for parameters at 1x1 resolution" 
             units_str = "see description attribute" 
             Fill_Value = -999.9 
             description(1) = "\n"
             description(2) = 
     $ "Statistics_1km:" 
             description(3) = 
     $ "  1. Successful Retrieval Rate (%)" 
             description(4) = 
     $ "  2. Land Cover Fraction (%)" 
             description(5) = 
     $ "  3. Water Cover Fraction (%)" 
             description(6) = 
     $ "  4. Snow Cover Fraction (%)" 
             description(7) = 
     $ "  5. Cloud Cover Fraction (%)" 
             description(8) = 
     $ "  6. Water Cloud Detected (%)" 
             description(9) = 
     $ "  7. Ice Cloud Detected (%)" 
             description(10) = 
     $ "  8. Mean of Water Cloud Optical Thickness" 
             description(11) = 
     $ "  9.Mean of Ice Cloud Optical Thickness" 
             description(12) = 
     $ "  10.Mean of Water Cloud Effective Particle Radius (microns)" 
             description(13) = 
     $ "  11.Mean of Ice Cloud Effective Diameter (microns)" 
             description(14) = 
     $ "  12.Mean Liquid Water Cloud Top Pressure (mb)" 
             description(15) = 
     $ "  13.Mean Ice Cloud Top Pressure (mb)" 
             description(16) = 
     $ "  14.Mean Undetermined Cloud Top Pressure (mb)" 
             description(17) = 
     $ "  15.Mean Liquid Water Cloud Top Temperature (K)" 
             description(18) = 
     $ "  16.Mean Ice Cloud Top Temperature (K)" 
             description(19) = 
     $ "  17.Mean Undetermined Cloud Top Temperature (K)" 
             DescCount = 19 
          End If

          if (InputName .EQ. "Scan_Start_Time") then
             NumberType = DFNT_FLOAT64
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "TAI time at start of scan replicated across the swath" 
             units_str = "seconds since 1993-1-1 00:00:00.0 0" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.1558e9 
             Fill_Value = -999.9 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Geolocation data not applicable" 
          End If

          if (InputName .EQ. "Latitude") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Geodetic Latitude" 
             units_str = "degrees_north" 
             valid_range(1) = -90.0 
             valid_range(2) = 90.0 
             Fill_Value = -999.9 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Longitude") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Geodetic Longitude" 
             units_str = "degrees_east" 
             valid_range(1) = -180.0 
             valid_range(2) = 180.0 
             Fill_Value = -999.9 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Solar_Zenith") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Solar Zenith Angle, Cell to Sun" 
             units_str = "degrees" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Solar_Zenith_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Solar Zenith Angle, Cell to Sun, Day Data Only" 
             units_str = "degrees" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Solar_Zenith_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Solar Zenith Angle, Cell to Sun, Night Data Only" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Solar_Azimuth") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Solar Azimuth Angle, Cell to Sun" 
             units_str = "degrees" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Solar_Azimuth_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Solar Azimuth Angle, Cell to Sun, Day Data Only" 
             units_str = "degrees" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Solar_Azimuth_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Solar Azimuth Angle, Cell to Sun, Night Data Only" 
             units_str = "degrees" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Sensor_Zenith") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Sensor Zenith Angle, Cell to Sensor" 
             units_str = "degrees" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Sensor_Zenith_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Sensor Zenith Angle, Cell to Sensor, Day Data Only" 
             units_str = "degrees" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Sensor_Zenith_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Sensor Zenith Angle, Cell to Sensor, Night Data Only" 
             units_str = "degrees" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Sensor_Azimuth") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Sensor Azimuth Angle, Cell to Sensor" 
             units_str = "degrees" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Sensor_Azimuth_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Sensor Azimuth Angle, Cell to Sensor, Day Data Only" 
             units_str = "degrees" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Sensor_Azimuth_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Sensor Azimuth Angle, Cell to Sensor, Night Data Only" 
             units_str = "degrees" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Brightness_Temperature") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km,Band_Number"
             long_name = 
     $ "Observed Brightness Temperature from Averaged Radiances in a 5x5 1-km " //
     $ "Pixel Region"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Surface_Temperature") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Surface Temperature from Ancillary Data" 
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Non MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Surface_Pressure") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Surface Pressure from Ancillary Data" 
             units_str = "hPa" 
             valid_range(1) = 8000 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Non MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Height_Method") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Index Indicating MODIS Bands Used for Cloud Top Pressure Retrieval" 
             units_str = "none" 
             valid_range(1) = 1 
             valid_range(2) = 6 
             Fill_Value = 127 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " The values in this SDS are set to mean the following:" 
             description(4) = 
     $ " 1 -- CO2-slicing retrieval, bands 36/35" 
             description(5) = 
     $ " 2 -- CO2-slicing retrieval, bands 35/34" 
             description(6) = 
     $ " 3 -- CO2-slicing retrieval, bands 35/33" 
             description(7) = 
     $ " 4 -- CO2-slicing retrieval, bands 34/33" 
             description(8) = 
     $ " 6 -- IR-window retrieval, band 31" 
             DescCount = 8 
          End If

          if (InputName .EQ. "Cloud_Top_Height") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Geopotential Height at Retrieved Cloud Top Pressure Level (rounded to nearest 50 m)" 
             units_str = "meters" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32767 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Height_Nadir") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Geopotential Height at Retrieved Cloud Top Pressure Level " //
     $ "for Sensor Zenith (View) Angles  <=32 Degrees (rounded to nearest 50 m)"
             units_str = "meters" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32767 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Height_Nadir_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Geopotential Height at Retrieved Cloud Top Pressure Level " //
     $ "for Sensor Zenith (View) Angles  <=32 Degrees, Day Data Only (rounded to nearest 50 m)"
             units_str = "meters" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32767 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Height_Nadir_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Geopotential Height at Retrieved Cloud Top Pressure Level " //
     $ "for Sensor Zenith (View) Angles  <=32 Degrees, Night Data Only (rounded to nearest 50 m)"
             units_str = "meters" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32767 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Pressure") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Top Pressure Level (rounded to nearest 5 mb)" 
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Pressure_Nadir") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Top Pressure Level for Sensor Zenith (View) " //
     $ "Angles <= 32 Degrees (rounded to nearest 5 mb)"
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Pressure_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Top Pressure Level, Night Data Only (rounded to nearest 5 mb)" 
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Pressure_Nadir_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Top Pressure Level for Sensor Zenith (View) " //
     $ "Angles <= 32 Degrees (rounded to nearest 5 mb), Night Data Only"
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Pressure_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Top Pressure Level, Day Only (rounded to nearest 5 mb)" 
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Pressure_Nadir_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Top Pressure Level for Sensor Zenith (View) " //
     $ "Angles <= 32 Degrees (rounded to nearest 5 mb), Day Data Only"
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Temperature") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level" 
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Temperature_Nadir") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Temperature from Ancillary Data at Retrieved Cloud Top " //
     $ "Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Temperature_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Temperature from Ancillary Data at Retrieved Cloud Top Pressure " //
     $ "Level, Night Only"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Temperature_Nadir_Night") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Temperature from Ancillary Data at Retrieved Cloud " //
     $ "Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees, Night Data Only"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Temperature_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Temperature from Ancillary Data at Retrieved Cloud Top Pressure " //
     $ "Level, Day Only"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Temperature_Nadir_Day") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Temperature from Ancillary Data at Retrieved Cloud " //
     $ "Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Tropopause_Height") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Tropopause Height from Ancillary Data" 
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Fraction") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Fraction_Nadir") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from " //
     $ "1-km Cloud Mask for Sensor Zenith (View) Angles <= 32 Degrees"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Fraction_Night") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud " //
     $ "Mask, Night Only"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Fraction_Nadir_Night") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from " //
     $ "1-km Cloud Mask for Sensor Zenith (View) Angles <= 32 Degrees, Night Data Only"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Fraction_Day") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud " //
     $ "Mask, Day Only"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Fraction_Nadir_Day") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from " //
     $ "1-km Cloud Mask for Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Effective_Emissivity") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Effective Emissivity from Cloud Top Pressure Retrieval" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Effective_Emissivity_Nadir") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Effective Emissivity from Cloud Top Pressure " //
     $ "Retrieval for Sensor Zenith (View) Angles <= 32 Degrees"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Effective_Emissivity_Night") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Effective Emissivity from Cloud Top Pressure " //
     $ "Retrieval, Night Only"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Effective_Emissivity_Nadir_Night") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Effective Emissivity from Cloud Top " //
     $ "Pressure Retrieval for Sensor Zenith (View) Angles <= 32 Degrees, Night Data Only"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Effective_Emissivity_Day") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Effective Emissivity from Cloud Top Pressure " //
     $ "Retrieval, Day Only"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Effective_Emissivity_Nadir_Day") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Effective Emissivity from Cloud Top " //
     $ "Pressure Retrieval for Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Pressure_Infrared") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Top Pressure from IR Window Retrieval" 
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Spectral_Cloud_Forcing") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km,Band_Forcing"
             long_name = 
     $ "Spectral Cloud Forcing (cloud minus clear radiance)" 
             units_str = "Watts/meter2/steradian/micron" 
             valid_range(1) = -2000 
             valid_range(2) = 2000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Top_Pressure_From_Ratios") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km,Band_Ratio"
             long_name = 
     $ "Cloud Top Pressure Levels from Ratios of Bands " //
     $ "36/35, 35/34, 35/33, 34/33 from the CO2-slicing Algorithm"
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " Cloud top pressure level solutions in the following array locations:" 
             description(4) = 
     $ " Band_Ratio:mod06 = 1:      MODIS bands 36/35" 
             description(5) = 
     $ " Band_Ratio:mod06 = 2:      MODIS bands 35/34" 
             description(6) = 
     $ " Band_Ratio:mod06 = 3:      MODIS bands 35/33" 
             description(7) = 
     $ " Band_Ratio:mod06 = 4:      MODIS bands 34/33" 
             description(8) = 
     $ " Band_Ratio:mod06 = 5:      Not used" 
             DescCount = 8 
          End If

          if (InputName .EQ. "Radiance_Variance") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Band 31 Radiance Standard Deviation " 
             units_str = "Watts/meter2/steradian/micron " 
             valid_range(1) = 0 
             valid_range(2) = 20 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Phase_Infrared") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Phase from 8.5 and 11 um Bands" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 6 
             Fill_Value = 127 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " The values in this SDS indicate the following cloud phase:" 
             description(4) = 
     $ " 0 -- cloud free" 
             description(5) = 
     $ " 1 -- water cloud" 
             description(6) = 
     $ " 2 -- ice cloud" 
             description(7) = 
     $ " 3 -- mixed phase cloud" 
             description(8) = 
     $ " 6 -- undetermined phase" 
             DescCount = 8 
          End If

          if (InputName .EQ. "Cloud_Phase_Infrared_Night") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Phase from 8.5 and 11 um Bands, Night Only" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 6 
             Fill_Value = 127 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " The values in this SDS indicate the following cloud phase:" 
             description(4) = 
     $ " 0 -- cloud free" 
             description(5) = 
     $ " 1 -- water cloud" 
             description(6) = 
     $ " 2 -- ice cloud" 
             description(7) = 
     $ " 3 -- mixed phase cloud" 
             description(8) = 
     $ " 6 -- undetermined phase" 
             DescCount = 8 
          End If

          if (InputName .EQ. "Cloud_Phase_Infrared_Day") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Cloud Phase from 8.5 and 11 um Bands, Day Only" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 6 
             Fill_Value = 127 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " The values in this SDS indicate the following cloud phase:" 
             description(4) = 
     $ " 0 -- cloud free" 
             description(5) = 
     $ " 1 -- water cloud" 
             description(6) = 
     $ " 2 -- ice cloud" 
             description(7) = 
     $ " 3 -- mixed phase cloud" 
             description(8) = 
     $ " 6 -- undetermined phase" 
             DescCount = 8 
          End If

          if (InputName .EQ. "Cloud_Phase_Infrared_1km") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Phase at 1-km resolution from 8.5-" //
     $ " 11 um BTDs and cloud emissivity ratios (12/11, 8.5/11, and 7.2/11 um)"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 6 
             Fill_Value = 127 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " The values in this SDS indicate the following cloud phase:" 
             description(4) = 
     $ " 0 -- cloud free" 
             description(5) = 
     $ " 1 -- water cloud" 
             description(6) = 
     $ " 2 -- ice cloud" 
             description(7) = 
     $ " 3 -- mixed phase cloud" 
             description(8) = 
     $ " 6 -- undetermined phase" 
             DescCount = 8 
          End If

          if (InputName .EQ. "os_top_flag_1km") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Upper Tropospheric/Lower Stratospheric (UTLS) Cloud Flag at 1-km " //
     $ " resolution - valid from -50 to +50 Degrees Latitude"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 2 
             Fill_Value = 127 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " The values in this SDS indicate the following:" 
             description(4) = 
     $ " 0 -- stratospheric cloud test not performed" 
             description(5) = 
     $ " 1 -- stratospheric cloud not indicated" 
             description(6) = 
     $ " 2 -- stratospheric cloud indicated (BTD35-33 > 0.5K)" 
             DescCount = 6 
          End If

          if (InputName .EQ. "cloud_top_pressure_1km") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Top Pressure at 1-km resolution from LEOCAT, " //
     $ "Cloud Top Pressure Level"
             units_str = "hPa" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -999 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "cloud_top_height_1km") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Top Height at 1-km resolution from LEOCAT, " //
     $ "Geopotential Height at Retrieved Cloud Top Pressure Level"
             units_str = "meters" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -999 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "cloud_top_temperature_1km") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Top Temperature at 1-km resolution from LEOCAT, " //
     $ "Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -999 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "cloud_emissivity_1km") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Emissivity at 1-km resolution from LEOCAT " //
     $ "Cloud Top Pressure Retrieval"
             units_str = "unitless" 
             valid_range(1) = 0 
             valid_range(2) = 100 
             Fill_Value = 127 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "cloud_top_method_1km") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Index Indicating the MODIS Band(s) Used to Produce the Cloud Top " //
     $ "Pressure Result"
             units_str = "unitless" 
             valid_range(1) = 0 
             valid_range(2) = 6 
             Fill_Value = 127 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " The values in this SDS are set to mean the following:" 
             description(4) = 
     $ " 1 -- CO2-slicing retrieval, bands 36/35" 
             description(5) = 
     $ " 2 -- CO2-slicing retrieval, bands 35/34" 
             description(6) = 
     $ " 3 -- CO2-slicing retrieval, bands 35/33" 
             description(7) = 
     $ " 4 -- CO2-slicing retrieval, bands 34/33" 
             description(8) = 
     $ " 6 -- IR-window retrieval, band 31" 
             DescCount = 8 
          End If

          if (InputName .EQ. "surface_temperature_1km") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Surface Temperature for Each 1-km MODIS Pixel " //
     $ "Interplated from Ancillary Data"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -999 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Non MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "cloud_emiss11_1km") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "11 micron Cloud Emissivity at 1-km resolution from LEOCAT " //
     $ "for All Clouds"
             units_str = "unitless" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -999 
             scale_factor = 0.0001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "cloud_emiss12_1km") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "12 micron Cloud Emissivity at 1-km resolution from LEOCAT " //
     $ "for All Clouds"
             units_str = "unitless" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -999 
             scale_factor = 0.0001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "cloud_emiss13_1km") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "13.3 micron Cloud Emissivity at 1-km resolution from LEOCAT " //
     $ "for All Clouds"
             units_str = "unitless" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -999 
             scale_factor = 0.0001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "cloud_emiss85_1km") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "8.5 micron Cloud Emissivity at 1-km resolution from LEOCAT " //
     $ "for All Clouds"
             units_str = "unitless" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -999 
             scale_factor = 0.0001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Effective_Radius") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Particle Effective Radius two-" //
     $ "channel retrieval using band 7 and either band 1, 2, or 5 (specified in Quality_Assurance_1km)"
             units_str = "micron" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -9999 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Effective_Radius_16") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Particle Effective Radius two-" //
     $ "channel retrieval using band 6 and either band 1, 2, or 5 (specified in Quality_Assurance_1km)"
             units_str = "micron" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -9999 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Effective_Radius_37") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Particle Effective Radius two-" //
     $ "channel retrieval using band 20 and either band 1, 2, or 5 (specified in Quality_Assurance_1km)"
             units_str = "micron" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -9999 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Optical_Thickness") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Optical Thickness two-" //
     $ "channel retrieval using band 7 and either band 1, 2, or 5 (specified in Quality_Assurance_1km)"
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -9999 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Effective_Radius_1621") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Particle Effective Radius two-" //
     $ "channel retrieval using band 7 and band 6"
             units_str = "micron" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -9999 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Optical_Thickness_1621") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Optical Thickness two-" //
     $ "channel retrieval using band 7 and band 6 "
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -9999 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Water_Path") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Column Water Path two-" //
     $ "band retrieval using band 7 and either band 1, 2, or 5 (specified in Quality_Assurance_1km)"
             units_str = "g/m^2" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -9999 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Water_Path_1621") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Column Water Path two-" //
     $ "band retrieval using band 7 and band 6"
             units_str = "g/m^2" 
             valid_range(1) = 0 
             valid_range(2) = 10000 
             Fill_Value = -9999 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Effective_Radius_Uncertainty") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Effective Particle Radius (from band 7) Relative Uncertainty (Percent)" 
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Effective_Radius_Uncertainty_16") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Effective Particle Radius (from band 6) Relative Uncertainty (Percent)" 
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Effective_Radius_Uncertainty_37") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Effective Particle Radius (from band 20) Relative Uncertainty (Percent)" 
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Optical_Thickness_Uncertainty") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Optical Thickness Relative Uncertainty (Percent)" 
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Water_Path_Uncertainty") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Water Path Relative Uncertainty (Percent)" 
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Effective_Radius_Uncertainty_1621") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Effective Particle Radius Relative Uncertainty (Percent) using band 7 and band 6" 
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Optical_Thickness_Uncertainty_1621") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Optical Thickness Relative Uncertainty (Percent) using band 7 and band 6" 
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cloud_Water_Path_Uncertainty_1621") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Water Path Relative Uncertainty (Percent) using band 7 and band 6" 
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Above_Cloud_Water_Vapor_094") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Above-cloud water vapor amount from 0.94um channel, ocean only, tau > 5." 
             units_str = "cm" 
             valid_range(1) = 0 
             valid_range(2) = 1500 
             Fill_Value = -9999 
             scale_factor = .01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "IRW_Low_Cloud_Temperature_From_COP") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Low Cloud Temperature from IR Window retrieval using cloud emissivity " //
     $ " based on cloud optical thickness"
             units_str = "K" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Cloud_Phase_Optical_Properties") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Phase Determination Used in Optical Thickness/Effective Radius Retrieval" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 4 
             Fill_Value = 0 
             scale_factor = 1. 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " The values in this SDS are set to mean the following:" 
             description(4) = 
     $ " 0 -- cloud mask undetermined" 
             description(5) = 
     $ " 1 -- clear sky" 
             description(6) = 
     $ " 2 -- liquid water cloud" 
             description(7) = 
     $ " 3 -- ice cloud" 
             description(8) = 
     $ " 4 -- undetermined phase cloud (but retrieval is attempted as " // 
     $ " liquid water)"
             DescCount = 8 
          End If

          if (InputName .EQ. "Cloud_Multi_Layer_Flag") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cloud Multi Layer Identification From MODIS Shortwave Observations" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 9 
             Fill_Value = 0 
             scale_factor = 1. 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " Flag for multi-layer multi-phase cloud situations. Values 2 through 9" 
             description(4) = 
     $ " indicate the success of various multi-layer cloud tests. Value of 0" 
             description(5) = 
     $ " indicates no retrieval, value of 1 indicates single layer cloud. The" 
             description(6) = 
     $ " other values are of increasing confidence level." 
             DescCount = 6 
          End If

          if (InputName .EQ. "Cirrus_Reflectance") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cirrus Reflectance" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 8000 
             Fill_Value = -9999 
             scale_factor = 0.0002 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Cirrus_Reflectance_Flag") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Cirrus Reflectance Flag" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 3 
             Fill_Value = 157 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
             description(1) = "0: bad data, 1: non-cirrus pixel, 2: cirrus pixel, 3: contrail pixel"
             DescCount = 1 
          End If

          if (InputName .EQ. "Cloud_Mask_5km") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cloud_Mask_5km_Num_Bytes,Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "First Byte of MODIS Cloud Mask Plus Additional Stats for L3 (2nd Byte)" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " Bit fields within each byte are numbered from the left:" 
             description(4) = 
     $ " 7, 6, 5, 4, 3, 2, 1, 0." 
             description(5) = 
     $ " The left-most bit (bit 7) is the most significant bit." 
             description(6) = 
     $ " The right-most bit (bit 0) is the least significant bit." 
             description(7) = 
     $ " " 
             description(8) = 
     $ " First Byte" 
             description(9) = 
     $ " " 
             description(10) = 
     $ " bit field       Description                             Key" 
             description(11) = 
     $ " ---------       -----------                             ---" 
             description(12) = 
     $ " " 
             description(13) = 
     $ " 0               Cloud Mask Flag                      0 = Not " // 
     $ " determined"
             description(14) = 
     $ "                                                      1 = Determined" 
             description(15) = 
     $ " " 
             description(16) = 
     $ " 2, 1            Unobstructed FOV Quality Flag        00 = Cloudy" 
             description(17) = 
     $ "                                                      01 = Uncertain" 
             description(18) = 
     $ "                                                      10 = Probably " // 
     $ " Clear"
             description(19) = 
     $ "                                                      11 = Confident " // 
     $ " Clear"
             description(20) = 
     $ "                 PROCESSING PATH" 
             description(21) = 
     $ "                 ---------------" 
             description(22) = 
     $ " 3               Day or Night Path                    0 = Night " // 
     $ " / 1 = Day"
             description(23) = 
     $ " 4               Sunglint Path                        0 = Yes " // 
     $ "   / 1 = No"
             description(24) = 
     $ " 5               Snow/Ice Background Path             0 = Yes " // 
     $ "   / 1 = No"
             description(25) = 
     $ " 7, 6            Land or Water Path                   00 = Water" 
             description(26) = 
     $ "                                                      01 = Coastal" 
             description(27) = 
     $ "                                                      10 = Desert" 
             description(28) = 
     $ "                                                      11 = Land" 
             description(29) = 
     $ " " 
             description(30) = 
     $ " Second Byte" 
             description(31) = 
     $ "--------------------------------------------------------------------------" 
             description(32) = 
     $ " " 
             description(33) = 
     $ " 1, 0            Sun-glint Under CTP Retrieval        00 = No CTP Ret." 
             description(34) = 
     $ "                                                      01 = No Sun-glint" 
             description(35) = 
     $ "                                                      10 = Sun-glint" 
             description(36) = 
     $ " " 
             description(37) = 
     $ " 3, 2            Snow/Ice Under CTP Retrieval         00 = No CTP Ret." 
             description(38) = 
     $ "                                                      01 = No Snow/Ice" 
             description(39) = 
     $ "                                                      10 = Snow/Ice" 
             description(40) = 
     $ " " 
             description(41) = 
     $ " 6, 5, 4         Surface Type Under CTP Retrieval    000 = No CTP Ret." 
             description(42) = 
     $ "                                                     001 = Water" 
             description(43) = 
     $ "                                                     010 = Coast" 
             description(44) = 
     $ "                                                     011 = Desert" 
             description(45) = 
     $ "                                                     100 = Land" 
             description(46) = 
     $ "                                                     101 = Other" 
             description(47) = 
     $ " 7               Spare" 
             description(48) = 
     $ " " 
             DescCount = 48 
          End If

          if (InputName .EQ. "Quality_Assurance_5km") then
             NumberType = DFNT_INT8
             DimList = 
     $ "QA_Parameter_5km,Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Quality Assurance at 5x5 Resolution" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             description(1) = "See MODIS atmosphere QA plan for details"
             Geolocation_Pointer = "Internal geolocation arrays" 
             DescCount = 1 
          End If

          if (InputName .EQ. "Cloud_Mask_1km") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cloud_Mask_1km_Num_Bytes,Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "MODIS Cloud Mask, L2 MOD06 QA Plan" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             description(1) = "See MODIS atmosphere QA plan for details"
             Geolocation_Pointer = "External MODIS geolocation product" 
             DescCount = 1 
          End If

          if (InputName .EQ. "Asymmetry_Parameter_Ice") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "RadTran_NWL,RadTran_NRE"
             long_name = 
     $ "Ice Asymmetry Parameter from the phase functions used to generate " //
     $ " the forward lookup tables"
             units_str = "none" 
             valid_range(1) = 0. 
             valid_range(2) = 1. 
             Fill_Value = 0. 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 7 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 15 
             Cell_Along_Swath_Sampling(3) = 1 
             description(1) = "Ice Asymmetry Parameter from the phase functions used to generate the forward lookup tables"
             Geolocation_Pointer = "Geolocation product not applicable" 
             DescCount = 1 
          End If

          if (InputName .EQ. "Single_Scatter_Albedo_Ice") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "RadTran_NWL,RadTran_NRE"
             long_name = 
     $ "Ice single scatter albedo from the phase functions used to generate " //
     $ " the forward lookup tables"
             units_str = "none" 
             valid_range(1) = 0. 
             valid_range(2) = 1. 
             Fill_Value = 0. 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 7 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 15 
             Cell_Along_Swath_Sampling(3) = 1 
             description(1) = "Ice single scatter albedo from the phase functions used to generate the forward lookup tables"
             Geolocation_Pointer = "Geolocation product not applicable" 
             DescCount = 1 
          End If

          if (InputName .EQ. "Cloud_Mask_SPI") then
             NumberType = DFNT_INT16
             DimList = 
     $ "SPI_nband,Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Dispersion in bands 1 (plane 1) and 2 (plane 2) from 250m reflectance " //
     $ " statistics of cloud mask"
             units_str = "percent" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
          End If

          if (InputName .EQ. "Quality_Assurance_1km") then
             NumberType = DFNT_INT8
             DimList = 
     $ "QA_Parameter_1km,Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Quality Assurance at 1x1 Resolution" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ "Quality Assurance 1km reports on Cloud Optical Properties" 
             description(4) = 
     $ "algorithm performance.  Refer to MOD_PR06OD User Documentation and the" 
             description(5) = 
     $ "MODIS atmosphere QA plan for complete descriptions and coding examples." 
             description(6) = 
     $ " " 
             description(7) = 
     $ " Bit fields within each byte are numbered from the left:" 
             description(8) = 
     $ " 7, 6, 5, 4, 3, 2, 1, 0." 
             description(9) = 
     $ " The left-most bit (bit 7) is the most significant bit." 
             description(10) = 
     $ " The right-most bit (bit 0) is the least significant bit." 
             description(11) = 
     $ " " 
             description(12) = 
     $ " " 
             description(13) = 
     $ " Byte 0" 
             description(14) = 
     $ " ------" 
             description(15) = 
     $ " " 
             description(16) = 
     $ " bit field       Description                             Key" 
             description(17) = 
     $ " ---------       -----------                             ---" 
             description(18) = 
     $ " " 
             description(19) = 
     $ " Byte 0 -----------------------------------------------------------------" 
             description(20) = 
     $ "  0              Optical Thickness General QA       0 = Not Useful" 
             description(21) = 
     $ "                                                    1 = Useful" 
             description(22) = 
     $ "  2,1            Optical Thickness Confidence QA    00 = No confidence" 
             description(23) = 
     $ "                                                    01 = Marginal" 
             description(24) = 
     $ "                                                    10 = Good" 
             description(25) = 
     $ "                                                    11 = Very Good" 
             description(26) = 
     $ "  4,3            Optical Thickness out-of-bounds    00 = OT < 100" 
             description(27) = 
     $ "                                                    01 = 100 < OT < 150" 
             description(28) = 
     $ "                                                    10 = OT > 150" 
             description(29) = 
     $ "                                                    11 = Albedo " // 
     $ " too high"
             description(30) = 
     $ "  5              Effective Radius General QA        0 = Not Useful" 
             description(31) = 
     $ "                                                    1 = Useful" 
             description(32) = 
     $ "  7,6            Effective Radius Confidence QA     00 = No confidence" 
             description(33) = 
     $ "                                                    01 = Marginal" 
             description(34) = 
     $ "                                                    10 = Good" 
             description(35) = 
     $ "                                                    11 = Very Good" 
             description(36) = 
     $ " Byte 1 -----------------------------------------------------------------" 
             description(37) = 
     $ "  0              Liquid Water Path General QA       0 = Not Useful" 
             description(38) = 
     $ "                                                    1 = Useful" 
             description(39) = 
     $ "  2,1            Liquid Water Path Confidence QA    00 = No confidence" 
             description(40) = 
     $ "                                                    01 = Marginal" 
             description(41) = 
     $ "                                                    10 = Good" 
             description(42) = 
     $ "                                                    11 = Very Good" 
             description(43) = 
     $ "  5,4,3          1621 Retrieval processing path     000 = No Cloud Mask" 
             description(44) = 
     $ "                                                    001 = No Cloud" 
             description(45) = 
     $ "                                                    010 = Water Cloud" 
             description(46) = 
     $ "                                                    011 = Ice Cloud" 
             description(47) = 
     $ "                                                    100 = Unknown Cloud" 
             description(48) = 
     $ "  6              1621 Retrieval Outcome             0 = Failed/No " // 
     $ " attempt"
             description(49) = 
     $ "                                                    1 = Successful" 
             description(50) = 
     $ " Byte 2 -----------------------------------------------------------------" 
             description(51) = 
     $ "  2,1,0          Primary retrieval processing path  000 = No Cloud Mask" 
             description(52) = 
     $ "                                                    001 = No Cloud" 
             description(53) = 
     $ "                                                    010 = Water Cloud" 
             description(54) = 
     $ "                                                    011 = Ice Cloud" 
             description(55) = 
     $ "                                                    100 = Unknown Cloud" 
             description(56) = 
     $ "  3              Retrieval Outcome                  0 = Failed/No " // 
     $ " attempt"
             description(57) = 
     $ "                                                    1 = Successful" 
             description(58) = 
     $ "  4              Rayleigh Correction                0 = No Correction" 
             description(59) = 
     $ "                                                    1 = Correction" 
             description(60) = 
     $ "  5              Water Vapor Correction             0 = No Correction" 
             description(61) = 
     $ "                                                    1 = Correction" 
             description(62) = 
     $ "  7,6            Band Used for Optical Thickness Retrieval" 
             description(63) = 
     $ "                                                    00 = No attempt" 
             description(64) = 
     $ "                                                    01 = .645 micron" 
             description(65) = 
     $ "                                                    10 = .858 micron" 
             description(66) = 
     $ "                                                    11 = 1.24 micron" 
             description(67) = 
     $ " Byte 3 -----------------------------------------------------------------" 
             description(68) = 
     $ "  0              Optical Thickness 1621 General QA  0 = Not Useful" 
             description(69) = 
     $ "                                                    1 = Useful" 
             description(70) = 
     $ "  2,1            Optical Thickness 1621 Condifence QA" 
             description(71) = 
     $ "                                                    00 = No confidence" 
             description(72) = 
     $ "                                                    01 = Marginal" 
             description(73) = 
     $ "                                                    10 = Good" 
             description(74) = 
     $ "                                                    11 = Very Good" 
             description(75) = 
     $ "  3              Effective Radius 1621 General QA   0 = Not Useful" 
             description(76) = 
     $ "                                                    1 = Useful" 
             description(77) = 
     $ "  5,4            Effective Radius 1621 Confidence QA" 
             description(78) = 
     $ "                                                    00 = No confidence" 
             description(79) = 
     $ "                                                    01 = Marginal" 
             description(80) = 
     $ "                                                    10 = Good" 
             description(81) = 
     $ "                                                    11 = Very Good" 
             description(82) = 
     $ "  6,7            Clear Sky Restoral Type QA" 
             description(83) = 
     $ "                                       00 = Not Restored" 
             description(84) = 
     $ "                                       01 = Restored Via Edge detection" 
             description(85) = 
     $ "                                       10 = Restored Via Spatial " // 
     $ " Variance"
             description(86) = 
     $ "                                       11 = Restored Via 250m Tests" 
             description(87) = 
     $ " Byte 4 -----------------------------------------------------------------" 
             description(88) = 
     $ "  0              Water Path 1621 General QA         0 = Not Useful" 
             description(89) = 
     $ "                                                    1 = Useful" 
             description(90) = 
     $ "  2,1            Water Path 1621 Confidence QA      00 = No confidence" 
             description(91) = 
     $ "                                                    01 = Marginal" 
             description(92) = 
     $ "                                                    10 = Good" 
             description(93) = 
     $ "                                                    11 = Very Good" 
             description(94) = 
     $ "  5,4,3          Multi Layer Cloud Flag      000 = Cloud Mask Undet" 
             description(95) = 
     $ "                                             001 = Decision tree stop" 
             description(96) = 
     $ "                                             010 = single layer: water" 
             description(97) = 
     $ "                                             011 = multi layer: water" 
             description(98) = 
     $ "                                             100 = single layer: ice" 
             description(99) = 
     $ "                                             101 = multi layer: ice" 
             description(100) = 
     $ "                                             110 = single layer: " // 
     $ " unknown"
             description(101) = 
     $ "                                             111 = multi layer: unknown" 
             description(102) = 
     $ "  7,6            spare" 
             description(103) = 
     $ " Byte 5 -----------------------------------------------------------------" 
             description(104) = 
     $ "  0              Phase difference multilayer test result" 
             description(105) = 
     $ "  1              Delta precipitable water multilayer test result" 
             description(106) = 
     $ "  2              Delta precipitable water at 900mb multilayer " // 
     $ " test result"
             description(107) = 
     $ "  3              Tau difference VIS-NIR multilayer test result" 
             description(108) = 
     $ "  4              Pavolonis-Heidinger multilayer test result" 
             description(109) = 
     $ "  7,6,5            spare" 
             DescCount = 109 
          End If
       End If

       if (ProgramName .EQ. "MOD07") then
          
          Number_of_Instrument_Scans = 203
          Maximum_Number_of_1km_Frames = 1354
          TitleCount = 1
          title = "MODIS Level 2 Atmospheric Profiles"
          ScaleFactor_AddOffset_Application = "Value=scale_factor*(stored integer - add_offset)"
          Band_Number = "24, 25, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36"
          Pressure_Levels = "5, 10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 620, 700, 780, 850, 920, 950, 1000 hPa"
          HistoryCount = 1
          history = "$Id: MOD07.V2.CDL,v 1.1 2005/12/14 16:44:04 kathys Exp $"
          
          if (InputName .EQ. "Band_Number") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Band_Number"
             long_name = 
     $ "MODIS Band Number" 
             units_str = "none" 
          End If

          if (InputName .EQ. "Pressure_Level") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Pressure_Level"
             long_name = 
     $ "Pressure Levels of Profiles" 
             units_str = "hPa" 
          End If

          if (InputName .EQ. "Pressure_Levels") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Pressure_Level"
             long_name = 
     $ "Pressure Levels" 
             units_str = "hPa" 
             scale_factor = 1.0 
             add_offset = 0.0 
          End If

          if (InputName .EQ. "Scan_Start_Time") then
             NumberType = DFNT_FLOAT64
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "TAI time at start of scan replicated across the swath" 
             units_str = "seconds since 1993-1-1 00:00:00.0 0" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.1558e9 
             Fill_Value = -999.9 
          End If

          if (InputName .EQ. "Latitude") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Geodetic Latitude" 
             units_str = "degrees_north" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = -90.0 
             valid_range(2) = 90.0 
             Fill_Value = -999.9 
          End If

          if (InputName .EQ. "Longitude") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Geodetic Longitude" 
             units_str = "degrees_east" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = -180.0 
             valid_range(2) = 180.0 
             Fill_Value = -999.9 
          End If

          if (InputName .EQ. "Solar_Zenith") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Solar Zenith Angle, Cell to Sun" 
             units_str = "degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Solar_Azimuth") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Solar Azimuth Angle, Cell to Sun" 
             units_str = "degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Sensor_Zenith") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Sensor Zenith Angle, Cell to Sensor" 
             units_str = "degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Sensor_Azimuth") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Sensor Azimuth Angle, Cell to Sensor" 
             units_str = "degrees" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Brightness_Temperature") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Band_Number"
             long_name = 
     $ "Brightness Temperature" 
             units_str = "K" 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Cloud_Mask") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "MODIS Cloud Mask, First Byte" 
             units_str = "none" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " Bit fields within each byte are numbered from the left:" 
             description(4) = 
     $ " 7, 6, 5, 4, 3, 2, 1, 0." 
             description(5) = 
     $ " The left-most bit (bit 7) is the most significant bit." 
             description(6) = 
     $ " The right-most bit (bit 0) is the least significant bit." 
             description(7) = 
     $ " " 
             description(8) = 
     $ " bit field       Description                             Key" 
             description(9) = 
     $ " ---------       -----------                             ---" 
             description(10) = 
     $ " " 
             description(11) = 
     $ " 0               Cloud Mask Flag                      0 = Not " // 
     $ " determined"
             description(12) = 
     $ "                                                      1 = Determined" 
             description(13) = 
     $ " " 
             description(14) = 
     $ " 2, 1            Unobstructed FOV Quality Flag        00 = Cloudy" 
             description(15) = 
     $ "                                                      01 = Uncertain" 
             description(16) = 
     $ "                                                      10 = Probably " // 
     $ " Clear"
             description(17) = 
     $ "                                                      11 = Confident " // 
     $ " Clear"
             description(18) = 
     $ "                 PROCESSING PATH" 
             description(19) = 
     $ "                 ---------------" 
             description(20) = 
     $ " 3               Day or Night Path                    0 = Night " // 
     $ " / 1 = Day"
             description(21) = 
     $ " 4               Sunglint Path                        0 = Yes " // 
     $ "   / 1 = No"
             description(22) = 
     $ " 5               Snow/Ice Background Path             0 = Yes " // 
     $ "   / 1 = No"
             description(23) = 
     $ " 7, 6            Land or Water Path                   00 = Water" 
             description(24) = 
     $ "                                                      01 = Coastal" 
             description(25) = 
     $ "                                                      10 = Desert" 
             description(26) = 
     $ "                                                      11 = Land" 
             DescCount = 26 
          End If

          if (InputName .EQ. "Skin_Temperature") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Skin Temperature" 
             units_str = "K" 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Surface_Pressure") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Surface Pressure" 
             units_str = "hPa" 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Non MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 8000 
             valid_range(2) = 11000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Surface_Elevation") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Surface Elevation" 
             units_str = "m" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Non MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = -400 
             valid_range(2) = 8840 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Processing_Flag") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Processing Flag" 
             units_str = "none" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 1 
             Fill_Value = 127 
          End If

          if (InputName .EQ. "Tropopause_Height") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Tropopause Height" 
             units_str = "hPa" 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 10 
             valid_range(2) = 11000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Guess_Temperature_Profile") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Pressure_Level"
             long_name = 
     $ "Guess Temperature Profile" 
             units_str = "K" 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Non MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Guess_Moisture_Profile") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Pressure_Level"
             long_name = 
     $ "Guess Mixing ratio Profile" 
             units_str = "g/kg" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Non MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Retrieved_Temperature_Profile") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Pressure_Level"
             long_name = 
     $ "Retrieved Temperature Profile" 
             units_str = "K" 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "MODIS Ouput" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Retrieved_Moisture_Profile") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Pressure_Level"
             long_name = 
     $ "Retrieved Dew Point Temperature Profile" 
             units_str = "K" 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "MODIS Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Retrieved_WV_Mixing_Ratio_Profile") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Pressure_Level"
             long_name = 
     $ "Retrieved Water Vapour Mixing Ratio Profile" 
             units_str = "g/kg" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Retrieved_Height_Profile") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath,Pressure_Level"
             long_name = 
     $ "Retrieved Geopotential Height Profile" 
             units_str = "m" 
             scale_factor = 1.0 
             add_offset = -32500.0 
             Param_Type = "MODIS Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = -32500 
             valid_range(2) = 32500 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Total_Ozone") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Total Ozone Burden" 
             units_str = "Dobson" 
             scale_factor = 0.1 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 5000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Total_Totals") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Total Totals" 
             units_str = "K" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 8000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Lifted_Index") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Lifted Index" 
             units_str = "K" 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = -2000 
             valid_range(2) = 4000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "K_Index") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "K_Index" 
             units_str = "K" 
             scale_factor = 0.01 
             add_offset = -15000.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 11500 
             valid_range(2) = 20000 
             Fill_Value = -32768 
          End If

          if (InputName .EQ. "Water_Vapor") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Total Column Precipitable Water Vapor - IR Retrieval" 
             units_str = "cm" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
          End If

          if (InputName .EQ. "Water_Vapor_Direct") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Total Column Precipitable Water Vapor - Direct IR Retrieval" 
             units_str = "cm" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
          End If

          if (InputName .EQ. "Water_Vapor_Low") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Precipitable Water Vapor Low - IR Retrieval" 
             units_str = "cm" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
          End If

          if (InputName .EQ. "Water_Vapor_High") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Precipitable Water Vapor High - IR Retrieval" 
             units_str = "cm" 
             scale_factor = 0.001 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             valid_range(1) = 0 
             valid_range(2) = 20000 
             Fill_Value = -9999 
          End If

          if (InputName .EQ. "Quality_Assurance") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Output_Parameter,Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Quality Assurance Parameters" 
             units_str = "none" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ "                Product and Optional run time QA flags" 
             description(3) = 
     $ " " 
             description(4) = 
     $ "                   Product run time QA flags" 
             description(5) = 
     $ " QA Flag Name     Number of     Bit Value   Description" 
             description(6) = 
     $ "                    Bits" 
             description(7) = 
     $ " -------------------------------------------------------------" 
             description(8) = 
     $ " Retrieved Temperature 1            0       not useful" 
             description(9) = 
     $ "   Profile QA                       1       useful" 
             description(10) = 
     $ " " 
             description(11) = 
     $ " Retrieved Temperature 2           0-3      4 confidence" 
             description(12) = 
     $ "   Profile                                    levels" 
             description(13) = 
     $ "   Confidence QA" 
             description(14) = 
     $ " " 
             description(15) = 
     $ " Spare                1" 
             description(16) = 
     $ " " 
             description(17) = 
     $ " Retrieved Moisture   1             0       not useful" 
             description(18) = 
     $ "   Profile                          1       useful" 
             description(19) = 
     $ "   QA" 
             description(20) = 
     $ " " 
             description(21) = 
     $ " Retrieved Moisture   2            0-3      4 confidence" 
             description(22) = 
     $ "   Profile                                    levels" 
             description(23) = 
     $ "   Confidence QA" 
             description(24) = 
     $ " " 
             description(25) = 
     $ " Spare                1" 
             description(26) = 
     $ " " 
             description(27) = 
     $ " Total Ozone          1             0       not useful" 
             description(28) = 
     $ "   Burden QA                        1       useful" 
             description(29) = 
     $ " " 
             description(30) = 
     $ " Total Ozone          2            0-3      4 confidence" 
             description(31) = 
     $ "   Burden                                     levels" 
             description(32) = 
     $ "   Confidence QA" 
             description(33) = 
     $ " " 
             description(34) = 
     $ " Spare                1" 
             description(35) = 
     $ " " 
             description(36) = 
     $ " Lifted Index         1             0       not useful" 
             description(37) = 
     $ "   Stability                        1       useful" 
             description(38) = 
     $ "   Index QA" 
             description(39) = 
     $ " " 
             description(40) = 
     $ " Lifted Index         2            0-3      4 confidence" 
             description(41) = 
     $ "   Stability                                  levels" 
             description(42) = 
     $ "   Confidence QA" 
             description(43) = 
     $ " " 
             description(44) = 
     $ " Spare                1" 
             description(45) = 
     $ " " 
             description(46) = 
     $ " K Index              1             0       not useful" 
             description(47) = 
     $ "   Stability                        1       useful" 
             description(48) = 
     $ "   Index QA" 
             description(49) = 
     $ " " 
             description(50) = 
     $ " K Index              2            0-3      4 confidence" 
             description(51) = 
     $ "   Stability                                  levels" 
             description(52) = 
     $ "   Confidence QA" 
             description(53) = 
     $ " " 
             description(54) = 
     $ " Spare                1" 
             description(55) = 
     $ " " 
             description(56) = 
     $ " Total Totals         1             0       not useful" 
             description(57) = 
     $ "   Stability                        1       useful" 
             description(58) = 
     $ "   Index QA" 
             description(59) = 
     $ " " 
             description(60) = 
     $ " Total Totals         2            0-3      4 confidence" 
             description(61) = 
     $ "   Stability                                  levels" 
             description(62) = 
     $ "   Confidence QA" 
             description(63) = 
     $ " " 
             description(64) = 
     $ " Spare                1" 
             description(65) = 
     $ " ---------------------- 3 bytes total -----------------------" 
             description(66) = 
     $ " " 
             description(67) = 
     $ "   Optional run time QA flags - processing path flags" 
             description(68) = 
     $ " QA Flag Name     Number of     Bit Value   Description" 
             description(69) = 
     $ "                    Bits" 
             description(70) = 
     $ " -------------------------------------------------------------" 
             description(71) = 
     $ " Number of            8            0-25" 
             description(72) = 
     $ "   Cloudy Pixels" 
             description(73) = 
     $ "   within 5x5 km" 
             description(74) = 
     $ "   box" 
             description(75) = 
     $ " " 
             description(76) = 
     $ " Number of            8            0-25" 
             description(77) = 
     $ "   Clear Pixels" 
             description(78) = 
     $ "   within 5x5 km" 
             description(79) = 
     $ "   box" 
             description(80) = 
     $ " " 
             description(81) = 
     $ " Number of            8            0-25" 
             description(82) = 
     $ "   Missing" 
             description(83) = 
     $ "   Pixels within" 
             description(84) = 
     $ "   5x5 km box" 
             description(85) = 
     $ " " 
             description(86) = 
     $ " Method               2             0       Statistical" 
             description(87) = 
     $ "  of Profiles                       1       Physical" 
             description(88) = 
     $ "  Retrieval                         2       Other" 
             description(89) = 
     $ "                                    3       No retrieval" 
             description(90) = 
     $ " " 
             description(91) = 
     $ " Method               2             0       RTE Perturbation" 
             description(92) = 
     $ "  of Ozone                          1       Upper and Lower" 
             description(93) = 
     $ "  Retrieval                                 Stratospheric Ozone Method" 
             description(94) = 
     $ "                                    2       Other" 
             description(95) = 
     $ "                                    3       No retrieval" 
             description(96) = 
     $ " " 
             description(97) = 
     $ " Spares               4" 
             description(98) = 
     $ " ---------------   4 bytes total -----------------------------" 
             description(99) = 
     $ " " 
             description(100) = 
     $ "   Optional run time QA flags - data resource flags" 
             description(101) = 
     $ " " 
             description(102) = 
     $ " QA Flag Name     Number of     Bit Value   Description" 
             description(103) = 
     $ "                    Bits" 
             description(104) = 
     $ " -------------------------------------------------------------" 
             description(105) = 
     $ " Guess Moisture       2             0       NCEP" 
             description(106) = 
     $ "   profile                          1       DAO" 
             description(107) = 
     $ "   source                           2       AIRS/AMSU" 
             description(108) = 
     $ "                                    3       Not used" 
             description(109) = 
     $ " " 
             description(110) = 
     $ " Guess Temperature    2             0       NCEP" 
             description(111) = 
     $ "   Profile                          1       DAO" 
             description(112) = 
     $ "   source                           2       AIRS/AMSU" 
             description(113) = 
     $ "                                    3       Not used" 
             description(114) = 
     $ " " 
             description(115) = 
     $ " Surface              2             0       NCEP" 
             description(116) = 
     $ "   Temperature                      1       DAO" 
             description(117) = 
     $ "   Over Land                        2       Other" 
             description(118) = 
     $ "                                    3       Not used" 
             description(119) = 
     $ " " 
             description(120) = 
     $ " Surface              2             0       Reynolds blended" 
             description(121) = 
     $ "   Temperature                      1       DAO" 
             description(122) = 
     $ "   Over Ocean                       2       Other" 
             description(123) = 
     $ "                                    3       Not used" 
             description(124) = 
     $ " " 
             description(125) = 
     $ " Surface              2             0       NCEP" 
             description(126) = 
     $ "   Pressure                         1       DAO" 
             description(127) = 
     $ "                                    2       Other" 
             description(128) = 
     $ "                                    3       Not used" 
             description(129) = 
     $ " " 
             description(130) = 
     $ " Ozone                2             0       TOMS" 
             description(131) = 
     $ "   Profile                          1       TOVS" 
             description(132) = 
     $ "   First                            2       DAO" 
             description(133) = 
     $ "   Guess                            3       Other" 
             description(134) = 
     $ " " 
             description(135) = 
     $ " Spares              12" 
             description(136) = 
     $ " ------------------ 2 bytes total ----------------------------" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             DescCount = 136 
          End If

          if (InputName .EQ. "Quality_Assurance_Infrared") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Water_Vapor_QA_Bytes,Cell_Across_Swath,Cell_Along_Swath"
             long_name = 
     $ "Run time QA flags" 
             units_str = "none" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ "  Water Vapor IR product and Optional run time QA flags" 
             description(4) = 
     $ " " 
             description(5) = 
     $ "                   Product run time QA flags" 
             description(6) = 
     $ " QA Flag Name     Number of     Bit Value   Description" 
             description(7) = 
     $ "                    Bits" 
             description(8) = 
     $ " -------------------------------------------------------------" 
             description(9) = 
     $ " IR Water Vapor       1             0       not useful" 
             description(10) = 
     $ "   QA                               1       useful" 
             description(11) = 
     $ " " 
             description(12) = 
     $ " IR Water Vapor       2            0-3      4 confidence" 
             description(13) = 
     $ "   Confidence QA                              levels" 
             description(14) = 
     $ " " 
             description(15) = 
     $ " Spares               5" 
             description(16) = 
     $ " ---------------------- 1 byte total -------------------------" 
             description(17) = 
     $ " " 
             description(18) = 
     $ "   Optional run time QA flags - processing path flags" 
             description(19) = 
     $ " QA Flag Name     Number of     Bit Value   Description" 
             description(20) = 
     $ "                    Bits" 
             description(21) = 
     $ " -------------------------------------------------------------" 
             description(22) = 
     $ " Number of            8            0-25" 
             description(23) = 
     $ "   Cloudy Pixels" 
             description(24) = 
     $ "   within 5x5 km" 
             description(25) = 
     $ "   box" 
             description(26) = 
     $ " " 
             description(27) = 
     $ " Number of            8            0-25" 
             description(28) = 
     $ "   Clear Pixels" 
             description(29) = 
     $ "   within 5x5 km" 
             description(30) = 
     $ "   box" 
             description(31) = 
     $ " " 
             description(32) = 
     $ " Number of            8            0-25" 
             description(33) = 
     $ "   Missing" 
             description(34) = 
     $ "   Pixels within" 
             description(35) = 
     $ "   5x5 km box" 
             description(36) = 
     $ " " 
             description(37) = 
     $ " IR Water Vapor       2             0      Split Window (11-" // 
     $ " 12) technique"
             description(38) = 
     $ "   Retrieval                        1      Integration of moisture " // 
     $ " profile"
             description(39) = 
     $ "   Method Used                      2      Other" 
             description(40) = 
     $ "                                    3      No Retrieval" 
             description(41) = 
     $ " " 
             description(42) = 
     $ " Spares               6" 
             description(43) = 
     $ " ---------------------- 4 bytes total -----------------------" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             DescCount = 43 
          End If
       End If

       if (ProgramName .EQ. "MOD35") then
          
          Number_of_Instrument_Scans = 203
          Maximum_Number_of_1km_Frames = 1354
          TitleCount = 1
          title = "MODIS Level 2 Cloud Mask"
          HistoryCount = 1
          history = "$Id: MOD35.V2.CDL,v 1.1.2.5 2002/09/16 18:24:50 raf Exp $"
          
          if (InputName .EQ. "Byte_Segment") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Byte_Segment"
          End If

          if (InputName .EQ. "Scan_Start_Time") then
             NumberType = DFNT_FLOAT64
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "TAI time at start of scan replicated across the swath" 
             units_str = "seconds since 1993-1-1 00:00:00.0 0" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.1558e9 
             Fill_Value = -999.9 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Latitude") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Geodetic Latitude" 
             units_str = "degrees_north" 
             valid_range(1) = -90.0 
             valid_range(2) = 90.0 
             Fill_Value = -999.99 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Longitude") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Geodetic Longitude" 
             units_str = "degrees_east" 
             valid_range(1) = -180.0 
             valid_range(2) = 180.0 
             Fill_Value = -999.99 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Solar_Zenith") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Solar Zenith Angle, Cell to Sun" 
             units_str = "degrees" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32767 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Solar_Azimuth") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Solar Azimuth Angle, Cell to Sun" 
             units_str = "degrees" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32767 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Sensor_Zenith") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Sensor Zenith Angle, Cell to Sensor" 
             units_str = "degrees" 
             valid_range(1) = 0 
             valid_range(2) = 18000 
             Fill_Value = -32767 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Sensor_Azimuth") then
             NumberType = DFNT_INT16
             DimList = 
     $ "Cell_Across_Swath_5km,Cell_Along_Swath_5km"
             long_name = 
     $ "Sensor Azimuth Angle, Cell to Sensor" 
             units_str = "degrees" 
             valid_range(1) = -18000 
             valid_range(2) = 18000 
             Fill_Value = -32767 
             scale_factor = 0.01 
             add_offset = 0.0 
             Param_Type = "MODIS Input" 
             Cell_Across_Swath_Sampling(1) = 3 
             Cell_Across_Swath_Sampling(2) = 1348 
             Cell_Across_Swath_Sampling(3) = 5 
             Cell_Along_Swath_Sampling(1) = 3 
             Cell_Along_Swath_Sampling(2) = 2028 
             Cell_Along_Swath_Sampling(3) = 5 
             Geolocation_Pointer = "Internal geolocation arrays" 
          End If

          if (InputName .EQ. "Ref_250m_stats") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Stats_250m_Dim,Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "MODIS 4x4 250-m reflectance statistics" 
             units_str = "reflectance" 
             valid_range(1) = 0.0 
             valid_range(2) = 2.0 
             Fill_Value = -999.0 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " Statistics are (3rd dimension):" 
             description(4) = 
     $ "          index 1:   band 1 mean" 
             description(5) = 
     $ "          index 2:   band 1 standard deviation" 
             description(6) = 
     $ "          index 3:   band 2 mean" 
             description(7) = 
     $ "          index 4:   band 2 standard deviation" 
             description(8) = 
     $ " " 
             DescCount = 8 
          End If

          if (InputName .EQ. "Cloud_Mask") then
             NumberType = DFNT_INT8
             DimList = 
     $ "Cell_Across_Swath_1km,Cell_Along_Swath_1km,Byte_Segment"
             long_name = 
     $ "MODIS Cloud Mask and Spectral Test Results" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
             description(1) = "\n"
             description(2) = 
     $ " " 
             description(3) = 
     $ " Bit fields within each byte are numbered from the left:" 
             description(4) = 
     $ " 7, 6, 5, 4, 3, 2, 1, 0." 
             description(5) = 
     $ " The left-most bit (bit 7) is the most significant bit." 
             description(6) = 
     $ " The right-most bit (bit 0) is the least significant bit." 
             description(7) = 
     $ " " 
             description(8) = 
     $ " bit field       Description                             Key" 
             description(9) = 
     $ " ---------       -----------                             ---" 
             description(10) = 
     $ " " 
             description(11) = 
     $ " 0               Cloud Mask Flag                      0 = Not " // 
     $ " determined"
             description(12) = 
     $ "                                                      1 = Determined" 
             description(13) = 
     $ " " 
             description(14) = 
     $ " 2, 1            Unobstructed FOV Quality Flag        00 = Cloudy" 
             description(15) = 
     $ "                                                      01 = Uncertain" 
             description(16) = 
     $ "                                                      10 = Probably " // 
     $ " Clear"
             description(17) = 
     $ "                                                      11 = Confident " // 
     $ " Clear"
             description(18) = 
     $ "                 PROCESSING PATH" 
             description(19) = 
     $ "                 ---------------" 
             description(20) = 
     $ " 3               Day or Night Path                    0 = Night " // 
     $ " / 1 = Day"
             description(21) = 
     $ " 4               Sunglint Path                        0 = Yes " // 
     $ "   / 1 = No"
             description(22) = 
     $ " 5               Snow/Ice Background Path             0 = Yes " // 
     $ "   / 1 = No"
             description(23) = 
     $ " 7, 6            Land or Water Path                   00 = Water" 
             description(24) = 
     $ "                                                      01 = Coastal" 
             description(25) = 
     $ "                                                      10 = Desert" 
             description(26) = 
     $ "                                                      11 = Land" 
             description(27) = 
     $ " ____ END BYTE 1 ______________" // 
     $ " ___________________________________________"
             description(28) = 
     $ " " 
             description(29) = 
     $ " bit field       Description                             Key" 
             description(30) = 
     $ " ---------       -----------                             ---" 
             description(31) = 
     $ " " 
             description(32) = 
     $ "                 1-KM FLAGS" 
             description(33) = 
     $ "                 ----------------------" 
             description(34) = 
     $ " 0               Non-" // 
     $ " cloud obstruction Flag              0 = Yes / 1 = No"
             description(35) = 
     $ " 1               Thin Cirrus Detected  (Solar)           0 = " // 
     $ " Yes / 1 = No"
             description(36) = 
     $ " 2               Snow cover from ancillary map           0 = " // 
     $ " Yes / 1 = No"
             description(37) = 
     $ " 3               Thin Cirrus Detected  (Infrared)        0 = " // 
     $ " Yes / 1 = No"
             description(38) = 
     $ " 4               Cloud Adjacency (cloudy, probably       0 = " // 
     $ " Yes / 1 = No"
             description(39) = 
     $ "                 cloudy plus 1-pixel adjacent)" 
             description(40) = 
     $ " 5               Cloud Flag -" // 
     $ "  IR Threshold               0 = Yes / 1 = No"
             description(41) = 
     $ " 6               High Cloud Flag -" // 
     $ "  CO2 Test              0 = Yes / 1 = No"
             description(42) = 
     $ " 7               High Cloud Flag -" // 
     $ "  6.7 micron Test       0 = Yes / 1 = No"
             description(43) = 
     $ " ____ END BYTE 2 ______________" // 
     $ " ___________________________________________"
             description(44) = 
     $ " " 
             description(45) = 
     $ " bit field       Description                             Key" 
             description(46) = 
     $ " ---------       -----------                             ---" 
             description(47) = 
     $ " " 
             description(48) = 
     $ " 0               High Cloud Flag -" // 
     $ "  1.38 micron Test      0 = Yes / 1 = No"
             description(49) = 
     $ " 1               High Cloud Flag - 3.9-" // 
     $ " 12 micron Test    0 = Yes / 1 = No"
             description(50) = 
     $ " 2               Cloud Flag -" // 
     $ "  IR Temperature             0 = Yes / 1 = No"
             description(51) = 
     $ "                              Difference" 
             description(52) = 
     $ " 3               Cloud Flag - 3.9-" // 
     $ " 11 micron Test         0 = Yes / 1 = No"
             description(53) = 
     $ " 4               Cloud Flag -" // 
     $ "  Visible Reflectance Test   0 = Yes / 1 = No"
             description(54) = 
     $ " 5               Cloud Flag -" // 
     $ "  Visible/NIR Reflectance    0 = Yes / 1 = No"
             description(55) = 
     $ "                              Ratio Test" 
             description(56) = 
     $ " 6               Cloud Flag -" // 
     $ "  NDVI Clear Sky Restoral    0 = Yes / 1 = No"
             description(57) = 
     $ "                              Test" 
             description(58) = 
     $ " 7               Cloud Flag -" // 
     $ "  Night Land and Polar       0 = Yes / 1 = No"
             description(59) = 
     $ "                              7.3-11 micron Test" 
             description(60) = 
     $ " ____ END BYTE 3 ______________" // 
     $ " ___________________________________________"
             description(61) = 
     $ " " 
             description(62) = 
     $ " bit field       Description                             Key" 
             description(63) = 
     $ " ---------       -----------                             ---" 
             description(64) = 
     $ " " 
             description(65) = 
     $ " 0               Cloud Flag - Ocean 8.6-" // 
     $ " 11 micron Test   0 = Yes / 1 = No"
             description(66) = 
     $ " 1               Cloud Flag -" // 
     $ "  Clear Sky Restoral Test    0 = Yes / 1 = No"
             description(67) = 
     $ "                 Spatial Consistency" 
             description(68) = 
     $ " 2               Cloud Flag -" // 
     $ "  Clear Sky Restoral Test    0 = Yes / 1 = No"
             description(69) = 
     $ "                 Polar Night, Land, Sun-glint" 
             description(70) = 
     $ " 3               Cloud Flag -" // 
     $ "  Surface Temperature Test   0 = Yes / 1 = No"
             description(71) = 
     $ " 4               Suspended Dust Flag                     0 = " // 
     $ " Yes / 1 = No"
             description(72) = 
     $ " 5               Cloud Flag - Night Ocean 8.6-" // 
     $ " 7.3 micron 0 = Yes / 1 = No"
             description(73) = 
     $ "                 Test" 
             description(74) = 
     $ " 6               Cloud Flag -" // 
     $ "  Night Ocean 11 micon       0 = Yes / 1 = No"
             description(75) = 
     $ "                 Spatial Variability Test" 
             description(76) = 
     $ " 7               Cloud Flag -" // 
     $ "  Night Ocean Low Emissivity 0 = Yes / 1 = No"
             description(77) = 
     $ "                 Low Cloud 3.9-11 micron Test" 
             description(78) = 
     $ " ____ END BYTE 4 ______________" // 
     $ " ___________________________________________"
             description(79) = 
     $ " " 
             description(80) = 
     $ " bit field       Description                             Key" 
             description(81) = 
     $ " ---------       -----------                             ---" 
             description(82) = 
     $ " " 
             description(83) = 
     $ "                 250-m Cloud Flag - Visible Tests" 
             description(84) = 
     $ "                 --------------------------------" 
             description(85) = 
     $ " 0               Element(1,1)                            0 = " // 
     $ " Yes / 1 = No"
             description(86) = 
     $ " 1               Element(1,2)                            0 = " // 
     $ " Yes / 1 = No"
             description(87) = 
     $ " 2               Element(1,3)                            0 = " // 
     $ " Yes / 1 = No"
             description(88) = 
     $ " 3               Element(1,4)                            0 = " // 
     $ " Yes / 1 = No"
             description(89) = 
     $ " 4               Element(2,1)                            0 = " // 
     $ " Yes / 1 = No"
             description(90) = 
     $ " 5               Element(2,2)                            0 = " // 
     $ " Yes / 1 = No"
             description(91) = 
     $ " 6               Element(2,3)                            0 = " // 
     $ " Yes / 1 = No"
             description(92) = 
     $ " 7               Element(2,4)                            0 = " // 
     $ " Yes / 1 = No"
             description(93) = 
     $ " ____ END BYTE 5 ______________" // 
     $ " ___________________________________________"
             description(94) = 
     $ " " 
             description(95) = 
     $ " bit field       Description                             Key" 
             description(96) = 
     $ " ----------      -----------                             ---" 
             description(97) = 
     $ " " 
             description(98) = 
     $ " 0               Element(3,1)                            0 = " // 
     $ " Yes / 1 = No"
             description(99) = 
     $ " 1               Element(3,2)                            0 = " // 
     $ " Yes / 1 = No"
             description(100) = 
     $ " 2               Element(3,3)                            0 = " // 
     $ " Yes / 1 = No"
             description(101) = 
     $ " 3               Element(3,4)                            0 = " // 
     $ " Yes / 1 = No"
             description(102) = 
     $ " 4               Element(4,1)                            0 = " // 
     $ " Yes / 1 = No"
             description(103) = 
     $ " 5               Element(4,2)                            0 = " // 
     $ " Yes / 1 = No"
             description(104) = 
     $ " 6               Element(4,3)                            0 = " // 
     $ " Yes / 1 = No"
             description(105) = 
     $ " 7               Element(4,4)                            0 = " // 
     $ " Yes / 1 = No"
             description(106) = 
     $ " ____ END BYTE 6 ______________" // 
     $ " ___________________________________________"
             DescCount = 106 
          End If

          if (InputName .EQ. "Quality_Assurance") then
             NumberType = DFNT_INT8
             DimList = 
     $ "QA_Dimension,Cell_Across_Swath_1km,Cell_Along_Swath_1km"
             long_name = 
     $ "Quality Assurance for Cloud Mask" 
             units_str = "none" 
             valid_range(1) = 0 
             valid_range(2) = 255 
             Fill_Value = 0 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "Output" 
             Cell_Across_Swath_Sampling(1) = 1 
             Cell_Across_Swath_Sampling(2) = 1354 
             Cell_Across_Swath_Sampling(3) = 1 
             Cell_Along_Swath_Sampling(1) = 1 
             Cell_Along_Swath_Sampling(2) = 2030 
             Cell_Along_Swath_Sampling(3) = 1 
             Geolocation_Pointer = "External MODIS geolocation product" 
             description(1) = "\n"
             description(2) = 
     $ "        Product Run Time QA Flags" 
             description(3) = 
     $ " QA Flag Name    Number of Bits  Bit Value   Description" 
             description(4) = 
     $ " -------------------------------------------------------------" 
             description(5) = 
     $ " Cloud Mask QA        1             0       not useful" 
             description(6) = 
     $ "   (1km)                            1       useful" 
             description(7) = 
     $ " " 
             description(8) = 
     $ " Cloud Mask           3            0-7      8 confidence levels" 
             description(9) = 
     $ "   Confidence QA" 
             description(10) = 
     $ "   (1km)" 
             description(11) = 
     $ " " 
             description(12) = 
     $ " Spares               4" 
             description(13) = 
     $ "  -------- ( 1 byte total ) --------------------------------" 
             description(14) = 
     $ " " 
             description(15) = 
     $ " Optional run time QA flags - Individual test application" 
             description(16) = 
     $ " " 
             description(17) = 
     $ " QA Flag Name    Number of Bits  Bit Value   Description" 
             description(18) = 
     $ " -------------------------------------------------------------" 
             description(19) = 
     $ " NCO test             1             0       Not Applied" 
             description(20) = 
     $ "                                    1       Applied" 
             description(21) = 
     $ " " 
             description(22) = 
     $ " Thin Cirrus          1             0       Not Applied" 
             description(23) = 
     $ "   test (Solar)                     1       Applied" 
             description(24) = 
     $ " " 
             description(25) = 
     $ " Ancillary Snow Map   1             NA      NA" 
             description(26) = 
     $ " " 
             description(27) = 
     $ " Thin Cirrus          1             0       Not Applied" 
             description(28) = 
     $ "   test(IR)                         1       Applied" 
             description(29) = 
     $ " " 
             description(30) = 
     $ " Cloud Adjancency     1             0       Not Applied" 
             description(31) = 
     $ "   test(IR)                         1       Applied" 
             description(32) = 
     $ " " 
             description(33) = 
     $ " IR Threshold         1             0       Not Applied" 
             description(34) = 
     $ "   test                             1       Applied" 
             description(35) = 
     $ " " 
             description(36) = 
     $ " High Cloud           1             0       Not Applied" 
             description(37) = 
     $ "   Test (CO2)                       1       Applied" 
             description(38) = 
     $ " " 
             description(39) = 
     $ " High Cloud           1             0       Not Applied" 
             description(40) = 
     $ "   Test (6.7)                       1       Applied" 
             description(41) = 
     $ " " 
             description(42) = 
     $ " High Cloud           1             0       Not Applied" 
             description(43) = 
     $ "   Test (1.38)                      1       Applied" 
             description(44) = 
     $ " " 
             description(45) = 
     $ " High Cloud           1             0       Not Applied" 
             description(46) = 
     $ "   Test (3.9-12um)                  1       Applied" 
             description(47) = 
     $ " " 
             description(48) = 
     $ " IR Temperature       1             0       Not Applied" 
             description(49) = 
     $ "   Difference Tests                 1       Applied" 
             description(50) = 
     $ " " 
             description(51) = 
     $ " 3.9-11um Test        1             0       Not Applied" 
             description(52) = 
     $ "                                    1       Applied" 
             description(53) = 
     $ " " 
             description(54) = 
     $ " .68 Reflectance      1             0       Not Applied" 
             description(55) = 
     $ "   Test                             1       Applied" 
             description(56) = 
     $ " " 
             description(57) = 
     $ " Visible Ratio        1             0       Not Applied" 
             description(58) = 
     $ "   Test                             1       Applied" 
             description(59) = 
     $ " " 
             description(60) = 
     $ " NDVI Clear Sky       1             0       Not Applied" 
             description(61) = 
     $ "   Restoral Test                    1       Applied" 
             description(62) = 
     $ " " 
             description(63) = 
     $ " 7.3-11 um Test       1             0       Not Applied" 
             description(64) = 
     $ "                                    1       Applied" 
             description(65) = 
     $ " " 
             description(66) = 
     $ " 8.6-11 um Test       1             0       Not Applied" 
             description(67) = 
     $ "                                    1       Applied" 
             description(68) = 
     $ " " 
             description(69) = 
     $ " Spatial Variability  1             0       Not Applied" 
             description(70) = 
     $ "   Clear Sky Restoral               1       Applied" 
             description(71) = 
     $ "   Test" 
             description(72) = 
     $ " " 
             description(73) = 
     $ " Clear Sky Restoral   1             0       Not Applied" 
             description(74) = 
     $ "   Tests                            1       Applied" 
             description(75) = 
     $ " " 
             description(76) = 
     $ " Night Water          1             0       Not Applied" 
             description(77) = 
     $ "   Spatial Variability              1       Applied" 
             description(78) = 
     $ " " 
             description(79) = 
     $ " Suspended Dust Test  1             0       Not Applied" 
             description(80) = 
     $ "                                    1       Applied" 
             description(81) = 
     $ " " 
             description(82) = 
     $ " Night Water          1             0       Not Applied" 
             description(83) = 
     $ "   8.6-7.3 um Test                  1       Applied" 
             description(84) = 
     $ " " 
             description(85) = 
     $ " Night Water 11 um    1             0       Not Applied" 
             description(86) = 
     $ "   Variability Test                 1       Applied" 
             description(87) = 
     $ " " 
             description(88) = 
     $ " Night Water Low      1             0       Not Applied" 
             description(89) = 
     $ "    Emissivity Low                  1       Applied" 
             description(90) = 
     $ "    Cloud Test" 
             description(91) = 
     $ " " 
             description(92) = 
     $ " 250 m Visible      1(16)           0       Not Applied" 
             description(93) = 
     $ "   Tests                            1       Applied" 
             description(94) = 
     $ "   (16 times as (line1,element1),(line1, element2),..." 
             description(95) = 
     $ "   (line4,element4))" 
             description(96) = 
     $ " " 
             description(97) = 
     $ "  -------- ( 5 bytes total ) --------------------------------" 
             description(98) = 
     $ " " 
             description(99) = 
     $ " Optional run time QA flags - data information flags" 
             description(100) = 
     $ " " 
             description(101) = 
     $ " QA Flag Name    Number of Bits  Bit Value   Description" 
             description(102) = 
     $ "  -------------------------------------------------------------" 
             description(103) = 
     $ " Number of            2             0       None" 
             description(104) = 
     $ "   bands used to                    1       1-7" 
             description(105) = 
     $ "   generate                         2       8-14" 
             description(106) = 
     $ "   cloud mask                       3       15-21" 
             description(107) = 
     $ " " 
             description(108) = 
     $ " Number of            2             0       None" 
             description(109) = 
     $ "   spectral                         1       1-3" 
             description(110) = 
     $ "   tests used to                    2       4-6" 
             description(111) = 
     $ "   generate                         3       7-9" 
             description(112) = 
     $ "   cloud mask" 
             description(113) = 
     $ " " 
             description(114) = 
     $ " Spares               4" 
             description(115) = 
     $ " " 
             description(116) = 
     $ "  -------- ( 1 byte total ) --------------------------------" 
             description(117) = 
     $ " " 
             description(118) = 
     $ " Optional run time QA flags - data resource flags" 
             description(119) = 
     $ " " 
             description(120) = 
     $ " QA Flag Name    Number of Bits  Bit Value   Description" 
             description(121) = 
     $ " -------------------------------------------------------------" 
             description(122) = 
     $ " Clear                2             0       MOD35" 
             description(123) = 
     $ "   Radiance                         1       Model forward calculation" 
             description(124) = 
     $ "   Origin                           2       Other" 
             description(125) = 
     $ "                                    3       Not used" 
             description(126) = 
     $ " " 
             description(127) = 
     $ " Surface              2             0       NCEP GDAS" 
             description(128) = 
     $ "   Temperature                      1       DAO" 
             description(129) = 
     $ "   Over Land                        2       MOD11" 
             description(130) = 
     $ "                                    3       Other" 
             description(131) = 
     $ " " 
             description(132) = 
     $ " Surface              2             0       Reynolds blended" 
             description(133) = 
     $ "   Temperature                      1       DAO" 
             description(134) = 
     $ "   Over Ocean                       2       MOD28" 
             description(135) = 
     $ "                                    3       Other" 
             description(136) = 
     $ " " 
             description(137) = 
     $ " Surface Winds        2             0       NCEP GDAS" 
             description(138) = 
     $ "                                    1       DAO" 
             description(139) = 
     $ "                                    2       Other" 
             description(140) = 
     $ "                                    3       Not used" 
             description(141) = 
     $ " " 
             description(142) = 
     $ " Ecosystem Map        2             0       Loveland NA 1km" 
             description(143) = 
     $ "                                    1       Olson Ecosystem" 
             description(144) = 
     $ "                                    2       MOD12" 
             description(145) = 
     $ "                                    3       Other" 
             description(146) = 
     $ " " 
             description(147) = 
     $ " Snow mask            2             0       MOD33" 
             description(148) = 
     $ "                                    1       SSMI Product" 
             description(149) = 
     $ "                                    2       Other" 
             description(150) = 
     $ "                                    3       Not used" 
             description(151) = 
     $ " " 
             description(152) = 
     $ " Ice cover            2             0       MOD42" 
             description(153) = 
     $ "                                    1       SSMI product" 
             description(154) = 
     $ "                                    2       Other" 
             description(155) = 
     $ "                                    3       Not used" 
             description(156) = 
     $ " " 
             description(157) = 
     $ " Land/Sea Mask        2             0       USGS 1 km 6 level" 
             description(158) = 
     $ "                                    1       USGS 1 km binary" 
             description(159) = 
     $ "                                    2       Other" 
             description(160) = 
     $ "                                    3       Not used" 
             description(161) = 
     $ " " 
             description(162) = 
     $ " Digital              1             0       EOS DEM" 
             description(163) = 
     $ "   Elevation                        1       Not used" 
             description(164) = 
     $ "   Model" 
             description(165) = 
     $ " " 
             description(166) = 
     $ " Precipitable         2             0       NCEP GDAS" 
             description(167) = 
     $ "   Water                            1       DAO" 
             description(168) = 
     $ "                                    2       MOD07" 
             description(169) = 
     $ "                                    3       Other" 
             description(170) = 
     $ " " 
             description(171) = 
     $ " Spare                5" 
             description(172) = 
     $ " " 
             description(173) = 
     $ "  -------- ( 3 bytes total ) --------------------------------" 
             description(174) = 
     $ " " 
             DescCount = 174 
          End If
       End If

       if (ProgramName .EQ. "MOD_PRCSRD") then
          
          Satellite_Name = "TERRA"
          Instrument_Name = "MODIS"
          TitleCount = 1
          title = "MODIS Daily Clear Sky Grids"
          Grid_Type = "Global Equal Area"
          Grid_Resolution = "25 km"
          Start_of_Compositing_Period = "unlimited"
          End_of_Compositing_Period = "unlimited"
          Files_Used_In_Compositing = unlimited
          File_Creation_Date = "unlimited"
          HistoryCount = 1
          history = "$Id: MOD_PRCSRD.CDL,v 1.1 2005/12/14 16:44:04 kathys Exp $"
          
          if (InputName .EQ. "Band_Number") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Band_Number"
             long_name = 
     $ "MODIS Band Number" 
             units_str = "none" 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_01") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_02") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_03") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_04") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_05") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_06") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_07") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_17") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_18") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_19") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_26") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_20") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_21") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_22") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_23") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_24") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_25") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_27") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_28") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_29") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_30") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_31") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_32") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_33") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_34") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_35") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day_Band_36") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_26") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_20") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_21") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_22") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_23") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_24") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_25") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_27") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_28") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_29") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_30") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_31") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_32") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_33") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_34") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_35") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night_Band_36") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Land_Water_Composite_Day") then
             NumberType = DFNT_INT32
             DimList = 
     $ "LW_Statistics,Grid_size"
             long_name = 
     $ "Land/Water Statistics for daily grid boxes. 1 - Total Observations, " //
     $ " 2 - Number of Land pixels 3 - Number of Water pixels"
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 50000.0 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Land_Water_Composite_Night") then
             NumberType = DFNT_INT32
             DimList = 
     $ "LW_Statistics,Grid_size"
             long_name = 
     $ "Land/Water Statistics for nightly grid boxes. 1 - Total Observations, " //
     $ " 2 - Number of Land pixels 3 - Number of Water pixels"
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 50000.0 
             Fill_Value = 0.0 
          End If
       End If

       if (ProgramName .EQ. "MOD_PRCSR8") then
          
          Satellite_Name = "TERRA"
          Instrument_Name = "MODIS"
          TitleCount = 1
          title = "MODIS Eight Day Clear Sky Grids"
          Grid_Type = "Global Equal Area"
          Grid_Resolution = "25 km"
          Start_of_Compositing_Period = "unlimited"
          End_of_Compositing_Period = "unlimited"
          Files_Used_In_Compositing = unlimited
          File_Creation_Date = "unlimited"
          HistoryCount = 1
          history = "$Id: MOD_PRCSR8.CDL,v 1.2 2005/12/14 16:45:42 kathys Exp $"
          
          if (InputName .EQ. "Band_Number") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Band_Number"
             long_name = 
     $ "MODIS Band Number" 
             units_str = "none" 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_01") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_02") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_03") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_04") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_05") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_06") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_07") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_17") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_18") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_19") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_26") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_20") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_21") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_22") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_23") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_24") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_25") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_27") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_28") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_29") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_30") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_31") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_32") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_33") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_34") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_35") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Clear_Radiance_Band_36") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_26") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_20") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_21") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_22") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_23") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_24") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_25") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_27") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_28") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_29") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_30") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_31") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_32") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_33") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_34") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_35") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Clear_Radiance_Band_36") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Bands 1-7, 17-19, 26 Reflectance, Bands 20-25, 27-36 Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Day_Land_Water_Composite") then
             NumberType = DFNT_INT32
             DimList = 
     $ "LW_Statistics,Grid_size"
             long_name = 
     $ "Statistics for daily grid boxes. 1 - Total Observations, 2 - " //
     $ " Number of Land pixels 3 - Number of Water pixels"
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 50000.0 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Eight_Night_Land_Water_Composite") then
             NumberType = DFNT_INT32
             DimList = 
     $ "LW_Statistics,Grid_size"
             long_name = 
     $ "Statistics for nightly grid boxes. 1 - Total Observations, 2 " //
     $ " - Number of Land pixels 3 - Number of Water pixels"
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 50000.0 
             Fill_Value = 0.0 
          End If
       End If

       if (ProgramName .EQ. "MOD_PRCSRB") then
          
          Satellite_Name = "TERRA"
          Instrument_Name = "MODIS"
          TitleCount = 1
          title = "MODIS Daily Clear Sky Radiance Latitude Biases"
          Bands = "31, 33, 34, 35, 36"
          Lat_Resolution = "1 degree latitude"
          File_Creation_Date = "unlimited"
          HistoryCount = 1
          history = "$Id: MOD_PRCSRB.CDL,v 1.1 2005/10/20 21:52:02 kathys Exp $"
          
          if (InputName .EQ. "Band_Number") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Band_Number"
             long_name = 
     $ "MODIS Band Number" 
             units_str = "none" 
          End If

          if (InputName .EQ. "Ocean_Clear_Sky_Bias") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Latitude_bin,Band_Number"
             long_name = 
     $ "Ocean Latitudinally averaged eight day observed-calculated clear sky radiances" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = -50.0 
             valid_range(2) = 50.0 
             Fill_Value = -999.99 
          End If

          if (InputName .EQ. "Land_Day_Clear_Sky_Bias") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Latitude_bin,Band_Number"
             long_name = 
     $ "Land Day Latitudinally averaged eight day observed-calculated clear sky radiances" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = -50.0 
             valid_range(2) = 50.0 
             Fill_Value = -999.99 
          End If

          if (InputName .EQ. "Land_Night_Clear_Sky_Bias") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Latitude_bin,Band_Number"
             long_name = 
     $ "Land Night Latitudinal averaged eight day observed-calculated clear sky radiances" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = -50.0 
             valid_range(2) = 50.0 
             Fill_Value = -999.99 
          End If

          if (InputName .EQ. "Ocean_Day_Summations") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Summation_stats,Latitude_bin,Band_Number"
             long_name = 
     $ "Ocean Day eight day clear sky radiance summations" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Ocean_Night_Summations") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Summation_stats,Latitude_bin,Band_Number"
             long_name = 
     $ "Ocean Night eight day clear sky radiance summations" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Land_Day_Summations") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Summation_stats,Latitude_bin,Band_Number"
             long_name = 
     $ "Land Day eight day clear sky radiance summations" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Land_Night_Summations") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Summation_stats,Latitude_bin,Band_Number"
             long_name = 
     $ "Land Night eight day clear sky radiance summations" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = 0.0 
          End If

          if (InputName .EQ. "Ocean_Day_Diagnostics") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Diagnostic_stats,Latitude_bin,Band_Number"
             long_name = 
     $ "Ocean Day eight day clear sky radiance diagnostics (clear sum/total, square root of bias)" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = -999.99 
          End If

          if (InputName .EQ. "Ocean_Night_Diagnostics") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Diagnostic_stats,Latitude_bin,Band_Number"
             long_name = 
     $ "Ocean Night eight day clear sky radiance diagnostics (clear sum/total, " //
     $ " square root of bias)"
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = -999.99 
          End If

          if (InputName .EQ. "Land_Day_Diagnostics") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Diagnostic_stats,Latitude_bin,Band_Number"
             long_name = 
     $ "Land Day eight day clear sky radiance diagnostics (clear sum/total, square root of bias)" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = -999.99 
          End If

          if (InputName .EQ. "Land_Night_Diagnostics") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Diagnostic_stats,Latitude_bin,Band_Number"
             long_name = 
     $ "Land Night eight day clear sky radiance diagnostics (clear sum/total, square root of bias)" 
             units_str = "Watts/meter2/steradian/micron" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = -999.99 
          End If
       End If

       if (ProgramName .EQ. "MOD35_CLRRAD") then
          
          Granule_Start_TAI = unlimited
          Granule_End_TAI = unlimited
          TitleCount = 1
          title = "MODIS Level 2 Cloud Mask Clear Radiances"
          HistoryCount = 1
          history = "$Id: MOD35_CLRRAD.CDL,v 1.1 2005/12/14 16:44:04 kathys Exp $"
          
          if (InputName .EQ. "Band_Number") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Band_Number"
             long_name = 
     $ "MODIS Band Number" 
             units_str = "none" 
          End If

          if (InputName .EQ. "Grid_Index_Day") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Size_of_One,Number_Cells_Day"
             long_name = 
     $ "Index of daytime granule grid boxes of clear sky to global radiance grid map" 
             units_str = "none" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 1.0 
             valid_range(2) = 814880.0 
             Fill_Value = -999.0 
          End If

          if (InputName .EQ. "Grid_Index_Night") then
             NumberType = DFNT_INT32
             DimList = 
     $ "Size_of_One,Number_Cells_Night"
             long_name = 
     $ "Index of nighttime granule grid boxes of clear sky to global radiance grid map" 
             units_str = "none" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 1.0 
             valid_range(2) = 814880.0 
             Fill_Value = -999.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Day") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Band_Number,Number_Cells_Day"
             long_name = 
     $ "Statistics for daily global grid boxes" 
             units_str = "Statistic Dependent" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = -999.0 
          End If

          if (InputName .EQ. "Clear_Radiance_Night") then
             NumberType = DFNT_FLOAT32
             DimList = 
     $ "Statistics,Band_Number,Number_Cells_Night"
             long_name = 
     $ "Statistics for nightly global grid boxes" 
             units_str = "Statistic Dependent" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 3.40e38 
             Fill_Value = -999.0 
          End If

          if (InputName .EQ. "Land_Water_Day") then
             NumberType = DFNT_INT32
             DimList = 
     $ "LW_Statistics,Number_Cells_Day"
             long_name = 
     $ "Number of Land/Water observations per daytime grid box" 
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 50000.0 
             Fill_Value = -999.0 
          End If

          if (InputName .EQ. "Land_Water_Night") then
             NumberType = DFNT_INT32
             DimList = 
     $ "LW_Statistics,Number_Cells_Night"
             long_name = 
     $ "Number of Land/Water observations per nighttime grid box" 
             units_str = "None" 
             scale_factor = 1.0 
             add_offset = 0.0 
             Param_Type = "MODIS Output" 
             Geolocation_Pointer = "None" 
             valid_range(1) = 0.0 
             valid_range(2) = 50000.0 
             Fill_Value = -999.0 
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

       if (ProgramName .EQ. "MOD06") then
          DimNames(1) = "Cell_Across_Swath_5km"
          DimNumbers(1) = 270
          DimNames(2) = "Cell_Along_Swath_5km"
          DimNumbers(2) = 406
          DimNames(3) = "Cell_Across_Swath_1km"
          DimNumbers(3) = 1354
          DimNames(4) = "Cell_Along_Swath_1km"
          DimNumbers(4) = 2030
          DimNames(5) = "Cell_Across_Swath_hkm"
          DimNumbers(5) = 2708
          DimNames(6) = "Cell_Along_Swath_hkm"
          DimNumbers(6) = 4060
          DimNames(7) = "Band_Number"
          DimNumbers(7) = 7
          DimNames(8) = "Band_Ratio"
          DimNumbers(8) = 5
          DimNames(9) = "Band_Forcing"
          DimNumbers(9) = 5
          DimNames(10) = "Band_Difference"
          DimNumbers(10) = 2
          DimNames(11) = "QA_Parameter_5km"
          DimNumbers(11) = 10
          DimNames(12) = "QA_Parameter_1km"
          DimNumbers(12) = 9
          DimNames(13) = "Cloud_Mask_1km_Num_Bytes"
          DimNumbers(13) = 2
          DimNames(14) = "Cloud_Mask_5km_Num_Bytes"
          DimNumbers(14) = 2
          DimNames(15) = "RadTran_NWL"
          DimNumbers(15) = 7
          DimNames(16) = "RadTran_NRE"
          DimNumbers(16) = 15
          DimNames(17) = "SPI_nband"
          DimNumbers(17) = 2
          DimNames(18) = "Statistic_Parameter_1km"
          DimNumbers(18) = 17
          DimCount = 18 
          DataNames(1) = "Band_Number" 
          DataValues(1,1) = 29
          DataValues(1,2) = 31
          DataValues(1,3) = 32
          DataValues(1,4) = 33
          DataValues(1,5) = 34
          DataValues(1,6) = 35
          DataValues(1,7) = 36
          DataCount(1) = 7
          NumOfData = 1 
          FieldNames(1) = "Band_Number"
          FieldNames(2) = "Statistics_1km"
          FieldNames(3) = "Scan_Start_Time"
          FieldNames(4) = "Latitude"
          FieldNames(5) = "Longitude"
          FieldNames(6) = "Solar_Zenith"
          FieldNames(7) = "Solar_Zenith_Day"
          FieldNames(8) = "Solar_Zenith_Night"
          FieldNames(9) = "Solar_Azimuth"
          FieldNames(10) = "Solar_Azimuth_Day"
          FieldNames(11) = "Solar_Azimuth_Night"
          FieldNames(12) = "Sensor_Zenith"
          FieldNames(13) = "Sensor_Zenith_Day"
          FieldNames(14) = "Sensor_Zenith_Night"
          FieldNames(15) = "Sensor_Azimuth"
          FieldNames(16) = "Sensor_Azimuth_Day"
          FieldNames(17) = "Sensor_Azimuth_Night"
          FieldNames(18) = "Brightness_Temperature"
          FieldNames(19) = "Surface_Temperature"
          FieldNames(20) = "Surface_Pressure"
          FieldNames(21) = "Cloud_Height_Method"
          FieldNames(22) = "Cloud_Top_Height"
          FieldNames(23) = "Cloud_Top_Height_Nadir"
          FieldNames(24) = "Cloud_Top_Height_Nadir_Day"
          FieldNames(25) = "Cloud_Top_Height_Nadir_Night"
          FieldNames(26) = "Cloud_Top_Pressure"
          FieldNames(27) = "Cloud_Top_Pressure_Nadir"
          FieldNames(28) = "Cloud_Top_Pressure_Night"
          FieldNames(29) = "Cloud_Top_Pressure_Nadir_Night"
          FieldNames(30) = "Cloud_Top_Pressure_Day"
          FieldNames(31) = "Cloud_Top_Pressure_Nadir_Day"
          FieldNames(32) = "Cloud_Top_Temperature"
          FieldNames(33) = "Cloud_Top_Temperature_Nadir"
          FieldNames(34) = "Cloud_Top_Temperature_Night"
          FieldNames(35) = "Cloud_Top_Temperature_Nadir_Night"
          FieldNames(36) = "Cloud_Top_Temperature_Day"
          FieldNames(37) = "Cloud_Top_Temperature_Nadir_Day"
          FieldNames(38) = "Tropopause_Height"
          FieldNames(39) = "Cloud_Fraction"
          FieldNames(40) = "Cloud_Fraction_Nadir"
          FieldNames(41) = "Cloud_Fraction_Night"
          FieldNames(42) = "Cloud_Fraction_Nadir_Night"
          FieldNames(43) = "Cloud_Fraction_Day"
          FieldNames(44) = "Cloud_Fraction_Nadir_Day"
          FieldNames(45) = "Cloud_Effective_Emissivity"
          FieldNames(46) = "Cloud_Effective_Emissivity_Nadir"
          FieldNames(47) = "Cloud_Effective_Emissivity_Night"
          FieldNames(48) = "Cloud_Effective_Emissivity_Nadir_Night"
          FieldNames(49) = "Cloud_Effective_Emissivity_Day"
          FieldNames(50) = "Cloud_Effective_Emissivity_Nadir_Day"
          FieldNames(51) = "Cloud_Top_Pressure_Infrared"
          FieldNames(52) = "Spectral_Cloud_Forcing"
          FieldNames(53) = "Cloud_Top_Pressure_From_Ratios"
          FieldNames(54) = "Radiance_Variance"
          FieldNames(55) = "Cloud_Phase_Infrared"
          FieldNames(56) = "Cloud_Phase_Infrared_Night"
          FieldNames(57) = "Cloud_Phase_Infrared_Day"
          FieldNames(58) = "Cloud_Phase_Infrared_1km"
          FieldNames(59) = "os_top_flag_1km"
          FieldNames(60) = "cloud_top_pressure_1km"
          FieldNames(61) = "cloud_top_height_1km"
          FieldNames(62) = "cloud_top_temperature_1km"
          FieldNames(63) = "cloud_emissivity_1km"
          FieldNames(64) = "cloud_top_method_1km"
          FieldNames(65) = "surface_temperature_1km"
          FieldNames(66) = "cloud_emiss11_1km"
          FieldNames(67) = "cloud_emiss12_1km"
          FieldNames(68) = "cloud_emiss13_1km"
          FieldNames(69) = "cloud_emiss85_1km"
          FieldNames(70) = "Cloud_Effective_Radius"
          FieldNames(71) = "Cloud_Effective_Radius_16"
          FieldNames(72) = "Cloud_Effective_Radius_37"
          FieldNames(73) = "Cloud_Optical_Thickness"
          FieldNames(74) = "Cloud_Effective_Radius_1621"
          FieldNames(75) = "Cloud_Optical_Thickness_1621"
          FieldNames(76) = "Cloud_Water_Path"
          FieldNames(77) = "Cloud_Water_Path_1621"
          FieldNames(78) = "Cloud_Effective_Radius_Uncertainty"
          FieldNames(79) = "Cloud_Effective_Radius_Uncertainty_16"
          FieldNames(80) = "Cloud_Effective_Radius_Uncertainty_37"
          FieldNames(81) = "Cloud_Optical_Thickness_Uncertainty"
          FieldNames(82) = "Cloud_Water_Path_Uncertainty"
          FieldNames(83) = "Cloud_Effective_Radius_Uncertainty_1621"
          FieldNames(84) = "Cloud_Optical_Thickness_Uncertainty_1621"
          FieldNames(85) = "Cloud_Water_Path_Uncertainty_1621"
          FieldNames(86) = "Above_Cloud_Water_Vapor_094"
          FieldNames(87) = "IRW_Low_Cloud_Temperature_From_COP"
          FieldNames(88) = "Cloud_Phase_Optical_Properties"
          FieldNames(89) = "Cloud_Multi_Layer_Flag"
          FieldNames(90) = "Cirrus_Reflectance"
          FieldNames(91) = "Cirrus_Reflectance_Flag"
          FieldNames(92) = "Cloud_Mask_5km"
          FieldNames(93) = "Quality_Assurance_5km"
          FieldNames(94) = "Cloud_Mask_1km"
          FieldNames(95) = "Asymmetry_Parameter_Ice"
          FieldNames(96) = "Single_Scatter_Albedo_Ice"
          FieldNames(97) = "Cloud_Mask_SPI"
          FieldNames(98) = "Quality_Assurance_1km"

          FieldCount = 98

       End If

       if (ProgramName .EQ. "MOD07") then
          DimNames(1) = "Cell_Across_Swath"
          DimNumbers(1) = 270
          DimNames(2) = "Cell_Along_Swath"
          DimNumbers(2) = 406
          DimNames(3) = "Band_Number"
          DimNumbers(3) = 12
          DimNames(4) = "Pressure_Level"
          DimNumbers(4) = 20
          DimNames(5) = "Output_Parameter"
          DimNumbers(5) = 10
          DimNames(6) = "Water_Vapor_QA_Bytes"
          DimNumbers(6) = 5
          DimCount = 6 
          DataNames(1) = "Band_Number" 
          DataValues(1,1) = 24
          DataValues(1,2) = 25
          DataValues(1,3) = 27
          DataValues(1,4) = 28
          DataValues(1,5) = 29
          DataValues(1,6) = 30
          DataValues(1,7) = 31
          DataValues(1,8) = 32
          DataValues(1,9) = 33
          DataValues(1,10) = 34
          DataValues(1,11) = 35
          DataValues(1,12) = 36
          DataCount(1) = 12
          DataNames(2) = "Pressure_Level" 
          DataValues(2,1) = 05.
          DataValues(2,2) = 10.
          DataValues(2,3) = 20.
          DataValues(2,4) = 30.
          DataValues(2,5) = 50.
          DataValues(2,6) = 70.
          DataValues(2,7) = 100.
          DataValues(2,8) = 150.
          DataValues(2,9) = 200.
          DataValues(2,10) = 250.
          DataValues(2,11) = 300.
          DataValues(2,12) = 400.
          DataValues(2,13) = 500.
          DataValues(2,14) = 620.
          DataValues(2,15) = 700.
          DataValues(2,16) = 780.
          DataValues(2,17) = 850.
          DataValues(2,18) = 920.
          DataValues(2,19) = 950.
          DataValues(2,20) = 1000.
          DataCount(2) = 20
          NumOfData = 2 
          FieldNames(1) = "Band_Number"
          FieldNames(2) = "Pressure_Level"
          FieldNames(3) = "Pressure_Levels"
          FieldNames(4) = "Scan_Start_Time"
          FieldNames(5) = "Latitude"
          FieldNames(6) = "Longitude"
          FieldNames(7) = "Solar_Zenith"
          FieldNames(8) = "Solar_Azimuth"
          FieldNames(9) = "Sensor_Zenith"
          FieldNames(10) = "Sensor_Azimuth"
          FieldNames(11) = "Brightness_Temperature"
          FieldNames(12) = "Cloud_Mask"
          FieldNames(13) = "Skin_Temperature"
          FieldNames(14) = "Surface_Pressure"
          FieldNames(15) = "Surface_Elevation"
          FieldNames(16) = "Processing_Flag"
          FieldNames(17) = "Tropopause_Height"
          FieldNames(18) = "Guess_Temperature_Profile"
          FieldNames(19) = "Guess_Moisture_Profile"
          FieldNames(20) = "Retrieved_Temperature_Profile"
          FieldNames(21) = "Retrieved_Moisture_Profile"
          FieldNames(22) = "Retrieved_WV_Mixing_Ratio_Profile"
          FieldNames(23) = "Retrieved_Height_Profile"
          FieldNames(24) = "Total_Ozone"
          FieldNames(25) = "Total_Totals"
          FieldNames(26) = "Lifted_Index"
          FieldNames(27) = "K_Index"
          FieldNames(28) = "Water_Vapor"
          FieldNames(29) = "Water_Vapor_Direct"
          FieldNames(30) = "Water_Vapor_Low"
          FieldNames(31) = "Water_Vapor_High"
          FieldNames(32) = "Quality_Assurance"
          FieldNames(33) = "Quality_Assurance_Infrared"

          FieldCount = 33

       End If

       if (ProgramName .EQ. "MOD35") then
          DimNames(1) = "Cell_Across_Swath_1km"
          DimNumbers(1) = 1354
          DimNames(2) = "Cell_Across_Swath_5km"
          DimNumbers(2) = 270
          DimNames(3) = "Cell_Along_Swath_1km"
          DimNumbers(3) = 2030
          DimNames(4) = "Cell_Along_Swath_5km"
          DimNumbers(4) = 406
          DimNames(5) = "Byte_Segment"
          DimNumbers(5) = 6
          DimNames(6) = "QA_Dimension"
          DimNumbers(6) = 10
          DimNames(7) = "Stats_250m_Dim"
          DimNumbers(7) = 4
          DimCount = 7 
          DataNames(1) = "Byte_Segment" 
          DataValues(1,1) = 1
          DataValues(1,2) = 2
          DataValues(1,3) = 3
          DataValues(1,4) = 4
          DataValues(1,5) = 5
          DataValues(1,6) = 6
          DataCount(1) = 6
          NumOfData = 1 
          FieldNames(1) = "Byte_Segment"
          FieldNames(2) = "Scan_Start_Time"
          FieldNames(3) = "Latitude"
          FieldNames(4) = "Longitude"
          FieldNames(5) = "Solar_Zenith"
          FieldNames(6) = "Solar_Azimuth"
          FieldNames(7) = "Sensor_Zenith"
          FieldNames(8) = "Sensor_Azimuth"
          FieldNames(9) = "Ref_250m_stats"
          FieldNames(10) = "Cloud_Mask"
          FieldNames(11) = "Quality_Assurance"

          FieldCount = 11

       End If

       if (ProgramName .EQ. "MOD_PRCSRD") then
          DimNames(1) = "Grid_size"
          DimNumbers(1) = 814880
          DimNames(2) = "Band_Number"
          DimNumbers(2) = 27
          DimNames(3) = "Statistics"
          DimNumbers(3) = 9
          DimNames(4) = "LW_Statistics"
          DimNumbers(4) = 3
          DimCount = 4 
          DataNames(1) = "Band_Number" 
          DataValues(1,1) = 1
          DataValues(1,2) = 2
          DataValues(1,3) = 3
          DataValues(1,4) = 4
          DataValues(1,5) = 5
          DataValues(1,6) = 6
          DataValues(1,7) = 7
          DataValues(1,8) = 17
          DataValues(1,9) = 18
          DataValues(1,10) = 19
          DataValues(1,11) = 26
          DataValues(1,12) = 20
          DataValues(1,13) = 21
          DataValues(1,14) = 22
          DataValues(1,15) = 23
          DataValues(1,16) = 24
          DataValues(1,17) = 25
          DataValues(1,18) = 27
          DataValues(1,19) = 28
          DataValues(1,20) = 29
          DataValues(1,21) = 30
          DataValues(1,22) = 31
          DataValues(1,23) = 32
          DataValues(1,24) = 33
          DataValues(1,25) = 34
          DataValues(1,26) = 35
          DataValues(1,27) = 36
          DataCount(1) = 27
          NumOfData = 1 
          FieldNames(1) = "Band_Number"
          FieldNames(2) = "Clear_Radiance_Day_Band_01"
          FieldNames(3) = "Clear_Radiance_Day_Band_02"
          FieldNames(4) = "Clear_Radiance_Day_Band_03"
          FieldNames(5) = "Clear_Radiance_Day_Band_04"
          FieldNames(6) = "Clear_Radiance_Day_Band_05"
          FieldNames(7) = "Clear_Radiance_Day_Band_06"
          FieldNames(8) = "Clear_Radiance_Day_Band_07"
          FieldNames(9) = "Clear_Radiance_Day_Band_17"
          FieldNames(10) = "Clear_Radiance_Day_Band_18"
          FieldNames(11) = "Clear_Radiance_Day_Band_19"
          FieldNames(12) = "Clear_Radiance_Day_Band_26"
          FieldNames(13) = "Clear_Radiance_Day_Band_20"
          FieldNames(14) = "Clear_Radiance_Day_Band_21"
          FieldNames(15) = "Clear_Radiance_Day_Band_22"
          FieldNames(16) = "Clear_Radiance_Day_Band_23"
          FieldNames(17) = "Clear_Radiance_Day_Band_24"
          FieldNames(18) = "Clear_Radiance_Day_Band_25"
          FieldNames(19) = "Clear_Radiance_Day_Band_27"
          FieldNames(20) = "Clear_Radiance_Day_Band_28"
          FieldNames(21) = "Clear_Radiance_Day_Band_29"
          FieldNames(22) = "Clear_Radiance_Day_Band_30"
          FieldNames(23) = "Clear_Radiance_Day_Band_31"
          FieldNames(24) = "Clear_Radiance_Day_Band_32"
          FieldNames(25) = "Clear_Radiance_Day_Band_33"
          FieldNames(26) = "Clear_Radiance_Day_Band_34"
          FieldNames(27) = "Clear_Radiance_Day_Band_35"
          FieldNames(28) = "Clear_Radiance_Day_Band_36"
          FieldNames(29) = "Clear_Radiance_Night_Band_26"
          FieldNames(30) = "Clear_Radiance_Night_Band_20"
          FieldNames(31) = "Clear_Radiance_Night_Band_21"
          FieldNames(32) = "Clear_Radiance_Night_Band_22"
          FieldNames(33) = "Clear_Radiance_Night_Band_23"
          FieldNames(34) = "Clear_Radiance_Night_Band_24"
          FieldNames(35) = "Clear_Radiance_Night_Band_25"
          FieldNames(36) = "Clear_Radiance_Night_Band_27"
          FieldNames(37) = "Clear_Radiance_Night_Band_28"
          FieldNames(38) = "Clear_Radiance_Night_Band_29"
          FieldNames(39) = "Clear_Radiance_Night_Band_30"
          FieldNames(40) = "Clear_Radiance_Night_Band_31"
          FieldNames(41) = "Clear_Radiance_Night_Band_32"
          FieldNames(42) = "Clear_Radiance_Night_Band_33"
          FieldNames(43) = "Clear_Radiance_Night_Band_34"
          FieldNames(44) = "Clear_Radiance_Night_Band_35"
          FieldNames(45) = "Clear_Radiance_Night_Band_36"
          FieldNames(46) = "Land_Water_Composite_Day"
          FieldNames(47) = "Land_Water_Composite_Night"

          FieldCount = 47

       End If

       if (ProgramName .EQ. "MOD_PRCSR8") then
          DimNames(1) = "Grid_size"
          DimNumbers(1) = 814880
          DimNames(2) = "Band_Number"
          DimNumbers(2) = 27
          DimNames(3) = "Statistics"
          DimNumbers(3) = 9
          DimNames(4) = "LW_Statistics"
          DimNumbers(4) = 3
          DimCount = 4 
          DataNames(1) = "Band_Number" 
          DataValues(1,1) = 1
          DataValues(1,2) = 2
          DataValues(1,3) = 3
          DataValues(1,4) = 4
          DataValues(1,5) = 5
          DataValues(1,6) = 6
          DataValues(1,7) = 7
          DataValues(1,8) = 17
          DataValues(1,9) = 18
          DataValues(1,10) = 19
          DataValues(1,11) = 26
          DataValues(1,12) = 20
          DataValues(1,13) = 21
          DataValues(1,14) = 22
          DataValues(1,15) = 23
          DataValues(1,16) = 24
          DataValues(1,17) = 25
          DataValues(1,18) = 27
          DataValues(1,19) = 28
          DataValues(1,20) = 29
          DataValues(1,21) = 30
          DataValues(1,22) = 31
          DataValues(1,23) = 32
          DataValues(1,24) = 33
          DataValues(1,25) = 34
          DataValues(1,26) = 35
          DataValues(1,27) = 36
          DataCount(1) = 27
          NumOfData = 1 
          FieldNames(1) = "Band_Number"
          FieldNames(2) = "Eight_Day_Clear_Radiance_Band_01"
          FieldNames(3) = "Eight_Day_Clear_Radiance_Band_02"
          FieldNames(4) = "Eight_Day_Clear_Radiance_Band_03"
          FieldNames(5) = "Eight_Day_Clear_Radiance_Band_04"
          FieldNames(6) = "Eight_Day_Clear_Radiance_Band_05"
          FieldNames(7) = "Eight_Day_Clear_Radiance_Band_06"
          FieldNames(8) = "Eight_Day_Clear_Radiance_Band_07"
          FieldNames(9) = "Eight_Day_Clear_Radiance_Band_17"
          FieldNames(10) = "Eight_Day_Clear_Radiance_Band_18"
          FieldNames(11) = "Eight_Day_Clear_Radiance_Band_19"
          FieldNames(12) = "Eight_Day_Clear_Radiance_Band_26"
          FieldNames(13) = "Eight_Day_Clear_Radiance_Band_20"
          FieldNames(14) = "Eight_Day_Clear_Radiance_Band_21"
          FieldNames(15) = "Eight_Day_Clear_Radiance_Band_22"
          FieldNames(16) = "Eight_Day_Clear_Radiance_Band_23"
          FieldNames(17) = "Eight_Day_Clear_Radiance_Band_24"
          FieldNames(18) = "Eight_Day_Clear_Radiance_Band_25"
          FieldNames(19) = "Eight_Day_Clear_Radiance_Band_27"
          FieldNames(20) = "Eight_Day_Clear_Radiance_Band_28"
          FieldNames(21) = "Eight_Day_Clear_Radiance_Band_29"
          FieldNames(22) = "Eight_Day_Clear_Radiance_Band_30"
          FieldNames(23) = "Eight_Day_Clear_Radiance_Band_31"
          FieldNames(24) = "Eight_Day_Clear_Radiance_Band_32"
          FieldNames(25) = "Eight_Day_Clear_Radiance_Band_33"
          FieldNames(26) = "Eight_Day_Clear_Radiance_Band_34"
          FieldNames(27) = "Eight_Day_Clear_Radiance_Band_35"
          FieldNames(28) = "Eight_Day_Clear_Radiance_Band_36"
          FieldNames(29) = "Eight_Night_Clear_Radiance_Band_26"
          FieldNames(30) = "Eight_Night_Clear_Radiance_Band_20"
          FieldNames(31) = "Eight_Night_Clear_Radiance_Band_21"
          FieldNames(32) = "Eight_Night_Clear_Radiance_Band_22"
          FieldNames(33) = "Eight_Night_Clear_Radiance_Band_23"
          FieldNames(34) = "Eight_Night_Clear_Radiance_Band_24"
          FieldNames(35) = "Eight_Night_Clear_Radiance_Band_25"
          FieldNames(36) = "Eight_Night_Clear_Radiance_Band_27"
          FieldNames(37) = "Eight_Night_Clear_Radiance_Band_28"
          FieldNames(38) = "Eight_Night_Clear_Radiance_Band_29"
          FieldNames(39) = "Eight_Night_Clear_Radiance_Band_30"
          FieldNames(40) = "Eight_Night_Clear_Radiance_Band_31"
          FieldNames(41) = "Eight_Night_Clear_Radiance_Band_32"
          FieldNames(42) = "Eight_Night_Clear_Radiance_Band_33"
          FieldNames(43) = "Eight_Night_Clear_Radiance_Band_34"
          FieldNames(44) = "Eight_Night_Clear_Radiance_Band_35"
          FieldNames(45) = "Eight_Night_Clear_Radiance_Band_36"
          FieldNames(46) = "Eight_Day_Land_Water_Composite"
          FieldNames(47) = "Eight_Night_Land_Water_Composite"

          FieldCount = 47

       End If

       if (ProgramName .EQ. "MOD_PRCSRB") then
          DimNames(1) = "Latitude_bin"
          DimNumbers(1) = 181
          DimNames(2) = "Band_Number"
          DimNumbers(2) = 5
          DimNames(3) = "Summation_stats"
          DimNumbers(3) = 4
          DimNames(4) = "Diagnostic_stats"
          DimNumbers(4) = 2
          DimCount = 4 
          DataNames(1) = "Band_Number" 
          DataValues(1,1) = 31
          DataValues(1,2) = 33
          DataValues(1,3) = 34
          DataValues(1,4) = 35
          DataValues(1,5) = 36
          DataCount(1) = 5
          NumOfData = 1 
          FieldNames(1) = "Band_Number"
          FieldNames(2) = "Ocean_Clear_Sky_Bias"
          FieldNames(3) = "Land_Day_Clear_Sky_Bias"
          FieldNames(4) = "Land_Night_Clear_Sky_Bias"
          FieldNames(5) = "Ocean_Day_Summations"
          FieldNames(6) = "Ocean_Night_Summations"
          FieldNames(7) = "Land_Day_Summations"
          FieldNames(8) = "Land_Night_Summations"
          FieldNames(9) = "Ocean_Day_Diagnostics"
          FieldNames(10) = "Ocean_Night_Diagnostics"
          FieldNames(11) = "Land_Day_Diagnostics"
          FieldNames(12) = "Land_Night_Diagnostics"

          FieldCount = 12

       End If

       if (ProgramName .EQ. "MOD35_CLRRAD") then
          DimNames(1) = "Statistics"
          DimNumbers(1) = 9
          DimNames(2) = "LW_Statistics"
          DimNumbers(2) = 3
          DimNames(3) = "Band_Number"
          DimNumbers(3) = 27
          DimNames(4) = "Size_of_One"
          DimNumbers(4) = 1
          DimNames(5) = "Number_Cells_Day"
          DimNumbers(5) = unlimited
          DimNames(6) = "Number_Cells_Night"
          DimNumbers(6) = unlimited
          DimCount = 6 
          DataNames(1) = "Band_Number" 
          DataValues(1,1) = 1
          DataValues(1,2) = 2
          DataValues(1,3) = 3
          DataValues(1,4) = 4
          DataValues(1,5) = 5
          DataValues(1,6) = 6
          DataValues(1,7) = 7
          DataValues(1,8) = 17
          DataValues(1,9) = 18
          DataValues(1,10) = 19
          DataValues(1,11) = 26
          DataValues(1,12) = 20
          DataValues(1,13) = 21
          DataValues(1,14) = 22
          DataValues(1,15) = 23
          DataValues(1,16) = 24
          DataValues(1,17) = 25
          DataValues(1,18) = 27
          DataValues(1,19) = 28
          DataValues(1,20) = 29
          DataValues(1,21) = 30
          DataValues(1,22) = 31
          DataValues(1,23) = 32
          DataValues(1,24) = 33
          DataValues(1,25) = 34
          DataValues(1,26) = 35
          DataValues(1,27) = 36
          DataCount(1) = 27
          NumOfData = 1 
          FieldNames(1) = "Band_Number"
          FieldNames(2) = "Grid_Index_Day"
          FieldNames(3) = "Grid_Index_Night"
          FieldNames(4) = "Clear_Radiance_Day"
          FieldNames(5) = "Clear_Radiance_Night"
          FieldNames(6) = "Land_Water_Day"
          FieldNames(7) = "Land_Water_Night"

          FieldCount = 7

       End If

       return
       end
