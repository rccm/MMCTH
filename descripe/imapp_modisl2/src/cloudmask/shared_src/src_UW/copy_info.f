      INTEGER FUNCTION COPY_INFO(GEO_FILE, OUT_FILE, 
     &  START, STRIDE, EDGE, ERROR_TEXT)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C    This routine copies several MODIS L1B geolocation parameters
C    to an existing UW MODIS L2 product file. The parameters are
C
C    Latitude,
C    Longitude,
C    Solar zenith angle,
C    Solar azimuth angle,
C    Sensor zenith angle,
C    Sensor azimuth angle,
C    Scan start time.
C
C!INPUT PARAMETERS:
C    GEO_FILE    Name of MODIS L1B geolocation file
C    OUT_FILE    Name of existing UW MODIS L2 product file
C    START       2-element array specifying the starting location for
C                each dimension of the 2D arrays (ZERO BASED)
C    STRIDE      2-element array specifying the number of values to skip
C                along each dimension of the 2D arrays
C    EDGE        Array specifying the number of values to read
C                along each dimension of the 2D arrays
C
C!OUTPUT PARAMETERS:
C    COPY_INFO   Result flag (0=success, -1=failure)
C    ERROR_TEXT  String describing cause of failure
C
C!REVISION HISTORY:
C
C!TEAM-UNIQUE HEADER:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
C!END
C-----------------------------------------------------------------------
      
      implicit none

c ... Include files
      include 'hdf.inc'
      
c ... Arguments
      character*(*) geo_file, out_file, error_text
      integer start(2), stride(2), edge(2)

c ... Local variables
      integer rtn, rtn_inp, rtn_out
      integer geo_id, out_id
      character*100 sds_name 
      integer rank, dimsizes(2), data_type, num_attrs
      integer latval_id(2), lonval_id(2),
     &        solzen_id(2), solazm_id(2),
     &        senzen_id(2), senazm_id(2),
     &        timval_id(2)
      
      integer num_pixels, num_lines, num_scans  
      integer line_index, pixel_index, scan_index
      integer line_count, pixel_count
      
      integer max_pixel
      parameter (max_pixel = 1499)
      double precision double_data_inp(0 : max_pixel)
      double precision double_data_out(0 : max_pixel)
      real float_data_inp(0 : max_pixel)
      real float_data_out(0 :max_pixel)
      integer*2 short_data_inp(0 : max_pixel)
      integer*2 short_data_out(0 : max_pixel)
      
      integer geo_start(2), geo_stride(2), geo_edge(2)
      integer out_start(2), out_stride(2), out_edge(2)
      
c ... External functions
      integer sfstart, sfselect, sfn2index, sfginfo, sfrdata, sfwdata, 
     &  sfendacc, sfend
      external sfstart, sfselect, sfn2index, sfginfo, sfrdata, sfwdata, 
     &  sfendacc, sfend
            
C-----------------------------------------------------------------------
C     INITIALIZATION
C-----------------------------------------------------------------------

c ... Set default return values
      copy_info = -1
      error_text = 'Unidentified error'

c ... Open the input file (DFACC_READ is defined in hdf.inc)
      geo_id = sfstart(geo_file, DFACC_READ)
      if (geo_id .eq. -1) then
        error_text = 'Could not open input geolocation file for reading'
        return
      endif
      
c ... Open the output file (DFACC_WRITE is defined in hdf.inc)
      out_id = sfstart(out_file, DFACC_WRITE)
      if (out_id .eq. -1) then
        error_text = 'Could not open output product file for writing'
        return
      endif
      
c ... Get number of pixels, lines, and scans
      latval_id(1) = sfselect(geo_id, sfn2index(geo_id, 'Latitude'))
      rtn = sfginfo(latval_id(1), sds_name, rank, dimsizes, data_type,
     &  num_attrs) 
      rtn = sfendacc(latval_id(1))
      num_pixels = dimsizes(1)
      num_lines  = dimsizes(2)
      num_scans  = num_lines / 10

c ... Get the SDS ID for each of the input arrays
      latval_id(1) = sfselect(geo_id, sfn2index(geo_id, 'Latitude'))
      lonval_id(1) = sfselect(geo_id, sfn2index(geo_id, 'Longitude'))
      solzen_id(1) = sfselect(geo_id, sfn2index(geo_id, 'SolarZenith'))
      solazm_id(1) = sfselect(geo_id, sfn2index(geo_id, 'SolarAzimuth'))
      senzen_id(1) = sfselect(geo_id, sfn2index(geo_id, 'SensorZenith'))
      senazm_id(1) = sfselect(geo_id, sfn2index(geo_id, 'SensorAzimuth'))
      timval_id(1) = sfselect(geo_id, sfn2index(geo_id, 'EV start time'))

c ... Check input SDS IDs
      if (latval_id(1) .eq. -1 .or.
     &    lonval_id(1) .eq. -1 .or.
     &    solzen_id(1) .eq. -1 .or.
     &    solazm_id(1) .eq. -1 .or.
     &    senzen_id(1) .eq. -1 .or.
     &    senazm_id(1) .eq. -1 .or.
     &    timval_id(1) .eq. -1) then
        error_text = 'Could not select input SDS array'
        return
      endif
                 
c ... Get the SDS ID for each of the output arrays
      latval_id(2) = sfselect(out_id, sfn2index(out_id, 'Latitude'))
      lonval_id(2) = sfselect(out_id, sfn2index(out_id, 'Longitude'))
      solzen_id(2) = sfselect(out_id, sfn2index(out_id, 'Solar_Zenith'))
      solazm_id(2) = sfselect(out_id, sfn2index(out_id, 'Solar_Azimuth'))
      senzen_id(2) = sfselect(out_id, sfn2index(out_id, 'Sensor_Zenith'))
      senazm_id(2) = sfselect(out_id, sfn2index(out_id, 'Sensor_Azimuth'))
      timval_id(2) = sfselect(out_id, sfn2index(out_id, 'Scan_Start_Time'))

c ... Check output SDS IDs
      if (latval_id(2) .eq. -1 .or.
     &    lonval_id(2) .eq. -1 .or.
     &    solzen_id(2) .eq. -1 .or.
     &    solazm_id(2) .eq. -1 .or.
     &    senzen_id(2) .eq. -1 .or.
     &    senazm_id(2) .eq. -1 .or.
     &    timval_id(2) .eq. -1) then
        error_text = 'Could not select output SDS array'
        return
      endif

C-----------------------------------------------------------------------
C     COPY GEOLOCATION ARRAYS
C-----------------------------------------------------------------------

c ... Initialize output line counter
      line_count = 0
      
c ... Loop through each line in the input file
      do line_index = start(2), (num_lines - 1), stride(2)

c ...   Set the input dimensions (1st dim = pixel, 2nd dim = line)
        geo_start(1)  = 0
        geo_start(2)  = line_index
        geo_stride(1) = 1
        geo_stride(2) = 1
        geo_edge(1)   = num_pixels
        geo_edge(2)   = 1

c ...   Set the output dimensions
        out_start(1)  = 0
        out_start(2)  = line_count
        out_stride(1) = 1
        out_stride(2) = 1
        out_edge(1)   = edge(1)
        out_edge(2)   = 1

c ...   Read and write the latitude data (type: float)
        rtn_inp = sfrdata(latval_id(1), geo_start, geo_stride, geo_edge,
     &    float_data_inp)
        pixel_count = 0
        do pixel_index = start(1), edge(1) * stride(1), stride(1)
          float_data_out(pixel_count) = float_data_inp(pixel_index)
          pixel_count = pixel_count + 1
        end do  
        rtn_out = sfwdata(latval_id(2), out_start, out_stride, out_edge,
     &    float_data_out)
        
c ...   Check read/write status
        if (rtn_inp .eq. -1 .or. rtn_out .eq. -1) then
          error_text = 'Error reading or writing latitude data'
          return
        endif
                
c ...   Read and write the longitude data (type: float)
        rtn_inp = sfrdata(lonval_id(1), geo_start, geo_stride, geo_edge,
     &    float_data_inp)
        pixel_count = 0
        do pixel_index = start(1), edge(1) * stride(1), stride(1)
          float_data_out(pixel_count) = float_data_inp(pixel_index)
          pixel_count = pixel_count + 1
        end do  
        rtn_out = sfwdata(lonval_id(2), out_start, out_stride, out_edge,
     &    float_data_out)

c ...   Check read/write status
        if (rtn_inp .eq. -1 .or. rtn_out .eq. -1) then
          error_text = 'Error reading or writing longitude data'
          return
        endif

c ...   Read and write the solar zenith data (type: short integer)
        rtn_inp = sfrdata(solzen_id(1), geo_start, geo_stride, geo_edge,
     &    short_data_inp)
        pixel_count = 0
        do pixel_index = start(1), edge(1) * stride(1), stride(1)
          short_data_out(pixel_count) = short_data_inp(pixel_index)
          pixel_count = pixel_count + 1
        end do  
        rtn_out = sfwdata(solzen_id(2), out_start, out_stride, out_edge,
     &    short_data_out)

c ...   Check read/write status
        if (rtn_inp .eq. -1 .or. rtn_out .eq. -1) then
          error_text = 'Error reading or writing solar zenith data'
          return
        endif

c ...   Read and write the solar azimuth data (type: short integer)
        rtn_inp = sfrdata(solazm_id(1), geo_start, geo_stride, geo_edge,
     &    short_data_inp)
        pixel_count = 0
        do pixel_index = start(1), edge(1) * stride(1), stride(1)
          short_data_out(pixel_count) = short_data_inp(pixel_index)
          pixel_count = pixel_count + 1
        end do  
        rtn_out = sfwdata(solazm_id(2), out_start, out_stride, out_edge,
     &    short_data_out)

c ...   Check read/write status
        if (rtn_inp .eq. -1 .or. rtn_out .eq. -1) then
          error_text = 'Error reading or writing solar azimuth data'
          return
        endif

c ...   Read and write the sensor zenith data (type: short integer)
        rtn_inp = sfrdata(senzen_id(1), geo_start, geo_stride, geo_edge,
     &    short_data_inp)
        pixel_count = 0
        do pixel_index = start(1), edge(1) * stride(1), stride(1)
          short_data_out(pixel_count) = short_data_inp(pixel_index)
          pixel_count = pixel_count + 1
        end do  
        rtn_out = sfwdata(senzen_id(2), out_start, out_stride, out_edge,
     &    short_data_out)

c ...   Check read/write status
        if (rtn_inp .eq. -1 .or. rtn_out .eq. -1) then
          error_text = 'Error reading or writing sensor zenith data'
          return
        endif

c ...   Read and write the sensor azimuth data (type: short integer)
        rtn_inp = sfrdata(senazm_id(1), geo_start, geo_stride, geo_edge,
     &    short_data_inp)
        pixel_count = 0
        do pixel_index = start(1), edge(1) * stride(1), stride(1)
          short_data_out(pixel_count) = short_data_inp(pixel_index)
          pixel_count = pixel_count + 1
        end do  
        rtn_out = sfwdata(senazm_id(2), out_start, out_stride, out_edge,
     &    short_data_out)

c ...   Check read/write status
        if (rtn_inp .eq. -1 .or. rtn_out .eq. -1) then
          error_text = 'Error reading or writing sensor azimuth data'
          return
        endif

c ...   Increment the output line counter
        line_count = line_count + 1
        
      end do

C-----------------------------------------------------------------------
C     SPECIAL CASE FOR TIME DATA
C-----------------------------------------------------------------------

c ... Initialize output line counter
      line_count = 0

c ... Loop through each scan in the input file 
c ... (geolocation file contains one time value per earth scan)
      do scan_index = 0, (num_scans - 1), 1

c ...   Set the input dimensions (1st dim = scan)
        geo_start(1)  = scan_index
        geo_stride(1) = 1
        geo_edge(1)   = 1

c ...   Set the output dimensions
        out_start(1)  = 0
        out_start(2)  = line_count
        out_stride(1) = 1
        out_stride(2) = 1
        out_edge(1)   = edge(1)
        out_edge(2)   = 1
        
c ...   Read and write the time data (type: double)
c ...   (note that two lines of output are written for every earth scan)
        rtn_inp = sfrdata(timval_id(1), geo_start, geo_stride, geo_edge,
     &    double_data_inp)
        do pixel_count = 0, (edge(1) - 1), 1
          double_data_out(pixel_count) = double_data_inp(0)
        end do  
        rtn_out = sfwdata(timval_id(2), out_start, out_stride, out_edge,
     &    double_data_out)
        out_start(2) = line_count + 1
        rtn_out = sfwdata(timval_id(2), out_start, out_stride, out_edge,
     &    double_data_out)

c ...   Check read/write status
        if (rtn_inp .eq. -1 .or. rtn_out .eq. -1) then
          error_text = 'Error reading or writing time data'
          return
        endif

c ...   Increment the output line counter
        line_count = line_count + 2

      end do      

C-----------------------------------------------------------------------
C     CLEANUP AND EXIT
C-----------------------------------------------------------------------

c ... End access to the input SDS IDs
      rtn = sfendacc(latval_id(1))
      rtn = sfendacc(lonval_id(1))
      rtn = sfendacc(solzen_id(1))
      rtn = sfendacc(solazm_id(1))
      rtn = sfendacc(senzen_id(1))
      rtn = sfendacc(senazm_id(1))
      rtn = sfendacc(timval_id(1))

c ... End access to the output SDS IDs
      rtn = sfendacc(latval_id(2))
      rtn = sfendacc(lonval_id(2))
      rtn = sfendacc(solzen_id(2))
      rtn = sfendacc(solazm_id(2))
      rtn = sfendacc(senzen_id(2))
      rtn = sfendacc(senazm_id(2))
      rtn = sfendacc(timval_id(2))

c ... Close the input file
      rtn = sfend(geo_id)

c ... Close the output file
      rtn = sfend(out_id)

c ... Set return values
      copy_info = 0
      error_text = 'Successful completion'
      
      END
