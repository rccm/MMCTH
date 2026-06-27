      SUBROUTINE MOD07_GET_DATA( CUBE_HANDLE, CUBE_FILE, GEO_HANDLE,
     &  CLD_HANDLE, SCAN, BEGIN_DATE, BEGIN_TIME )

c-----------------------------------------------------------------------
c
c !F77
c
c !DESCRIPTION:
c     Get all input data at 1km resolution needed by MOD07 for one scan.
c
c !INPUT PARAMETERS:
c     CUBE_HANDLE    Integer array returned by OPMFIL for L1B file
c     CUBE_FILE      Full name and path of L1B file
c     GEO_HANDLE     Integer array returned by OPMFIL for geolocation
c     CLD_HANDLE     Integer array returned by OPMFIL for cloudmask
c     SCAN           Scan number within L1B granule
c
c !OUTPUT PARAMETERS:
c     The following arrays in COMMON /MOD07_DATA/ are filled:
c     RADIANCE1      Radiances for IR bands
c     LAT1           Latitude
c     LON1           Longitude
c     ELV1           Surface elevation
c     ZEN1           Sensor zenith angle
c     TEM1           Surface temperature
c     WAT1           Surface water vapor mixing ratio
c     PRE1           Surface pressure
c     LAND1          Land/water mask
c     MASK1          Cloudmask
c     MASK1_QA       Cloudmask QA
c     PWAT1          GDAS Total Precipitable Water 
c
c !REVISION HISTORY:
c     March 2003: SWS added PWAT1
c     March 2004: SWS defined v_code array and added to call list for read_l1b
c         (to accommodate new level 1b reader v2.2)
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE
      
      include 'debug.inc'
      include 'mod07_data.inc'

c ... Arguments

      INTEGER cube_handle(*), geo_handle(*), cld_handle(*), scan
      CHARACTER*(*) cube_file, begin_date, begin_time

c ... Local variables

      INTEGER resol, data_size( 2 )
      CHARACTER*1 gain
      CHARACTER*6 cal_type
      LOGICAL error_flag
      
      INTEGER v_code( max_pixel, max_line, LAST_MODIS_BAND )
      REAL buffer( max_pixel, max_line )
      BYTE buf_un( max_pixel, max_line ), buf_sa( max_pixel, max_line ),
     &  v_flag( max_pixel, max_line )

      INTEGER i, j, k, ier, icode, mask_byte

      REAL pres( 16 ), temp( 16 ), mixr( 16 ), land, sfctmp, sfcpres,
     &  pwat, ugrd, vgrd, ozone, icec, sst
      INTEGER nise
      INTEGER met_date(4), ozn_date(4), ice_date(4), sst_date(4),
     &  nise_date(4)

      INTEGER ds_dim1_cm, ds_dim1_qa, ds_dim2, ds_dim3

      INTEGER rtn
      DOUBLE PRECISION tai_time, delta

      LOGICAL init

      INTEGER year, month, day, hour, minute
      DOUBLE PRECISION second
      CHARACTER*27 ccsds_time

C-----SSTG member lipo, /05/27/99, added two character variables:
C-----               word_begin_time, word_begin_date
      character*100 word_begin_time, word_begin_date

      
c ... External functions

      INTEGER read_geo_v2
      EXTERNAL read_geo_v2

      INTEGER pgs_td_utctotai
      EXTERNAL pgs_td_utctotai

      INTEGER compare_times
      EXTERNAL compare_times
      
      INTEGER date_to_ccsds
      EXTERNAL date_to_ccsds
      
c ... Save statements

      SAVE init, ccsds_time
            
c ... Data statements

      DATA init / .true. /
      
c-----------------------------------------------------------------------
c     READ LEVEL 1B 1KM RADIANCE DATA FOR ALL REQUIRED BANDS
c-----------------------------------------------------------------------

c ... Initialize 1km radiance storage array with the bad value
       call message( 'mod07_get_data',
     &  'before initialize ', 0, 3 )


      do k = 1, nbands
        do j = 1, max_line
          do i = 1, max_pixel
            radiance1( i, j, k ) = bad_value
          end do
        end do
      end do

c ... Start loop over 1km IR bands

      do k = 1, nbands

c ...   Initialize buffer arrays

        do j = 1, max_line
          do i = 1, max_pixel
            buffer( i, j ) = bad_value
            v_flag( i, j ) = -1
          end do
        end do
            
c ...   Get single band of radiance data

        gain = 'L'
        resol = 1
        cal_type = 'rad'
        error_flag = .false.

        call read_l1b( cube_handle, cube_file, scan, bands( k ), 
     &    gain, resol, cal_type, max_pixel, max_line,
     &    buffer, buf_un, buf_sa, v_flag, v_code, data_size, 
     &    error_flag )

        if ( error_flag ) call message( 'mod07_get_data',
     &      'Error reading MODIS L1B radiance ' //
     &      '[OPERATOR ACTION: Verify format of staged file. ' //
     &      'If error persists, contact SDST]', 0, 2 )
        
c ...   Copy good data into 1km radiance storage array

        do j = 1, data_size( 2 )
          do i = 1, data_size( 1 )
            if ( v_flag( i, j ) .eq. 0 .and. buffer( i, j ) .gt. 0.0 )
     &        radiance1( i, j, k ) = buffer( i, j )
          end do
        end do
                  
      end do

c-----------------------------------------------------------------------
c     READ LEVEL 1B 1KM GEOLOCATION DATA
c-----------------------------------------------------------------------

c ... Get latitude data for 1km bands

      do j = 1, max_line
        do i = 1, max_pixel
          lat1( i, j ) = bad_value
          buffer( i, j ) = bad_value
        end do
      end do

      icode = 1
      ier = read_geo_v2( geo_handle, icode, scan, max_pixel,
     &  max_line, buffer, data_size )
      if ( ier .eq. -1 ) call message( 'mod07_get_data',
     &    'Error reading MODIS L1B geolocation latitude ' //
     &    '[OPERATOR ACTION: Verify format of staged file. ' //
     &    'If error persists, contact SDST]', 0, 2 )

      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( i, j ) .ge. -90.0 .and.
     &         buffer( i, j ) .le.  90.0 )
     &      lat1( i, j ) = buffer( i, j )
        end do
      end do

c ... Get longitude data for 1km bands

      do j = 1, max_line
        do i = 1, max_pixel
          lon1( i, j ) = bad_value
          buffer( i, j ) = bad_value
        end do
      end do

      icode = 2
      ier = read_geo_v2( geo_handle, icode, scan, max_pixel,
     &  max_line, buffer, data_size )
      if ( ier .eq. -1 ) call message( 'mod07_get_data',
     &    'Error reading MODIS L1B geolocation longitude ' //
     &    '[OPERATOR ACTION: Verify format of staged file. ' //
     &    'If error persists, contact SDST]', 0, 2 )

      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( i, j ) .ge. -180.0 .and.
     &         buffer( i, j ) .le.  180.0 )
     &      lon1( i, j ) = buffer( i, j )
        end do
      end do

c ... Get surface elevation data for 1km bands

      do j = 1, max_line
        do i = 1, max_pixel
          elv1( i, j ) = bad_value
          buffer( i, j ) = bad_value
        end do
      end do

      icode = 3
      ier = read_geo_v2( geo_handle, icode, scan, max_pixel,
     &  max_line, buffer, data_size )
      if ( ier .eq. -1 ) call message( 'mod07_get_data',
     &    'Error reading MODIS L1B geolocation surface elevation ' //
     &    '[OPERATOR ACTION: Verify format of staged file. ' //
     &    'If error persists, contact SDST]', 0, 2 )

c ...  (The Dead Sea is the lowest place on the Earth's surface
c ...  at about 400 meters (1,300 feet) below sea level.
c ...  Mt. Everest is the highest at about 8840 meters (29,000 ft).)

      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( i, j ) .ge.  -400.0 .and.
     &         buffer( i, j ) .le.  8840.0 )
     &      elv1( i, j ) = buffer( i, j )
        end do
      end do

c ... Get sensor zenith data for 1km bands

      do j = 1, max_line
        do i = 1, max_pixel
          zen1( i, j ) = bad_value
          buffer( i, j ) = bad_value
        end do
      end do

      icode = 4
      ier = read_geo_v2( geo_handle, icode, scan, max_pixel,
     &  max_line, buffer, data_size )
      if ( ier .eq. -1 ) call message( 'mod07_get_data',
     &    'Error reading MODIS L1B geolocation sensor zenith ' //
     &    '[OPERATOR ACTION: Verify format of staged file. ' //
     &    'If error persists, contact SDST]', 0, 2 )

      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( i, j ) .ge.   0.0 .and.
     &         buffer( i, j ) .le.  70.0 )
     &      zen1( i, j ) = buffer( i, j )
        end do
      end do

c ... Get land/water mask for 1km bands

      do j = 1, max_line
        do i = 1, max_pixel
          land1( i, j ) = -1
          buffer( i, j ) = bad_value
        end do
      end do

      icode = 9
      ier = read_geo_v2( geo_handle, icode, scan, max_pixel,
     &  max_line, buffer, data_size )
      if ( ier .eq. -1 ) call message( 'mod07_get_data',
     &    'Error reading MODIS L1B geolocation landmask ' //
     &    '[OPERATOR ACTION: Verify format of staged file. ' //
     &    'If error persists, contact SDST]', 0, 2 )

      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( i, j ) .ge. 0.0 .and.
     &      buffer( i, j ) .le. 10.0 ) land1( i, j ) =
     &      int( buffer( i, j ) )
        end do
      end do

c-----------------------------------------------------------------------
c     READ LEVEL 1B CLOUD MASK DATA
c-----------------------------------------------------------------------

c ... Get cloud mask result array and QA array for 1km bands

      do k = 1, max_line

        do j = 1, max_pixel

          do i = 1, nmask
            mask1( i, j, k ) = 0
          end do

          do i = 1, nmask_qa
            mask1_qa( i, j, k ) = 0
          end do

        end do

      end do

      error_flag = .false.
      call read_cldmsk( cld_handle, scan,
     &  nmask, nmask_qa, max_pixel, max_line,
     &  ds_dim1_cm, ds_dim1_qa, ds_dim2, ds_dim3,
     &  mask1, mask1_qa, error_flag )

      if ( error_flag ) call message( 'mod07_get_data',
     &  'Error reading MODIS L1B cloudmask ' //
     &  '[OPERATOR ACTION: Verify format of staged file. ' //
     &  'If error persists, contact SDST]', 0, 2 )

c-----------------------------------------------------------------------
c     CHECK TIME OF NON-MODIS ANCILLARY DATA
c-----------------------------------------------------------------------

c ... On first pass, check time of ancillary data

C-----SSTG member lipo, /05/27/99, added two assignment statement
      word_begin_date = begin_date(1:10)
      word_begin_time = begin_time(1:15)

      if ( init ) then
      
c ...   Convert granule start date/time to TAI seconds

        rtn = pgs_td_utctotai( word_begin_date(1:10)//'T'//
     &    word_begin_time(1:15)//'Z', tai_time )
        if ( rtn .ne. 0 ) call message( 'mod07_get_data',
     &    'Error converting granule start date/time to TAI seconds' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2 )
        
c ...   Get ancillary data to obtain date arrays

        call get_ancillary( lat1( 1, 1 ), lon1( 1, 1 ),
     &    pres, temp, mixr, land, sfctmp, sfcpres, pwat, ugrd, vgrd,
     &    ozone, icec, sst, nise,
     &    met_date, ozn_date, ice_date, sst_date, nise_date )

c ...   Check if NCEP GDAS1 data was read correctly by
c ...   testing MET_DATE year field for validity (year must be GE 0)

        if ( met_date( 1 ) .lt. 0 ) call message( 'mod07_get_data',
     &    'Failed to read NCEP GDAS1 ancillary data' // char(10) //
     &    '[OPERATOR ACTION: Stage correct NCEP GDAS1 data and ' //
     &    're-run PGE. If error persists, contact SDST]', 0, 1 )
       
c ...   Compare granule start time (TAI seconds) to NCEP GDAS1
c ...   date/time only (MOD07 does not use OZN, ICE, SST, NISE data).
c ...   NCEP GDAS1 analysis files are produced at 0Z, 6Z, 12Z, 18Z
c ...   daily. Thus if we specify the closest file (in time) must be
c ...   used, the DELTA value (largest acceptable difference between
c ...   granule start time and ancillary data time) should be no more
c ...   than 3 hours, or 3*3600 seconds.

        delta = 3.0d+0 * 3600.0d+0
        rtn = compare_times( met_date, tai_time, delta )

c ...   Print warning message if difference is greater than 3 hours

        if ( rtn .ne. 0 ) call message( 'mod07_get_data',
     &    'L1B granule start time and NCEP GDAS1 analysis time ' //
     &    'differ by more than 3 hours ' // char(10) //
     &    '[OPERATOR ACTION: Stage correct NCEP GDAS1 data and ' //
     &    're-run PGE. If error persists, contact SDST]', 0, 1 )

c ...   Create NCEP GDAS1 CCSDS time string for use in debug output

        year   = met_date( 1 )
        month  = met_date( 2 )
        day    = met_date( 3 )
        hour   = met_date( 4 )
        minute = 0
        second = 0.0d+0
        rtn = date_to_ccsds( year, month, day, hour, minute, second,
     &    ccsds_time )
        if ( rtn .ne. 0 ) call message( 'mod07_get_data',
     &    'Error converting NCEP GDAS1 analysis time to CCSDS format' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 1 )

c ...   Unset initialization flag

        init = .false.
        
      endif

c-----------------------------------------------------------------------
c     READ NON-MODIS ANCILLARY DATA
c-----------------------------------------------------------------------

c ... Get ancillary data for 1km bands

      do j = 1, max_line
      
        do i = 1, max_pixel

c ...     Set output values to missing

          tem1( i, j ) = bad_value
          wat1( i, j ) = bad_value
          pre1( i, j ) = bad_value
          pwat1( i, j ) = bad_value

c ...     If lat/lon is in range, get ancillary data

          if ( lat1( i, j ) .ge.  -90.0 .and.
     &         lat1( i, j ) .le.   90.0 .and.
     &         lon1( i, j ) .ge. -180.0 .and.
     &         lon1( i, j ) .le.  180.0 ) then

            call get_ancillary( lat1( i, j ), lon1( i, j ), 
     &        pres, temp, mixr, land, sfctmp, sfcpres, pwat, ugrd, vgrd,
     &        ozone, icec, sst, nise,
     &        met_date, ozn_date, ice_date, sst_date, nise_date )
              
c ...       Get interpolated surface pressure and TPW for this lat/lon

            call get_sfcpres_pwat( lat1( i, j ), lon1( i, j ), sfcpres, pwat )           

c ...       If data values are good, save input values of
c ...       surface temperature, water vapor, and pressure
c ...       (note: only surface pressure is used in retrieval)

            if ( sfctmp .ge. 0.0 )    tem1( i, j ) = sfctmp
            if ( mixr( 1 ) .ge. 0.0 ) wat1( i, j ) = mixr( 1 )
            if ( sfcpres .ge. 0.0 )   pre1( i, j ) = sfcpres
            if ( pwat .ge. 0.0 )   pwat1( i, j ) = pwat

          endif

        end do

      end do
     
c-----------------------------------------------------------------------
c     WRITE DEBUG INFORMATION
c-----------------------------------------------------------------------

      if ( debug .eq. 1 ) then

        write( h_output, '(/,a,/)' ) 'MOD07_GET_DATA DEBUG INFO'

        i = 676
        j = 6
        write( h_output, '( ''Data for Nadir Pixel at '', ' //
     &    ' ''Scan'',i4,'', Pixel'',i4,'', Line'',i2)' ) scan, i, j

        mask_byte = int( mask1( 1, i, j ) )
        write( h_output,
     &    ' (''Lat ='',f9.3,'', Lon ='',f10.3, ' //
     &    ' '', Elev = '',f7.1, '', Landmask ='',i3, ' //
     &    ' '', Cloudmask Bit 0 ='',i2,'', Cloudmask Bits 1:2 ='',i3)' )
     &    lat1( i, j ), lon1( i, j ), elv1( i, j ), land1( i, j ),
     &    iand( mask_byte, 1 ),
     &    iand( mask_byte, 6 )

        write( h_output, '(''Ancillary Data: '', a27, ' //
     &    ' '', Sfc Temp ='',f7.2,'' K, Mixing Ratio ='',f7.2, ' //
     &    ' '' g/kg, Sfc Pressure ='',f7.2,'' hPa'')' )
     &    ccsds_time, tem1( i, j ), wat1( i, j ), pre1( i, j )

        write( h_output, '(''IR Bands: '',20i8)' )
     &    ( bands( k ), k = 1, nbands )

        write( h_output, '(''Radiance: '',20f8.2)' )
     &    ( radiance1( i, j, k ), k = 1, nbands )

      endif

      END
