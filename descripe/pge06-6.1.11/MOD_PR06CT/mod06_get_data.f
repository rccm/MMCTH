      subroutine mod06_get_data( cube_handle, cube_file, geo_handle,
     &  cld_handle, h_eco1, scan , begin_date, begin_time, met_date )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Get all input data at 1km resolution needed by the University
c    of Wisconsin MOD06 product for one scan.
c
c!Input Parameters:
c    CUBE_HANDLE    Integer array returned by OPMFIL for L1B file
c    CUBE_FILE      Full name and path of L1B file
c    GEO_HANDLE     Integer array returned by OPMFIL for geolocation
c    CLD_HANDLE     Integer array returned by OPMFIL for cloudmask
c    H_ECO1         1km north American exosystem file handle
c    SCAN           Scan number within L1B granule
c    BEGIN_DATE     Beginning date of data granule
c    BEGIN_TIME     Beginning time of data granule
c
c    MET_DATE       Year, month, day, and hour of gridded met data
c
c!Output Parameters:
c    The following arrays in COMMON /MOD06_DATA/ are filled:
c    RADIANCE1      Radiances for IR bands
c    LAT1           Latitude
c    LON1           Longitude
c    TEM1           Surface temperature
c    WAT1           Surface water vapor mixing ratio
c    OZN1           Total column ozone from NCEP
c    PRE1           Surface pressure
c    TPROF1         Temperatures at 101 standard levels (see array PSTD)
c    WPROF1         Water vapor mixing ratio at 101 standard levels
c    LAND1          Land/water mask
c    ECO            Ecosystem value based on Olson categorization
c    VIEW           Geolocation viewing zenith angle value
c    SOLZ           Geolocation solar zenith angle value
c    DN_FLAG        Pixel Day/night flag taken from instrument mode
c                    and solar zenith angle threshold
c    MASK1          Cloudmask
c    MASK1_QA       Cloudmask QA
c    CS_BIAS_CORR   Clear-sky radiance bias correction for CO2-slicing
c    MET_DATA_I     Met data grid cell for each 1-km MODIS pixel
c                    (coordinate in longitudinal direstion)
c    MET_DATA_J     Met data grid cell for each 1-km MODIS pixel
c                    (coordinate in latitudinal direction)
c
c!Revision History:
c    05-20-2021     R. Frey  Changed definition of night/day to depend 
c                            only on solar zenith angle.
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      save
      
      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'

c ... parameters

      real       EPS_ZERO
      parameter (EPS_ZERO = 0.00001)

c ... arguments

      integer cube_handle(*), geo_handle(*), cld_handle(*), 
     &        h_eco1, scan
      character*(*) cube_file, begin_date, begin_time

c ... local variables

      double precision delta, tai_time
      integer nise, resol, data_size( 2 ), met_date(4), ozn_date(4),
     &        ice_date(4), sst_date(4), nise_date(4),rtn,
     &        V_code(nxx, nyy, modis_bands)
      character*1 gain, scan_type
      character*6 cal_type
c rhucek 05/04/01: added work buffers charwork1 and charwork2 to hold substrings 
c                  of assumed-size character arguments begin_date and begin_time. 
      character*25 charwork1,charwork2
      logical error_flag, init, lsf(max_pixel, max_line)
      
      real buffer( max_pixel, max_line )
      byte buf_un( max_pixel, max_line ), buf_sa( max_pixel, max_line ),
     &  v_flag( max_pixel, max_line )

      integer i, j, k, ier, icode, mask_byte, ilat, grn_line
      integer icec_sst

      real pres( 26 ), temp( 26 ), mixr( 26 ), land, sfctmp, prmsl,
     &  prsfc, pwat, ugrd, vgrd, ozone, icec, sst, tozone, o3mr(6)
      real tlat(1), tlon(1),
     &     zn_biasdl(181,nbct), zn_biasnl(181,nbct),
     &     zn_biasoc(181,nbct)
      byte ecotest(1)

      integer ds_dim1_cm, ds_dim1_qa, ds_dim2, ds_dim3

      integer ni
      real pstd( 101 ), tstd( 101 ), wstd( 101 ), 
     &  latitude
      
      INTEGER year, month, day, hour, minute
      DOUBLE PRECISION second
      CHARACTER*27 ccsds_time

c ... set program name for error messaging
      character*32 routine_name
      parameter ( routine_name = 'mod06_get_data' )
      
c ... external functions

      integer read_geo_v2,pgs_td_utctotai,compare_times,date_to_ccsds
      external read_geo_v2,pgs_td_utctotai,compare_times,date_to_ccsds

c ... Save Statements
c whuang 05/09/01: Added common / mod06_data / and / mod06_debug / to 
c                  the SAVE statement
       SAVE init, ccsds_time, / mod06_data /, / mod06_debug /

c ... data statements

      DATA init / .true. /

c 03/10/98 fhliang commented out next 10 lines.
cc ... Define 1km IR band numbers
c
c      data bands/ 20, 22, 23, 24, 25, 27, 28, 29,
c     &  30, 31, 32, 33, 34, 35, 36 /
c
cc ... Defind 1km IR band numbers used by irphase
c      data ir_mband / 29, 31, 32 /
c
cc ... Defind 1km IR band array indices in band array
c      data ir_aband / 8, 10, 11 /

c ... Define 101 standard pressure levels (hPa)

      data pstd    / 0.0050,    0.0161,    0.0384,    0.0769,    0.1370,
     +    0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972,
     +    1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204,
     +    5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492,
     +   14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829,
     +   29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882,
     +   51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396,
     +   83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775,
     +  125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784,
     +  180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338,
     +  247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369,
     +  328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738,
     +  424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200,
     +  535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398,
     +  661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897,
     +  802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236,
     +  958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000/

c ... Initialize 1km radiance storage array with the bad value

      do k = 1, nbands
        do j = 1, max_line
          do i = 1, max_pixel
            radiance1( i, j, k ) = bad_value
          end do
        end do
      end do

c ... Extract the scan type indicating day or night mode from L1B file

      scan_type = ' '
      call mod06_Get_Swath_Metadata(cube_handle,scan,scan_type)

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
     &                 gain, resol, cal_type, max_pixel, max_line,
     &                 buffer, buf_un, buf_sa, v_flag, V_code,
     &                 data_size, error_flag )
        if ( error_flag ) then
          call message( routine_name,
     &      'Failed to extract 1km L1B data.' //
     &      ' [OPERATOR ACTION: Check input L1B, rerun PGE]',
     &      0, 1 )

        if ( debug .gt. 0 )
     &    write( h_output, '(''Error reading L1B radiance'')' )
        endif
        
c ...   Copy good data into 1km radiance storage array

        do j = 1, data_size( 2 )
          do i = 1, data_size( 1 )
            if ( v_flag( i, j ) .eq. 0 .and. buffer( i, j ) .gt. 0.0 )
     &        radiance1( i, j, k ) = buffer( i, j )
          end do
        end do
                  
      end do

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
      if ( ier .eq. -1 ) then
        call message( routine_name,
     &    'Failed to extract latitude from Geo. file' //
     &    ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     &    0, 1 )

        if ( debug .gt. 0 )
     &    write( h_output, '(''Error reading latitude'')' )
      endif
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( i, j ) .ge. -90.0 .and.
     &      buffer( i, j ) .le. 90.0 ) lat1( i, j ) = buffer( i, j )
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
      if ( ier .eq. -1 ) then
        call message( routine_name,
     &    'Failed to extract longitude from Geo. file' //
     &    ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     &    0, 1 )

        if ( debug .gt. 0 )
     &    write( h_output, '(''Error reading longitude'')' )
      endif
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( i, j ) .ge. -180.0 .and.
     &      buffer( i, j ) .le. 180.0 ) lon1( i, j ) = buffer( i, j )
        end do
      end do

c ... Get viewing zenith data 

      do j = 1, max_line
        do i = 1, max_pixel
          view( i, j ) = bad_value
          buffer( i, j ) = bad_value
        end do
      end do
      icode = 4
      ier = read_geo_v2( geo_handle, icode, scan, max_pixel,
     &  max_line, buffer, data_size )
      if ( ier .eq. -1 ) then
        call message( routine_name,
     &    'Failed to extract sensor zenith from Geo. file' //
     &    ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     &    0, 1 )

        if ( debug .gt. 0 )
     &    write( h_output, '(''Error reading viewing zenith '')' )
      endif
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if (abs(buffer(i,j)) .lt. 67.0 .and. abs(buffer(i,j)).gt. 0.0)
     &      view( i, j ) = buffer( i, j )
        end do
      end do

c ... Get solar zenith angle and day night flag

      do j = 1, max_line
        do i = 1, max_pixel
          solz( i, j ) = bad_value
          buffer( i, j ) = bad_value
          dn_flag( i, j ) = -1
        end do
      end do
      icode = 7
      ier = read_geo_v2( geo_handle, icode, scan, max_pixel,
     &  max_line, buffer, data_size )
      if ( ier .eq. -1 ) then
        call message( routine_name,
     &    'Failed to extract solar zenith from Geo. file' //
     &    ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     &    0, 1 )

        if ( debug .gt. 0 )
     &    write( h_output, '(''Error reading solar zenith '')' )
      endif
      do j = 1, data_size( 2 )
        grn_line = ((scan - 1) * 10) + j
        do i = 1, data_size( 1 )
            solz( i, j ) = buffer( i, j )
            if (solz(i,j) .gt. 0.0) then
c              if (scan_type .ne. 'D' .or. solz(i,j) .gt. 85.0) then
               if (solz(i,j) .gt. 85.0) then
                 dn_flag(i,j) = 0
               else 
                 dn_flag(i,j) = 1
               endif
               grn_dn_flag(i, grn_line) = dn_flag(i,j)
            endif
        end do
      end do

c ... Get land/water mask for 1km bands

      do j = 1, max_line
        do i = 1, max_pixel
          land1( i, j ) = -1
          eco( i, j ) = -1
          buffer( i, j ) = bad_value
        end do
      end do
      icode = 9
      ier = read_geo_v2( geo_handle, icode, scan, max_pixel,
     &  max_line, buffer, data_size )
      if ( ier .eq. -1 ) then
        call message( routine_name,
     &    'Error reading landmask from Geolocation file' //
     &    ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     &    0, 1 )
        if ( debug .gt. 0 )
     &    write( h_output, '(''Error reading landmask'')' )
      endif        
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( i, j ) .ge. 0.0 .and.
     &      buffer( i, j ) .le. 10.0 ) land1( i, j ) =
     &      int( buffer( i, j ) )
        end do
      end do

c-----------------------------------------------------------------------
c     CHECK TIME OF NON-MODIS ANCILLARY DATA
c-----------------------------------------------------------------------

c ... On first pass, check time of ancillary data

      if ( init ) then

c ...   Convert granule start date/time to TAI seconds

        tai_time = 0.0d0
        
c rhucek 05/04/01: replaced references to substrings of assumed-size 
c character variables begin_date and begin_time with local fixed-size
c variables charwork1 and charwork2. 
c       rtn = pgs_td_utctotai( begin_date(1:10)//'T'//
c    &    begin_time(1:15)//'Z', tai_time )

        charwork1 = begin_date(1:10)
        charwork2 = begin_time(1:15)
        rtn = pgs_td_utctotai(charwork1(1:10)//'T'//
     &    charwork2(1:15)//'Z', tai_time )
        if ( rtn .ne. 0 ) call message( 'routine_name',
     &    'Error converting granule start date/time to TAI seconds' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2 )

c ...   Get ancillary data to obtain date arrays

        call get_ancillary( lat1( 1, 1 ), lon1( 1, 1 ),
     &    pres, temp, mixr, land, sfctmp, prmsl, prsfc, pwat, tozone, o3mr, ugrd, vgrd,
     &    ozone, icec, sst, icec_sst, nise,met_date, ozn_date, ice_date, 
     &    sst_date, nise_date, zn_biasoc, zn_biasdl, zn_biasnl)

c ...   Check if NCEP GDAS1 data was read correctly by
c ...   testing MET_DATE year field for validity (year must be GE 0)

        if ( met_date( 1 ) .lt. 0 ) call message( 'routine_name',
     &    'Failed to read NCEP GDAS1 ancillary data' // char(10) //
     &    '[OPERATOR ACTION: Stage correct NCEP GDAS1 data and ' //
     &    're-run PGE. If error persists, contact SDST]', 0, 2 )

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

        if ( rtn .ne. 0 ) call message( 'routine_name',
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
        ccsds_time = '  '
        rtn = date_to_ccsds( year, month, day, hour, minute, second,
     &    ccsds_time )
        if ( rtn .ne. 0 ) call message( 'routine_name',
     &    'Error converting NCEP GDAS1 analysis time to CCSDS format' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 1 )

c----------------------------------------------------------------------
c       Okay, now read the ecosystem information

c ...   First, let's make sure we are reading from the correct
c ...   file.  Extract Madison Wisconsin eco_type from file
c       call read_goode( 1,1,43.083,-89.305,h_eco1,ecotest )
        tlat(1) = 43.083
        tlon(1) = -89.305
        call read_goode( 1,1,tlat,tlon,h_eco1,ecotest )
c whuang 05/11/01: changed value 14 to 22
        if ( ecotest(1) .ne. 22 ) call message( 'get_anc_data',
     &    'Extracted incorrect ecosystem value from 1 km file. ' //
     &    '[OPERATOR ACTION: Verify size and format of file. ' //
     &    'If error persists, contact SDST.]', 0, 2 )

c----------------------------------------------------------------------
c ...   Unset initialization flag

        init = .false.

      endif

c-----------------------------------------------------------------------
c     READ NON-MODIS ANCILLARY DATA
c-----------------------------------------------------------------------

c ... Read all ecosystem values for this scan cube

       call read_goode( data_size( 1 ), data_size( 2 ),
     &   lat1, lon1, h_eco1, eco )


c ... Get ancillary data for 1km bands

      do j = 1, max_line
        do i = 1, max_pixel

c ...     Set default output values

          tem1( i, j ) = bad_value
          wat1( i, j ) = bad_value
          ozn1( i, j ) = bad_value
          pre1( i, j ) = bad_value
          do k = 1, 101
            tprof1( k, i, j ) = bad_value
            wprof1( k, i, j ) = bad_value
          end do
          do k = 1, 6
            o3prof1( k, i, j ) = bad_value
          end do
          
c ...     Fill clear-sky bias correction array with missing values.
          do k = 1,nbct
            cs_bias_corr( i, j, k ) = bad_value
          end do

c ...     If lat/lon data is ok, get ancillary data

          if ( lat1( i, j ) .ge. -90.0 .and.
     &      lat1( i, j ) .le.  90.0 .and.
     &      lon1( i, j ) .ge. -180.0 .and.
     &      lon1( i, j ) .le.  180.0 ) then

            call get_ancillary( lat1( i, j ), lon1( i, j ),
     &        pres, temp, mixr, land, sfctmp, prmsl, prsfc, pwat, tozone, o3mr, ugrd, vgrd,
     &        ozone, icec, sst , icec_sst, nise, met_date, ozn_date, ice_date,
     &        sst_date, nise_date, zn_biasoc, zn_biasdl, zn_biasnl)

c ...       Fill land/sea flag array.
            if(land1( i, j ) .eq. 1 .or. land1( i, j ) .eq. 2 .or.
     &         land1( i, j ) .eq. 4 ) then
              lsf( i, j ) = .true.
            else
              lsf( i, j ) = .false.
            end if

c ...       Set surface data items if values are ok

            if (sfctmp .ge. 0.0) then
                tem1( i, j ) = sfctmp
            endif
            if (land1(i,j) .ne. 1 .and. land1(i,j) .ne. 2 .and.
     &          land1(i,j) .ne. 4 .and. icec_sst .le. 50) then
                tem1( i, j ) = sst
            endif

            if ( mixr( 1 ) .ge. 0.0 ) wat1( i, j ) = mixr( 1 )
            if ( tozone .ge. 0.0 ) ozn1( i, j ) = tozone
            if ( prsfc .ge. 0.0 ) pre1( i, j ) = prsfc
            if ( prmsl .ge. 0.0 ) premsl( i, j ) = prmsl

c ...       Interpolate temperature and water vapor mixing ratio profiles
c ...       to the 101 standard levels.
            ni = 26
            latitude = lat1( i, j )
            call extend_profile( ni, pres, temp, mixr, latitude,
     &                           tstd, wstd )

c ...       Copy temperature and moisture profiles into output array
            
            do k = 1, 101

              if ( tstd( k ) .ge. 0.0 ) tprof1( k, i, j ) = tstd( k )

c rhucek 10/23/01: changed logic for resetting wprof1
c              if ( abs( wstd( k ) ) .lt. EPS_ZERO) then
               if ( wstd( k ) .lt. EPS_ZERO .and. wstd( k) .ge. -0.50 ) then
                  wprof1( k, i, j ) = EPS_ZERO
              elseif ( wstd( k ) .ge. EPS_ZERO ) then
                  wprof1( k, i, j ) = wstd( k )
              endif

            end do

c ...       Copy ozone profile data to output array.
            do k = 1, 6
              o3prof1(k,i,j) = o3mr(k)
            enddo

c ...       Get clear-sky radiance biases.

            ilat = (90.5 - lat1( i, j ) ) + 1

            if( lsf(i,j) ) then

              if( dn_flag( i , j ) .eq. 1) then
                do k = 1, nbct
                  if(zn_biasdl(ilat,k) .ne. -999.99) then
                    cs_bias_corr( i, j, k ) = zn_biasdl( ilat , k )
                  end if
                end do
              else 
                do k = 1, nbct
                  if(zn_biasnl(ilat,k) .ne. -999.99) then
                    cs_bias_corr( i, j, k ) = zn_biasnl( ilat , k )
                  end if
                end do
              end if

            else

              do k = 1, nbct
                if(zn_biasoc(ilat,k) .ne. -999.99) then
                  cs_bias_corr( i, j, k ) = zn_biasoc( ilat , k )
                end if
              end do

              if(lat1(i,j) .le. -60.0 .and. icec_sst .gt. 0.5) then
                do k = 1, nbct
                  if(zn_biasdl(ilat,k) .ne. -999.99) then
                    cs_bias_corr( i, j, k ) = zn_biasdl(ilat,k)
                  end if
                end do
              end if

            end if

c ...       If this is a water pixel, then make sure ecosystem is
c ...       is not saying land - set equal to 14
            if (land1(i,j).eq.0.or.land1(i,j).eq.3.or.land1(i,j).eq.5 
     +        .or. land1(i,j).eq.6 .or. land1(i,j).eq.7) 
     +        eco(i,j) = 14
          endif

        end do
      end do
     
c ... Interpolate clear-sky radiance bias grid.

c     call csr_interp ( cs_bias_corr, lsf )

c ... Get cloud mask result array and QA array for 1km bands

      do k = 1, max_line
        do j = 1, max_pixel
          do i = 1, nmask
            mask1( i, j, k ) = 0
          end do
        end do
      end do
      do k = 1, max_line
        do j = 1, max_pixel
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
      if ( error_flag ) then
        call message( routine_name, 
     &    'Error reading cloudmask file' //
     &    ' [OPERATOR ACTION: Check cloud mask file, rerun PGE]',
     &    0, 1 )
        if ( debug .gt. 0 )
     &    write( h_output, '(''Error reading cloudmask'')' )
      endif

c ... Write debug information

      if ( debug .gt. 0 ) then

        write( h_output, '(/,a,/)' ) 'mod06_get_data debug'

        i = 676
        j = 6

        write( h_output, '( ''Data for Pixel at '', ' //
     &    ' ''Scan'',i4,'', Pixel'',i4,'', Line'',i2)' ) scan, i, j

        write( h_output, '(''IR Bands: '',20i7)' )
     &    ( bands( k ), k = 1, nbands )
        write( h_output, '(''Radiance: '',20f7.3)' )
     &    ( radiance1( i, j, k ), k = 1, nbands )

        mask_byte = int( mask1( 1, i, j ) )
        write( h_output,
     &    ' (''Lat ='',f9.3,'', Lon ='',f10.3,'', Landmask ='',i3, ' //
     &    ' '', Cloudmask Bit 0 ='',i2,'', Cloudmask Bits 1:2 ='',i3)')
     &    lat1( i, j ), lon1( i, j ), land1( i, j ),
     &    iand( mask_byte, 1 ),
     &    iand( mask_byte, 6 )

        write( h_output, '(''Ancillary Data: '', ' //
     &    ' ''Sfc Temp ='',f7.2,'' K, Mixing Ratio ='',f7.2, ' //
     &    ' '' g/kg, Total Ozone ='',f7.2,'' Dobsons,'', ' //
     &    ' '' Sfc Pressure ='',f7.2,'' hPa'', '' time '',a27)' )
     &    tem1( i, j ), wat1( i, j ), ozn1( i, j ), pre1( i, j ),
     &    ccsds_time

        write( h_output, '(''Auxillary Data: '', ' //
     &    ' ''Viewing Zenith ='',f7.2,'' Ecosystem Type ='',i10 )')
     &    view( i, j ), eco( i, j )

        write( h_output, '(''Auxillary Data: '', ' //
     &    ' ''Solar Zenith ='',f7.2,'' Day/Night Flag  ='',i10 )')
     &    solz( i, j ), dn_flag( i, j )

         write(h_output,'(1x)')
         write(h_output,'(''Interpolated and extrapolated 101-level profiles:'')')
         write( h_output, '(''Pres(hPa)  Temp(K)  Mixr(g/kg)'')' )
         do k = 1, 101
            write( h_output, '(f9.1,2x,f7.2,2x,f10.3)' )
     &      pstd(k), tprof1(k,i,j), wprof1(k,i,j)
         end do

      endif

      end
