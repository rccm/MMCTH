      SUBROUTINE GET_ANCILLARY( LAT, LON, PRES, TEMP, MIXR, LAND,
     &  SFCTMP, PRMSL, PWAT, UGRD, VGRD, OZONE, ICEC, SST, NISE,
     &  MET_DATE, OZN_DATE, ICE_DATE, SST_DATE, NISE_DATE )

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c      Retrieve ancillary data items for a given latitude and longitude.
c
c !INPUT PARAMETERS:
c      LAT       Latitude (degrees, -90S to +90.0N)
c      LON       Longitude (degrees, -180W to +180E, Greenwich=0)
c
c !OUTPUT PARAMETERS:
c      PRES      Array of pressure levels (hPa)
c      TEMP      Array of atmospheric temperatures (K) at PRES(0:15)
c      MIXR      Array of water vapor mixing ratios (g/kg) at PRES(0:15)
c      LAND      Land mask (0=water, 1=land)
c      SFCTMP    Surface temperature (K)
c      PRMSL     Pressure (hPa) at mean sea level
c      PWAT      Precipitable water (g/cm**2)
c      UGRD      Surface wind u component (m/s)
c      VGRD      Surface wind v component (m/s)
c      OZONE     Total ozone (Dobsons)
c      ICEC      Ice concentration (fraction)
c      SST       Sea surface temperature (K) - valid over ocean only
c      NISE      NSIDC NISE snow/ice extent (see read_nise.f)
c      MET_DATE  UTC date for parameters PRES-VGRD (year,month,day,hour)
c      OZN_DATE  UTC date for parameter  OZONE     (year,month,day,hour)
c      ICE_DATE  UTC date for parameter  ICEC      (year,month,day,hour)
c      SST_DATE  UTC date for parameter  SST       (year,month,day,hour)
c      NISE_DATE UTC date for parameter  NISE      (year,month,day,hour)
c    
c      The missing value for output parameters is MISSING (see below).
c
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c      Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !DESIGN NOTES:
c      (1) On the first call, this subroutine unpacks and reads 4 files:
c          NCEP GDAS1 meteorological analysis data,
c          NCEP TOVS total ozone data,
c          NCEP SSMI ice concentration data,
c          NCEP sea surface temperature data,
c          NSIDC NISE snow/ice extent data.
c          On subsequent calls, data is obtained from SAVEd arrays.
c
c      (2) This subroutine will not cause an exit. If errors are
c          detected (e.g. missing or bad input file), the subroutine
c          will write a 'Recoverable error' message to the LogStatus
c          file, and will return missing value(s) for the parameter(s)
c          it failed to read.
c
c      (3) No checking of data validity is done within this routine.
c          Missing data values are used only where the input file was
c          either missing or bad. The user is responsible for checking
c          that ancillary data values (e.g. SST) are within an
c          acceptable range for user's application.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Input arguments
      
      REAL lat, lon

c ... Output arguments

      REAL pres(0:15), temp(0:15), mixr(0:15),
     &  land, sfctmp, prmsl, pwat, ugrd, vgrd, ozone, icec, sst
      INTEGER nise
      INTEGER met_date(4), ozn_date(4), ice_date(4), sst_date(4),
     &  nise_date(4)
      
c ... Parameters

      REAL missing
      PARAMETER ( missing = -999.0 )
      
c ... Local variables

      LOGICAL INIT
      
      INTEGER lun, i, j, k, level, ios, header( 8 ),
     &  pcfnum, reclen, status, version
         
      REAL x, x0, dx, y, y0, dy, p( 0:15 ), satmix,
     &  met_grid( 0:359, 0:180, 0:29 ),
     &  ozn_grid( 0:359, 0:180 ),
     &  ice_grid( 0:719, 0:359 ),
     &  sst_grid( 0:359, 0:179 )

      CHARACTER*160 errmsg

      LOGICAL met_success, ozn_success, ice_success, sst_success,
     &  nise_success      
      
      CHARACTER*255 nise_file
      
      INTEGER gridsize
      PARAMETER ( gridsize = 721 )      
      CHARACTER*1 nise_north( gridsize, gridsize )
      CHARACTER*1 nise_south( gridsize, gridsize )

      INTEGER met_year, met_month, met_day, met_hour
      INTEGER ozn_year, ozn_month, ozn_day, ozn_hour
      INTEGER ice_year, ice_month, ice_day, ice_hour
      INTEGER sst_year, sst_month, sst_day, sst_hour
      INTEGER nise_year, nise_month, nise_day, nise_hour
      INTEGER nise_minute, nise_second, nise_fsecond, rtn_code
      CHARACTER*255 ECS_DateTimeGroup, HDF_AttributeName


c ... Include files for PGS toolkit

      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_SMF.f'

c ... Include file for PCF numbers

      INCLUDE 'Atmos_AncData.inc'
                  
c ... External functions

      INTEGER modis_grib_driver
      EXTERNAL modis_grib_driver
      
      REAL svpwat
      EXTERNAL svpwat
      
      INTEGER openr_temp
      EXTERNAL openr_temp
      
      INTEGER pgs_io_gen_openf, pgs_pc_getreference
      EXTERNAL pgs_io_gen_openf, pgs_pc_getreference

      INTEGER read_nise
      EXTERNAL read_nise

      INTEGER ezlh_convert
      EXTERNAL ezlh_convert
                  
c ... Save statement

      SAVE

c ... Temperature and moisture profile pressure levels (hPa)

      DATA p / 1000.0, 925.0, 850.0, 700.0, 500.0, 400.0, 300.0,
     &   250.0, 200.0, 150.0, 100.0,  70.0,  50.0,  30.0,  20.0,  10.0 /

c ... Initialization flag

      DATA init / .true. /

c-----------------------------------------------------------------------
c     INITIALIZATION
c-----------------------------------------------------------------------

c ... Open and read input data files if this is the first call

      if ( init ) then

c ...   Set data ingest success/fail flags

        met_success = .false.
        ozn_success = .false.
        ice_success = .false.
        sst_success = .false.
        nise_success = .false.

c ...   Unpack GRIB files

        errmsg = ' '
        status = modis_grib_driver( 
     &    LUN_GDAS_0ZF, LUN_TEMP_GDAS_0ZF,
     &    LUN_OZ_DAILY, LUN_TEMP_OZ_DAILY,
     &    LUN_SEA_ICE,  LUN_TEMP_SEA_ICE,
     &    met_year, met_month, met_day, met_hour,
     &    ozn_year, ozn_month, ozn_day, ozn_hour,
     &    ice_year, ice_month, ice_day, ice_hour,
     &    errmsg )
        if ( status .ne. 0 ) then
          level = 1
          call message( 'get_ancillary', errmsg, status, level )
        endif

c ...   Open the unpacked met file

        pcfnum = LUN_TEMP_GDAS_0ZF
        reclen = 360*181*30*4
        status = openr_temp( pcfnum, reclen, lun )
        if ( status .ne. 0 ) then
          level = 1
          write( errmsg,'(''Error opening met PCF#'',i12)') pcfnum
          call message( 'get_ancillary', errmsg //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )
        endif

c ...   Read the unpacked met file

        if ( status .eq. 0 ) then
          read( lun, rec = 1, iostat = ios ) met_grid
          if ( ios .ne. 0 ) then
            level = 1
            write( errmsg,'(''Error opening met PCF#'',i12)') pcfnum
            call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', ios, level )
          else
            met_success = .true.
          endif
          call modis_io_gen_closef( pcfnum, lun )
        endif
        
c ...   Open the unpacked ozone file

        pcfnum = LUN_TEMP_OZ_DAILY
        reclen = 360*181*4
        status = openr_temp( pcfnum, reclen, lun )
        if ( status .ne. 0 ) then
          level = 1
          write( errmsg,'(''Error opening ozone PCF#'',i12)') pcfnum
          call message( 'get_ancillary', errmsg //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )
        endif

c ...   Read the unpacked ozone file

        if ( status .eq. 0 ) then
          read( lun, rec = 1, iostat = ios ) ozn_grid
          if ( ios .ne. 0 ) then
            level = 1
            write( errmsg,'(''Error reading ozone PCF#'',i12)') pcfnum
            call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', ios, level )
          else
            ozn_success = .true.
          endif
          call modis_io_gen_closef( pcfnum, lun )
        endif

c ...   Open the unpacked ice file

        pcfnum = LUN_TEMP_SEA_ICE
        reclen = 720*360*4
        status = openr_temp( pcfnum, reclen, lun )
        if ( status .ne. 0 ) then
          level = 1
          write( errmsg,'(''Error opening ice PCF#'',i12)') pcfnum
          call message( 'get_ancillary', errmsg //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )
        endif

c ...   Read the unpacked ice file

        if ( status .eq. 0 ) then
          read( lun, rec = 1, iostat = ios ) ice_grid
          if ( ios .ne. 0 ) then
            level = 1
            write( errmsg,'(''Error reading ice PCF#'',i12)') pcfnum
            call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', ios, level )
          else
            ice_success = .true.
          endif  
          call modis_io_gen_closef( pcfnum, lun )
        endif

c ...   Open the ASCII SST file (this is NOT a temporary file)

        pcfnum = LUN_REYNSST
        reclen = 1
        version = 1
        status = pgs_io_gen_openf( pcfnum, PGSd_IO_Gen_RSeqFrm, reclen,
     &    lun, version )
        if ( status .ne. PGS_S_SUCCESS ) then
          level = 1
          write( errmsg,'(''Error opening sst PCF#'',i12)') pcfnum
          call message( 'get_ancillary', errmsg //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )
        endif

c ...   Read the ASCII SST file

        if ( status .eq. PGS_S_SUCCESS ) then

          read( lun, '( 8i5 )', iostat = ios ) header
          if ( ios .ne. 0 ) then
            level = 1
            write( errmsg, 
     &        '(''Error reading sst header PCF#'',i12)') pcfnum
            call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', ios, level )
          else
            sst_year  = header( 1 )
            sst_month = header( 2 )
            sst_day   = header( 3 )
            sst_hour  = 0
          endif            

          read( lun, '( 20f4.2 )', iostat = ios ) sst_grid
          if ( ios .ne. 0 ) then
            level = 1
            write( errmsg,
     &        '(''Error reading sst data PCF#'',i12)') pcfnum
            call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', ios, level )
          else
            sst_success = .true.
          endif
          call modis_io_gen_closef( pcfnum, lun )

        endif

c ...   Get the name of the NISE snow/ice file

        pcfnum = LUN_NISE
        version = 1
        status = pgs_pc_getreference( pcfnum, version, nise_file )
        if ( status .ne. PGS_S_SUCCESS ) then
          level = 1
          write( errmsg,'(''Error opening NISE PCF#'',i12)') pcfnum
          call message( 'get_ancillary', errmsg //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )
        endif
        
c ...   Read the NISE snow/ice file (no date or time in this file yet!)

        status = read_nise( nise_file, gridsize, 
     &    nise_north, nise_south )
        if ( status .ne. 0 ) then
          level = 1
          write( errmsg,'(''Error reading NISE PCF#'',i12)') pcfnum
          call message( 'get_ancillary', errmsg //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )
        else
          version = 1
          HDF_AttributeName = 'coremetadata.0'
          ECS_DateTimeGroup = 'EndingDateTime'

          call Parse_ECS_DateTime(pcfnum, version,
     &                            HDF_AttributeName,
     &                            ECS_DateTimeGroup,
     &                            nise_year, nise_month, nise_day,
     &                            nise_hour, nise_minute, nise_second,
     &                            nise_fsecond, rtn_code)

          if ( rtn_code .ne. 0 ) then
            nise_success = .false.
            level = 1
            write( errmsg,'(''Error reading NISE metadata on PCF#'',i12)') pcfnum
            call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
          else
            nise_success = .true.
          endif !if ( rtn_code .ne. 0 )

        endif !if ( status .ne. 0 )
            
c ...   Unset initialization flag
        
        init = .false.
        
      endif

c-----------------------------------------------------------------------
c     SET MISSING VALUES
c-----------------------------------------------------------------------

c ... Set return values to missing

      do k = 0, 15
        pres( k ) = missing
        temp( k ) = missing
        mixr( k ) = missing
      end do
      land   = missing
      sfctmp = missing
      prmsl  = missing
      pwat   = missing
      ugrd   = missing
      vgrd   = missing
      ozone  = missing
      icec   = missing
      sst    = missing
      nise   = int( missing )
      do i = 1, 4
        met_date( i )  = int( missing )
        ozn_date( i )  = int( missing )
        ice_date( i )  = int( missing )
        nise_date( i ) = int( missing )
      end do
      
c-----------------------------------------------------------------------
c     GET MET AND OZONE DATA
c-----------------------------------------------------------------------

      if ( met_success .or. ozn_success ) then
      
c ...   Compute cell coordinates in met and ozn grids

        x = min( max( lon,  -179.99 ), 179.99 )
        if( x .lt. 0.0 ) x = x + 360.0
        x0 = 0.0
        dx = 1.0
        i = int( ( x - x0 + 0.5*dx ) / dx )
        if( i .eq. 360 ) i = 0
 
        y = min( max( lat, -89.99 ), 89.99 )
        y0 = 90.0
        dy = -1.0
        j = int( ( y - y0 + 0.5*dy ) / dy )

      endif

      if ( met_success ) then
      
c ...   Save output pressure levels

        do k = 0, 15
          pres( k ) = p( k )
        end do

c ...   Save output met data
c ...   (note that water vapor profile is relative humidity (%))

        do k = 0, 15
          temp( k ) = met_grid( i, j, k )
        end do
        do k = 0, 6
          mixr( k ) = met_grid( i, j, k + 16 )
        end do
        land   = met_grid( i, j, 23 )
        sfctmp = met_grid( i, j, 24 )
        prmsl  = met_grid( i, j, 25 ) * 0.01
        pwat   = met_grid( i, j, 26 )
        ugrd   = met_grid( i, j, 27 )
        vgrd   = met_grid( i, j, 28 )

c ...   Convert relative humidity profile (%) to mixing ratio (g/kg)

        do k = 0, 6
 
c ...     Compute mixing ratio at 100% relative humidity
 
          satmix = 622.0 * svpwat( temp( k ) ) / pres( k )

c ...     Convert relative humidity to mixing ratio

          mixr( k ) = satmix * 0.01 * mixr( k )

        end do

c ...   Extrapolate mixing ratio profile from 300 hPa to 10 hPa

        do k = 7, 15
          mixr( k ) = max( mixr( 6 ), 0.003 ) * ( pres( k ) / 300.0 )**3
          mixr( k ) = max( mixr( k ), 0.003 )
        end do

c ...   Save date

        met_date( 1 ) = met_year
        met_date( 2 ) = met_month
        met_date( 3 ) = met_day
        met_date( 4 ) = met_hour
        
      endif
      
      if ( ozn_success ) then

c ...   Save output ozone data

        ozone  = ozn_grid( i, j )
        
c ...   Save date

        ozn_date( 1 ) = ozn_year
        ozn_date( 2 ) = ozn_month
        ozn_date( 3 ) = ozn_day
        ozn_date( 4 ) = ozn_hour

      endif

c-----------------------------------------------------------------------
c     GET ICE DATA
c-----------------------------------------------------------------------

      if ( ice_success ) then
      
c ...   Compute cell coordinates in ice grid

        x = min( max( lon, -179.99 ), 179.99 )
        if( x .lt. 0.0 ) x = x + 360.0
        x0 = 0.25
        dx = 0.5
        i = int( ( x - x0 + 0.5*dx ) / dx )
        if( i .eq. 720 ) i = 0
 
        y = min( max( lat, -89.99 ), 89.99 )
        y0 = 89.75
        dy = -0.5
        j = int( ( y - y0 + 0.5*dy ) / dy )

c ...   Save output ice data

        icec = ice_grid( i, j )

c ...   Save date

        ice_date( 1 ) = ice_year
        ice_date( 2 ) = ice_month
        ice_date( 3 ) = ice_day
        ice_date( 4 ) = ice_hour

      endif

c-----------------------------------------------------------------------
c     GET SST DATA
c-----------------------------------------------------------------------

      if ( sst_success ) then
                        
c ...   Compute cell coordinates in sst grid

        x = min( max( lon, -179.99 ),  179.99 )
        x0 = -179.5
        dx = 1.0
        i = int( ( x - x0 + 0.5*dx ) / dx )
 
        y = min( max( lat, -89.99 ), 89.99 )
        y0 = -89.5
        dy = 1.0
        j = int( ( y - y0 + 0.5*dy ) / dy )

c ...   Save output sst data

        sst = sst_grid( i, j ) + 273.15

c ...   Save date

        sst_date( 1 ) = sst_year
        sst_date( 2 ) = sst_month
        sst_date( 3 ) = sst_day
        sst_date( 4 ) = sst_hour

c ...   Correct Y2K problem with SST year

        if ( sst_date( 1 ) .le. 99 ) then
          if ( sst_date( 1 ) .ge. 70 ) then
            sst_date( 1 ) = sst_date( 1 ) + 1900
          else
            sst_date( 1 ) = sst_date( 1 ) + 2000
          endif
        endif
        
      endif

c-----------------------------------------------------------------------
c     GET NISE DATA
c-----------------------------------------------------------------------

      if ( nise_success ) then

c ...   Get cell coordinates for southern or northern hemisphere
c ...   (Note the grid name strings are 'Sl' and 'Nl', L not 1)
      
        x = min( max( lon, -179.99 ),  179.99 )
        y = min( max( lat, -89.99 ), 89.99 )

        if ( y .lt. 0.0 ) then
          status = ezlh_convert( 'Sl', y, x, i, j )
        else
          status = ezlh_convert( 'Nl', y, x, i, j )
        endif

        if ( status .ne. 0 ) then

          level = 1
          call message( 'get_ancillary', 
     &      'Error converting lat,lon to NISE col,row' //
     &      ' [OPERATOR ACTION: CONTACT SDST]', status, level )

        else
        
c ...     Save output NISE data for southern or northern hemisphere
        
          if ( y .lt. 0.0 ) then
            nise = ichar( nise_south( i, j ) )
          else
            nise = ichar( nise_north( i, j ) )
          endif

        endif
        
c ...   Save date

        nise_date( 1 ) = nise_year
        nise_date( 2 ) = nise_month
        nise_date( 3 ) = nise_day
        nise_date( 4 ) = nise_hour

      endif
              
      END
