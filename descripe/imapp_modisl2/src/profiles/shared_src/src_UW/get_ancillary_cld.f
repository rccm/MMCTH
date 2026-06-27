      SUBROUTINE GET_ANCILLARY( LAT, LON, PRES, TEMP, MIXR, LAND,
     &  SFCTMP, PRMSL, PWAT, UGRD, VGRD, OZONE, ICEC, SST, NISE,
     &  MET_DATE, OZN_DATE, ICE_DATE, SST_DATE, NISE_DATE, PXL_BIAS )

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
c      PXL_BIAS  Mean difference between observed and modeled clear-sky 
c                radiances (currently not used)
c    
c      The missing value for output parameters is MISSING (see below).
c
c !REVISION HISTORY:
c
c  09/12/03 R. Frey:  Added calls to pgs_pc_getreference to check for
c    the presence of grib met and sea ice files, and also for NISE files.
c    Changed error levels to "2" if any files cannot be opened or read.
c
c  12/11/02 R. Hucek: Commented code statements that access the NCEP 
c    Reynolds sst and TOVS ozone products.  
c
c  11/21/02 R. Frey:  Added subroutine read_reynsst to enable reading of
c    either formatted or unformatted Reynolds SST files.
c
c  06/04 Collection 5  R. Frey:  Added logic to acquire clear-sky
c    radiance biases.
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

      REAL pres(0:25), temp(0:25), mixr(0:25),
     &  land, sfctmp, prmsl, pwat, ugrd, vgrd, ozone, icec, sst
      INTEGER nise
      INTEGER met_date(4), ozn_date(4), ice_date(4), sst_date(4),
     &  nise_date(4)
      
c ... Parameters

      REAL missing
      PARAMETER ( missing = -999.0 )
      
      INTEGER npoints_x, npoints_y
      PARAMETER ( npoints_x = 360 )
      PARAMETER ( npoints_y = 180 )

      INTEGER modnum, csr_nmbins, csr_nmstats, csr_nmbnds
      PARAMETER ( modnum = 5 )
      PARAMETER ( csr_nmbins = 814880 )
      PARAMETER ( csr_nmstats = 9 )
      PARAMETER ( csr_nmbnds = 5 )

c ... Local variables

      LOGICAL INIT
      
      INTEGER lun, i, j, k, level, ios, ret, grid_index,
     &  pcfnum, reclen, status, version, modfil(modnum),
     &  mcsbnds(csr_nmbnds)
         
      REAL x, x0, dx, y, y0, dy, p( 0:25 ), satmix, xlon, 
     &  met_grid( 0:359, 0:180, 0:53 ),
     &  ozn_grid( 0:359, 0:180 ),
     &  ice_grid( 0:719, 0:359 ),
     &  sst_grid( 0:npoints_x-1, 0:npoints_y-1 ), sst_bl
      real csr_data(csr_nmstats,csr_nmbins,csr_nmbnds),
     &     pxl_bias(csr_nmbnds),num_csr,sum_csr

      CHARACTER*3   sst_file_fmt
      CHARACTER*8   ESDT_name
      character*10 value_text
      CHARACTER*160 errmsg

      LOGICAL met_success, ozn_success, ice_success, sst_success,
     &  nise_success, csr_success
      
      CHARACTER*255 nise_file, grib_name, CSR_file
      
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
      INCLUDE 'PGS_SMF.f'

c ... Include file for PCF numbers
      INCLUDE 'Atmos_AncData.inc'
                  
c ... External subroutines
      EXTERNAL bl_int
      EXTERNAL blint_met
      EXTERNAL read_csr

c ... External functions

      INTEGER modis_grib_driver
      EXTERNAL modis_grib_driver
      
      REAL modis_bright
      EXTERNAL modis_bright
      
      REAL ppv
      EXTERNAL ppv
      
      INTEGER openr_temp
      EXTERNAL openr_temp
      
      INTEGER pgs_pc_getreference
      EXTERNAL pgs_pc_getreference

      INTEGER read_nise
      EXTERNAL read_nise

      INTEGER ezlh_convert
      EXTERNAL ezlh_convert
                  
      INTEGER get_grid_index
      EXTERNAL get_grid_index
                  
c ... Save statement.
      SAVE

c ... Temperature and moisture profile pressure levels (hPa)

      DATA p / 1000.0, 975.0, 950.0, 925.0, 900.0, 850.0, 800.0,
     &   750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0, 400.0,
     &   350.0, 300.0, 250.0, 200.0, 150.0, 100.0,  70.0,  50.0,
     &    30.0,  20.0,  10.0 /

c ... Initializations.

      DATA init / .true. /
      DATA mcsbnds / 31, 33, 34, 35, 36 /

c-----------------------------------------------------------------------
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
        csr_success = .false.


c-----------------------------------------------------------------------
c      Get NCEP Meteorological Data
c-----------------------------------------------------------------------

c       Check presence of grib met file.
        version = 1
        errmsg    = ' '
        status  = pgs_pc_getreference( LUN_GDAS_0ZF, version, grib_name )
        if ( status .eq. PGS_S_SUCCESS ) then

c ...     Unpack grib met file and write to binary file
          ESDT_name = 'GDAS_0ZF'
          status = modis_grib_driver( LUN_GDAS_0ZF, LUN_TEMP_GDAS_0ZF, 
     1                                ESDT_name, errmsg,
     2                                met_year, met_month, 
     3                                met_day, met_hour ) 

          if ( status .ne. 0 ) then
            level = 2
            call message( 'get_ancillary', errmsg, status, level )

c ....    Open unpacked met file
          else
            pcfnum = LUN_TEMP_GDAS_0ZF
            reclen = 360*181*54*4
            status = openr_temp( pcfnum, reclen, lun )

            if ( status .ne. 0 ) then
              level = 2
              write( errmsg,'(''Error opening GDAS on PCF#'',i12)') 
     &           pcfnum
              call message( 'get_ancillary', errmsg //
     &           ' [OPERATOR ACTION: Contact SDST]', status, level )

c ......    Read the unpacked met file
            else 
              if ( status .eq. 0 ) then
                read( lun, rec = 1, iostat = ios ) met_grid

                if ( ios .ne. 0 ) then
                  level = 2
                  write( errmsg,'(''Error reading GDAS on PCF#'',i12)') 
     &               pcfnum
                  call message( 'get_ancillary', errmsg //
     &               ' [OPERATOR ACTION: Contact SDST]', ios, level )
                else
                  met_success = .true.
                endif

                call modis_io_gen_closef( pcfnum, lun )
              endif   
            endif   
          endif   

        else

          level = 1
          write( errmsg,'(''No entry found for GDAS on PCF#'',
     &       i12)') LUN_GDAS_0ZF
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )

        end if

cc-----------------------------------------------------------------------
cc     Get TOVS Ozone Data
cc-----------------------------------------------------------------------
cc ...   Unpack grib ozone file and write to binary file
c        errmsg    = ' '
c        ESDT_name = 'OZ_DAILY'
c        status    = modis_grib_driver( LUN_OZ_DAILY, LUN_TEMP_OZ_DAILY, 
c     1                                 ESDT_name, errmsg, 
c     2                                 ozn_year, ozn_month, 
c     3                                 ozn_day, ozn_hour ) 
c
c        if ( status .ne. 0 ) then
c          level = 1
c          call message( 'get_ancillary', errmsg, status, level ) 
c
cc ....  Open unpacked ozone file
c        else    
c          pcfnum = LUN_TEMP_OZ_DAILY
c          reclen = 360*181*4
c          status = openr_temp( pcfnum, reclen, lun )
c
c          if ( status .ne. 0 ) then
c            level = 1
c            write( errmsg,'(''Error opening ozone PCF#'',i12)') pcfnum
c            call message( 'get_ancillary', errmsg //
c     &           ' [OPERATOR ACTION: Contact SDST]', status, level ) 
c
cc ......  Read the unpacked ozone file
c          else    
c            if ( status .eq. 0 ) then
c              read( lun, rec = 1, iostat = ios ) ozn_grid
c
c              if ( ios .ne. 0 ) then
c                level = 1
c                write( errmsg,'(''Error opening ozone PCF#'',i12)') 
c     &             pcfnum
c                call message( 'get_ancillary', errmsg //
c     &             ' [OPERATOR ACTION: Contact SDST]', ios, level )
c              else
c                ozn_success = .true.
c              endif
c
c              call modis_io_gen_closef( pcfnum, lun )
c            endif   
c          endif   
c        endif   


c-----------------------------------------------------------------------
c     Get SSM/I sea ice concentration 
c-----------------------------------------------------------------------
c
c       Check presence of grib sea ice concentration file.
        version = 1
        ESDT_name = 'SEA_ICE'
        errmsg    = ' '
        status  = pgs_pc_getreference( LUN_SEA_ICE, version, grib_name )
        if ( status .eq. PGS_S_SUCCESS ) then

c ...     Unpack grib sea ice file and write to binary file
          status = modis_grib_driver( LUN_SEA_ICE, LUN_TEMP_SEA_ICE, 
     1                                ESDT_name, errmsg, 
     2                                ice_year, ice_month, 
     3                                ice_day, ice_hour ) 

          if ( status .ne. 0 ) then
            level = 2
            call message( 'get_ancillary', errmsg, status, level ) 

          else    
c ......    Open unpacked ice file
            pcfnum = LUN_TEMP_SEA_ICE
            reclen = 720*360*4 
            status = openr_temp( pcfnum, reclen, lun )

            if ( status .ne. 0 ) then
              level = 2
              write( errmsg,'(''Error opening sea ice on PCF#'',i12)') 
     &           pcfnum
              call message( 'get_ancillary', errmsg //
     &           ' [OPERATOR ACTION: Contact SDST]', status, level ) 

            else    
c ........    Read the unpacked ice file
              if ( status .eq. 0 ) then
                read( lun, rec = 1, iostat = ios ) ice_grid

                if ( ios .ne. 0 ) then
                  level = 2
                  write( errmsg,'(''Error reading sea ice on PCF#'',
     &               i12)') pcfnum
                  call message( 'get_ancillary', errmsg //
     &              ' [OPERATOR ACTION: Contact SDST]', ios, level ) 
                else    
                  ice_success = .true.
                endif   

                call modis_io_gen_closef( pcfnum, lun )
              endif   
            endif  
          endif   

        else

          level = 1
          write( errmsg,'(''No entry found for sea ice on PCF#'',i12)')
     &       LUN_SEA_ICE
          call message( 'get_ancillary', errmsg //
     &       ' [OPERATOR ACTION: Contact SDST]', status, level )

        end if


cc----------------------------------------------------------------------
cc     Get Reynolds sea surface temperature 
cc----------------------------------------------------------------------
 
         call read_reynsst ( npoints_x, npoints_y, sst_year, sst_month, 
     &                       sst_day, sst_hour, sst_grid, sst_file_fmt, 
     &                       sst_success)
 
        if( (.not. sst_success) ) then
          level  =  1
          status = -1
          write( errmsg,'(''Problem in read_reynsst '')')
          call message( 'get_ancillary', errmsg //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )
        end if

c----------------------------------------------------------------------
c     Get NISE snow extent 
c----------------------------------------------------------------------
c ...   Get the name of the NISE snow/ice file
        pcfnum  = LUN_NISE
        version = 1
        status  = pgs_pc_getreference( pcfnum, version, nise_file )
        if ( status .eq. PGS_S_SUCCESS ) then

c ...     Read the NISE snow/ice file 
          status = read_nise( nise_file, gridsize, 
     &      nise_north, nise_south )

          if ( status .ne. 0 ) then
            level = 2
            write( errmsg,'(''Error opening NISE on PCF#'',i12)') 
     &         pcfnum
            if ( status .lt. -1)  
     &      write( errmsg, '(''Error reading NISE on PCF#'',i12)') 
     &         pcfnum
            call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )

          else
            version = 1
            HDF_AttributeName = 'coremetadata.0'
            ECS_DateTimeGroup = 'EndingDateTime'

            call Parse_ECS_DateTime(pcfnum, version,
     &                              HDF_AttributeName,
     &                              ECS_DateTimeGroup,
     &                              nise_year, nise_month, nise_day,
     &                              nise_hour, nise_minute, nise_second,
     &                              nise_fsecond, rtn_code)

            if ( rtn_code .ne. 0 ) then
              nise_success = .false.
              level = 2
              write( errmsg,'(''Error parsing NISE metadata on PCF#'',
     &           i12)') pcfnum
              call message( 'get_ancillary', errmsg //
     &           ' [OPERATOR ACTION: Contact SDST]', status, level )
            else
              nise_success = .true.
            endif

          endif
            
        else

          level = 1
          write( errmsg,'(''No entry found for NISE on PCF#'',i12)') 
     &       pcfnum
          call message( 'get_ancillary', errmsg //
     &       ' [OPERATOR ACTION: Contact SDST]', status, level )

        end if

c----------------------------------------------------------------------
c     Get clear-sky radiance data. 
c----------------------------------------------------------------------

c       pcfnum = LUN_CSR8d
c       status = pgs_pc_getreference( pcfnum, version, CSR_file )

c       if( status .ne. PGS_S_SUCCESS ) then
          
c         level = 1
c         call message( 'get_ancillary',
c    &      'Error getting CSR filename from PCF ' //
c    &      ' [OPERATOR ACTION: Contact SDST]', status, level )

c       else

c         Open and read CSR data.

c         call read_csr(CSR_file,csr_nmbins,csr_nmstats,csr_nmbnds,
c    &                  modnum,csr_data,csr_success)

c         if( (.not. csr_success) ) then
c           level = 1
c           status = -1
c           write( errmsg,'(''Problem in read_csr '')')
c           call message( 'get_ancillary', errmsg //
c    &        ' [OPERATOR ACTION: Contact SDST]', status, level )
c         end if

c       endif

c----------------------------------------------------------------------

c ...   Unset initialization flag
        
        init = .false.
        
      endif

c-----------------------------------------------------------------------
c     SET MISSING VALUES
c-----------------------------------------------------------------------

c ... Set return values to missing

      ozone  = missing
      icec   = missing
      sst    = missing
      nise   = int( missing )
      do i = 1, 4
        ozn_date( i )  = int( missing )
        ice_date( i )  = int( missing )
        nise_date( i ) = int( missing )
      end do
      do i = 1, csr_nmbnds
        pxl_bias(i) = missing
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
      
c ...     Set met return values to missing
          do k = 0, 25
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
          do k = 1, 4
            met_date( k )  = int( missing )
          end do

c ...     Save output pressure levels

          do k = 0, 25
            pres( k ) = p( k )
          end do

c ...     Save output met data
c ...     (note that water vapor profile is relative humidity (%))

          do k = 0, 25
            temp( k ) = met_grid( i, j, k )
          end do
          do k = 0, 20
            mixr( k ) = met_grid( i, j, k + 26 )
          end do
          land   = met_grid( i, j, 47 )
               
          call blint_met( met_grid(0,0,48), i, j, y, x, sfctmp, ret )
          if(ret .lt. 0) then
          level = 2
          status = -1
          write( errmsg,'(''Problem in blint_met '')')
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
          end if

          call blint_met( met_grid(0,0,49), i, j, y, x, prmsl, ret )
          if(ret .lt. 0) then
          level = 2
          status = -1
          write( errmsg,'(''Problem in blint_met '')')
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
          end if
          prmsl = prmsl * 0.01

          pwat   = met_grid( i, j, 50 )
          ugrd   = met_grid( i, j, 51 )
          vgrd   = met_grid( i, j, 52 )

c ...     Convert relative humidity profile (%) to mixing ratio (g/kg)

          do k = 0, 20
 
c ...       Compute mixing ratio at 100% relative humidity
 
            satmix = 622.0 * ppv( temp( k ) ) / pres( k )

c ...       Convert relative humidity to mixing ratio

            mixr( k ) = satmix * 0.01 * mixr( k )

          end do

c ...     Extrapolate mixing ratio profile from 100 hPa to 10 hPa

          do k = 20, 25
            mixr( k ) = max( mixr( 20 ), 0.003 ) * ( pres( k ) / 100.0 )**3
            mixr( k ) = max( mixr( k ), 0.003 )
          end do

c ...     Save date

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

c       Binary SST data is shifted 180 degrees relative to ASCII data.
        if(sst_file_fmt .eq. 'fmt') then
          x = min( max( lon, -179.99 ),  179.99 )
          x0 = -179.5
          xlon = lon + 180.0
        else if(sst_file_fmt .eq. 'unf') then
          if(lon .lt. 0.0) then
            xlon = lon + 360.0
          else
            xlon = lon
          end if
          x = min( max( xlon, 0.00 ),  359.99 )
          x0 = 0.5
        end if

        dx = 1.0
        i = int( ( x - x0 + 0.5*dx ) / dx )
 
        y = min( max( lat, -89.99 ), 89.99 )
        y0 = -89.5
        dy = 1.0
        j = int( ( y - y0 + 0.5*dy ) / dy )

c ...   Bi-linearly interpolate SST

        call bl_int(sst_grid, i, j, lat, xlon, sst_bl, ret)
        if(ret .lt. 0) then
          level = 2
          status = -1
          write( errmsg,'(''Problem in bl_int '')')
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
        end if

c ...   Save output sst data

        sst = sst_bl + 273.15

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

c-----------------------------------------------------------------------
c     GET CLEAR RADIANCE BIAS DATA
c-----------------------------------------------------------------------

c     if ( csr_success ) then

c       Get CSR grid index for the latitude and logitude of current
c       pixel.

c       status = get_grid_index( lat, lon, grid_index )
c       if ( status .ne. 0 ) then
c         write( value_text, '( i10 )' ) status
c         call message( 'clear_radiance',
c    &      'Bad result code returned from get_grid_index (' //
c    &      value_text // '); returning to caller [OPERATOR ACTION:' //
c    &      'Contact SDST]', 0, 1 )
c         return
c       endif

c       Get number of clear sky observations and sum of observed clear
c       minus calculated clear radiances for 5 bands.  Then calculate
c       mean bias in the 25-km bin.
c       do i = 1, csr_nmbnds

c         num_csr = csr_data(2,grid_index,i)
c         if(num_csr .ge. 1.0) then
c           sum_csr = csr_data(8,grid_index,i)
c           pxl_bias(i) = sum_csr / num_csr
c         end if

c       end do

c     end if

c-----------------------------------------------------------------------

      return
      end
