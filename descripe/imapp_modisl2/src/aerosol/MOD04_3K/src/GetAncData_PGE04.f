      SUBROUTINE GetAncData_PGE04(lat,lon,sfctmp,ugrd,vgrd,pwat,ozone)

      implicit none

      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'Atmos_AncData.inc'

C-----------------------------------------------------------------------
C
C !F77
C
C !DESCRIPTION:
C      Retrieve ancillary data items for a given latitude and longitude.
C
C !INPUT PARAMETERS:
C      LAT       Latitude (degrees, -90S to +90.0N)
C      LON       Longitude (degrees, -180W to +180E, Greenwich=0)
C
C !OUTPUT PARAMETERS:
C      SFCTMP    Surface temperature (K)
C      UGRD      Surface wind u component (m/s)
C      VGRD      Surface wind v component (m/s)
C      PWAT      Precipitable water (g/cm**2)
C      OZONE     Total ozone (Dobsons)
C
C !REVISION HISTORY:
c 01/23/98 fhliang
c fixed prolog: added colons and '!Team-Unique Header:'.
c
C      11-DEC-1997 Liam Gumley CIMSS/SSEC
C          Changed open for SST ASCII file to MODIS_IO_GEN_OPENF.
C      10-DEC-1997 Liam Gumley CIMSS/SSEC
C          Added call to MODIS_GRIB_DRIVER which unpacks the GRIB input
C          files. Updated the PCF numbers of the resulting unpacked
C          files and the ASCII format SST file.
C      06-AUG-1997 Liam Gumley CIMSS/SSEC
C          Changed water vapor profile to mixing ratio (g/kg).
C          Interpolated water vapor profile up to 10 hPa via power law.
C      04-AUG-1997 Liam Gumley CIMSS/SSEC
C          Created.
C
c!Team-Unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
C !DESIGN NOTES:
C      This program reads four input files:
C      met.dat - unpacked binary NCEP GDAS1 meteorological analysis data
C      ozn.dat - unpacked binary NCEP TOVS total ozone data
C      ice.dat - unpacked binary NCEP SSMI ice concentration data
C
C      These input files are created prior to the PGE that calls this
C      subroutine.
C
C      Selected local variables for possible future use
C      ICEC      Ice concentration (fraction)
C      LAND      Land mask (0=water, 1=land)
C      MIXR      Array of water vapor mixing ratios (g/kg) at PRES(0:15)
C      TEMP      Array of atmospheric temperatures (K) at PRES(0:15)
C      PRES      Array of pressure levels (hPa)
C      PRMSL     Pressure (hPa) at mean sea level
C
C !END
C
C-----------------------------------------------------------------------

C-----Parameters
      REAL        missing
      PARAMETER ( missing = -999.0 )

C-----Function arguments
      REAL lat, lon, ozone, pwat, sfctmp, ugrd, vgrd

C-----Local variables
      CHARACTER*8   ESDT_name
      CHARACTER*160 errmsg
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) file_name

      INTEGER met_date(4), ozn_date(4), ice_date(4)
      INTEGER met_year, met_month, met_day, met_hour
      INTEGER ozn_year, ozn_month, ozn_day, ozn_hour
      INTEGER ice_year, ice_month, ice_day, ice_hour
c rh
      INTEGER n
      INTEGER i, ios, j, k, level, lun, pcfnum, reclen,
     1        status, version

      LOGICAL INIT
      LOGICAL met_success, ozn_success, ice_success

      REAL icec, land, prmsl, satmix, xlon,
     1     x, x0, dx, y, y0, dy
      REAL ice_grid( 0:719, 0:359 ),
     1     met_grid( 0:359, 0:180, 0:30 ),
     2     ozn_grid( 0:359, 0:180 )
      REAL p(0:15), pres(0:15), temp(0:15), mixr(0:15)


c ... External functions
      REAL     svpwat
      EXTERNAL svpwat

      INTEGER  ezlh_convert, modis_grib_driver, openr_temp,
     &         pgs_pc_getreference
      EXTERNAL ezlh_convert, modis_grib_driver, openr_temp,
     &         pgs_pc_getreference


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

c-----Open and read input data files if this is the first call and
c     set data ingest success/fail flags

      if ( init ) then
        met_success = .false.
        ozn_success = .false.
        ice_success = .false.


c-----------------------------------------------------------------------
c     Get NCEP Meteorological Data
c-----------------------------------------------------------------------
         version   = 1
         file_name = ' '
         errmsg    = ' '

         status = pgs_pc_getreference( LUN_GDAS_0ZF, version, file_name )

c--------met file entry exists in PCF
         if ( status .eq. PGS_S_SUCCESS ) then
            ESDT_name = 'GDAS_0ZF'
            status    = modis_grib_driver( LUN_GDAS_0ZF, LUN_TEMP_GDAS_0ZF,
     1                                     ESDT_name, errmsg,
     2                                     met_year, met_month,
     3                                     met_day, met_hour )

c-----------failed to open or read met grib product: system problem, exit immediately
            if ( status .ne. 0 ) then
               level = 2
               call message( 'GetAncData_PGE04', errmsg, status, level )

c-----------open gdas temp file
            else
               pcfnum = LUN_TEMP_GDAS_0ZF
               reclen = 360*181*31*4

               status = openr_temp( pcfnum, reclen, lun )

c--------------failed to open gdas temp file: system problem, exit immediately
               if ( status .ne. 0 ) then
                  level = 2
                  write( errmsg,'(''Error opening GDAS temp file on LUN'',i12)') pcfnum

                  call message( 'GetAncData_PGE04', errmsg //
     1            ' [OPERATOR ACTION: Contact SDST]', status, level )

c--------------read gdas temp file
               else
                  read( lun, rec = 1, iostat = ios ) met_grid

c-----------------failed to read gdas temp file: system problem, exit immediately
                  if ( ios .ne. 0 ) then
                     level = 2
                     write( errmsg,'(''Error reading GDAS temp file on LUN'',i12)')
     1               pcfnum
                     call message( 'GetAncData_PGE04', errmsg //
     1               ' [OPERATOR ACTION: Contact SDST]', ios, level )
                  else
                     met_success = .true.
                  endif ! failed to read gdas temp file

                  call modis_io_gen_closef( pcfnum, lun )
               endif ! failed to open gdas temp file

            endif ! failed to open grib file

c--------no met file entry in pcf; write warning message and run without
         else
            level = 1
            write( errmsg, '(''No GDAS file entry found on PCF LUN'',i12)') LUN_GDAS_0ZF
            call message( 'GetAncData_PGE04', errmsg, status, level )

         endif ! no met file entry in PCF 

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

      do i = 1, 4
        met_date( i )  = int( missing )
        ozn_date( i )  = int( missing )
        ice_date( i )  = int( missing )
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
        Ozone  = met_grid( i, j, 30 )

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

c      if ( ozn_success ) then

c ...   Save output ozone data

c        ozone  = ozn_grid( i, j )

c rh
 


c ...   Save date

c        ozn_date( 1 ) = ozn_year
c        ozn_date( 2 ) = ozn_month
c        ozn_date( 3 ) = ozn_day
c        ozn_date( 4 ) = ozn_hour

c      endif

cc-----------------------------------------------------------------------
cc     GET ICE DATA
cc-----------------------------------------------------------------------
c
c      if ( ice_success ) then
c
cc ...   Compute cell coordinates in ice grid
c
c        x = min( max( lon, -179.99 ), 179.99 )
c        if( x .lt. 0.0 ) x = x + 360.0
c        x0 = 0.25
c        dx = 0.5
c        i = int( ( x - x0 + 0.5*dx ) / dx )
c        if( i .eq. 720 ) i = 0
c
c        y = min( max( lat, -89.99 ), 89.99 )
c        y0 = 89.75
c        dy = -0.5
c        j = int( ( y - y0 + 0.5*dy ) / dy )
c
cc ...   Save output ice data
c
c        icec = ice_grid( i, j )
c
cc ...   Save date
c
c        ice_date( 1 ) = ice_year
c        ice_date( 2 ) = ice_month
c        ice_date( 3 ) = ice_day
c        ice_date( 4 ) = ice_hour
c
c      endif

      END
