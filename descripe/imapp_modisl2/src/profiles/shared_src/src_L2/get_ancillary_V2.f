      SUBROUTINE GET_ANCILLARY( LAT, LON, PRES, TEMP, MIXR, LAND,
     &  SFCTMP, PRMSL, PWAT, UGRD, VGRD, OZONE, ICEC, SST )

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
C      PRES      Array of pressure levels (hPa)
C      TEMP      Array of atmospheric temperatures (K) at PRES(0:15)
C      MIXR      Array of water vapor mixing ratios (g/kg) at PRES(0:15)
C      LAND      Land mask (0=water, 1=land)
C      SFCTMP    Surface temperature (K)
C      PRMSL     Pressure (hPa) at mean sea level
C      PWAT      Precipitable water (g/cm**2)
C      UGRD      Surface wind u component (m/s)
C      VGRD      Surface wind v component (m/s)
C      OZONE     Total ozone (Dobsons)
C      ICEC      Ice concentration (fraction)
C      SST       Sea surface temperature (K) - valid over ocean only
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
C      sst.dat - native ASCII format NCEP sea surface temperature data
C
C      These input files are created prior to the PGE that calls this
C      subroutine.
C
C !END
C
C-----------------------------------------------------------------------

      implicit none

      include 'PGS_IO.f'
            
c ... input arguments
      
      real lat, lon

c ... output arguments

      real pres( 0:15), temp( 0:15 ), mixr( 0:15 ),
     &     land, sfctmp, prmsl, pwat, ugrd, vgrd, ozone, icec, sst
      
c ... local variables

      integer lun, i, j, k, level, ios, header( 0:7 ), init,
     &  pcfnum, reclen, version
         
      real x, x0, dx, y, y0, dy, p( 0:15 ), satmix,
     &  met_grid( 0:359, 0:180, 0:29 ),
     &  ozn_grid( 0:359, 0:180 ),
     &  ice_grid( 0:719, 0:359 ),
     &  sst_grid( 0:359, 0:179 )

c ... variables for grib converter

      integer status
      integer modis_grib_driver
      character*150 errmsg

c ... functions

      real svpwat
      external svpwat
      
c ... save statement

      save

c ... define temperature and moisture profile pressure levels (hPa)

      data p / 1000.0, 925.0, 850.0, 700.0, 500.0, 400.0, 300.0,
     &   250.0, 200.0, 150.0, 100.0,  70.0,  50.0,  30.0,  20.0,  10.0 /

c ... define initialization flag (1=initialize)

      data init / 1 /

c ... define logical unit number

      data lun / 21 /

c ... define error level

      data level / 2 /
      
c ... read input data files if this is the first call

      if ( init .eq. 1 ) then

c ...   Unpack GRIB files

        errmsg = ' '
        status = 0
        status = modis_grib_driver( 900000, 497000, 900020, 497020,
     &    900040, 497040, errmsg )
        if ( status .lt. 0 ) call message( 'get_ancillary',
     &    errmsg, status, level )

c ...   Read the unpacked met data

        pcfnum = 497000
        reclen = 360*181*30*4
        call modis_io_gen_temp_openf( pcfnum, PGSD_IO_GEN_RDIRUNF,
     &    reclen, lun )
        read( lun, rec = 1, iostat = ios ) met_grid
        if ( ios .ne. 0 ) call message( 'get_ancillary',
     &    'Meteorological ancillary data read failed' //
     &    ' [OPERATOR action: Contact SDST]', ios, level )
        call modis_io_gen_closef( pcfnum, lun )

c ...   Read the unpacked ozone data

        pcfnum = 497020
        reclen = 360*181*4
        call modis_io_gen_temp_openf( pcfnum, PGSD_IO_GEN_RDIRUNF,
     &    reclen, lun )
        read( lun, rec = 1, iostat = ios ) ozn_grid
        if ( ios .ne. 0 ) call message( 'get_ancillary',
     &    'Ozone ancillary data read failed' //
     &    ' [OPERATOR action: Contact SDST]', ios, level )
        call modis_io_gen_closef( pcfnum, lun )

c ...   Read the unpacked ice data
        
        pcfnum = 497040
        reclen = 720*360*4
        call modis_io_gen_temp_openf( pcfnum, PGSD_IO_GEN_RDIRUNF,
     &    reclen, lun )
        read( lun, rec = 1, iostat = ios ) ice_grid
        if ( ios .ne. 0 ) call message( 'get_ancillary',
     &    'Ice ancillary data read failed' //
     &    ' [OPERATOR action: Contact SDST]', ios, level )
        call modis_io_gen_closef( pcfnum, lun )

c ...   Read the ASCII SST data
        
        pcfnum = 900030
        reclen = 1
        version = 1
        call modis_io_gen_openf( pcfnum, PGSD_IO_GEN_RSEQFRM,
     &    reclen, lun, version )
        read( lun, '( 8i5 )', iostat = ios ) header
        if ( ios .ne. 0 ) call message( 'get_ancillary',
     &    'SST ancillary header read failed' //
     &    ' [OPERATOR action: Contact SDST]', ios, level )
        read( lun, '( 20f4.2 )', iostat = ios ) sst_grid
        if ( ios .ne. 0 ) call message( 'get_ancillary',
     &    'SST ancillary data read failed' //
     &    ' [OPERATOR action: Contact SDST]', ios, level )
        call modis_io_gen_closef( pcfnum, lun )
        
        init = 0
        
      endif
      
c ... save output pressure levels

      do k = 0, 15
        pres( k ) = p( k )
      end do

c ... compute cell coordinates in met and ozn grids

      x = min( max( lon,  -179.99 ), 179.99 )
      if( x .lt. 0.0 ) x = lon + 360.0
      x0 = 0.0
      dx = 1.0
      i = int( ( x - x0 + 0.5*dx ) / dx )
      if( i .eq. 360 ) i = 0
      
      y = min( max( lat, -89.99 ), 89.99 )     
      y0 = 90.0
      dy = -1.0
      j = int( ( y - y0 + 0.5*dy ) / dy )

c ... save output met and ozone data
c ... (note that water vapor profile is relative humidity (%))

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
      ozone  = ozn_grid( i, j )

c ... convert relative humidity profile (%) to mixing ratio (g/kg)

      do k = 0, 6
      
c ...   compute mixing ratio at 100% relative humidity
     
        satmix = 622.0 * svpwat( temp( k ) ) / pres( k )

c ...   convert relative humidity to mixing ratio

        mixr( k ) = satmix * 0.01 * mixr( k )

      end do

c ... extrapolate mixing ratio profile from 300 hPa to 10 hPa

      do k = 7, 15
        mixr( k ) = max( mixr( 6 ), 0.003 ) * ( pres( k ) / 300.0 )**3
        mixr( k ) = max( mixr( k ), 0.003 )
      end do

c ... compute cell coordinates in ice grid

      x = min( max( lon, -179.99 ), 179.99 )
      if( x .lt. 0.0 ) x = lon + 360.0
      x0 = 0.25
      dx = 0.5
      i = int( ( x - x0 + 0.5*dx ) / dx )
      if( i .eq. 720 ) i = 0
      
      y = min( max( lat, -89.99 ), 89.99 )     
      y0 = 89.75
      dy = -0.5
      j = int( ( y - y0 + 0.5*dy ) / dy )

c ... save output ice data

      icec = ice_grid( i, j )
            
c ... compute cell coordinates in sst grid

      x = min( max( lon, -179.99 ),  179.99 )
      x0 = -179.5
      dx = 1.0
      i = int( ( x - x0 + 0.5*dx ) / dx )
      
      y = min( max( lat, -89.99 ), 89.99 )     
      y0 = -89.5
      dy = 1.0
      j = int( ( y - y0 + 0.5*dy ) / dy )

c ... save output sst data

      sst = sst_grid( i, j ) + 273.15

      end
