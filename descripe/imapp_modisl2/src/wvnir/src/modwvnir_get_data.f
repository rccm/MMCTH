C--------------------------------------------------------------------
C  Copyright (C) 2002,  Space Science and Engineering Center, 
C  University C  of Wisconsin-Madison, Madison WI.
C      
C  This program is free software; you can redistribute it 
C  and/or modify it under the terms of the GNU General 
C  Public License as published by the Free Software Foundation; 
C  either version 2 of the License, or (at your option) any 
C  later version.
C
C  This program is distributed in the hope that it will be 
C  useful, but WITHOUT ANY WARRANTY; without even the implied 
C  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
C  See the  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public 
C  License along with this program; if not, write to the Free 
C  Software Foundation, Inc., 59 Temple Place, Suite 330, 
C  Boston, MA  02111-1307 USA
C--------------------------------------------------------------------
C
C
      subroutine modwvnir_get_data(l1b_1km_lun,geo_1km_lun,datatype_1km,
     &           datatype_geo,interleave_1km,interleave_geo,
     &           resolution_1km,resolution_geo,
     &           offset_1km,offset_geo,samples_1km,
     &           samples_geo,lines_1km,lines_geo,
     &           error_1km,error_geo,bands_1km,bands_geo,
     &           bandnames_1km,bandnames_geo,bandunits_1km,
     &           bandunits_geo,scan, debug)

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Get all input data at 1km resolution needed by the 
C    Direct Broadcast NIR water vapor product.
c
c!Input Parameters:
C    l1b_1km_lun     LUN for open 1km L1B data file
C    geo_1km_lun     LUN for open geolocation 1km file
C    datatype_1km    Format of data in L1b file
C    datatype_geo    Format of data in geolocation file
C    interleave_1km  Order data is saved in L1b file
C    interleave_geo  Order data is saved in geolocation file
C    resolution_1km  L1b file data resolution
C    resolution_geo  Geolocation file data resolution
C    offset_1km      Offset, if any, of 1km L1B data
C    offset_geo      Offset, if any, of geolocation file data
C    samples_1km     Number of elements in 1km L1b file
C    samples_geo     Number of elements in geolocation file
C    lines_1km       Number of lines of data in 1km L1b file
C    lines_geo       Number of lines of data in geolocation file
C    error_1km       Bad data value of data in 1km L1b file
C    error_geo       Bad data value of data in geolocation file
C    bands_1km       Number of bands in 1km L1b file
C    bands_geo       Number of bands in geolocation file
C    bandnames_1km   Band numbers in 1km L1b file
C    bandnames_geo   Band numbers in geolocation file
C    bandunits_1km   Units for each band in the 1km L1b data file
C    bandunits_geo   Units for each band in the geolocation file
c    scan            Scan number within L1B granule
C    debug           0: no debug messages / 1: debug messages
c
c!Output Parameters:
c    The following arrays in COMMON /MODWVNIR_DATA/ are filled:
c    RADIANCE1      Radiances for NIR bands
c    VIEW           Viewing zenith angle value
c    SOLZ           Solar zenith angle value
c    SENA           Viewing azimuth angle value
c    SOLZ           Solar azimuth angle value
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      
      include 'modwvnir_data.inc'

c ... arguments

      integer l1b_1km_lun,geo_1km_lun,
     +        scan,datatype_1km,datatype_geo,
     +        resolution_1km,resolution_geo,offset_1km,
     +        offset_geo,samples_1km,samples_geo,
     +        lines_1km,lines_geo,bands_1km,bands_geo

      real error_1km,error_geo
 
      character*(*) interleave_1km,bandnames_1km(*),bandunits_1km(*),
     +              interleave_geo,bandnames_geo(*),bandunits_geo(*)

c ... local variables

      integer nise, resol, data_size( 2 )
      character*6 cal_type,band_no
      character*80 req_bandname
      
      real buffer( max_pixel *  max_line )
      byte v_flag( max_pixel *  max_line ), 
     &        m_buf (max_pixel *  max_line)

      integer i, j, k, m, n, scan_flag, status


      integer ni
      
      integer imet, jmet, last_imet, last_jmet
      real x, x0, dx, y, y0, dy

c ... set program name for error messaging
      character*32 routine_name
      parameter ( routine_name = 'modwvnir_get_data' )

c ... external functions

      integer modwvnir_read_flat_file
      external modwvnir_read_flat_file

c ... Save Statements
      SAVE / MODWVNIR_DATA /

c ... Initialize 1km radiance storage array with the bad value

      do k = 1, NBANDS
        do j = 1, MAX_LINE
          do i = 1, MAX_PIXEL
            RADIANCE1( i, j, k ) = BAD_VALUE
          end do
        end do

c ...   Initialize buffer arrays

        m = 1 
        do j = 1, max_line
          do i = 1, max_pixel
            buffer( m ) = bad_value
            v_flag( m ) = -1
            m = m + 1
          end do
        end do

c ...   Get single band of radiance data
c ...   Whe are in the NIR, so get data in terms of reflectance

        resol = 1
        cal_type = 'ref'
 
        write(req_bandname,*) bands(k)

        status = modwvnir_read_flat_file(scan,resol,req_bandname,cal_type,
     &                             l1b_1km_lun,datatype_1km,
     &                             interleave_1km,resolution_1km,
     &                             error_1km,offset_1km,samples_1km,
     &                             lines_1km,bands_1km,bandnames_1km,
     &                             bandunits_1km,max_pixel,max_line,
     &                             buffer,v_flag,data_size(1),
     &                             data_size(2))
        if (status .lt. 0) then
          if ( debug .gt. 0 ) then
            write(*, '(''Error reading L1B radiance'')' )
            write(*, '(5x,''read status = '',I9)') status
          endif
          write(band_no,'(I6)' ) bands(k)
          call message( routine_name,
     &      'Failed to extract 1km L1B data for band ' // band_no
     &      ,0, 2 )
        endif

c ...   Copy good data into 1km radiance storage array

        n = 1
        do j = 1, data_size( 2 )
          do i = 1, data_size( 1 )
            radiance1( i, j, k ) = buffer( n )
            n = n + 1
          end do
        end do

      end do

c ... Get geolocation data from file


c ... Extract viewing zenith information
 
      m = 1
      do j = 1, max_line
        do i = 1, max_pixel
          view( i, j ) = bad_value
          buffer( m ) = bad_value
          v_flag( m ) = -1
          m = m + 1
        end do
      end do
 
      resol = 1
      req_bandname = 'SensorZenith'
      cal_type = 'deg'

      status = modwvnir_read_flat_file(scan,resol,req_bandname,cal_type,
     &                           geo_1km_lun,datatype_geo,
     &                           interleave_geo,resolution_geo,
     &                           error_geo,offset_geo,samples_geo,
     &                           lines_geo,bands_geo,bandnames_geo,
     &                           bandunits_geo,max_pixel,max_line,
     &                           buffer,v_flag,data_size(1),
     &                           data_size(2))
      if (status .lt. 0) then
        if ( debug .gt. 0 ) then
          write(*, '(''Error reading Viewing Zenith. '')' )
          write(*, '(5x,''read status = '',I9)') status
        endif
        call message( routine_name,
     &    'Failed to extract Viewing Zenith information ' 
     &     ,0, 2 )
       endif

c ... Fill array
      n = 1
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if (abs(buffer(n)) .lt. 67.0 .and. abs(buffer(n)).gt. 0.0)
     &      view( i, j ) = buffer( n )
            n = n + 1
        end do
      end do


c ... Extract solar zenith information
 
      m = 1
      do j = 1, max_line
        do i = 1, max_pixel
          solz( i, j ) = bad_value
          buffer( m ) = bad_value
          m = m + 1
        end do
      end do
 
      resol = 1
      req_bandname = 'SolarZenith'
      cal_type = 'deg'

      status = modwvnir_read_flat_file(scan,resol,req_bandname,cal_type,
     &                           geo_1km_lun,datatype_geo,
     &                           interleave_geo,resolution_geo,
     &                           error_geo,offset_geo,samples_geo,
     &                           lines_geo,bands_geo,bandnames_geo,
     &                           bandunits_geo,max_pixel,max_line,
     &                           buffer,v_flag,data_size(1),
     &                           data_size(2))
      if (status .lt. 0) then
        if ( debug .gt. 0 ) then
          write(*, '(''Error reading Solar Zenith. '')' )
          write(*, '(5x,''read status = '',I9)') status
        endif
        call message( routine_name,
     &    'Failed to extract Solar Zenith information ' 
     &     ,0, 2 )
       endif

c ... Fill array
      n = 1
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
            solz( i, j ) = buffer( n )
            n = n + 1
        end do
      end do

c ... Extract sensor azimuth information
 
      m = 1
      do j = 1, max_line
        do i = 1, max_pixel
          sena( i, j ) = bad_value
          buffer( m ) = bad_value
          m = m + 1
        end do
      end do
 
      resol = 1
      req_bandname = 'SensorAzimuth'
      cal_type = 'deg'

      status = modwvnir_read_flat_file(scan,resol,req_bandname,cal_type,
     &                           geo_1km_lun,datatype_geo,
     &                           interleave_geo,resolution_geo,
     &                           error_geo,offset_geo,samples_geo,
     &                           lines_geo,bands_geo,bandnames_geo,
     &                           bandunits_geo,max_pixel,max_line,
     &                           buffer,v_flag,data_size(1),
     &                           data_size(2))
      if (status .lt. 0) then
        if ( debug .gt. 0 ) then
          write(*, '(''Error reading Sensor Azimuth. '')' )
          write(*, '(5x,''read status = '',I9)') status
        endif
        call message( routine_name,
     &    'Failed to extract Sensor Azimuth information ' 
     &     ,0, 2 )
       endif

c ... Fill array
      n = 1
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
            sena( i, j ) = buffer( n )
            n = n + 1
        end do
      end do

c ... Extract solar azimuth information
 
      m = 1
      do j = 1, max_line
        do i = 1, max_pixel
          sola( i, j ) = bad_value
          buffer( m ) = bad_value
          m = m + 1
        end do
      end do
 
      resol = 1
      req_bandname = 'SolarAzimuth'
      cal_type = 'deg'

      status = modwvnir_read_flat_file(scan,resol,req_bandname,cal_type,
     &                           geo_1km_lun,datatype_geo,
     &                           interleave_geo,resolution_geo,
     &                           error_geo,offset_geo,samples_geo,
     &                           lines_geo,bands_geo,bandnames_geo,
     &                           bandunits_geo,max_pixel,max_line,
     &                           buffer,v_flag,data_size(1),
     &                           data_size(2))
      if (status .lt. 0) then
        if ( debug .gt. 0 ) then
          write(*, '(''Error reading Solar Azimuth. '')' )
          write(*, '(5x,''read status = '',I9)') status
        endif
        call message( routine_name,
     &    'Failed to extract Solar Azimuth information ' 
     &     ,0, 2 )
       endif

c ... Fill array
      n = 1
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
            sola( i, j ) = buffer( n )
            n = n + 1
        end do
      end do

c ... Extract latitude information
 
      m = 1
      do j = 1, max_line
        do i = 1, max_pixel
          lat1( i, j ) = bad_value
          buffer( m ) = bad_value
          m = m + 1
        end do
      end do
 
      resol = 1
      req_bandname = 'Latitude'
      cal_type = 'deg'

      status = modwvnir_read_flat_file(scan,resol,req_bandname,cal_type,
     &                           geo_1km_lun,datatype_geo,
     &                           interleave_geo,resolution_geo,
     &                           error_geo,offset_geo,samples_geo,
     &                           lines_geo,bands_geo,bandnames_geo,
     &                           bandunits_geo,max_pixel,max_line,
     &                           buffer,v_flag,data_size(1),
     &                           data_size(2))
      if (status .lt. 0) then
        if ( debug .gt. 0 ) then
          write(*, '(''Error reading Latitude. '')' )
          write(*, '(5x,''read status = '',I9)') status
        endif
        call message( routine_name,
     &    'Failed to extract Latitude information '
     &     ,0, 2 )
      endif

c ... Fill array
      n = 1
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( n ) .ge. -90.0 .and.
     &      buffer( n ) .le. 90.0 ) lat1( i, j ) = buffer( n )
            n = n + 1
        end do
      end do


c ... Extract longitude information

      m = 1
      do j = 1, max_line
        do i = 1, max_pixel
          lon1( i, j ) = bad_value
          buffer( m ) = bad_value
          m = m + 1
        end do
      end do

      resol = 1
      req_bandname = 'Longitude'
      cal_type = 'deg'

      status = modwvnir_read_flat_file(scan,resol,req_bandname,cal_type,
     &                           geo_1km_lun,datatype_geo,
     &                           interleave_geo,resolution_geo,
     &                           error_geo,offset_geo,samples_geo,
     &                           lines_geo,bands_geo,bandnames_geo,
     &                           bandunits_geo,max_pixel,max_line,
     &                           buffer,v_flag,data_size(1),
     &                           data_size(2))
      if (status .lt. 0) then
        if ( debug .gt. 0 ) then
          write(*, '(''Error reading Longitude. '')' )
          write(*, '(5x,''read status = '',I9)') status
        endif
        call message( routine_name,
     &    'Failed to extract Longitude information '
     &     ,0, 2 )
       endif

c ... Fill array
      n = 1
      do j = 1, data_size( 2 )
        do i = 1, data_size( 1 )
          if ( buffer( n ) .ge. -180.0 .and.
     &      buffer( n ) .le. 180.0 ) lon1( i, j ) = buffer( n )
            n = n + 1
        end do
      end do

c ...............................................................
c ... Write debug information ...................................

      if ( debug .gt. 1 ) then

        write( *, '(/,a,/)' ) 'modwvnir_get_data debug'

        i = 676
        j = 6

        write( *, '( ''Data for Nadir Pixel at '', ' //
     &    ' ''Scan'',i4,'', Pixel'',i4,'', Line'',i2)' ) scan, i, j

        write( *, '(''IR Bands: '',20i7)' )
     &    ( bands( k ), k = 1, nbands )

        write( *, '(''Radiance: '',20f7.3)' )
     &    ( radiance1( i, j, k ), k = 1, nbands )
        do j = 1, MAX_LINE
          write( *, '(i4, 20f7.3)' )
     &    j, ( radiance1( i, j, k ), k = 1, nbands )
        end do
        j=6
        write( *, '(''Auxillary Data: '', ' //
     &    ' ''Solar Zenith ='',f7.2,'' Sensor Aimuth = '', f7.2,'' Day/Night Flag  ='',i10 )')
     &    solz( i, j ), sena(i, j), -1

        write( *,
     &    ' (''Lat ='',f9.3,'', Lon ='',f10.3)')
     &    lat1( i, j ), lon1( i, j )

c ...............................................................

      endif

      end
