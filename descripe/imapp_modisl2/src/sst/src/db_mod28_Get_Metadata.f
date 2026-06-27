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
      subroutine db_mod28_Get_Metadata(header_1km,header_geo,
     &           nscans,npixels,datatype_1km,datatype_geo,
     &           interleave_1km,interleave_geo,
     &           resolution_1km,resolution_geo,offset_1km,
     &           offset_geo,samples_1km,samples_geo,
     &           lines_1km,lines_geo,error_1km,
     &           error_geo,bands_1km,bands_geo,bandnames_1km,
     &           bandnames_geo,bandunits_1km,bandunits_geo,
     &           platform_1km)
      implicit none
      save

C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C    Extract the information needed to begin processing the 
c     data for the Sea Surface Temperature (MOD28)
c     MODIS Direct Broadcast algorithm.
C
C !DESCRIPTION:
C    Extracts information from both the the Geolocation and
C    L1b header.
C
C !INPUT PARAMETERS:
C    header_1km      L1b 1km file header information.
C    header_geo      Geolocation file header information.
C  
C !OUTPUT PARAMETERS:
C    nscans          Total number of 10 line scans in the l1b data set
C    npixels         Total number of pixels in a line of data
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
C
C !REVISION HISTORY
C
C !TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C----------------------------------------------------------------------

      include           'db_mod28uw_data.inc'
      include           'db_mod28uw_debug.inc'

c ... array arguments 
      character*(*) header_1km,header_geo,
     +              interleave_1km,interleave_geo,
     +              bandnames_1km(*),bandnames_geo(*),
     +              bandunits_1km(*),
     +              bandunits_geo(*),platform_1km

c ... scalar arguments
      integer nscans,npixels,datatype_1km,datatype_geo,
     +        resolution_1km,resolution_geo,
     +        offset_1km,offset_geo,samples_1km,
     +        samples_geo,lines_1km,lines_geo,
     +        bands_1km,bands_geo

      real error_1km,error_geo

c --- external functions
      integer hdrgetkeydbl
      integer hdrgetkeystr
      integer hdrgetkeyint

c --- internal variables
      character*(PATH_MAX) filename
      character*80 keyname,function
      double precision d_value
      integer keyindex,status,len
      logical remove_all


c --- initalize strings
      function = 'db_mod6_get_metadata'

      remove_all = .TRUE.

c ... initialize variables
      d_value = 0.0d0


C ---------------- 1KM -------------------

c --- copy file name to internal string
      call strcompress( header_1km, remove_all, len )
      filename = header_1km(1:len)

c --- extract platform name from header file
      platform_1km = '  '
      keyname = 'platform'
      keyindex = 1
      status=hdrgetkeystr(filename,keyname,keyindex,
     &                    platform_1km)
      if( status.lt.0 ) then
         call message( function, 'Could not get platform name from
     &                 header file ', 0, 2)
         return
      endif

c --- datatype  
      keyname = "data type"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    datatype_1km) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - 1km datatype', 0, 3) 
         return
      endif

c --- interleave  
      keyname = "interleave"
      keyindex = 1 
      status=hdrgetkeystr(filename,keyname,keyindex,
     &                    interleave_1km) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - 1km interleave', 0, 3) 
         return
      endif

c --- resolution  
      keyname = "resolution"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    resolution_1km) 
      if( status.lt.0 ) resolution_1km = 1

c --- offset  
      keyname = "offset"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    offset_1km) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - 1km offset', 0, 3) 
         return
      endif

c --- error  
      keyname = "bad value"
      keyindex = 1 
      status=hdrgetkeydbl(filename,keyname,keyindex,
     &                    d_value) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - 1km error', 0, 3) 
         return
      endif
      error_1km = REAL(d_value)

c --- get the number of scans in the image  
      keyname = "lines"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    lines_1km) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - 1km lines', 0, 3) 
         return
      endif
      nscans = lines_1km / MAX_LINE

c --- get the maximum number of pixels per scan
      keyname = "samples"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    samples_1km) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - 1km samples', 0, 3) 
         return
      endif
      npixels = samples_1km

c --- get the number of bands
      keyname = "bands"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    bands_1km) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - 1km bands', 0, 3) 
         return
      endif

c --- get the band names and units
      do keyindex = 1, bands_1km
      keyname = "band names"
      status=hdrgetkeystr(filename,keyname,keyindex,
     &                    bandnames_1km(keyindex)) 
      if( status.lt.0 ) then
         call message(function,'FAILED - 1km band names',0,3) 
         return
      endif
      keyname = "band units"
      status=hdrgetkeystr(filename,keyname,keyindex,
     &                    bandunits_1km(keyindex)) 
      if( status.lt.0 ) then
         call message(function,'FAILED - 1km band units',0,3) 
         return
      endif
      enddo

C ---------------- GEO -------------------

c --- copy file name to internal string
      call strcompress( header_geo, remove_all, len )
      filename = header_geo(1:len)

c --- datatype  
      keyname = "data type"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    datatype_geo) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - Geo datatype', 0, 3) 
         return
      endif

c --- interleave  
      keyname = "interleave"
      keyindex = 1 
      status=hdrgetkeystr(filename,keyname,keyindex,
     &                    interleave_geo) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - Geo interleave', 0, 3) 
         return
      endif

c --- resolution  
      keyname = "resolution"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    resolution_geo) 
      if( status.lt.0 ) resolution_geo = 1

c --- offset  
      keyname = "offset"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    offset_geo) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - Geo offset', 0, 3) 
         return
      endif

c --- error  
      keyname = "bad value"
      keyindex = 1 
      status=hdrgetkeydbl(filename,keyname,keyindex,
     &                    d_value) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - Geo error', 0, 3) 
         return
      endif
      error_geo = REAL(d_value)

c --- get the number of scans in the image  
      keyname = "lines"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    lines_geo) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - Geo lines', 0, 3) 
         return
      endif

c --- get the maximum number of pixels per scan
      keyname = "samples"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    samples_geo) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - Geo samples', 0, 3) 
         return
      endif

c --- get the number of bands
      keyname = "bands"
      keyindex = 1 
      status=hdrgetkeyint(filename,keyname,keyindex,
     &                    bands_geo) 
      if( status.lt.0 ) then
         call message( function, 'FAILED - Geo bands', 0, 3) 
         return
      endif

c --- get the band names and units
      do keyindex = 1, bands_geo
        keyname = "band names"
        status=hdrgetkeystr(filename,keyname,keyindex,
     &                    bandnames_geo(keyindex)) 
        if( status.lt.0 ) then
          call message(function,'FAILED - Geo band names',status,3) 
          return
        endif
        keyname = "band units"
        status=hdrgetkeystr(filename,keyname,keyindex,
     &                    bandunits_geo(keyindex)) 
        if( status.lt.0 ) then
          bandunits_geo(keyindex) = ' '
        endif
      enddo


c ---------------------------------------------------------------------
      if (debug .gt. 0) then
        WRITE(h_output,'(10x,'' Processing dimensions '')')
        WRITE(h_output,'(10x,'' Scans in Image    Pixels per Scan'','//
     &     '/,5x,I15,5x,I15/)') nscans, npixels
      endif
c ---------------------------------------------------------------------

      return
      end
