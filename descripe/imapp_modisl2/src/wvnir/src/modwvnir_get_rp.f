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
      subroutine modwvnir_get_rp(  
     &                           cfgname,
     &                           rp_debug,
     &                           l1b_1km_hdr,
     &                           l1b_geo_hdr,
     &                           beg_lin,
     &                           nlins,
     &                           beg_ele,
     &                           neles
     &                           )

      implicit none
      save

C!F77 ************************************************************
C ...
C!Description:
C     Routine which gets the runtime debug parameter 
C      and data header file names for the MODWVNIR direct
C      broadcast code.
C      debug values: 0 = no debug writes
C                    1 = debug information
C
C!Input parameters:
C cfgname      Configuration file name
C
C!Output Parameters:
C rp_debug     Debug value - 0 or 1
C l1b_1km_hdr  Header file name for 1km data
C l1b_geo_hdr  Header file name for geolocation data
C beg_lin      Beginning line number to process
C nlins        Number of lines to process
C beg_ele      Beginning element number to process
C neles        Number of elements to process
C
C!Revision History:
C
C!Team-Unique Header:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C!References and Credits:
C See Cloud Mask ATBD-MOD-25.
C          external subroutines:  message.f
C
C!Design Notes:
C
C!END****************************************************************

      include 'modwvnir_data.inc'

C --- constants
      character*80 routine_name
      parameter (routine_name =  'modwvnir_get_rp')

C --- parameters
      character*(*) l1b_1km_hdr, l1b_geo_hdr, cfgname
      integer rp_debug, beg_lin, nlins, beg_ele, neles

C --- external functions
      integer hdrgetkeyint, hdrgetkeystr

C --- internal variables
      character*80  keyname
      character*(PATH_MAX) filename
      integer keyindex, status, len_name
      logical remove_all
      character*72 text

C --- initialize function status
      remove_all = .FALSE.

C --- initialize filename
      filename = ' '

C --- get debug value
      keyname = 'debug'
      keyindex = 1
      status = hdrgetkeyint( cfgname, keyname, keyindex, rp_debug )
      if( status.lt.0 ) rp_debug = 0

C --- get the name of the L1B geo header file
      keyname = 'GEO_1KM_HDR'
      keyindex = 1
      status = hdrgetkeystr( cfgname, keyname, keyindex, filename )
      if( status.lt.0 ) then
         call message( routine_name,
     &                 'FAILED: No L1B Geo header file name (' // keyname // ')',
     &                 0, 2 )
         return
      else
         call strcompress( filename, remove_all, len_name )
         l1b_geo_hdr = filename(1:len_name)
      endif

C --- get the name of the L1B 1km header file
      keyname = 'L1B_1KM_HDR'
      keyindex = 1
      status = hdrgetkeystr( cfgname, keyname, keyindex, filename )
      if( status.lt.0 ) then
         call message( routine_name,
     &                 'FAILED: No L1B 1KM header file name (' // keyname // ')',
     &                 0, 2 )
         return
      else
         call strcompress( filename, remove_all, len_name )
         l1b_1km_hdr = filename(1:len_name)
      endif

C --- Retrieve beginning line number
      keyname = 'bline'
      keyindex = 1
      status=hdrgetkeyint(filename(1:len_name),keyname,keyindex,beg_lin)
      if( status.lt.0 ) beg_lin = 0
      if( beg_lin.lt.0 ) then

         call message( routine_name,
     &                 'FAILED: Invalid beginning 1km line number '
     &                 // keyname // ' @ ' // filename(1:len_name) // ')', 
     &                 beg_lin, 2 )
         return
      endif

      write(text, '(''INFO: Beginning 1km line number = '',i4)') beg_lin
      call message( routine_name, text , 0, rp_debug )

C --- Retrieve number of lines to process
      keyname = 'lines'
      keyindex = 1
      status=hdrgetkeyint(filename(1:len_name),keyname,keyindex,nlins)
      if( status.lt.0 ) then
        call message( routine_name,
     &                 'FAILED: Invalid number of 1km lines (' 
     &                 // keyname // ' @ ' // filename(1:len_name) // ')', 
     &                 nlins, 2 )
          return
      endif
      write(text, '(''INFO: Number of 1km lines = '',i4)') nlins
      call message( routine_name, text , 0, rp_debug )

C --- Retrieve beginning element number
      keyname = 'belem'
      keyindex = 1
      status=hdrgetkeyint(filename(1:len_name),keyname,keyindex,beg_ele)
      if( status.lt.0 ) beg_ele = 0
      if( beg_ele.lt.0 .or. beg_ele.gt.1354) then
         call message( routine_name,
     &                 'FAILED: Invalid beginning 1km element '
     &                 // keyname // ' @ ' // filename(1:len_name) // ')', 
     &                 beg_ele, 2 )
         return
      endif
      write(text, '(''INFO: Beginning 1km element = '',i4)') beg_ele
      call message( routine_name, text , 0, rp_debug )

C --- Retrieve number of elements to process
      keyname = 'samples'
      keyindex = 1
      status=hdrgetkeyint(filename(1:len_name),keyname,keyindex,neles)
      if( status.lt.0 ) neles = -1
      if( neles.lt.0 .or. neles.gt.MAX_PIXEL) then
         call message( routine_name,
     &                 'FAILED: Invalid number of samples/line '
     &                 // keyname // ' @ ' // filename(1:len_name) // ')', 
     &                  neles, 2 )
         return
      endif
      write(text, '(''INFO: Number of samples/line = '',i4)') neles
      call message( routine_name, text , 0, rp_debug )

      return
      end
