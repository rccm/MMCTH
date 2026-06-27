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
      subroutine modwvnir_chk_input(nscans,npixels,beg_lin,nlins,
     &                              beg_ele,neles,beg_scan,ibes,
     &                              debug)

      implicit none
      save

c ... Include files Needed for time conversion
      include 'modwvnir_data.inc'

c-------------------------------------------------------------------
C!F77 
C 
C!DESCRIPTION:
C     Routine for checking consistency of input parameters.  
C
c!Input Parameters:
c nscans        Number of scans for this granule
c npixels       Number of pixels in this granule
c beg_lin       Beginning line number to process
c nlins         Number of lines to process
c beg_ele       Beginning element number to process
c neles         Number of elements to process
C debug         0: no debug messages / 1: debug messages
c
c!Output Parameters:
c beg_scan          Beginning scan number to process based upon
c                    user inputs
c nscans            Updated Number of scans to process based upon
c                     user inputs
c npixels           Updated number of pixesl to process based upon
c                     user inputs
c ibes              Beginning pixel number to process
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-25.
c Externals:  error messaging routine   -  message.f
c
c!END
c-------------------------------------------------------------------

c     scalar arguments
      integer nscans,beg_lin,nlins,beg_ele,neles,beg_scan,ibes,
     &        npixels

c     local scalars
      integer out_lines_1km,out_elements_1km,num_lines,num_eles
    
      character*72 text

c ... Define processing intervals
c ... nscans will be number of scans to process
c ... ibes is beginning pixel number, npixels will be number of pixesl
c ...  to process.  
      beg_scan = (beg_lin - 1) / max_line + 1
      out_lines_1km = nscans * max_line
      out_elements_1km = npixels
      num_lines = beg_lin + nlins - 1
      num_eles = beg_ele + neles - 1
      
      if (num_lines .gt. out_lines_1km .or. num_eles .gt. npixels)
     &  call message( 'modwvnir_chk_input',
     &    'ERROR: Number of lines or elements to process was .gt. max'
     &    // char(10) // 'Check runtime parameters '
     &    // 'which define processing interval. ', 0, 2 )

      if (beg_lin .eq. 0 .and. beg_ele .eq. 0) then
          ibes = 1
          beg_scan = 1
      else
          ibes = beg_ele
          nscans = (num_lines -1) / max_line + 1
          npixels = neles + ibes - 1
c ...     Check if at least one pixel can be processed
          if (nlins .lt. 1 .or. neles .lt. 1)
     &    call message( 'modwvnir_chk_input',
     &     'ERROR: Number of lines or elements to process was lt 1'
     &     // char(10) // ' Check runtime parameters '
     &     // 'which define processing interval. ', 0, 2 )
      endif

      if(nscans .le. 0 .or. npixels .le. 0 .or. ibes .le. 0 .or.
     &   ibes .le. 0) then
        call message( 'modwvnir_chk_input',
     &  'ERROR: Number of lines or elements to process was .le. 0.'
     &    // char(10) // ' Check runtime parameters '
     &    // 'which define processing interval. ',  0, 2 )
      endif

c ---------------------------------------------------------------------
      if(beg_lin .eq. 0) then
        write(text, '(''INFO: Processing dimensions: '',2i6)') nscans, npixels
      else
        write(text, '(''INFO: Processing dimensions: '',2i6)') nscans-beg_scan+1 , neles
      endif
      call message( 'modwvnir_chk_input', text , 0, debug )
c ---------------------------------------------------------------------

      return
      end
