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
      subroutine db_mod28_initialize_output()
      
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Initialize output product arrays for the MODIS Product MOD28
c    SST product.  This version is for the direct braodcast product.
c
c!Input Parameters:
c
c!Output Parameters:
c
c    Processing Variables:
c
c
c    The following arrays in COMMON /MOD28_DATA/ are initialized:
c    raw_radiance               Radiance
c    brightess_temp             Brightness temperature 
c    product_qa                 Output product qa array
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      save
          
      include 'db_mod28uw_data.inc'
      
c ... arguments 

c ... array arguments

c ... local variables

      integer i, j, k

c ... Brightness temperature
      do j = 1, max_line
        do i = 1, max_pixel
          sst1(i,j) = bad_value
          sst4(i,j) = bad_value
          do k = 1, n_output
            raw_radiance( i, j, k ) = bad_value
            brightness_temp( i, j, k ) = bad_value
          end do
        end do
      end do

      return
      end
