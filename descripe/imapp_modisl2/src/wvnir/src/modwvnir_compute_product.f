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
      subroutine modwvnir_compute_product(line,pixel)

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Compute NIR water vapor products.
c
c!Input Parameters:
c    LINE             Line number within a swath (1-10)
c    PIXEL            Pixel number within a 1km scan (1-1500)
c
c    The following array in COMMON /MODWVNIR_DATA/ is used:
c    NN_INPUT       Input vector for the neural network
c                   (logarithms of radiance ratios, observation and 
c                   illumination geometry)
c
c    The following arrays in COMMON /MODWVNIR_NN_DATA/ are used:
c    OFFSET1, OFFSET2, SCALE    Scaling coefficients for input and output data
c    INV, OUTV                  Internal network matrices
c
c!Output Parameters:
c    The following array in COMMON /MODWVNIR_DATA/ is filled:
c    WVNIR_OUT      NIR water vapor output product
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c-----------------------------------------------------------------------
      
      implicit none
      save
      
      include 'modwvnir_data.inc'
      
c ... arguments

      integer line, pixel, i, j

c ... local variables

      real zw( DIMH + 1 )
      real sig_in
      real  bad_input

c ... external functions

      real  sigmoide, modwvnir_check_nn_input

c ... Check whether input is within range of training data
      bad_input = modwvnir_check_nn_input()

c      write(*, '("Input flag: ", 2i5,": ",  f14.1)'), pixel, line, bad_input

      if ( bad_input .gt. -0.5 ) then
            
c ... Linear scaling of input
c      write(*, '(8f12.4)'), nn_input
        
	do i = 1, DIMI
	  nn_input(i) = OFFSET1(i) + SCALE(i) * (nn_input(i) - OFFSET2(i))
        end do
  
        nn_input(DIMI + 1) = BIAS_I
   
        
        do i = 1, DIMH
          sig_in = 0.0
          do j = 1, DIMI + 1
            sig_in = sig_in + nn_input(j) * INV(j, i)
          end do
          zw(i) = sigmoide(sig_in,TEMPERATURE / DIMI )
        end do
        zw(DIMH + 1) = BIAS_H
        
        sig_in = 0.0
        do j = 1, DIMH + 1
          sig_in = sig_in + zw(j) * OUTV(j, 1)
        end do
        wvnir_out(pixel, line) = sigmoide(sig_in,TEMPERATURE / DIMH )
        wvnir_out(pixel, line) =  
     &          (wvnir_out(pixel, line) - OFFSET1(DIMI + 1)) / 
     &          SCALE(DIMI + 1) + OFFSET2(DIMI + 1)

c ... g/cm**2 -> kg/m**2
        wvnir_out(pixel, line) =  wvnir_out(pixel, line) * 10.

      else
        wvnir_out(pixel, line) = bad_input
c      write(*, '("Input data outside training data: ", 2i5,": ",  f14.1)'), pixel, line, bad_input
      endif
      return
      end
