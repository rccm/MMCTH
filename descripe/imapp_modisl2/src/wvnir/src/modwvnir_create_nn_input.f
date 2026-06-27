C--------------------------------------------------------------------
C  Copyright (C) 2002,  Institut fuer Weltraumwissenschaften
C  Freie Universitaet Berlin, Berlin
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

      subroutine modwvnir_create_nn_input(line,pixel)

C-----------------------------------------------------------------------
C
C!F77
C
C!Description:
C    This program creates the input for the neural
C    network used  for the Direct Broadcast version of the NIR water vapor
C    product.
C
C!Input parameters:
c    LINE             Line number within a swath (1-10)
c    PIXEL            Pixel number within a 1km scan (1-1500)
C
c    The following arrays in COMMON /MODWVNIR_DATA/ are used:
c    RADIANCE1      Radiances for NIR bands
c    VIEW           Viewing zenith angle value
c    SOLZ           Solar zenith angle value
c    SENA           Viewing azimuth angle value
c    SOLZ           Solar azimuth angle value
C
C!Output Parameters:
c    The following arrays in COMMON /MODWVNIR_DATA/ are filled:
c    NN_INPUT       Input vector for the neural network
c                   (logarithms of radiance ratios, observation and 
c                   illumination geometry)
C!End
C
C-----------------------------------------------------------------------

      
      implicit none
      save
      
      include 'modwvnir_data.inc'
      
c ... arguments

      integer line, pixel

c ... local variables
      real d2r
      parameter(d2r = 0.0174533 )     
      real rel_azi
      
      
      rel_azi = abs(SENA(pixel, line) - SOLA(pixel, line))
      if( rel_azi.gt. 180 ) rel_azi = 360. - rel_azi
      rel_azi = abs(180. - rel_azi)
             
      nn_input(1) = alog(RADIANCE1(pixel, line, 2) / RADIANCE1(pixel, line, 1))
      nn_input(2) = alog(RADIANCE1(pixel, line, 3) / RADIANCE1(pixel, line, 1))
      nn_input(3) = alog(RADIANCE1(pixel, line, 4) / RADIANCE1(pixel, line, 1))
      nn_input(4) = alog(RADIANCE1(pixel, line, 3) / RADIANCE1(pixel, line, 2))
      nn_input(5) = cos(SOLZ(pixel,line)*d2r)
      nn_input(6) = cos(VIEW(pixel,line)*d2r)

c     multiply the relative azimuth with the sine of the viewing zenith in order
c     to get rid of the azimuth jump around nadir
      nn_input(7) = sin(VIEW(pixel,line)*d2r) * rel_azi
      return
      end
