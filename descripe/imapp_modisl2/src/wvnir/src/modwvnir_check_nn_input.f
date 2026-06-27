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
      real function modwvnir_check_nn_input()

      implicit none
      save
      
      include 'modwvnir_data.inc'
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    This program compares the actual input data for the neural network
c    with the appropriate minima and maxima of the training data. 
c
c!Input Parameters:
c    The following array in COMMON /MODWVNIR_DATA/ is used:
c    NN_INPUT       Input vector for the neural network
c                   (logarithms of radiance ratios, observation and 
c                   illumination geometry)
c    The following arrays in COMMON /MODWVNIR_NN_DATA/ are used:
c    MIN_TRAIN_IN   minima of input training data
c    MAX_TRAIN_IN   maxima of input training data
c
c!Return Parameters:
c    The routine returns zero if all values are within the allowed range;
c    if at least one value is too small or too large, respectively, 
c    a negative value is returned.
c    This value contains bitwise information on which input dimension(s)
c    are affected
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c-----------------------------------------------------------------------
      
c ... arguments

      real  i, value
      
      value = 0.0
      do i = 1, DIMI
c        if( nn_input(i).lt.MIN_TRAIN_IN(i) ) value = value + 2.0**i
        if( nn_input(i).lt.MIN_TRAIN_IN(i) ) value = 1.0
      end do
      do i = 1, DIMI
c        if( nn_input(i).gt.MAX_TRAIN_IN(i) ) value = value + 2.0**(i+DIMI)
        if( nn_input(i).gt.MAX_TRAIN_IN(i) ) value = 1.0
      end do
      
      modwvnir_check_nn_input = -1.0 * value
      return
      end

