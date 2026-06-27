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
      real function sigmoide(x, temperature)

      implicit none
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Compute the sigmoide for the activation of neural network neurons
c
c!Input Parameters:
c    x                input value
c    temperature      temperature
c
c!Return Parameters:
c    calculated sigmoide
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c-----------------------------------------------------------------------
      real      x, temperature
      
      sigmoide = temperature * x * (-1.)
      if( sigmoide .lt. -80. ) sigmoide = -80.
      if( sigmoide .gt.  80. ) sigmoide =  80.
      sigmoide = 1. / (1. + exp(sigmoide))
      
      return
      end
