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
        integer function siflip(ispin)   

c ...   Routine which will flip bytes for a 
c ...   integer*4 variable (ispin)
c ....  version of 10.01.02

        implicit none
        character*1 ch   
        character*4 carg   
        integer ispin,ispout
        equivalence (carg,ispout)

        save

        ispout = ispin
        ch = carg(1:1)
        carg(1:1) = carg(4:4) 
        carg(4:4) = ch
        ch = carg(2:2)
        carg(2:2) = carg(3:3) 
        carg(3:3) = ch
        siflip = ispout

        return
        end
