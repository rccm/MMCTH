C*----------------------------------------------------------------------
C    Copyright (C) 2011,  Space Science and Engineering Center,
C    University C  of Wisconsin-Madison, Madison WI.
C
C    This program is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program.  If not, see <http://www.gnu.org/licenses/>.
C------------------------------------------------------------------------*/
C
C
      INTEGER FUNCTION STRING_LENGTH(STRING)

c ... Returns the index of the last non-blank character in a string
c ... (Returns zero if the string contains all blanks)

      implicit none

      character*(*) string
      integer index

      string_length = 0
      do index = 1, len(string)
        if (string(index:index) .ne. ' ') string_length = index
      end do
      
      END
