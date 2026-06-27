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
      real function read_real(lun_sst,iword,iret)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C   Read one word (4 bytes) out of the Reynolds blended SST file
C   and flip the bytes if it is a little endian machine.  This
C   routine works for words of type REAL
C
C !INPUT PARAMETERS:
C  lun_sst        Reynolds blended sst binary file file handle
C  iword          Record number (word) to read out of the file
C
C !OUTPUT PARAMETERS:
C  iret           Return value. Success = 0, Fail = -1
C  read_real      Output real sst value read from file
C
C !REVISION HISTORY:
C
C!Team-Unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

c ...   arguments
        integer lun_sst, iword, iret

c ...   local variables
        integer ios
        real xval

c ...   external functions
        real sbflip
        logical big_endian

        read(lun_sst,rec=iword,iostat=ios) xval
        if( ios.eq.0 ) then
           if(.not.big_endian()) xval = sbflip(xval)
           read_real = xval
        else
           iret = -1
           read_real = -9999.0
        endif

        return
        end

