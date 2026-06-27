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
      subroutine modwvnir_file_close(l1b_1km_lun,geo_1km_lun,
     &                         modwvnir_lun,hdr_lun,output_type,hdf_lun)

      implicit none
      save

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
c Program which will close files used in the main DB WVNIR program. 
C
C!INPUT PARAMETERS:    None.
C
C!OUTPUT PARAMETERS:
C
C   l1b_1km_lun        LUN for open 1km L1B data file
C   geo_1km_lun        LUN for open geo file
C   modwvnir_lun       LUN for MODWVNIR output binary file
C   hdr_lun            LUN for MODWVNIR output binary header file
C   output_type        Flag defining type of output: 
C                         1 = binary only
C                         2 = hdf only 
C                         3 = binary and hdf
C   hdf_lun            LUN for MODWVNIR output HDF file
C
C!REVISION HISTORY:
c
C!TEAM-UNIQUE HEADER:
C
C!REFERENCES AND CREDITS
C
C EXTERNALS:
C
C       NAMED CONSTANTS:
C
C
C!END-------------------------------------------------------------------

      include 'hdf.inc'

c --- arguments
      integer l1b_1km_lun, geo_1km_lun,output_type,hdf_lun, 
     &        modwvnir_lun,hdr_lun

c ... local variables
      integer rtn

c ... HDF external routines
      integer sfend
      external sfend

c ... close files
      close(l1b_1km_lun)
      close(geo_1km_lun)
      if (output_type .eq. 1 .or. output_type .eq. 3) then
        close(modwvnir_lun)
        close(hdr_lun)
      endif

      if (output_type .eq. 2 .or. output_type .eq. 3) then
        rtn = sfend(hdf_lun)
      endif

      return
      end
