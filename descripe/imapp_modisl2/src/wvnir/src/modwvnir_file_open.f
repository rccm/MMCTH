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
      subroutine modwvnir_file_open( 
     &                                cfgname,
     &                                l1b_1km_neles,
     &                                output_type,
     &                                l1b_1km_lun,
     &                                geo_1km_lun,
     &                                modwvnir_lun,
     &                                hdr_lun,
     &                                hdf_lun,
     &                                debug
     &                              )
      implicit none
      save

      include 'modwvnir_file_IO.inc'
      include 'modwvnir_data.inc'
      include 'hdf.inc'

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
c Program which will open files and return unit numbers for 
c all files used in the main DB WVNIR program. 
C
C!INPUT PARAMETERS>
C   cfgname            Configuration file name
C   output_type        Flag defining type of output: 
C                         1 = binary only
C                         2 = hdf only 
C                         3 = binary and hdf
C   debug              0: no debug messages / 1: debug messages
C
C!OUTPUT PARAMETERS:
C
C   l1b_1km_neles      Number of elements per scan line of data
C   l1b_1km_lun        LUN for open 1km L1B data file
C   geo_1km_lun        LUN for open geo file
C   modwvnir_lun       LUN for MODWVNIR output file
C   hdr_lun            LUN for MODWVNIR output header file
C   hdf_lun            LUN for MOD28 output HDF file
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
C       FUNCTIONS AND SUBROUTINES:
C               hdrgetkeystr.f
c               message.f
c               strcompress.f
C
C!END-------------------------------------------------------------------

c --- constants
      character*18 routine_name
      parameter (routine_name = 'modwvnir_file_open')

c --- arguments
      integer l1b_1km_neles, l1b_1km_lun, geo_1km_lun, output_type,
     &        modwvnir_lun, hdr_lun, hdf_lun
      character*(*) cfgname

c --- external functions
      integer hdrgetkeystr, sfstart

c --- internal variables
      character*12 req
      character*80 keyname
      character*(PATH_MAX) file
      integer keyindex, level, status, io_err, len_file
      logical remove_all

c --- initialize the lun
      l1b_1km_lun  = -1
      geo_1km_lun  = -1
      modwvnir_lun = -1
      hdr_lun	= -1
      hdf_lun	= -1


c --- initialize dummy variables
      req = ' '
      file = ' '

C 1KM Level 1B MODIS DATA FILES

c --- get the file requirement
      keyname = 'L1B_1KM_REQ'
      keyindex = 1
      status = hdrgetkeystr( cfgname, keyname, keyindex, req )
      if( status.lt.0 ) req = 'OPTIONAL'
      level = 3
      if( req(1:3).eq.'MAN' ) level = 2

c --- get file access information 
      keyname = 'L1B_1KM_NAME'
      keyindex = 1
      status = hdrgetkeystr( cfgname, keyname, keyindex, file )

c --- check for error reading info file
      if( status.lt.0 ) then 
         call message( 
     &                 routine_name,
     &                'FAILED: Read of 1km file info (' // keyname // ')',
     &                 0, level 
     &               )

      else
                   
c ------ issue the OPEN request
         OPEN( 
     &         FILE=file,
     &         UNIT=L1B_1KM_UNIT,
     &         FORM=L1B_1KM_FORM,
     &         ACCESS=L1B_1KM_ACCESS,
     &         STATUS=L1B_1KM_STATUS,
     &         RECL=l1b_1km_neles*4,
     &         IOSTAT=io_err, 
     &         ERR=9005
     &       ) 

         call message(
     &                 routine_name,
     &                'INFO: Using 1km L1B file: '// file ,
     &                 0, debug
     &               )

         l1b_1km_lun = L1B_1KM_UNIT
         goto 140

c ------ Failed to open 
9005     continue 
         call message( 
     &                 routine_name,
     &                'FAILED: Could not open 1km L1B file:  '// file , 
     &                 0, level 
     &               )

      endif
140   continue
              

C GEOLOCATION FILES: LATITUDE, LONGITUDE, HEIGHT, ...

c --- get the file requirement
      keyname = 'GEO_1KM_REQ'
      keyindex = 1
      status = hdrgetkeystr( cfgname, keyname, keyindex, req )
      if( status.lt.0 ) req = 'OPTIONAL'
      level = 3
      if( req(1:3).eq.'MAN' ) level = 2

c --- get file access information 
      keyname = 'GEO_1KM_NAME'
      keyindex = 1
      status = hdrgetkeystr( cfgname, keyname, keyindex, file )

c --- check for error reading info file
      if( status.lt.0 ) then 
         call message( 
     &                 routine_name,
     &                'FAILED: Read of Geo file info (' // keyname // ')',
     &                 0, level 
     &               )

      else
                   
c ------ issue the OPEN request
         OPEN( 
     &         FILE=file,
     &         UNIT=GEO_1KM_UNIT,
     &         STATUS=GEO_1KM_STATUS,
     &         FORM=GEO_1KM_FORM,
     &         ACCESS=GEO_1KM_ACCESS,
     &         RECL=l1b_1km_neles*4,
     &         IOSTAT=io_err, 
     &         ERR=9006
     &       ) 

         call message(
     &                 routine_name,
     &                'INFO: Using 1km GEO file: '// file ,
     &                 0, debug
     &               )

         geo_1km_lun = GEO_1KM_UNIT
         goto 150

c ------ Failed to open 
9006     continue 
         call message( 
     &                 routine_name,
     &                'FAILED: Could not open GEO file: '// file ,
     &                 0, level 
     &               )

      endif
150   continue

C OUTPUT FILES: WVNIR
C ***  NEW RUNTIME ARGUMENT DECIDES WHAT KIND OF OUTPUT FILE TO CREATE
C ********************************************************************
c ... Binary File first

      IF (output_type .eq. 1 .or. output_type .eq. 3) then

c ---   get the file requirement
        keyname = 'MODWVNIR_REQ'
        keyindex = 1
        status = hdrgetkeystr( cfgname, keyname, keyindex, req )
        if( status.lt.0 ) req = 'OPTIONAL'
        level = 3
        if( req(1:3).eq.'MAN' ) level = 2

c ---   get file access information 
        keyname = 'MODWVNIR_NAME'
        keyindex = 1
        status = hdrgetkeystr(cfgname, keyname, keyindex,file )
        call strcompress( file, remove_all, len_file ) 

c ---   check for error reading info file
        if( status.lt.0 ) then 
           call message( 
     &                   routine_name,
     &                  'FAILED: Read of WVNIR output file info (' // keyname // ')',
     &                   0, level 
     &                 )

        else
                   
c ------   issue the OPEN request
           OPEN( 
     &           FILE=file(1:len_file),
     &           UNIT=MODWVNIR_UNIT,
     &           STATUS=MODWVNIR_STATUS,
     &           ACCESS=MODWVNIR_ACCESS,
     &           RECL=(MAX_PIXEL) * 4,
     &           IOSTAT=io_err, 
     &           ERR=9017
     &         )

           call message(
     &                   routine_name,
     &                  'INFO: Creating WVNIR Output file: '// file ,
     &                   0, debug
     &                 )

           modwvnir_lun = MODWVNIR_UNIT
           goto 260

c ------   Failed to open 
9017       continue 
           call message( 
     &                   routine_name,
     &                  'FAILED: Could not open output file: ' // file,
     &                   0, level 
     &                 )

        endif
260     continue

C OUTPUT FILES: WVNIR HEADER FILE

c ---   get the file requirement
        keyname = 'MODWVNIR_REQ'
        keyindex = 1
        status = hdrgetkeystr( cfgname, keyname, keyindex, req )
        if( status.lt.0 ) req = 'OPTIONAL'
        level = 3
        if( req(1:3).eq.'MAN' ) level = 2

c ---   get file access information 
        keyname = 'MODWVNIR_HDR'
        keyindex = 1
        status = hdrgetkeystr( cfgname, keyname, keyindex,file )
        call strcompress( file, remove_all, len_file ) 


c ---   check for error reading info file
        if( status.lt.0 ) then 
           call message( 
     &                   routine_name,
     &                  'FAILED: Read of WVNIR output header file info (' // keyname // ')',
     &                   0, level 
     &                 )

        else
                   
c ------   issue the OPEN request
           OPEN( 
     &           FILE=file(1:len_file),
     &           UNIT=MODWVNIR_HDR_UNIT,
     &           STATUS=MODWVNIR_HDR_STATUS,
     &           ACCESS=MODWVNIR_HDR_ACCESS,
     &           IOSTAT=io_err, 
     &           ERR=9018
     &         ) 

           call message(
     &                   routine_name,
     &                  'INFO: Creating WVNIR Output header file: '// file ,
     &                   0, debug
     &                 )
           hdr_lun = MODWVNIR_HDR_UNIT
           goto 261

c ------   Failed to open 
9018       continue 
           call message( 
     &                   routine_name,
     &                  'FAILED: Could not open WVNIR output header file:' // file,
     &                   0, level 
     &                 )

261       continue
          endif
C **** END OF BINARY FILE CREATION IF BLOCK
C **************************************************************
        endif

C ***  NEW RUNTIME ARGUMENT DECIDES WHAT KIND OF OUTPUT FILE TO CREATE
C ********************************************************************
C ---- HDF output file ---------------------------------------
      if (output_type .eq. 2 .or. output_type .eq. 3) then

c ---   get the file requirement
        keyname = 'MODWVNIR_HDF_REQ'
        keyindex = 1
        status = hdrgetkeystr( cfgname, keyname, keyindex, req )
        if( status.lt.0 ) level=2

c ---   get file access information
        keyname = 'MODWVNIR_HDF_NAME'
        keyindex = 1
        status = hdrgetkeystr( cfgname, keyname, keyindex,file )

c ---   check for error reading info file
        if( status.lt.0 ) then
           call message(
     &                   routine_name,
     &                  'FAILED: Opening of output WVNIR HDF file. ',
     &                   0, level
     &                 )

        else

c ------   issue the OPEN request - this is an HDF file)
c ...      Open the input HDF file for read/write) 
c ...      (DFACC_READ is defined in hdf.inc)
           hdf_lun = sfstart(file, DFACC_CREATE)
           if (hdf_lun .eq. -1) then
              call message(
     &                     routine_name,
     &                    'FAILED: Could not open output hdf WVNIR file. ',
     &                     0, level
     &                     )
           endif
           call message(
     &                   routine_name,
     &                  'INFO: Creating WVNIR HDF Output file: '// file ,
     &                   0, 4
     &                  )
        endif
C **** END OF HDF OUTPUT FILE CREATION IF BLOCK
C **************************************************************
      endif

      call message(routine_name, 'INFO: finished', 0, 0)
      return
      end
