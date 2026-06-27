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
      subroutine db_mod28_file_open( 
     &                                filename,
     &                                l1b_1km_neles,
     &                                output_type,
     &                                l1b_1km_lun,
     &                                geo_1km_lun,
     &                                mod28_lun,
     &                                hdr_lun,
     &                                hdf_lun,
     &                                debug_lun
     &                              )
      implicit none
      save

      include 'hdf.inc'
      include 'sst.inc'
      include 'db_mod28uw_data.inc'

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
c Program which will open files and return unit numbers for 
c all files used in the main UW SST retreival program.
C
C!INPUT PARAMETERS:
C   filename           Configuration file name
C   l1b_1km_neles      Number of elements per scan line of data
C   output_type        Flag defining type of output: 
C                         1 = binary only
C                         2 = hdf only 
C                         3 = binary and hdf
C
C!OUTPUT PARAMETERS:
C
C   l1b_1km_lun        LUN for open 1km L1B data file
C   geo_1km_lun        LUN for open geo file
C   mod28_lun          LUN for MOD28 output binary file
C   hdr_lun            LUN for MOD28 output binary header file
C   hdf_lun            LUN for MOD28 output HDF file
C   debug_lun          LUN for open debug messages file
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
c               db_message.f
C
C!END-------------------------------------------------------------------

c --- constants
      character*18 FUNCTION
      parameter (FUNCTION = 'db_mod28_file_open')

c --- arguments
      integer l1b_1km_neles, l1b_1km_lun, geo_1km_lun, debug_lun, 
     &mod28_lun, hdr_lun, hdf_lun, output_type
      character*(*) filename

c --- External routines
      integer read_int, hdrgetkeystr, sfstart
      real read_real
      external read_int, read_real, hdrgetkeystr, sfstart


c --- internal variables
      character*12 req
      character*80 keyname
      character*(PATH_MAX) file
      integer len_file, keyindex, status, io_err, level, len, i, j, iret
      integer iyrst,imst,idst,iyrnd,imnd,idnd,ndays,index,temp1,iword
      integer oisst_lun
      logical remove_all

c --- Common block used for writing debug output statements
      integer debug
      integer h_output
      common / bug / debug, h_output

c --- initialize the lun
      l1b_1km_lun  = -1
      geo_1km_lun  = -1
      debug_lun	= -1
      hdr_lun	= -1
      hdf_lun	= -1
      oisst_lun = -1

      remove_all = .true.

c --- initialize dummy variables
      req = ' '
      file = ' '


C 1KM Level 1B MODIS DATA FILES

c --- get the file requirement
      keyname = 'L1B_1KM_REQ'
      keyindex = 1
      status = hdrgetkeystr( filename, keyname, keyindex, req )
      if( status.lt.0 ) req = 'OPTIONAL'
      level = 3
      if( req(1:3).eq.'MAN' ) level = 2

c --- get file access information 
      keyname = 'L1B_1KM_NAME'
      keyindex = 1
      status = hdrgetkeystr( filename, keyname, keyindex, file )

c --- check for error reading info file
      if( status.lt.0 ) then 
         call message( 
     &                 FUNCTION,
     &                'FAILED: Read of 1km file info . ',
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
     &                 '--',
     &                'Using 1km L1B file: '// file ,
     &                 0, 4
     &               )

         l1b_1km_lun = L1B_1KM_UNIT
         goto 140

c ------ Failed to open 
9005     continue 
         call message( 
     &                 FUNCTION,
     &                'FAILED: Could not open 1km L1B file. ', 
     &                 0, level 
     &               )

      endif
140   continue
      if( debug.ne.0 ) then
         call strcompress( file, remove_all, len_file )
         write(*,*) '  L1B 1km file: ',file(1:len_file)
      endif
              


C GEOLOCATION FILES: LATITUDE, LONGITUDE, HEIGHT, ...

c --- get the file requirement
      keyname = 'GEO_1KM_REQ'
      keyindex = 1
      status = hdrgetkeystr( filename, keyname, keyindex, req )
      if( status.lt.0 ) req = 'OPTIONAL'
      level = 3
      if( req(1:3).eq.'MAN' ) level = 2

c --- get file access information 
      keyname = 'GEO_1KM_NAME'
      keyindex = 1
      status = hdrgetkeystr( filename, keyname, keyindex, file )

c --- check for error reading info file
      if( status.lt.0 ) then 
         call message( 
     &                 FUNCTION,
     &                'FAILED: Read of Geo file info . ',
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
     &                 '--',
     &                'Using 1km GEO file: '// file ,
     &                 0, 4
     &               )

         geo_1km_lun = GEO_1KM_UNIT
         goto 150

c ------ Failed to open 
9006     continue 
         call message( 
     &                 FUNCTION,
     &                'FAILED: Could not open GEO file. ',
     &                 0, level 
     &               )

      endif
150   continue
      if( debug.ne.0 ) then
         call strcompress( file, remove_all, len_file )
         write(*,*) '  GEO file: ',file(1:len_file)
      endif
            
C OISST FILE

      do i = 1,360
        do j = 1,180
          oisst(i,j) = -1.0
        end do
      end do
c --- get the file requirement
      keyname = 'OISST_REQ'
      keyindex = 1
      status = hdrgetkeystr( filename, keyname, keyindex, req )
      if( status.lt.0 ) req = 'OPTIONAL'
      level = 3
      if( req(1:3).eq.'MAN' ) level = 2

c --- get file access information 
      keyname = 'OISST_NAME'
      keyindex = 1
      status = hdrgetkeystr( filename, keyname, keyindex, file )

c --- check for error reading info file
      if( status.lt.0 ) then 
         call message( 
     &                 FUNCTION,
     &                'FAILED: Read of OISST file info . ',
     &                 0, level 
     &               )

      else

c ------ issue the OPEN request
         OPEN( 
     &         FILE=file,
     &         UNIT=OISST_UNIT,
     &         STATUS=OISST_STATUS,
     &         FORM=OISST_FORM,
     &         ACCESS=OISST_ACCESS,
     &         RECL=OISST_RECL,
     &         IOSTAT=io_err, 
     &         ERR=9009
     &       ) 

         call message(
     &                 '--',
     &                'Using OISST file: '// file ,
     &                 0, 4
     &               )

         oisst_lun = OISST_UNIT

c This is a bit naughty, but I'm going to do it anyway - read the OISST file
c here and close it, it saves passing the lun all over the place and the info
c on whether it is a mandatory input. -Jim Davies
                    

c ...      Although this was written as a sequential unformatted
c ...      file, it is open and read as a direct access unformatted
c ...      file.  Since the file was written on a big endian machine
c ...      in fortran, there are extra characters between each
c ...      of the different data sets contained in the file.  You
c ...      have no way of knowing this when working on a little
c ...      endian machine.

           iret = -1
           temp1 = read_int(oisst_lun,1,iret)
           iyrst = read_int(oisst_lun,2,iret)
           imst = read_int(oisst_lun,3,iret)
           idst = read_int(oisst_lun,4,iret)
           iyrnd = read_int(oisst_lun,5,iret)
           imnd = read_int(oisst_lun,6,iret)
           idnd = read_int(oisst_lun,7,iret)
           ndays = read_int(oisst_lun,8,iret)
           index = read_int(oisst_lun,9,iret)
           iword = 12
           do j = 1 , 180
             do i=1 , 360
               oisst(i,j) = read_real(oisst_lun,iword,iret)
               iword = iword+1
             enddo
           enddo

           if( iret.ne.0 ) then
               level = 1
               call message( 'get_ancillary_cld',
     &                       'Error reading sst file',
     &                       0, level )
           endif

           CLOSE ( UNIT=oisst_lun )

         goto 180

c ------ Failed to open 
9009     continue 
         call message( 
     &                 FUNCTION,
     &                'FAILED: Could not open OISST file. ',
     &                 0, level 
     &               )

      endif
180   continue
      if( debug.ne.0 ) then
         call strcompress( file, remove_all, len_file )
         write(*,*) '  OISST file: ',file(1:len_file)
      endif

      
C OUTPUT FILES: SST
C ***  NEW RUNTIME ARGUMENT DECIDES WHAT KIND OF OUTPUT FILE TO CREATE
C ********************************************************************

      IF (output_type .eq. 1 .or. output_type .eq. 3) then

c ---   get the file requirement
        keyname = 'MOD28_REQ'
        keyindex = 1
        status = hdrgetkeystr( filename, keyname, keyindex, req )
        if( status.lt.0 ) req = 'OPTIONAL'
        level = 3
        if( req(1:3).eq.'MAN' ) level = 2

c ---   get file access information 
        keyname = 'MOD28_NAME'
        keyindex = 1
        status = hdrgetkeystr(filename, keyname, keyindex,file )
        call strcompress( file, remove_all, len ) 

c ---   check for error reading info file
        if( status.lt.0 ) then 
           call message( 
     &                   FUNCTION,
     &                  'FAILED: Read of SST output file info . ',
     &                   0, level 
     &                 )

        else
                   
c ------   issue the OPEN request
           OPEN( 
     &           FILE=file(1:len),
     &           UNIT=MOD28_UNIT,
     &           STATUS=MOD28_STATUS,
     &           ACCESS=MOD28_ACCESS,
     &           RECL=(max_pixel) * 4,
     &           IOSTAT=io_err, 
     &           ERR=9017
     &         )

           call message(
     &                   '--',
     &                  'Creating SST Output file: '// file ,
     &                   0, 4
     &                 )

           mod28_lun = MOD28_UNIT
           goto 260

c ------   Failed to open 
9017       continue 
           call message( 
     &                   FUNCTION,
     &                  'FAILED: Could not open output file. ',
     &                   0, level 
     &                 )

        endif
260     continue
        if( debug.ne.0 ) then
           call strcompress( file, remove_all, len_file )
           write(*,*) '  SST Output binary file: ',file(1:len_file)
        endif

C ...   OUTPUT FILES: SST HEADER FILE

c ---   get the file requirement
        keyname = 'MOD28_REQ'
        keyindex = 1
        status = hdrgetkeystr( filename, keyname, keyindex, req )
        if( status.lt.0 ) req = 'OPTIONAL'
        level = 3
        if( req(1:3).eq.'MAN' ) level = 2

c ---   get file access information 
        keyname = 'MOD28_HDR'
        keyindex = 1
        status = hdrgetkeystr( filename, keyname, keyindex,file )
        call strcompress( file, remove_all, len ) 


c ---   check for error reading info file
        if( status.lt.0 ) then 
           call message( 
     &                   FUNCTION,
     &                  'FAILED: Read of SST output header file info . ',
     &                   0, level 
     &                 )

        else
                   
c ------   issue the OPEN request
           OPEN( 
     &           FILE=file(1:len),
     &           UNIT=MOD28_HDR_UNIT,
     &           STATUS=MOD28_HDR_STATUS,
     &           ACCESS=MOD28_HDR_ACCESS,
     &           IOSTAT=io_err, 
     &           ERR=9018
     &         ) 

           hdr_lun = MOD28_HDR_UNIT
           goto 261

c ------   Failed to open 
9018       continue 
           call message( 
     &                   FUNCTION,
     &                  'FAILED: Could not open output header file. ',
     &                   0, level 
     &                 )
  
        endif
261     continue
        if( debug.ne.0 ) then
           call strcompress( file, remove_all, len_file )
           write(*,*) '  SST Output binary header file: ',file(1:len_file)
        endif
C **** END OF BINARY FILE CREATION IF BLOCK
C **************************************************************
      endif

C ***  NEW RUNTIME ARGUMENT DECIDES WHAT KIND OF OUTPUT FILE TO CREATE
C ********************************************************************
C ---- HDF output file ---------------------------------------
      if (output_type .eq. 2 .or. output_type .eq. 3) then

c ---   get the file requirement
        keyname = 'MOD28_HDF_REQ'
        keyindex = 1
        status = hdrgetkeystr( filename, keyname, keyindex, req )
        if( status.lt.0 ) level=2

c ---   get file access information
        keyname = 'MOD28_HDF_NAME'
        keyindex = 1
        status = hdrgetkeystr( filename, keyname, keyindex,file )

c ---   check for error reading info file
        if( status.lt.0 ) then
           call message(
     &                   FUNCTION,
     &                  'FAILED: Opening of output SST HDF file. ',
     &                   0, level
     &                 )

        else

c ------   issue the OPEN request - this is an HDF file)
c ...      Open the input HDF file for read/write) 
c ...      (DFACC_READ is defined in hdf.inc)
           hdf_lun = sfstart(file, DFACC_CREATE)
           if (hdf_lun .eq. -1) then
              call message(
     &                     FUNCTION,
     &                    'FAILED: Could not open output hdf SST file. ',
     &                     0, level
     &                     )
           endif
           call message(
     &                   '--',
     &                  'Creating SST HDF Output file: '// file ,
     &                   0, 4
     &                  )
        endif
        if( debug.ne.0 ) then
           call strcompress( file, remove_all, len_file )
           write(*,*) ' Output SST HDF file: ',file(1:len_file)
        endif
C **** END OF HDF OUTPUT FILE CREATION IF BLOCK
C **************************************************************
      endif



C  RUNTIME DEBUG FILE

c --- get the file requirement
      keyname = 'DEBUG_REQ'
      keyindex = 1
      status = hdrgetkeystr( filename, keyname, keyindex, req )
      if( status.lt.0 ) req = 'OPTIONAL'
      level = 3
      if( req(1:3).eq.'MAN' ) level = 2

c --- get file access information 
      keyname = 'DEBUG_NAME'
      keyindex = 1
      status = hdrgetkeystr( filename, keyname, keyindex, file )
      call strcompress( file, remove_all, len ) 

c --- check for error reading info file
      if( status.lt.0 ) then 
         call message( 
     &                 FUNCTION,
     &                'FAILED: Read of RunTime Debug file info. ',
     &                 0, level 
     &               )

      else
                   
c ------ issue the OPEN request
         OPEN( 
     &         FILE=file(1:len),
     &         UNIT=DEBUG_UNIT,
     &         STATUS=DEBUG_STATUS,
     &         FORM=DEBUG_FORM,
     &         ACCESS=DEBUG_ACCESS,
     &         IOSTAT=io_err, 
     &         ERR=9030 
     &       ) 

         debug_lun = DEBUG_UNIT
         goto 270

c ------ Failed to open 
9030     continue 
         call message( 
     &                 FUNCTION,
     &                'FAILED: Open RunTime Debug file. ',
     &                 0, level 
     &               )

      endif
270   continue
      if( debug.ne.0 ) then
         call strcompress( file, remove_all, len_file )
         write(*,*) '  Debug file: ',file(1:len_file)
      endif

      return
      end

c ***************************************************************

        real function sbflip(spin)

c ...   Hal Wolf routine which will flip bytes for a
c ...   real*4 variable (spin)
c ....  version of 10.01.02

        implicit none
        character*1 ch
        character*4 carg
        real spin,spout
        equivalence (carg,spout)

        save

        spout = spin
        ch = carg(1:1)
        carg(1:1) = carg(4:4)
        carg(4:4) = ch
        ch = carg(2:2)
        carg(2:2) = carg(3:3)
        carg(3:3) = ch
        sbflip = spout

        return
        end
