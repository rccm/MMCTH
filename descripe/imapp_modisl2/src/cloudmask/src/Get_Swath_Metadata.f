        subroutine Get_Swath_Metadata(modfil,Cube,scan_number,
     &                                scan_pixels,scan_type,
     &                                Sstart_time,nadir)

        IMPLICIT NONE
        SAVE

c     Common Block for debugging code
      common / bug / debug, h_output

C!F77------------------------------------------------------------------
C
C!DESCRIPTION:  Extracts Swath metadata for the given scan.
C
C!INPUT PARAMETERS:   
C
C       INTEGER  modfil(3)  Passes vital stats about input L1B granule
C       INTEGER  Cube        Starting record number for extracting dta
c
C!OUTPUT PARAMETERS: 
C       INTEGER  scan_number  Scan number in granule
C       CHAR     scan_pixels  Number of pixels in this scan
C       INTEGER  scan_type    Scan type - Day or night
C       DBL      Sstart_time  Starting scan time (TAI)
C       INTEGER  nadir        Nadir pixel number
C
C!REVISION HISTORY:  
C
C!TEAM-UNIQUE HEADER: 
c This software is developed by the MODIS Science Data Support
c Team for the National Aeronautics and Space Administration,
c Goddard Space Flight Center, under contract NAS5-32373.
c
C!REFERENCES and CREDITS:
C
C!DESIGN NOTES:
C
C    Externals:
C       Functions:
C          GMTBL            (mapi.inc) 
c       Subroutines:
c          message          Utility which writes message to LogStatus
c                           file
C!END-----------------------------------------------------------------

c following line was commented by fhliang on 01/05/98
c     INCLUDE 'hdf.inc' 
      INCLUDE 'mapi.inc' 
      INCLUDE 'PGS_MODIS_39500.f' 
      INCLUDE 'PGS_SMF.f'

C     Declarations

      integer b_size
      parameter (b_size = 8)

c     arguments
      character*1 scan_type
      double precision Sstart_time
      integer modfil(MODFILLEN)
      integer Cube,scan_number,scan_pixels,nadir

c     local variables
      character*64  table_name, field_name
      integer   rtn, nrecords,buff_size,h_output,debug,Start

c     initialize variables needed for extraction of all metadata
      table_name = 'Level 1B Swath Metadata'
      nrecords = 1
      rtn = 0
      buff_size = b_size
c ... Actual starting record is not 1, it is zero
      Start = Cube - 1

c     use mapi routines to extract correct information
c     First call to get Scan Number in granule
      field_name = 'Scan Number'
c     call first to get correct buff_size filled in
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,scan_number)
c     second call actually gets data
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,scan_number)
      if (rtn .ne. 0) then
        call message( 'Get_Swath_Metadata',
     &    'Error getting scan number from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]',
     &    0, 2 )
      end if

c     use mapi routines to extract correct information
c     First call to get Number of pixels in this scan
      buff_size = b_size
      field_name = 'EV_Frames'
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,scan_pixels)
c     second call actually gets data
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,scan_pixels)
      if (rtn .ne. 0) then
        call message( 'Get_Swath_Metadata',
     &    'Error getting EV_Frames from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]',
     &    0, 2 )
      end if
c
c     use mapi routines to extract correct information
c     First call to get Scan type (day or night mode)
      buff_size = b_size
      field_name = 'Scan Type'
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,scan_type)
c     second call actually gets data
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,scan_type)
      if (rtn .ne. 0) then
        call message( 'Get_Swath_Metadata',
     &    'Error getting Scan Type from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]',
     &    0, 2 )
      end if

c     use mapi routines to extract correct information
c     First call to get the starting time for this particular scan
      buff_size = b_size
      field_name = 'EV Sector Start Time'
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,Sstart_time)
c     second call actually gets data
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,Sstart_time)
      if (rtn .ne. 0) then
        call message( 'Get_Swath_Metadata',
     &    'Error getting Start Time (TAI) from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]',
     &    0, 2 )
      end if

c     use mapi routines to extract correct information
c     First call to get the starting time for this particular scan
      buff_size = b_size
      field_name = 'Nadir_Frame_Number'
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,nadir)
c     second call actually gets data
      rtn = GMTBL(modfil,table_name,'  ',field_name,Start,nrecords,
     *            buff_size,nadir)
      if (rtn .ne. 0) then
        call message( 'Get_Swath_Metadata',
     &    'Error getting nadir pixel number from L1B file' // 
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]',
     &    0, 2 )
      end if

c ........debug statement ..........................................
      if (debug .gt. 0) then
        write(h_output,'(/,15x,'' Scan Level Metadata '')')
        write(h_output,'(/,2x,'' Record No.  Scan number  Last pix ''
     +        // '' Nadir  Stype     Start Time'')')
        write(h_output,'(2x,4I10,A5,3x,d17.10,/)') Start,
     +    scan_number, scan_pixels, nadir,scan_type,
     +    Sstart_time
      endif
c ....................................................................

      return 
      end
