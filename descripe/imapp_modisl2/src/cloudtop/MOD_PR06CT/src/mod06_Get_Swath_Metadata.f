        subroutine mod06_Get_Swath_Metadata(modfil,Cube,scan_type)

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
C       INTEGER  scan_type    Scan type - Day or night
C
C!REVISION HISTORY:  
C $Id: mod06_Get_Swath_Metadata.f,v 1.1 1999/05/24 15:49:50 kis Exp $
C
C!TEAM-UNIQUE HEADER: 
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

      INCLUDE 'mapi.inc' 
      INCLUDE 'PGS_MODIS_39500.f' 
      INCLUDE 'PGS_SMF.f'

C     Declarations

      integer b_size
      parameter (b_size = 8)

c     arguments
      character*1 scan_type
      integer modfil(MODFILLEN)
      integer Cube

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

c ........debug statement ..........................................
      if (debug .gt. 0) then
        write(h_output,'(/,15x,'' Scan Level Metadata '')')
        write(h_output,'(/,2x,'' Scan Type '')')
        write(h_output,'(2x,A10,/)') Scan_type
      endif
c ....................................................................

      return 
      end
