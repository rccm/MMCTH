      SUBROUTINE TEST_FOR_250M_FILE(NO_250)   

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Check for 250 m file - if MISSING, return true
c
c !INPUT PARAMETERS:
c     None
c
c !OUTPUT PARAMETERS:
c     NO_250  Logical variable to tell if 250 meter data is used or not
c
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE
      
      INCLUDE 'mod35.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_PC_9.f'
      
c ... Arguments
      LOGICAL no_250

c ... Local variables

      character*128 qkmfile
      INTEGER ier, file_id, file_version
      CHARACTER*32 routine_name
      PARAMETER ( routine_name = 'test_for_250m_file' )

c ... External functions

      INTEGER pgs_pc_getreference 
      EXTERNAL pgs_pc_getreference 
            
c ... Get QKM file name from process control file

      file_id = LRN_L1b_250
      file_version = 1
      qkmfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, qkmfile )
      if( (ier .ne. PGS_S_SUCCESS) .and. (ier .ne. PGSPC_W_NO_REFERENCE_FOUND) ) then

         call message( routine_name,
     &        'Error getting MODIS 250-meter filename from PCF' //
     &        ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )
      endif

c ... R. Cintineo 28 May 2014 - added ".or. (qkmfile .eq. './MISSING')" to handle case of QKM filename "MISSING"
      if ((qkmfile .eq. ' ') .or. (qkmfile .eq. './MISSING')) no_250=.true.

      END
