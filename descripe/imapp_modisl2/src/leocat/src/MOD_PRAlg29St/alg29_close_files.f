      SUBROUTINE ALG29_CLOSE_FILES( LEO_HANDLE, CLD_HANDLE, MET_GROUPS)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C    Close input and output files.
C
C !INPUT PARAMETERS:
C    CLD_HANDLE      MAPI hdf and SDS information on cloud mask file
C    MET_GROUPS      Contains group names for metadata filing
C
C !OUTPUT PARAMETERS:
C    None
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'mapi.inc'
      INCLUDE 'mod35.inc'
      
c ... Arguments

      INTEGER  leo_handle(modfillen),cld_handle( modfillen )
      CHARACTER*(pgsd_met_group_name_l)
     &           met_groups(pgsd_met_num_of_groups)

c ... Local variables

      INTEGER inventorymetadata, archivemetadata
      PARAMETER ( inventorymetadata = 2, archivemetadata = 3 )

      INTEGER rtn, num_of_hdfattrnms
      CHARACTER*32 routine_name
      PARAMETER ( routine_name = 'alg29_close_files' )
      LOGICAL opened
      
      CHARACTER*(max_ecs_name_l-1) hdfattnms( pgsd_met_num_of_groups )

c ... External functions

      integer  pgs_io_gen_closef
      external pgs_io_gen_closef, message

c ... Close the C6 file

      rtn = clmfil(leo_handle)
      if( rtn .eq. 0 ) then
        call message( 'file_close',
     &  'Success closing C6 file ', 0, 3 )
      else
        call message( 'file_close',
     &  'Error closing C6 file' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

c ... Close the output product file
c ... In this case, we need to use cpmfil to close the file in order
c ... to attach the metadata fields to it.

      num_of_hdfattrnms = 2
      hdfattnms( inventorymetadata ) = mecs_core
      hdfattnms( archivemetadata ) = mecs_archive

      rtn = cpmfil(cld_handle,met_groups,hdfattnms,num_of_hdfattrnms)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing output MOD_PR35 file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing output MOD_PR35 file ' //
     &  '[OPERATOR ACTION: Contact MODIS SDST]', 0, 2 )
      endif

      END
