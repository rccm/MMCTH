      SUBROUTINE MOD07_CLOSE_FILES( CUBE_HANDLE, GEO_HANDLE, CLD_HANDLE, 
     &  OUT_HANDLE, MET_GROUPS)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C    Close input and output files.
C
C !INPUT PARAMETERS:
C    CUBE_HANDLE     MAPI hdf and SDS information on L1B file
C    GEO_HANDLE      MAPI hdf and SDS information on geolocation file
C    CLD_HANDLE      MAPI hdf and SDS information on cloud mask file
C    OUT_HANDLE      MAPI hdf and SDS information on output  file
C    MET_GROUPS      Contains group names form metadata filing
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
      INCLUDE 'mod07_pcfnum.inc'
      
c ... Arguments

      INTEGER cube_handle( modfillen ), geo_handle( modfillen ), 
     &        cld_handle( modfillen ), out_handle( modfillen )
      CHARACTER*(pgsd_met_group_name_l)
     &           met_groups(pgsd_met_num_of_groups)

c ... Local variables

      INTEGER archivemetadata, inventorymetadata
      PARAMETER ( inventorymetadata = 2, archivemetadata = 3 )

      INTEGER rtn, num_of_hdfattrnms
      CHARACTER*32 routine_name
      PARAMETER ( routine_name = 'mod07_close_files' )
      LOGICAL opened
      
      CHARACTER*(max_ecs_name_l-1) hdfattnms( pgsd_met_num_of_groups )

c ... External functions

      integer  pgs_io_gen_closef
      external pgs_io_gen_closef, message
            
c ... Close the 1km L1B file

      rtn = clmfil(cube_handle)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing scan cube 1km data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing scan cube 1km data file ' //
     &  '[OPERATOR ACTION: Contact MODIS SDST]', 0, 2 )
      endif

c ... Close the geolocation file

      rtn = clmfil(geo_handle)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing geolocation data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing geolocation data file ' //
     &  '[OPERATOR ACTION: Contact MODIS SDST]', 0, 2 )
      endif

c ... Close the cloud mask input file

      rtn = clmfil(cld_handle)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing cloud mask data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing cloud mask data file ' //
     &  '[OPERATOR ACTION: Contact MODIS SDST]', 0, 2 )
      endif

c ... Close view angle and regresssion coefficient files

      inquire( iuctw, opened = opened )
      if( opened ) call modis_io_gen_closef( trc_pcfnum, iuctw )

      inquire( iuang, opened = opened )
      if( opened ) call modis_io_gen_closef( ang_pcfnum, iuang )

c ... Close the output product file
c ... In this case, we need to use cpmfil to close the file in order
c ... to attach the metadata fields to it.

      num_of_hdfattrnms = 2
      hdfattnms( inventorymetadata ) = mecs_core
      hdfattnms( archivemetadata ) = mecs_archive

      rtn = cpmfil(out_handle,met_groups,hdfattnms,num_of_hdfattrnms)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing output MOD_PR07 file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing output MOD_PR07 file ' //
     &  '[OPERATOR ACTION: Contact MODIS SDST]', 0, 2 )
      endif

      END
