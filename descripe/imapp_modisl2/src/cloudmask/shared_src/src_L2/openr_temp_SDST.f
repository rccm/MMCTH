      INTEGER FUNCTION OPENR_TEMP( PCFNUM, RECLEN, LUN )
      
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Open a TEMPORARY file for read, direct, unformatted access
C     using the PGS toolkit (use this routine for binary files).
C
C !INPUT PARAMETERS:
C     PCFNUM        File number in Process Control File
C     RECLEN        Record length (bytes)
C
C !OUTPUT PARAMETERS:
C     OPENR_TEMP    Return status flag (0=Success, -1=Failure)
C     LUN           FORTRAN logical unit number for the opened file.
C
C !REVISION HISTORY:
C 
C
C !Team-Unique Header:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Input arguments

      INTEGER pcfnum, reclen
      
c ... Output arguments

      INTEGER lun      

c ... Local variables

      INTEGER status, handle
      
c ... External functions

      INTEGER pgs_io_gen_temp_openf
      EXTERNAL pgs_io_gen_temp_openf
      
c ... Include files for toolkit calls

      include 'PGS_IO.f'
      include 'PGS_SMF.f'

c ... Set return values
      
      openr_temp = -1
      lun = -1
      
c ... Open the file for read, direct, unformatted access
      
      status = pgs_io_gen_temp_openf( PGSd_IO_Gen_NoEndurance, pcfnum,
     &  PGSd_IO_Gen_RDirUnf, reclen, handle )
      if ( status .eq. PGS_S_SUCCESS ) then
        openr_temp = 0
        lun = handle
      endif

      END
