      SUBROUTINE MODIS_IO_GEN_CLOSEF(FILE_PCF_INDEX, FILE_HANDLE,
     &                               FILE_VERSION, rtn_status)


      IMPLICIT NONE
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C       Subroutine to close a file using the SDP TK generic file
C       close function.
C
C !INPUT PARAMETERS:
C
C       FILE_HANDLE - FORTRAN LOGICAL UNIT NUMBER.  THIS IS THE
C                     NATIVE FORTRAN LOGICAL REFERENCE TO THE FILE.
C       FILE_PCF_INDEX - THIS IS THE INDEX NUMBER OF THE FILE AS
C                        LISTED IN THE PROCESS CONTROL FILE (PCF).
C
C       FILE_VERSION  VERSION NUMBER OF FILE
C
C !OUTPUT PARAMETERS:
C       rtn_status    Subroutine return status (-1=fail; 0=success)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C
C      THIS SOFTWARE WAS DEVELOPED BY THE MODIS SCIENCE DATA SUPPORT
C      TEAM FOR THE NATIONAL AERONAUTICS AND SPACE ADMINISTRATION,
C      GODDARD SPACE FLIGHT CENTER, UNDER CONTRACT NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C
C !DESIGN NOTES:
C
C  EXTERNALS:
C
C    FUNCTIONS AND SUBROUTINES:
C        MODIS_SMF_SETDYNAMICMSG (science code)
C        string_locm             (science code)
C        PGS_IO_GEN_CLOSEF       (libPGSTK.a)
C
C    NAMED CONSTANTS:
C        MODIS_W_GENERIC         (PGS_MODIS_39500.f)
C        MODIS_S_SUCCESS         (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS           (PGS_SMF.f)
C
C !END
C-----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'MODIS_IO_GEN_CLOSEF')


C DECLARATION OF VARIABLES

      INTEGER RETURNSTATUS
      INTEGER FILE_PCF_INDEX
      INTEGER FILE_HANDLE
      INTEGER FILE_VERSION
      integer rtn_status


c local variable declarations
      character*1024 msgbuf
      character*25   msg25_LUN,msg25_VRSN

      integer rtn
      integer fbyte_LUN,fbyte_VRSN,
     1        lbyte_LUN,lbyte_VRSN


c function declarations
      integer string_loc, PGS_IO_GEN_CLOSEF



c------------------------
c Initialization
c------------------------

      RETURNSTATUS = 0
      rtn_status   = 0

c.....construct status messsage components
      write(msg25_LUN,'(I25)')  FILE_PCF_INDEX
      write(msg25_VRSN,'(I25)') FILE_VERSION

      rtn = string_loc(msg25_LUN,  fbyte_LUN,  lbyte_LUN)
      rtn = string_loc(msg25_VRSN, fbyte_VRSN, lbyte_VRSN)


c-------------------------------------------------------------------------------
c Close generic file using SDPTK.  Report success or failure to LogStatus file 
c-------------------------------------------------------------------------------

      RETURNSTATUS = PGS_IO_GEN_CLOSEF(FILE_HANDLE)

      IF (RETURNSTATUS .NE. PGS_S_SUCCESS) THEN
         rtn_status = -1

         msgbuf = 'PGS_IO_GEN_CLOSEF unable to close generic file on LUN = '
     1            // msg25_LUN(fbyte_LUN:lbyte_LUN) // ' and VRSN = '
     2            // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC, msgbuf, FUNCNAME)

      ELSE
         msgbuf = 'Closed generic file on LUN = '
     1            // msg25_LUN(fbyte_LUN:lbyte_LUN) // ' and VRSN = '
     2            // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS, msgbuf, FUNCNAME)

      END IF

      RETURN
      END

