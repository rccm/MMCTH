      SUBROUTINE MODIS_IO_GEN_OPENF(FILE_PCF_INDEX,FILE_ACCESS,
     &           RECORD_LENGTH, FILE_HANDLE, FILE_VERSION, rtn_status)

      IMPLICIT NONE

      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_MODIS_39500.f'

C---------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C      SUBROUTINE TO OPEN A FILE USING THE SDP TK GNERIC FILE
C      OPEN FUNCTION.
C
C !INPUT PARAMETERS:
C
C      FILE_PCF_INDEX  LOGICAL IDENTIFIER OF THE FILE IN THE PROCESS
C                      CONTROL FILE (PCF).
C      FILE_ACCESS     SDPTK FILE ACCESS TYPE
C      RECORD_LENGTH   LENGTH OF DATA RECORD (IGNORED, EXCEPT FOR
C                      DIRECT ACCESS FILES)
C      FILE_VERSION    VERSION NUMBER OF FILE
C
C !OUTPUT PARAMETERS:
C
C      FILE_HANDLE    FORTRAN LOGICAL FILE UNIT NUMBER.  THIS IS THE
C                     REFERENCE TO THE FILE THAT IS USED FOR NATIVE
C                     FORTRAN READ AND WRITE STATEMENTS.
C
C      rtn_status     Subroutine return status (-1=fail; 0=success)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C
C      THIS SOFTWARE WAS DEVELOPED BY THE MODIS SCIENCE DATA SUPPORT
C      TEAM FOR THE NATIONAL AERONAUTICS AND SPACE ADMINISTRATION,
C      GODDARD SPACE FLIGHT CENTER, UNDER CONTRACT NAS5-32373
C
C !REFERENCES AND CREDITS:
C
C
C !DESIGN NOTES:
C
C  EXTERNALS:
C
C      FUNCTIONS AND SUBROUTINES:
C          MODIS_SMF_SETDYNAMICMSG (science code)
C          string_locm             (science code)
C          PGS_IO_GEN_OPENF        (libPGSTK.a)
C
C      NAMED CONSTANTS:
C          MODIS_W_GENERIC         (PGS_MODIS_39500.f)
C          MODIS_S_SUCCESS         (PGS_MODIS_39500.f)
C          PGS_S_SUCCESS           (PGS_SMF.f)
C
C !END
C---------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'MODIS_IO_GEN_OPENF')


C DECLARATION OF VARIABLES

      INTEGER RETURNSTATUS
      INTEGER FILE_PCF_INDEX
      INTEGER FILE_ACCESS
      INTEGER RECORD_LENGTH
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
      integer string_loc, PGS_IO_GEN_OPENF



c------------------------
c Initialization
c------------------------

      RETURNSTATUS = 0
      rtn_status = 0

c.....set status message components 
      write(msg25_LUN,'(I25)') FILE_PCF_INDEX
      write(msg25_VRSN,'(I25)') FILE_VERSION

      rtn = string_loc(msg25_LUN, fbyte_LUN, lbyte_LUN)
      rtn = string_loc(msg25_VRSN, fbyte_VRSN, lbyte_VRSN)
      

c-------------------------------------------------------------------------------
c Open generic file using SDPTK.  Report success or failure to LogStatus file.  
c-------------------------------------------------------------------------------

      RETURNSTATUS = PGS_IO_GEN_OPENF(FILE_PCF_INDEX, FILE_ACCESS,
     &               RECORD_LENGTH, FILE_HANDLE, FILE_VERSION)


      IF (RETURNSTATUS .NE. PGS_S_SUCCESS) THEN
         rtn_status = -1

         msgbuf = 'PGS_IO_GEN_OPENF unable to open generic file on LUN = '
     1            // msg25_LUN(fbyte_LUN:lbyte_LUN) // ' and VRSN = '
     2            // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC, msgbuf, FUNCNAME)

      ELSE
         msgbuf = 'Opened generic file on LUN = '
     1            // msg25_LUN(fbyte_LUN:lbyte_LUN) // ' and VRSN = '
     2            // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS, msgbuf, FUNCNAME)

      END IF

      RETURN
      END
