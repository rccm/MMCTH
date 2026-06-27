      INTEGER FUNCTION Delete_CommonTempFiles_Atmos()

      IMPLICIT NONE

      INCLUDE 'PGS_IO_1.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'Atmos_AncData.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Delete PCF temporary I/O files containing processed atmosphere
C		ancillary data
C
C
C !INPUT PARAMETERS:
C  None
C
C
C !OUTPUT PARAMETERS:
C  None
C
C !REVISION HISTORY:
c Revision 1.3  1998/11/18  20:12:10  rhucek
c
C
C
C !REFERENCES AND CREDITS:
C
C    Written by Rich Hucek/Po Li  Nov. 1998
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    lipo@ltpmail.gsfc.nasa.gov
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C
C !DESIGN NOTES:
C
C  SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C    Named Constants:
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        MODIS_M_GENERIC           (PGS_MODIS_39500.f)
C	 NUM_ATMOS_COMMON_ANC	   (Atmos_AncData.inc)
C        PGSIO_E_GEN_FILE_NODEL    (PGS_IO_1.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        Delete_TempFiles          (science code)
C        String_Loc		   (science code)
C
C    Variables:
C	 LUN_Temp_AtmosCommonAnc   (Atmos_AncData.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER Statements
      Character*(*)  FUNCNAME
      Parameter     (FUNCNAME = 'Delete_CommonTempFiles_Atmos')

      Integer        SUCCEED,     FAIL
      Parameter     (SUCCEED = 0, FAIL = -1)


C Argument Declarations - None


C Local variable/array declarations
      character*25   msg25_LUN
      character*1024 msgbuf

      Integer PGS_IO_Gen_Temp_Delete
      Integer rtn, rtn2, i, LUN, string_loc, fbyte, lbyte

      Logical error_flag


c-------------------------------------------------------------------------------
c Begin executable code
c-------------------------------------------------------------------------------
      error_flag = .FALSE.
      Delete_CommonTempFiles_Atmos = FAIL

C.....Loop over number of temporary I/O files
      Do 100 i = 1, NUM_ATMOS_COMMON_ANC
         LUN = LUN_Temp_AtmosCommonAnc(i)

         write(msg25_LUN,'(I25)') LUN
	 rtn2 = string_loc(msg25_LUN,fbyte,lbyte)

         rtn = PGS_IO_Gen_Temp_Delete( LUN ) 

         If (rtn .EQ. PGSIO_E_GEN_FILE_NODEL) then
            msgbuf = 'no temporary file opened on LUN ' // msg25_LUN(fbyte:lbyte)
     1      // ' ... no file to delete' 

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,msgbuf,FUNCNAME)
      
         Else if (rtn .NE. PGS_S_SUCCESS) then
            error_flag = .TRUE.

            msgbuf = 'deleting temporary file on PCF LUN ' //  msg25_LUN(fbyte:lbyte)
     1      // char(10) // 'Operator Action: Check for corrupt PCF file and rerun '
     2      // char(10) // 'as necessary.  If problem persists, notify SDST.'
 
            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         Endif

100   Continue

      if ( .NOT.error_flag) Delete_CommonTempFiles_Atmos = SUCCEED

      return

      end
