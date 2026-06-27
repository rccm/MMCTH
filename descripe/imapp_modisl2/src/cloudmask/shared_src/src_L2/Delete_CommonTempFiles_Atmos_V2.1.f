      INTEGER FUNCTION Delete_CommonTempFiles_Atmos(LUN_List,NUM_Temp_file)

      IMPLICIT NONE

      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Delete PCF temporary I/O files containing processed atmosphere
C               ancillary data
C
C
C !INPUT PARAMETERS:
C  Integer LUN_List       Array containing a list of Logical Unit Number(LUN) to be deleted
C  Integer NUM_Temp_file  The number of (LUN) to be deleted
C
C
C !OUTPUT PARAMETERS: None
C
C !REVISION HISTORY:
C
C !REFERENCES AND CREDITS:
C
C    Written by Rich Hucek/Po Li  Nov. 1998
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    rhucek@ltpmail.gsfc.nasa.gov
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
C
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        Delete_TempFiles          (science code)
C        String_Loc                (science code)
C
C !END
C----------------------------------------------------------------------

C PARAMETER Statements
      Character*(*)  FUNCNAME
      Parameter     (FUNCNAME = 'Delete_CommonTempFiles_Atmos')

      Integer        SUCCEED,     FAIL
      Parameter     (SUCCEED = 0, FAIL = -1)

      Integer        LUN_List(*), NUM_Temp_file

C Argument Declarations - None


C Local variable/array declarations
      character*25   msg25_LUN
      character*1024 msgbuf

      Integer Delete_PCF_TempFile
      Integer rtn, rtn2, i, LUN, string_loc, fbyte, lbyte

      Logical error_flag


c-------------------------------------------------------------------------------
c Begin executable code
c-------------------------------------------------------------------------------

      error_flag = .FALSE.
      Delete_CommonTempFiles_Atmos = FAIL

C.....Loop over number of temporary I/O files
      Do 100 i = 1, NUM_Temp_file

         LUN = LUN_List(i)

         if ((LUN.LT.1).OR.(LUN.GE.10000.AND.LUN.LE.10999)) THEN

            error_flag = .TRUE.
            write(msg25_LUN,'(I25)') LUN
            rtn2 = string_loc(msg25_LUN,fbyte,lbyte)

            msgbuf =
     1      'PCF LUN number ' // msg25_LUN(fbyte:lbyte) // ' is not vaild range.'
     2      // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else
            rtn = Delete_PCF_TempFile(LUN)

            if (rtn.EQ.FAIL) THEN
                error_flag = .TRUE.
                write(msg25_LUN,'(I25)') LUN
                rtn2 = string_loc(msg25_LUN,fbyte,lbyte)

                msgbuf =
     1          'Delete_TempFiles returned error removing atmos ancillary data file '
     2          // char(10) // 'on PCF Temporary I/O LUN ' // msg25_LUN(fbyte:lbyte)
     3          // char(10) // 'Operator Action:  Refer to prior low level LogStatus error'
     4          // char(10) // 'messages produced by routine Delete_CommonTempFiles.'

                Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            endif

         endif

100   continue

      if ( .NOT.error_flag) Delete_CommonTempFiles_Atmos = SUCCEED

      return

      end
