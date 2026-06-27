      INTEGER FUNCTION Delete_PCF_TempFile(LUN_Tempfile)

      IMPLICIT NONE

      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Delete a file entry in the temporary I/O section of pcf.
C
C
C !INPUT PARAMETERS:
C  Integer LUN_Tempfile  Logical Unit Number(LUN) of temporary I/O file
C
C
C !OUTPUT PARAMETERS: None

C
C !REVISION HISTORY:
C
C !REFERENCES AND CREDITS:
C
C    Written by        Po Li  Nov. 1998
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
C
C
C    Functions and Subroutines:
C        PGS_IO_Gen_Temp_Delete    (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C
C
C
C    Named Constant:
C	 MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C
C
C !END
C----------------------------------------------------------------------

C.....Function argument
      character*(*)    FUNCNAME
      Parameter       (FUNCNAME = 'Delete_PCF_TempFile')

      integer          SUCCEED,     FAIL
      Parameter       (SUCCEED = 0, FAIL    = -1)

      integer          LUN_TempFile

C.....local variables
      character*25     msg25
      character*1024   msgbuf

      integer rtn, rtn2, fbyte, lbyte, PGS_IO_Gen_Temp_Delete, string_loc


c-------------------------------------------------------------------------------
c Begin executable code
c-------------------------------------------------------------------------------

      rtn = PGS_IO_Gen_Temp_Delete(LUN_Tempfile)

      if (rtn.ne.PGS_S_SUCCESS) then
         write(msg25,'(I25)') LUN_Tempfile
         rtn2 = String_Loc(msg25,fbyte,lbyte)

         msgbuf =
     1   'PGS_IO_Gen_Temp_Delete detected error removing '
     2   // char(10) // 'pcf temporary I/O file on LUN ' // msg25(fbyte:lbyte)
     3   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     4   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'
     5   // char(10) // 'messages produced by routine Delete_PCF_TempFile.'

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC, msgbuf, FUNCNAME)

         Delete_PCF_TempFile = FAIL

      else

         Delete_PCF_TempFile = SUCCEED
      end if

      return

      end
