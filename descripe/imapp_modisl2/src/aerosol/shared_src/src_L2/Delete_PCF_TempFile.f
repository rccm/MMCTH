      INTEGER FUNCTION Delete_PCF_TempFile(LUN_Tempfile)
      
      IMPLICIT NONE

      INCLUDE 'PGS_IO_1.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Delete a file entry in the temporary I/O section
C	 	of pcf.
C                
C !INPUT PARAMETERS:
C  Integer LUN_Tempfile		Logical Unit Number(LUN) of temporary I/O file
C
C !OUTPUT PARAMETERS:
C  None  
C
C !REVISION HISTORY:
C Revision 1.3  1998/11/18  20:14:30  rhucek
C
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
      character*(*)   		FUNCNAME
      Parameter    	       (FUNCNAME = 'Delete_PCF_TempFile')
      
      integer 			SUCCEED,     FAIL
      Parameter 	       (SUCCEED = 0, FAIL    = -1)

      integer 			LUN_TempFile

C.....local variables
      character*25 		msg25
      character*1024 		msgbuf

      integer rtn, rtn2, fbyte, lbyte, PGS_IO_Gen_Temp_Delete, string_loc

      logical error_flag

c-------------------------------------------------------------------------------
c Begin executable code
c-------------------------------------------------------------------------------
      error_flag          = .FALSE.
      Delete_PCF_TempFile = FAIL

      write(msg25,'(I25)') LUN_Tempfile 
      rtn2 = string_loc(msg25,fbyte,lbyte)

      rtn = PGS_IO_Gen_Temp_Delete( LUN_Tempfile ) 

      If (rtn .EQ. PGSIO_E_GEN_FILE_NODEL) then
         msgbuf = 'No file opened on LUN ' // msg25(fbyte:lbyte)
     1   // ', and none to delete' 

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,msgbuf,FUNCNAME)
      
      Else if (rtn .NE. PGS_S_SUCCESS) then
         error_flag = .TRUE.

         msgbuf = 'deleting temporary file on PCF LUN ' //  msg25(fbyte:lbyte)
     1   // char(10) // 'Operator Action: Check for corrupt PCF file and rerun '
     2   // char(10) // 'as necessary.  If problem persists, notify SDST.'
 
         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      Endif


      if ( .NOT. error_flag) Delete_PCF_TempFile = SUCCEED

      return

      end
