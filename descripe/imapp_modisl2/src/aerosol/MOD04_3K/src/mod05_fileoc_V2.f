*/  Modified by JC Guu  10/11/96
*/  The variables groups, HDFAttNms, and NumHandles
*/  are added to the argument list.

      SUBROUTINE FILEOC(ioflag,modis_flag,handle_trans,
     &                  handle_weight,handle_QC,FN_L1B,
     &                  modfil_L1B,modfil_Geo,modfil_CldMsk,
     &                  modfil_Anc,modfil_mod05,groups, 
     &                  HDFAttNms,NumHandles)
      IMPLICIT NONE
      SAVE
      INCLUDE 'mapi.inc'

*/  Modified by JC Guu 08/23/96
*/  The file name for mod05.inc was mod05.inc.

      INCLUDE 'mod05.inc'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_MODIS_39500.f'

*/  CCR
*/  Modified by JC Guu  12/30/96
*/  PGS_PC.f is included.

      INCLUDE 'PGS_PC.f'

C----------------------------------------------------------------------
C!F77
C
C!Description:
C USING TWO SEPARATE CALLS, THIS SUBROUTINE IS USED BOTH TO OPEN AND 
C CLOSE ALL FILES REQUIRED DURING PRODUCT GENERATION.  TO OPEN ALL 
C FILES (AT THE START OF PROCESSING),  THE SUBROUTINE IS CALLED WITH 
C ARGUMENT IOFLAG SET TO OPEN.  WHEN PROCESSING IS COMPLETE, A SECOND 
C CALL TO FILE_OPN_CLS IS MADE.  THIS TIME, IOFLAG IS SET TO CLOSE.  
C ALL FILES WILL NOW BE CLOSED TO THE PROCESS.
C
C!Input Parameters:
C CHARACTER*5 IOFLAG     Used for either file open or file close
C LOGICAL     modis_flag Used for processing MODIS data
C
C!Output Parameters:
C INTEGER     handle_*   File handle for generic files.
C INTEGER     modfil*    File handle structure for HDF files.
C
C!Revision History:
C $Log: mod05_fileoc.f,v $
c Revision 1.5  1996/12/31  16:31:17  jguu
c The size of HDF_FILENAME is changed to PGSd_PC_FILE_PATH_MAX.
c
c Revision 1.4  1996/11/01  19:32:23  jguu
c *** empty log message ***
c
c Revision 1.2  1996/08/29  19:59:43  jguu
c The file name for the include file mod05.inc is changed to
c mod05.inc. The type of FN_L1B in FILEOC() and the argument
c list for calling CPMFIL() are changed.
c
c Revision 1.1  1996/07/08  18:32:08  vlin
c Initial revision
c
c Revision 1.7  1995/11/21  00:10:58  rhucek
c Added comment within body of code describing alternate actions associated
c with value of "modis_flag" switch.
c
c Revision 1.6  1995/11/07  18:19:53  vlin
c Updated prolog for input & output parameters
c
c Revision 1.5  1995/11/03  19:37:36  vlin
c Moved PCF logical identifiers to "mod05.inc"
c
c Revision 1.4  1995/11/01  15:24:38  vlin
c Before code walkthrough
c
c Revision 1.3  1995/10/18  19:30:11  vlin
c Move all calls to subroutine "MODIS_IO_GEN_OPEN" & "MODIS_IO_GEN_CLOSE"
c from vapor.f to mod05_fileoc.f
c
c Revision 1.1  1995/10/06  14:33:39  vlin
c Initial revision
c
C
C!Team-unique Header:
C
C     This software was developed by the MODIS Science Data Support Team 
C     (SDST) for the National Aeronautics and Space Administration, 
C     Goddard Space Flight Center, under contract NAS5-32373.
C
C REFERENCES and CREDITS:
C
C     Written by 
C     Vicky Lin                 10/06/95
C     Research and Data Systems Corporation
C     SAIC/GSC MODIS Science Data Support Office
C     7501 Forbes Blvd,     Seabrook  MD 20706
C
C     vlin@ltpmail.gsfc.nasa.gov
C
C DESIGN NOTES:
C
C Externals: 
C
C     Functions and Subroutines:
C     PGS_IO_GEN_OPENF                (libPGSTK.a)
C     PGS_IO_GEN_CLOSEF               (libPGSTK.a)
C     PGS_SMF_SetStaticMsg            (libPGSTK.a)
C     PGS_SMF_GenerateStatusReport    (libPGSTK.a)
C     OPMFIL                          (libmapi.a)
C     CLMFIL                          (libmapi.a)
C     CPMFIL                          (libmapi.a)
C
C     Named Constants:
C     MODIS_S_GENERIC               (PGS_MODIS_39500.f)
C     MODIS_F_GENERIC               (PGS_MODIS_39500.f)
C     PGS_S_SUCCESS                 (PGS_SMF.f)
C     PGSd_IO_Gen_RSeqUnf           (PGS_IO.f)
C     PGSd_IO_Gen_WSeqUnf           (PGS_IO.f)
C
C Internals:
C
C     Functions and Subroutines:
C     MODIS_IO_GEN_OPENF
C     MODIS_IO_GEN_CLOSEF
C     CONCATENATE
C     MODIS_SMF_SETDYNAMICMSG
C
C     Variables:
C     MSG*       Temporary variable for holding character value.
C     RTN        Temporary variable for holding integer value.
C     TRANS*,weight,INPUT,EXT,S*,L*,HDFFILE  PCF logical identifier
C     handle_*     Logical unit number
C     FILE_ACCESS  File access type 
C     FILE_VERSION Version of file (minimum value = 1)
C
C!End      
C----------------------------------------------------------------------
C DECLARATION OF VARIABLES

      CHARACTER*5  IOFLAG
      CHARACTER*1024 MSG

*/ Modified by JC Guu 7/22/96
*/ The type of FN_L1B is changed from CHARACTER * 100 TO CHARACTER* (*)

*/  CCR
*/  Modified by JC Guu  12/30/96
*/  The size of HDF_FILENAME is changed to PGSd_PC_FILE_PATH_MAX.

      CHARACTER*(PGSd_PC_FILE_PATH_MAX) HDF_FILENAME
      CHARACTER*100  MSG_BUF1,MSG_BUF2
      CHARACTER* (*) FN_L1B

C by DAC, 9/20/96
 
*/  Modified by JC Guu  10/16/96
*/  The type for the elements of groups and HDFAttNms
*/  is changed to CHARACTER*(*).

      CHARACTER*(*) groups(*)
      CHARACTER*(*) HDFAttNms(*)

      INTEGER  NumHandles
C
      INTEGER  I,RTN,handle_weight,FILE_ACCESS,FILE_VERSION,RECLEN,
     &         modfil_L1B(MODFILLEN),modfil_Geo(MODFILLEN),
     &         modfil_CldMsk(MODFILLEN),modfil_Anc(MODFILLEN),
     &         modfil_mod05(MODFILLEN),modfil_QC(MODFILLEN),
     &         handle_trans(6),handle_QC
      LOGICAL  modis_flag

C DECLARATION OF FUNCTIONS
      INTEGER pgs_pc_getreference, pgs_io_gen_openf

      integer num_args
      integer FlagRA
      character FlagBuff*10
      integer iargc

      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
C      This is the default value
         FlagRA = 0
      endif

C INITIALIZE VARIABLES
      RTN = -5555
      FILE_ACCESS = 0
      HDF_FILENAME = ' '
      FILE_VERSION = 1

      IF (IOFLAG.EQ.'OPEN') THEN
         FILE_ACCESS = PGSd_IO_Gen_RSeqFrm

         CALL MODIS_IO_GEN_OPENF(LRN_TRANS1,FILE_ACCESS,RECLEN,
     2   handle_trans(1),FILE_VERSION)

         CALL MODIS_IO_GEN_OPENF(LRN_TRANS2,FILE_ACCESS,RECLEN,
     2   handle_trans(2),FILE_VERSION)

         CALL MODIS_IO_GEN_OPENF(LRN_TRANS3,FILE_ACCESS,RECLEN,
     2   handle_trans(3),FILE_VERSION)

         CALL MODIS_IO_GEN_OPENF(LRN_TRANS4,FILE_ACCESS,RECLEN,
     2   handle_trans(4),FILE_VERSION)

         CALL MODIS_IO_GEN_OPENF(LRN_TRANS5,FILE_ACCESS,RECLEN,
     2   handle_trans(5),FILE_VERSION)

         CALL MODIS_IO_GEN_OPENF(LRN_TRANS6,FILE_ACCESS,RECLEN,
     2   handle_trans(6),FILE_VERSION)

C         WRITE(*,*) 'HDF_FILENAME, FILE_VERSION'
C         WRITE(*,*) HDF_FILENAME, FILE_VERSION
C         HDF_FILENAME = ' '
C         FILE_VERSION = 1

         CALL MODIS_IO_GEN_OPENF(LRN_weight,FILE_ACCESS,RECLEN,
     2   handle_weight,FILE_VERSION)
      
C Open MOD05 QC file.  First try to open assuming file exists. 
C IF unsuccessful, open file for create (PGSd_IO_Gen_WSeqFrm).

         FILE_ACCESS = PGSd_IO_Gen_USeqFrm
         RTN = pgs_io_gen_openf(LRN_QC,FILE_ACCESS,RECLEN,
     2         handle_QC,FILE_VERSION)
      
c rhucek 02/09/98:  Supplemented status messaging when opening QC file.
         If (rtn .NE. PGS_S_SUCCESS) Then
            msg =  'PGS_IO_GEN_OpenF unable to open MOD05_QC file in update mode, '
     1      // char(10) // 'i.e. file does not exit.  Next try opening file in "create" mode.' 
  
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,MSG,'FILEOC')

            FILE_ACCESS = PGSd_IO_Gen_WSeqFrm 
            rtn = pgs_io_gen_openf(LRN_QC,FILE_ACCESS,RECLEN,handle_QC,FILE_VERSION)

            If (rtn .NE. PGS_S_SUCCESS) Then
               msg = 'PGS_IO_GEN_OpenF able to open MOD05_QC in "create" mode. '
     1         // char(10) // 'Continue processing without access to MOD05_QC.'

               CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,MSG,'FILEOC')
            Else
               msg = 'PGS_IO_GEN_OpenF successfully opened MOD05_QC in "create" mode. '
               CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,'FILEOC')
            End If
         END IF
      
      ELSE IF (IOFLAG.EQ.'CLOSE') THEN
         CALL MODIS_IO_GEN_CLOSEF(LRN_TRANS1,handle_trans(1))

         CALL MODIS_IO_GEN_CLOSEF(LRN_TRANS2,handle_trans(2))

         CALL MODIS_IO_GEN_CLOSEF(LRN_TRANS3,handle_trans(3))

         CALL MODIS_IO_GEN_CLOSEF(LRN_TRANS4,handle_trans(4))

         CALL MODIS_IO_GEN_CLOSEF(LRN_TRANS5,handle_trans(5))

         CALL MODIS_IO_GEN_CLOSEF(LRN_TRANS6,handle_trans(6))

         CALL MODIS_IO_GEN_CLOSEF(LRN_weight,handle_weight)

         CALL MODIS_IO_GEN_CLOSEF(LRN_QC,handle_QC)

      ELSE
         MSG = 
     1   'INVALID VALUE OF IOFLAG'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &   'mod05_fileoc_V2.f')

      END IF

C-----------------------------------------------------------------------
C If MODIS_FLAG is .TRUE., open data files for MODIS processing; 
C otherwise, open data files for "non-MODIS" processing.  Only MODIS
C processing option is presently available.
C-----------------------------------------------------------------------

      IF (IOFLAG.EQ.'OPEN'.and. modis_flag) THEN

C READ L1B FILENAME FROM PCF
       FILE_VERSION = 1
       if( FlagRA .eq. 1) then
          RTN = pgs_pc_getreference(LRN_L1B_RA,FILE_VERSION,FN_L1B)
C CHECK IF L1B FILENAME WAS RETRIEVED SUCCESSFULLY
          WRITE(MSG_BUF1,'(I10)') LRN_L1B_RA 
       else
          RTN = pgs_pc_getreference(LRN_L1B,FILE_VERSION,FN_L1B)
C CHECK IF L1B FILENAME WAS RETRIEVED SUCCESSFULLY
          WRITE(MSG_BUF1,'(I10)') LRN_L1B 
       endif

       MSG_BUF2 = 'RETRIEVED NAME OF Level 1B FILE'
       CALL CONCATENATE(MSG_BUF2, MSG_BUF1, MSG)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
         MSG=
     1   'L1B filename was not retrieved successfully'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       ELSE
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       ENDIF

C OPEN L1B FILE
       RTN = OPMFIL(FN_L1B,'r',modfil_L1B)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Successfully opened L1B file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,
     &    'mod05_fileoc_V2.f')
       ELSE 
          MSG = 
     1   'failed to open L1B file using MAPI OPMFIL'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       END IF

C READ Geolocation FILENAME FROM PCF
       FILE_VERSION = 1
       HDF_FILENAME = ' '
       RTN = pgs_pc_getreference(LRN_Geo,FILE_VERSION,HDF_FILENAME)

C CHECK IF Geolocation FILENAME WAS RETRIEVED SUCCESSFULLY
       WRITE(MSG_BUF1,'(I10)') LRN_Geo
       MSG_BUF2 = 'RETRIEVED NAME OF Geolocation FILE'
       CALL CONCATENATE(MSG_BUF2, MSG_BUF1, MSG)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
         MSG=
     1   'Geolocation filename was not retrieved successfully'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       ELSE
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       ENDIF

C OPEN Geolocation file 
       RTN = OPMFIL(HDF_FILENAME,'r',modfil_Geo)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Successfully opened Geolocation file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,
     &    'mod05_fileoc_V2.f')
       ELSE
         MSG=
     1   'failed to open Geolocation file using MAPI OPMFIL'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       END IF

C READ Cloud Mask FILENAME FROM PCF
       FILE_VERSION = 1
       HDF_FILENAME = ' '
       RTN = pgs_pc_getreference(LRN_CldMsk,FILE_VERSION,HDF_FILENAME)

C CHECK IF Cloud Mask FILENAME WAS RETRIEVED SUCCESSFULLY
       WRITE(MSG_BUF1,'(I10)') LRN_CldMsk
       MSG_BUF2 = 'RETRIEVED NAME OF Cloud Mask FILE'
       CALL CONCATENATE(MSG_BUF2, MSG_BUF1, MSG)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
         MSG=
     1   'Cloud mask filename was not retrieved successfully'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       ELSE
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       ENDIF

C OPEN Cloud Mask FILE
       RTN = OPMFIL(HDF_FILENAME,'r',modfil_CldMsk)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Successfully opened Cloud Mask file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,
     &    'mod05_fileoc_V2.f')
       ELSE
          MSG =
     1   'failed to open cloud mask file using MAPI OPMFIL'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       END IF

C OPEN Ancillary Data File
C       FILE_VERSION = 1
C       HDF_FILENAME = ' '
C       RTN = pgs_pc_getreference(LRN_Anc,FILE_VERSION,HDF_FILENAME)
C
C       WRITE(MSG_BUF1,'(I10)') LRN_Anc
C       MSG_BUF2 = 'RETRIEVED NAME OF ANCILLARY DATA FILE'
C       CALL CONCATENATE(MSG_BUF2, MSG_BUF1, MSG)
C       IF (RTN.NE.PGS_S_SUCCESS) THEN
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
C       ELSE
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,'FILEOC')
C       ENDIF
C       RTN = OPMFIL(HDF_FILENAME,'r',modfil_Anc)
C       IF (RTN.EQ.mapiok) THEN
C          MSG = '... Successfully opened Ancillay data file'
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
C       ELSE
C          MSG = 'MAPI OPMFIL FOR Ancillay data FILE'
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
C       END IF

C READ OUTPUT HDF FILENAME FROM PCF
       FILE_VERSION = 1
       HDF_FILENAME = ' '
       RTN = pgs_pc_getreference(LRN_MOD05,FILE_VERSION,HDF_FILENAME)

C CHECK IF OUTPUT HDF FILENAME WAS RETRIEVED SUCCESSFULLY
       WRITE(MSG_BUF1,'(I10)') LRN_MOD05
       MSG_BUF2 = 'RETRIEVED NAME OF OUTPUT HDF FILE'
       CALL CONCATENATE(MSG_BUF2, MSG_BUF1, MSG)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
         MSG=
     1   'MOD05 HDF filename was not retrieved successfully'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       ELSE
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       ENDIF
CDAC
C
C OPEN MOD05 HDF OUTPUT FILE
C
C       RTN = OPMFIL(HDF_FILENAME,'w',modfil_mod05)
       RTN = OPMFIL(HDF_FILENAME,'a',modfil_mod05)
CDAC

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Successfully opened "mod05.hdf"'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,
     &    'mod05_fileoc_V2.f')
       ELSE
         MSG=
     1   'failed to open MOD05 HDF file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &    'mod05_fileoc_V2.f')
       END IF

      ELSE IF (IOFLAG.EQ.'CLOSE'.and.modis_flag) THEN

C CLMFIL terminates the access of MODIS API routines to a MODIS-HDF
C file created using OPMFIL.  Only pre-existing files should be closed 
C by CLMFIL.

      RTN = CLMFIL(modfil_L1B)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Successfully closed L1B file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,
     &   'mod05_fileoc_V2.f')
      ELSE
         MSG=
     1   'failed to close L1B file using MAPI CLMFIL'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         MSG = 'MAPI CLMFIL FOR L1B FILE'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &   'mod05_fileoc_V2.f')
      END IF

      RTN = CLMFIL(modfil_Geo)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Successfully closed Geolocation file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,
     &   'mod05_fileoc_V2.f')
      ELSE
         MSG=
     1   'failed to close geolocation file using MAPI CLMFIL'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &   'mod05_fileoc_V2.f')
      END IF

      RTN = CLMFIL(modfil_CldMsk)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Successfully closed Cloud Mask file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,
     &   'mod05_fileoc_V2.f')
      ELSE
         MSG=
     1   'failed to close cloud mask file using MAPI CLMFIL'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &   'mod05_fileoc_V2.f')
      END IF

C      RTN = CLMFIL(modfil_Anc)
C
C      IF (RTN.EQ.mapiok) THEN
C         MSG = '... Successfully closed ancillary data file'
C         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
C      ELSE
C         MSG = 'MAPI CLMFIL FOR ancillary data FILE'
C         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
C      END IF

*/  Modified by JC Guu 08/08/96
*/  The argument list for CPMFIL in M-API2.0
*/  is changed.
C      RTN = CPMFIL(modfil_mod05,' ',' ',0)
C
C DAC, 10/7/97
C
      RTN = CPMFIL(modfil_mod05,groups,HDFAttNms,NumHandles)

C      IF(RTN.EQ.0) THEN
C        CALL EXIT(0)
C      ELSE
C        CALL EXIT(-1)
C      END IF

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Successfully closed "mod05.hdf"'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,
     &   'mod05_fileoc_V2.f')
      ELSE
         MSG=
     1   'failed to close mod05 HDF file using MAPI CPMFIL'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,
     &   'mod05_fileoc_V2.f')
      END IF

C------------------------------------------------------
CEND if statement is here
C------------------------------------------------------
      END IF

      RETURN
      END
