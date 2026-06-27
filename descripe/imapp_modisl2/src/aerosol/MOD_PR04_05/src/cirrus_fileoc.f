      SUBROUTINE FILEOC_MOD06CD(ioflag,handle_trans,handle_QC,FN_L1B,
     &modfil_L1B,modfil_Geo,modfil_CldMsk,groups,HDFAttNms,NumHandles)
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
C
C!Output Parameters:
C INTEGER     handle_*   File handle for generic files.
C INTEGER     modfil*    File handle structure for HDF files.
C
C!Revision History:
C Revision for cirrus algorithm  1996/11/27 by Dr. Allen Chu
C Change some array names such as modfil_mod05 to modfil_cirrus.
C Include file mod05.inc is changed to cirrus.inc.
C Some other name changes are subject to cirrus output, such as
C cirrus.hdf.
C
C $Log: cirrus_fileoc.f,v $
C  Modified by Liqun Ma   Feb. 1998
C  Modified some error message buf.
C  Modified open file mode from creat to append
C
c Revision 1.2  1996/08/29  19:59:43  jguu
c The file name for the include file cirrus.inc is changed to
c cirrus.inc. The type of FN_L1B in FILEOC() and the argument
c list for calling CPMFIL() are changed.
c
c Revision 1.1  1996/07/08  18:32:08  vlin
c Initial revision
c
c Revision 1.7  1995/11/21  00:10:58  rhucek
c Added comment within body of code describing alternate actions associated
c
c Revision 1.6  1995/11/07  18:19:53  vlin
c Updated prolog for input & output parameters
c
c Revision 1.5  1995/11/03  19:37:36  vlin
c Moved PCF logical identifiers to "cirrus.inc"
c
c Revision 1.4  1995/11/01  15:24:38  vlin
c Before code walkthrough
c
c Revision 1.3  1995/10/18  19:30:11  vlin
c Move all calls to subroutine "MODIS_IO_GEN_OPEN" & "MODIS_IO_GEN_CLOSE"
c from vapor.f to cirrus_fileoc.f
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
C----------------------------------------------------------------------
C
*/  Modified by JC Guu  10/11/96
*/  The variables groups, HDFAttNms, and NumHandles
*/  are added to the argument list.

      IMPLICIT NONE
      SAVE
      INCLUDE 'mapi.inc'
      INCLUDE 'cirrus.inc'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_MODIS_39500.f'

C DECLARATION OF VARIABLES

      CHARACTER*5  IOFLAG
      CHARACTER*1024 MSG

*/ fhliang 07/15/98: added following 2 lines
      integer rtn_string_loc, fbyte, lbyte
      integer string_loc

*/ Modified by JC Guu 7/22/96
*/ The type of FN_L1B is changed from CHARACTER * 100 TO CHARACTER* (*)

      CHARACTER*100  HDF_FILENAME,MSG_BUF1,MSG_BUF2
      CHARACTER* (*) FN_L1B

C by DAC, 9/20/96

*/  Modified by JC Guu  10/16/96
*/  The type for the elements of groups and HDFAttNms
*/  is changed to CHARACTER*(*).

      CHARACTER*(*) groups(*)
      CHARACTER*(*) HDFAttNms(*)

      INTEGER  NumHandles

      INTEGER  RTN,FILE_ACCESS,FILE_VERSION,RECLEN,
     &         modfil_L1B(MODFILLEN),modfil_Geo(MODFILLEN),
     &         modfil_CldMsk(MODFILLEN),
     &         handle_trans,handle_QC

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

c---------------han----------------
         CALL MODIS_IO_GEN_OPENF(LRN_TRANS,FILE_ACCESS,RECLEN,
     2   handle_trans,FILE_VERSION)


C Open Cirrus QC file.  First try to open assuming file does
C exist;  IF unsuccessful, open file for writing.

         FILE_ACCESS = PGSd_IO_Gen_USeqFrm
C         RTN = pgs_io_gen_openf(LRN_QC,FILE_ACCESS,RECLEN,
C     2         handle_QC,FILE_VERSION)

C         IF (rtn .NE. PGS_S_SUCCESS) THEN
C            FILE_ACCESS = PGSd_IO_Gen_WSeqFrm
C            rtn = pgs_io_gen_openf(LRN_QC,FILE_ACCESS,RECLEN,
C     2            handle_QC,FILE_VERSION)
C         END IF

      ELSE IF (IOFLAG.EQ.'CLOSE') THEN

         CALL MODIS_IO_GEN_CLOSEF(LRN_TRANS,handle_trans)

C         CALL MODIS_IO_GEN_CLOSEF(LRN_QC,handle_QC)

      ELSE
         MSG = 'INVALID VALUE OF IOFLAG'
     2   // char(10) // 'Operator Action: Set IOFLAG to either CLOSE or OPEN'

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')

      END IF

C-----------------------------------------------------------------------
C If MODIS_FLAG is .TRUE., open data files for MODIS processing;
C otherwise, open data files for "non-MODIS" processing.  Only MODIS
C processing option is presently available.
C-----------------------------------------------------------------------

      IF (IOFLAG.EQ.'OPEN') THEN

C READ L1B FILENAME FROM PCF
       FILE_VERSION = 1
       if( FlagRA .eq. 1) then
          RTN = pgs_pc_getreference(LRN_L1B_RA,FILE_VERSION,FN_L1B)
C CHECK IF L1B FILENAME WAS RETRIEVED SUCCESSFULLY
          WRITE(MSG_BUF1,'(I10)') LRN_L1B_RA
       else
          RTN = pgs_pc_getreference(LRN_L1B,FILE_VERSION,FN_L1B)
C     CHECK IF L1B FILENAME WAS RETRIEVED SUCCESSFULLY
          WRITE(MSG_BUF1,'(I10)') LRN_L1B
       endif

       rtn_string_loc = string_loc(MSG_BUF1,fbyte,lbyte)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
          MSG = 'pgs_pc_getreference unble to retrieve L1B FILE.'
     2    // char(10) // 'Operator Action: Correct PCF file name'
     3    // ' reference on LUN ' // MSG_BUF1(fbyte:lbyte)
     4    // char(10) // 'and rerun code.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
       ELSE
          MSG = 'Retrieved L1B file on LUN ' // MSG_BUF1(fbyte:lbyte)
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,'FILEOC')
       ENDIF

C OPEN L1B FILE
       RTN = OPMFIL(FN_L1B,'r',modfil_L1B)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Opened L1B file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
       ELSE
          MSG = 'MAPI OPMFIL FOR L1B FILE'
     2            // char(10) // 'Operator Action: Stage non-corrupt'
     3            // char(10) // 'version or correct PCF reference to'
     4            // char(10) // 'required file and return PGE'

          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
       END IF

C READ Geolocation FILENAME FROM PCF
       FILE_VERSION = 1
       HDF_FILENAME = ' '
       RTN = pgs_pc_getreference(LRN_Geo,FILE_VERSION,HDF_FILENAME)

C CHECK IF Geolocation FILENAME WAS RETRIEVED SUCCESSFULLY
       WRITE(MSG_BUF1,'(I10)') LRN_Geo
       rtn_string_loc = string_loc(MSG_BUF1,fbyte,lbyte)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
          MSG = 'pgs_pc_getreference unble to retrieve Geolocation FILE'
     2    // char(10) // 'Operator Action: Correct PCF file name'
     3    // ' reference on LUN ' // MSG_BUF1(fbyte:lbyte)
     4    // char(10) // 'and rerun code.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
       ELSE
          MSG = 'Retrieved Geolocation file on LUN '
     2    // MSG_BUF1(fbyte:lbyte)
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,'FILEOC')
       ENDIF

C OPEN Geolocation file
       RTN = OPMFIL(HDF_FILENAME,'r',modfil_Geo)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Opened Geolocation file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
       ELSE
          MSG = 'MAPI OPMFIL FOR Geolocation FILE'
     2            // char(10) // 'Operator Action: Stage non-corrupt'
     3            // char(10) // 'version or correct PCF reference to'
     4            // char(10) // 'required file and return PGE'

          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
       END IF

C READ Cloud Mask FILENAME FROM PCF
       FILE_VERSION = 1
       HDF_FILENAME = ' '
       RTN = pgs_pc_getreference(LRN_CldMsk,FILE_VERSION,HDF_FILENAME)

C CHECK IF Cloud Mask FILENAME WAS RETRIEVED SUCCESSFULLY
       WRITE(MSG_BUF1,'(I10)') LRN_CldMsk
       rtn_string_loc = string_loc(MSG_BUF1,fbyte,lbyte)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
          MSG = 'pgs_pc_getreference unble to retrieve Cloud Mask FILE.'
     2    // char(10) // 'Operator Action: Correct PCF file name'
     3    // ' reference on LUN ' // MSG_BUF1(fbyte:lbyte)
     4    // char(10) // 'and rerun code.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
       ELSE
          MSG = 'Retrieved Cloud Mask file on LUN '
     2    // MSG_BUF1(fbyte:lbyte)
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,'FILEOC')
       ENDIF

C OPEN Cloud Mask FILE
       RTN = OPMFIL(HDF_FILENAME,'r',modfil_CldMsk)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Opened Cloud Mask file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
       ELSE
          MSG = 'MAPI OPMFIL FOR Cloud Mask FILE'
     2            // char(10) // 'Operator Action: Stage non-corrupt'
     3            // char(10) // 'version or correct PCF reference to'
     4            // char(10) // 'required file and return PGE'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
       END IF

C----han----------
C OPEN Ancillary Data File
C       FILE_VERSION = 1
C       HDF_FILENAME = ' '
C       RTN = pgs_pc_getreference(LRN_Anc,FILE_VERSION,HDF_FILENAME)
C
C       WRITE(MSG_BUF1,'(I10)') LRN_Anc
C       MSG_BUF2 = 'RETRIEVED NAME OF ANCILLARY DATA FILE'
C       CALL CONCATENATE(MSG_BUF2, MSG_BUF1, MSG)
C
C       IF (RTN.NE.PGS_S_SUCCESS) THEN
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
C       ELSE
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,'FILEOC')
C       ENDIF
C
C       RTN = OPMFIL(HDF_FILENAME,'r',modfil_Anc)
C       IF (RTN.EQ.mapiok) THEN
C          MSG = '... Successfully opened Ancillay data file'
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
C       ELSE
C          MSG = 'MAPI OPMFIL FOR Ancillay data FILE'
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
C       END IF
C
C READ OUTPUT HDF FILENAME FROM PCF
C       FILE_VERSION = 1
C       HDF_FILENAME = ' '
C       RTN = pgs_pc_getreference(LRN_mod06cd,FILE_VERSION,HDF_FILENAME)
C
C CHECK IF OUTPUT HDF FILENAME WAS RETRIEVED SUCCESSFULLY
C       WRITE(MSG_BUF1,'(I10)') LRN_mod06cd
C       rtn_string_loc = string_loc(MSG_BUF1,fbyte,lbyte)
C
C       IF (RTN.NE.PGS_S_SUCCESS) THEN
C          MSG = 'pgs_pc_getreference unble to retrieve OUTPUT HDF FILE.'
C     2    // char(10) // 'Operator Action: Correct PCF file name'
C     3    // ' reference on LUN ' // MSG_BUF1(fbyte:lbyte)
C     4    // char(10) // 'and rerun code.'
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
C       ELSE
C          MSG = 'Retrieved OUTPUT HDF file on LUN '
C     2    // MSG_BUF1(fbyte:lbyte)
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,'FILEOC')
C       ENDIF
C
C OPEN OUTPUT HDF FILE
C      RTN = OPMFIL(HDF_FILENAME,'w',modfil_MOD06cd)
C       RTN = OPMFIL(HDF_FILENAME,'a',modfil_MOD06cd)
C
C  Liqun Ma  Modified the error msg buf   02/18/98
C       IF (RTN.EQ.mapiok) THEN
C          MSG = '... Opened "mod06cd output file"'
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
C       ELSE
C          MSG = 'MAPI OPMFIL FOR "mod06cd output file"'
C     2            // char(10) // 'Operator Action: Stage non-corrupt'
C     3            // char(10) // 'version or correct PCF reference to'
C     4            // char(10) // 'required file and return PGE'
C          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
C       END IF

      ELSE IF (IOFLAG.EQ.'CLOSE') THEN

C CLMFIL terminates the access of MODIS API routines to a MODIS-HDF
C file created using OPMFIL.  Only pre-existing files should be closed
C by CLMFIL.

      RTN = CLMFIL(modfil_L1B)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Closed L1B file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
      ELSE
         MSG = 'MAPI CLMFIL FOR L1B FILE'
     &         // char(10) // 'Operator Action: Check system'
     &         // char(10) // 'disk resources; If adequate,'
     &         // char(10) // 'contact SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
      END IF

      RTN = CLMFIL(modfil_Geo)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Closed Geolocation file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
      ELSE
         MSG = 'MAPI CLMFIL FOR Geolocation FILE'
     &         // char(10) // 'Operator Action: Check system'
     &         // char(10) // 'disk resources; If adequate,'
     &         // char(10) // 'contact SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
      END IF

      RTN = CLMFIL(modfil_CldMsk)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Closed Cloud Mask file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
      ELSE
         MSG = 'MAPI CLMFIL FOR Cloud Mask FILE'
     &         // char(10) // 'Operator Action: Check system'
     &         // char(10) // 'disk resources; If adequate,'
     &         // char(10) // 'contact SDST.'

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
      END IF

c-----------han--------------
c      RTN = CLMFIL(modfil_Anc)
C
c      IF (RTN.EQ.mapiok) THEN
c         MSG = '... Successfully closed ancillary data file'
c         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
c      ELSE
c         MSG = 'MAPI CLMFIL FOR ancillary data FILE'
c         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
c      END IF
C
*/  Modified by JC Guu 08/08/96
*/  The argument list for CPMFIL in M-API2.0
*/  is changed.
C      RTN = CPMFIL(modfil_cirrus,' ',' ',0)
C
C DAC, 9/19/96
C
C      RTN = CPMFIL(modfil_MOD06cd,groups,HDFAttNms,NumHandles)
C
C  Liqun Ma  Modified the error msg buf   02/18/98
C      IF (RTN.EQ.mapiok) THEN
C         MSG = '... Closed "mod06cd output file"'
C         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,'FILEOC')
C      ELSE
C         MSG = 'MAPI CPMFIL FOR "mod06cd output file"'
C     &         // char(10) // 'Operator Action: Check system'
C     &         // char(10) // 'disk resources; If adequate,'
C     &         // char(10) // 'contact SDST.'
C         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,'FILEOC')
C      END IF

C------------------------------------------------------
C!END if statement is here
C------------------------------------------------------
      END IF

      RETURN
      END

