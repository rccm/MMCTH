      SUBROUTINE FILEOC(ioflag, idebug, 
     &                  handle_trans, handle_QC, FN_L1B,
     &                  modfil_L1B, modfil_Geo,
     &                  modfil_CldMsk, modfil_MOD06cd)

      IMPLICIT NONE
      SAVE

      INCLUDE 'cirrus.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_SMF.f'

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
C  rhucek 09/13/2005
C  Replaced parameter LRN_L1BDS (destriped L1B product) with LRN_L1B 
C  throughout.  LRN_L1B references the native MODIS L1B product.  
C  The destriping algorithm sometimes produced unrealistic amplitudes
C  in band 26. 
C  
C  
C  rhucek 04/14/2005
C  Replaced LRN_L1B with LRN_L1BDS throughout.  LRN_L1BDS references
C  the destriped L1B product rather than the original MCST file. 
C
C  rhucek 05/11/2000
C  Added code and logic to open QC file only if debug option is set.
C  If in debug mode, first query the operating system to identify if
C  QC file exists (via Inquire statement).  If QC file exists, open it 
C  in update mode.  Otherwise, open it as a new file.
C 
c  If it does not exist, open a new 
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
c rhucek 4/30/01:  defined parameter FUNCNAME
C DEFINE PARAMETERS 
      character*(*) FUNCNAME
      parameter    (FUNCNAME = 'FILEOC')


C DECLARATION OF VARIABLES
*/ Modified by JC Guu 7/22/96
*/ The type of FN_L1B is changed from CHARACTER * 100 TO CHARACTER* (*)
c  rhucek 05/07/01: allocate full-path file names  using SDPTK PARAMETER for 
c  max file name size
      CHARACTER* (*) FN_L1B
      CHARACTER*5    IOFLAG
      CHARACTER*25   MSG25
      CHARACTER*100  MSG100
      CHARACTER*1024 MSG
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) FN_QC,HDF_FILENAME  

      INTEGER  IDEBUG, RTN,FILE_ACCESS,FILE_VERSION,RECLEN,
     &         modfil_L1B(MODFILLEN),modfil_Geo(MODFILLEN),
     &         modfil_CldMsk(MODFILLEN),
     &         modfil_MOD06cd(MODFILLEN),
     &         handle_trans,handle_QC

*/ fhliang 07/15/98: added following 2 lines
      INTEGER rtn_string_loc, fbyte, fbyte_n, lbyte, lbyte_n,
     &        string_loc

C DECLARATION OF FUNCTIONS
      INTEGER pgs_pc_getreference, pgs_io_gen_openf

      LOGICAL exists

      integer num_args
      integer FlagRA
      character FlagBuff*10
      integer iargc
      
      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
C     This is the default value
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
C   Ping comments out the following two lines because we do not 
C          use water vapor trans table in the new algorithm (4/19/2001)
c         CALL MODIS_IO_GEN_OPENF(LRN_TRANS,FILE_ACCESS,RECLEN,
c     2   handle_trans,FILE_VERSION)

C--------running in debug mode, open QC file
         IF (IDEBUG .eq. 1) THEN 
            write(msg25,'(I25)') LRN_QC
            rtn_string_loc = String_Loc(msg25,fbyte,lbyte)

C-----------check for valid PCF file reference
            FILE_VERSION = 1
            RTN = pgs_pc_getreference(LRN_QC,FILE_VERSION,FN_QC)

            IF (RTN .NE. PGS_S_SUCCESS) THEN
               MSG = 'Cirrus Detection code running in debug mode, but '
     1         // char(10) // 'pgs_pc_getreference unable to retrieve QC file path ' 
     2         // 'name from PCF LUN ' // msg25(fbyte:lbyte) // '.'
     3         // char(10) // 'Operator Action:  Check for corrupted PCF file. ' 
     4         // char(10) // 'If a fault is identified, stage correct PCF file '
     5         // char(10) // 'and rerun PGE.  Otherwise, notify SDST.'
 
               CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)  
     
C-----------open QC file
            ELSE
               rtn_string_loc = String_Loc(FN_QC,fbyte_n,lbyte_n)
               INQUIRE(FILE=FN_QC,EXIST=exists)

C--------------open append mode
               IF (exists) THEN
                  FILE_ACCESS  = PGSd_IO_Gen_USeqFrm
                  FILE_VERSION = 1
                  RTN = pgs_io_gen_openf(LRN_QC,FILE_ACCESS,RECLEN,handle_QC,FILE_VERSION)

                  IF (rtn .NE. PGS_S_SUCCESS) THEN
                     MSG = 'pgs_io_gen_open unable to open existing MOD06QC file on LUN ' 
     1               // char(10) // msg25(fbyte:lbyte) // ': ' // FN_QC(fbyte_n:lbyte_n) // '.'
     2               // char(10) // 'Operator Action:  Check for corrupted PCF file. ' 
     3               // char(10) // 'If a fault is identified, stage correct PCF file '
     4               // char(10) // 'and rerun PGE.  Otherwise, notify SDST.'

                    CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)  
                  END IF

C--------open as new file
               ELSE
                  FILE_ACCESS  = PGSd_IO_Gen_WSeqFrm
                  FILE_VERSION = 1
                  rtn = pgs_io_gen_openf(LRN_QC,FILE_ACCESS,RECLEN,handle_QC,FILE_VERSION)

                  IF (rtn .NE. PGS_S_SUCCESS) THEN
                     MSG = 'pgs_io_gen_openf able to open new MOD06QC file on LUN '
     1               // char(10) // msg25(fbyte:lbyte) // ': ' // FN_QC(fbyte_n:lbyte_n) // '.'
     2               // char(10) // 'Operator Action:  Check for corrupted PCF file or system. ' 
     3               // char(10) // 'problem.  If a fault is identified, stage correct PCF file '
     4               // char(10) // 'and rerun PGE.  Otherwise, notify SDST.'

                    CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
                  END IF 

               END IF   !open new or append

            END IF   !check for valid PCF reference 

         END IF   !check for debug mode

      ELSE IF (IOFLAG.EQ.'CLOSE') THEN

C Ping comments out the following one line because the file is not openned
c        CALL MODIS_IO_GEN_CLOSEF(LRN_TRANS,handle_trans)

         IF (IDEBUG .eq. 1) CALL MODIS_IO_GEN_CLOSEF(LRN_QC,handle_QC)

      ELSE
         MSG = 'IOFLAG ("' // IOFLAG // '") invalid.'
     1   // char(10) // 'Operator Action: Notify SDST.'

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
      END IF  !end open/close check

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
       else
          RTN = pgs_pc_getreference(LRN_L1B,FILE_VERSION,FN_L1B)
       endif

C CHECK IF L1B FILENAME WAS RETRIEVED SUCCESSFULLY
       if( FlagRA .eq. 1) then
          WRITE(MSG100,'(I10)') LRN_L1B_RA
       else
          WRITE(MSG100,'(I10)') LRN_L1B
       endif
       rtn_string_loc = string_loc(MSG100,fbyte,lbyte)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
          MSG = 'pgs_pc_getreference unble to retrieve L1B FILE.'
     2    // char(10) // 'Operator Action: Correct PCF file name'
     3    // ' reference on LUN ' // MSG100(fbyte:lbyte)
     4    // char(10) // 'and rerun code.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
       ELSE
          MSG = 'Retrieved L1B file on LUN ' // MSG100(fbyte:lbyte)
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,FUNCNAME)
       ENDIF

C OPEN L1B FILE
       RTN = OPMFIL(FN_L1B,'r',modfil_L1B)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Opened L1B file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
       ELSE
          MSG = 'MAPI OPMFIL FOR L1B FILE'
     2            // char(10) // 'Operator Action: Stage non-corrupt'
     3            // char(10) // 'version or correct PCF reference to'
     4            // char(10) // 'required file and return PGE'

          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
       END IF

C READ Geolocation FILENAME FROM PCF
       FILE_VERSION = 1
       HDF_FILENAME = ' '
       RTN = pgs_pc_getreference(LRN_Geo,FILE_VERSION,HDF_FILENAME)

C CHECK IF Geolocation FILENAME WAS RETRIEVED SUCCESSFULLY
       WRITE(MSG100,'(I10)') LRN_Geo
       rtn_string_loc = string_loc(MSG100,fbyte,lbyte)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
          MSG = 'pgs_pc_getreference unble to retrieve Geolocation FILE'
     2    // char(10) // 'Operator Action: Correct PCF file name'
     3    // ' reference on LUN ' // MSG100(fbyte:lbyte)
     4    // char(10) // 'and rerun code.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
       ELSE
          MSG = 'Retrieved Geolocation file on LUN '
     2    // MSG100(fbyte:lbyte)
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,FUNCNAME)
       ENDIF

C OPEN Geolocation file
       RTN = OPMFIL(HDF_FILENAME,'r',modfil_Geo)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Opened Geolocation file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
       ELSE
          MSG = 'MAPI OPMFIL FOR Geolocation FILE'
     2            // char(10) // 'Operator Action: Stage non-corrupt'
     3            // char(10) // 'version or correct PCF reference to'
     4            // char(10) // 'required file and return PGE'

          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
       END IF

C READ Cloud Mask FILENAME FROM PCF
       FILE_VERSION = 1
       HDF_FILENAME = ' '
       RTN = pgs_pc_getreference(LRN_CldMsk,FILE_VERSION,HDF_FILENAME)


C CHECK IF Cloud Mask FILENAME WAS RETRIEVED SUCCESSFULLY
       WRITE(MSG100,'(I10)') LRN_CldMsk
       rtn_string_loc = string_loc(MSG100,fbyte,lbyte)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
          MSG = 'pgs_pc_getreference unble to retrieve Cloud Mask FILE.'
     2    // char(10) // 'Operator Action: Correct PCF file name'
     3    // ' reference on LUN ' // MSG100(fbyte:lbyte)
     4    // char(10) // 'and rerun code.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
       ELSE
          MSG = 'Retrieved Cloud Mask file on LUN '
     2    // MSG100(fbyte:lbyte)
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,FUNCNAME)
       ENDIF

C OPEN Cloud Mask FILE
       RTN = OPMFIL(HDF_FILENAME,'r',modfil_CldMsk)

       IF (RTN.EQ.mapiok) THEN
          MSG = '... Opened Cloud Mask file'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
       ELSE
          MSG = 'MAPI OPMFIL FOR Cloud Mask FILE'
     2            // char(10) // 'Operator Action: Stage non-corrupt'
     3            // char(10) // 'version or correct PCF reference to'
     4            // char(10) // 'required file and return PGE'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
       END IF


C READ OUTPUT HDF FILENAME FROM PCF
       FILE_VERSION = 1
       HDF_FILENAME = ' '
       RTN = pgs_pc_getreference(LRN_mod06cd,FILE_VERSION,HDF_FILENAME)

C CHECK IF OUTPUT HDF FILENAME WAS RETRIEVED SUCCESSFULLY
       WRITE(MSG100,'(I10)') LRN_mod06cd
       rtn_string_loc = string_loc(MSG100,fbyte,lbyte)

       IF (RTN.NE.PGS_S_SUCCESS) THEN
          MSG = 'pgs_pc_getreference unble to retrieve OUTPUT HDF FILE.'
     2    // char(10) // 'Operator Action: Correct PCF file name'
     3    // ' reference on LUN ' // MSG100(fbyte:lbyte)
     4    // char(10) // 'and rerun code.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
       ELSE
          MSG = 'Retrieved OUTPUT HDF file on LUN '
     2    // MSG100(fbyte:lbyte)
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,FUNCNAME)
       ENDIF

C OPEN OUTPUT HDF FILE
c      RTN = OPMFIL(HDF_FILENAME,'w',modfil_MOD06cd)
       RTN = OPMFIL(HDF_FILENAME,'a',modfil_MOD06cd)

C  Liqun Ma  Modified the error msg buf   02/18/98
       IF (RTN.EQ.mapiok) THEN
          MSG = '... Opened "mod06cd output file"'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
       ELSE
          MSG = 'MAPI OPMFIL FOR "mod06cd output file"'
     2            // char(10) // 'Operator Action: Stage non-corrupt'
     3            // char(10) // 'version or correct PCF reference to'
     4            // char(10) // 'required file and return PGE'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
       END IF

      ELSE IF (IOFLAG.EQ.'CLOSE') THEN

C CLMFIL terminates the access of MODIS API routines to a MODIS-HDF
C file created using OPMFIL.  Only pre-existing files should be closed
C by CLMFIL.

      RTN = CLMFIL(modfil_L1B)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Closed L1B file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
      ELSE
         MSG = 'MAPI CLMFIL FOR L1B FILE'
     &         // char(10) // 'Operator Action: Check system'
     &         // char(10) // 'disk resources; If adequate,'
     &         // char(10) // 'contact SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
      END IF

      RTN = CLMFIL(modfil_Geo)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Closed Geolocation file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
      ELSE
         MSG = 'MAPI CLMFIL FOR Geolocation FILE'
     &         // char(10) // 'Operator Action: Check system'
     &         // char(10) // 'disk resources; If adequate,'
     &         // char(10) // 'contact SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
      END IF

      RTN = CLMFIL(modfil_CldMsk)

      IF (RTN.EQ.mapiok) THEN
         MSG = '... Closed Cloud Mask file'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
      ELSE
         MSG = 'MAPI CLMFIL FOR Cloud Mask FILE'
     &         // char(10) // 'Operator Action: Check system'
     &         // char(10) // 'disk resources; If adequate,'
     &         // char(10) // 'contact SDST.'

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
      END IF

c-----------han--------------
c      RTN = CLMFIL(modfil_Anc)

c      IF (RTN.EQ.mapiok) THEN
c         MSG = '... Successfully closed ancillary data file'
c         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
c      ELSE
c         MSG = 'MAPI CLMFIL FOR ancillary data FILE'
c         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
c      END IF


c rhucek 7/24/99: replaced call to cpmfil with clmfil
      RTN = CLMFIL(modfil_MOD06cd)

C  Liqun Ma  Modified the error msg buf   02/18/98
      IF (RTN.EQ.mapiok) THEN
         MSG = '... Closed "mod06cd output file"'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,MSG,FUNCNAME)
      ELSE
         MSG = 'MAPI CPMFIL FOR "mod06cd output file"'
     &         // char(10) // 'Operator Action: Check system'
     &         // char(10) // 'disk resources; If adequate,'
     &         // char(10) // 'contact SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,MSG,FUNCNAME)
      END IF

C------------------------------------------------------
C!END if statement is here
C------------------------------------------------------
      END IF

      RETURN
      END

