      PROGRAM Driver_MOD06CD

      IMPLICIT NONE
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'cirrus.inc'

C-----------------------------------------------------------------------
C !F77
C
C !Description: This program is a driver for MODIS MOD06CD products.
C
C
C !Input Parameters: none
C
C !Output Parameters: none
C
C !Revision History:
c 05/19/1999 fhliang
c changed the declaration for LRN_PCF_Attr.
c
C  Modified by Liqun Ma   03/11/98
C  Updated the logical and some error msgs.
C
C  Modified by Liqun Ma   02/18/98
C  Some unreferenced variables and unused include files are moved out.
C  Some temporary output were moved out
C  Prolog was updated
C  Replaced calling Set_CoreMetadata and Set_ArchiveMetadata
C  to calling Update_InvMet_MOD06 and Update_ArchMet_MOD06
C  Added input file check section.
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Science Team
C   for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center.
C
C !REFERENCES and CREDITS:
C
C   Written by
C   Dr. Wei Han                                    11/24/97
C   SFA, Inc.
C   Naval Research Laboratory
C   Code 7212
C   Washington, DC 20375
C
C !DESIGN NOTES:
C
C   At present, this program is set up to process one granule of MODIS
C   sensor data.
C
C
C   Functions and Subroutines:
c
C       Process_Mod06CD
C       Set_InvMet_MOD06
C       Set_ArchMet_MOD06
C       MODIS_SMF_SETDYNAMICMSG
C       pgs_pc_getreference
C       pgs_met_getpcattr_s
C       OPMFIL
C       CPMFIL
C
C !END
C-----------------------------------------------------------------------

C.....Declare and define FORTRAN PARAMETERs
      CHARACTER*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'Driver_MOD06CD.f')

      INTEGER        ARCHIVEMETADATA
      PARAMETER     (ARCHIVEMETADATA = 3)

      INTEGER        INVENTORYMETADATA
      PARAMETER     (INVENTORYMETADATA = 2)

      INTEGER        MAX_NUM_PSA,               NUM_PSA
      PARAMETER     (MAX_NUM_PSA = 16,          NUM_PSA = 0)

      INTEGER        MAX_NUM_MeasParm,          NUM_MeasParm
      PARAMETER     (MAX_NUM_MeasParm = 2,      NUM_MeasParm = 0)

      INTEGER        MAX_NUM_MODIS_Inputs,      NUM_MODIS_InputFiles
      PARAMETER     (MAX_NUM_MODIS_Inputs = 30, NUM_MODIS_InputFiles = 3)

      INTEGER        MAX_NUM_ArchivePSA_SC,     NUM_ArchivePSA_SC
      PARAMETER     (MAX_NUM_ArchivePSA_SC = 8, NUM_ArchivePSA_SC = 0)

      INTEGER        NUM_Input_LRN
C......Ping Yang made 1 modification (4/19/2001)
C      PARAMETER     (NUM_Input_LRN = 4)
      PARAMETER     (NUM_Input_LRN = 3)

C.....function declarations
      INTEGER pgs_met_getpcattr_s,
     1        pgs_pc_getreference

      INTEGER Set_ArchMet_MOD06,
     1        Set_InvMet_MOD06,
     2        chk_input_L2

c.....other declarations
      CHARACTER*512 char_buf, msgbuf
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) HDF_FILENAME

      INTEGER i, rtn, FILE_VERSION, NumHandles, Vrsn_No
      INTEGER Modfil_mod06cd(MODFILLEN)

      LOGICAL error_flag, error_flag2, ExtGeoPntr_Flag


C---For metadata setup---------------------------
      CHARACTER*256                     HDFAttrName, ECSParmName
      CHARACTER*(MAX_ECS_NAME_L-1)      HDFAttNms(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(PGSd_MET_GROUP_NAME_L) MET_Handles(PGSd_MET_NUM_OF_GROUPS)

c.....Measured Parameters
      CHARACTER*30 Name_MeasParm(MAX_NUM_MeasParm),
     *             Auto_QA_Flag (MAX_NUM_MeasParm),
     *             Auto_QA_Expl (MAX_NUM_MeasParm)
      INTEGER      QA_Miss(2)

c.....Input Pointers
      INTEGER Input_LRN(NUM_Input_LRN),Input_Vrsn(NUM_Input_LRN)

c.....MODIS Input Files
      INTEGER VRSN_MODIS_InputFiles(MAX_NUM_MODIS_Inputs)
      INTEGER  LRN_MODIS_InputFiles(MAX_NUM_MODIS_Inputs)

c.....PSAs
      CHARACTER*30 PSAName (MAX_NUM_PSA)
      REAL         PSAValue(MAX_NUM_PSA)

c.....Product Specific (or Archive) PSAs
      CHARACTER*100  Name_ArchivePSA_SC(MAX_NUM_ArchivePSA_SC)
      REAL          Value_ArchivePSA_SC(MAX_NUM_ArchivePSA_SC)

      integer num_args
      integer FlagRA
      character FlagBuff*10
      integer iargc

      data ExtGeoPntr_Flag /.true. /
      data QA_Miss / 20,20 /
      
      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
C     This is the default value
         FlagRA = 0
      endif
      
C Added by Liqun Ma
C.....Opening message
      msgbuf = 'MODIS process MOD_PR06CD initiated.'
      Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,msgbuf,FUNCNAME)

C Initialization
      error_flag = .FALSE.
      error_flag2 = .FALSE.

C Check Start and Stop time from PCF and MODIS input files
C (L1B, geolocation, cloud mask )
C modified by Liqun Ma 02/18/98: write the erorr message to LogStates.
C
      Vrsn_No = 1

      if( FlagRA .eq. 1) then
         rtn=chk_input_L2(LRN_L1B_RA,Vrsn_No)
      else
         rtn=chk_input_L2(LRN_L1B,Vrsn_No)
      endif


      If (rtn.ne.MAPIOK) Then
         error_flag = .true.
         msgbuf = 'L1B file does not pass time checking'
     2            // char(10) // 'Operator Action: Stage MODIS L1B'
     3            // char(10) // 'granule that matches collection'
     4            // char(10) // 'period, and rerun code. If L1B'
     5            // char(10) // 'granule and collection period'
     6            // char(10) // 'match, contact SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf

C rhucek 9/12/05: L1B destriped product no longer used as input 
C
C rhucek 4/14/05: check that destriped L1B product is available on 
C LUN 430000 
C     rtn=chk_input_L2(LRN_L1BDS,Vrsn_No)
C
C     If (rtn.ne.MAPIOK) Then
C        error_flag = .true.
C        msgbuf = 'Destriped L1B file does not pass time checking'
C    2            // char(10) // 'Operator Action: Check for destriped' 
C    3            // char(10) // 'file on LUN 430000. Notify SDST.' 
C
C        Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
C     EndIf

      rtn=chk_input_L2(LRN_GEO,Vrsn_No)

      If (rtn.ne.MAPIOK) Then
         error_flag = .true.
         msgbuf = 'Geolocation file does not pass time checking'
     2            // char(10) // 'Operator Action: Stage MODIS'
     3            // char(10) // 'geolocation granule that matces'
     4            // char(10) // 'collection period, and rerun code.'
     5            // char(10) // 'If geolocation granule and collection'
     6            // char(10) // 'period match, contact SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf

      rtn=chk_input_L2(LRN_CldMsk,Vrsn_No)
      If (rtn.ne.MAPIOK) Then
         error_flag = .true.
         msgbuf = 'Cloud mask file does not pass time checking'
     2            // char(10) // 'Operator Action: Stage MODIS cloud'
     3            // char(10) // 'mask granule that matches collection'
     4            // char(10) // 'period, and rerun code. If cloud'
     5            // char(10) // 'mask granule and collection period'
     6            // char(10) // 'match, contact SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


C     IF anyone or all required inputs are not available, exit now.
      If (error_flag)  Then
         msgbuf = 'Required inputs not available, MOD_PR06CD exit(1)'
         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         Call exit(1)
      EndIf




C--Retrieve value of ECS Parameter 'DAYNIGHTFLAG' from HDF attribute
C--'CoreMetadata.0' in Geolocation Product.  DAYNIGHTFLAG contains
C--a string value.

       Vrsn_No     = 1
       HDFAttrName = MECS_CORE
       ECSParmName = MCORE_DAYNIGHTFLAG

       rtn=pgs_met_getpcattr_s(LRN_GEO,Vrsn_No,HDFAttrName,ECSParmName,char_buf)

      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .true.
         msgbuf = 'pgs_met_getpcattr_s unable to rtrieve'
     2       // char(10) // 'DAYNIGHTFLAG from Geolocation product'
     3       // char(10) // 'Operator Action: Stage valid'
     4       // char(10) // 'geolocation file with ECS metadata field'
     5       // char(10) // 'DAYNIGHTFLAG set.  Rerun code.'

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      Else

C--------If (char_buf .ne. 'Night') Then
C--------Comment from lipo, add two more night flags: "NIGHT" and "night"

         If (char_buf .ne. 'Night'
     1      .and. char_buf .ne. 'night'
     2      .and. char_buf .ne. 'NIGHT') Then

C-----------Process MOD_PR06CD
            CALL Process_Mod06CD(error_flag2)

            If (error_flag2)  Then     !never happens
               error_flag = .true.

               msgbuf =
     1         'Process_Mod06CD failed.'
     2         // char(10) // 'Operator Action: Notify SDST. '

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            EndIf


C-----------Set Core Metadata parameters
            do 10 I = 1, NUM_Input_LRN
               Input_Vrsn(i) = 1
   10       continue

            do 20 I = 1, PGSd_MET_NUM_OF_GROUPS
               MET_Handles(I)= ' '
   20       continue

            Input_LRN(1) = LRN_Geo
            if( FlagRA .eq. 1) then
               Input_LRN(2) = LRN_L1B_RA
            else
               Input_LRN(2) = LRN_L1B
            endif

            Input_LRN(3) = LRN_CldMsk

C            Input_LRN(4) = LRN_TRANS   ! changed by ping (4/19/20001)

            do 40 I = 1, NUM_MODIS_InputFiles
               LRN_MODIS_InputFiles(I)  = Input_LRN(I)
               VRSN_MODIS_InputFiles(I) = 1
   40       continue
C-----------End of metadata setup---------------------------


C-----------Set ECS Inventory Metadata
            rtn = Set_InvMet_MOD06( ExtGeoPntr_Flag,
     1                              NUM_Input_LRN,Input_LRN,Input_Vrsn,
     2                              NUM_MeasParm,Name_MeasParm,QA_Miss,
     3                              Auto_QA_Flag,Auto_QA_Expl,
     4                              NUM_PSA,PSAName,PSAValue,
     5                              MET_Handles)

            If (rtn.ne.mapiok) Then
               error_flag = .true.

               msgbuf = 'Set_InvMet_MOD06 failed. '
     1         // char(10) // 'Operator Action: Notify SDST'

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC, msgbuf,FUNCNAME)
            Else
               msgbuf = 'Set_InvMet_MOD06 passed without error'
               call modis_smf_setdynamicmsg(MODIS_S_SUCCESS,msgbuf,FUNCNAME)
            EndIf


C-----------Set ECS Archive Metadata
            rtn = Set_ArchMet_MOD06( MET_Handles,
     1                               NUM_MODIS_InputFiles,
     2                               LRN_MODIS_InputFiles,
     3                               VRSN_MODIS_InputFiles,
     4                               NUM_ArchivePSA_SC,
     5                               Name_ArchivePSA_SC,
     6                               Value_ArchivePSA_SC )

            If (rtn.ne.mapiok) Then
               error_flag = .true.

               msgbuf = 'Set_ArchMet_MOD06 failed'
     2         // char(10) // 'Operator Action: Notify SDST'

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
     2              msgbuf,FUNCNAME)
            Else
               msgbuf = 'Set_ArchMet_MOD06 passed without error'
               call modis_smf_setdynamicmsg(MODIS_S_SUCCESS,msgbuf,FUNCNAME)
            EndIf


C----Open HDF product file
            FILE_VERSION=1
            rtn = pgs_pc_getreference(LRN_MOD06cd,FILE_VERSION,HDF_FILENAME)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf = 'pgs_p c_getreference unable to get MOD06 product'
     &         // char(10) // 'file name for ECS metadata update.'
     &         // char(10) // 'Operator Action: Correct PCF file'
     &         // char(10) // 'name reference and rerun code.'

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            Else
               rtn = OPMFIL(hdf_filename,'a',modfil_mod06cd)

               If (rtn.ne.mapiok) Then
                  error_flag = .true.

                  msgbuf = 'OPMFIL Failed to open output HDF file'
     &            // 'for ECS metadata update '
     &            // char(10) // 'Operator Action: Stage'
     &            // char(10) // 'non-corrupt version or correct'
     &            // char(10) // 'PCF reference to required file'
     &            // char(10) // 'and rerun PGE.'

                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               Else


C---Write HDF attribute "CoreMetadata.0" to HDF product and close file
C---with single call to M-API routine CPMFIL.
                  NumHandles = 2

                  HDFAttNms(INVENTORYMETADATA) = MECS_CORE
                  HDFAttNms(ARCHIVEMETADATA)   = MECS_ARCHIVE

                  rtn = CPMFIL(modfil_mod06cd,MET_Handles,HDFAttNms,NumHandles)

                  If (RTN.EQ.mapiok) Then
                     msgbuf = '... Closed output HDF file after ECS metadata update.'
                     CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS, msgbuf,FUNCNAME)
                  Else
                     error_flag = .true.
                     msgbuf = 'CPMFIL returned error closing Cloud Product (MOD06_L2) '
     2                    // 'after ECS metadata update.'
     3                    // char(10) // 'Operator Action: Check system'
     4                    // char(10) // 'disk resources; If adequate,'
     5                    // char(10) // 'contact SDST'
                     CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf,FUNCNAME)
                  EndIf

               EndIf
            EndIf

         Else
            msgbuf =
     1      'DayNightFlag is "Night".  MOD_PR06CD exits before '
     2      // char(10) // 'processing any data and without '
     3      // 'updating the MOD06 product file in any way.'

            CALL MODIS_SMF_SETDYNAMICMSG( MODIS_M_GENERIC,msgbuf,FUNCNAME )
         EndIf

      EndIf


C.....Added by Liqun Ma 1998/03/11
C.....Check process error_flag.  Write terminating message and exit.
      If (error_flag) Then
         msgbuf = 'MODIS process MOD_PR06CD failed, exiting code 1'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,msgbuf,FUNCNAME)

         Call exit(1)
      Else
         msgbuf = 'MODIS process MOD_PR06CD completed successfully, exiting code 0.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)

         Call exit(0)
      EndIf

      END
