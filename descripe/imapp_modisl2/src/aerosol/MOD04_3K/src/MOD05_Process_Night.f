      Subroutine MOD05_Process_Night()

C-----------------------------------------------------------------------
C !F77
C
C !Description:  Open MOD05 nighttime product file and insert ECS
C                metadata.
C
C !Input PARAMETERs:
c    None
C
C !Output PARAMETERs:
C    None
C
c !REVISION HISTORY:
c 10/07/1999 fhliang
c removed multiple specification of data type; put non-DATA specification
c statement preceding DATA statements.
c
c 01/28/1998 fhliang
c fixed prolog.
c
c !TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C
C !REFERENCES and CREDITS:
C
C
C !DESIGN NOTES:
C
C
C Externals:
C   Named Constants:   MAPIOK                   (mapi.inc)
c                      MAX_ECS_NAME_L           (mapi.inc)
c                      MODFILLEN                (mapi.inc)
c                      PGS_S_SUCCESS            (PGS_SMF.f)
c                      MODIS_F_GENERIC          (PGS_MODIS_39500.f)
c                      MODIS_A_GENERIC          (PGS_MODIS_39500.f)
c                      PGSd_PC_FILE_PATH_MAX    (PGS_PC.f)
c                      PGSd_MET_GROUP_NAME_L    (PGS_MET.f, in mapi.inc)
c                      PGSd_MET_NUM_OF_GROUPS   (PGS_MET.f, in mapi.inc)
c
C
C Internals:
C
C   Functions and Subroutines:
C          MODIS_SMF_SETDYNAMICMSG
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT None

      INCLUDE 'mod05.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_SMF.f'


c-----Declare PARAMETERs
      CHARACTER*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'MOD05_Process_Night')

      INTEGER    No_PSA
      PARAMETER (No_PSA = 15)


c-----Declare local variables
      CHARACTER*25   msg25
      CHARACTER*1024 msgbuf
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) HDF_FILENAME
      CHARACTER*(PGSd_MET_GROUP_NAME_L) groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(MAX_ECS_NAME_L-1)   HDFAttNms(PGSd_MET_NUM_OF_GROUPS)

      INTEGER file_version,i,MODFIL_MOD05(MODFILLEN),
     1        NumHandles,fbyte,lbyte,rtn,rtn2

      INTEGER pgs_pc_getreference,String_Loc

      REAL QA_Metadata_MOD05(No_PSA)


c-----------------------------------------------------------------------
c Begin executable code
c-----------------------------------------------------------------------

      CALL MODIS_SMF_SETDYNAMICMSG(MODIS_A_GENERIC,
     1     'Setting Metadata for MOD05 night time process',FUNCNAME)

      FILE_VERSION = 1
      HDF_FILENAME = ' '
      write(msg25,'(I25)') LRN_MOD05
      rtn2 = String_Loc(msg25,fbyte,lbyte)

c-----retrieve name of MOD05 product file from PCF
      rtn = pgs_pc_getreference(LRN_MOD05,FILE_VERSION,HDF_FILENAME)

      If (rtn .ne. PGS_S_SUCCESS) Then
         msgbuf =
     1   'pgs_pc_getreference unable to retrieve MOD05 product file '
     2   // char(10) // 'name from PCF LUN ' // msg25(fbyte:lbyte)
     3   // char(10) // 'Operatior Action: Check system enviroment, PCF and '
     4   // char(10) // 'SDPTK configuration. If a fault is identified, correct '
     5   // char(10) // 'and rerun PGE. Otherwise, notify SDST'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf,FUNCNAME)
      Else

c--------Open mod05 product file
         rtn = OPMFIL(HDF_FILENAME,'a',modfil_mod05)

         If (rtn .ne. MAPIOK) Then
         msgbuf =
     1   'MAPI function opmfil unable to open MOD05 product file '
     2   // char(10) // 'on PCF LUN ' // msg25(fbyte:lbyte)
     3   // char(10) // 'Operatior Action: Check system enviroment, PCF and '
     4   // char(10) // 'SDPTK configuration. If a fault is identified, correct '
     5   // char(10) // 'and rerun PGE. Otherwise, notify SDST'

            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf,
     1                                   FUNCNAME)
         Else
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,
     1      '... Successfully opened mod05 HDF product file"',FUNCNAME)

            Do I=1,No_PSA
               QA_Metadata_MOD05(I)=0.0
            End Do

c-----------Write ECS metadata to product file
            CALL MOD05_Night_ECSMeta_Data(QA_Metadata_MOD05,groups,
     1                                    HDFAttNms,NumHandles)

c-----------Close MOD05 nighttime product file
            rtn = CPMFIL(modfil_mod05,groups,HDFAttNms,NumHandles)

            If (rtn .EQ. MAPIOK) THEN
               CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,
     1         '... Successfully closed "mod05.hdf"',FUNCNAME)
            Else

         msgbuf =
     1   'MAPI function opmfil unable to close MOD05 product file '
     2   // char(10) // 'on PCF LUN ' // msg25(fbyte:lbyte)
     3   // char(10) // 'Operatior Action: Check related error messages '
     4   // char(10) // 'possibly due to  the failure of setting metadata and etc '
     5   // char(10) // 'and rerun PGE. Otherwise, notify SDST'

               CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     1         '... MAPI CPMFIL FOR "mod05.hdf"',FUNCNAME)

            End If

         End If  ! Opened product file

      End If  ! Retrieved product file name


      Return

      End


C*********************************************************************


       SUBROUTINE MOD05_Night_ECSMeta_Data(QA_Metadata_MOD05,
     * MET_groups,HDFAttNms,NumHandles)

C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:   This subroutine sets up varaibles for ecs core meta data
C                 and product specific attributes
C
C !INPUT PARAMETERS:
C
C
C !OUTPUT PARAMETERS:
C              groups
C              HDFAttNms
C              NumHandles
C
C !REVISION HISTORY:
c 01/29/98 fhliang
c added NCSA acknowledgement;
c fixed prolog.
c
C*/  Modified by JC Guu  11/05/96
C*/  One variables numScanline and pixPerScanline
C*/  are added to the argument list.
C
C*/  Modified by Dr. Allen Chu  7/24/97
C*/  For V2 software, including inventory and archive metadata
C*/  Details can be seen in MODIS atmosphere QA plan
C
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !END
C-----------------------------------------------------------------------
C
      IMPLICIT  NONE

      INCLUDE 'mod05.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'mapi.inc'

c declare PARAMETERs
      CHARACTER*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'MOD05_Night_ecsMeta_data')

C------- V2 setup (9/25/97)
      INTEGER No_PCF_Attr,No_PSA,No_ParaName
C rhucek 12/18/97: changed No_PCF_Attr from 6 to 7
c rhucek 06/04/99:  No_PCF_Attr no longer used
      PARAMETER (No_ParaName = 2, No_PSA = 15)
c     PARAMETER (No_PCF_Attr = 6, No_ParaName = 2, No_PSA = 15)
      CHARACTER*1024 msgbuf

      CHARACTER*(PGSd_MET_GROUP_NAME_L) groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(MAX_ECS_NAME_L-1) HDFAttNms(PGSd_MET_NUM_OF_GROUPS)
c rhucek 06/04/99:  Name_PCF_Attr(No_PCF_Attr) no longer used
C      CHARACTER*60 Name_PCF_Attr(No_PCF_Attr),PSAName(No_PSA),
      CHARACTER*60 PSAName(No_PSA),
     *             Value_Para_Name(No_ParaName),
     *             Auto_QA_Flag(No_ParaName),Auto_QA_Expl(No_ParaName)

c rhucek 06/04/99:  LRN_PCF_Attr(No_PCF_Attr) no longer used
      INTEGER I,rtn,No_Input_LRN,INVENTORYMETADATA,
     *        ARCHIVEMETADATA,
C     *        LRN_Prod,LRN_PCF_Attr(No_PCF_Attr),ARCHIVEMETADATA,
     *        QA_Miss(2),
     *        NumHandles

C\* No_Input_LRN=3 for at nighttime
      PARAMETER (INVENTORYMETADATA=2,ARCHIVEMETADATA=3,No_Input_LRN=2)
      INTEGER Input_LRN(No_Input_LRN),Input_Vrsn(No_Input_LRN)

      INTEGER Set_CoreMetadata,Set_ArchiveMetadata
      INTEGER Vrsn_No

      REAL QA_Metadata_MOD05(No_PSA)

      LOGICAL ExtGeoPntr_Flag
      REAL PSAValue(No_PSA)
C------- Archive metadata
      INTEGER LRN_NUM_Of_ArchiveRP_Pairs
      PARAMETER (LRN_NUM_Of_ArchiveRP_Pairs = 412250)

      INTEGER MAX_NUM_MODIS_Inputs
      PARAMETER (MAX_NUM_MODIS_Inputs = 30)

      INTEGER NUM_ArchivePSA_SC
C      PARAMETER (NUM_ArchivePSA_SC = 8)
      PARAMETER (NUM_ArchivePSA_SC = 0)

      INTEGER NUM_MODIS_InputFiles
      PARAMETER (NUM_MODIS_InputFiles = 2)

C      CHARACTER*100 Name_ArchivePSA_SC(NUM_ArchivePSA_SC)
      CHARACTER*100 Name_ArchivePSA_SC(1)

      INTEGER LRN_MODIS_InputFiles(MAX_NUM_MODIS_Inputs)
      INTEGER VRSN_MODIS_InputFiles(MAX_NUM_MODIS_Inputs)

C      real Value_ArchivePSA_SC(NUM_ArchivePSA_SC)
      real Value_ArchivePSA_SC(1)

c-----lipo 06/03/99:  Define LRN_OF_NUM_INV_RP_PAIRS
c-----functionally replaces variables and arrays: No_PCF_Attr, Name_PCF_Attr, LRN_PCF_Attr
      integer    LRN_OF_NUM_INV_RP_PAIRS
      parameter (LRN_OF_NUM_INV_RP_PAIRS = 412050)

C------- ECS core metadata and PSA
      CHARACTER*(PGSd_MET_GROUP_NAME_L)
     *          MET_Groups(PGSd_MET_NUM_OF_GROUPS)
      data Auto_QA_Flag / 'Failed','Passed' /,
     *     Auto_QA_Expl / 'No test done','Tests Passed'/

c------lipo 06/03/99:  No_PCF_Attr, Name_PCF_Attr, LRN_PCF_Attr nolonger used;
c------replaced functionally by use of PARAMETER LRN_OF_NUM_INV_RP_PAIRS
CC-----rhucek 12/18/97: added ASSOCIATED* attributes
C      data Name_PCF_Attr / 'REPROCESSINGACTUAL',
C     1                     'REPROCESSINGPLANNED',
C     2                     'LOCALVERSIONID',
C     3                     'PGEVERSION',
C     4                     'ASSOCIATEDPLATFORMSHORTNAME.1',
C     5                     'ASSOCIATEDINSTRUMENTSHORTNAME.1',
C     6                     'ASSOCIATEDSENSORSHORTNAME.1' /
C
c     data Name_PCF_Attr / 'REPROCESSINGACTUAL', 'REPROCESSINGPLANNED',
c    *                     'LOCALVERSIONID', 'PGEVERSION',
c    *                     'PLATFORMSHORTNAME.1',
c    *                     'INSTRUMENTSHORTNAME.1'/
      data ExtGeoPntr_Flag /.false. /
C rhucek 12/18/97: added 412056 to specification of array LRN_PCF_Attr
C      data LRN_PCF_Attr /412050,412051,412052,412053,412054,412055,
C     *                   412056/

      data PSAName / 'SuccessfulRetrievalPct_NIR',
     *               'SuccessfulRetrievalPct_IR',
     *               'LowConfidentClearPct',
     *               'DayProcessedPct','NightProcessedPct',
     *               'SunglintProcessedPct',
     *               'Snow_IceSurfaceProcessedPct',
     *               'LandProcessedPct','WaterProcessedPct',
     *               'ShadowFoundPct',
     *               'ThinCirrusSolarFoundPct','ThinCirrusIR_FoundPct',
     *               'NonCloudObstructionFoundPct',
     *               'MaxSolarZenithAngle','MinSolarZenithAngle' /
      data PSAValue / 50.00,50.00,60.00,90.00,10.00,10.00,10.00,
     *                65.00,35.00,10.00,20.00,20.00,10.00,30.00,10.00/
      data QA_Miss / 20,20 /

      data Name_ArchivePSA_SC /' '/
      data Value_ArchivePSA_SC / 0.00/

C      data Name_ArchivePSA_SC /'VeryGoodQualityDataPct_NIR',
C     +                         'GoodQualityDataPct_NIR',
C     +                         'MarginalQualityDataPct_NIR',
C     +                         'BadQualityDataPct_NIR',
C     +                         'VeryGoodQualityDataPct_IR',
C     +                         'GoodQualityDataPct_IR',
C     +                         'MarginalQualityDataPct_IR',
C     +                         'BadQualityDataPct_IR'/

C      data Value_ArchivePSA_SC / 0.00,100.00,0.00,0.00,
C     +                          64.65,3.80,15.96,17.0/

      Value_Para_Name(1) = 'Water_Vapor_Near_Infrared'
      Value_Para_Name(2) = 'Water_Vapor_Infrared'
C------- V2 setup (9/25/97)

      NumHandles = 2

      HDFAttNms(INVENTORYMETADATA) = MECS_CORE
      HDFAttNms(ARCHIVEMETADATA) = MECS_ARCHIVE

      do 10 I = 1, No_Input_LRN
         Input_Vrsn(i) = 1
   10 continue
      do 20 I = 1, PGSd_MET_NUM_OF_GROUPS
         groups(I)= ' '
   20 continue
      do 30 I = 1, PGSd_MET_NUM_OF_GROUPS
         MET_Groups(I)= ' '
   30 continue

      Input_LRN(1) = LRN_Geo
      Input_LRN(2) = LRN_MOD07

      do 40 I = 1, NUM_MODIS_InputFiles
         LRN_MODIS_InputFiles(I)=Input_LRN(I)
         VRSN_MODIS_InputFiles(I)=1
   40 continue

c\* Set PSA metadata values

      PSAValue(1)=0.0
      QA_Miss(1)=100

      Vrsn_No=1
      Call  Additional_Attr_MOD05(LRN_CldMsk,LRN_MOD07,
     &         Vrsn_No,PSAValue,QA_Miss(2),No_PSA,PSAName)

      Auto_QA_Flag(1)='Failed'
      Auto_QA_Expl(1)='NoSolarBandWaterRetrieval'

c\* Adding ECS core metadata to HDF file

C------lipo 06/03/99: Using latest version Set_CoreMetadata calling interface(V2.1)
C      rtn = Set_CoreMetadata(LRN_MCF,ExtGeoPntr_Flag,
C     *      No_PCF_Attr,LRN_PCF_Attr,Name_PCF_Attr,
C     *      No_Input_LRN,Input_LRN,Input_Vrsn,
C     *      No_ParaName,Value_Para_Name,QA_Miss,Auto_QA_Flag,
C     *      Auto_QA_Expl,No_PSA,PSAName,PSAValue,MET_Groups)

      rtn = Set_CoreMetadata(LRN_MCF,ExtGeoPntr_Flag,
     *      No_Input_LRN,Input_LRN,Input_Vrsn,
     *      No_ParaName,Value_Para_Name,QA_Miss,Auto_QA_Flag,
     *      Auto_QA_Expl,No_PSA,PSAName,PSAValue,
     *      LRN_OF_NUM_INV_RP_PAIRS,MET_Groups)


      If (rtn.ne.MAPIOK) then
         msgbuf =
     1   'Set_CoreMetadata failed'
     2   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)
      Else
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,
     2        'Set_CoreMetadata passed without error', FUNCNAME)
      End If

c\* Adding archive metadata to HDF file

       rtn = -1
       rtn = Set_ArchiveMetadata(LRN_NUM_Of_ArchiveRP_Pairs,
     +                           NUM_MODIS_InputFiles,
     +                           LRN_MODIS_InputFiles,
     +                           VRSN_MODIS_InputFiles,
     +                           NUM_ArchivePSA_SC,
     +                           Name_ArchivePSA_SC,
     +                           Value_ArchivePSA_SC,
     +                           MET_Groups)

      If (rtn.ne.MAPIOK) then
         msgbuf =
     1   'Set_ArchiveMetadata failed'
     2   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

      Else
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,
     2        'Set_ArchiveMetadata passed without error', FUNCNAME)
      End If


      RETURN

      END
