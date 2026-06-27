       SUBROUTINE MOD04_Metadata(QA_Metadata_MOD04,Modfil_L1B,RTN_NCEP,
     * MET_groups,HDFAttNms,NumHandles)
C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:   This subroutine sets up varaibles for ecs core meta data
C                 and product specific attributes
C
C !INPUT PARAMETERS:
C         modfil_L1B
C
C !OUTPUT PARAMETERS:
C              groups
C              HDFAttNms
C              NumHandles
C
C !REVISION HISTORY:
c 10/08/1999 fhliang
c put non-DATA specification statement preceding DATA statements.
c
c 01/29/1998 fhliang
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
C !END
C-----------------------------------------------------------------------
C
      IMPLICIT  NONE
      SAVE

c rhucek 01/16/98: SAVE statement added to retain values of named
c                COMMON /ONE/ (in mod04.inc).

      INCLUDE   'mod04.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'mapi.inc'

c rhucek 01/16/98: Commented 5 used include files
c     INCLUDE 'PGS_SMF.f'
c     INCLUDE 'PGS_PC.f'
c     INCLUDE 'PGS_PC_9.f'
c     INCLUDE 'PGS_IO.f'
c     INCLUDE 'PGS_IO_1.f'

      INTEGER RTN_NCEP

C------- Archive metadata
      integer LRN_NUM_Of_ArchiveRP_Pairs
      parameter (LRN_NUM_Of_ArchiveRP_Pairs = 407250)

      integer MAX_NUM_MODIS_Inputs
      parameter (MAX_NUM_MODIS_Inputs = 30)

      integer NUM_ArchivePSA_SC
      parameter (NUM_ArchivePSA_SC = 8)

      integer NUM_MODIS_InputFiles
      parameter (NUM_MODIS_InputFiles = 5)

      integer ODL_IN_MEMORY
      PARAMETER (ODL_IN_MEMORY=1)
      
      integer        FAIL,    SUCCEED
      PARAMETER     (FAIL=-1, SUCCEED=0)

      character*100 Name_ArchivePSA_SC(NUM_ArchivePSA_SC)

      integer LRN_MODIS_InputFiles(MAX_NUM_MODIS_Inputs)
      integer VRSN_MODIS_InputFiles(MAX_NUM_MODIS_Inputs)

      real Value_ArchivePSA_SC(NUM_ArchivePSA_SC)

      integer pgs_met_init,pgs_met_write

C------- V2 setup (9/25/97)
      INTEGER No_PSA,No_ParaName,Vrsn_No

      PARAMETER (No_ParaName = 2, No_PSA = 15)

      CHARACTER*(PGSd_MET_GROUP_NAME_L) groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(MAX_ECS_NAME_L-1) HDFAttNms(PGSd_MET_NUM_OF_GROUPS)

      CHARACTER*60 PSAName(No_PSA),
     *             Value_Para_Name(No_ParaName),
     *             Auto_QA_Flag(No_ParaName),Auto_QA_Expl(No_ParaName)

      INTEGER I,rtn,No_Input_LRN,INVENTORYMETADATA,
     *        ARCHIVEMETADATA,
     *        QA_Miss(2),
     *        NUM_of_HDFAttrNms,
     *        Modfil_L1B(MODFILLEN),NumHandles

      PARAMETER (INVENTORYMETADATA=2,ARCHIVEMETADATA=3,No_Input_LRN=18)
      INTEGER Input_LRN(No_Input_LRN),Input_Vrsn(No_Input_LRN)

c rhucek 12/18/97: added Set_CoreMetadata_QC declaration
      INTEGER Set_CoreMetadata,Set_CoreMetadata_QC,Set_ArchiveMetadata
      REAL QA_Metadata_MOD04(No_PSA)

      LOGICAL ExtGeoPntr_Flag
      REAL PSAValue(No_PSA)

C------- ECS core metadata and PSA
      CHARACTER*(PGSd_MET_GROUP_NAME_L)
     *          MET_Groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*1024 msgbuf

      integer  Set_PSA, pgs_pc_getconfigdata
      character*4 pcf_satid
      character*26 doi_char
      integer LUN_Sat_Instrument
      parameter (LUN_Sat_Instrument = 800510)

      integer    LRN_OF_NUM_INV_RP_PAIRS
      parameter (LRN_OF_NUM_INV_RP_PAIRS = 407050)

c    all data statements......
      data Auto_QA_Flag / 'Failed','Passed' /,
     *     Auto_QA_Expl / 'No land found','No test done'/

      data ExtGeoPntr_Flag /.false. /

      data PSAName / 'SuccessfulRetrievalPct_Land',
     *               'SuccessfulRetrievalPct_Ocean',
     *               'LowConfidentClearPct',
     *               'DayProcessedPct',
     *               'NightProcessedPct',
     *               'SunglintProcessedPct',
     *               'Snow_IceSurfaceProcessedPct',
     *               'LandProcessedPct',
     *               'WaterProcessedPct',
     *               'ShadowFoundPct',
     *               'ThinCirrusSolarFoundPct',
     *               'ThinCirrusIR_FoundPct',
     *               'NonCloudObstructionFoundPct',
     *               'MaxSolarZenithAngle',
     *               'MinSolarZenithAngle' /

      data PSAValue / 50.00,50.00,60.00,90.00,10.00,10.00,10.00,
     *                65.00,35.00,10.00,20.00,20.00,10.00,30.00,10.00/

      data QA_Miss / 20,20 /

      data Name_ArchivePSA_SC /'VeryGoodQualityDataPct_Land',
     +                         'GoodQualityDataPct_Land',
     +                         'MarginalQualityDataPct_Land',
     +                         'BadQualityDataPct_Land',
     +                         'VeryGoodQualityDataPct_Ocean',
     +                         'GoodQualityDataPct_Ocean',
     +                         'MarginalQualityDataPct_Ocean',
     +                         'BadQualityDataPct_Ocean'/

      data Value_ArchivePSA_SC / 0.00,100.00,0.00,0.00,
     +                           64.65,3.80,15.96,17.0/

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

c
      Value_Para_Name(1) = 'Corrected_Optical_Depth_Land'
      Value_Para_Name(2) = 'Effective_Optical_Depth_Average_Ocean'
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
      Input_LRN(2) = LRN_L1B_250
      Input_LRN(3) = LRN_L1B_500
      if( FlagRA .eq. 1) then
         Input_LRN(4) = LRN_L1B_1km_RA
      else
         Input_LRN(4) = LRN_L1B_1km
      endif
      Input_LRN(5) = LRN_CldMsk
      Input_LRN(6) = S1
      Input_LRN(7) = S2
      Input_LRN(8) = S3
      Input_LRN(9) = L1
      Input_LRN(10) = L2
      Input_LRN(11) = L3
      Input_LRN(12) = LRN_LUT466
      Input_LRN(13) = LRN_LUT553
      Input_LRN(14) = LRN_LUT644
      Input_LRN(15) = LRN_LUT213
      Input_LRN(16) = LRN_LUTMAP
      IF(RTN_NCEP.EQ.0) THEN
        Input_LRN(17) = LRN_WISC_ANC_met   
      ENDIF
       Input_LRN(18) = INSCI
      do 40 I = 1, NUM_MODIS_InputFiles
         LRN_MODIS_InputFiles(I)=Input_LRN(I)
         VRSN_MODIS_InputFiles(I)=1
   40 continue

      Vrsn_No=1

      Call  Additional_Attr_MOD04(LRN_CldMsk,Vrsn_No,
     &                      PSAValue,No_PSA,PSAName)

C
C QA_Miss over land is defined as 100 - successful rate of
C retrievals
C

      PSAValue(1)=QA_Metadata_MOD04(1)
      QA_Miss(1)=100-NINT(PSAValue(1))


      IF(PSAValue(8).GT.0.0) THEN
        IF(PSAValue(1).GT.25.0) THEN
          Auto_QA_Flag(1)='Passed'
          Auto_QA_Expl(1)='SuccessfulRetrievalPct>25%'
        ELSE
          Auto_QA_Flag(1)='Suspect'
          Auto_QA_Expl(1)='FurtherInvestigationNeeded'
        ENDIF
      ENDIF
C
C QA_Miss over Ocean is defined as 100 - successful rate of
C retrievals
C

      PSAValue(2)=QA_Metadata_MOD04(2)
      QA_Miss(2)=100-NINT(PSAValue(2))


      IF(PSAValue(9).GT.0.0) THEN
        IF(PSAValue(2).GT.25.0) THEN
          Auto_QA_Flag(2)='Passed'
          Auto_QA_Expl(2)='SuccessfulRetrievalPct>25%'
        ELSE
          Auto_QA_Flag(2)='Suspect'
          Auto_QA_Expl(2)='FurtherInvestigationNeeded'
        ENDIF
      ENDIF

C  Adding ECS core metadata to HDF file

      rtn = Set_CoreMetadata(LRN_MCF,ExtGeoPntr_Flag,
     *      No_Input_LRN,Input_LRN,Input_Vrsn,
     *      No_ParaName,Value_Para_Name,QA_Miss,Auto_QA_Flag,
     *      Auto_QA_Expl,No_PSA,PSAName,PSAValue,
     *      LRN_OF_NUM_INV_RP_PAIRS,MET_Groups)

      if (rtn.ne.mapiok) then
         msgbuf =
     1   'Set_CoreMetadata failed in MOD04 HDF file'
     2   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,
     &   'Metadata_MOD04_V2')

      else
         call modis_smf_setdynamicmsg(MODIS_S_SUCCESS,
     2   'Set_CoreMetadata passed without error',
     3   'Metadata_MOD04_V2')
      end if

C-----------------------------------------------------------------------
C  Set the doi Attributes (PSAs).
C-----------------------------------------------------------------------

c           Get satellite instrument name.
      rtn = pgs_pc_getconfigdata(LUN_Sat_Instrument,pcf_satid)
      if (rtn .ne. 0) then
         call message( 'Metadata_MOD04_V2',
     &        'Error reading instrument name from pcf LUN 800510.' //
     &        ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &        0, 1 )
      endif

      if (pcf_satid .eq. 'AM1M') then
         doi_char = '10.5067/MODIS/MOD04_L2.006'
      else
         doi_char = '10.5067/MODIS/MYD04_L2.006'
      endif

      rtn = -1
      rtn = Set_PSA(MET_Groups, 'identifier_product_doi', No_PSA+1, doi_char)
      If (rtn .eq. FAIL) Then
         
         msgbuf = 'Set_PSA detected error setting DOI PSAs.'
     1        // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2        // char(10) // 'messages originating in function Set_InvPSA_Atmos.'
         
         
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,'Metadata_MOD04_V2')
      EndIf

      rtn = -1
      rtn = Set_PSA(MET_Groups, 'identifier_product_doi_authority', No_PSA+2, 'http://dx.doi.org')
      If (rtn .eq. FAIL) Then
         
         msgbuf = 'Set_PSA detected error setting DOI PSAs.'
     1        // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2        // char(10) // 'messages originating in function Set_InvPSA_Atmos.'
         
         
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,'Metadata_MOD04_V2')
      EndIf


C  Adding archive metadata to HDF file

       rtn = -1

       rtn = Set_ArchiveMetadata(LRN_NUM_Of_ArchiveRP_Pairs,
     +                           NUM_MODIS_InputFiles,
     +                           LRN_MODIS_InputFiles,
     +                           VRSN_MODIS_InputFiles,
     +                           NUM_ArchivePSA_SC,
     +                           Name_ArchivePSA_SC,
     +                           Value_ArchivePSA_SC,
     +                           MET_Groups)

      if (rtn.ne.mapiok) then
         msgbuf =
     1   'Set_ArchiveMetadata failed in MOD04 HDF file'
     2   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,
     &   'Metadata_MOD04_V2')

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,
     2   'Set_ArchiveMetadata passed without error',
     3   'Metadata_MOD04_V2')
      end if

C  Adding metadata to mod05.qc file

      rtn = Set_CoreMetadata_QC(LRN_MCFQC,LRN_QCMET)

      if (rtn.ne.mapiok) then
         msgbuf =
     1   'Set_CoreMetadata failed in MOD04 QC file'
     2   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,
     &   'Metadata_MOD04_V2')

      else
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,
     2   'Create metadata file to non-HDF (mod04_QC) without error',
     3   'Metadata_MOD04_V2')
      end if

      RETURN
      END




C*********************************************************************
       SUBROUTINE Additional_Attr_MOD04(LRN_CldMsk,Vrsn_No,
     &   PSAValue,No_PSA,PSAName)

C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:   This subroutine reads additional attributes from
C                 cloud mask file
C
C !INPUT PARAMETERS:
C
C     LRN_CldMsk                      Cloud mask logical unit number
C     Vrsn_No                         Version number (=1)
C     No_PSA                          Number of PSAs (=15)
C     PSAName                         Names of PSAs (see below)
C
C !OUTPUT PARAMETERS:
C
C     PSAValue contains values of the followings:
C
C     LowConfidentClearPct            Percent of cloudy pixels found
C     DayProcessedPct                 Percent of daytime coverage
C     NightProcessedPct               Percent of nighttime coverage
C     SunglintProcessedPct            Percent of sun glint coverage
C     Snow_IceSurfaceProcessedPct     Percent of snow/ice coverage
C     LandProcessedPct                Percent of land coverage
C     WaterProcessedPct               Percent of water coverage
C     ShadowFoundPct                  Percent of shadow found
C     ThinCirrusSolarFoundPct         Percent of thin cirrus solar detecetd
C     ThinCirrusIR_FoundPct           Percent of thin cirrus IR detected
C     NonCloudObstructionPct          Percent of non-obstruction found
C     MaxSolarZenithAngle             Maximum solar zenith angle
C     MinSolarZenithAngle             Minimum solar zenith angle
C
C !REVISION HISTORY:
c 01/28/98 fhliang
c fixed prolog.
c
C   Written by Dr. Allen Chu  12/23/97
C
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !END
C-----------------------------------------------------------------------

      IMPLICIT  NONE
      SAVE

      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'mapi.inc'


      INTEGER No_PSA,rtn,I
      INTEGER LRN_CldMsk, Vrsn_No, pgs_met_getpcattr_s

      Real PSAValue(*)

      CHARACTER*8 char_buf
      CHARACTER*256 HDFAttrName,ECSParmName
      CHARACTER*1024 msgbuf

      CHARACTER*(*) PSAName(*)

C
C Extract additional attritubes from cloud mask file
C

      HDFAttrName = 'CoreMetadata.0'

      DO I=3,No_PSA

      ECSParmName = PSAName(I)

      rtn=pgs_met_getpcattr_s(LRN_CldMsk,Vrsn_No,HDFAttrName,
     1                        ECSParmName,char_buf)
      IF(rtn.ne.mapiok) THEN
        msgbuf =
     1  'Retrieve PSA attributes from cloud mask failed'
        CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,
     &  'Metadata_MOD04_V2')
      ENDIF

      read(char_buf,'(f8.2)') PSAValue(I)

      ENDDO

      RETURN
      END
