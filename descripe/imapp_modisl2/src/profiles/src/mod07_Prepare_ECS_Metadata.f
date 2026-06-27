      SUBROUTINE MOD07_PREPARE_ECS_METADATA( PCBAD, PCN4S, PCNND,
     &  PCNNT, PCNNG, PCNNI, PCNNL, PCNNW, PCNNS, PCNNV, PCNNR,
     &  PCNNC, MAX_SOL, MIN_SOL, MET_GROUPS )

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     This program is designed to drive the Core, Archive and
C     Inventory Metadata subroutines which place the values
C     into storage until they are attached to the product file
C     using mapi CPMFIL.
C
C !INPUT PARAMETERS:
C     PCBAD      Percentage of bad pixels
C     PCN4S      Percentage of pixels in category 4 per granule
C     PCNND      Percentage of day processed pixels per granule
C     PCNNT      Percentage of night processed pixels per granule
C     PCNNG      Percentage of sunglint found pixles per granule
C     PCNNI      Percentage of snow/ice processed pixels per granule
C     PCNNL      Percentage of land processed pixels per granule
C     PCNNW      Percentage of water processed pixels per granule
C     PCNNS      Percentage of shadow pixels per granule
C     PCNNV      Percentage of thin cirrus (vis) pixels per granule
C     PCNNR      Percentage of thin cirrus (ir) pixels per granule
C     PNCNC      Percentage of non-cloud obstruction pixels per granule
C     MAX_SOL    Maximum solar zenith angle for this granule
C     MIN_SOL    Minimum solar zenith angle for this granule
C
C !OUTPUT PARAMETERS:
C     MET_GROUPS Array containing the names of the "MASTER" groups
C                as defined in the MCF.
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C---------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_MODIS_39500.f'
      include 'Atmos_ECSMET.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'mod07_meta.inc'
      INCLUDE 'mod07_pcfnum.inc'

c ... Arguments

      REAL pcbad, pcn4s, pcnnd, pcnnt, pcnng, pcnni, pcnnl,
     &  pcnnw, pcnns, pcnnv, pcnnr, pcnnc, max_sol, min_sol
      CHARACTER*(PGSd_MET_GROUP_NAME_L) 
     &  MET_Groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER i,rtn,LRN_Prod_MCF

      CHARACTER*60 Auto_QA_Flag(No_ParaName),Auto_QA_Expl(No_ParaName),
     &  Value_Para_Name(No_ParaName)

      REAL PSAValue(No_PSA),Value_ArchivePSA_SC(MAX_NUM_ARCHIVE_PSA_SC)

      INTEGER QA_Miss(No_ParaName)

c ... External functions

      INTEGER Set_ArchiveMetadata, Set_CoreMetadata
      EXTERNAL Set_ArchiveMetadata, Set_CoreMetadata

      integer  Set_PSA, pgs_pc_getconfigdata
      character*4 pcf_satid
      character*255 doi_char
      integer LUN_Sat_Instrument
      parameter (LUN_Sat_Instrument = 800510)
      character*512 msgbuf


      integer num_args
      integer FlagRA
      character FlagBuff*10
      num_args = command_argument_count()
      
      if(num_args == 1) then
         call get_command_argument(1,FlagBuff)
         read (FlagBuff,*) FlagRA
      else
      !This is the default value
         FlagRA = 0
      endif

c ... Initializate of the output holder of the ODL metadata

      do i = 1, PGSd_MET_NUM_OF_GROUPS
         MET_Groups(i) = ' '
      end do

c ... Place metadata values in correct variables

      PSAValue(1)  = 100. - pcbad
      PSAValue(2)  = pcn4s
      PSAValue(3)  = pcnnd
      PSAValue(4)  = pcnnt
      PSAValue(5)  = pcnng
      PSAValue(6)  = pcnni
      PSAValue(7)  = pcnnl
      PSAValue(8)  = pcnnw
      PSAValue(9)  = pcnns
      PSAValue(10) = pcnnv
      PSAValue(11) = pcnnr
      PSAValue(12) = pcnnc
      PSAValue(13) = max_sol
      PSAValue(14) = min_sol
      
c ... Place metadata measured parameters in correct locations

      Value_Para_Name(1) = 'Water_Vapor_Infrared'

c ... Save percent missing for each parameter

      do i = 1 , No_ParaName
        QA_Miss(i) = nint(pcbad)
      enddo

c ... Save automatic QA flag

      do i = 1 , No_ParaName
        if (pcbad .gt. 90.) then
          Auto_QA_Flag(i) =  'Failed'
        else
          Auto_QA_Flag(i) =  'Passed'
        endif
        Auto_QA_Expl(i) =  'Passed: >10% usable; Failed: <10% usable'
      enddo

c ... Set our archive PSA code version number

      Value_ArchivePSA_SC(1) = mod07_ver(1)
      Value_ArchivePSA_SC(2) = mod07_ver(2)
      Value_ArchivePSA_SC(3) = mod07_ver(3)

C ... Assign LRN of MCF

      LRN_Prod_MCF = mcf_pcfnum

C---------------------------------------------------------------------
C     Set CoreMetadata
C---------------------------------------------------------------------

      rtn = -1

c rhucek 05/24/99:  Using latest version Set_CoreMetadata calling interface 
c     rtn = Set_CoreMetadata(LRN_Prod_MCF,ExtGeoPntr_Flag,
c    &      No_PCF_Attr,LRN_PCF_Attr,Name_PCF_Attr,
c    &      No_Input_LRN,Input_LRN,Input_Vrsn,
c    &      No_ParaName,Value_Para_Name,QA_Miss,Auto_QA_Flag,
c    &      Auto_QA_Expl,No_PSA,PSAName,PSAValue,MET_Groups)

      if( FlagRA == 1) then
         rtn = Set_CoreMetadata(LRN_Prod_MCF,ExtGeoPntr_Flag,
     &        No_Input_LRN,Input_LRN_RA,Input_Vrsn,
     &        No_ParaName,Value_Para_Name,QA_Miss,Auto_QA_Flag,
     &        Auto_QA_Expl,No_PSA,PSAName,PSAValue,
     &        LRN_OF_NUM_INV_RP_PAIRS,MET_Groups)
      else
         rtn = Set_CoreMetadata(LRN_Prod_MCF,ExtGeoPntr_Flag,
     &        No_Input_LRN,Input_LRN,Input_Vrsn,
     &        No_ParaName,Value_Para_Name,QA_Miss,Auto_QA_Flag,
     &        Auto_QA_Expl,No_PSA,PSAName,PSAValue,
     &        LRN_OF_NUM_INV_RP_PAIRS,MET_Groups)
      endif

      if ( rtn .eq. mapiok ) then
         call message('mod07_Prepare_ECS_Metadata',
     &  'Successfully set ECS Core (Inventory) metadata',
     &  0, 3 )
      else
        call message( 'mod07_Prepare_ECS_Metadata',
     &  'Error setting ECS Core (Inventory) metadata [OPERATOR ACTION: Notify SDST]',
     &  0, 1 )
      endif

C-----------------------------------------------------------------------
C  Set the doi Attributes (PSAs).
C-----------------------------------------------------------------------

c           Get satellite instrument name.
      rtn = pgs_pc_getconfigdata(LUN_Sat_Instrument,pcf_satid)
      if (rtn .ne. 0) then
         call message( 'mod07_prepare_ecs_metadata',
     &        'Error reading instrument name from pcf LUN 800510.' //
     &        ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &        0, 1 )
      endif

      if (pcf_satid .eq. 'AM1M') then
         doi_char = '10.5067/MODIS/MOD07_L2.006'
      else
         doi_char = '10.5067/MODIS/MYD07_L2.006'
      endif

      rtn = Set_PSA(MET_Groups, 'identifier_product_doi', No_PSA+1, doi_char)
      If (rtn .eq. FAIL) Then
         
         msgbuf = 'Set_PSA detected error setting DOI PSAs.'
     1        // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2        // char(10) // 'messages originating in function Set_InvPSA_Atmos.'
         
         
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,'mod07_Prepare_ECS_Metadata')
      EndIf

      rtn = Set_PSA(MET_Groups, 'identifier_product_doi_authority', No_PSA+2, 'http://dx.doi.org')
      If (rtn .eq. FAIL) Then
         
         msgbuf = 'Set_PSA detected error setting DOI PSAs.'
     1        // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2        // char(10) // 'messages originating in function Set_InvPSA_Atmos.'
         
         
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,'mod07_Prepare_ECS_Metadata')
      EndIf

C---------------------------------------------------------------------
C      Set Archive Metadata
C
C Note:  Set_ArchiveMetadata uses the same three arrays (No_Input_LRN,
C        Input_LRN and Input_Vrsn) to identify input file entries in 
C        the PCF as does Set_CoreMetadata.  Unlike Set_CoreMetadata, 
C        it retrieves and saves the file name, rather than the "UR", 
C        (Universal Reference) to the HDF product metadata. 
C---------------------------------------------------------------------

C rhucek 08/06/02:  Replaced NUM_MODIS_InputFiles, LRN_MODIS_InputFiles
C                   and VRSN_MODIS_InputFiles with No_Input_LRN, 
C                   Input_LRN and Input_Vrsn.  Declarations and 
C                   definitions of the former were removed; the latter
C                   are defined in mod07_meta.inc.

       rtn = -1

      if( FlagRA == 1) then
         rtn = Set_ArchiveMetadata(LRN_NUM_Of_ArchiveRP_Pairs,
     &        No_Input_LRN,
     &        Input_LRN_RA,
     &        Input_Vrsn,
     &        NUM_ArchivePSA_SC,
     &        Name_ArchivePSA_SC,
     &        Value_ArchivePSA_SC,
     &        MET_Groups)
      else
         rtn = Set_ArchiveMetadata(LRN_NUM_Of_ArchiveRP_Pairs,
     &        No_Input_LRN,
     &        Input_LRN,
     &        Input_Vrsn,
     &        NUM_ArchivePSA_SC,
     &        Name_ArchivePSA_SC,
     &        Value_ArchivePSA_SC,
     &        MET_Groups)
      endif

      if ( rtn .eq. mapiok ) then
         call message('mod07_Prepare_ECS_Metadata',
     &  'Successfully set ECS Archive metadata',
     &  0, 3 )
      else
        call message( 'mod07_Prepare_ECS_Metadata',
     &  'Error setting ECS Archive metadata [OPERATOR ACTION: Notify SDST]',
     &  0, 1 )
      endif

      END
