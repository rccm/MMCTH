      Integer Function Set_ArchMet_MOD06( MET_Handles,
     1                                    NUM_Of_MODIS_InputFiles,
     2                                    LUN_MODIS_InputFiles,
     3                                    VRSN_MODIS_InputFiles,
     4                                    NUM_Of_ArchivePSA_SC,
     5                                    Name_ArchivePSA_SC,
     6                                    Value_ArchivePSA_SC )

      implicit none

      include 'mod06_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Set MOD06_L2 (MODIS L2 cloud product) ECS "Archive" metadata objects
C   to internal memory.  In the context used here, set means to associate
C   a value with an ODL (Object Description Language) object defined in
C   the MOD06_L2 metadata configuration file (MCF).  Function calls to the
C   ECS Science Data Processing Toolkit (SDPTK) are used for this purpose.
C   No metadata fields are actually written to the product file by
C   Set_ArchMet_MOD06 - only the association of values to ODL objects is
C   made at this time.  Actual insertion of Archive metadata into the HDF
C   product file takes place later with a call to the SDPTK routine
C   PGS_MET_Write.  This call is made from outside Set_ArchMet_MOD06, after
C   all other ECS metadata relevant to the process have been set.
C
C   The specific metadata fields set by Set_ArchMet_MOD06, their
C   origin, and the value of the Metadata Configuration File (MCF)
C   "Data Location" parameter are listed below.
C      GEO = Geolocation product granule
C      PCF = Process Control File
C      RP  = PCF runtime Parameter (RP)
C      PGE = Science code
C
C                                          Source of    Data Location
C         ECS Metadata Objects               Value      Value in MCF
C   -------------------------------        ---------    -------------
C      1    AlgorithmPackageAcceptanceDate    PCF (RP)      PGE
C      2    AlgorithmPackageMaturityCode      PCF (RP)      PGE
C      3    AlgorithmPackageName              PCF (RP)      PGE
C      4    AlgorithmPackageVersion           PCF (RP)      PGE
C      5    LongName (optional)               PCF (RP)      PGE
C      6    InstrumentName                    PCF (RP)      PGE
C      7    ExclusionGringFlag                GEO           PGE
C      8    GringPointLatitude                GEO           PGE
C      9    GringPointLongitude               GEO           PGE
C     10    GringPointSequenceNo              GEO           PGE
C
C     Product Specific Attributes (PSAs) Objects
C   -----------------------------------------------
C     11    LocalInputGranuleID        MODIS input products PGE
C
C     12    Algorithm_Version_Cloud_Top_Property_IR
C                                             PCF (RP)      PGE

C     13    Algorithm_Version_Cloud_Phase_IR  PCF (RP)      PGE
C     14    Algorithm_Version_Cloud_Property_VIS
C                                             PCF (RP)      PGE
C
C
C !INPUT PARAMETERS:
C
C  integer NUM_Of_MODIS_InputFiles
C                Variable containing the number of MODIS input
C                products used by process.
C
C  integer LUN_MODIS_InputFiles(*)
C                Array containing MODIS input product LUNs used by
C                process (See also VRSN_MODIS_InputFiles below).
C
C                A one-to-one association between the elements of
C                the arrays LUN_MODIS_InputFiles and
C                VRSN_MODIS_InputFiles is assumed.
C
C  integer VRSN_MODIS_InputFiles(*)
C                Array containing file version numbers of MODIS input
C                products used by currrent process (See
C                LUN_MODIS_InputFiles above).
C
C                A one-to-one association between the elements of
C                the arrays LUN_MODIS_InputFiles and
C                VRSN_MODIS_InputFiles is assumed.
C
C  integer NUM_Of_ArchivePSA_SC
C                Variable containing the number of archive Product
C                Specific Attributes (PSA) with granule dependent
C                values determined by the science code (SC) to be set
C                by process.  If none, pass NUM_Of_ArchivePSA_SC=0.
C
C                NOTE1: Currently MOD06_L2 contains no ArchivePSA_SC.
C
C                NOTE2: Some archive PSAs are known prior to running
C                       PGE06 and defined as RPs in the PCF.  (See
C                       LUN_OF_NUM_ARCH_RP_MOD06 in mod06_ECSMET.inc).
C                       Still others are granule dependent, but are set
C                       internally by Set_ArchMet_MOD06.  The
C                       LOCALINPUTGRANULEID is one.
C
C                NOTE3: Although Archive PSAs must appear as ODL data
C                       objects in the product MCF, they are not stored
C                       in the ECS "inventory" database for search and
C                       query purposes as are the AdditionalAttribute
C                       fields (also referred to as PSAs) of ECS
C                       Inventory metadata.
C
C
C  character*(*) Name_ArchivePSA_SC(*)
C                Array containing the names of archive PSAs with granule
C                dependent values determined by the SC (See
C                NUM_Of_ArchivePSA_SC above and Value_ArchivePSA_SC below).
C
C                A one-to-one association between the elements of
C                arrays Name_ArchivePSA_SC and Value_ArchivePSA_SC
C                is assumed.
C
C  Real Value_ArchivePSA_SC(*)
C                Array containing the values of archive PSAs with granule
C                dependent values determined by the SC (See
C                NUM_Of_ArchivePSA_SC and Name_ArchivePSA_SC above).
C
C                A one-to-one association between the elements of
C                arrays Name_ArchivePSA_SC and Value_ArchivePSA_SC
C                is assumed.
C
C                NOTE: Passed values are restricted to range from -9999.99
C                to +99999.99 and are stored in the product file (by
C                Set_ArchMet_MOD06 as formatted, F8.2 ASCII character strings.
C
C
C  character MET_Handles(20)
C                 Array containing the names of the "MASTER" groups as
C                 defined in the MCF (There may be up to 20 "MASTER"
C                 groups in the MCF.).  "MET_Handles" must be the same
C                 array as initialized previously in a call to function
C                 Update_InvMeta_MOD06.  Update_InvMeta_MOD06 returns the
C                 initialized array back to the calling routine.
C
C                 In FORTRAN, element 1 of the MET_Handles array refers
C                 to the MCF file, element 2 to ECS inventory metadata,
C                 and element 3 to archive metadata.
C
C
C !OUTPUT PARAMETERS:  None
C
C
C !REVISION HISTORY:
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by Rich Hucek
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        Set_ArchPSA_MOD06         (science code)
C        Set_GRing                 (science code)
C        Set_LclInputGranID_MOD06  (science code)
C        Set_RP_Data_Atmos         (science code)
C        String_Loc                (science code)
C
C    Named Constant:
C        ARCHIVEMETADATA           (mod06_ECSMET.inc)
C        FAIL                      (mod06_ECSMET.inc)
C        LUN_OF_NUM_ARCH_RP_MOD06  (mod06_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSd_MET_GROUP_NAME_L     (PGS_MET.f: included in mapi.inc)
C        SUCCEED                   (mod06_ECSMET.inc)
C
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Set_ArchMet_MOD06')


c-----Declaration of function arguments
      character*(*) MET_Handles(*),Name_ArchivePSA_SC(*)

      integer NUM_Of_MODIS_InputFiles
      integer LUN_MODIS_InputFiles(*)
      integer VRSN_MODIS_InputFiles(*)
      integer NUM_Of_ArchivePSA_SC

      real Value_ArchivePSA_SC(*)


c-----Declaration of local variables
      character*25   msg25_MET
      character*1024 msgbuf

      integer string_loc,
     1        Set_RP_Data_Atmos, Set_LclInputGranID_MOD06,
     2        Set_GRing, Set_ArchPSA_MOD06

      integer fbyte_MET,
     1        lbyte_MET,
     2        rtn,rtn_loc,rtn_MET

      logical error_flag


C------------------------
C Initialization
C------------------------

      Set_ArchMet_MOD06 = FAIL
      error_flag = .false.


c-------------------------------------------------------------------------------
c Input Argument Checks
c 1 - Check string len of MET_Handles array
c
c-------------------------------------------------------------------------------
      rtn_MET = LEN( MET_Handles(1) )

      write(msg25_MET,'(I25)')  rtn_MET
      rtn_loc = string_loc(msg25_MET,fbyte_MET,lbyte_MET)

      If (rtn_MET .lt. PGSd_MET_GROUP_NAME_L) Then
         msgbuf =
     1   'MET_Handles array element size (' // msg25_MET(fbyte_MET:lbyte_MET)
     2   // ' characters) declared < PGSd_MET_GROUP_NAME_L.'
     3   // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         return
      End If


C----------------------------------------------------------------------
C Set ECS metadata objects that are specified as name/value pairs in
C the USER DEFINED RUNTIME PARAMETERS section of the PCF.
C----------------------------------------------------------------------

      rtn = Set_RP_Data_Atmos( Met_Handles,
     1                         LUN_OF_NUM_ARCH_RP_MOD06,
     2                         ARCHIVEMETADATA )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_RP_Data_Atmos detected error setting MOD06 archive RP data.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_RP_Data_Atmos.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


C----------------------------------------------------------------------
C Set ECS GPolygon group attributes.
C----------------------------------------------------------------------

      rtn = Set_GRing(MET_Handles, ARCHIVEMETADATA)

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_GRing detected error setting ECS GPolygon group attributes.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_GRing.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C-----------------------------------------------------------------
C Update MOD06 LocalInputGranuleID
C-----------------------------------------------------------------

      rtn = Set_LclInputGranID_MOD06( MET_Handles,
     1                                NUM_Of_MODIS_InputFiles,
     2                                LUN_MODIS_InputFiles,
     3                                VRSN_MODIS_InputFiles )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_LclInputGranID_MOD06 detected error incrementing ECS '
     1            // char(10) // 'LocalInputGranuleID metadata attribute.'
     2            // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3            // char(10) // 'messages originating from call to routine '
     4            // 'Set_LclInputGranID_MOD06.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


C----------------------------------------------------------------------
C Set archive Product Specific Attributes (PSAs).  PSA values are
C stored in the product as formatted, F8.2 ASCII character strings
C
C----------------------------------------------------------------------

      rtn = Set_ArchPSA_MOD06( MET_Handles,
     1                         NUM_Of_ArchivePSA_SC,
     2                         Name_ArchivePSA_SC,
     3                         Value_ArchivePSA_SC )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_ArchPSA_MOD06 detected error setting archive PSA attribute.'
     1            // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     2            // char(10) // 'messages originating from call to function '
     3            // 'Set_ArchPSA_MOD06.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


      If (.not.error_flag) Set_ArchMet_MOD06 = SUCCEED

      Return

      End



c---------------------------------------------------------------------------------------
      subroutine Copy_LclInputGranID_MOD06( NUM_MODIS_FileNames,
     1                                      MODIS_FileNames,
     2                                      rtn_flag)

      implicit none

      include 'mod06_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_MET_13.f'
      include 'PGS_MODIS_39500.f'
      include 'PGS_SMF.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Copy existing LOCALINPUTGRANULEID archive metadata 
C                from the MOD06_L2 product file to array storage.  
C                Return array to calling routine.
C
C !INPUT PARAMETERS:
C
C !OUTPUT PARAMETERS:       
C  Integer NUM_MODIS_FileNames    Variable containing the number of
C                                 MODIS file names in MOD06
C                                 LOCALINPUTGRANULEID.
C
C  Character*(*) MODIS_FileNames(*)
C                                 Array containing the set of MODIS 
C                                 file names retrieved from MOD06
C                                 LOCALINPUTGRANULEID.
C
C  Integer rtn_flag               Subroutine return flag 0 = SUCCESS
C                                                       -1 = FAIL
C 
C
C !REVISION HISTORY:
C
C
C !REFERENCES AND CREDITS:
C
C    Written by        Liqun Ma  Feb. 1998 
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    lma@ltpmail.gsfc.nasa.gov
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
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error 
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        string_loc                (science code)
C
C
C    Named Constant:
C        FAIL                      (mod06_ECSMET.inc)
C        LUN_MOD06                 (mod06_ECSMET.inc)
C        MAX_NUM_LCLPNTRS          (mod06_ECSMET.inc)
C        MCORE_LOCALINPUTGRANULEID (mapi.inc)
C        MECS_ARCHIVE              (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSMET_E_NULL_PARAMETER   (PGS_MET_13.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (mod06_ECSMET.inc)
C        VRSN_MOD06                (mod06_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)   FUNCNAME
      PARAMETER      (FUNCNAME = 'Copy_LclInputGranID_MOD06')


C input argument declarations
      character*(*)  MODIS_FileNames(*)
      integer        NUM_MODIS_FileNames, rtn_flag


C local variable declarations
      character*25     msg25_LUN,msg25_VRSN
      character*255    HDF_FileAttrName, ECS_AttrName
      character*(1024) msgbuf

      integer String_Loc, PGS_MET_GetPCAttr_s
      integer rtn_loc, rtn,
     1        fbyte,fbyte_LUN,fbyte_VRSN,
     2        lbyte,lbyte_LUN,lbyte_VRSN

      logical error_flag


c initialization
      NUM_MODIS_FileNames =  0
      MODIS_FileNames(1)  =  PGSd_MET_STR_END 
      rtn_flag            =  FAIL
      error_flag          = .FALSE.
      
c-----------------------------------------------------------------------
c Set up status message variables 
c-----------------------------------------------------------------------

      write(msg25_LUN,'(I25)') LUN_MOD06
      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

      write(msg25_VRSN,'(I25)') VRSN_MOD06
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)


c-------------------------------------------------------------------------------
c Copy existing MOD06 LocalInputGranuleID to internal memory, and count the
c number of file names found. 
c-------------------------------------------------------------------------------

      HDF_FileAttrName  = MECS_ARCHIVE
      ECS_AttrName      = MCORE_LOCALINPUTGRANULEID 
      rtn_loc           = String_Loc(ECS_AttrName,fbyte,lbyte)

c-----read LocalInputGranuleID
      rtn = PGS_MET_GetPCAttr_s( LUN_MOD06, 
     1                           VRSN_MOD06, 
     2                           HDF_FileAttrName, 
     3                           ECS_AttrName, 
     4                           MODIS_FileNames )  

      If (rtn .EQ. PGS_S_SUCCESS) Then
         NUM_MODIS_FileNames = 1
   
c--------count the number of MODIS file names retrieved from MOD06 LocalInputGranuleID 
         Do WHILE ( (MODIS_FileNames(NUM_MODIS_FileNames) .NE. PGSd_MET_STR_END) .AND.
     1              (NUM_MODIS_FileNames .LE. MAX_NUM_LCLPNTRS) ) 
   
            NUM_MODIS_FileNames = NUM_MODIS_FileNames + 1
         End Do
   
         NUM_MODIS_FileNames = NUM_MODIS_FileNames - 1


c-----LocalInputGranuleID not yet set 
      Else If ( rtn.eq.PGSMET_E_NULL_PARAMETER .OR. rtn.eq.PGSMET_E_PARAMETER_NOT_SET) Then
         error_flag = .TRUE.

         msgbuf = 
     1   MECS_ARCHIVE // ' in MOD06 product file, yet the required object ' 
     2   // MCORE_LOCALINPUTGRANULEID // ' is not set - improbable PGE06 processing sequence. '
     3   // char(10) // 'Set ' // MCORE_LOCALINPUTGRANULEID // ' now and complete current process. '
     4   // char(10) // 'Operator Action:  Notify SDST.' 
     
         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 


c-----Failed to read LocalInputGranuleID 
      Else  
         error_flag = .TRUE.

         msgbuf = 'PGS_MET_GetPCAttr_s detected error reading ECS attribute '
     1   // ECS_Attrname(fbyte:lbyte)
     2   // char(10) // 'from MODIS cloud product (LUN = ' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     3   // ' and VRSN number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '). '
     5   // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     6   // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7   // char(10) // 'fault is identified, stage correct PCF/input file and '
     8   // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

      End If   !  INPUTPOINTER read check

c-----set return flag
      If (.NOT. error_flag) rtn_flag = SUCCEED 

      Return

      End



c---------------------------------------------------------------------------------------
      Integer Function Set_ArchPSA_MOD06( MET_Handles,
     1                                    NUM_NewArchPSA_SC,
     2                                    Name_NewArchPSA_SC,
     3                                    Value_NewArchPSA_SC )

      implicit none

      include 'mod06_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_MET_13.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C
C !INPUT PARAMETERS:
C
C  integer NUM_NewArchPSA_SC
C                Variable containing the number of archive Product
C                Specific Attributes (PSA) with granule dependent
C                values determined by the science code (SC) to be set
C                by process.  If none, pass NUM_NewArchPSA_SC=0.
C
C                NOTE1: Currently MOD06_L2 contains no ArchivePSA_SC.
C
C                NOTE2: Some archive PSAs are known prior to running
C                       PGE06 and defined as RPs in the PCF.  (See
C                       LUN_OF_NUM_ARCH_RP_MOD06 in mod06_ECSMET.inc).
C                       Still others are granule dependent, but are set
C                       internally by Set_ArchMet_MOD06.  The
C                       LOCALINPUTGRANULEID is one.
C
C                NOTE3: Although Archive PSAs must appear as ODL data
C                       objects in the product MCF, they are not stored
C                       in the ECS "inventory" database for search and
C                       query purposes as are the AdditionalAttribute
C                       fields (also referred to as PSAs) of ECS
C                       Inventory metadata.
C
C
C  character*(*) Name_NewArchPSA_SC(*)
C                Array containing the names of archive PSAs with granule
C                dependent values determined by the SC (See
C                NUM_NewArchPSA_SC above and Value_NewArchPSA_SC below).
C
C                A one-to-one association between the elements of
C                arrays Name_NewArchPSA_SC and Value_NewArchPSA_SC
C                is assumed.
C
C  Real Value_NewArchPSA_SC(*)
C                Array containing the values of archive PSAs with granule
C                dependent values determined by the SC (See
C                NUM_NewArchPSA_SC and Name_NewArchPSA_SC above).
C
C                A one-to-one association between the elements of
C                arrays Name_NewArchPSA_SC and Value_NewArchPSA_SC
C                is assumed.
C
C                NOTE: Passed values are restricted to range from -9999.99
C                to +99999.99 and are stored in the product file (by
C                Set_ArchMet_MOD06 as formatted, F8.2 ASCII character strings.
C
C
C  character MET_Handles(20)
C                 Array containing the names of the "MASTER" groups as
C                 defined in the MCF (There may be up to 20 "MASTER"
C                 groups in the MCF.).  "MET_Handles" must be the same
C                 array as initialized previously in a call to function
C                 Update_InvMeta_MOD06.  Update_InvMeta_MOD06 returns the
C                 initialized array back to the calling routine.
C
C                 In FORTRAN, element 1 of the MET_Handles array refers
C                 to the MCF file, element 2 to ECS inventory metadata,
C                 and element 3 to archive metadata.
C
C
C !OUTPUT PARAMETERS:  None
C
C
C !REVISION HISTORY:
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by Rich Hucek
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        PGS_MET_GetPCAttr_s       (libPGSTK.a)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        Query_HDF_Attr            (science code)
C        String_Loc                (science code)
C
C    Named Constant:
C        ARCHIVEMETADATA           (mod06_ECSMET.inc)
C        FAIL                      (mod06_ECSMET.inc)
C        LUN_MOD06                 (mod06_ECSMET.inc)
C        LUN_MCF_MOD06             (mod06_ECSMET.inc)
C        MECS_ARCHIVE              (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        MODIS_M_GENERIC           (PGS_MODIS_39500.f)
C        MODIS_W_GENERIC           (PGS_MODIS_39500.f)
C        NUM_ARCH_PSA_SC_MOD06     (mod06_ECSMET.inc)
C        PGSMET_E_NULL_PARAMETER   (PGS_MET_13.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (mod06_ECSMET.inc)
C        VRSN_MOD06                (mod06_ECSMET.inc)
C
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Set_ArchPSA_MOD06')

      real           REL_EQUALITY_EPS
      parameter     (REL_EQUALITY_EPS=0.000001)


c-----Declaration of function arguments
      character*(*) MET_Handles(*),Name_NewArchPSA_SC(*)
      integer       NUM_NewArchPSA_SC
      real          Value_NewArchPSA_SC(*)


c-----Declaration of local variables
      Character*8    msg8
      Character*25   msg25_a, msg25_b, msg25_c, msg25_d, msg25_LUN,
     1               msg25_MCF, msg25_VRSN
      Character*150  PSA_Name,PSA_Value
      Character*1024 msgbuf

      integer fbyte,fbyte_a,fbyte_b,fbyte_c,fbyte_d,fbyte_LUN,fbyte_MCF,fbyte_VRSN,
     1        lbyte,lbyte_a,lbyte_b,lbyte_c,lbyte_d,lbyte_LUN,lbyte_MCF,lbyte_VRSN,
     2        rtn1,rtn2,rtn_loc,rtn_status,i0,i1,i2,i3

      integer      NUM_Blank_PSA
      real         Value

      integer pgs_met_setattr_s,
     1        pgs_met_getpcattr_s,
     2        string_loc

      integer Lcl_NUM_NewPSA

      logical Name_Is_Valid,error_flag, arch_exists



C------------------------
C Initialization
C------------------------

      Set_ArchPSA_MOD06  = SUCCEED
      arch_exists        = .FALSE.
      rtn_status         = 0
      error_flag         = .FALSE.

      write(msg25_a,'(i25)')   NUM_NewArchPSA_SC
      write(msg25_b,'(i25)')   NUM_ARCH_PSA_SC_MOD06
      write(msg25_MCF,'(i25)') LUN_MCF_MOD06

      rtn_loc = String_Loc(msg25_a,fbyte_a,lbyte_a)
      rtn_loc = String_Loc(msg25_b,fbyte_b,lbyte_b)
      rtn_loc = String_Loc(msg25_MCF,fbyte_MCF,lbyte_MCF)


c-----------------------------------------------------------------------
c Perform input argument checks
c-----------------------------------------------------------------------

c-----check for positive number of PSAs
      If (NUM_NewArchPSA_SC .LE. 0) Then

         msgbuf =
     1   'The number of new PSAs (' // msg25_a(fbyte_a:lbyte_a) // ') on input '
     2   // 'argument list is < 1.'
     3   // char(10) // 'No new archive PSAs to be set by current process, but '
     4   // 'existing ones will be copied.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)
      End If


c-----check for blank PSAs
      NUM_Blank_PSA = 0

      Do 50 i0 = 1,NUM_NewArchPSA_SC
         If (Name_NewArchPSA_SC(i0) .eq. BLANK) NUM_Blank_PSA = NUM_Blank_PSA + 1
50    Continue

      If (NUM_Blank_PSA .GT. 0) Then
         write(msg25_c,'(i25)') NUM_Blank_PSA
         rtn_loc = String_Loc(msg25_c,fbyte_c,lbyte_c)

         msgbuf =
     1   msg25_c(fbyte_c:lbyte_c) // ' PSAs on input argument list have blank names.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)
      End If


c-----------------------------------------------------------------------
c Check if ECS Archive metadata is in product output file
c-----------------------------------------------------------------------

      call Query_HDF_Attr( LUN_MOD06,
     1                     VRSN_MOD06,
     2                     MECS_ARCHIVE,
     3                     rtn_status,
     4                     arch_exists )


c--------------------------------------------------------------------------------
c Return if Error, or if there are no new PSAs to set and no existing PSA to copy
c--------------------------------------------------------------------------------

      If (rtn_status .eq. FAIL) then
         Set_ArchPSA_MOD06 = FAIL

         msgbuf =
     1   'Query_HDF_Attr detected error while searching MOD06 product for '
     2   // 'presence of ' // MECS_ARCHIVE // '.'
     3   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4   // char(10) // 'messages generated by routine Query_HDF_Attr.'

         Return
      Else If ( (NUM_NewArchPSA_SC .le. 0) .and. (.NOT. arch_exists) ) Then
         Return
      End If


c--------------------------------------------------------------------------------
c If archive metadata are in file, loop over PSAs in MODIS Cloud Product and copy
c existing PSAs to internal memory.
c--------------------------------------------------------------------------------

      If (arch_exists) Then

         Do 100 i1 = 1, NUM_ARCH_PSA_SC_MOD06
            PSA_Name  = Name_ArchPSA_SC_MOD06(i1)
            rtn_loc   = string_loc(PSA_Name,fbyte,lbyte)

            rtn1 = PGS_MET_GetPCAttr_s( LUN_MOD06,
     1                                  VRSN_MOD06,
     2                                  MECS_ARCHIVE,
     3                                  PSA_Name,
     4                                  PSA_Value )


c-----------PSA read successful. Now copy PSA to memory
            If (rtn1 .EQ. PGS_S_SUCCESS) Then

               rtn2 = PGS_MET_SetAttr_s( MET_Handles,
     1                                   PSA_NAME,
     2                                   PSA_Value )

               If (rtn2 .EQ. FAIL) Then
                  error_flag = .TRUE.

                  msgbuf =
     1            'PGS_MET_SetAttr_s detected error setting MOD06 Archive PSA '
     2            // PSA_Name(fbyte:lbyte)
     3            // char(10) // 'Operator Action:  Check for correct MCF on LUN '
     4            // msg25_MCF(fbyte_MCF:lbyte_MCF) // '.  If wrong, '
     5            // char(10) // 'stage correct MCF and rerun PGE.  Otherwise, notify SDST.'



                  Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If


c-----------PSA is not in archive metadata.
            Else If (rtn1 .EQ. PGSMET_E_NULL_PARAMETER) Then

                  msgbuf =
     1            'Query for Archive PSA ' // PSA_Name(fbyte:lbyte) // ' in ' // MECS_ARCHIVE 
     2            // ' unsuccessful - '
     3            // char(10) // PSA_Name(fbyte:lbyte) // ' not yet in MOD06 ' // MECS_ARCHIVE // '! '
     4            // char(10) // 'Operator Action:  Disregard prior SDPTK error message '
     5            // 'PGSMET_E_NULL_PARAMETER.'

                  Call MODIS_SMF_SetDynamicMsg(MODIS_M_GENERIC,msgbuf,FUNCNAME)

            Else
               error_flag = .TRUE.

               write(msg25_LUN, '(I25)') LUN_MOD06
               write(msg25_VRSN,'(I25)') VRSN_MOD06

               rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)
               rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)

               msgbuf = 'PGS_MET_GetPCAttr_d unable to retrieve archive '// PSA_Name(fbyte:lbyte)
     1         // char(10) // 'from MOD06 product on LUN = '// msg25_LUN(fbyte_LUN:lbyte_LUN)
     2         // ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     4         // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     5         // char(10) // 'fault is identified, stage correct PCF/input file and '
     6         // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            End IF   ! End PSA existence check

100      Continue   ! End loop on MOD06 PSAs

      End If   ! check on arch_exists


c-------------------------------------------------------------------------------
c Return if Error during PSA copy
c-------------------------------------------------------------------------------

      If (error_flag) then
         Set_ArchPSA_MOD06 = FAIL
         Return
      End If


C-----------------------------------------------------------------
c Loop on list of new PSAs.  Set PSAs with valid names and values
c PSA values are stored in the product as formatted, F8.2 ASCII
c character strings
C-----------------------------------------------------------------

c-----test for too many archive PSAs
      Lcl_NUM_NewPSA = NUM_NewArchPSA_SC

      If (NUM_NewArchPSA_SC .gt. NUM_ARCH_PSA_SC_MOD06) Then
         error_flag = .true.
         Lcl_NUM_NewPSA = NUM_ARCH_PSA_SC_MOD06

         msgbuf =
     1   'Number of new archive PSAs (=' // msg25_a(fbyte_a:lbyte_a) // ') exceeds total '
     2   // char(10) // 'number in MOD06 product (' // msg25_b(fbyte_b:lbyte_b) // ').'
     3   // 'Only first ' // msg25_b(fbyte_b:lbyte_b) // ' archive PSAs to be examined.'
     4   // char(10) // 'Operator Action:  Notify SDST.'

         Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


      Do 200 i2 = 1, Lcl_NUM_NewPSA
         PSA_Name = Name_NewArchPSA_SC(i2)
         rtn_loc  = string_loc(PSA_Name,fbyte,lbyte)

c--------check validity of new PSA name
         Name_Is_Valid = .FALSE.

         Do 300 i3 = 1, NUM_ARCH_PSA_SC_MOD06
            If (Name_ArchPSA_SC_MOD06(i3) .eq. PSA_Name) Name_Is_Valid = .TRUE.
300      Continue


c--------invalid PSA name
         If (.NOT. Name_Is_Valid) Then
            error_flag = .TRUE.

            msgbuf = 'Archive PSA name ' // PSA_Name(fbyte:lbyte) // ' not on MOD06 metadata list.'
     1      // char(10) // PSA_Name(fbyte:lbyte) // ' Name/Value pair not set.'
     2      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------valid PSA name
         Else

c-----------check for "No Set" option
            If ( (abs((Value-NO_SET_FLAG)/NO_SET_FLAG) .gt. REL_EQUALITY_EPS) ) Then

c--------------PSA value out of range; cannot set Name/Value pair
               If (Abs(Value) .gt. MAX_ABS_VALUE_PSA) Then
                  error_flag = .TRUE.

                  write(msg25_d,'(E25.6)') Value
                  rtn_loc = String_Loc(msg25_d,fbyte_d,lbyte_d)

                  msgbuf =
     1            'Archive PSA value ' // msg25_d(fbyte_d:lbyte_d) // ' is out of bounds!'
     2            // char(10) // 'PSA ' // PSA_Name(fbyte:lbyte)
     3            // ' Name/Value pair not updated.'
     4            // char(10) // 'Operator Action:  Notify SDST.'


                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------------Set PSA Name/Value pair
               Else
                  write(msg8,'(f8.2)') Value

                  rtn2 = PGS_MET_SetAttr_s( MET_Handles,
     1                                      PSA_Name,
     2                                      Value )

                  If (rtn2 .EQ. FAIL) Then
                     error_flag = .TRUE.

                     msgbuf =
     1               'PGS_MET_SetAttr_s detected error setting archive PSA '
     2               // PSA_Name(fbyte:lbyte) // ' Name/Value pair.'
     3               // char(10) // 'Operator Action:  Check for correct MCF on LUN '
     4               // msg25_MCF(fbyte_MCF:lbyte_MCF) // '.  If wrong, '
     5               // char(10) // 'stage correct MCF and rerun PGE.  Otherwise, notify SDST.'

                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  Else
                     msgbuf =
     1               'Archive PSA ' // PSA_Name(fbyte:lbyte) // ' Name/Value pair successfully set.'

                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
                  End If

               End If   ! PSA valid range check

            End If   ! PSA "ok to set" test

         End If   ! PSA name check

200   Continue   ! Loop on Input PSAs


      If (.not.error_flag) Set_ArchPSA_MOD06 = SUCCEED

      Return

      End



c---------------------------------------------------------------------------------------
      Integer Function Set_LclInputGranID_MOD06( MET_Handles,
     +                                           NUM_MODIS_InputFiles,
     +                                           LUN_MODIS_InputFiles,
     +                                           VRSN_MODIS_InputFiles)

      implicit none

      include 'mod06_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Copy, increment and set an updated list of PGE06 MODIS 
C                input product file pointers.  The file pointers are 
C                names in MODIS convention that are stored as ECS Archive 
C                metadata (under attribute LOCALINPUTGRANULEID) in the 
C                MOD06 product file.  
C
C
C !INPUT PARAMETERS:
C  character*(*) MET_Handles(*)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           In FORTRAN, element 1 of the MET_Handles
C                           array is reserved for the MCF file.   ECS
C                           inventory metadata is referenced as as
C                           element 2 and archive metadata as element 3.
C
C  integer  NUM_MODIS_InputFiles
C                           Variable containing the number of MODIS input
C                           files used by current process.
C
C  integer  LUN_MODIS_InputFiles(*)
C                           Array containing PCF LUN numbers of the MODIS input
C                           files used by current process.
C                           (See also VRSN_MODIS_InputFiles below.)
C
C  integer  VRSN_MODIS_InputFiles(*)
C                           Array containing the file version numbers of MODIS
C                           input products used by current process.
C                           (See LUN_MODIS_InputFiles above.)
C
C                           A one-to-one correspondence between the elements of
C                           the arrays LUN_MODIS_InputFiles and
C                           VRSN_MODIS_InputFiles is assumed.
C
C
C !OUTPUT PARAMETERS:       None
C
C !REVISION HISTORY:
C
C
C !REFERENCES AND CREDITS:
C
C    Written by        Liqun Ma  Feb. 1998
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    lma@ltpmail.gsfc.nasa.gov
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_SetAttr_s           (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG     (science code)
C        Query_HDF_Attr              (science code)
C        Copy_LclInputGranID_MOD06   (science code)
C        Update_LclInputGranID_MOD06 (science code)
C
C
C    Named Constant:
C        ARCHIVEMETADATA            (mod06_ECSMET.inc)
C        FAIL                       (mod06_ECSMET.inc)
C        LUN_MOD06                  (mod06_ECSMET.inc)
C        MAX_NUM_LCLPNTRS           (mod06_ECSMET.inc)
C        MCORE_LOCALINPUTGRANULEID  (mapi.inc)
C        MECS_ARCHIVE               (mapi.inc)
C        MODIS_E_GENERIC            (PGS_MODIS_39500.f)
C        PGSd_PC_VALUE_LENGTH_MAX   (PGS_PC.f)
C        PGS_S_SUCCESS              (PGS_SMF.f)
C        PGSd_MET_STR_END           (PGS_MET.f, included in mapi.inc)
C        SUCCEED                    (mod06_ECSMET.inc)
C        VRSN_MOD06                 (mod06_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Set_LclInputGranID_MOD06')


C input argument declarations
      character*(*) MET_Handles(*)
      integer NUM_MODIS_InputFiles
      integer LUN_MODIS_InputFiles(*)
      integer VRSN_MODIS_InputFiles(*)


C local variable declarations
      character*255    HDF_FileAttrName, ECS_AttrName
      character*(1024) msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) LCLGRANID(MAX_NUM_LCLPNTRS) 
     1                                     /MAX_NUM_LCLPNTRS*PGSd_MET_STR_END/

      integer NUM_Old_LCLGRANID, NUM_Total_LCLGRANID
      integer fbyte,lbyte
      integer rtn, rtn_loc, rtn_status, rtn_status_query 
      integer PGS_MET_SetAttr_s, String_Loc

      logical archive_exists, error_flag


c-----Initialization
      Set_LclInputGranID_MOD06 = FAIL
      error_flag               = .FALSE.
      NUM_Old_LCLGRANID        = 0
      NUM_Total_LCLGRANID      = 0

c-----------------------------------------------------------------------
c Check if ECS archive metadata are in MOD06 product file
c Exit routine if Subroutine Query_HDF_Attr returns FAIL.
c-----------------------------------------------------------------------

      HDF_FileAttrName  = MECS_ARCHIVE

      call Query_HDF_Attr( LUN_MOD06,
     1                     VRSN_MOD06,
     2                     HDF_FileAttrName,
     3                     rtn_status_query,
     4                     archive_exists )

      If (rtn_status_query .eq. FAIL) then

         msgbuf = 
     1   'Query_HDF_Attr detected error searching MOD06 cloud product file '
     2   // char(10) // 'for ' // MECS_ARCHIVE // '.'
     3   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4   // char(10) // 'messages generated by routine Query_HDF_Attr.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Return
      End If


c-------------------------------------------------------------------------------
c If archive metadata are in file, copy existing LocalInputGranuleID to 
c LCLGRANID array.  
c-------------------------------------------------------------------------------

      If (archive_exists) Then

         Call Copy_LclInputGranID_MOD06( NUM_Old_LCLGRANID,
     1                                   LCLGRANID,
     2                                   rtn_status )

         If (rtn_status .eq. FAIL) Then
            error_flag = .TRUE.

            msgbuf =
     1      'Copy_LclInputGranID_MOD06 detected error copying LocalInputGranuleID '
     2      // char(10) // 'from MOD06 product file. '
     3      // char(10) // 'Operator Action:  Refer to prior low level error LogStatus '
     4      // char(10) // 'messages generated by routine Copy_LclInputGranID_MOD06.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If
      End If


c-------------------------------------------------------------------------------
c Aggregate new MODIS input file names to LCLGRANID array.  Track the total 
c number (this process + previous processes) of unique MODIS files used as 
c input to PGE06. 
c-------------------------------------------------------------------------------

c-----if unable to read existing LocalInputGranuleID, set only MODIS file
c-----names used by current process.

      If (error_flag) NUM_Old_LCLGRANID = 0

      Call Update_LclInputGranID_MOD06( NUM_MODIS_InputFiles,
     1                                  LUN_MODIS_InputFiles,
     2                                  VRSN_MODIS_InputFiles,
     3                                  NUM_Old_LCLGRANID,
     4                                  LCLGRANID,
     5                                  NUM_Total_LCLGRANID,
     6                                  rtn_status )

      If (rtn_status .eq. FAIL) Then
         error_flag = .TRUE.

         msgbuf = 
     1   'Update_LclInputGranID_MOD06 detected error setting LocalInputGranuleID '
     2   // char(10) // 'from MOD06 product file. '
     3   // char(10) // 'Operator Action:  Refer to prior low level error LogStatus '
     4   // char(10) // 'messages generated by routine Update_LclInputGranID_MOD06.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


c-------------------------------------------------------------------------------
c set LOCALINPUTGRANULEID
c-------------------------------------------------------------------------------

c-----If sufficient storage, add end-of-string marker to LCLGRANID array 
      If (NUM_Total_LCLGRANID .LT. MAX_NUM_LCLPNTRS) 
     1   LCLGRANID(NUM_Total_LCLGRANID + 1) = PGSd_MET_STR_END

      ECS_AttrName = MCORE_LOCALINPUTGRANULEID

      rtn = PGS_MET_SetAttr_s(MET_Handles(ARCHIVEMETADATA),ECS_AttrName,LCLGRANID)

      If (rtn .NE. PGS_S_SUCCESS) Then
          error_flag = .TRUE.
          rtn_loc = String_Loc(ECS_AttrName,fbyte,lbyte)

          msgbuf = 'PGS_MET_SetAttr_s detected error setting attribute '
     1    // ECS_AttrName(fbyte:lbyte)
     2    // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3    // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4    // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


      If (.NOT. error_flag) Set_LclInputGranID_MOD06 = SUCCEED

      Return

      End



c---------------------------------------------------------------------------------------
      Subroutine Update_LclInputGranID_MOD06( NUM_MODIS_InputFiles,
     1                                        LUN_MODIS_InputFiles,
     2                                        VRSN_MODIS_InputFiles,
     3                                        NUM_Old_LCLGRANID,
     4                                        LCLGRANID,
     5                                        NUM_Total_LCLGRANID,
     6                                        rtn_status)

      implicit none

      include 'mod06_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'
      include 'PGS_PC.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Increment the existing list of PGE06 MODIS input product
C                files names with additional names used by the current
C                process.  Only unique new are aggregated to the list
C                which is maintained and passed as a function array
C                argument.
C
C
C !INPUT PARAMETERS:
C  integer  NUM_MODIS_InputFiles
C                           Variable containing the number of MODIS product
C                           files used as input by current process.
C
C  integer  LUN_MODIS_InputFiles(*)
C                           Array containing PCF LUN numbers of MODIS product
C                           files used as input by current process.
C                           (See also VRSN_MODIS_InputFiles below.)
C
C  integer  VRSN_MODIS_InputFiles(*)
C                           Array containing the PCF version numbers of MODIS
C                           product files used as input by current process.
C                           (See LUN_MODIS_InputFiles above.)
C
C                           A one-to-one correspondence between the elements of
C                           the arrays LUN_MODIS_InputFiles and
C                           VRSN_MODIS_InputFiles is assumed.
C
C  integer  NUM_Old_LCLGRANID
C                           Variable containing the number of unique MODIS
C                           product files used as input to PGE06 prior to the
C                           start of the current process.
C
C
C
C !OUTPUT PARAMETERS:
C  NUM_Total_LCLGRANID      Total number of LocalGranuleIDs stored in LCLGRANID
C                           array.
C
C  integer rtn_status       Subroutine return flag: 0 = SUCCESS, -1 = FAIL.
C
C
C !INPUT/OUTPUT PARAMETERS:
C  character*(*) LCLGRANID(*)
C                           Array that on input contains the names of unique MODIS
C                           product files (in MODIS convention) used as input to PGE06
C                           prior to the start of the current process.  On output,
C                           the array is updated to include the names of the unique
C                           MODIS products used by the current process.
C
C
C !REVISION HISTORY:
C
C
C !REFERENCES AND CREDITS:
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_GetPCAttr_s       (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        Retrieve_LclGranID_Set    (Met_Common_Atmos.f)
C        string_loc                (science code)
C
C
C    Named Constant:
C        FAIL                       (mod06_ECSMET.inc)
C        MAX_NUM_LCLPNTRS           (mod06_ECSMET.inc)
C        MCORE_LOCALINPUTGRANULEID  (mapi.inc)
C        MODIS_E_GENERIC            (PGS_MODIS_39500.f)
C        MODIS_M_GENERIC            (PGS_MODIS_39500.f)
C        PGSd_PC_VALUE_LENGTH_MAX   (PGS_PC.f)
C        SUCCEED                    (mod06_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Update_LclInputGranID_MOD06')


C input argument declarations
      integer NUM_MODIS_InputFiles
      integer LUN_MODIS_InputFiles(*)
      integer VRSN_MODIS_InputFiles(*)
      integer NUM_Old_LCLGRANID
      character*(*) LCLGRANID(*)
      integer rtn_status


C local variable declarations
      character*25     msg25_2,msg25_3,msg25_4
      character*(1024) msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) New_LCLGRANID(MAX_NUM_LCLPNTRS)

      integer Array_Index
      integer i1, i2, rtn_flag, rtn_loc
      integer fbyte2,fbyte3,fbyte4,
     1        lbyte2,lbyte3,lbyte4
      integer NUM_Duplicates,NUM_New_LCLGRANID,NUM_Pntrs,
     1        NUM_Total_LCLGRANID,NUM_Unique
      integer String_Loc

      logical error_flag, In_LCLGRANID


c-----Initialization
      rtn_status          = FAIL
      error_flag          = .FALSE.
      NUM_Duplicates      = 0
      NUM_Unique          = 0
      NUM_New_LCLGRANID   = 0

      NUM_Pntrs           = NUM_MODIS_InputFiles
      If (NUM_MODIS_InputFiles .gt. MAX_NUM_LCLPNTRS) NUM_Pntrs = MAX_NUM_LCLPNTRS

      If (NUM_Old_LCLGRANID .gt. 0) Then
         NUM_Total_LCLGRANID = NUM_Old_LCLGRANID
         Array_Index         = NUM_Old_LCLGRANID
      Else
         NUM_Total_LCLGRANID = 0
         Array_Index         = 0
      End If


c-------------------------------------------------------------------------------
c Retrieve LocalGranuleIDs of MODIS input files.  Track the number of MODIS
c file names (i.e. LocalGranuleIDs) already in file, and total number of
c MODIS file names (old + new).
c-------------------------------------------------------------------------------

      If (NUM_MODIS_InputFiles .GT. 0) Then

c--------retrieve LocalInputGranuleID set for current process
         Call Retrieve_LclGranID_Set( NUM_Pntrs,
     1                                LUN_MODIS_InputFiles,
     2                                VRSN_MODIS_InputFiles,
     3                                NUM_New_LCLGRANID,
     4                                New_LCLGRANID,
     5                                rtn_flag )

         If (rtn_flag .EQ. SUCCEED) Then

            Do 100 i1 = 1, NUM_New_LCLGRANID
               In_LCLGRANID = .FALSE.

c--------------search against existing file name set for duplicate name
               Do 200 i2 = 1, Array_Index 

                  If (LCLGRANID(i2) .EQ. New_LCLGRANID(i1)) Then
                     In_LCLGRANID = .TRUE.
                     NUM_Duplicates   =  NUM_Duplicates + 1
                  End If

200            Continue


c--------------new name is unique
               If ( .NOT. In_LCLGRANID ) Then

c-----------------new name is non-blank - add to LOCALINPUTGRANULEID array if
c-----------------sufficient storage

                  If (New_LCLGRANID(i1) .ne. BLANK) Then
                     NUM_Total_LCLGRANID = NUM_Total_LCLGRANID + 1
                     NUM_Unique          = NUM_Unique + 1

                     If (NUM_Total_LCLGRANID .LE. MAX_NUM_LCLPNTRS) Then
                        Array_Index = NUM_Total_LCLGRANID
                        LCLGRANID(Array_Index) = New_LCLGRANID(i1)
                     Else
                        Array_Index = MAX_NUM_LCLPNTRS

                        write(msg25_2,'(i25)') NUM_Total_LCLGRANID
                        rtn_loc = String_Loc(msg25_2,fbyte2,lbyte2)

                        msgbuf = 'Total number of MODIS file inputs ' // msg25_2(fbyte2:lbyte2)
     1                  // ' exceeds '// FUNCNAME
     2                  // char(10) // 'internal storage buffer size.'
     3                  // char(10) // 'New MODIS_File ' //  New_LCLGRANID(i1)
     4                  // ' not aggregated to ECS '// 'LocalInputGranuleID array.'
     5                  // char(10) // 'Operator Action:  Notify SDST.'

                        Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                     End If

                  Else
                     msgbuf =
     1               'LocalGranuleID retrieved from MODIS input product file is blank - '
     2               // char(10) // 'not aggregated to LocalInputGranuleID metadata array.'
     3               // char(10) // 'Operator Action:  Notify SDST.'

                     Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

               End If

100         Continue

c-----------report the number of new and duplicate file references
            write(msg25_3,'(i25)') NUM_Unique
            rtn_loc = String_Loc(msg25_3,fbyte3,lbyte3)

            write(msg25_4,'(i25)') NUM_Duplicates
            rtn_loc = String_Loc(msg25_4,fbyte4,lbyte4)


c--------report the number of unique versus duplicate input MODIS file references
c--------use different messages depending on whether core metadata exists

            If (NUM_Old_LCLGRANID .LE. 0) Then  ! implies core does not exist.

               msgbuf =
     1         'No preexisting ' // MCORE_LOCALINPUTGRANULEID // ' array in MOD06 product file - '
     2         // char(10) // msg25_3(fbyte3:lbyte3) // ' unique MODIS input file names to be set '
     3         // 'by current process.'
     4         // char(10) // msg25_4(fbyte4:lbyte4) // ' repeat file names in input stream not '
     5         // 'duplicated to ' // MCORE_LOCALINPUTGRANULEID

            Else

               msgbuf =
     1         msg25_3(fbyte3:lbyte3) // ' unique MODIS input file names added to preexisting '
     2         // char(10) // MCORE_LOCALINPUTGRANULEID // ' array.  ' // msg25_4(fbyte4:lbyte4) //
     3         ' input file names nonunique and not duplicated.'

            End If

            Call MODIS_SMF_SetDynamicMsg(MODIS_M_GENERIC,msgbuf,FUNCNAME)

         Else
            error_flag = .TRUE.

            msgbuf = 'Retrieve_LclGranID_Set detected error reading local name '
     1      // char(10) // 'of input product.'
     2      // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3      // char(10) // 'messages originating from call to routine Retrieve_LclGranID_Set.'


            Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         End If   ! Retrieve LocalGranuleId from MODIS inputs

      End If   !  If (NUM_MODIS_InputFiles .GT. 0)


      If (.NOT. error_flag) rtn_status = SUCCEED

      Return

      End



C-----------------------------------------------------------------------

      Subroutine Retr_LocalGranID(NUM_Of_MODIS_InputFiles,
     1                            LUN_MODIS_InputFiles,
     2                            VRSN_MODIS_InputFiles,
     3                            NUM_LocalGranID,
     4                            LocalGranID,
     5                            rtn_flag)

      implicit none

      include 'mod06_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Retrieve ECS metadata attribute LOCALGRANULEID from 
C                a set of MODIS files identified by PCF LUN and Version
C                numbers.
C
C
C !INPUT PARAMETERS:
C
C  integer NUM_Of_MODIS_InputFiles
C                Variable containing the number of MODIS input
C                products used by process.
C
C  integer LUN_MODIS_InputFiles(*)
C                Array containing MODIS input product LUNs used by
C                process (See also VRSN_MODIS_InputFiles below).
C
C  integer VRSN_MODIS_InputFiles(*)
C                Array containing file version numbers of MODIS input
C                products (See LUN_MODIS_InputFiles above).
C
C                A one-to-one correspondence between the elements of
C                the arrays LUN_MODIS_InputFiles and
C                VRSN_MODIS_InputFiles is assumed.
C
C !OUTPUT PARAMETERS:  See INPUT/OUTPUT PARAMETERS
C
C !INPUT/OUTPUT PARAMETERS:
C
C  character*(*) LocalGranID(*)
C                Array of LOCALGRANULEIDs retrieved from array of input
C                file references (LUN_MODIS_InputFiles and VRSN_MODIS_InputFiles)
C
C  integer NUM_LocalGranID
C                Number of LOCALGRANULEID attributes successfully retrieved
C                from input file references array (LUN_MODIS_InputFiles and
C                VRSN_MODIS_InputFiles)
C
C  integer rtn_flag
C                Procedure return flag (SUCCEED=0, FAIL=-1)
C
C
C
C !REVISION HISTORY:
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by Rich Hucek, February 1998
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if an error occurs
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_GetPCAttr_s        (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG    (science code)
C        String_Loc                 (science code)
C
C    Named Constant:
C        FAIL                       (mod06_ECSMET.inc)
C        MCORE_LOCALINPUTGRANULEID  (mapi.inc)
C        MECS_CORE                  (mapi.inc)
C        MODIS_E_GENERIC            (PGS_MODIS_39500.f)
C        PGSd_PC_VALUE_LENGTH_MAX   (PGS_PC.f)
C        PGS_S_SUCCESS              (PGS_SMF.f)
C        SUCCEED                    (mod06_ECSMET.inc)
C
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*) FUNCNAME
      parameter (FUNCNAME = 'Retr_LocalGranID')

c-----Declaration of function arguments
      character*(*) LocalGranID(*)

      integer NUM_Of_MODIS_InputFiles
      integer LUN_MODIS_InputFiles(*)
      integer VRSN_MODIS_InputFiles(*)
      integer NUM_LocalGranID,rtn_flag


c-----Declaration of local variables
      character*25  msg25a, msg25b, msg25c
      character*60  AttrN
      character*512 msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) AttrV

      integer pgs_met_getpcattr_s,string_loc
      integer fbytea,fbyteb,fbytec,
     2        lbytea,lbyteb,lbytec,
     3        LUN,VRSN,icounter,rtn


C------------------------
C Initialization
C------------------------
      rtn_flag        = SUCCEED
      AttrN           = MCORE_LOCALGRANULEID
      NUM_LocalGranID = 0


C-----------------------------------------------------------------------
C Aggregate metadata attribute LOCALGRANULEID over list of input file
C references
C-----------------------------------------------------------------------

c-----loop over the number of MODIS input files
      Do 100 icounter = 1, NUM_Of_MODIS_InputFiles
         VRSN = VRSN_MODIS_InputFiles(icounter)
         LUN  = LUN_MODIS_InputFiles(icounter)

c-----set message buffers for error reporting
         write(msg25a,'(i25)') LUN
         write(msg25b,'(i25)') VRSN
         write(msg25c,'(i25)') icounter

         rtn = String_Loc(msg25a,fbytea,lbytea)
         rtn = String_Loc(msg25b,fbyteb,lbyteb)
         rtn = String_Loc(msg25c,fbytec,lbytec)

c--------check file version number
         If (VRSN .LT. 1)  Then
            rtn_flag = FAIL

            msgbuf = 'Unexpected PCF version number (' // msg25b(fbyteb:lbyteb)
     1      // ') of MODIS input file passed to ' // FUNCNAME
     2      // char(10) // MCORE_LOCALGRANULEID // ' will not be set for element '
     3      // msg25c(fbytec:lbytec) // ' of input arrays.'
     4      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------check file LUN
         Else If ( (LUN .LE. 0) .OR. (LUN.GE.10000 .AND. LUN.LE.10999) ) Then
            rtn_flag = FAIL

            msgbuf = 'MODIS input file LUN number ' // msg25a(fbytea:lbytea) //
     1      ' passed to ' // FUNCNAME // ' is out of bounds.'
     2      // char(10) // MCORE_LOCALGRANULEID // ' cannot be retrieved '
     3      // 'for element ' // msg25c(fbytec:lbytec) // ' of input arrays.'
     4      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------retrieve LOCALGRANULEID
         Else
            rtn = pgs_met_getpcattr_s(LUN,VRSN,MECS_CORE,AttrN,AttrV)

            If (rtn .ne. PGS_S_SUCCESS) Then
               rtn_flag = FAIL

               msgbuf =
     1         'pgs_met_getpcattr_s detected error retrieving ' // MCORE_LOCALGRANULEID
     2         // char(10) // 'for input array element ' // msg25c(fbytec:lbytec)
     3         // char(10) // 'File PCF LUN = ' // msg25a(fbytea:lbytea)
     4         // ';  Version number = ' // msg25b(fbyteb:lbyteb)
     5         // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     6         // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7         // char(10) // 'fault is identified, stage correct PCF/input file and '
     8         // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            Else
               NUM_LocalGranID = NUM_LocalGranID + 1
               LocalGranID(NUM_LocalGranID) = AttrV
            End If

         End If

100   Continue

      Return

      End
