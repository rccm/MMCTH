      Integer Function Set_InvMet_MOD04( ExtGeoPntr_Flag,
     1                                   NUM_InputPntr,LUN_InputPntr,VRSN_InputPntr,
     2                                   NUM_MeasParm,Name_MeasParm,PctMissing_MeasParm,
     3                                   AutoFlag_MeasParm,AutoFlagExp_MeasParm,
     4                                   NUM_PSA,Name_PSA,Value_PSA,
     5                                   MET_Handles )

      implicit none

      include 'mod04_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-------------------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Set MOD04_L2 (MODIS L2 Aersol product) ECS "Inventory" metadata objects to 
C   internal memory.  In the context used here, set means to associate a value 
C   with an ODL (Object Description Language) object defined in the MOD04_L2 
C   metadata configuration file (MCF).  Function calls to the ECS Science Data 
C   Processing Toolkit (SDPTK) are used for this purpose. No metadata fields are
C   actually written to the product file by Set_InvMet_MOD04 - only the 
C   association of values to ODL objects is made at this time.  Actual insertion
C   of Archive metadata into the HDF product file takes place later with a call 
C   to the SDPTK routine PGS_MET_Write.  This call is made from outside 
C   Set_InvMet_MOD04, after all other ECS metadata relevant to the process have 
C    been set.
C
C   The specific metadata fields set by Set_InvMet_MOD04, their
C   origin, and the value of the Metadata Configuration File (MCF)
C   "Data Location" parameter are listed below.
C      GEO = Geolocation product granule
C      MCF = Metadata Configuration File
C      PCF = Process Control File
C      RP  = PCF runtime Parameter (RP)
C      PGE = Science code
C      TK  = SDP Toolkit
C
C                                          Source of    Data Location
C     ECS Core Metadata Objects              Value      Value in MCF
C   -------------------------------        ---------    -------------
C      1    ShortName                         MCF           MCF
C      2    VersionID                         MCF           MCF
C      4    ReprocessingActual                PCF (RP)      PCF
C      5    ReprocessingPlanned               PCF (RP)      PCF
C      6    LocalGranuleID                  MCF,GEO         PGE
C      7    DayNightFlag                      GEO           PGE
C      8    ProductionDateTime                TK            TK
C      9    LocalVersionID                    PCF (RP)      PCF
C     10    PGEVersion                        PCF (RP)      PCF
C     11    InputPointer                      PCF           PGE
C     12    RangeBeginningTime                GEO           PGE
C     13    RangeEndingTime                   GEO           PGE
C     14    RangeBeginningDate                GEO           PGE
C     15    RangeEndingDate                   GEO           PGE
C     16    EastBoundingCoordinate            GEO           PGE
C     17    WestBoundingCoordinate            GEO           PGE
C     18    NorthBoundingCoordinate           GEO           PGE
C     19    SouthBoundingCoordinate           GEO           PGE
C     20    OrbitNumber                       GEO           PGE
C     21    EquatorCrossingLongitude.1        GEO           PGE
C     22    EquatorCrossingDate.1             GEO           PGE
C     23    EquatorCrossingTime.1             GEO           PGE
C     24    ParameterName.n                   PGE           PGE
C     25    AutomaticQualityFlag.n            PGE           PGE
C     26    AutomaticQualityFlagExplanation.n PGE           PGE
C     31    QAPercentMissingData.n            PGE           PGE
C     32    AdditionalAttributeName.n (PSA)   PGE           PGE
C     33    ParameterValue.n (PSA value)      PGE           PGE
C     34    AncillaryInputType.1       Set_InvMet_MOD04     PGE
C     35    AncillaryInputPointer.1           PCF           PGE
C     36    AssociatedPlatformShortName.1     PCF (RP)      PCF
C     37    AssociatedInstrumentShortName.1   PCF (RP)      PCF
C     38    AssociatedSensorShortName.1       PCF (RP)      PCF
C
C !INPUT PARAMETERS:
C  logical ExtGeoPntr_Flag  Variable (.TRUE. or .FALSE. ) indicating
C                           whether the current process produces 1-km data
C                           fields and thus requires a metadata reference
C                           to the 1-km Geolocation (MOD03) product
C
C                           MOD04CT sets this flag to .FALSE.; MOD04CD
C                           and MOD04OD set it to .TRUE.
C
C  integer  NUM_InputPntr   Variable containing the number of input
C                           files used by the current process.
C
C                           Included in the count are ancillary data files,
C                           look up tables, and MODIS product files.  System
C                           files, such as the MCF and PCF, are not counted.
C
C  integer  LUN_InputPntr(*)
C                           Array containing the PCF LUNS of the input
C                           files used by the current process
C                           (See NUM_InputPntr above).
C
C                           A one-to-one correspondence between the elements
C                           of the arrays LUN_InputPntr and VRSN_InputPntr
C                           is assumed.
C
C  integer  VRSN_InputPntr(*)
C                           Array containing PCF version numbers of the
C                           input files used by the current process.
C                           (See NUM_InputPntr and LUN_InputPntr above.)
C
C  integer NUM_MeasParm     Variable containing the number of
C                           MeasuredParameters to be set by the current
C                           process.  If none, pass NUM_MeasParm=0.
C
C  character*(*) Name_MeasParm(*)
C                           Array containing the names of the
C                           MeasuredParameters to be set by the
C                           current process.  If none, pass
C                           NUM_MeasParm=0.
C
C                           A one-to-one association between the
C                           elements of arrays PctMissing_MeasParm,
C                           Name_MeasParm, AutoFlag_MeasParm and
C                           AutoFlagExp_MeasParm is assumed.
C
C  integer PctMissing_MeasParm(*)
C                           Array containing the percentages of missing
C                           data for the MeasuredParameters to be set
C                           by the current process.
C                           (See Name_MeasParm above).
C
C                           Note that the percentage missing data are
C                           to be passed as integers values.
C
C  character*(*) AutoFlag_MeasParm(*)
C                           Array of flags referring to the quality
C                           of the MeasuredParameters to be set by the
C                           current process.
C                           (See Name_MeasParm above).
C
C                           Valid quality values are:
C                           "Passed"/"Failed"/"Suspect".
C
C  character*(*) AutoFlagExp_MeasParm(*)
C                           Array containing the text explanations of the
C                           "criteria" used to define the AutomaticQualityFlags
C                           to be set by the current process.
C                           (See Name_MeasParm and AutoFlag_MeasParm
C                           above).
C
C  integer NUM_PSA          Variable containing the number of
C                           Product-Specific Attributes (PSAs) to be
C                           set by current process.
C
C  character*(*) Name_PSA(*)
C                           Array containing the names of PSAs to be set
C                           by current process.
C
C                           A one-to-one correspondence between the
C                           elements of arrays Name_PSA and Value_PSA
C                           is assumed.
C
C  Real Value_PSA(*)        Array containing the values of the PSAs to
C                           be set by the current process (See Name_PSA
C                           above).  Values are restricted to the range
C                           +99999.99 and -9999.99 and are set to the
C                           product in F8.2 format;  A value of -99999.0
C                           is interpreted as a flag meaning the PSA is
C                           not to be set.
C
C !OUTPUT PARAMETERS:
C  character*(*) MET_Handles(*)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           MET_Handles is defined during initialization
C                           of the ECS Metadata Configuration File (MCF).
C                           In FORTRAN, element 1 of the MET_Handles
C                           array refers to the MCF file.  ECS inventory
C                           metadata are referenced as element 2 and
C                           archive metadata as element 3.
C
C
C !REVISION HISTORY:
C 2001/06/03  rhucek
C Revised messaging within function Set_InputPntr_MOD04.  Trailing blanks removed
C from file names and references to "MODIS input file" replaced by "input file".
C
C 2001/03/30  rhucek
C Renamed subprogram "Get_LUN_Of_LocalVrsnID_Atmos" to "Get_LUN_Of_LclVrsnID_Atmos"
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
C    Written by         Richard Hucek 2/1998
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        Get_LUN_Of_LclVrsnID_Atmos (science code)
C        PGS_MET_Init               (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG    (science code)
C        Set_AncInputPntr_MOD04     (science code)
C        Set_GeoData_Atmos          (science code)
C        Set_InputPntr_MOD04        (science code)
C        Set_InvPSA_MOD04           (science code)
C        Set_LclGranID              (science code)
C        Set_MeasParm_MOD04         (science code)
C        Set_RP_Data_Atmos          (science code)
C        string_loc                 (science code)
C
C
C    Named Constant:
C        FAIL                      (mod04_ECSMET.inc)
C        INVENTORYMETADATA         (mod04_ECSMET.inc)
C        LUN_MCF_MOD04             (mod04_ECSMET.inc)
C        LUN_OF_NUM_INV_RP_MOD04   (mod04_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        PGSd_MET_GROUP_NAME_L     (PGS_MET.f: included in mapi.inc)
C        SUCCEED                   (mod04_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Set_InvMet_MOD04')

C Function argument declarations
      character*(*) MET_Handles(*)
      character*(*) AutoFlag_MeasParm(*),AutoFlagExp_MeasParm(*),Name_MeasParm(*)
      character*(*) Name_PSA(*)

      integer       NUM_InputPntr,LUN_InputPntr(*),VRSN_InputPntr(*)
      integer       NUM_MeasParm,NUM_PSA
      integer       PctMissing_MeasParm(*)
      integer       LUN_Of_LclVrsnID

      logical       ExtGeoPntr_Flag

      real          Value_PSA(*)


C other variable declarations
      character*25   msg25_MET,msg25_MCF,msg25_LUN
      character*2048 msgbuf

      integer i, len_MET,
     1        rtn,rtn_flag,rtn_loc

      integer Set_InvPSA_MOD04,   Set_InputPntr_MOD04,
     1        Set_MeasParm_MOD04, Set_AncInputPntr_MOD04,
     2        Set_RP_Data_Atmos,  Set_GeoData_Atmos,
     3        Set_LclGranID,      String_Loc,
     4        pgs_met_init

      integer fbyte_MET, fbyte_MCF, fbyte_LUN,
     1        lbyte_MET, lbyte_MCF, lbyte_LUN

      logical error_flag


C------------------------
C Initialization
C------------------------

      Set_InvMet_MOD04 = FAIL
      error_flag = .false.

      len_MET = LEN( MET_Handles(1) )

      write(msg25_MCF,'(I25)') LUN_MCF_MOD04
      write(msg25_LUN,'(I25)') LUN_OF_NUM_INV_RP_MOD04  
      write(msg25_MET,'(I25)') len_MET

      rtn_loc = string_loc(msg25_MCF,fbyte_MCF,lbyte_MCF)
      rtn_loc = string_loc(msg25_LUN,fbyte_LUN,lbyte_LUN)
      rtn_loc = string_loc(msg25_MET,fbyte_MET,lbyte_MET)

      Do i = 1, PGSd_MET_NUM_OF_GROUPS
         MET_Handles(i) = ' '
      EndDo


c-------------------------------------------------------------------------------
c Check string len of MET_Handles array
c-------------------------------------------------------------------------------

      If (len_MET .lt. PGSd_MET_GROUP_NAME_L) Then

         msgbuf =
     1   'MET_Handles array element size (' // msg25_MET(fbyte_MET:lbyte_MET)
     2   // ' characters) < PGSd_MET_GROUP_NAME_L.'
     3   // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )

         return
      End If


c-------------------------------------------------------------------------------
c Initialize MOD04 MCF file which also sets inventory attribute fields whose
c values are defined in the MCF.
c-------------------------------------------------------------------------------

      rtn = pgs_met_init( LUN_MCF_MOD04, MET_Handles )

      If (rtn.ne.PGS_S_SUCCESS) Then
         msgbuf = 'pgs_met_init unable to initialize MOD04 MCF on LUN='
     1   // msg25_MCF(fbyte_MCF:lbyte_MCF) // '.'
     2   // char(10) // 'Operator Action:  Check for valid MCF file. If wrong or corrupted,'
     3   // char(10) // 'stage correct MCF and rerun PGE. Otherwise, notify SDST.'


         Call MODIS_SMF_SETDYNAMICMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )

         Return
      End If

C-----------------------------------------------------------------------
C Set MOD04 inventory attribute fields whose values are to be retrieved
C from the Geolocation product.  These fields include:
C
C "DayNightFlag"
C "EastBoundingCoordinate"
C "NorthBoundingCoordinate"
C "SouthBoundingCoordinate"
C "WestBoundingCoordinate"
C "OrbitNumber"
C "RangeBeginningTime"
C "RangeEndingTime"
C "RangeBeginningDate"
C "RangeEndingDate"
C "EquatorCrossingLongitude"
C "EquatorCrossingDate"
C "EquatorCrossingTime".
C---------------------------------------------------------------------

      rtn = Set_GeoData_Atmos( Met_Handles )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_GeoData_Atmos detected error setting MOD04 inventory '
     1   // char(10) // 'metadata retrieved from MODIS Geolocation product.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_GeoData_Atmos.'


         Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
      EndIf


C----------------------------------------------------------------------
C Set MOD04 inventory attributes whose values are retrieved as USER
C DEFINED RUNTIME PARAMETERs from PCF.  The following fields are
C included:
C
C "ReprocessingActual"
C "ReprocessingPlanned"
C "LocalVersionID"
C "PGEVersion"
C "AssociatedPlatformShortName.1"
C "AssociatedInstrumentShortName.1"
C "AssociatedSensorShortName.1"
C----------------------------------------------------------------------

      rtn = Set_RP_Data_Atmos( Met_Handles,
     1                         LUN_OF_NUM_INV_RP_MOD04,
     2                         INVENTORYMETADATA )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_RP_Data_Atmos detected error setting MOD04 inventory RP data.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_RP_Data_Atmos.'


         Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
      EndIf


C-----------------------------------------------------------------
C Set MOD04 InputPointer
C-----------------------------------------------------------------

      rtn = Set_InputPntr_MOD04( MET_Handles,
     1                           NUM_InputPntr,
     2                           LUN_InputPntr,
     3                           VRSN_InputPntr )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_InputPntr_MOD04 detected error incrementing MOD04 '
     1    // 'INPUTPOINTER.'
     2    // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     3    // char(10) // 'messages originating in function Set_InputPntr_MOD04.'


         Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
      EndIf


C-----------------------------------------------------------------
C Set MOD04_L2 Measured Parameters
C-----------------------------------------------------------------

      rtn = Set_MeasParm_MOD04( MET_Handles,
     1                          NUM_MeasParm,
     2                          Name_MeasParm,
     3                          PctMissing_MeasParm,
     4                          AutoFlag_MeasParm,
     5                          AutoFlagExp_MeasParm )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_MeasParm_MOD04 detected error setting MOD04 MeasuredParameters.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_MeasParm_MOD04.'


         Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
      EndIf


C-----------------------------------------------------------------
C Set MOD04 Product Specific Metadata
C-----------------------------------------------------------------

      rtn = Set_InvPSA_MOD04( MET_Handles,
     1                        NUM_PSA,
     2                        Name_PSA,
     3                        Value_PSA )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf =
     1    'Set_InvPSA_MOD04 detected error setting MOD04 Product Specific Attributes (PSAs).'
     2    // char(10) // 'Operator Action:  Refer to prior low level LogStatus error messages '
     3    // char(10) // 'originating in function Set_InvPSA_MOD04.'

         Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
      EndIf


C-----------------------------------------------------------------------
C Set MOD04 "LocalGranuleID"
C-----------------------------------------------------------------------

      Call Get_LUN_Of_LclVrsnID_Atmos( LUN_OF_NUM_INV_RP_MOD04,
     1                                 LUN_Of_LclVrsnID,
     2                                 rtn_flag)

      If (rtn_flag .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Get_LUN_Of_LclVrsnID_Atmos detected error identifying MOD04 '
     1   // char(10) // 'LOCALVERSIONID LUN in RP group beginning at LUN='
     2   // msg25_LUN(fbyte_LUN:lbyte_LUN) // '.'
     3   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4   // char(10) // 'messages originating in routine Get_LUN_Of_LclVrsnID_Atmos.'


         Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )

      Else
         rtn = Set_LclGranID( Met_Handles, LUN_Of_LclVrsnID )

         If (rtn .eq. FAIL) Then
            error_flag = .true.

            msgbuf = 'Set_LclGranID detected error setting MOD04 LocalGranuleID.'
     1      // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2      // char(10) // 'messages originating in function Set_LclGranID.'


            Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
         EndIf

      EndIf


C-----------------------------------------------------------------
C Set MOD04 "AncillaryInputPointer" group
C-----------------------------------------------------------------
C
C      rtn = Set_AncInputPntr_MOD04( Met_Handles,
C     1                              ExtGeoPntr_Flag )
C
C      If (rtn .eq. FAIL) Then
C         error_flag = .TRUE.
C
C         msgbuf = 'Set_AncInputPntr_MOD04 detected error setting MOD04 '
C     1    // char(10) // 'AncillaryInputGranule metadata group.'
C     2    // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
C     3    // char(10) // 'messages originating in function Set_AncInputPntr_MOD04.'
C
C
C         Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
C      EndIf
 
 
 
C--------------------------------------------------------------------------------------
C If error_flag is false, exit SUCCEED; otherwise exit -1 
C--------------------------------------------------------------------------------------

       If (.not.error_flag) Set_InvMet_MOD04 = SUCCEED
 
       Return
 
       End



C--------------------------------------------------------------------------------------
      Integer Function Copy_MeasParm_MOD04(MET_Handles)

      implicit none

      include 'mod04_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MET_13.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Copy existing ECS MeasuredParameter inventory metadata
C                groups from the MOD04_L2 product file to internal memory.
C                MOD04_L2 ECS Measured Parameter groups consist of the
C                following 4 required ECS ODL objects:
C
C      "PARAMETERNAME"
C      "QAPERCENTMISSINGDATA"
C      "AUTOMATICQUALITYFLAG"
C      "AUTOMATICQUALITYFLAGEXPLANATION"
C
C
C !INPUT PARAMETERS:
C
C  character*(*) MET_Handles(*)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           MET_Handles is defined during initialization
C                           of the ECS Metadata Configuration File (MCF).
C                           In FORTRAN, element 1 of the MET_Handles
C                           array refers to the MCF file.  ECS inventor 
C                           metadata are referenced as element 2 and
C                           archive metadata as element 3.
C
C
C !OUTPUT PARAMETERS:       None
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by         Richard Hucek
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
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
C        FAIL                    (mod04_ECSMET.inc)
C        FUNCNAME_PGS_MET_SET_i  (mod04_ECSMET.inc)
C        FUNCNAME_PGS_MET_SET_s  (mod04_ECSMET.inc)
C        INVENTORYMETADATA       (mod04_ECSMET.inc)
C        LUN_MOD04               (mod04_ECSMET.inc)
C        MCORE_*                 (mapi.inc)
C        MECS_CORE               (mapi.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        MODIS_M_GENERIC         (PGS_MODIS_39500.f)
C        NUM_MEAS_PARM_MOD04     (mod04_ECSMET.inc)
C        PGSMET_E_NULL_PARAMETER (PGS_MET_13.f)
C        PGS_S_SUCCESS           (PGS_SMF.f)
C        SUCCEED                 (mod04_ECSMET.inc)
C        VRSN_MOD04              (mod04_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'Copy_MeasParm_MOD04')


C function arguments declarations
      character*(*) MET_Handles(*)

C other variable declarations
      character*25   msg25d,msg25_LUN,msg25_VRSN
      character*60   Name
      character*128  AttrN,AttrV_s,Field_Name(4),HDF_AttrName
      character*512  msgbuf_GET_err,msgbuf_SET_err,msgbuf_GET_NULL
      character*1024 msgbuf

      integer i,i1,rtn,rtn_loc
      integer AttrV_i
      integer pgs_met_setattr_s, pgs_met_setattr_i,String_Loc,
     2        PGS_MET_GetPCAttr_i,PGS_MET_GetPCAttr_s
      integer fbyte,fbytea,fbyteb,lbyte,lbytea,lbyteb
      integer fbyte_LUN,fbyte_VRSN,lbyte_LUN,lbyte_VRSN

      logical error_flag


C------------------------
C Initialization
C------------------------
      Copy_MeasParm_MOD04 = SUCCEED
      error_flag          = .FALSE.
      HDF_AttrName        = MECS_CORE

c-----------------------------------------------------------------------
c Set up status message variables
c-----------------------------------------------------------------------

      write(msg25_LUN, '(I25)') LUN_MOD04
      write(msg25_VRSN,'(I25)') VRSN_MOD04

      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)


      Field_Name(1) = MCORE_PARAMETERNAME
      Field_Name(2) = MCORE_PERCENT_MISSING
      Field_Name(3) = MCORE_AUTO_QUALITY
      Field_Name(4) = MCORE_AQUAL_FLG


c-----loop on MeasuredParameters in MOD04 product set

      Do i1 = 1, NUM_MEAS_PARM_MOD04
         Name = Name_MeasParm_MOD04(i1)
         write(msg25d,'(i25)') i1 
         rtn_loc = String_Loc(msg25d,fbyte,lbyte)
         rtn_loc = String_Loc(Name,fbytea,lbytea)


c--------loop on 4 required ODL objects

         Do i=1,4
            rtn_loc = String_Loc(Field_Name(i),fbyteb,lbyteb)
            AttrN   = Field_Name(i)(fbyteb:lbyteb) // '.' // msg25d(fbyte:lbyte)
            rtn_loc = String_Loc(AttrN,fbyteb,lbyteb)


            msgbuf_GET_NULL =
     1      'MeasuredParameter ' // AttrN(fbyteb:lbyteb) // ' not yet written to ' 
     2      // MECS_CORE // '. ' 
     3      // char(10) // 'Operator Action:  Disregard prior SDPTK error message '
     4      // 'PGSMET_E_NULL_PARAMETER.'


            msgbuf_GET_err =
     1      ' unable to retrieve ECS attribute '// AttrN(fbyteb:lbyteb)
     2      // char(10) // 'in MeasuredParameter group ' // Name(fbytea:lbytea) // '.'
     3      // char(10) // 'MOD04 target file is on LUN = '// msg25_LUN(fbyte_LUN:lbyte_LUN)
     4      // ' and VRSN = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.'
     5      // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     6      // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     7      // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'


            msgbuf_SET_err =
     1      ' unable to set ECS attribute '// AttrN(fbyteb:lbyteb) // ' '
     2      // char(10) // 'in MeasuredParameter group ' //  Name(fbytea:lbytea) // '.'
     3      // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     4      // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     5      // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'


c-----------data type of ODL object is integer
            If (i.eq.2) Then
               rtn = PGS_MET_GetPCAttr_i(LUN_MOD04, VRSN_MOD04, HDF_AttrName, AttrN, AttrV_i)

c--------------read and set integer ODL object
               If (rtn .eq. PGS_S_SUCCESS) Then
                  rtn = PGS_MET_SetAttr_i(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_i)

                  If (rtn .ne. PGS_S_SUCCESS) Then
                     error_flag = .true.
                     msgbuf = FUNCNAME_PGS_MET_SET_i // msgbuf_SET_err
                     call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

               Else If (rtn .eq. PGSMET_E_NULL_PARAMETER) Then
                  msgbuf = msgbuf_GET_NULL
                  Call MODIS_SMF_SetDynamicMsg(MODIS_M_GENERIC,msgbuf,FUNCNAME)
               Else
                  error_flag = .true.
                  msgbuf = FUNCNAME_PGS_MET_GET_i // msgbuf_GET_err
                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If

c-----------data type of ODL object is character
            Else
               rtn = PGS_MET_GetPCAttr_s(LUN_MOD04, VRSN_MOD04, HDF_AttrName, AttrN, AttrV_s)

c--------------read and set character ODL object
               If (rtn .eq. PGS_S_SUCCESS) Then
                  rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_s)

                  If (rtn .ne. PGS_S_SUCCESS) Then
                     error_flag = .true.
                     msgbuf = FUNCNAME_PGS_MET_SET_s // msgbuf_SET_err
                     call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

               Else If (rtn .eq. PGSMET_E_NULL_PARAMETER) Then
                  msgbuf = msgbuf_GET_NULL
                  Call MODIS_SMF_SetDynamicMsg(MODIS_M_GENERIC,msgbuf,FUNCNAME)
               Else
                  error_flag = .true.
                  msgbuf = FUNCNAME_PGS_MET_GET_s // msgbuf_SET_err
                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If

            End If   ! data type of ODL object (integer or string)

         END DO   ! number of ODL objects per MearuredParameter

      END DO    ! number of MeasuredParameters

      If (error_flag) Copy_MeasParm_MOD04 = FAIL

      Return
      END



c--------------------------------------------------------------------------------------
      Integer Function Set_AncInputPntr_MOD04(MET_Handles,ExtGeoPntr_Flag)

      implicit none

      include 'mod04_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_MET_13.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Update ECS and user-defined product specific "Inventory" metadata
C   previously written to the MODIS L2 aerosol product (MOD04_L2) by one or
C   more prior processes.  In the context used here, set means to associate
C   a value with a metadata parameter name via function calls to the ECS
C   Science Data Processing Toolkit (SDPTK).
C
C
C !INPUT PARAMETERS:
C  logical ExtGeoPntr_Flag  Variable (.TRUE. or .FALSE. ) indicating
C                           whether the current process produces 1-km data
C                           fields and thus requires a metadata reference
C                           to the 1-km Geolocation (MOD03) product
C
C                           MOD04CT sets this flag to .FALSE.; MOD04CD
C                           and MOD04OD set it to .TRUE.
C
C
C
C  character*(*) MET_Handles(*)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           MET_Handles is defined during initialization
C                           of the ECS Metadata Configuration File (MCF).
C                           In FORTRAN, element 1 of the MET_Handles
C                           array refers to the MCF file.  ECS inventory
C                           metadata are referenced as element 2 and
C                           archive metadata as element 3.
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
C    Written by         Richard Hucek 2/1998
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    rhucek@ltpmail.gsfc.nasa.gov
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
C        Set_AncInputPntr          (science code)
C        string_loc                (science code)
C
C
C    Named Constant:
C        FAIL                      (mod04_ECSMET.inc)
C        LUN_GEO                   (mod04_ECSMET.inc)
C        LUN_MCF_MOD04             (mod04_ECSMET.inc)
C        MCORE_ANCIL_POINTER       (mapi.inc)
C        MECS_CORE                 (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSd_PC_UREF_LENGTH_MAX   (PGS_PC.f)
C        PGSMET_E_NULL_PARAMETER   (PGS_MET_13.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (mod04_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Set_InvMet_MOD04')

C Function argument declarations
      character*(*) MET_Handles(*)

      logical       ExtGeoPntr_Flag


C other variable declarations
      character*25  msg25_LUN,msg25_VRSN
      character*128 AttrN
      character*512 msgbuf
      character*(PGSd_PC_UREF_LENGTH_MAX) UR_AncPntr

      integer    rtn,rtn_status,rtn_loc

      integer    PGS_MET_GetPCAttr_s,
     1           String_Loc,
     2           Set_AncInputPntr

      integer    fbyte,fbyte_LUN,fbyte_VRSN,
     1           lbyte,lbyte_LUN,lbyte_VRSN

      logical    core_exists,error_flag


C------------------------
C Initialization
C------------------------

      Set_AncInputPntr_MOD04 = FAIL
      error_flag = .false.


C-----------------------------------------------------------------
C Set ECS "AncillaryInputPointer" group
C-----------------------------------------------------------------

c-----set reference to external geolocation file
      If (ExtGeoPntr_Flag) Then
         rtn = Set_AncInputPntr(Met_Handles)

         If (rtn .ne. SUCCEED) Then
            error_flag = .TRUE.

            msgbuf =
     1      'Set_AncInputPntr detected error setting AncillaryInputGranule group. '
     2      // char(10) // 'Operator Action:  Refer to prior low level error LogStatus '
     3      // char(10) // 'messages generated by routine Set_AncInputPntr.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If

c-----do not set reference to external geolocation file unless it already exists 
      Else

         call Query_HDF_Attr( LUN_MOD04,
     1                        VRSN_MOD04,
     2                        MECS_CORE,
     3                        rtn_status,
     4                        core_exists )


         If (rtn_status .eq. FAIL) then
            error_flag = .TRUE.

            msgbuf =
     1      'Query_HDF_Attr detected error searching MOD04 aerosol product '
     2      // char(10) // 'for HDF file attribute ' // MECS_CORE // '.'
     3      // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4      // char(10) // 'messages generated by routine Query_HDF_Attr.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------core_exists - copy existing GEO reference if any
         Else if (core_exists) Then
            AttrN = MCORE_ANCIL_POINTER // '.1'
            rtn   = String_Loc(AttrN,fbyte,lbyte)

            rtn = PGS_MET_GetPCAttr_s( LUN_MOD04,
     1                                 VRSN_MOD04,
     2                                 MECS_CORE,
     3                                 AttrN,
     4                                 UR_AncPntr )

c-----------AncInputPntr already in file - copy to memory.
            If (rtn .eq. PGS_S_SUCCESS) Then
               rtn = Set_AncInputPntr(Met_Handles)

               If (rtn .ne. SUCCEED) Then
                  error_flag = .TRUE.

                  msgbuf =
     1            'Set_AncInputPntr detected error setting AncillaryInputGranule group.'
     2            // char(10) // 'Operator Action:  Refer to prior low level error LogStatus '
     3            // char(10) // 'messages generated by routine Set_AncInputPntr.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If

c-----------AncInputPntr not in core metadata and not to be set current process. 
            Else If ( rtn .eq. PGSMET_E_NULL_PARAMETER ) Then

               msgbuf =
     1         AttrN(fbyte:lbyte) // ' not yet written to ' // MECS_CORE // '. '
     2         // char(10) // 'Operator Action:  Disregard prior SDPTK error message '
     3         // 'PGSMET_E_NULL_PARAMETER.'

               Call MODIS_SMF_SetDynamicMsg(MODIS_M_GENERIC,msgbuf,FUNCNAME)

c-----------unknow error probing core metadata
            Else
               error_flag = .TRUE.

               write(msg25_LUN,'(I25)') LUN_GEO
               rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

               write(msg25_VRSN,'(I25)') VRSN_GEO
               rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)

               msgbuf = 'PGS_MET_GetPCAttr_s unable to retrieve '// AttrN(fbyte:lbyte)
     1         // char(10) // 'from MOD04 product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2         // ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     4         // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     5         // char(10) // 'fault is identified, stage correct PCF/input file and '
     6         // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            End If   ! check return from PGS_MET_GetPCAttr_s

         End If   !  check return from Query_HDF_Attr

      End If   ! check on ExtGeoPntr_Flag


      If (.not.error_flag) Set_AncInputPntr_MOD04 = SUCCEED

      Return

      End



c-------------------------------------------------------------------------------
      Integer Function Set_InputPntr_MOD04( MET_Handles,
     1                                      NUM_InputPntr,
     2                                      LUN_InputPntr, 
     3                                      VRSN_InputPntr )

      implicit none

      include 'mod04_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_MET_13.f'
      include 'PGS_MODIS_39500.f'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'

C-------------------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Increment the Universal Reference (UR) list of input
C                files used to generate the MODIS Aerosol Product (MOD04).
C                The specific inputs of the current process are added to
C                the existing UR list, if any, and the updated list is
C                set to memory.  In a subsequent processing step, these
C                data are written to the product as ECS Inventory
C                metadata (in attribute INPUTPOINTER)
C
C
C !INPUT PARAMETERS:
C  character*(*) MET_Handles(*)
C                           Array containing the names of "MASTER"
C                           groups defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           MET_Handles is defined during initialization
C                           of the ECS Metadata Configuration File (MCF).
C                           In FORTRAN, element 1 of the MET_Handles
C                           array refers to the MCF file.  ECS inventory
C                           metadata are referenced as element 2 and
C                           archive metadata as element 3.
C
C  integer  NUM_InputPntr   Variable containing the number of input
C                           files used by the current process.  Included
C                           in the count are ancillary data files, look up
C                           tables, and MODIS product files.  System files,
C                           such as the MCF and PCF, are not counted.
C
C  integer  LUN_InputPntr(*)
C                           Array containing the PCF LUNS of the input
C                           files used by the current process
C                           (See NUM_InputPntr above).
C
C                           A one-to-one association between the elements
C                           of the arrays LUN_InputPntr and VRSN_InputPntr
C                           is assumed.
C
C  integer  VRSN_InputPntr(*)
C                           Array containing PCF version numbers of the
C                           input files used by the current process.
C                           (See NUM_InputPntr and LUN_InputPntr above.)
C
C                           A one-to-one association between the elements
C                           of the arrays LUN_InputPntr and VRSN_InputPntr
C                           is assumed.
C
C
C !OUTPUT PARAMETERS:       None
C
C !REVISION HISTORY:
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
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        Query_HDF_Attr            (science code)
C        string_loc                (science code)
C
C
C    Named Constant:
C        FAIL                       (mod04_ECSMET.inc)
C        INVENTORYMETADATA          (mod04_ECSMET.inc)
C        LUN_MOD04                  (mod04_ECSMET.inc)
C        MAX_NUM_PNTRS              (mod04_ECSMET.inc)
C        MODIS_E_GENERIC            (PGS_MODIS_39500.f)
C        MODIS_M_GENERIC            (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS              (PGS_SMF.f)
C        PGSMET_E_PARAMETER_NOT_SET (PGS_MET_13.f)
C        PGSd_MET_STR_END           (PGS_MET.f)
C        PGSd_PC_UREF_LENGTH_MAX    (PGS_PC.f)
C        SUCCEED                    (mod04_ECSMET.inc)
C        VRSN_MOD04                 (mod04_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Set_InputPntr_MOD04')


C input argument declarations
      character*(*)    MET_Handles(*)
      integer          NUM_InputPntr,LUN_InputPntr(*),VRSN_InputPntr(*)


C local variable declarations
      character*25     msg25_2,msg25_3,msg25_4,msg25_LUN,msg25_VRSN
      character*256    HDF_FileAttrName, ECS_AttrName
      character*(1024) msgbuf
      character*(PGSd_PC_UREF_LENGTH_MAX) New_UR(MAX_NUM_PNTRS),
     1                                    UR(MAX_NUM_PNTRS) /MAX_NUM_PNTRS*PGSd_MET_STR_END/

      integer  String_Loc, PGS_MET_GetPCAttr_s, PGS_MET_SetAttr_s, PGS_MET_SetMultiAttr_s

      integer  i1, i2, rtn, rtn_loc, rtn_flag,
     1         fbyte,fbyte2,fbyte3,fbyte4,fbyte_LUN,fbyte_VRSN,fbyte_Name,
     2         lbyte,lbyte2,lbyte3,lbyte4,lbyte_LUN,lbyte_VRSN,lbyte_Name

      integer  NUM_Duplicates,NUM_UR,NUM_New_UR,NUM_Pntr,NUM_Total_UR,NUM_Unique,
     1         Array_Index, MCF_MAX_Pntr

      logical  core_exists, error_flag, In_UR_Array, RTRV_new_UR_flag

      integer  rtn_status


c-----Initialization
      Set_InputPntr_MOD04 = FAIL
      core_exists         = .FALSE.
      error_flag          = .FALSE.
      NUM_Unique          = 0
      NUM_Duplicates      = 0
      NUM_UR              = 0
      NUM_Total_UR        = 0

      rtn_status          = 0
      RTRV_new_UR_flag    = .FALSE.


c-----------------------------------------------------------------------
c Check whether ECS Core (inventory) metadata exist in product output 
c file.  If unable to make query (e.g. no file staged to LUN, file 
c corrupted, etc), exit routine with status 'FAIL'. 
c-----------------------------------------------------------------------

      HDF_FileAttrName = MECS_CORE

      call Query_HDF_Attr( LUN_MOD04,
     1                     VRSN_MOD04,
     2                     HDF_FileAttrName,
     3                     rtn_status,
     4                     core_exists)

      If (rtn_status .eq. FAIL) then

         msgbuf =
     1   'Query_HDF_Attr detected error in attempt to query ' // MECS_CORE 
     2   // ' in MODIS aerosol product.'
     3   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4   // char(10) // 'messages generated by routine Query_HDF_Attr.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Return
      End If


c-------------------------------------------------------------------------------
c If Core metadata are in file, copy existing InputPointer to internal memory
c (assumes ECS attribute "InputPointer" exists if Core metadata exist). Track
c the number of URs already in file, and total number of URs (old + new).
c-------------------------------------------------------------------------------

      HDF_FileAttrName = MECS_CORE
      ECS_AttrName     = MCORE_INPUT_POINTER
      rtn_loc          = String_Loc(ECS_AttrName,fbyte,lbyte)

      If (core_exists) Then

c--------read InputPointer
         rtn = PGS_MET_GetPCAttr_s(LUN_MOD04,VRSN_MOD04,HDF_FileAttrName,ECS_AttrName,UR)

c--------Successfully read inputpointer
         If (rtn .EQ. PGS_S_SUCCESS) Then
            NUM_UR = 1

c-----------search existing UR array for PGSd_MET_STR_END.  Note MAX_NUM_PNTRS set
c-----------equal to ECS NUM_VAL field of MCF INPUTPOINTER object.  NUM_VAL is the
c-----------maximum number of URs SDPTK will write to ECS product metadata.

            Do WHILE ( (UR(NUM_UR) .NE. PGSd_MET_STR_END) .AND.
     1                 (NUM_UR .LE. MAX_NUM_PNTRS) )
               NUM_UR = NUM_UR + 1
            End Do

            NUM_UR       = NUM_UR - 1
            NUM_Total_UR = NUM_UR


c--------inputpointer not yet set
         Else If (rtn .eq. PGSMET_E_PARAMETER_NOT_SET) Then
            error_flag = .TRUE.

            msgbuf =
     1      MCORE_INPUT_POINTER // ' not in existing ' // MECS_CORE // '. '
     2      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c--------Failed to read inputpointer
         Else
            error_flag = .TRUE.

c-----------Set up status message variables
            write(msg25_LUN,'(I25)') LUN_MOD04
            rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

            write(msg25_VRSN,'(I25)') VRSN_MOD04
            rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)

            msgbuf = 'PGS_MET_GetPCAttr_s detected error while reading ECS attribute '
     1      // ECS_Attrname(fbyte:lbyte)
     2      // char(10) // 'from MODIS aerosol product on LUN = ' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     3      // ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     4      // char(10) // ECS_Attrname(fbyte:lbyte)//' not set to HDF product file.'
     5      // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     6      // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7      // char(10) // 'fault is identified, stage correct PCF/input file and '
     8      // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If   !  INPUTPOINTER read check

      End If   ! if Core is present


c-------------------------------------------------------------------------------
c Return if error detected.
c-------------------------------------------------------------------------------

      If (error_flag) Return


c-------------------------------------------------------------------------------
c Retrieve new URs (i.e. UR set for current process), copy to existing UR
c array and eliminate duplicates.
c-------------------------------------------------------------------------------
      Array_Index = NUM_UR

      If (NUM_InputPntr .gt. 0) Then
         NUM_Pntr = NUM_InputPntr

         If (NUM_InputPntr .GT. MAX_NUM_PNTRS) NUM_Pntr = MAX_NUM_PNTRS

         Call Retrieve_UR_Set(NUM_Pntr, LUN_InputPntr, VRSN_InputPntr,
     1                        NUM_New_UR, New_UR, rtn_flag)


c--------copy new UR set to existing URs.  Eliminate duplicates
         If (rtn_flag .EQ. SUCCEED) Then
            RTRV_new_UR_flag = .TRUE.

            Do 100 i1 = 1, NUM_New_UR
               In_UR_Array = .FALSE.
               rtn_loc = String_Loc( New_UR(i1), fbyte_Name, lbyte_Name )

c--------------search against existing UR set for duplicate 
               Do 200 i2 = 1, Array_Index

                  If (UR(i2) .EQ. New_UR(i1)) Then
                     NUM_Duplicates      =  NUM_Duplicates + 1
                     In_UR_Array = .TRUE.
                  End If

200            Continue

c--------------new name is unique
               If ( .NOT. In_UR_Array ) Then

c-----------------new name is non-blank - add to InputPointer array if
c-----------------sufficient storage

                  If (New_UR(i1) .ne. BLANK) Then
                     NUM_Total_UR = NUM_Total_UR + 1
                     NUM_Unique   = NUM_Unique + 1

                     If (NUM_Total_UR .LE. MAX_NUM_PNTRS) Then
                        Array_Index     = NUM_Total_UR
                        UR(Array_Index) = New_UR(i1)

                     Else
                        error_flag  = .TRUE.
                        Array_Index = MAX_NUM_PNTRS

                        write(msg25_2,'(i25)') NUM_Total_UR
                        rtn_loc = String_Loc(msg25_2,fbyte2,lbyte2)

                        msgbuf = 'Total number of input files ' // msg25_2(fbyte2:lbyte2)
     1                  // ' exceeds '// FUNCNAME
     2                  // char(10) // 'internal buffer size.'
     3                  // char(10) // 'Unique input file ' // New_UR(i1)(fbyte_Name:lbyte_Name)
     4                  // char(10) // ' cannot be aggregated to '// 'InputPointer array.'
     5                  // char(10) // 'Operator Action:  Notify SDST.'

                        Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                     End If

                  Else
                     error_flag = .TRUE.

                     msgbuf =
     1               'UR for input product file is blank - not aggregated to InputPointer array.'
     3               // char(10) // 'Operator Action:  Notify SDST.'

                     Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

               End If   ! UR check for uniqueness 

100         Continue 

         Else   ! failed to retrieve UR set for current process
            error_flag = .TRUE.

            msgbuf = 
     1      'Retrieve_UR_Set detected error acquiring input file URs. '
     2      // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     3      // char(10) // 'messages originating within function Retrieve_UR_Set.'

            Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         End If   ! check return from Retrieve_UR_Set

      End If   ! check NUM_InputPntr


c-----If sufficient storage in UR buffer, add end-of-string marker to UR array
      If ( Array_Index .lt. MAX_NUM_PNTRS) Then
         Array_Index     = Array_Index + 1
         UR(Array_Index) = PGSd_MET_STR_END
      End If


c-----------------------------------------------------------------------
c Set revised InputPointer to SDPTK internal memory
c-----------------------------------------------------------------------

c     rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),ECS_AttrName,UR)
      rtn = PGS_MET_SetMultiAttr_s(MET_Handles(INVENTORYMETADATA),
     &                             ECS_AttrName,MAX_NUM_PNTRS,UR)


      If (rtn .EQ. PGS_S_SUCCESS) Then

c--------report the number of unique versus duplicate file references
c--------use different messages depending on whether core metadata exists

         write(msg25_3,'(i25)') NUM_Unique
         write(msg25_4,'(i25)') NUM_Duplicates

         rtn_loc = String_Loc(msg25_3,fbyte3,lbyte3)
         rtn_loc = String_Loc(msg25_4,fbyte4,lbyte4)

         If (core_exists .AND. RTRV_new_UR_flag) Then
            msgbuf =
     1      msg25_3(fbyte3:lbyte3) // ' unique input file URs added to existing INPUTPOINTER array.'
     2      // char(10) // msg25_4(fbyte4:lbyte4) // ' input URs in existing '
     3      // 'INPUTPOINTER array - not duplicated here.'

         Else If (.NOT.core_exists .AND. RTRV_new_UR_flag) Then
            msgbuf =
     1      'No preexisting ' // MCORE_INPUT_POINTER // ' array in MOD04 product file - '
     2      // char(10) // msg25_3(fbyte3:lbyte3) // ' unique input file pointers (URs) set by '
     3      // 'current process. '
     4      // char(10) // msg25_4(fbyte4:lbyte4) // ' duplicate URs in input stream not '
     5      // 'repeated in ' // MCORE_INPUT_POINTER // '.'

         Else
            msgbuf = 
     1      'No current process URs aggregated to ' // MCORE_INPUT_POINTER // ' array. '
         End If

         Call MODIS_SMF_SetDynamicMsg(MODIS_M_GENERIC,msgbuf,FUNCNAME)

      Else
         error_flag = .TRUE.

         msgbuf =
     1   'PGS_MET_SetMultiAttr_s detected error setting attribute ' 
     2   // ECS_AttrName(fbyte:lbyte) // '. '
     3   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     4   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     5   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If

      If (.NOT.error_flag) Set_InputPntr_MOD04 = SUCCEED

      Return

      End



c--------------------------------------------------------------------------------------
      Integer Function Set_InvPSA_MOD04(MET_Handles,
     1                                  NUM_NewPSA,
     2                                  Name_NewPSA,
     3                                  Value_NewPSA)


      Implicit None

      Include 'mod04_ECSMET.inc'
      Include 'mapi.inc'
      Include 'PGS_MET_13.f'
      Include 'PGS_MODIS_39500.f'
      Include 'PGS_SMF.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Increment and set the number ECS PSA attributes to be
C                stored as ECS metadata in the MOD04 product file
C
C !INPUT PARAMETERS:
C
C  integer NUM_NewPSA       Variable containing the number of
C                           additional Product-Specific Attributes (PSAs)
C                           to be set by process.
C
C  character*(*) Name_NewPSA(*)
C                           Array containing the additional PSA names
C                           to be set by process
C
C                           A one-to-one correspondence between the
C                           elements of arrays Name_NewPSA and Value_NewPSA
C                           is assumed.
C
C  Real Value_NewPSA(*)     Array containing the additional PSA values
C                           to be set by process (See Name_NewPSA
C                           above).  Values are restricted to range
C                           plus (+) and minus (-) 99999.99 and are
C                           set to the product in F8.2 format;  A
C                           value of -99999.0 is interpreted as a
C                           flag meaning the PSA is not to be set.
C
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
C !OUTPUT PARAMETERS:       None
C
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by         Richard Hucek
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        Query_HDF_Attr            (science code)
C        Set_PSA                   (science code)
C        string_loc                (science code)
C
C
C    Named Constant:
C        FAIL                    (mod04_ECSMET.inc)
C        LUN_MOD04               (mod04_ECSMET.inc)
C        MCORE_*                 (mapi.inc)
C        MECS_CORE               (mapi.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        MODIS_M_GENERIC         (PGS_MODIS_39500.f)
C        MODIS_W_GENERIC         (PGS_MODIS_39500.f)
C        PGSMET_E_NULL_PARAMETER (PGS_MET_13.f)
C        PGS_S_SUCCESS           (PGS_SMF.f)
C        SUCCEED                 (mod04_ECSMET.inc)
C        VRSN_MOD04              (mod04_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------


C-----Parameter declarations
      Character*(*)  FUNCNAME
      Parameter     (FUNCNAME = 'Set_InvPSA_MOD04')

      real           REL_EQUALITY_EPS
      parameter     (REL_EQUALITY_EPS=0.000001)


C-----function argument declarations
      Character*(*)  Name_NewPSA(*), MET_Handles(*)
      Integer        NUM_NewPSA
      Real           Value_NewPSA(*)


C-----Local variable declarations
      Character*8    msg8
      Character*25   msg25_1,msg25_2,msg25_4,msg25_LUN,msg25_VRSN
      Character*30   Name,PSA_Name
      Character*60   PSA_Value
      Character*1024 msgbuf

      Integer fbyte1,fbyte2,fbyte4,fbyte_LUN,fbyte_VRSN,
     1        lbyte1,lbyte2,lbyte4,lbyte_LUN,lbyte_VRSN,
     2        rtn1,rtn2,rtn_loc,rtn_status,
     3        i0,i1,i2,i3,NUM_Blank_PSA,PSA_Class

      Logical core_exists, error_flag, Name_Is_Valid

      Real    Value

      integer Lcl_NUM_NewPSA

c-----Functions declarations
      Integer pgs_met_getpcattr_s     ! SDPTK
      Integer Set_PSA,string_loc      ! Other


c-----Initialization
      error_flag       = .FALSE.
      Set_InvPSA_MOD04 = SUCCEED
      core_exists      = .FALSE.
      rtn_status       = 0


c-----------------------------------------------------------------------
c Perform input argument checks
c-----------------------------------------------------------------------

c-----check for positive number of PSAs
      If (NUM_NewPSA .LE. 0) Then
         write(msg25_1,'(i25)') NUM_NewPSA
         rtn_loc = String_Loc(msg25_1,fbyte1,lbyte1)

         msgbuf =
     1   'The number of new PSAs (' // msg25_1(fbyte1:lbyte1) // ') on input '
     2   // 'argument list is < 1.'
     3   // char(10) // 'No new PSAs to be set by current process, but '
     4   // 'existing ones will be copied. '

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)
      End If

c-----check for blank PSAs
      NUM_Blank_PSA = 0

      Do 50 i0 = 1,NUM_NewPSA
         If (Name_NewPSA(i0) .eq. BLANK) NUM_Blank_PSA = NUM_Blank_PSA + 1
50    Continue

      If (NUM_Blank_PSA .GT. 0) Then
         write(msg25_1,'(i25)') NUM_Blank_PSA
         rtn_loc = String_Loc(msg25_1,fbyte1,lbyte1)

         msgbuf =
     1   msg25_1(fbyte1:lbyte1) // ' PSAs on input argument list have blank names.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)
      End If


c-----------------------------------------------------------------------
c Check if ECS Core (inventory) metadata is in product output file
c-----------------------------------------------------------------------

      call Query_HDF_Attr(LUN_MOD04,
     1                    VRSN_MOD04,
     2                    MECS_CORE, 
     3                    rtn_status,
     4                    core_exists)


c--------------------------------------------------------------------------------
c Return if Error, or if there are no new PSAs to set and no existing PSA to copy
c--------------------------------------------------------------------------------

      If (rtn_status .eq. FAIL) then
         Set_InvPSA_MOD04 = FAIL

         msgbuf =
     1   'Query_HDF_Attr detected error while searching MOD04 product for '
     2   // 'presence of ' // MECS_CORE // '.'
     3   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4   // char(10) // 'messages generated by routine Query_HDF_Attr.'

         Return
      Else If ( (NUM_NewPSA .le. 0) .and. (.NOT. core_exists) ) Then
         Return
      End If


c-------------------------------------------------------------------------------
c If core metadata are in file, loop over PSAs in MODIS AEROSOL Product and copy
c existing PSAs to internal memory.
c-------------------------------------------------------------------------------

      If (core_exists) Then

         Do 100 i1 = 1, NUM_PSA_MOD04
            PSA_Name     = Name_PSA_MOD04(i1)
            PSA_Class    = i1

            rtn1 = pgs_met_getpcattr_s(LUN_MOD04,VRSN_MOD04,MECS_CORE,PSA_Name,PSA_Value)
            rtn_loc = string_loc(PSA_Name,fbyte1,lbyte1)


c-----------PSA read was successful. Copy it to memory
            If (rtn1 .EQ. PGS_S_SUCCESS) Then

c--------------retrieved string has leading blanks removed; reformat to original "f8.2" format. 
               read(PSA_Value,'(f8.2)') Value
               write(msg8,'(f8.2)') Value

               rtn2 = Set_PSA(MET_Handles,PSA_NAME,PSA_Class,msg8)

               If (rtn2 .EQ. FAIL) Then
                  error_flag = .TRUE.

                  msgbuf =
     1            'Set_PSA detected error setting MOD04 PSA to ECS internal memory.'
     2            // char(10) // 'PSA ' // PSA_Name(fbyte1:lbyte1) // ' not set.'
     3            // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4            // char(10) // 'messages from routine Set_PSA.'

                  Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If


c-----------PSA is not in core metadata.
            Else If (rtn1 .EQ. PGSMET_E_NULL_PARAMETER) Then

               If (rtn1 .EQ. PGSMET_E_NULL_PARAMETER) Then

                  msgbuf =
     1            'PSA ' // PSA_Name(fbyte1:lbyte1) // ' not yet written to ' // MECS_CORE // '.'
     2            // char(10) // 'Operator Action:  Disregard prior SDPTK error message '
     3            // 'PGSMET_E_NULL_PARAMETER.'

                  Call MODIS_SMF_SetDynamicMsg(MODIS_M_GENERIC,msgbuf,FUNCNAME)
               End If

            Else
               error_flag = .TRUE.

               write(msg25_LUN, '(I25)') LUN_MOD04
               write(msg25_VRSN,'(I25)') VRSN_MOD04

               rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)
               rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)

               msgbuf = 'PGS_MET_GetPCAttr_d unable to retrieve '// PSA_Name(fbyte1:lbyte1)
     1         // char(10) // 'from MOD04 product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2         // ' and VRSN = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3         // char(10) // PSA_Name(fbyte1:lbyte1) // ' not set.'
     4         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5         // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6         // char(10) // 'fault is identified, stage correct PCF/input file and '
     7         // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            End If   ! End PSA existence check

100      Continue   ! End loop on MOD04 PSAs

      End If


c-------------------------------------------------------------------------------
c Return if Error during PSA copy, return
c-------------------------------------------------------------------------------

      If (error_flag) then
         Set_InvPSA_MOD04 = FAIL
         Return
      End If


c-------------------------------------------------------------------------------
c Loop on list of new PSAs.  Set PSAs with valid names and values
c-------------------------------------------------------------------------------
      Lcl_NUM_NewPSA = NUM_NewPSA

c-----check for too many new PSAs
      If (NUM_NewPSA .gt. NUM_PSA_MOD04) Then
         error_flag = .true.
         Lcl_NUM_NewPSA = NUM_PSA_MOD04

         write(msg25_1,'(I25)') NUM_NewPSA
         rtn_loc = string_loc(msg25_1,fbyte1,lbyte1)
         write(msg25_2,'(I25)') NUM_PSA_MOD04
         rtn_loc = string_loc(msg25_2,fbyte2,lbyte2)

         msgbuf = 'Number of PSAs (=' // msg25_1(fbyte1:lbyte1) // ') exceeds '
     1   // char(10) // 'number on NUM_PSA_MOD04 (= ' // msg25_2(fbyte2:lbyte2) // ').'
     2   // char(10) // 'Operator Action:  Notify SDST.'

         Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      Else
         Do 200 i2 = 1, Lcl_NUM_NewPSA
            Name  = Name_NewPSA(i2)
            Value = Value_NewPSA(i2)
            rtn_loc = string_loc(Name,fbyte1,lbyte1)

c-----------check validity of new PSA name and get ECS "Class" attribute
            Name_Is_Valid = .FALSE.

            Do 300 i3 = 1, NUM_PSA_MOD04
               If (Name_PSA_MOD04(i3) .eq. Name) Then
                  Name_Is_Valid = .TRUE.
                  PSA_Class = i3
               End If
 300        Continue

c-----------invalid PSA name
            If (.NOT. Name_Is_Valid) Then
               error_flag = .TRUE.

               msgbuf = 'PSA name ' // Name(fbyte1:lbyte1) // ' not on MOD04 metadata list.'
     1         // char(10) // Name(fbyte1:lbyte1) // ' Name/Value pair not set.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----------valid PSA name
            Else

c--------------check for "No Set" option
               If ( (abs((Value-NO_SET_FLAG)/NO_SET_FLAG) .gt. REL_EQUALITY_EPS) ) Then

c-----------------PSA value out of range check; cannot set Name/Value pair
                  If (Abs(Value) .gt. MAX_ABS_VALUE_PSA) Then
                     error_flag = .TRUE.

                     write(msg25_4,'(E25.6)') Value
                     rtn_loc = String_Loc(msg25_4,fbyte4,lbyte4)

                     msgbuf =
     1               'PSA value ' // msg25_4(fbyte4:lbyte4) // ' is out of bounds! '
     2               // char(10) // 'PSA ' // Name(fbyte1:lbyte1)
     3               // ' Name/Value pair not updated.'
     4               // char(10) // 'Operator Action:  Notify SDST. '


                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----------------Set PSA Name/Value pair
                  Else
                     write(msg8,'(f8.2)') Value

                     rtn2 = Set_PSA(MET_Handles,Name,PSA_Class,msg8)

                     If (rtn2 .EQ. FAIL) Then
                        error_flag = .TRUE.

                        msgbuf =
     1                  'Set_PSA detected error setting PSA '
     2                  // Name(fbyte1:lbyte1) // ' Name/Value pair. '
     3                  // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4                  // char(10) // 'messages originating within routine Set_PSA.'

                        Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                     Else
                        msgbuf = 'PSA ' // Name(fbyte1:lbyte1) // ' Name/Value pair successfully set.'

                        Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
                     End If

                  End If   ! PSA valid range check

               End If   ! PSA "ok to set" test 

            End If   ! PSA name check

200      Continue   ! Loop on Input PSAs

      End If   ! check on NUM_NewPSA


      If (error_flag) Set_InvPSA_MOD04 = FAIL


      Return

      End



c--------------------------------------------------------------------------------------
      Integer Function Set_MeasParm_MOD04( MET_Handles,
     1                                     NUM_NewMeasParm,
     2                                     Name_NewMeasParm,
     3                                     PctMissing_NewMeasParm,
     4                                     AutoFlag_NewMeasParm,
     5                                     AutoFlagExp_NewMeasParm )

      implicit none

      include 'mod04_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Set new ECS MeasuredParameter inventory metadata groups
C                and copy existing ones from the MOD04_L2 product file to 
C                internal memory.  MOD04_L2 ECS Measured Parameter groups
C                consists of the following 4 required ECS ODL objects:
C
C      "PARAMETERNAME"                    ==> Name_NewMeasParm
C      "QAPERCENTMISSINGDATA"             ==> PctMissing_NewMeasParm
C      "AUTOMATICQUALITYFLAG"             ==> AutoFlag_NewMeasParm
C      "AUTOMATICQUALITYFLAGEXPLANATION"  ==> AutoFlagExp_NewMeasParm
C
C
C !INPUT PARAMETERS:
C
C  character*(*) MET_Handles(*)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           MET_Handles is defined during initialization 
C                           of the ECS Metadata Configuration File (MCF).
C                           In FORTRAN, element 1 of the MET_Handles
C                           array refers to the MCF file.  ECS inventory 
C                           metadata are referenced as element 2 and 
C                           archive metadata as element 3.
C
C  integer NUM_NewMeasParm  Variable containing the number of 
C                           MeasuredParameters to be set by the current 
C                           process.  If none, NUM_NewMeasParm=0.
C  
C  character*(*) Name_NewMeasParm(*)  
C                           Array containing the names of the 
C                           MeasuredParameters to be set by the 
C                           current process.  If none, be sure to pass
C                           NUM_NewMeasParm=0. 
C
C                           A one-to-one association between the 
C                           elements of arrays PctMissing_NewMeasParm,
C                           Name_NewMeasParm, AutoFlag_NewMeasParm and
C                           AutoFlagExp_NewMeasParm is assumed.
C
C  integer PctMissing_NewMeasParm(*)
C                           Array containing the percentages of missing 
C                           data for the MeasuredParameters to be set 
C                           by the current process.  
C                           (See Name_NewMeasParm above). 
C        
C                           Note that the percentage missing data are 
C                           to be passed as integers values. 
C        
C  character*(*) AutoFlag_NewMeasParm(*)  
C                           Array of flags referring to the quality 
C                           of the MeasuredParameters to be set by the 
C                           current process.  
C                           (See Name_NewMeasParm above). 
C
C                           Valid quality values are: 
C                           "Passed"/"Failed"/"Suspect".  
C
C  character*(*) AutoFlagExp_NewMeasParm(*)
C                           Array containing the text explanations of the
C                           "criteria" used to define the AutomaticQualityFlags
C                           to be set by the current process.  
C                           (See Name_NewMeasParm and AutoFlag_NewMeasParm 
C                           above). 
C
C
C !OUTPUT PARAMETERS:       None
C
C !REVISION HISTORY:
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
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        Copy_MeasParm_MOD04       (science code)
C        Query_HDF_Attr            (science code)
C        Set_NewMeasParm_MOD04     (science code)
C
C
C    Named Constant:
C        FAIL                    (mod04_ECSMET.inc)
C        LUN_MOD04               (mod04_ECSMET.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        MECS_CORE               (mapi.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        SUCCEED                 (mod04_ECSMET.inc)
C        VRSN_MOD04              (mod04_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'Set_MeasParm_MOD04')

C function arguments declarations
      character*(*) MET_Handles(*),
     1              Name_NewMeasParm(*),
     2              AutoFlag_NewMeasParm(*),
     3              AutoFlagExp_NewMeasParm(*)

      integer NUM_NewMeasParm,
     1        PctMissing_NewMeasParm(*)


C other variable declarations

      character*1024 msgbuf

      integer rtn
      integer Copy_MeasParm_MOD04, Set_NewMeasParm_MOD04 

      logical core_exists, error_flag
      integer rtn_status_query


C------------------------
C Initialization
C------------------------
      Set_MeasParm_MOD04 =  FAIL 
      error_flag         = .FALSE.
      core_exists        = .FALSE.
      rtn_status_query         =  0


c----------------------------------------------------------------------
c Check if ECS Core (inventory) metadata is in product output file 
C If Query_HDF_Attr returns FAIL, exit routine.
c-----------------------------------------------------------------------

      call Query_HDF_Attr( LUN_MOD04,
     1                     VRSN_MOD04,
     2                     MECS_CORE, 
     3                     rtn_status_query,
     4                     core_exists )

      
      If (rtn_status_query .eq. FAIL) then

         msgbuf = 
     1   'Query_HDF_Attr detected error while searching MOD04 aerosol product ' 
     2   // 'for presence of ' // MECS_CORE // '.'
     3   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error ' 
     4   // char(10) // 'messages generated by routine Query_HDF_Attr.' 

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Return
      End If 
      
      
c-------------------------------------------------------------------------------
c If core metadata exist in MODIS Aerosol Product, copy existing 
c MeasuredParameters to internal memory.
c-------------------------------------------------------------------------------

      If (core_exists) Then 
         rtn = Copy_MeasParm_MOD04(MET_Handles) 

         If (rtn .eq. FAIL) Then
            error_flag = .TRUE.

            msgbuf = 
     1      'Copy_MeasParm_MOD04 detected error copying MeasuredParameter data ' 
     2      // char(10) // 'from MOD04 product file. ' 
     3      // char(10) // 'Operator Action:  Refer to prior low level error LogStatus ' 
     4      // char(10) // 'messages generated by routine Copy_MeasParm_MOD04.' 

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If

      End If


c-------------------------------------------------------------------------------
c Set MeasuredParameters provided by current process 
c-------------------------------------------------------------------------------

      If ( NUM_NewMeasParm .gt. 0 ) Then

         rtn = Set_NewMeasParm_MOD04( MET_Handles, 
     1                                NUM_NewMeasParm,
     2                                Name_NewMeasParm,
     3                                PctMissing_NewMeasParm,
     4                                AutoFlag_NewMeasParm,
     5                                AutoFlagExp_NewMeasParm )

         If (rtn .eq. FAIL) Then
            error_flag = .TRUE.

            msgbuf =
     1      'Set_NewMeasParm_MOD04 detected error setting MeasuredParameters '
     2      // char(10) // 'by current process. '
     3      // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4      // char(10) // 'messages generated by routine Set_NewMeasParm_MOD04.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If

      End If


c-----If no error found, return SUCCEED = 0
      If (.not.error_flag) Set_MeasParm_MOD04 = SUCCEED


      Return

      END



c--------------------------------------------------------------------------------------
      Integer Function Set_NewMeasParm_MOD04( MET_Handles,
     1                                        NUM_NewMeasParm,
     2                                        Name_NewMeasParm,
     3                                        PctMissing_NewMeasParm,
     4                                        AutoFlag_NewMeasParm,
     5                                        AutoFlagExp_NewMeasParm )


      implicit none

      include 'mod04_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Set ECS MeasuredParameter groups computed and passed
C                by current process to internal memory.  A MOD04_L2 ECS
C                MeasuredParameter group contains the following 4 ECS
C                data model ODL (Object Description Language) objects:
C
C      "PARAMETERNAME"
C      "QAPERCENTMISSINGDATA"
C      "AUTOMATICQUALITYFLAG"
C      "AUTOMATICQUALITYFLAGEXPLANATION"
C
C
C !INPUT PARAMETERS:
C
C  character*(*) MET_Handles(*)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           MET_Handles is defined during initialization
C                           of the ECS Metadata Configuration File (MCF).
C                           In FORTRAN, element 1 of the MET_Handles
C                           array refers to the MCF file.  ECS inventory
C                           metadata are referenced as element 2 and
C                           archive metadata as element 3.
C
C  integer NUM_NewMeasParm  Variable containing the number of
C                           MeasuredParameters to be set by the current
C                           process.  If none, NUM_NewMeasParm=0.
C
C  character*(*) Name_NewMeasParm(*)
C                           Array containing the names of the
C                           MeasuredParameters to be set by the
C                           current process.  If none, be sure to pass
C                           NUM_NewMeasParm=0.
C
C                           A one-to-one association between the
C                           elements of arrays PctMissing_NewMeasParm,
C                           Name_NewMeasParm, AutoFlag_NewMeasParm and
C                           AutoFlagExp_NewMeasParm is assumed.
C
C  integer PctMissing_NewMeasParm(*)
C                           Array containing the percentages of missing
C                           data for the MeasuredParameters to be set
C                           by the current process.
C                           (See Name_NewMeasParm above).
C
C                           Note that percentage missing data are
C                           passed as integers values between 0-100.
C
C  character*(*) AutoFlag_NewMeasParm(*)
C                           Array of flags referring to the quality
C                           of the MeasuredParameters to be set by the
C                           current process.
C                           (See Name_NewMeasParm above).
C
C                           Valid quality values are:
C                           "Passed"/"Failed"/"Suspect".
C
C  character*(*) AutoFlagExp_NewMeasParm(*)
C                           Array containing the text explanations of the
C                           "criteria" used to define the AutomaticQualityFlags
C                           to be set by the current process.
C                           (See Name_NewMeasParm and AutoFlag_NewMeasParm
C                           above).
C
C
C !OUTPUT PARAMETERS:  None
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by         Richard Hucek
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        string_loc                (science code)
C
C
C    Named Constant:
C        BLANK                   (mod04_ECSMET.inc)
C        FAIL                    (mod04_ECSMET.inc)
C        FUNCNAME_PGS_MET_SET_i  (mod04_ECSMET.inc)
C        FUNCNAME_PGS_MET_SET_s  (mod04_ECSMET.inc)
C        INVENTORYMETADATA       (mod04_ECSMET.inc)
C        MCORE_*                 (mapi.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        NUM_MEAS_PARM_MOD04     (mod04_ECSMET.inc)
C        PGS_S_SUCCESS           (PGS_SMF.f)
C        SUCCEED                 (mod04_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

c PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'Set_NewMeasParm_MOD04')

c function argument declarations
      character*(*) MET_Handles(*),
     1              Name_NewMeasParm(*),
     2              AutoFlag_NewMeasParm(*),
     3              AutoFlagExp_NewMeasParm(*)

      integer NUM_NewMeasParm,
     1        PctMissing_NewMeasParm(*)


c other variable declarations
      character*25   Class,msg25
      character*60   Name
      character*128  AttrN,AttrV_s
      character*1024 msgbuf_SET_err
      character*2048 msgbuf

      integer i1,i2,i3,rtn,rtn_loc
      integer AttrV_i
      integer Class_Index
      integer pgs_met_setattr_s, pgs_met_setattr_i,String_Loc
      integer fbytea,fbyteb,fbytec,fbyted,fbytee,
     1        lbytea,lbyteb,lbytec,lbyted,lbytee

      logical Name_Is_Valid, error_flag, error_flag_group


C------------------------
C Initialization
C------------------------

      Set_NewMeasParm_MOD04 =  FAIL
      error_flag            = .FALSE.


c-------------------------------------------------------------------------------
c Loop on list of new Measured Parameters
c-------------------------------------------------------------------------------

      Do 100 i1 = 1, NUM_NewMeasParm
         Name    = Name_NewMeasParm(i1)
         rtn_loc = string_loc(Name,fbytea,lbytea)

c--------check validity of Measured Parameters name
         Name_Is_Valid = .FALSE.

         Do 200 i2 = 1, NUM_MEAS_PARM_MOD04

            If (Name_MeasParm_MOD04(i2) .eq. Name) Then
               Name_Is_Valid = .TRUE.
               Class_Index = i2
            End If

200      Continue

c--------invalid Measure Parameter Name
         If (.NOT. Name_Is_Valid) Then
            error_flag = .TRUE.

            If (Name .eq. BLANK) Then
               msgbuf = 
     1         'MeasuredParameter group name is blank - cannot be set.'
     2         // char(10) // 'Operator Action:  Notify SDST.'
            Else

               msgbuf =
     1         'MeasuredParameter group ' // Name(fbytea:lbytea) // ' is not on MOD04 list.'
     2         // char(10) // 'MeasuredParameter group ' // Name(fbytea:lbytea) // ' not set.'
     3         // char(10) // 'Operator Action:  Notify SDST.'

            End If


            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------valid Measure Parameter Name
         Else
            write(Class, '(I25)') Class_Index
            rtn_loc = string_loc(Class, fbyteb, lbyteb)


c-----------loop over 4 MeasuredParameter ODL objects used by MODIS atmosphere
            error_flag_group = .FALSE.

            Do 300 i3 = 1, 4

               If (i3 .eq. 1) Then   ! is PARAMETERNAME
                  AttrN = MCORE_PARAMETERNAME // '.' // Class(fbyteb:lbyteb)
                  AttrV_s = Name_NewMeasParm(i1)
                  rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_s)


               Else If (i3 .eq. 2) Then   ! is QAPERCENTMISSINGDATA
                  AttrV_i = PctMissing_NewMeasParm(i1)
                  AttrN   = MCORE_PERCENT_MISSING // '.' // Class(fbyteb:lbyteb)

                  write(msg25, '(I25)') AttrV_i

                  rtn_loc = string_loc(msg25,fbyted,lbyted)
                  rtn_loc = string_loc(AttrN,fbytee,lbytee)

                  If ( (AttrV_i.GE.0) .AND. (AttrV_i.LE.100) )Then
                     continue
                  Else
                     error_flag = .TRUE.

                     msgbuf =
     1               AttrN(fbytee:lbytee) // ' (' // msg25(fbyted:lbyted) // ') in MeasuredParameter '
     2               // 'group ' // Name(fbytea:lbytea)
     3               // char(10) // 'is out of range 0-100.  To be set anyway.'
     4               // char(10) // 'Operator Action:  Notify SDST.'

                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

                  rtn = PGS_MET_SetAttr_i(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_i)


               Else If (i3 .eq. 3) Then   ! is AUTOMATICQUALITYFLAG
                  AttrN   = MCORE_AUTO_QUALITY // '.' // Class(fbyteb:lbyteb)
                  rtn_loc = string_loc(AttrN,fbytee,lbytee)
                  AttrV_s = AutoFlag_NewMeasParm(i1)
                  rtn_loc = string_loc(AttrV_s,fbyted,lbyted)

                  If (AttrV_s.EQ.'Passed' .OR. AttrV_s.EQ.'Failed' .OR. AttrV_s.EQ.'Suspect') Then
                     Continue
                  Else
                     error_flag = .TRUE.

                     If (AttrV_s .EQ. BLANK) Then
                        msgbuf =
     1                  AttrN(fbytee:lbytee) // ' in MeasuredParameter group ' // Name(fbytea:lbytea)
     2                  // char(10) // 'is blank - not one of valid values:  Passed/Failed/Suspect.'
     3                  // 'To be set anyway. '
     4                  // char(10) // 'Operator Action:  Notify SDST.'
                     Else

                        msgbuf =
     1                  AttrN(fbytee:lbytee) // ' in MeasuredParameter group ' // Name(fbytea:lbytea)
     2                  // char(10) // 'is:  ' // AttrV_s(fbyted:lbyted)
     3                  // char(10) // '- not one of valid values:  Passed/Failed/Suspect.'
     4                  // 'To be set anyway.'
     5                  // char(10) // 'Operator Action:  Notify SDST.'
                     End If

                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

                  rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_s)


               Else   ! is AUTOMATICQUALITYFLAGEXPLANATION
                  AttrN = MCORE_AQUAL_FLG // '.' // Class(fbyteb:lbyteb)
                  rtn_loc = string_loc(AttrN,fbytee,lbytee)
                  AttrV_s = AutoFlagExp_NewMeasParm(i1)

                  If (AttrV_s .EQ. BLANK) Then
                     msgbuf =
     1               AttrN(fbytee:lbytee) // ' in MeasuredParameter group ' // Name(fbytea:lbytea)
     2               // char(10) // 'is blank - not a valid value.  To be set anyway.'
     3               // char(10) // 'Operator Action:  Notify SDST.'

                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

                  rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_s)

               End If

c--------------unable to set MeasuredParameter object - report to LogStatus
               If (rtn .ne. PGS_S_SUCCESS) Then
                  error_flag       = .TRUE.
                  error_flag_group = .TRUE.
                  rtn_loc          = string_loc(AttrN,fbytec,lbytec)

                  msgbuf_SET_err =
     1            ' unable to set ECS attribute '// AttrN(fbytec:lbytec)
     2            // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3            // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4            // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                                 msgbuf = FUNCNAME_PGS_MET_SET_s // msgbuf_SET_err
                  If (i3 .eq. 2) msgbuf = FUNCNAME_PGS_MET_SET_i // msgbuf_SET_err

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If

300         Continue

c-----------report success message to LogStatus
            If (.not. error_flag_group) Then
                msgbuf =
     1          'MeasuredParameter group ' // Name(fbytea:lbytea) // ' successfully set.'

                Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
            End If

         End If

100   Continue


c-----return SUCCEED = 0 if no errors found
      If (.not.error_flag) Set_NewMeasParm_MOD04 = SUCCEED


      Return

      END
