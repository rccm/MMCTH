      Integer Function Set_CoreMetadata(LRN_MCF,ExtGeoPntr_Flag,
     1                 No_InputPntr,LRN_InputPntr,Vrsn_InputPntr,
     2                 NUM_MeasParm,Name_MeasParm,PctMissing_MeasParm,
     3                 AutoFlag_MeasParm,AutoFlagExp_MeasParm,
     4                 NUM_PSA,Name_PSA,Value_PSA,
     5                 LUN_Of_NUM_INV_RP,MET_Handles)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C   Set values to the EOSDIS Core System (ECS) inventory metadata
C   fields that comprise ECS-HDF file attribute 'CoreMetadata.0'.  In
C   this context, set means to associate a value with an ECS attribute
C   name via function calls to the ECS Science Data Processing Toolkit.
C
C !DESCRIPTION:
C
C   The inventory metadata fields set by Set_CoreMetadata, the data
C   source for the field values, and Metadata Configuration File (MCF)
C   "Data Location" parameter are listed below.  Only a small number of
C   field values must actually be computed/determined by the science
C   code (referred to here as the PGE).  Most can be read in at run
C   time from the MODIS Geolocation product (many of the same fields
C   required by L2 products are carried in the upstream MOD03 product
C   granule) or as USER DEFINE RUNTIME PARAMETERs (RP) in the Process
C   Control File (PCF).  A few fields are static and can be predefined
C   within the MCF or even within the Set_CoreMetadata module itself.
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
C     34    AncillaryInputType.1       Set_CoreMetadata     PGE
C     35    AncillaryInputPointer.1           PCF           PGE
C     36    AssociatedPlatformShortName.1     PCF (RP)      PCF
C     37    AssociatedInstrumentShortName.1   PCF (RP)      PCF
C     38    AssociatedSensorShortName.1       PCF (RP)      PCF
C
C !INPUT PARAMETERS:
C  integer  LRN_MCF         Variable containing the PCF logical
C                           reference number (LRN) to the metadata
C                           configuration file (MCF)
C
C  logical ExtGeoPntr_Flag  Variable (.TRUE. or .FALSE. ) indicating
C                           whether a reference to the 1-km Geolocation
C                           (MOD03) product is to be written to product.
C                           If the product contains 1-km science
C                           parameters but not 1-km Geolocation data,
C                           then set ExtGeoPntr_Flag to .TRUE..
C                           Otherwise set it to .FALSE.
C
C  integer  No_InputPntr    Number of input datasets for granule.  This
C                           includes ancillary data files, look up
C                           tables, and other MODIS product files.
C                           System files, such as MCF, are not included.
C
C  integer  LRN_InputPntr(*)
C                           Array containing the PCF LRNs of all input
C                           datasets for granule (See No_InputPntr above).
C
C                           A one-to-one correspondence between the elements
C                           of the arrays LRN_InputPntr and Vrsn_InputPntr
C                           is assumed.
C
C  integer  Vrsn_InputPntr(*)
C                           Array containing the PCF version numbers of all
C                           input datasets for granule.  (See No_InputPntr
C                           and LRN_InputPntr above)
C
C  integer NUM_MeasParm     Variable containing the number of Measured
C                           Parameters to be set for granule.  If none,
C                           set NUM_MeasParm=0.
C
C  character*(*) Name_MeasParm(*)
C                           Array containing the names of all Measured
C                           Parameters to be set for granule.
C
C                           A one-to-one correspondence between the
C                           elements of arrays PctMissing_MeasParm,
C                           Name_MeasParm, AutoFlag_MeasParm and
C                           AutoFlagExp_MeasParm is assumed.
C
C  integer PctMissing_MeasParm(*)
C                           Array containing the percentage of missing
C                           data for each Measured Parameter to be set
C                           for granule (See Name_MeasParm above).
C
C                           Note that the percentage missing data are
C                           passed as integers to Set_CoreMetadata.
C
C  character*(*) AutoFlag_MeasParm(*)
C                           Array containing a "simple" assessment of
C                           the quality of each Measured Parameter
C                           value to be set for granule.
C
C                           Valid values for the assessment are
C                           "Passed"/"Suspect"/"Failed".
C                           (See also Name_MeasParm above).
C
C  character*(*) AutoFlagExp_MeasParm(*)
C                           Array containing explanations of the "criteria"
C                           used to define the AutomaticQualityFlags
C                           set by the current process.
C                           (See Name_MeasParm and AutoFlag_MeasParm
C                           above).
C
C  integer NUM_PSA          Variable containing the number of
C                           Product-Specific Attributes (PSAs) to be
C                           set for granule.
C
C  character*(*) Name_PSA(*)
C                           Array containing the names of all PSAs
C                           to be set for granule
C
C                           A one-to-one correspondence between the
C                           elements of arrays Name_PSA and Value_PSA
C                           is assumed.
C
C  Real Value_PSA(*)        Array containing the values of all PSAs
C                           to be set for granule (See Name_PSA
C                           above).  Values are restricted to range
C                           plus (+) and minus (-) 99999.99 and are
C                           set to the product in F8.2 format;  A
C                           value of -99999.0 is interpreted as a
C                           flag meaning the PSA is not to be set.
C
C  integer LUN_Of_NUM_INV_RP
C                           PCF RP LUN containing the number of
C                           "Inventory" RP metadata objects to be set
C                           by current process.  A list of these
C                           metadata names and values is expected to
C                           follow immediately in the PCF.
C
C !OUTPUT PARAMETERS:
C  character MET_Handles(20)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           In FORTRAN, element 1 of the MET_Handles
C                           array is reserved for the MCF file.   ECS
C                           inventory metadata is referenced as as
C                           element 2 and archive metadata as element 3.
C
C
C !REVISION HISTORY:
c revision 1.5  2001/03/30  rhucek
c renamed subprogram "Get_LUN_Of_LocalVrsnID_Atmos" to "Get_LUN_Of_LclVrsnID_Atmos" 
c
c Revision 1.4  1999/02/17  20:06:30  rhucek
c updated prolog.
c
c Revision 1.3  1999/02/03  21:51:11  fhliang
c added a new argument LUN_Of_NUM_INV_RP.
c
c Revision 1.2  1999/02/03  15:56:16  fhliang
c used routine Get_LUN_Of_LclGranID_Atmos to get LUN_Of_LclVrsnID used in
c routine Set_LclGranID.
c
c Revision 1.1  1999/01/29  16:58:11  fhliang
c Initial revision
c
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C
C !DESIGN NOTES:
C
C   Function Set_CoreMetadata sets (to internal memory) ECS metadata
C   required for archival and queries on HDF product files stored
C   at the DAACs.  No metadata fields are actually written to the product
C   file during the call to Set_CoreMetadata;  only the association of
C   values to fields is made at this time.  Actual insertion of attribute
C   'CoreMetadata.0' into the HDF product file takes place later with a
C   call to the SDPTK routine PGS_MET_Write.  This call is made from
C   outside of Set_CoreMetadata, and only after any and all other core
C   metadata fields have been set by the PGE.  If the MODIS HDF
C   Application Programming Interface (M-API) Utility is being used,
C   PGS_MET_Write is called automatically when the HDF product file is
C   closed with CPMFIL.  If one does not use M-API, then a direct call
C   to PGS_MET_Write is required to attach 'CoreMetadata.0' to the file
C   before it is closed.
C
C  Returns:     SUCCEED if successful, FAIL if an error occurs
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_Init                 (libPGSTK.a)
C        Get_LUN_Of_LclVrsnID_Atmos   (science code)
C        MODIS_SMF_SETDYNAMICMSG      (science code)
C        Set_AncInputPntr             (science code)
C        Set_GeoData_Atmos            (science code)
C        Set_InputPntr_Atmos          (science code)
C        Set_InvPSA_Atmos             (science code)
C        Set_LclGranID                (science code)
C        Set_MeasParm_Atmos           (science code)
C        Set_RP_Data_Atmos            (science code)
C
C    Named Constant:
C        FAIL                         (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC              (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS                (PGS_SMF.f)
C        SUCCEED                      (Atmos_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Set_CoreMetadata')


C other variable declarations
      integer  LRN_MCF, NUM_MeasParm, No_InputPntr, NUM_PSA
      integer  LRN_InputPntr(*),Vrsn_InputPntr(*), PctMissing_MeasParm(*)

      character*(*) AutoFlag_MeasParm(*),AutoFlagExp_MeasParm(*),
     1              Name_MeasParm(*),Name_PSA(*)

      integer       LUN_Of_NUM_INV_RP
      character*(*) MET_Handles(*)
      logical       ExtGeoPntr_Flag
      real          Value_PSA(*)

      integer       LUN_Of_LclVrsnID

      integer  pgs_met_init

      integer  Set_AncInputPntr,
     1         Set_GeoData_Atmos,
     2         Set_InputPntr_Atmos,
     3         Set_InvPSA_Atmos,
     4         Set_LclGranID,
     5         Set_MeasParm_Atmos,
     6         Set_RP_Data_Atmos,
     7         string_loc

      integer  rtn, rtn_flag, rtn_loc
      integer  fbyte_LUN, lbyte_LUN

      character*512 msgbuf
      character*25  msg25_LUN

      logical error_flag


C------------------------
C Initialization
C------------------------

      Set_CoreMetadata = FAIL
      error_flag = .false.

C-------------------------------------------------------------------
C Initialize metadata tool defining array MET_Handles, and set
C inventory attribute fields whose values are provided in the MCF.
C-------------------------------------------------------------------

      rtn = pgs_met_init(LRN_MCF,MET_Handles)

      If (rtn.ne.PGS_S_SUCCESS) Then
         msgbuf = 'pgs_met_init unable to initialize MCF.'
     1     // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2     // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3     // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Return
      End If

C---------------------------------------------------------------------
C Set ECS inventory attribute fields whose values are to be retrieved
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

      rtn = Set_GeoData_Atmos(Met_Handles)

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_GeoData_Atmos detected error setting ECS Geo product attribtutes.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_GeoData_Atmos.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C----------------------------------------------------------------------
C Set ECS attributes whose values are retrieved as USER DEFINED
C RUNTIME PARAMETERs from PCF.  The following fields are included:
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
     1                         LUN_Of_NUM_INV_RP,
     2                         INVENTORYMETADATA )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_RP_Data_Atmos detected error setting ECS inventory RP data.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_RP_Data_Atmos.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C--------------------------------------------------------------------
C Set ECS inventory metadata attributes whose values are passed in
C as actual function arguments from PGE.  These include the Measured
C Parameter fields:
C
C "ParameterName"
C "AutomaticQualityFlag" (optional)
C "AutomaticQualityFlagExplanation" (optional)
C "QAPercentMissingData"
C--------------------------------------------------------------------

      rtn = Set_MeasParm_Atmos( MET_Handles,
     1                          NUM_MeasParm,
     2                          Name_MeasParm,
     3                          PctMissing_MeasParm,
     4                          AutoFlag_MeasParm,
     5                          AutoFlagExp_MeasParm )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_MeasParm_Atmos detected error setting ECS MeasuredParameters.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_MeasParm_Atmos.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C-----------------------------------------------------------------------
C  Set Product Specific Attributes (PSAs).
C-----------------------------------------------------------------------

      rtn = Set_InvPSA_Atmos( NUM_PSA,
     1                        Name_PSA,
     2                        Value_PSA,
     3                        MET_Handles )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_InvPSA_Atmos detected error setting ECS PSAs.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_InvPSA_Atmos.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C-----------------------------------------------------------------
C Set "LocalGranuleID".  This is the product file name in the
C MODIS file naming convention.
C-----------------------------------------------------------------

      Call Get_LUN_Of_LclVrsnID_Atmos( LUN_Of_NUM_INV_RP,
     1                                 LUN_Of_LclVrsnID,
     2                                 rtn_flag)

      If (rtn_flag .eq. FAIL) Then
         error_flag = .true.

         write(msg25_LUN,'(I25)') LUN_Of_NUM_INV_RP
         rtn_loc = string_loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

         msgbuf = 'Get_LUN_Of_LclVrsnID_Atmos detected error identifying '
     1   // char(10) // 'LOCALVERSIONID LUN in RP group beginning at LUN='
     2   // msg25_LUN(fbyte_LUN:lbyte_LUN) // '.'
     3   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4   // char(10) // 'messages originating in routine Get_LUN_Of_LclVrsnID_Atmos.'

         Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )

      Else
         rtn = Set_LclGranID( Met_Handles, LUN_Of_LclVrsnID )

         If (rtn .eq. FAIL) Then
            error_flag = .true.

            msgbuf = 'Set_LclGranID detected error setting ECS LocalGranuleID.'
     1      // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2      // char(10) // 'messages originating in function Set_LclGranID.'

            Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
         EndIf

      EndIf



C---------------------------------------------------------------------
C Set Universal References (URs) or "InputPointers" to input data set.
C---------------------------------------------------------------------

      rtn = Set_InputPntr_Atmos( No_InputPntr,
     1                           LRN_InputPntr,
     2                           Vrsn_InputPntr,
     3                           MET_Handles )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_InputPntr_Atmos detected error setting ECS INPUTPOINTER attribute.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating in function Set_InputPntr_Atmos.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C-----------------------------------------------------------------
C Set attribute field "AncillaryInputPointer"
C-----------------------------------------------------------------

      If (ExtGeoPntr_Flag) rtn = Set_AncInputPntr(Met_Handles)

      If (rtn .eq. FAIL) Then
         error_flag = .TRUE.

         msgbuf = 'Set_AncInputPntr detected error setting ECS '
     1      // char(10) // 'AncillaryInputGranule metadata group. '
     2      // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     3      // char(10) // 'messages originating in function Set_AncInputPntr. '


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



      If (.not.error_flag) Set_CoreMetadata = SUCCEED

      Return

      End
