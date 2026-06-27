      Integer Function Update_InvMet_MOD06(LUN_MCF,ExtGeoPntr_Flag,
     1                 No_PCF_Attr,LUN_PCF_Attr,Name_PCF_Attr,
     2                 No_InputPntr,LUN_InputPntr,Vrsn_InputPntr,
     3                 No_MeasParm,Name_MeasParm,PctMissing_MeasParm,
     4                 AutoFlag_MeasParm,AutoFlagExp_MeasParm,
     5                 No_PSA,Name_PSA,Value_PSA,
     6                 MET_Handles)

      implicit none
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Update ECS and user-defined product specific "Inventory" metadata 
C   previously written to the MODIS L2 cloud product (MOD06_L2) by one or 
C   more prior processes.  In the context used here, set means to associate 
C   a value with a metadata parameter name via function calls to the ECS 
C   Science Data Processing Toolkit (SDPTK).
C
C   The specific metadata fields updated by Update_InvMet_MOD06, their
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
C     34    AncillaryInputType.1       Update_InvMet_MOD06  PGE
C     35    AncillaryInputPointer.1           PCF           PGE
C     36    AssociatedPlatformShortName.1     PCF (RP)      PCF
C     37    AssociatedInstrumentShortName.1   PCF (RP)      PCF
C     38    AssociatedSensorShortName.1       PCF (RP)      PCF
C
C !INPUT PARAMETERS:
C  integer  LUN_MCF         Variable containing the PCF logical 
C                           reference number (LUN) to the metadata 
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
C  integer No_PCF_Attr      The number of ECS inventory fields with 
C                           values defined as USER DEFINED RUNTIME 
C                           PARAMETERS in the PCF file.  Typical fields 
C                           that may be set in the PCF are:  
C                           ReprocessingActual, 
C                           ReprocessingPlanned, 
C                           LocalVersionID, 
C                           PGEVersionID, 
C                           AssociatedPlatformShortName,
C                           AssociatedInstrumentShortName,
C                           AssociatedSensorShortName.
C
C                           The LocalVersionID is used to distinguish
C                           different processing versions of the same
C                           data granule.  By MODIS convention, the 
C                           version ID is to be represented by a 
C                           3-character string consisting only of 
C                           numeric characters (e.g. '101' or
C                           '001').    
C                           
C
C  integer LUN_PCF_Attr(*)  Array containing the PCF logical reference
C                           numbers (LUN) of ECS inventory field values 
C                           listed as USER DEFINED RUNTIME PARAMETERS. 
C
C                           A one-to-one correspondence between the 
C                           elements of the arrays LUN_PCF_Attr and
C                           Name_PCF_Attr is assumed.
C 
C  character*(*) Name_PCF_Attr(*)
C                           Array containing the names of ECS inventory 
C                           metadata fields whose values are listed as 
C                           USER DEFINED RUNTIME PARAMETERS in the PCF.
C                           Note:  these names are not retrieved from 
C                           the PCF as are their associated values; they 
C                           must be assigned by the science code. (See
C                           LUN_PCF_Attr above) 
C                     
C  integer  No_InputPntr    Number of input datasets required by process.  
C                           This includes ancillary data files, look up 
C                           tables, and other MODIS product files.  
C                           System files, such as MCF, are not included.
C
C  integer  LUN_InputPntr(*)
C                           Array containing the PCF LUNs of input 
C                           datasets used by process (See No_InputPntr above). 
C
C                           A one-to-one correspondence between the elements 
C                           of the arrays LUN_InputPntr and Vrsn_InputPntr
C                           is assumed.
C
C  integer  Vrsn_InputPntr(*)  
C                           Array containing the PCF version numbers of 
C                           input datasets used by process.  (See No_InputPntr 
C                           and LUN_InputPntr above) 
C
C  integer No_MeasParm      Variable containing the number of additional 
C                           Measured Parameters to be set by process.  If 
C                           none, No_MeasParm=0.
C  
C  character*(*) Name_MeasParm(*)  
C                           Array containing the names of additional
C                           Measured Parameters to be set by process. 
C
C                           A one-to-one correspondence between the 
C                           elements of arrays PctMissing_MeasParm,
C                           Name_MeasParm, AutoFlag_MeasParm and
C                           AutoFlagExp_MeasParm is assumed.
C
C  integer PctMissing_MeasParm(*)
C                           Array containing the percentage of missing 
C                           data for each additional Measured Parameter to 
C                           be set by process (See Name_MeasParm above). 
C                         
C                           Note that the percentage missing data are 
C                           passed to Update_InvMet_MOD06 as integers. 
C                            
C  character*(*) AutoFlag_MeasParm(*)  
C                           Array containing a "simple" assessment of 
C                           the quality of each additional Measured 
C                           Parameter value to be set by process 
C
C                           Valid values are Passed/Failed/Suspect.  
C                           If the value of AutoFlag_MeasParm is blank 
C                           (' '), neither AutoFlag_MeasParm nor 
C                           AutoFlagExp_MeasParm will be set by function 
C                           Update_InvMet_MOD06 (See Name_MeasParm above).  
C
C  character*(*) AutoFlagExp_MeasParm(*)
C                           Array containing the "criterion" for the
C                           quality judgement made on each Measured 
C                           Parameter value to be set by process.  If 
C                           the value of the AutoFlag_MeasParm is set to 
C                           blank (' '), the AutoFlagExp_MeasParm will not 
C                           be set (see AutoFlag_MeasParm and 
C                           Name_MeasParm above).
C
C  integer No_PSA           Variable containing the number of additional 
C                           Product-Specific Attributes (PSAs) to be
C                           set by process. 
C 
C  character*(*) Name_PSA(*)   
C                           Array containing the names of additional PSAs
C                           to be set by process. 
C
C                           A one-to-one correspondence between the 
C                           elements of arrays Name_PSA and Value_PSA 
C                           is assumed.
C
C  Real Value_PSA(*)        Array containing the values of additional 
C                           PSAs to be set by process (See Name_PSA
C                           above).  Values are restricted to range
C                           plus (+) and minus (-) 99999.99 and are
C                           set to the product in F8.2 format;  A 
C                           value of -99999.0 is interpreted as a 
C                           flag meaning the PSA is not to be set. 
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
C        PGS_MET_Init              (libPGSTK.a)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        PGS_MET_GetSetAttr_*      (libPGSTK.a)
C        PGS_PC_GetConfigData      (libPGSTK.a)
C        PGS_PC_GetUniversalRef    (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        Gen_Modis_FileName        (science code)
C        string_loc                (science code)
C
C
C    Named Constant:
C        MCORE_*               (mapi.inc)
C        MECS_CORE             (mapi.inc)
C        MODIS_E_GENERIC       (PGS_MODIS_39500.f)
C        MODIS_W_GENERIC       (PGS_MODIS_39500.f)
C        PGSd_PC_*             (PGS_PC.f)
C        PGS_S_SUCCESS         (PGS_SMF.f)
C    
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Update_InvMet_MOD06')

      integer       NUM_MEAS_PARM_MOD06
      parameter    (NUM_MEAS_PARM_MOD06 = 1)

      integer       MAX_NO_PNTRS
      parameter    (MAX_NO_PNTRS=1000)

      integer       FAIL
      parameter    (FAIL=-1)

      integer       SUCCEED
      parameter    (SUCCEED = 0)

      integer       LUN_GEO
      parameter    (LUN_GEO=600000)

      integer       INVENTORYMETADATA
      parameter    (INVENTORYMETADATA = 2 )

C other variable declarations
      character*8   ESDT_Name

      character*25  msg25,msg25b,msg25c,msg25_LUN,msg25_VRSN
      character*128  fname
      character*128 AttrN,AttrV,Field_Name(8),buf_char
      character*512 msgbuf, runtimepar,parmname
      character*(*) MET_Handles(*)
      character*(PGSd_PC_UREF_LENGTH_MAX) UR_Ref(MAX_NO_PNTRS)
      character*(PGSd_PC_VALUE_LENGTH_MAX) RunTimeParm
      character*(PGSd_PC_VALUE_LENGTH_MAX) LOCAL_VRSNID


      character*(*) AutoFlag_MeasParm(*),AutoFlagExp_MeasParm(*),
     2              Name_MeasParm(*),Name_PSA(*),Name_PCF_Attr(*)

      integer    i,No_MeasParm,No_PCF_Attr,No_InputPntr,No_PSA,
     2           icounter,rtn,rtn_loc,LUN_MCF

      integer    buf_int(4),
     2           LUN_PCF_Attr(*),LUN_InputPntr(*),Vrsn_InputPntr(*),
     3           PctMissing_MeasParm(*),LUN_LCLVRNSID,VRSN_GEO

      integer    pgs_met_init,pgs_met_setattr_d,
     2           pgs_met_setattr_s, pgs_met_setattr_i,
     3           PGS_MET_GetPCAttr_i,PGS_MET_GetPCAttr_d,
     4           PGS_MET_GetPCAttr_s,PGS_PC_GetConfigData,
     6           PGS_PC_GetUniversalRef,pgs_met_getsetattr_s,
     7           Update_PSA_MOD06,Update_InputPntr_MOD06,
     8           Update_MeasParm_MOD06,Gen_Modis_FileName,String_Loc

      integer    fbyte,fbyteb,fbytec,fbyte_LUN,fbyte_VRSN,
     1           lbyte,lbyteb,lbytec,lbyte_LUN,lbyte_VRSN

      real   Value_PSA(*)

      double precision buf_dbl(4)

      logical    error_flag,ExtGeoPntr_Flag,loopexit
      
C------------------------
C Initialization
C------------------------

      Update_InvMet_MOD06 = FAIL 
      error_flag = .false.
      fname = ' '

      Do 10 i = 1, 4
         buf_int(i) = 0
         buf_dbl(i) = 0.0D0
   10 continue

      VRSN_GEO = 1
c-----Set up variables for status messaging
      write(msg25_LUN,'(I25)') LUN_GEO
      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

      write(msg25_VRSN,'(I25)') VRSN_GEO
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)


C-------------------------------------------------------------------
C Initialize metadata tool defining array MET_Handles, and set 
C inventory attribute fields whose values are provided in the MCF.
C-------------------------------------------------------------------

      rtn = pgs_met_init(LUN_MCF,MET_Handles)

      If (rtn.ne.PGS_S_SUCCESS) Then
         msgbuf = 'pgs_met_init unable to initialize MCF'
     2      // char(10) // 'Operator Action:  Check for valid MCF file. If wrong or corrupted,'
     3      // char(10) // 'stage correct MCF and rerun PGE. Otherwise, notify SDST.'


         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'FUNCNAME')

         Return
      End If

C-----------------------------------------------------------------------
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

      Field_Name(1) = MCORE_EAST_BOUND
      Field_Name(2) = MCORE_NORTH_BOUND
      Field_Name(3) = MCORE_SOUTH_BOUND
      Field_Name(4) = MCORE_WEST_BOUND

      Do 20 i = 1, 4
         AttrN = Field_Name(i)

         rtn = PGS_MET_GetPCAttr_d(LUN_GEO,VRSN_GEO,MECS_ARCHIVE,AttrN,buf_dbl)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.
            rtn = String_Loc(AttrN,fbyte,lbyte)

            msgbuf = 'PGS_MET_GetPCAttr_d unable to retrieve ECS attribute '//AttrN(fbyte:lbyte) 
     2         // char(10) // 'from Geolocation product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     3         //  ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     4         // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     5         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     6         // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7         // char(10) // 'fault is identified, stage correct PCF/input file and '
     8         // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

            call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

         Else
            rtn = PGS_MET_SetAttr_d(MET_Handles(INVENTORYMETADATA),AttrN,buf_dbl)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)

               msgbuf = 'PGS_MET_SetAttr_d unable to set ECS attribute ' //AttrN(fbyte:lbyte) 
     1          // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2          // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3          // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'


               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
            End If

         End If

   20 Continue

      Field_Name(1) = MCORE_DAYNIGHTFLAG
      Field_Name(2) = MCORE_RANGE_BEG_TIME
      Field_Name(3) = MCORE_RANGE_ENDING_TIME
      Field_Name(4) = MCORE_RANGE_BEG_DATE
      Field_Name(5) = MCORE_RANGE_ENDING_DATE
      Field_Name(6) = MCORE_EQUATCROSSINGDATE//'.1'
      Field_Name(7) = MCORE_EQUATCROSSINGTIME//'.1'

      Do 25 i = 1, 7
         buf_char = ' '
         AttrN = Field_Name(i)

         rtn = PGS_MET_GetPCAttr_s(LUN_GEO,1,MECS_CORE,AttrN,buf_char)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.
            rtn = String_Loc(AttrN,fbyte,lbyte)

            msgbuf = 'PGS_MET_GetPCAttr_s unable to retrieve ECS attribute '// AttrN(fbyte:lbyte)
     1         // char(10) // 'from Geolocation product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2         //  ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3         // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     4         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5         // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6         // char(10) // 'fault is identified, stage correct PCF/input file and '
     7         // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

            call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

         Else
            rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,buf_char)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)

               msgbuf = 'PGS_MET_SetAttr_s unable to set ECS attribute ' //AttrN(fbyte:lbyte)
     2            // char(10) // 'Operator Action:  Check for valid MCF file. If wrong or corrupted,'
     3            // char(10) // 'stage correct MCF and rerun PGE. Otherwise, notify SDST.'

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
            End If

         End If

   25 Continue


      AttrN = MCORE_ORBIT_NUM//'.1'

      rtn = PGS_MET_GetPCAttr_i(LUN_GEO,1,MECS_CORE,AttrN,buf_int)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.
         rtn = String_Loc(AttrN,fbyte,lbyte)

         msgbuf = 'PGS_MET_GetPCAttr_i unable to retrieve ECS attribute '//AttrN(fbyte:lbyte)
     1     // char(10) // 'from Geolocation product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2     //  ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3     // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     4     // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5     // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6     // char(10) // 'fault is identified, stage correct PCF/input file and '
     7     // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

      Else
         rtn = PGS_MET_SetAttr_i(MET_Handles(INVENTORYMETADATA),AttrN,buf_int)

         If (rtn.ne.PGS_S_SUCCESS) Then
             error_flag = .true.
             rtn = String_Loc(AttrN,fbyte,lbyte)

             msgbuf = 'PGS_MET_SetAttr_i unable to set ECS attribute ' //AttrN(fbyte:lbyte)
     2            // char(10) // 'Operator Action:  Check for valid MCF file. If wrong or corrupted,'
     3            // char(10) // 'stage correct MCF and rerun PGE.  Otherwise, notify SDST.'

             call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
         End If

      End If


      AttrN = MCORE_EQUATCROSSINGLONG//'.1'

      rtn = PGS_MET_GetPCAttr_d(LUN_GEO,1,MECS_CORE,AttrN,buf_dbl)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.
         rtn = String_Loc(AttrN,fbyte,lbyte)

         msgbuf = 'PGS_MET_GetPCAttr_d unable to retrieve ECS attribute '//AttrN(fbyte:lbyte) // char(10)
     1       // 'from Geolocation product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2       //  ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3       // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     4       // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5       // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6       // char(10) // 'fault is identified, stage correct PCF/input file and '
     7       // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
      Else
         rtn = PGS_MET_SetAttr_d(MET_Handles(INVENTORYMETADATA),AttrN,buf_dbl)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.
            rtn = String_Loc(AttrN,fbyte,lbyte)

            msgbuf = 'PGS_MET_SetAttr_d unable to set ECS attribute ' // AttrN(fbyte:lbyte)
     1            // char(10) // 'Operator Action:  Check for valid MCF file. If wrong or corrupted,'
     2            // char(10) // 'stage correct MCF and rerun PGE. Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
         End If

      End If

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

      Do 30 icounter = 1, No_PCF_Attr
         AttrN = Name_PCF_Attr(icounter)

C........Next 6 lines for messaging support 
         rtn = String_Loc(AttrN,fbyte,lbyte)
         write(msg25b,'(I25)') LUN_PCF_Attr(icounter)
         rtn = string_loc(msg25b,fbyteb,lbyteb)

         write(msg25c,'(I25)') icounter 
         rtn = string_loc(msg25c,fbytec,lbytec)

         If (AttrN .eq. ' ') Then
            error_flag = .true.

            msgbuf = 'Input array element ' // msg25c(fbytec:lbytec) 
     1       // ' containing ECS attribute name is blank.' 
     2       // char(10) // 'Update_InvMet_MOD06 unable to set unknown attribute '
     3       // char(10) // 'Operator Action:  Notify SDST.'
                     
            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
         Else

            rtn = PGS_PC_GetConfigData(LUN_PCF_Attr(icounter),RunTimeParm)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf = 'PGS_PC_GetConfigData unable to retrieve value of ECS attribute ' 
     1          // char(10) // AttrN(fbyte:lbyte) // ' on Runtime Parameter LUN = ' 
     2          // msg25b(fbyteb:lbyteb)
     3          // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     4          // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     5          // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     6          // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'


               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
            Else

C..............ECS tools allow setting a "blank" value
               If (RunTimeParm .eq. ' ') Then
                  error_flag = .true.

                  msgbuf = 'Run time parameter value on LUN = '//msg25b(fbyteb:lbyteb)//' is blank ' 
     1                  // char(10) // 'Attribute name is ' // AttrN(fbyte:lbyte) 
     2                  // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter '
     3                  // char(10) // 'value.  If blank, stage correct PCF and rerun PGE. '
     4                  // char(10) // 'Otherwise, notify SDST.'


                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If
        
               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,RunTimeParm)

               If (rtn.ne.PGS_S_SUCCESS) Then
                  error_flag = .true.

                  msgbuf = 'PGS_MET_SetAttr_s unable to set ECS attribute ' // AttrN(fbyte:lbyte) 
     1             // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2             // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3             // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'


                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
               End If

            End If 

         End If

   30 Continue

C--------------------------------------------------------------------
C Set ECS inventory metadata attributes whose values are passed in 
C as function arguments from PGE.  These include the Measured 
C Parameter fields:
C
C "ParameterName"
C "AutomaticQualityFlag" (optional)
C "AutomaticQualityFlagExplanation" (optional)
C "QAPercentMissingData"
C
C and the Product Specific Attribute (PSA) fields:
C
C "AdditionalAttributeName"
C "ParameterValue",
C--------------------------------------------------------------------

C-----------------------------------------------------------------
C Reset MOD06_L2 Measured Parameters 
C-----------------------------------------------------------------

      rtn = Update_MeasParm_MOD06(MET_Handles,NUM_MEAS_PARM_MOD06)

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Update_MeasParm_MOD06 detected error setting ECS '
     1    // char(10) // 'Measured Parameter attribute group.'
     2    // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3    // char(10) // 'messages originating within function Update_MeasParm_MOD06.'

      
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


C-----------------------------------------------------------------
C Update MOD06 Product Specific Metadata
C-----------------------------------------------------------------

      rtn = Update_PSA_MOD06(MET_Handles,No_PSA,Name_PSA,Value_PSA)

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 
     1    'Update_PSA_MOD06 detected error setting ECS Product Specific Attributes (PSAs). '
     2    // char(10) // 'Operator Action:  Refer to prior low level LogStatus messages '
     3    // char(10) // 'originating within function Update_PSA_MOD06.'
         
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


C-----------------------------------------------------------------------
C Set "LocalGranuleID".  This is the product file name in the
C MODIS file naming convention.
C-----------------------------------------------------------------------

      i=1
      loopexit=.FALSE.
      LUN_LCLVRNSID=-1

      Do While(i.le.No_PCF_Attr .and. loopexit.eqv..FALSE.)
          runtimepar=Name_PCF_Attr(i)

          If(runtimepar(1:14).eq.'LOCALVERSIONID') Then
              LUN_LCLVRNSID=LUN_PCF_Attr(i)
              loopexit=.TRUE.
          End If

          i=i+1
      End Do


*/ Retrieve the value of SHORTNAME from the MCF.

      parmname=MCORE_SHORT_NAME
      rtn=pgs_met_getsetattr_s(MET_Handles(INVENTORYMETADATA),parmname,ESDT_Name)
 
      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf =
     1   'pgs_met_getsetattr_s unable to read the ESDT shortname from MCF.'
     2     // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3     // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4     // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      else


C .......Retrieve value of LOCALVERSIONID from PCF runtime parameter.
         rtn=pgs_pc_getconfigdata(LUN_LCLVRNSID,LOCAL_VRSNID)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .true.
            write(msg25,'(I25)') LUN_LCLVRNSID
            rtn = string_loc(msg25,fbyte,lbyte)

            msgbuf =
     1      'pgs_pc_getconfigdata unable to read LOCALVERSIONID'//' on LUN '// msg25(fbyte:lbyte)
     2      // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     3      // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     4      // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'



            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
         Else 
c ..........Construct LOCALGRANULEID 
            rtn = Gen_Modis_FileName(ESDT_Name,LOCAL_VRSNID,fname)

            If (rtn.ne.SUCCEED) Then
               error_flag = .true.
               msgbuf = 'Function Gen_Modis_FileName failed.'
     1                  // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     2                  // char(10) // 'messages originating within function Gen_Modis_FileName.'


               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            Else

c .............Set LOCALGRANULEID
               AttrN = MCORE_LOCALGRANULEID

               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,fname)

               If (rtn .ne. PGS_S_SUCCESS) Then
                  error_flag = .true.
                  rtn = String_Loc(AttrN,fbyte,lbyte)

                  msgbuf = 'PGS_MET_SetAttr_s unable to set '// AttrN(fbyte:lbyte)
     1             // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2             // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3             // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'



                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If
            End If
         End If
      End If


C-----------------------------------------------------------------
C Update MOD06 InputPointer 
C-----------------------------------------------------------------

      rtn = Update_InputPntr_MOD06(MET_Handles,No_InputPntr,LUN_InputPntr,VRSN_InputPntr)

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Update_InputPntr_MOD06 detected error incrementing ECS '
     1            // char(10) // 'INPUTPOINTER metadata attribute. ' 
     2            // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3            // char(10) // 'messages originating within function Update_InputPntr_MOD06.'

         
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C-----------------------------------------------------------------
C Set attribute field "AncillaryInputPointer" 
C-----------------------------------------------------------------

      If (ExtGeoPntr_Flag) Then
          AttrN = MCORE_ANCIL_INPUT_TYPE//'.1'
          AttrV = 'Geolocation'

          rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV)

          If (rtn.ne.PGS_S_SUCCESS) Then
              error_flag = .true. 
              rtn = String_Loc(AttrN,fbyte,lbyte)

              msgbuf = 'PGS_MET_SetAttr_s unable to set '// AttrN(fbyte:lbyte) 
     1         // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2         // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3         // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

              Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
          End If

c ---     Only 1 geolocation file expected 

          rtn = PGS_PC_GetUniversalRef(LUN_GEO,VRSN_GEO,UR_Ref(1))

          If (rtn .ne. PGS_S_SUCCESS) Then
             error_flag = .true.

             msgbuf = 'PGS_PC_GetUniversalRef unable to retrieve '
     1          // 'Geolocation UR.' // char(10) // 'LUN = ' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     2          // ' Version_No = '//msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3          // char(10) // 'Operator Action:  Check PCF Universal Reference Identifier. '
     4          // char(10) // 'If invalid, stage correct PCF and rerun PGE, Otherwise, notify SDST.'


             Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

          Else
             AttrN = MCORE_ANCIL_POINTER//'.1'

             rtn = PGS_MET_SetAttr_s (MET_Handles(INVENTORYMETADATA),AttrN,UR_Ref(1))

             If (rtn.ne.PGS_S_SUCCESS) Then
                 error_flag = .true. 
                 rtn = String_Loc(AttrN,fbyte,lbyte)

                 msgbuf = 'PGS_MET_SetAttr_s unable to set '// AttrN(fbyte:lbyte)
     1            // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2            // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3            // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                 call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
             End If

          End If
      End If

      If (.not.error_flag) Update_InvMet_MOD06 = SUCCEED 

      Return

      End

C-----------------------------------------------------------------------

      Integer Function Update_InputPntr_MOD06(MET_Handles, No_InputPntr,
     +                                        LUN_InputPntr, VRSN_InputPntr)

      implicit none
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Increment and set an updated list of process input files
C                to be stored as MOD06 product Inventory metadata (in 
C                attribute INPUTPOINTER) 
C
C
C !INPUT PARAMETERS:
C  character MET_Handles(20)
C                           Array containing the names of "MASTER"
C                           groups defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           In FORTRAN, element 1 of the MET_Handles
C                           array is reserved for the MCF file.   ECS
C                           inventory metadata is referenced as as
C                           element 2 and archive metadata as element 3.
C
C  integer  No_InputPntr    Variable containing the number of additional 
C                           input datasets to be set for granule.  
C                           It includes ancillary data files, look up
C                           tables, and other MODIS product files.
C                           System files, such as MCF, are not included.
C
C  integer  LUN_InputPntr(*)
C                           Array containing the PCF LUNS of additional files
C                           used during granule processing (See No_InputPntr 
C                           above).
C
C                           A one-to-one correspondence between the elements
C                           of the arrays LUN_InputPntr and VRSN_InputPntr
C                           is assumed.
C
C  integer  VRSN_InputPntr(*)
C                           Array containing PCF version numbers of additional
C                           input files used during granule processing.  
C                           (See No_InputPntr and LUN_InputPntr above)
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
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Update_InputPntr_MOD06')

      integer INVENTORYMETADATA
      Parameter (INVENTORYMETADATA = 2)

      integer MAX_NO_PNTRS
      parameter (MAX_NO_PNTRS=50)

      integer LUN_MOD06
      Parameter (LUN_MOD06 = 412500)

      integer    FAIL,    SUCCEED
      parameter (FAIL=-1, SUCCEED=0)


C input argument declarations
      character*(*) MET_Handles(*)
      integer No_InputPntr,LUN_InputPntr(*),VRSN_InputPntr(*)


C local variable declarations
      character*25     msg25_2,msg25_3,msg25_4,msg25_LUN,msg25_VRSN
      character*255    HDF_AttrName, ECS_AttrName
      character*(1024) msgbuf
      character*(PGSd_PC_UREF_LENGTH_MAX) New_UR(MAX_NO_PNTRS),
     1                                    UR(MAX_NO_PNTRS) /MAX_NO_PNTRS*PGSd_MET_STR_END/ 

      integer rtn_flag, String_Loc, PGS_MET_GetPCAttr_s, PGS_MET_SetAttr_s
      integer i1, i2, rtn_loc, rtn, 
     1        fbyte,fbyte1,fbyte2,fbyte3,fbyte4,fbyte_LUN,fbyte_VRSN,
     2        lbyte,lbyte1,lbyte2,lbyte3,lbyte4,lbyte_LUN,lbyte_VRSN
      integer No_Duplicates,No_UR,No_New_UR,No_Pntr,No_Total_UR,No_Unique,
     1        UR_Array_Index,VRSN_MOD06

      logical error_flag, In_Existing_UR_Array 


c-----Initialization
      Update_InputPntr_MOD06 = SUCCEED
      error_flag = .FALSE.
      No_Unique = 0 
      No_Duplicates = 0

c-----Read existing InputPointer from MOD06 product file 
      VRSN_MOD06 = 1
      HDF_AttrName = MECS_CORE
      ECS_AttrName = MCORE_INPUT_POINTER
      rtn_loc = String_Loc(ECS_Attrname,fbyte1,lbyte1)

      write(msg25_LUN,'(I25)') LUN_MOD06
      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

      write(msg25_VRSN,'(I25)') VRSN_MOD06
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)


      rtn = PGS_MET_GetPCAttr_s(LUN_MOD06, VRSN_MOD06, HDF_AttrName, ECS_AttrName, UR)  

c-----If read successful, append new URs to existing UR array
      If (rtn .EQ. PGS_S_SUCCESS) Then
         No_UR = 1

c--------search existing UR array for PGSd_MET_STR_END 
         Do WHILE ( (UR(No_UR) .NE. PGSd_MET_STR_END) .AND.
     1              (No_UR .LE. MAX_NO_PNTRS) ) 

            No_UR = No_UR + 1
         End Do

         No_UR = No_UR - 1
         No_Total_UR = No_UR

c--------retrieve UR set for current process 
         No_Pntr = No_InputPntr
         If (No_InputPntr .GT. MAX_NO_PNTRS) No_Pntr = MAX_NO_PNTRS  

         Call Retrieve_UR_Set(No_Pntr,LUN_InputPntr, VRSN_InputPntr,
     1                        No_New_UR,New_UR,rtn_flag)

         If (rtn_flag .EQ. SUCCEED) Then
 
            Do 100 i1 = 1, No_New_UR
               In_Existing_UR_Array = .FALSE.

c--------------search UR array set for duplicate name  
               Do 200 i2 = 1, No_UR

                  If (New_UR(i1) .EQ. UR(i2)) Then
                     In_Existing_UR_Array = .TRUE.
                     No_Duplicates = No_Duplicates + 1
                  End If

 200           Continue


c--------------check results of duplicate name search 
               If ( .NOT. In_Existing_UR_Array ) Then
                  No_Total_UR = No_Total_UR + 1
                  No_Unique   = No_Unique + 1

                  If (No_Total_UR .LE. MAX_NO_PNTRS) Then 
                     UR_Array_Index = No_Total_UR
                     UR(UR_Array_Index) = New_UR(i1)
                  Else
                     write(msg25_2,'(i25)') No_Total_UR
                     rtn_loc = String_Loc(msg25_2,fbyte2,lbyte2) 

                     msgbuf = 'Total number of URs ' // msg25_2(fbyte2:lbyte2) // ' exceeds ' 
     1               // FUNCNAME // ' internal storage buffer size.'  
     2               // char(10) // 'New UR ' //  New_UR(i1) // ' not aggregated to ECS '
     3               // 'INPUTPOINTER array.'

                     Call MODIS_SMF_SetDynamicMsg(MODIS_W_GENERIC,msgbuf,FUNCNAME)
                  End If  
               End If


100         Continue 

c-----------report the number of new, unique duplicate file references
            write(msg25_3,'(i25)') No_Unique
            rtn_loc = String_Loc(msg25_3,fbyte3,lbyte3)

            write(msg25_4,'(i25)') No_Duplicates
            rtn_loc = String_Loc(msg25_4,fbyte4,lbyte4)

            msgbuf = msg25_3(fbyte3:lbyte3) // ' unique file names added to '
     1      // 'INPUTPOINTER array.  '
     2      // char(10) // msg25_4(fbyte4:lbyte4) // ' duplicate file names found.'
 
            Call MODIS_SMF_SetDynamicMsg(MODIS_A_GENERIC,msgbuf,FUNCNAME) 

         Else
            error_flag = .TRUE.

            msgbuf = 'Retrieve_UR_Set detected error acquiring input file URs. '
     1      // char(10) // 'Only MOD06 URs already in product are returned.  '
     2      // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3      // char(10) // 'messages originating within function Retrieve_UR_Set.'

 
            Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If
            

C--------If sufficient storage in UR buffer, add end-of-string marker to UR array 
         If ( No_Total_UR .lt. MAX_NO_PNTRS) Then
            UR_Array_Index = No_Total_UR + 1
            UR(UR_Array_Index) = PGSd_MET_STR_END
         End If


         ECS_AttrName = MCORE_INPUT_POINTER

         rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),ECS_AttrName,UR)

         If (rtn .NE. PGS_S_SUCCESS) Then
             error_flag = .TRUE. 
             rtn_loc = String_Loc(ECS_AttrName,fbyte,lbyte)

             msgbuf = 'PGS_MET_SetAttr_s detected error while setting attribute ' 
     1        // ECS_AttrName(fbyte:lbyte)
     2        // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3        // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4        // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
         End If

      Else  
         error_flag = .TRUE.

         msgbuf = 'PGS_MET_GetPCAttr_s detected error while reading ECS attribute '
     1   // ECS_Attrname(fbyte1:lbyte1)
     2   // char(10) // 'from MODIS cloud product on LUN = ' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     3   // ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     4   // char(10) // ECS_Attrname(fbyte1:lbyte1)//' not set to HDF product file. '
     5   // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     6   // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7   // char(10) // 'fault is identified, stage correct PCF/input file and '
     8   // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
      End If   !  INPUTPOINTER read check

      If (error_flag) Update_InputPntr_MOD06 = FAIL 

      Return

      End

C-----------------------------------------------------------------------

      Subroutine Retrieve_UR_Set(No_InputPntr,LUN_InputPntr,
     1                           VRSN_InputPntr,No_UR,UR,rtn_flag)

      implicit none
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Retrieve ECS Universal References (UR) for a specified
C                file set (identified by PCF LUN and version numbers) 
C
C !INPUT PARAMETERS:
C  integer  No_InputPntr    Number of valid file references in LUN_InputPntr
C                           and VRSN_InputPntr array arguments.  Ancillary 
C                           data files, look up tables, and MODIS product 
C                           files are included.  System files, such as the 
C                           MCF, are not included.
C
C  integer  LUN_InputPntr(*)
C                           Array containing a list of PCF file LUNs (See 
C                           No_InputPntr above).
C
C                           A one-to-one correspondence between the elements
C                           of the arrays LUN_InputPntr and VRSN_InputPntr
C                           is assumed.
C
C  integer  VRSN_InputPntr(*)
C                           Array containing a list of PCF file version 
C                           numbers  (See No_InputPntr and LUN_InputPntr 
C                           above).
C
C !OUTPUT PARAMETERS:
C  character UR            Array containing non-blank URs retrieved from PCF 
C
C  integer No_UR           Variable containing number of non-blank URs retrieved 
C                          from PCF 
C
C  integer rtn_flag        Procedure return status: SUCCEED (0) or FAIL (-1) 
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Retrieve_UR_Set')

      integer MAX_NO_PNTRS
      parameter (MAX_NO_PNTRS=50)

      integer FAIL, SUCCEED
      parameter (FAIL=-1, SUCCEED=0)

C passed-in argument declarations
      integer No_InputPntr,LUN_InputPntr(*), VRSN_InputPntr(*),No_UR
      character*(PGSd_PC_UREF_LENGTH_MAX) UR(MAX_NO_PNTRS)
      integer rtn_flag

C local variable declarations
      integer Version_No,LUN
      integer icounter,rtn
      integer fbyteb,fbytec,fbyted,lbyteb,lbytec,lbyted

      character*25 msg25b,msg25c,msg25d
      character*512 msgbuf
      character*(PGSd_PC_UREF_LENGTH_MAX) work_buf

      integer PGS_PC_GetUniversalRef
      integer String_Loc

c-----initialization
      rtn_flag = SUCCEED
      No_UR = 0


c-----loop over number of file references
      Do 100 icounter = 1, No_InputPntr
         Version_No = VRSN_InputPntr(icounter)
         LUN = LUN_InputPntr(icounter)

c-----set message buffers for error reporting 
         write(msg25b,'(i25)') LUN 
         write(msg25c,'(i25)') Version_No 
         write(msg25d,'(i25)') icounter

         rtn = String_Loc(msg25b,fbyteb,lbyteb)
         rtn = String_Loc(msg25c,fbytec,lbytec)
         rtn = String_Loc(msg25d,fbyted,lbyted)

c--------set error flag if file version number is out of bounds
         If (Version_No .LT. 1) Then
            rtn_flag = FAIL

            msgbuf = 'PCF file version number ' // msg25c(fbytec:lbytec) // ' on LUN '
     2       // msg25b(fbyteb:lbyteb) // ' is out of bounds.'
     3       // char(10) // 'UR of element ' // msg25d(fbyted:lbyted) 
     4       // ' of input pointer array cannot be retrieved.' 
     5       // char(10) // 'Operator Action:  Notify SDST.'

       
            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------set error flag if file LUN is out of bounds
         Else If (LUN.LT.1 .OR. (LUN.GE.10000 .AND. LUN.LE.10999)) Then
            rtn_flag = FAIL

            msgbuf = 'PCF file LUN number ' // msg25b(fbyteb:lbyteb)
     1       // ' is out of bounds.  ' 
     2       // char(10) // 'UR of element ' // msg25d(fbyted:lbyted) 
     3       // ' of input pointer array cannot be retrieved.' 
     4       // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------retrieve PCF UR
         Else 
            rtn = PGS_PC_GetUniversalRef(LUN, Version_No, work_buf)

            IF (rtn .ne. PGS_S_SUCCESS) Then
               rtn_flag = FAIL

               msgbuf = 'PGS_PC_GetUniversalRef unable to retrieve '
     1          // 'Universal Reference Identifier '//char(10)
     2          // 'on PCF LUN = '// msg25b(fbyteb:lbyteb) // ';  Version Number = '
     3          // msg25c(fbytec:lbytec) // char(10)
     4          // char(10) // 'Operator Action:  Check PCF Universal Reference Identifier. '
     5          // char(10) // 'If invalid, stage correct PCF and rerun PGE.  Otherwise, notify SDST.'


               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

c-----------if UR successfully read, insert it into output array
            Else
               No_UR = No_UR + 1
               UR(No_UR) = work_buf
            End If !(rtn .ne. PGS_S_SUCCESS)

         End If !(Version_No .LT. 1)

100   Continue

      Return
      End

C-----------------------------------------------------------------------

      Integer Function Update_PSA_MOD06(MET_Handles,NUM_PSA,Name_PSA,Value_PSA)
 
      Implicit None

C-----Insert Include files
      Include 'mapi.inc'
      Include 'PGS_MET_13.f'
      Include 'PGS_SMF.f' 
      Include 'PGS_MODIS_39500.f' 

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Increment and set the number ECS PSA attributes to be 
C                stored as MOD06 product metadata. 
C
C !INPUT PARAMETERS:
C
C  integer NUM_PSA          Variable containing the number of
C                           additional Product-Specific Attributes (PSAs) 
C                           to be set by process.
C
C  character*(*) Name_PSA(*)
C                           Array containing the additional PSA names 
C                           to be set by process 
C
C                           A one-to-one correspondence between the
C                           elements of arrays Name_PSA and Value_PSA
C                           is assumed.
C
C  Real Value_PSA(*)        Array containing the additional PSA values 
C                           to be set by process (See Name_PSA
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
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        string_loc                (science code)
C
C
C    Named Constant:
C        MECS_CORE               (mapi.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        MODIS_W_GENERIC         (PGS_MODIS_39500.f)
C        PGSMET_E_NULL_PARAMETER (PGS_MET_13.f
C        PGS_S_SUCCESS           (PGS_SMF.f)
C
C
C !END
C----------------------------------------------------------------------


C-----Parameter declarations
      Character*(*)  FUNCNAME
      Parameter     (FUNCNAME = 'Update_PSA_MOD06')

      Character*(*)  BLANK
      Parameter     (BLANK = ' ')

      integer    FAIL,    SUCCEED
      parameter (FAIL=-1, SUCCEED=0)

      integer    LUN_MOD06
      parameter (LUN_MOD06 = 412500)

      integer    NUM_PSA_MOD06
      parameter (NUM_PSA_MOD06 = 20)

      real        MAX_ABS_VALUE_PSA
      parameter ( MAX_ABS_VALUE_PSA = 9999.99 )

      real        REL_EQUALITY_EPS,          NO_SET_FLAG
      parameter ( REL_EQUALITY_EPS=0.000001, NO_SET_FLAG = -99999.0 )


C-----Product argument declarations
      Character*(*) Name_PSA(*), MET_Handles(*)
      Integer       NUM_PSA
      Real          Value_PSA(*)


C-----Local variable declarations
      Character*8    msg8
      Character*25   msg25_1,msg25_2,msg25_4,msg25_LUN,msg25_VRSN
      Character*30   PSA_Name
      Character*60   PSA_Value
      Character*30   Name_PSA_MOD06(NUM_PSA_MOD06)
     1               / 'SuccessCloudTopPropRtrPct_IR',     'SuccessCloudPhaseRtrPct_IR',
     2                 'SuccessCloudOpticalPropRtr_VIS',   'LowCloudDetectedPct_IR',
     3                 'MidCloudDetectedPct_IR',           'HighCloudDetectedPct_IR',
     4                 'ThinCloudDetectedPct_IR',          'ThickCloudDetectedPct_IR',
     5                 'OpaqueCloudDetectedPct_IR',        'CirrusCloudDetectedPct_IR',
     6                 'IceCloudDetectedPct_IR',           'WaterCloudDetectedPct_IR',
     7                 'MixedCloudDetectedPct_IR',         'CloudPhaseUncertainPct_IR',
     8                 'OceanCoverFractionPct',            'LandCoverFractionPct',
     9                 'SnowCoverFractionPct',             'CloudCoverFractionPct_VIS',
     1                 'WaterCloudDetectedPct_VIS',        'IceCloudDetectedPct_VIS'       / 
      character*100  HDF_FileAttrName
      Character*1024 msgbuf

      Integer fbyte1,fbyte2,fbyte4,fbyte_LUN,fbyte_VRSN,
     1        lbyte1,lbyte2,lbyte4,lbyte_LUN,lbyte_VRSN,
     2        rtn_loc,rtn1,rtn2
      Integer i0,i1,i2,NUM_Blank_PSA,PSA_CLass,VRSN_MOD06
      Integer pgs_met_getpcattr_s,Set_PSA,string_loc

      Logical error_flag, On_PSA_Update_List

c-----Initialization
      error_flag = .FALSE.
      Update_PSA_MOD06 = SUCCEED 
       

c-----------------------------------------------------------------------
c Perform input argument checks
c-----------------------------------------------------------------------

      VRSN_MOD06 = 1
c-----Set up variables for status messaging
      write(msg25_LUN,'(I25)') LUN_MOD06
      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

      write(msg25_VRSN,'(I25)') VRSN_MOD06
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)


C --- check for positive number of PSAs 
      If (NUM_PSA .LE. 0) Then
         write(msg25_1,'(i25)') NUM_PSA
         rtn_loc = String_Loc(msg25_1,fbyte1,lbyte1)

         msgbuf = 'The number of additional PSAs is ' // msg25_1(fbyte1:lbyte1) 
     1   // char(10) // 'No new PSAs to be set by ' // FUNCNAME 
     2   // ', but existing ones will be copied.  ' 

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)
      End If
 
c-----check for blank PSAs
      NUM_Blank_PSA = 0

      Do 50 i0 = 1,NUM_PSA
         If (Name_PSA(i0) .eq. BLANK) NUM_Blank_PSA = NUM_Blank_PSA + 1 
 50   Continue

      If (NUM_Blank_PSA .GT. 0) Then
         write(msg25_1,'(i25)') NUM_Blank_PSA
         rtn_loc = String_Loc(msg25_1,fbyte1,lbyte1)

         msgbuf = msg25_1(fbyte1:lbyte1) // ' PSAs on input argument '
     1            // 'list have blank names.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME) 
      End If
   

c-----------------------------------------------------------------------
c Loop over PSAs in MODIS Cloud Product.  Copy existing PSAs and set 
c any new ones to internal memory.
c-----------------------------------------------------------------------

      Do 100 i1 = 1, NUM_PSA_MOD06
         PSA_Name     = Name_PSA_MOD06(i1)
         PSA_Class    = i1
         HDF_FileAttrName = MECS_CORE
         
         rtn1 = pgs_met_getpcattr_s(LUN_MOD06,VRSN_MOD06,HDF_FileAttrName,PSA_Name,PSA_Value)

         rtn_loc = string_loc(PSA_Name,fbyte1,lbyte1)

c--------PSA read was successful. Copy it to memory
         If (rtn1 .EQ. PGS_S_SUCCESS) Then

            rtn2 = Set_PSA(MET_Handles,PSA_NAME,PSA_Class,PSA_Value)

            If (rtn2 .EQ. FAIL) Then
               error_flag = .true.
               msgbuf = 
     2         'Set_PSA detected error setting MOD06 PSA to ECS internal memory. '
     3         // char(10) // 'PSA ' // PSA_Name(fbyte1:lbyte1) // ' not set.' 
     4         // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     5         // char(10) // 'messages originating within routine Set_PSA. '

 
               Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            End If


c--------PSA is not in product file.  
         Else If (rtn1 .EQ. PGSMET_E_NULL_PARAMETER) Then
 
c--------code for possible future use
c        Else If ( (rtn1 .EQ. PGSMET_E_NULL_PARAMETER) .OR. 
c    1             (rtn1 .EQ. PGSMET_E_FILETOODL_ERR) ) Then

            If (rtn1 .EQ. PGSMET_E_NULL_PARAMETER) Then
               msgbuf = 'Oops, PSA ' // PSA_Name(fbyte1:lbyte1) 
     1         // ' is not in product file! ' // char(10) 
     2         // 'Set it now if it is on update PSA list.' 

               Call MODIS_SMF_SetDynamicMsg(MODIS_W_GENERIC,msgbuf,FUNCNAME)
            End If


c-----------code for possible future use
c           If (rtn1 .EQ. PGSMET_E_FILETOODL_ERR) Then 
c              msgbuf = ' HDF attribute ' // MECS_CORE 
c    1         // ' is not in product file. ' // char(10) 
c    2         // 'Set only PSAs on update list.' 
c           End If


            On_PSA_Update_List = .FALSE.

c-----------Loop on PSA input argument list 
            Do 200 i2 = 1, NUM_PSA
               write(msg25_2,'(i25)') i2 
               rtn_loc = String_Loc(msg25_2,fbyte2,lbyte2)

c--------------PSA on update list 
               If ( (Name_PSA(i2) .EQ. PSA_Name) .AND. 
     1              (abs((Value_PSA(i2)-NO_SET_FLAG)/NO_SET_FLAG) .gt. REL_EQUALITY_EPS) ) Then

                 On_PSA_Update_List = .TRUE.

c-----------------PSA value out of range check; cannot set Name/Value pair
                  If ( Abs(Value_PSA(i2)) .gt. MAX_ABS_VALUE_PSA) Then
                     error_flag = .TRUE.
        
                     write(msg25_4,'(E25.6)') Value_PSA(i2)
                     rtn_loc = String_Loc(msg25_4,fbyte4,lbyte4)

                     msgbuf = 
     1                'PSA value ' // msg25_4(fbyte4:lbyte4) // ' is out of bounds! ' 
     2                // char(10) // 'PSA ' // PSA_Name(fbyte1:lbyte1) 
     3                // ' Name/Value pair not updated.'
     4                // char(10) // 'Operator Action:  Notify SDST. '


                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               
c-----------------Set PSA Name/Value pair 
                  Else
                     write(msg8,'(f8.2)') Value_PSA(i2) 

                     rtn2 = Set_PSA(MET_Handles,PSA_NAME,PSA_Class,msg8)

                     If (rtn2 .EQ. FAIL) Then
                        error_flag = .TRUE.
                        msgbuf = 'Set_PSA detected error setting PSA ' 
     1                  // PSA_Name(fbyte1:lbyte1) // ' Name/Value pair. ' 
     2                  // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3                  // char(10) // 'messages originating within routine Set_PSA. '


                        Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                     Else
                        msgbuf = 'PSA ' //  PSA_Name(fbyte1:lbyte1) 
     1                   // ' Name/Value pair successfully updated.'

                        Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
                     End If
    

                  End If   ! PSA valid value check

               End If   ! PSA on update list check 

 200        Continue   ! Loop on Input PSAs
 
            If (.NOT. On_PSA_Update_List) Then
               msgbuf = 'PSA ' // PSA_Name(fbyte1:lbyte1) 
     1         // ' not set; not on update PSA list.' 

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_A_GENERIC,msgbuf,FUNCNAME)
            End If

         Else 
            error_flag = .TRUE.

            msgbuf = 'PGS_MET_GetPCAttr_d unable to retrieve '// PSA_Name(fbyte1:lbyte1) 
     1       // char(10) // 'from MOD06 product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2       // ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3       // char(10) // PSA_Name(fbyte1:lbyte1) // ' not set.'
     4       // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5       // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6       // char(10) // 'fault is identified, stage correct PCF/input file and '
     7       // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End IF   ! End PSA existence check 

 100  Continue   ! End loop on MOD06 PSAs


      If (error_flag) Update_PSA_MOD06 = FAIL

      Return
      End

C-----------------------------------------------------------------------

      integer function Set_PSA(MET_Handles,PSA_Name,PSA_Class,PSA_Value)
      Implicit None

C-----Insert Include files
      Include 'PGS_SMF.f' 
      Include 'PGS_MODIS_39500.f' 
      Include 'mapi.inc' 

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Set an ECS AdditionalAttribute name/value metadata 
C                pair.
C
C
C !INPUT PARAMETERS:
C
C  character*(*) PSA_Name   Variable containing the name of PSA
C                           to be set for granule
C
C  character*(*) PSA_Value  Variable containing the value of PSA
C                           to be set for granule.
C
C  integer       PSA_Class  Variable containing the PSA class.
C
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
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by         Liqun Ma                   Feb. 1998
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    lma@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C
C
C    Named Constant:
C        MODIS_E_GENERIC       (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS         (PGS_SMF.f)
C
C  Internals:
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG
C        string_loc
C
C
C !END
C----------------------------------------------------------------------


C-----Parameter declarations
      Character*(*)  FUNCNAME
      Parameter     (FUNCNAME = 'Set_PSA')
     
      Integer    INVENTORYMETADATA
      Parameter (INVENTORYMETADATA = 2)

      integer    FAIL,    SUCCEED
      Parameter (FAIL=-1, SUCCEED=0)

C-----function argument declarations
      Character*(*) MET_Handles(*),PSA_Name,PSA_Value
      Integer       PSA_Class

C-----Local variable declarations
      character*25 msg25
      character*190 msgbuf
      character*40 AttrN

      Integer rtn,pgs_met_setattr_s,String_Loc
      Integer fbyte,lbyte
      Integer fbyteb,lbyteb

      logical error_flag


c-----Initialization
      Set_PSA=FAIL
      error_flag = .false.
      write(msg25,'(i25)') PSA_Class
      rtn = String_Loc(msg25,fbyte,lbyte)

c-----Set PSA Name
      AttrN = MCORE_ADD_ATTRIBUTENAME// '.' // msg25(fbyte:lbyte) 

      rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,PSA_Name)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true. 
         rtn = String_Loc( AttrN,fbyteb,lbyteb)

         msgbuf = 'PGS_MET_SetAttr_s unable to set ' // AttrN(fbyteb:lbyteb)
     1    // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2    // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3    // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
      End If


c-----Set PSA value 
      AttrN = MCORE_PARAMETERVALUE // '.' // msg25(fbyte:lbyte) 

      rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,PSA_Value)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true. 
         rtn = String_Loc( AttrN,fbyteb,lbyteb)

         msgbuf = 'PGS_MET_SetAttr_s unable to set '// AttrN(fbyteb:lbyteb)
     1    // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2    // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3    // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
      End If       
 
      If(.not.error_flag) Set_PSA=SUCCEED

      Return
      END

C-----------------------------------------------------------------------

      Integer Function Update_MeasParm_MOD06(MET_Handles,No_MeasParm)
      implicit none
      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Read and reset MOD06 ECS Measured Parameter inventory 
C                metadata group as previously written to product file by 
C                prior process. The MOD06 ECS Measured Parameter group 
C                consists of the following ODL objects:
C
C                "PARAMETERNAME"
C                "QAPERCENTMISSINGDATA" 
C                "AUTOMATICQUALITYFLAG" 
C                "AUTOMATICQUALITYFLAGEXPLANATION"
C
C
C !INPUT PARAMETERS:
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
C  integer No_MeasParm      Variable containing the number of Measured
C                           Parameters contained in MOD06. At launch, this
C                           is 1. 
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
C !END
C----------------------------------------------------------------------
C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Update_MeasParm_MOD06')

      integer LUN_MOD06
      Parameter (LUN_MOD06 = 412500)

      integer INVENTORYMETADATA
      Parameter (INVENTORYMETADATA = 2)

      integer    FAIL,    SUCCEED
      Parameter (FAIL=-1, SUCCEED=0)


C function arguments declarations
      integer    No_MeasParm
      character*(*) MET_Handles(*)

C other variable declarations
      character*25  msg25d,msg25_LUN,msg25_VRSN
      character*128 AttrN,AttrV_s,Field_Name(4),HDF_AttrName
      character*512 msgbuf

      integer    i,icounter,rtn,rtn_loc,VRSN_MOD06,AttrV_i
      integer    pgs_met_setattr_s, pgs_met_setattr_i,String_Loc,
     2           PGS_MET_GetPCAttr_i,PGS_MET_GetPCAttr_s
      integer    fbyte,fbyteb,lbyte,lbyteb
      integer    fbyte_LUN,fbyte_VRSN,lbyte_LUN,lbyte_VRSN

      logical error_flag


C------------------------
C Initialization
C------------------------
      Update_MeasParm_MOD06 = FAIL
      error_flag = .false.
      VRSN_MOD06 = 1
      HDF_AttrName = MECS_CORE

c-----Set up variables for status messaging
      write(msg25_LUN,'(I25)') LUN_MOD06
      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

      write(msg25_VRSN,'(I25)') VRSN_MOD06
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)



      Field_Name(1) = MCORE_PARAMETERNAME
      Field_Name(2) = MCORE_PERCENT_MISSING
      Field_Name(3) = MCORE_AUTO_QUALITY
      Field_Name(4) = MCORE_AQUAL_FLG


      Do icounter = 1, No_MeasParm
         write(msg25d,'(i25)') icounter
         rtn = String_Loc(msg25d,fbyte,lbyte)

         do i=1,4

            rtn = String_Loc(Field_Name(i),fbyteb,lbyteb)
            AttrN = Field_Name(i)(fbyteb:lbyteb) // '.' // msg25d(fbyte:lbyte) 

            if(i.eq.2) then
               rtn = PGS_MET_GetPCAttr_i(LUN_MOD06, VRSN_MOD06, HDF_AttrName, AttrN, AttrV_i)

               If (rtn .ne. PGS_S_SUCCESS) Then
                  error_flag = .true.
                  rtn = String_Loc(AttrN,fbyteb,lbyteb)

                  msgbuf = 'PGS_MET_GetPCAttr_i unable to retrieve '// AttrN(fbyteb:lbyteb)
     1            // char(10) // 'from MOD06 product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2            //  ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3            // char(10) // AttrN(fbyteb:lbyteb) // ' not set.'
     4            // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5            // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6            // char(10) // 'fault is identified, stage correct PCF/input file and '
     7            // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               Else
 
                  rtn = PGS_MET_SetAttr_i(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_i)

                  If (rtn .ne. PGS_S_SUCCESS) Then
                     error_flag = .true.
                     rtn = String_Loc(AttrN,fbyteb,lbyteb)

                     msgbuf = 'PGS_MET_SetAttr_i unable to set field '// AttrN(fbyteb:lbyteb)
     1                // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2                // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3                // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                     call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
                  End If

               End If
            else
               rtn = PGS_MET_GetPCAttr_s(LUN_MOD06, VRSN_MOD06, HDF_AttrName, AttrN, AttrV_s)

               If (rtn .ne. PGS_S_SUCCESS) Then
                  error_flag = .true.
                  rtn = String_Loc(AttrN,fbyteb,lbyteb)

                  msgbuf = 'PGS_MET_GetPCAttr_s unable to retrieve '//AttrN(fbyteb:lbyteb)
     1            // char(10) // 'from MOD06 product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2            //  ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3            // char(10) // AttrN(fbyteb:lbyteb) // ' not set.'
     4            // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5            // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6            // char(10) // 'fault is identified, stage correct PCF/input file and '
     7            // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               Else

                  rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_s)

                  If (rtn .ne. PGS_S_SUCCESS) Then
                     error_flag = .true.
                     rtn = String_Loc(AttrN,fbyteb,lbyteb)

                     msgbuf = 'PGS_MET_SetAttr_s unable to set field '// AttrN(fbyteb:lbyteb)
     1                // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2                // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3                // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                     call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

               End If

         end if

         END DO
      END DO

      If (.not.error_flag) Update_MeasParm_MOD06 = SUCCEED

      Return
      END
