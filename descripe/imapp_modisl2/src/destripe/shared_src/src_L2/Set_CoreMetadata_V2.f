      Integer Function Set_CoreMetadata(LRN_MCF,ExtGeoPntr_Flag,
     1                 No_PCF_Attr,LRN_PCF_Attr,Name_PCF_Attr,
     2                 No_InputPntr,LRN_InputPntr,Vrsn_InputPntr,
     3                 No_MeasParm,Name_MeasParm,PctMissing_MeasParm,
     4                 AutoFlag_MeasParm,AutoFlagExp_MeasParm,
     5                 No_PSA,Name_PSA,Value_PSA,
     6                 MET_Handles)

      implicit none
      include 'hdf.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C   Sets the value of EOSDIS Core System (ECS) inventory metadata 
C   fields that comprise the 'CoreMetadata.0' ECS-HDF attribute.  In 
C   this context, set means to associate a value with an ECS attribute 
C   name via function calls to the ECS Science Data Processing 
C   Toolkit (SDPTK V5.2).
C
C !DESCRIPTION:
C
C   Function Set_CoreMetadata sets the metadata field values within 
C   'CoreMetadata.0' which satisfy ECS's requirement for mandatory 
C   inventory metadata on HDF product files.  No metadata fields are 
C   actually written to the product file during the call to 
C   Set_CoreMetadata;  only the association of values to fields is made 
C   at this time.  Actual insertion of attribute 'CoreMetadata.0' into 
C   the HDF product file takes place later with a call to the SDPTK 
C   routine PGS_MET_Write.  This call is made from outside of 
C   Set_CoreMetadata, and only after any and all other core metadata 
C   fields have been set by the PGE.  If the MODIS HDF Application 
C   Programming Interface (M-API) Utility is being used, PGS_MET_Write 
C   is called automatically when the HDF product file is closed with 
C   CPMFIL.  If one does not use M-API, then a direct call to 
C   PGS_MET_Write is required to attach 'CoreMetadata.0' to the file 
C   before it is closed.
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
C                           alphanumeric characters (e.g. 'a01' or
C                           '001').    
C
C
C  integer LRN_PCF_Attr(*)  Array containing the PCF logical reference
C                           numbers (LRN) of ECS inventory field values 
C                           listed as USER DEFINED RUNTIME PARAMETERS. 
C
C                           A one-to-one correspondence between the 
C                           elements of the arrays LRN_PCF_Attr and
C                           Name_PCF_Attr is assumed.
C
C  character*(*) Name_PCF_Attr(*)
C                           Array containing the names of ECS inventory 
C                           metadata fields whose values are listed as 
C                           USER DEFINED RUNTIME PARAMETERS in the PCF.
C                           Note:  these names are not retrieved from 
C                           the PCF as are their associated values; they 
C                           must be assigned by the science code. (See
C                           LRN_PCF_Attr above) 
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
C  integer No_MeasParm      Variable containing the number of Measured 
C                           Parameters to be set for granule.  If 
C                           none, No_MeasParm=0.
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
C                           passed to Set_CoreMetadata as integers. 
C                            
C  character*(*) AutoFlag_MeasParm(*)  
C                           Array containing a "simple" assessment of 
C                           the quality of each Measured Parameter 
C                           value to be set for granule.  
C
C                           Recommended values for the assessment are 
C                           "passed" or "failed".  If the value of 
C                           AutoFlag_MeasParm is blank (' '), neither 
C                           AutoFlag_MeasParm nor AutoFlagExp_MeasParm
C                           will be set by function Set_CoreMetadata
C                           (See Name_MeasParm above).  
C
C  character*(*) AutoFlagExp_MeasParm(*)
C                           Array containing the "criterion" for the
C                           quality judgement made on each Measured 
C                           Parameter value.  If the value of the 
C                           AutoFlag_MeasParm is set to blank (' '), 
C                           the AutoFlagExp_MeasParm will not be set
C                           (see AutoFlag_MeasParm and Name_MeasParm
C                           above).
C
C  integer No_PSA           Variable containing the number of 
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
C Revision 1.10  1998/07/09  14:34:17  fhliang
C changed 'MECS_CORE' to 'MECS_ARCHIVE' on L.441.
C
C Revision 1.9  1998/06/08  20:22:52  fhliang
C Updated error messages with "Operator Action" strings.
C
C Revision 1.8  1997/12/03  17:36:45  rhucek
C changed lun reference to geolocation from 6000000 to 600000.
C
C Revision 1.7  1997/11/26  18:10:48  rhucek
C Added AssociatedPlatformInstrumentSensor Group and removed
C SensorCharacteristic Group.
C
C Revision 1.6  1997/11/19  20:43:00  lma
C Integrated "Gen_Modis_FileName" to generate file name
C (localgranuleID), move out "v2fname" subroutine.
C
C Revision 1.4  1997/07/27  21:26:52  rhucek
C Numbered continuation lines 2,3,4,etc instead of *,*,*.  Updated
C prolog to distinguish parameters read in as RUNTIME PARAMETERS.
C
C Revision 1.3  1997/07/15  01:05:52  rhucek
C Updated error messages
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by         Vicky Lin                  June 1997
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    vlin@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C  Returns:     MAPIOK if successful, MFAIL if an error occurs
C
C  Externals:
C
C    Functions and Subroutines:
C        sfstart, sfend            (HDF library)
C        PGS_MET_Init              (libPGSTK.a)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        PGS_MET_GetSetAttr_*      (libPGSTK.a)
C        PGS_PC_GetConfigData      (libPGSTK.a)
C        PGS_PC_GetUniversalRef    (libPGSTK.a)
C
C
C    Named Constant:
C        MAPIOK                (mapi.inc)
C        MCORE_*               (mapi.inc)
C        MECS_CORE             (mapi.inc)
C        MAX_ECS_NAME_L        (mapi.inc)
C        MODFILLEN             (mapi.inc)
C        MFAIL                 (mapi.inc)
C        MODIS_E_GENERIC       (PGS_MODIS_39500.f)
C        MODIS_W_GENERIC       (PGS_MODIS_39500.f)
C        PGSd_PC_*             (PGS_PC.f)
C        PGSd_MET_*            (PGS_MET.f: included in mapi.inc)
C        PGS_S_SUCCESS         (PGS_SMF.f)
C        SUCCEED               (hdf.inc)
C
C  Internals:
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG
C        Gen_Modis_FileName
C        string_loc
C
C    Variables:
C        Field_Name   Core metadata field name
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'Set_CoreMetadata')

      character*(*)  BLANK
      parameter     (BLANK = ' ')

      integer        MAX_NO_PNTRS
      parameter     (MAX_NO_PNTRS=1000)

      integer        NO_SET_UR
      parameter     (NO_SET_UR = -99999)

      real           MAX_ABS_VALUE_PSA
      parameter     (MAX_ABS_VALUE_PSA = 9999.99 )

      real           Rel_equality_EPS,            NoSet_Flag
      parameter     (Rel_equality_EPS = 0.000001, NoSet_Flag = -99999.0)

      integer        LRN_GEO
      parameter     (LRN_GEO = 600000)

      integer        VRSN_GEO
      parameter     (VRSN_GEO = 1)

      integer        INVENTORYMETADATA
      parameter     (INVENTORYMETADATA = 2)

C other variable declarations
      character*8   ESDT_Name

      character*8   msg8
      character*25  msg25,msg25b,msg25c,msg25d
      character*128  fname
      character*128 AttrN,AttrV,Field_Name(8),buf_char
      character*512 msgbuf, runtimepar,parmname
      character*(*) MET_Handles(*)
      character*(PGSd_PC_UREF_LENGTH_MAX)  UR_Ref(MAX_NO_PNTRS)
      character*(PGSd_PC_VALUE_LENGTH_MAX) RunTimeParm
      character*(PGSd_PC_VALUE_LENGTH_MAX) LOCAL_VRSNID


      character*(*) AutoFlag_MeasParm(*),AutoFlagExp_MeasParm(*),
     2              Name_MeasParm(*),Name_PSA(*),Name_PCF_Attr(*)

      integer    i,No_MeasParm,No_PCF_Attr,No_InputPntr,No_PSA,
     2           icounter,rtn,Version_No,LRN_MCF

      integer    LRN,buf_int(4),
     2           LRN_PCF_Attr(*),LRN_InputPntr(*),Vrsn_InputPntr(*),
     3           PctMissing_MeasParm(*),LRN_VRN

      integer    pgs_met_init,pgs_met_setattr_d,
     2           pgs_met_setattr_s, pgs_met_setattr_i,
     3           PGS_MET_GetPCAttr_i,PGS_MET_GetPCAttr_d,
     4           PGS_MET_GetPCAttr_s,PGS_PC_GetConfigData,
     6           PGS_PC_GetUniversalRef,String_Loc,
     7           Gen_Modis_FileName,pgs_met_getsetattr_s

      integer    NUM_nonBlank_UR,NUM_Pointers
      integer    fbyte,fbyteb,fbytec,fbyted,
     7           lbyte,lbyteb,lbytec,lbyted

      real   Value_PSA(*)

      double precision buf_dbl(4)

      logical error_flag,ExtGeoPntr_Flag,loopexit

C------------------------
C Initialization
C------------------------

      Set_CoreMetadata = MFAIL
      error_flag = .false.
      fname = ' '

      Do 10 i = 1, 4
         buf_int(i) = 0
         buf_dbl(i) = 0.0D0
   10 continue

      write(msg25b,'(I25)') LRN_GEO
      rtn = string_loc(msg25b,fbyteb,lbyteb)

      write(msg25c,'(I25)') VRSN_GEO
      rtn = string_loc(msg25c,fbytec,lbytec)

C-------------------------------------------------------------------
C Initialize metadata tool defining array MET_Handles, and set 
C inventory attribute fields whose values are provided in the MCF.
C-------------------------------------------------------------------

      rtn = pgs_met_init(LRN_MCF,MET_Handles)

      If (rtn.ne.PGS_S_SUCCESS) Then
         msgbuf = 'pgs_met_init unable to initialize MCF.'
     2     // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3     // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4     // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

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

      Field_Name(1) = MCORE_EAST_BOUND
      Field_Name(2) = MCORE_NORTH_BOUND
      Field_Name(3) = MCORE_SOUTH_BOUND
      Field_Name(4) = MCORE_WEST_BOUND

      Do 20 i = 1, 4
         AttrN = Field_Name(i)

         rtn = PGS_MET_GetPCAttr_d(LRN_GEO,VRSN_GEO,MECS_ARCHIVE, AttrN,buf_dbl)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.
            rtn = String_Loc(AttrN,fbyte,lbyte)

            msgbuf = 'PGS_MET_GetPCAttr_d unable to retrieve '// AttrN(fbyte:lbyte) 
     2      // char(10) // 'from Geolocation product (LUN = ' // msg25b(fbyteb:lbyteb)
     3      // ', Version = ' // msg25c(fbytec:lbytec) // ').'
     3      // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     3      // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     4      // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     5      // char(10) // 'fault is identified, stage correct PCF/input file and '
     6      // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

            call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

         Else
            rtn = PGS_MET_SetAttr_d(MET_Handles(INVENTORYMETADATA),AttrN,buf_dbl)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)
               msgbuf = 'PGS_MET_SetAttr_d unable to set ' //AttrN(fbyte:lbyte) 
     2          // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3          // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4          // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

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

         rtn = PGS_MET_GetPCAttr_s(LRN_GEO,VRSN_GEO,MECS_CORE,AttrN,buf_char)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.
            rtn = String_Loc(AttrN,fbyte,lbyte)

            msgbuf = 'PGS_MET_GetPCAttr_s unable to retrieve '// AttrN(fbyte:lbyte) // '.'
     2      // char(10) // 'from Geolocation product (LUN = ' // msg25b(fbyteb:lbyteb)
     3      // ', Version = ' // msg25c(fbytec:lbytec) // ').'
     4      // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     5      // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     6      // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7      // char(10) // 'fault is identified, stage correct PCF/input file and '
     8      // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

            call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

         Else
            rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,buf_char)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)

               msgbuf = 'PGS_MET_SetAttr_s unable to set ' // AttrN(fbyte:lbyte) // '.'
     2          // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3          // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4          // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

            End If
         End If

   25 Continue


      AttrN = MCORE_ORBIT_NUM//'.1'

      rtn = PGS_MET_GetPCAttr_i(LRN_GEO,VRSN_GEO,MECS_CORE,AttrN,buf_int)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.
         rtn = String_Loc(AttrN,fbyte,lbyte)

         msgbuf = 'PGS_MET_GetPCAttr_i unable to retrieve '// AttrN(fbyte:lbyte) 
     2   // char(10) // 'from Geolocation product (LUN = ' // msg25b(fbyteb:lbyteb)
     3   // ', Version = ' // msg25c(fbytec:lbytec) // ').'
     4   // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     5   // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     6   // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7   // char(10) // 'fault is identified, stage correct PCF/input file and '
     8   // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

      Else
         rtn = PGS_MET_SetAttr_i(MET_Handles(INVENTORYMETADATA),AttrN,buf_int)

         If (rtn.ne.PGS_S_SUCCESS) Then
             error_flag = .true.
             rtn = String_Loc(AttrN,fbyte,lbyte)

             msgbuf = 'PGS_MET_SetAttr_i unable to set ' // AttrN(fbyte:lbyte) // '.'
     2       // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3       // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4       // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

             call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

         End If
      End If


      AttrN = MCORE_EQUATCROSSINGLONG//'.1'

      rtn = PGS_MET_GetPCAttr_d(LRN_GEO,VRSN_GEO,MECS_CORE,AttrN,buf_dbl)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.
         rtn = String_Loc(AttrN,fbyte,lbyte)

         msgbuf = 'PGS_MET_GetPCAttr_d unable to retrieve '// AttrN(fbyte:lbyte) 
     2   // char(10) // 'from Geolocation product (LUN = ' // msg25b(fbyteb:lbyteb)
     3   // ', Version = ' // msg25c(fbytec:lbytec) // ').'
     4   // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     5   // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     6   // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7   // char(10) // 'fault is identified, stage correct PCF/input file and '
     8   // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

      Else
         rtn = PGS_MET_SetAttr_d(MET_Handles(INVENTORYMETADATA),AttrN,buf_dbl)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.
            rtn = String_Loc(AttrN,fbyte,lbyte)

            msgbuf = 'PGS_MET_SetAttr_d unable to set ' // AttrN(fbyte:lbyte) // '.'
     2      // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3      // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4      // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

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

         If (AttrN .eq. ' ') Then
            error_flag = .true.

            msgbuf = 'PCF attribute name is blank.' // char(10) //
     2      'Set_CoreMetadata unable to set unknown attribute '
     3       //char(10)//'Operator Action: Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

         Else

            rtn = PGS_PC_GetConfigData(LRN_PCF_Attr(icounter),RunTimeParm)

            write(msg25b,'(I25)') LRN_PCF_Attr(icounter)
            rtn = string_loc(msg25b,fbyteb,lbyteb)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)

               msgbuf = 'PGS_PC_GetConfigData unable to retrieve ' // AttrN(fbyte:lbyte) // '.'
     2         // char(10) // 'runtime parameter value (LUN = ' // msg25b(fbyteb:lbyteb) // '). '
     3         // AttrN(fbyte:lbyte) // ' not set.'
     4         // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter value.'
     5         // char(10) // 'If incorrect, stage correct PCF and rerun PGE. '
     6         // char(10) // 'Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

            Else

C..............ECS tools allow setting a "blank" value
               If (RunTimeParm .eq. ' ') Then
                  error_flag = .true.
                  rtn = String_Loc(AttrN,fbyte,lbyte)

                  msgbuf = 'Run time parameter value (LUN = ' // msg25b(fbyteb:lbyteb) // ') is blank.'
     2            // char(10) // 'Attribute name is ' // AttrN(fbyte:lbyte) // '.'
     3            // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter '
     4            // char(10) // 'value, If empty, stage correct PCF and rerun PGE. '
     5            // char(10) // 'Otherwise, notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               End If
        
               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,RunTimeParm)

               If (rtn.ne.PGS_S_SUCCESS) Then
                  error_flag = .true.
                  rtn = String_Loc(AttrN,fbyte,lbyte)

                  msgbuf = 'PGS_MET_SetAttr_s unable to set ' // AttrN(fbyte:lbyte) 
     2            // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3            // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4            // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

               End If
            End If 
         End If

   30 Continue

C--------------------------------------------------------------------
C Set ECS inventory metadata attributes whose values are passed in 
C as actual function arguments from PGE.  These include the Measured 
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


C - Let's set Measured Parameter attributes.  There are four 
C - fields to set, some of which are optional (see above notes).

      Do 40 icounter = 1, No_MeasParm
         AttrV = Name_MeasParm(icounter)
         write(msg25d,'(i25)') icounter
         rtn = String_Loc(msg25d,fbyted,lbyted)

c - If Measured Parameter name is blank, report error
         If (AttrV .eq. ' ') Then
            error_flag = .true.

            msgbuf = 'Measured Parameter name is blank.' // char(10) 
     2      // 'Set_CoreMetadata unable to set unknown attribute '
     3      // char(10) // 'Measured Parameter number is ' // msg25d(fbyted:lbyted) // '.'
     4      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

         Else

c - Set Measured Parameter attribute PARAMETERNAME 
            AttrN = MCORE_PARAMETERNAME// '.' // msg25d(fbyted:lbyted) 

            rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)

               msgbuf =
     1         'pgs_met_setattr_s unable to set attribute ' // AttrN(fbyte:lbyte) // '.'
     2         // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3         // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4         // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

            End If


c ---       Set Measured Parameter QAPERCENTMISSINGDATA 
            AttrN = MCORE_PERCENT_MISSING // '.' // msg25d(fbyted:lbyted) 

            rtn = PGS_MET_SetAttr_i(MET_Handles(INVENTORYMETADATA),AttrN,PctMissing_MeasParm(icounter))

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)

               msgbuf =
     1         'pgs_met_setattr_i unable to set attribute ' // AttrN(fbyte:lbyte) // '.'
     2               // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3               // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4               // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

            End If


C ---       Set Measured Parameters AutomaticQualityFlag and
C ---       AutomaticQualityFlagExplanation if requested.

            If (AutoFlag_MeasParm(icounter).ne.' ') Then
               AttrN = MCORE_AUTO_QUALITY // '.' // msg25d(fbyted:lbyted)
  
               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AutoFlag_MeasParm(icounter))

               If (rtn.ne.PGS_S_SUCCESS) Then
                  error_flag = .true.
                  rtn = String_Loc(AttrN,fbyte,lbyte)

                  msgbuf =
     1            'pgs_met_setattr_s unable to set attribute ' // AttrN(fbyte:lbyte) // '.'
     2            // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3            // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4            // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
               End If

               AttrN = MCORE_AQUAL_FLG // '.' // msg25d(fbyted:lbyted)  

               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA), AttrN,AutoFlagExp_MeasParm(icounter))

               If (rtn.ne.PGS_S_SUCCESS) Then
                  error_flag = .true.
                  rtn = String_Loc(AttrN,fbyte,lbyte)

                  msgbuf =
     1            'pgs_met_setattr_s unable to set attribute ' // AttrN(fbyte:lbyte) // '.'
     2            // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3            // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4            // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'


                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

               End If
            End If
         End If

   40 Continue


C-----------------------------------------------------------------------
C  Set Product Specific Attributes.  There are two fields: 
C  "AdditionalAttributeName"  & "ParameterValue"
C-----------------------------------------------------------------------

C --- check for positive number of PSAs 
      If (No_PSA .LE. 0) Then
         write(msg25,'(i25)') No_PSA
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'The number of PSAs (' // msg25(fbyte:lbyte) // ') is less than 1.'
     2   // char(10) // 'No PSAs set by Set_CoreMetadata.' 


         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)
 
C --- number of PSAs is positive
      Else
 
         Do 50 icounter = 1, No_PSA
            write(msg25,'(i25)') icounter
            rtn = String_Loc(msg25,fbyte,lbyte)

            write(msg25b,'(E25.3)') Value_PSA(icounter)
            rtn = String_Loc(msg25b,fbyteb,lbyteb)

            If ( abs( (Value_PSA(icounter) - NoSet_Flag)/NoSet_Flag) .gt. Rel_equality_EPS) Then

C --- check for blank PSA name
                If (Name_PSA(icounter) .eq. ' ') Then
                   error_flag = .true.

                   msgbuf = 'PSA name is blank.' // char(10) // 'Set_CoreMetadata unable to set ' 
     2             // msg25(fbyte:lbyte) // 'th element of PSA array.' 
     3             // char(10) // 'Operator Action:  Notify SDST.'

                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


C --- PSA value out of range
                Else If ( Abs(Value_PSA(icounter)) .gt. MAX_ABS_VALUE_PSA) Then
                   error_flag = .true.
           
                   rtn = String_Loc(Name_PSA(icounter),fbytec,lbytec)

                   msgbuf = 'PSA value ' // msg25b(fbyteb:lbyteb) // ' is out of bounds.' // char(10) 
     2             // 'PSA ' // Name_PSA(icounter)(fbytec:lbytec) // ' not set.'
     3             // char(10) // 'Operator Action:  Notify SDST.'


                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  
C --- set PSA name 
                Else
                   AttrN = MCORE_ADD_ATTRIBUTENAME // '.' // msg25(fbyte:lbyte) 

                   rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,Name_PSA(icounter))

                   If (rtn.ne.PGS_S_SUCCESS) Then
                       error_flag = .true. 
                       rtn = String_Loc( AttrN,fbyte,lbyte)

                       msgbuf =
     1                 'pgs_met_setattr_s unable to set attribute ' // AttrN(fbyte:lbyte) // '.'
     2                 // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3                 // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4                 // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                       call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

                   End If


C --- Set PSA value 
                   AttrN = MCORE_PARAMETERVALUE // '.' // msg25(fbyte:lbyte) 
                   write(msg8,'(f8.2)') Value_PSA(icounter)

                   rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,msg8)

                   If (rtn.ne.PGS_S_SUCCESS) Then
                       error_flag = .true. 
                       rtn = String_Loc( AttrN,fbyte,lbyte)

                       msgbuf =
     1                 'pgs_met_setattr_s unable to set attribute ' // AttrN(fbyte:lbyte) // '.'
     2                 // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3                 // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4                 // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                       call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

                   End If       ! pgs_met_setAttr_s check
               End If        ! check for blank PSA name
            End IF       ! check for PSA NoSet_Flag

   50    continue    ! Loop on PSAs

      End If  ! check for PSA < 0



C-----------------------------------------------------------------
C Set "LocalGranuleID".  This is the product file name in the
C MODIS file naming convention.
C-----------------------------------------------------------------

      i=1
      loopexit=.FALSE.
      LRN_VRN=-1

      Do While(i.le.No_PCF_Attr .and. loopexit.eqv..FALSE.)
          runtimepar=Name_PCF_Attr(i)

          If(runtimepar(1:14).eq.'LOCALVERSIONID') Then
              LRN_VRN=LRN_PCF_Attr(i)
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
     2   'pgs_met_getsetattr_s unable to read the ESDT shortname from MCF.'
     3   // char(10) // 'Operator Action:  Check for wrong MCF file. If corrupted,'
     4   // char(10) // 'stage correct MCF and rerun PGE. Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      else


C .......Retrieve value of LOCALVERSIONID from PCF runtime parameter.
         rtn=pgs_pc_getconfigdata(LRN_VRN,LOCAL_VRSNID)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .true.
            write(msg25,'(I25)') LRN_VRN
            rtn = string_loc(msg25,fbyte,lbyte)

            msgbuf =
     2      'pgs_pc_getconfigdata unable to read localversionid'//' on LUN '// msg25(fbyte:lbyte) // '.'
     3      // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     4      // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     5      // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         else



c ..........Construct LOCALGRANULEID
            rtn = Gen_Modis_FileName(ESDT_Name,LOCAL_VRSNID,fname)

            If (rtn.ne.SUCCEED) Then
               error_flag = .true.
               msgbuf = 'Function Gen_Modis_FileName failed.'
     2         // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3         // char(10) // 'messages generated by function Gen_Modis_FileName.'

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else

c .............Set LOCALGRANULEID
               AttrN = MCORE_LOCALGRANULEID

               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,fname)

               If (rtn .ne. PGS_S_SUCCESS) Then
                  error_flag = .true.
                  rtn = String_Loc(AttrN,fbyte,lbyte)

                  msgbuf = 'PGS_MET_SetAttr_s unable to set '// AttrN(fbyte:lbyte) // '.'
     2            // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3            // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4            // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               End If
            End If
         End If
      End If



C-----------------------------------------------------------------
C Set Universal References (URs) or "InputPointers" to input data
C set.  URs are retrieved by reading field 5 of PCF file line 
C entries.
C-----------------------------------------------------------------

      NUM_nonBlank_UR = 0
      NUM_Pointers    = No_InputPntr

c - first check if number of input pointers exceed 
c - Set_CoreMetadata internal limit.

      If (NUM_Pointers .gt. MAX_NO_PNTRS) Then
         NUM_Pointers = MAX_NO_PNTRS
         error_flag = .true.
         write(msg25,'(i25)') NUM_Pointers 

         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'Number of Universal References exceeds ' // char(10) 
     2   // 'number of elements that can be set by ' // 'Set_CoreMetadata.' 
     3   // char(10) // 'Only ' // msg25(fbyte:lbyte) // ' URs set.'
     4   // char(10) // 'Operator Action: Notify SDST. '

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If

c-----loop over array passed from application
      Do 60 icounter = 1, NUM_Pointers 
         Version_No = Vrsn_InputPntr(icounter)
         LRN        = LRN_InputPntr(icounter)

c --- set message buffers for error reporting 

         write(msg25b,'(i25)') LRN 
         write(msg25c,'(i25)') Version_No 
         write(msg25d,'(i25)') icounter

         rtn = String_Loc(msg25b,fbyteb,lbyteb)
         rtn = String_Loc(msg25c,fbytec,lbytec)
         rtn = String_Loc(msg25d,fbyted,lbyted)

c--------ignore InputPointers with LUN entries equal to NO_SET_UR
         If (LRN .EQ. NO_SET_UR) Then
            continue

         Else If (Version_No .LT. 1)  Then
            error_flag = .true.

            msgbuf = 
     1       'Input Pointer PCF file version number is '  // msg25c(fbytec:lbytec) // '.' 
     2       // char(10) // 'Set_CoreMetadata unable to ' // 'set element ' 
     3       // msg25d(fbyted:lbyted) // ' of input UR array.' 
     4       // char(10) // 'Operator Action: Notify SDST. '

       
            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else If (LRN.LT.1 .OR. (LRN.GE.10000 .AND. LRN.LE.10999)) Then
            error_flag = .true.

            msgbuf = 'Input Pointer LRN number out of bounds.'
     1      // char(10) // 'Set_CoreMetadata unable to set element ' // msg25d(fbyted:lbyted) 
     2      // char(10) // ' of input UR array (LUN = ' // msg25b(fbyteb:lbyteb) // ').'
     3      // char(10) // 'Operator Action: Notify SDST. '

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else 
            rtn = PGS_PC_GetUniversalRef(LRN, Version_No, UR_Ref(icounter) )

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf = 'PGS_PC_GetUniversalRef unable to retrieve element '
     1         // msg25d(fbyted:lbyted) // ' of input UR array. '
     2         // char(10) // '(LUN = ' // msg25b(fbyteb:lbyteb)
     3         // ';  Version = ' // msg25c(fbytec:lbytec) // ').'
     4         // char(10) // 'Operator Action:  Check PCF Universal Reference identifier. '
     5         // char(10) // 'If invalid, stage correct PCF and rerun PGE.  Otherwise, '
     6         // char(10) // 'notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

            Else If (UR_Ref(icounter) .EQ. BLANK) Then
               error_flag = .TRUE.

               msgbuf =
     1         'UR on LUN ' // msg25b(fbyteb:lbyteb) // ' and VRSN ' // msg25c(fbytec:lbytec)
     2         // ' is blank. '
     3         // char(10) // 'Operator Action:  Check PCF Universal Reference identifier. '
     4         // char(10) // 'If invalid, stage correct PCF and rerun PGE.  Otherwise, '
     5         // char(10) // 'notify SDST.'

               Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            Else
               NUM_nonBlank_UR         = NUM_nonBlank_UR + 1
               UR_Ref(NUM_nonBlank_UR) = UR_Ref(icounter) 
            End If

         End If

   60 Continue

C - If sufficient space in UR_Ref work buffer, place end-of-data
C   marker into element NUM_nonBlank_U+1.

      If (NUM_nonBlank_UR .lt. MAX_NO_PNTRS)
     2   UR_Ref(NUM_nonBlank_UR+1) = PGSd_MET_STR_END


      AttrN = MCORE_INPUT_POINTER
      rtn   = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,UR_Ref)

      If (rtn.ne.PGS_S_SUCCESS) Then
          error_flag = .true. 
          rtn = String_Loc(AttrN,fbyte,lbyte)

          msgbuf = 'PGS_MET_SetAttr_s unable to set ' // AttrN(fbyte:lbyte) // '.'
     2    // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3    // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4    // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

          call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
      End If


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

             msgbuf = 'PGS_MET_SetAttr_s unable to set ' // AttrN(fbyte:lbyte) // '.'
     2       // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3       // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4       // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

              Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
          End If

c ---     Only 1 geolocation file expected 
          Version_No = 1

          rtn = PGS_PC_GetUniversalRef(LRN_GEO,Version_No,UR_Ref(1))

          If (rtn .ne. PGS_S_SUCCESS) Then
             error_flag = .true.

             write(msg25b,'(i25)') LRN_GEO 
             rtn = String_Loc(msg25b,fbyteb,lbyteb)

             write(msg25c,'(i25)') Version_No 
             rtn = String_Loc(msg25c,fbytec,lbytec)

             msgbuf = 'PGS_PC_GetUniversalRef unable to retrieve Geolocation UR.'
     2       // char(10) // 'LUN = ' // msg25b(fbyteb:lbyteb) // ', '
     3       // ' Version_No = '//msg25c(fbytec:lbytec)
     4       // char(10) // 'Operator Action:  Check PCF for valid Universal Reference '
     5       // char(10) // 'identifier.  If incorrect, stage correct PCF and rerun PGE.'
     6       // char(10) // 'Otherwise notify SDST.'

             Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

          Else
             AttrN = MCORE_ANCIL_POINTER//'.1'

             rtn = PGS_MET_SetAttr_s (MET_Handles(INVENTORYMETADATA),AttrN,UR_Ref(1))

             If (rtn.ne.PGS_S_SUCCESS) Then
                 error_flag = .true. 
                 rtn = String_Loc(AttrN,fbyte,lbyte)

                 msgbuf = 'PGS_MET_SetAttr_s unable to set ' // AttrN(fbyte:lbyte) // '.'
     2            // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3            // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4            // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                 call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

             End If
          End If
      End If


      If (.not.error_flag) Set_CoreMetadata = MAPIOK 

      Return

      End
