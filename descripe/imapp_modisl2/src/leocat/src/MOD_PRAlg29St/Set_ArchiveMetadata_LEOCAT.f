      Integer Function Set_ArchiveMetadata( LUN_NUM_Of_ArchiveRP_Pairs,
     &                                      NUM_Of_MODIS_InputFiles,
     &                                      LUN_MODIS_InputFiles,
     &                                      VRSN_MODIS_InputFiles,
     &                                      NUM_Of_ArchivePSA_SC,
     &                                      Name_ArchivePSA_SC,
     &                                      Value_ArchivePSA_SC,
     &                                      MET_Handles, NDVIfile,
     &                                      thresfile)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C   Sets EOSDIS Core System (ECS) and user-defined product specific
C   metadata parameter values that are to appear in MODIS Level 2 HDF
C   products as elements of the file attribute 'ArchiveMetadata.0'.  In
C   the context used here,  set means to associate a value with a
C   metadata parameter name via function calls to the ECS Science Data
C   Processing Toolkit (SDPTK).
C
C !DESCRIPTION:
C
C   Function Set_ArchiveMetadata sets ECS and product specific attribute
C   (PSA) values that are to be written to the HDF file attribute
C   ArchiveMetadata.0.  Note that archive PSAs are defined as needed
C   by data producers, and they are not registered as formal ECS ODL
C   (Object Description Language) data model objects as are the
C   AdditionalAttribute fields of the product inventory metadata,
C   which are also referred to as PSAs.
C
C   No metadata fields are actually written to the product file during
C   the call to Set_ArchiveMetadata; only the association of values to
C   fields is made at this time.  Actual insertion of HDF attribute
C   'ArchiveMetadata.0' into the HDF product file takes place later
C   with a call to the SDPTK routine PGS_MET_Write (by passing it
C   element 3 of array 'MET_Handles' as actual first argument).  This
C   call is made from outside Set_ArchiveMetadata after all fields
C   comprising 'ArchiveMetadata.0' have been set by the PGE.  If the
C   MODIS HDF Application Programming Interface (M-API) Utility is being
C   used, PGS_MET_Write is called automatically when the HDF product
C   file is closed with CPMFIL.  If one does not use M-API, then a
C   direct call to PGS_MET_Write is required to attach
C   'ArchiveMetadata.0' to the file before it is closed.
C
C   The specific metadata fields set by Set_ArchiveMetadata, their
C   origin, and the value of the Metadata Configuration File (MCF)
C   "Data Location" parameter are listed below.  Only a small number
C   of field values must actually be computed/determined by the
C   science code during granule processing (referred to here as the
C   PGE).  Most fields can be read in at run time from the MODIS
C   Geolocation (GEO) product granule) or preset as USER DEFINED
C   RUNTIME PARAMETERs (RPs) in the Process Control File (PCF).
C
C                                          Source of    Data Location
C     ECS Core Metadata Objects              Value      Value in MCF
C   -------------------------------        ---------    -------------
C      1    AlgorithmPackageAcceptanceDate    PCF (RP)      PGE
C      2    AlgorithmPackageMaturityCode      PCF (RP)      PGE
C      3    AlgorithmPackageName              PCF (RP)      PGE
C      4    AlgorithmPackageVersion           PCF (RP)      PGE
C      5    LongName (optional)               PCF (RP)      PGE
C      6    InstrumentName                    PCF (RP)      PGE
C      7    LocalInputGranuleID        MODIS Input Products PGE
C      8    ExclusionGringFlag                GEO           PGE
C      9    GringPointLatitude                GEO           PGE
C     10    GringPointLongitude               GEO           PGE
C     11    GringPointSequenceNo              GEO           PGE
C
C     Product Specific Attributes (PSAs) or Objects
C        (Examples for V2 MOD04 product follow)
C   -----------------------------------------------
C     12    AlgorithmSoftwareVersionLand      PCF (RP)      PGE
C     13    AlgorithmSoftwareVersionOcean     PCF (RP)      PGE
C     14    VeryGoodQualityDataPct_Land       PGE           PGE
C     15    GoodQualityDataPct_Land           PGE           PGE
C     16    MarginalQualityDataPct_Land       PGE           PGE
C     17    BadQualityDataPct_Land            PGE           PGE
C     18    VeryGoodQualityDataPct_Ocean      PGE           PGE
C     19    GoodQualityDataPct_Ocean          PGE           PGE
C     20    MarginalQualityDataPct_Ocean      PGE           PGE
C     21    BadQualityDataPct_Ocean           PGE           PGE
C
C
C !INPUT PARAMETERS:
C
C  integer LUN_NUM_Of_ArchiveRP_Pairs
C                Variable containing the PCF logical reference number
C                of USER DEFINED RUNTIME PARAMETER (RP),
C                "NUM_Of_ArchiveRP_Pairs".  "NUM_Of_ArchiveRP_Pairs"
C                is to be set to the INTEGER number of archive metadata
C                object name/value pairs listed as RPs in the PCF.
C                (Note that NUM_Of_ArchiveRP_Pairs is not itself an
C                archive metadata object.  It is used internally by
C                Set_ArchiveMetadata to determine the actual number of
C                archive metadata objects appearing as PCF RPs.)
C
C                If a process lists no archive metadata object
C                name/value pairs in the PCF, set
C                LUN_NUM_Of_ArchiveRP_Pairs to zero or a negative
C                number.  Set_ArchiveMetadata then will read no
C                RPs from the PCF.
C
C                By convention, function Set_ArchiveMetadata will
C                assume that LUNs in the range from
C                LUN_NUM_Of_ArchiveRP_Pairs + 1 to
C                LUN_NUM_Of_ArchiveRP_Pairs + NUM_Of_ArchiveRP_Pairs
C                refer to RPs containing the -names- of "archive"
C                metadata objects.
C
C                Similarly, Set_ArchiveMetadata will assume that LUNs
C                in the range from
C                NUM_Of_ArchiveRP_Pairs + NUM_Of_ArchiveRP_Pairs + 1
C                to NUM_Of_ArchiveRP_Pairs + 2*NUM_Of_ArchiveRP_Pairs
C                refer to RPs containing the -values- of "archive"
C                metadata objects.
C
C                A one-to-one correspondence between the "names" in
C                LUN range from LUN_NUM_Of_ArchiveRP_Pairs + 1 to
C                LUN_NUM_Of_ArchiveRP_Pairs + NUM_Of_ArchiveRP_Pairs
C                and the "values" in the range from
C                LUN_NUM_Of_ArchiveRP_Pairs + NUM_Of_ArchiveRP_Pairs + 1
C                to
C                LUN_NUM_Of_ArchiveRP_Pairs + 2*NUM_Of_ArchiveRP_Pairs
C                is assumed.
C
C                Archive metadata objects that should be listed as
C                PCF RPs are:
C                  AlgorithmPackageAcceptanceDate
C                  AlgorithmPackageMaturityCode
C                  AlgorithmPackageName
C                  AlgorithmPackageVersion
C                  LongName
C                  InstrumentName
C                  AlgorithmSoftwareVersionLand
C                  AlgorithmSoftwareVersionOcean
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
C  integer NUM_Of_ArchivePSA_SC
C                Variable containing the number of archive Product
C                Specific Attributes (PSA) set by science code (i.e.
C                SC).  Note that other archive PSAs are known prior
C                to running PGE and they may be defined as PCF RPs
C                (See LUN_NUM_Of_ArchiveRP_Pairs above).
C
C  character*(*) Name_ArchivePSA_SC(*)
C                Array containing the names of archive PSAs to be set
C                for granule by science code (See Value_ArchivePSA_SC below).
C
C                Note that archive PSAs are defined as needed by data
C                producers, and they are not registered as formal ECS
C                ODL (Object Description Language) data model objects
C                as are the AdditionalAttribute fields of the product
C                inventory metadata which are also referred to as PSAs.
C
C  real Value_ArchivePSA_SC(*)
C                Array containing the values of archive PSAs to be set
C                for granule by science code (See Name_ArchivePSA above).
C                Values are restricted to range from -9999.99 to
C                +99999.99 and are to be stored in the product as
C                a formatted, F8.2 ASCII character string.
C
C                A one-to-one relationship between the elements of
C                arrays Name_ArchivePSA_SC and Value_ArchivePSA_SC
C                is assumed.
C
C
C  character MET_Handles(20)
C                 Array containing the names of the "MASTER" groups as
C                 defined in the MCF (There may be up to 20 "MASTER"
C                 groups in the MCF.).  "Met_Handles" must be the same
C                 array as initialized previously in a call to function
C                 Set_CoreMetadata.  Set_CoreMetadata returns the
C                 initialized array back to the calling routine.
C
C                 In FORTRAN, element 1 of the Met_Handles array refers
C                 to the MCF file, element 2 to ECS inventory metadata,
C                 and element 3 to archive metadata.
C
C
C !OUTPUT PARAMETERS: None
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
C    Written by Rich Hucek, Fay Liang, and Liqun Ma     September 1997
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C    fhliang@ltpmail.gsfc.nasa.gov
C    lma@ltpmail.gsfc.nasa.gov
C
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED if successful, FAIL if an error occurs
C
C  Externals:
C
C    Functions and Subroutines:
C     MODIS_SMF_SETDYNAMICMSG  (atmos shared code)
C     Set_ArchPSA              (atmos shared code)
C     Set_GRing                (atmos shared code)
C     Set_LclInputGranID       (atmos shared code)
C     Set_RP_Data_Atmos        (atmos shared code)
C
C    Named Constant:
C     ARCHIVEMETADATA          (Atmos_ECSMET.inc)
C     FAIL                     (Atmos_ECSMET.inc)
C     MODIS_E_GENERIC          (PGS_MODIS_39500.f)
C     SUCCEED                  (Atmos_ECSMET.inc)
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*) FUNCNAME
      parameter (FUNCNAME = 'Set_ArchiveMetadata')


c-----Declaration of function arguments
      character*(*) MET_Handles(*),Name_ArchivePSA_SC(*)
      character*34 ndvifile
      character*23 thresfile

      integer LUN_NUM_Of_ArchiveRP_Pairs
      integer NUM_Of_MODIS_InputFiles
      integer LUN_MODIS_InputFiles(*)
      integer VRSN_MODIS_InputFiles(*)
      integer NUM_Of_ArchivePSA_SC

      real Value_ArchivePSA_SC(*)


c-----Declaration of local variables
      character*512 msgbuf

      integer rtn

      integer Set_ArchPSA,
     2        Set_GRing,
     3        Set_LclInputGranID,
     4        Set_RP_Data_Atmos

      logical error_flag



C------------------------
C Initialization
C------------------------

      Set_ArchiveMetadata = FAIL
      error_flag = .false.


C----------------------------------------------------------------------
C Set ECS metadata objects that are specified as name/value pairs in
C the USER DEFINED RUNTIME PARAMETERS section of the PCF.
C----------------------------------------------------------------------

      rtn = Set_RP_Data_Atmos( Met_Handles,
     2                         LUN_NUM_Of_ArchiveRP_Pairs,
     3                         ARCHIVEMETADATA )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_RP_Data_Atmos detected error setting ECS archive RP data.'
     2   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     3   // char(10) // 'messages originating in function Set_RP_Data_Atmos.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C----------------------------------------------------------------------
C Set ECS GPolygon group attributes.
C----------------------------------------------------------------------

      rtn = Set_GRing(MET_Handles, ARCHIVEMETADATA)

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_GRing detected error setting ECS GPolygon group attributes.'
     2   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     3   // char(10) // 'messages originating in function Set_GRing.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C----------------------------------------------------------------------
C Set LOCALINPUTGRANULEID.  This is an array containing the names of
C all MODIS product files (i.e. product's ECS LOCALGRANULEID metadata
C field) used by the process.  File names are stored in the MODIS
C naming convention.
C----------------------------------------------------------------------

      rtn = Set_LclInputGranID( MET_Handles,
     2                          NUM_Of_MODIS_InputFiles,
     3                          LUN_MODIS_InputFiles,
     4                          VRSN_MODIS_InputFiles, ndvifile, 
     5                          thresfile)

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_LclInputGranID detected error incrementing ECS '
     2            // char(10) // 'LocalInputGranuleID metadata attribute.'
     3            // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     4            // char(10) // 'messages originating from call to routine '
     5            // 'Set_LclInputGranID.'


         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



C----------------------------------------------------------------------
C Set archive Product Specific Attributes (PSAs).  PSA values are
C stored in the product as formatted, F8.2 ASCII character strings
C----------------------------------------------------------------------

      rtn = Set_ArchPSA( MET_Handles,
     2                   NUM_Of_ArchivePSA_SC,
     3                   Name_ArchivePSA_SC,
     4                   Value_ArchivePSA_SC )

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Set_ArchPSA detected error setting archive PSA attribute.'
     2            // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3            // char(10) // 'messages originating from call to function '
     4            // 'Set_ArchPSA.'

         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf



      If (.not.error_flag) Set_ArchiveMetadata = SUCCEED

      Return

      End
