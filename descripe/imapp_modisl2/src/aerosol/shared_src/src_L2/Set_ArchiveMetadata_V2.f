      Integer Function Set_ArchiveMetadata( LUN_NUM_Of_ArchiveRP_Pairs,
     &                                      NUM_Of_MODIS_InputFiles,
     &                                      LUN_MODIS_InputFiles,
     &                                      VRSN_MODIS_InputFiles,
     &                                      NUM_Of_ArchivePSA_SC,
     &                                      Name_ArchivePSA_SC,
     &                                      Value_ArchivePSA_SC,
     &                                      MET_Handles )

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
C                object name/value pairs listed as as RPs in the PCF.
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
C                for granule by science code (See
C                Value_ArchivePSA_SC below).  
C
C                Note that archive PSAs are defined as needed by data 
C                producers, and they are not registered as formal ECS 
C                ODL (Object Description Language) data model objects 
C                as are the AdditionalAttribute fields of the product 
C                inventory metadata which are also referred to as PSAs. 
C
C  Real Value_ArchivePSA_SC(*)
C                Array containing the values of archive PSAs to be set
C                for granule by science code (See Name_PSA above).
C                Values are restricted to range from -9999.99 to
C                +99999.99 and are to be stored in the product as 
C                a formatted, F8.2 ASCII character string.  
C
C                A one-to-one relationship between the elements of
C                arrays Name_ArchivePSA_SC and Value_ArchivePSA_SC
C                is assumed.
C
C
C  character Met_Handles(20)
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
C !OUTPUT PARAMETERS: 
C
C  None 
C
C
C !REVISION HISTORY:
c Revision 1.4  1998/07/23  13:22:37  rhucek
c Changed prolog to list array Met_Handles(20) under Input Parameters.  There
c are no Output Parameters.
c
c Revision 1.3  1998/06/08  20:42:52  fhliang
c Update error messages with "Operator Action" strings;
c Modified error messages.
c
c Revision 1.2  1997/09/24  11:28:36  rhucek
c Changed the declaration of variables buf_dbl and buf_i to
c 4 element arrays, buf_dbl(4) and buf_i(4).
c
c Revision 1.1  1997/09/11  21:03:29  rhucek
c Initial revision
c
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
C        PGS_MET_Init              (libPGSTK.a)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        PGS_PC_GetConfigData      (libPGSTK.a)
C        PGS_PC_GetReference       (libPGSTK.a)
C
C    Named Constant:
C        FAIL                  (hdf.inc)
C        MCORE_*               (mapi.inc)
C        MECS_CORE             (mapi.inc)
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
C        String_Loc
C        CONCATENATE.f
C
C    Variables:
C
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*) FUNCNAME
      parameter (FUNCNAME = 'Set_ArchiveMetadata')

      character*1 BLANK
      parameter (BLANK = ' ')

      integer MAX_NUM_PNTRS, MAX_NUM_RP_Pairs, MAX_NUM_PSA 
      parameter (MAX_NUM_PNTRS=1000, MAX_NUM_RP_Pairs=50,
     2           MAX_NUM_PSA = 50)

      integer LUN_GEO
      parameter (LUN_GEO = 600000)

      integer ARCHIVE
      parameter (ARCHIVE = 3)


c-----Declaration of function arguments
      character*(*) MET_Handles(*),Name_ArchivePSA_SC(*)

      integer LUN_NUM_Of_ArchiveRP_Pairs
      integer NUM_Of_MODIS_InputFiles
      integer LUN_MODIS_InputFiles(*)
      integer VRSN_MODIS_InputFiles(*)
      integer NUM_Of_ArchivePSA_SC

      real Value_ArchivePSA_SC(*)


c-----Declaration of local variables
      character*8 msg8
      character*25 FORM, msg25, msg25a, msg25b, msg25c
      character*60 AttrN, Field_Name(4)
      character*512 msgbuf,msgbuf1,msgbuf2
      character*100 MODIS_FileName(MAX_NUM_PNTRS)
      character*(PGSd_PC_VALUE_LENGTH_MAX) NUM_Of_ArchiveRP_Pairs 
      character*(PGSd_PC_VALUE_LENGTH_MAX) RP_Name, RP_Value
      character*(PGSd_PC_VALUE_LENGTH_MAX) buf_s

      double precision buf_dbl(4)

      integer LUN,LUN_N,LUN_V,NUM_RP,VRSN
      integer pgs_pc_getconfigdata, 
     2        pgs_met_setattr_d, pgs_met_setattr_s, pgs_met_setattr_i,
     3        pgs_met_getpcattr_d, pgs_met_getpcattr_s, 
     4        pgs_met_getpcattr_i, string_loc
      integer fbyte,fbytea,fbyteb,fbytec,
     2        lbyte,lbytea,lbyteb,lbytec,
     3        buf_i(4),i,icounter,IOS,NUM_Pointers,rtn,rtn2

      real    PSA_V
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

      If (LUN_NUM_Of_ArchiveRP_Pairs .gt. 0) Then

c-----read the number of archive RP name/value pairs from the PCF 
         rtn = pgs_pc_getconfigdata(LUN_NUM_Of_ArchiveRP_Pairs,NUM_Of_ArchiveRP_Pairs)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .true.
            write(msg25,'(I25)') LUN_NUM_Of_ArchiveRP_Pairs
            rtn = string_loc(msg25,fbyte,lbyte)

            msgbuf = 'pgs_pc_getconfigdata unable to read the number'
     2      // char(10) // 'of PCF RP name/value pairs on LUN ' // msg25(fbyte:lbyte) // '.'
     3      // char(10) // 'No archive metadata RP name/value pairs set by Set_ArchiveMetadata.' 
     4      // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     5      // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     6      // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'

            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else

c-----construct integer FORMAT code to convert RP string to integer value 
            write(msg25,'(I25)') PGSd_PC_VALUE_LENGTH_MAX 
            rtn = string_loc(msg25,fbyte,lbyte)
            FORM = '(I' // msg25(fbyte:lbyte) // ')'

            read(NUM_Of_ArchiveRP_Pairs,FORM,IOSTAT=IOS) NUM_RP

            rtn = string_loc(NUM_Of_ArchiveRP_Pairs,fbyte,lbyte)

c-----If blank string, reset first and last non-blank bytes to 1
            If (rtn .eq. FAIL) Then
               fbyte = 1
               lbyte = 1
            End If

c-----test result of internal read for data conversion error 
            If (IOS .ne. 0) then
               error_flag = .true.
  
               msgbuf = 'FORTRAN internal read error converting integer value of '
     2         // char(10) // 'RP string NUM_Of_ArchiveRP_Pairs.'
     3         // char(10) // 'NUM_Of_ArchiveRP_Pairs=' // NUM_Of_ArchiveRP_Pairs(fbyte:lbyte)
     4         // char(10) // 'No archive metadata RP name/value pairs set by Set_ArchiveMetadata.' 
     5         // char(10) // 'Operator Action:  Notify SDST.'

               Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----test for out of bounds NUM_Of_ArchiveRP_Pairs. 
            Else If (NUM_RP .le. 0) Then
               error_flag = .true.
               write(msg25a,'(I25)') NUM_RP
               rtn = string_loc(msg25a,fbytea,lbytea)

               msgbuf = 'NUM_Of_ArchiveRP_Pairs (' // msg25a(fbytea:lbytea) // ') is out of bounds.'
     2         // char(10) // 'No archive metadata RP name/value pairs set by Set_ArchiveMetadata.' 
     3         // char(10) // 'Operator Action:  Notify SDST.'

               Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----test for too many RPs 
            Else If (NUM_RP .gt. MAX_NUM_RP_Pairs) Then 
               error_flag = .true.  
               write(msg25a,'(I25)') LUN_NUM_Of_ArchiveRP_Pairs
               rtn = string_loc(msg25a,fbytea,lbytea)

               msgbuf = 'Number of archive metadata RP name/value pairs exceeds '
     2         // char(10) // 'Set_ArchiveMetadata internal limit.'
     3         // char(10) // 'RP string on LUN ' // msg25a(fbytea:lbytea) // ' has value of '
     4         // NUM_Of_ArchiveRP_Pairs(fbyte:lbyte) // '.'
     5         // char(10) // 'No archive metadata RP name/value pairs set by Set_ArchiveMetadata.'
     6         // char(10) // 'Operator Action:  Check for valid PCF file.  If wrong or '
     7         // char(10) // 'corrupted, stage correct PCF and rerun PGE.  Otherwise, '
     8         // char(10) // 'notify SDST. '

               Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----Valid number of RPs...Set RP pairs
            Else

c-----loop over and set archive RP name/value pairs.
               DO 100 i = 1, NUM_RP
                  LUN_N = LUN_NUM_Of_ArchiveRP_Pairs + i
                  LUN_V = LUN_NUM_Of_ArchiveRP_Pairs + i + NUM_RP

                  rtn = pgs_pc_getconfigdata(LUN_N, RP_Name)

                  write(msg25a,'(I25)') LUN_N 
                  rtn2 = String_Loc(msg25a,fbytea,lbytea)
                  write(msg25b,'(I25)') i
                  rtn2 = String_Loc(msg25b,fbyteb,lbyteb)
                  write(msg25c,'(I25)') LUN_V 
                  rtn2 = String_Loc(msg25c,fbytec,lbytec)
      

c-----set error flag if unable to retrieve RP archive metadata name 
                  If (rtn .ne. PGS_S_SUCCESS) Then 
                     error_flag = .true.

                     msgbuf = 'pgs_pc_getconfigdata unable to retrieve archive metadata field'
     2               // char(10) // 'name on RP LUN '// msg25a(fbytea:lbytea) 
     3               // '. Archive metadata RP pair '// msg25b(fbyteb:lbyteb) // ' not set.'  
     4               // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     5               // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     6               // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'
   
                     Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

                  Else If (RP_Name .eq. BLANK) Then
                     error_flag = .true.

                     msgbuf = 'Archive metadata name on RP LUN ' // msg25a(fbytea:lbytea) 
     2               // ' is blank.' // char(10)// 'Set_ArchiveMetadata unable to set '
     3               // 'an unknown attribute.' // char(10) // 'Archive metadata RP pair ' 
     4               // msg25b(fbyteb:lbyteb) // ' not set.'
     5               // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter '
     6               // char(10) // 'value.  If blank, stage correct PCF and rerun PGE.'
     7               // char(10) // 'Otherwise, notify SDST.'

                     Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

                  Else
                     rtn = pgs_pc_getconfigdata(LUN_V, RP_Value)

c-----set error flag if unable to retrieve RP archive metadata value
                     If (rtn .ne. PGS_S_SUCCESS) Then
                        error_flag = .true.

                        msgbuf = 'pgs_pc_getconfigdata unable to retrieve archive metadata value '
     2                  // char(10) // 'on RP LUN ' // msg25c(fbytec:lbytec) // '.'
     3                  // char(10) // 'Archive metadata RP pair ' // msg25b(fbyteb:lbyteb) // ' not set.'
     4                  // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     5                  // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     6                  // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'

                        Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

                     Else
                        rtn = pgs_met_setattr_s(MET_Handles(ARCHIVE),RP_Name, RP_Value)

c-----set error flag if unable to set archive metadata name/value pair
                        If (rtn .ne. PGS_S_SUCCESS) Then
                           error_flag = .true.

                           rtn = String_Loc(RP_Name,fbyte,lbyte)
                           rtn = String_Loc(RP_Value,fbytec,lbytec)

                           msgbuf = 'pgs_met_setattr_s unable to set archive metadata RP name/value '
     2                     // 'pair ' // msg25b(fbyteb:lbyteb) // '.'
     3                     // char(10) // 'Metadata object name is ' // RP_Name(fbyte:lbyte) // '.'
     4                     // char(10) // 'Metadata object value is ' // RP_Value(fbytec:lbytec) // '.'
     5                     // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     6                     // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     7                     // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                           Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

                        End If !---pgs_met_setattr_s check
                     End If !---pgs_pc_getconfigdata check: RP Value
                  End If !---pgs_pc_getconfigdata check: RP name

 100           Continue

            End If !---check for blank NUM_Of_ArchiveRP_Pairs
         End If !-----pgs_pc_getconfigdata check  

c-----Issue warning if LUN_NUM_Of_ArchiveRP_Pairs is < 1
      Else 
         write(msg25,'(I25)') LUN_NUM_Of_ArchiveRP_Pairs
         rtn = string_loc(msg25,fbyte,lbyte)
     
         msgbuf = 'LUN_NUM_Of_ArchiveRP_Pairs (LUN = ' // msg25(fbyte:lbyte) 
     2   // ') is not a positive integer.' 
     3   // char(10) // 'No archive metadata RP name/value pairs set by Set_ArchiveMetadata.' 

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)

      End IF !-----end bounds checking on LUN_NUM_Of_ArchiveRP_Pairs



C----------------------------------------------------------------------
C Set ECS GPolygon group attributes.  These include the four metadata
C objects:
C
C 1 - ExclusionGRingFlag 
C 2 - GRingPointLatitude
C 3 - GRingPointLongitude
C 4 - GRingPointSequenceNo
C----------------------------------------------------------------------

      Field_Name(1) = MCORE_EXCLUS_GRING_FLG
      Field_Name(2) = MCORE_GRING_POINT_LAT
      Field_Name(3) = MCORE_GRING_POINT_LON
      Field_Name(4) = MCORE_GRING_POINT_NUM

      VRSN = 1 

c-----set message buffers for error reporting
      write(msg25a,'(I25)') LUN_GEO
      rtn = String_Loc(msg25a,fbytea,lbytea)
      write(msg25b,'(I25)') VRSN
      rtn = String_Loc(msg25b,fbyteb,lbyteb)

      Do 200 i = 1, 4
         AttrN = Field_Name(i)

         rtn = string_loc(AttrN,fbyte,lbyte)

         msgbuf1 = 'pgs_met_getpcattr unable to retrieve '// AttrN(fbyte:lbyte) // char(10)
     2   // 'from MODIS geolocation product (LUN = ' // msg25a(fbytea:lbytea) 
     3   // ',  Version = ' // msg25b(fbyteb:lbyteb) // ').'
     4   // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     5   // char(10) // 'Operator Action:  Check for corrupted Geolocation or PCF '
     6   // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     7   // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

         msgbuf2 =  'pgs_met_setattr unable to set G-Ring attribute ' // AttrN(fbyte:lbyte)
     2   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'


c-----Retrieve and set ExclusionGRingFlag (Y=inner, N=outer)
         If (i .eq. 1) Then
            rtn = pgs_met_getpcattr_s(LUN_GEO, VRSN, MECS_CORE,AttrN, buf_s)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf1,FUNCNAME)
            Else
               rtn = pgs_met_setattr_s(MET_Handles(ARCHIVE), AttrN,buf_s)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf2,FUNCNAME)
               End If
            End If


c-----Retrieve and set GRingPointLatitude and GRingPointLongitude
         Else If ( (i.eq.2) .or. (i.eq.3) ) Then
            rtn = pgs_met_getpcattr_d(LUN_GEO, VRSN, MECS_CORE,AttrN, buf_dbl)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf1,FUNCNAME)
            Else
               rtn = pgs_met_setattr_d(MET_Handles(ARCHIVE), AttrN,buf_dbl)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf2,FUNCNAME)
               End If
            End If

c-----Retrieve and set GRingPointSequenceNo
         Else If (i.eq.4) Then

            rtn = pgs_met_getpcattr_i(LUN_GEO, VRSN, MECS_CORE,AttrN, buf_i)

c--------set error flag if unable to retrieve metadata value
            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf1,FUNCNAME)
            Else
               rtn = pgs_met_setattr_i(MET_Handles(Archive), AttrN,buf_i)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf2,FUNCNAME)
               End If
            End If
         End If

 200  Continue



C----------------------------------------------------------------------
C Set LOCALINPUTGRANULEID.  This is an array containing the names of 
C all MODIS product files (i.e. product's ECS LOCALGRANULEID metadata 
C field) used by the process.  File names are stored in the MODIS 
C naming convention.
C----------------------------------------------------------------------

      AttrN = MCORE_LOCALGRANULEID
      NUM_Pointers = NUM_Of_MODIS_InputFiles 

c-----first check if number of MODIS input files exceeds
c-----Set_ArchiveMetadata internal limit.

      If (NUM_Pointers .gt. MAX_NUM_PNTRS) Then
         NUM_Pointers = MAX_NUM_PNTRS
         error_flag = .true.
         write(msg25,'(i25)') NUM_Pointers

         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'Number of MODIS input files exceeds ' // char(10) 
     2   // 'Set_ArchiveMetadata internal file limit.' // char(10) // 'Only '
     3   // msg25(fbyte:lbyte)  // ' element of LOCALINPUTGRANULEID array set.'  
     4   // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If

c-----Loop over the number of MODIS input files
      Do 300 icounter = 1, NUM_Pointers
         VRSN = VRSN_MODIS_InputFiles(icounter)
         LUN = LUN_MODIS_InputFiles(icounter)

c-----set message buffers for error reporting 

         write(msg25a,'(i25)') LUN 
         write(msg25b,'(i25)') VRSN 
         write(msg25c,'(i25)') icounter

         rtn = String_Loc(msg25c,fbytec,lbytec)

         If (VRSN .LT. 1)  Then
            error_flag = .true.
            MODIS_FileName(icounter) = BLANK 

            msgbuf = 'MODIS input file PCF version number (' // msg25b(fbyteb:lbyteb) // ')' 
     2      // ' is less than one.'
     3      // char(10) // 'LOCALINPUTGRANULEID array element ' // msg25c(fbytec:lbytec)
     4      // ' is set to blank.'
     5      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else If (LUN .LE. 0) Then 
            error_flag = .true.
            MODIS_FileName(icounter) = BLANK 

            msgbuf = 
     1      'MODIS input file LUN number ' // msg25a(fbytea:lbytea) // ' is out of bounds.' 
     2      // char(10) // 'LOCALINPUTGRANULEID array element ' // msg25c(fbytec:lbytec) 
     3      // ' is set to blank.'
     4      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else 
            rtn = pgs_met_getpcattr_s(LUN, VRSN, MECS_CORE,AttrN, MODIS_FileName(icounter))

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.
               MODIS_FileName(icounter) = BLANK 

               msgbuf = 'pgs_met_getpcattr_s unable to retrieve LOCALINPUTGRANULEID' 
     2         // char(10) // 'from MODIS product (LUN = ' // msg25a(fbytea:lbytea)
     3         // ', Version = ' // msg25b(fbyteb:lbyteb) // ').' // char(10) 
     4         // 'LOCALINPUTGRANULEID array element ' // msg25c(fbytec:lbytec) // ' is set to blank.'
     5         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     6         // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     7         // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

            End If
         End If

  300 Continue

c-----If sufficient space in MODIS_FileName work buffer, place end-of-data
c-----marker into element NUM_Pointers+1.

      If (NUM_Pointers .lt. MAX_NUM_PNTRS)
     2   MODIS_FileName(NUM_Pointers+1) = PGSd_MET_STR_END

      AttrN = MCORE_LOCALINPUTGRANULEID
      rtn = PGS_MET_SetAttr_s(MET_Handles(ARCHIVE),AttrN,MODIS_FileName)

      If (rtn.ne.PGS_S_SUCCESS) Then
          error_flag = .true. 
          rtn = String_Loc(AttrN,fbyte,lbyte)

     
          msgbuf = 'PGS_MET_SetAttr_s unable to set ' // AttrN(fbyte:lbyte) // '.'
     2    // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3    // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4    // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

          Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
      End If



C----------------------------------------------------------------------
C Set archive Product Specific Attributes (PSAs).  PSA values are
C stored in the product as formatted, F8.2 ASCII character strings 
C                                                                
C----------------------------------------------------------------------

c-----test for too many archive PSAs 

      If (NUM_Of_ArchivePSA_SC .gt. MAX_NUM_PSA) Then
         error_flag = .true.
         write(msg25a,'(I25)') NUM_Of_ArchivePSA_SC 
         rtn = string_loc(msg25a,fbytea,lbytea)

         msgbuf = 'Number of archive PSAs (' // msg25a(fbytea:lbytea) 
     2   // ') exceeds Set_ArchiveMetadata' // char(10) // 'internal limit.  Error assumed.'
     3   // char(10) // 'No PSAs set by Set_ArchiveMetadata.'
     4   // char(10) // 'Operator Action:  Notify SDST.'


         Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      Else

         DO 400 i = 1, NUM_Of_ArchivePSA_SC
            PSA_V = Value_ArchivePSA_SC(i)

            write(msg25a,'(E25.8)') PSA_V
            write(msg25b,'(I25)') i
            rtn = string_loc(msg25a,fbytea,lbytea)
            rtn = string_loc(msg25b,fbyteb,lbyteb)


c-----test for in-range PSA value 
            If (PSA_V.gt.99999.99 .OR. PSA_V.lt.-9999.99) Then 
               error_flag = .true.

               msgbuf = 'Archive PSA value ' // msg25a(fbytea:lbytea) // ' is out of range.' 
     2         // char(10) // 'PSA ' // msg25b(fbyteb:lbyteb) // ' not set.'
     3         // char(10) // 'Operator Action:  Notify SDST.'

               Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else
               write(msg8, '(F8.2)') PSA_V 
               AttrN = Name_ArchivePSA_SC(i)

c-----test for blank PSA name
               If (AttrN .eq. BLANK) Then
                  write(msg25,'(i25)') i
                  rtn = String_Loc(msg25,fbyte,lbyte)
                  error_flag = .true.

                  msgbuf = 'Archive PSA name is blank.' // char(10)
     2            // 'Set_ArchiveMetadata unable to set unknown PSA.'
     3            // char(10) // 'PSA number is ' // msg25(fbyte:lbyte) // '.'
     4            // char(10) // 'Operator Action:  Notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

C-----Set PSA value
               Else
                  rtn = pgs_met_setattr_s(MET_Handles(ARCHIVE),AttrN, msg8)

c-----set error flag if unable to set archive metadata
                  If (rtn.ne.PGS_S_SUCCESS) Then
                     error_flag = .true.
                     rtn = String_Loc(AttrN,fbyte,lbyte)

                     msgbuf = 'pgs_met_setattr_s unable to set archive PSA field '
     2               // char(10) // AttrN(fbyte:lbyte)
     3               // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     4               // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     5               // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If
               End If
            End If

 400     Continue

      End If


      If (.not.error_flag) Set_ArchiveMetadata = SUCCEED 

      Return

      End
