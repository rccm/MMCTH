      Integer Function Update_ArchMet_MOD06( LUN_NUM_Of_ArchiveRP_Pairs,
     &                                       NUM_Of_MODIS_InputFiles,
     &                                       LUN_MODIS_InputFiles,
     &                                       VRSN_MODIS_InputFiles,
     &                                       NUM_Of_ArchivePSA_SC,
     &                                       Name_ArchivePSA_SC,
     &                                       Value_ArchivePSA_SC,
     &                                       MET_Handles )

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
C   Update ECS and user-defined product specific "Archive" metadata 
C   previously written to the MODIS L2 cloud product (MOD06_L2) by one or 
C   more preceding processes.  In the context used here, set means 
C   to associate a value with a metadata parameter name via function calls 
C   to the ECS Science Data Processing Toolkit (SDPTK).
C
C   The specific metadata fields updated by Update_ArchMet_MOD06, their
C   origin, and the value of the Metadata Configuration File (MCF)
C   "Data Location" parameter are listed below.  
C      GEO = Geolocation product granule
C      PCF = Process Control File
C      RP  = PCF runtime Parameter (RP) 
C      PGE = Science code
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
C      7    ExclusionGringFlag                GEO           PGE
C      8    GringPointLatitude                GEO           PGE
C      9    GringPointLongitude               GEO           PGE
C     10    GringPointSequenceNo              GEO           PGE
C
C     Product Specific Attributes (PSAs) or Objects
C        (Examples for V2 MOD04 product follow)
C   -----------------------------------------------
C     11    LocalInputGranuleID        MODIS Input Products PGE
C     12    Algorithm_Version_Cloud_Top_Property_IR
C                                             PCF (RP)      PGE
C     13    Algorithm_Version_Cloud_Phase_IR  PCF (RP)      PGE
C     14    Algorithm_Version_Cloud_Property_VIS
C                                             PCF (RP)      PGE
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
C                Update_ArchMet_MOD06 to determine the actual number of
C                archive metadata objects appearing as PCF RPs.)
C     
C                If a process lists no archive metadata object
C                name/value pairs in the PCF, set
C                LUN_NUM_Of_ArchiveRP_Pairs to zero or a negative
C                number.  Update_ArchMet_MOD06 then will read no
C                RPs from the PCF.
C     
C                By convention, function Update_ArchMet_MOD06 will
C                assume that LUNs in the range from
C                LUN_NUM_Of_ArchiveRP_Pairs + 1 to
C                LUN_NUM_Of_ArchiveRP_Pairs + NUM_Of_ArchiveRP_Pairs
C                refer to RPs containing the -names- of "archive"
C                metadata objects.
C
C                Similarly, Update_ArchMet_MOD06 will assume that LUNs
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
C                  Algorithm_Version_Cloud_Top_Property_IR
C                  Algorithm_Version_Cloud_Phase_IR
C                  Algorithm_Version_Cloud_Property_VIS
C
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
C                 Update_InvMeta_MOD06.  Update_InvMeta_MOD06 returns the
C                 initialized array back to the calling routine.
C
C                 In FORTRAN, element 1 of the Met_Handles array refers
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
C        PGS_MET_Init              (libPGSTK.a)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        PGS_PC_GetConfigData      (libPGSTK.a)
C        PGS_PC_GetReference       (libPGSTK.a)
C
C    Named Constant:
C        MCORE_*               (mapi.inc)
C        MECS_CORE             (mapi.inc)
C        MODIS_E_GENERIC       (PGS_MODIS_39500.f)
C        MODIS_W_GENERIC       (PGS_MODIS_39500.f)
C        PGSd_PC_*             (PGS_PC.f)
C        PGSd_MET_*            (PGS_MET.f: included in mapi.inc)
C        PGS_S_SUCCESS         (PGS_SMF.f)
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
      parameter (FUNCNAME = 'Update_ArchMet_MOD06')

      character*1  BLANK
      parameter   (BLANK = ' ')

      integer    MAX_NUM_PNTRS,      MAX_NUM_RP_Pairs,    MAX_NUM_PSA 
      parameter (MAX_NUM_PNTRS=1000, MAX_NUM_RP_Pairs=50, MAX_NUM_PSA = 50)

      integer    LUN_GEO
      parameter (LUN_GEO=600000)

      integer    ARCHIVE
      parameter (ARCHIVE = 3)

      integer    FAIL,    SUCCEED
      parameter (FAIL=-1, SUCCEED=0)

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
      character*25 FORM, msg25, msg25a, msg25b, msg25c, msg25d
      character*60 AttrN, Field_Name(4)
      character*1012 msgbuf,msgbuf1,msgbuf2
      character*(PGSd_PC_VALUE_LENGTH_MAX) NUM_Of_ArchiveRP_Pairs 
      character*(PGSd_PC_VALUE_LENGTH_MAX) RP_Name, RP_Value
      character*(PGSd_PC_VALUE_LENGTH_MAX) buf_s

      double precision buf_dbl(4)

      integer LUN_N,LUN_V,NUM_RP,VRSN
      integer pgs_pc_getconfigdata, 
     2        pgs_met_setattr_d, pgs_met_setattr_s, pgs_met_setattr_i,
     3        pgs_met_getpcattr_d, pgs_met_getpcattr_s, 
     4        pgs_met_getpcattr_i, string_loc, Update_LclInputGranID_MOD06
      integer fbyte,fbytea,fbyteb,fbytec,fbyted,
     2        lbyte,lbytea,lbyteb,lbytec,lbyted,
     3        buf_i(4),i,IOS,rtn,
     4        rtn2

      real    PSA_V
      logical error_flag


C------------------------
C Initialization
C------------------------

      Update_ArchMet_MOD06 = FAIL
      error_flag = .false.

C----------------------------------------------------------------------
C Set ECS metadata objects that are specified as name/value pairs in 
C the USER DEFINED RUNTIME PARAMETERS section of the PCF.  
C----------------------------------------------------------------------

      If (LUN_NUM_Of_ArchiveRP_Pairs .gt. 0) Then

c-----read the number of archive RP name/value pairs from the PCF 
         write(msg25d,'(I25)') LUN_NUM_Of_ArchiveRP_Pairs
         rtn = string_loc(msg25d,fbyted,lbyted)

         rtn = pgs_pc_getconfigdata(LUN_NUM_Of_ArchiveRP_Pairs,
     2                              NUM_Of_ArchiveRP_Pairs)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .true.

            msgbuf = 
     1      'pgs_pc_getconfigdata unable to read the number' 
     2      // char(10) // 'of PCF RP name/value pairs on LUN ' // msg25d(fbyted:lbyted) 
     3      // char(10) // 'No archive metadata RP name/value pairs set by Update_ArchMet_MOD06.' 
     4      // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     5      // char(10) // 'entry.  If LUN is nonexistent or value is blank, stage correct '
     6      // char(10) // 'PCF and rerun PGE.  Otherwise, notify SDST.'

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
  
               msgbuf = 
     1         'FORTRAN internal read error converting NUM_Of_ArchiveRP_Pairs RP string '
     2         // char(10) // 'value (' // NUM_Of_ArchiveRP_Pairs(fbyte:lbyte) // ') to integer. '
     3         // char(10) // 'No archive metadata RP name/value pairs set by Update_ArchMet_MOD06.' 
     4         // char(10) // 'Operator Action:  Notify SDST.'

               Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----test for out of bounds NUM_Of_ArchiveRP_Pairs. 
            Else If (NUM_RP .le. 0) Then
               error_flag = .true.
               write(msg25a,'(I25)') NUM_RP
               rtn = string_loc(msg25a,fbytea,lbytea)

               msgbuf = 
     1         'NUM_Of_ArchiveRP_Pairs on LUN = ' // msg25d(fbyted:lbyted) // ' (' 
     2         // msg25a(fbytea:lbytea) // ') is less than or equal to 0.'
     3         // char(10) // 'No archive metadata RP name/value pairs set by Update_ArchMet_MOD06.' 
     4         // char(10) // 'Operator Action:  Check for valid PCF file.  If wrong or '
     5         // char(10) // 'corrupted, stage correct PCF and rerun PGE.  Otherwise, '
     6         // char(10) // 'notify SDST. '

               Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----test for too many RPs 
            Else If (NUM_RP .gt. MAX_NUM_RP_Pairs) Then 
               error_flag = .true.  
               write(msg25a,'(I25)') LUN_NUM_Of_ArchiveRP_Pairs
               rtn = string_loc(msg25a,fbytea,lbytea)

               msgbuf = 
     1         'Number of archive metadata RP name/value pairs '
     2         // 'exceeds Update_ArchMet_MOD06 internal limit.'
     3         // char(10) // 'RP string on LUN ' // msg25a(fbytea:lbytea) // ' has value of '
     4         // NUM_Of_ArchiveRP_Pairs(fbyte:lbyte)
     5         // char(10) // 'No archive metadata RP name/value pairs set by Update_ArchMet_MOD06.' 
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
                  write(msg25b,'(I25)') i
                  write(msg25c,'(I25)') LUN_V 
                  rtn2 = String_Loc(msg25a,fbytea,lbytea)
                  rtn2 = String_Loc(msg25b,fbyteb,lbyteb)
                  rtn2 = String_Loc(msg25c,fbytec,lbytec)
      

c-----set error flag if unable to retrieve RP archive metadata name 
                  If (rtn .ne. PGS_S_SUCCESS) Then 
                     error_flag = .true.

                     msgbuf = 
     1               'pgs_pc_getconfigdata unable to retrieve '
     2               // 'archive metadata attribute name on RP LUN ' // msg25a(fbytea:lbytea) 
     3               // char(10) // 'Archive metadata RP pair ' // msg25b(fbyteb:lbyteb) // ' not set.'  
     4               // char(10) // 'Operator Action:  Check for valid PCF file.  If wrong or '
     5               // char(10) // 'corrupted, stage correct PCF and rerun PGE.  Otherwise, '
     6               // char(10) // 'notify SDST. '

                     Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

                  Else If (RP_Name .eq. BLANK) Then
                     error_flag = .true.

                     msgbuf =
     1               'Archive metadata name on RP LUN ' // msg25a(fbytea:lbytea) // ' is blank.'
     2               // char(10) // 'Update_ArchMet_MOD06 unable to set an unknown attribute.'
     3               // char(10) // 'Archive metadata RP pair #' // msg25b(fbyteb:lbyteb) // ' not set.'
     4               // char(10) // 'Operator Action:  Check for valid PCF file.  If wrong or '
     5               // char(10) // 'corrupted, stage correct PCF and rerun PGE.  Otherwise, '
     6               // char(10) // 'notify SDST. '

                     Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

                  Else
                     rtn = pgs_pc_getconfigdata(LUN_V, RP_Value)

c-----set error flag if unable to retrieve RP archive metadata value
                     If (rtn .ne. PGS_S_SUCCESS) Then
                        error_flag = .true.

                        msgbuf =
     1                  'pgs_pc_getconfigdata unable to retrieve archive metadata value on RP LUN '
     2                  // msg25c(fbytec:lbytec)
     3                  // char(10) // 'Archive metadata RP pair #' // msg25b(fbyteb:lbyteb) // ' not set'
     4                  // char(10) // 'Operator Action:  Check for valid PCF file.  If wrong or '
     5                  // char(10) // 'corrupted, stage correct PCF and rerun PGE.  Otherwise, '
     6                  // char(10) // 'notify SDST. '

                        Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

                     Else
                        rtn = pgs_met_setattr_s(MET_Handles(ARCHIVE),RP_Name,RP_Value)

c-----set error flag if unable to set archive metadata name/value pair
                        If (rtn .ne. PGS_S_SUCCESS) Then
                           error_flag = .true.

                           rtn = String_Loc(RP_Name,fbyte,lbyte)
                           rtn = String_Loc(RP_Value,fbytec,lbytec)

                           msgbuf =
     1                     'pgs_met_setattr_s unable to set archive metadata RP name/value ' 
     2                     // 'pair #' // msg25b(fbyteb:lbyteb)
     3                     // char(10) // 'RP name is '// RP_Name(fbyte:lbyte) 
     4                     // char(10) // 'RP value is ' // RP_Value(fbytec:lbytec)
     5                     // char(10) // 'Operator Action:  Check for valid MCF file and correct PCF '
     6                     // char(10) // 'reference to MCF file.  If MCF/PCF files are wrong or '
     7                     // char(10) // 'corrupted, stage correct files and rerun PGE.  Otherwise, '
     8                     // char(10) // 'notify SDST.'


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
     
         msgbuf = 
     2   'Variable LUN_NUM_Of_ArchiveRP_Pairs (' // msg25(fbyte:lbyte) 
     3   // ') is not a positive integer.' 
     4   // char(10) // 'No archive metadata RP name/value pairs set '
     5   // char(10) // 'by Update_ArchMet_MOD06.' 

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,
     2        FUNCNAME)

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
      write(msg25b,'(I25)') VRSN
      rtn = String_Loc(msg25a,fbytea,lbytea)
      rtn = String_Loc(msg25b,fbyteb,lbyteb)

      Do 200 i = 1, 4
         AttrN = Field_Name(i)

         rtn = string_loc(AttrN,fbyte,lbyte)

c--------
         msgbuf1 = 'pgs_met_getpcattr_s unable to retrieve ECS attribute ' // AttrN(fbyte:lbyte)
     1   // char(10) // 'from MODIS geolocation product on LUN = ' // msg25a(fbytea:lbytea) 
     2   // ';  Version number = ' // msg25b(fbyteb:lbyteb)
     3   // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     4   // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     5   // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6   // char(10) // 'fault is identified, stage correct PCF/input file and '
     7   // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


         msgbuf2 = 'pgs_met_setattr_s unable to set G-Ring attribute ' // AttrN(fbyte:lbyte)
     1   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'


c-----Retrieve and set ExclusionGRingFlag (Y=inner, N=outer)
         If (i .eq. 1) Then
            rtn = pgs_met_getpcattr_s(LUN_GEO, VRSN, MECS_CORE, AttrN, buf_s)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf1,FUNCNAME)
            Else
               rtn = pgs_met_setattr_s(MET_Handles(ARCHIVE), AttrN, buf_s)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf2,FUNCNAME)
               End If
            End If


c-----Retrieve and set GRingPointLatitude and GRingPointLongitude
         Else If ( (i.eq.2) .or. (i.eq.3) ) Then
            rtn = pgs_met_getpcattr_d(LUN_GEO, VRSN, MECS_CORE, AttrN, buf_dbl)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf1,FUNCNAME)

            Else
               rtn = pgs_met_setattr_d(MET_Handles(ARCHIVE), AttrN, buf_dbl)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf2,FUNCNAME)
               End If
            End If

c-----Retrieve and set GRingPointSequenceNo
         Else If (i.eq.4) Then

            rtn = pgs_met_getpcattr_i(LUN_GEO, VRSN, MECS_CORE, AttrN, buf_i)

c-----------set error flag if unable to retrieve metadata value
            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf1,FUNCNAME)

            Else
               rtn = pgs_met_setattr_i(MET_Handles(Archive), AttrN, buf_i)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf2,FUNCNAME)
               End If

            End If

         End If

 200  Continue


C-----------------------------------------------------------------
C Update MOD06 LocalInputGranuleID 
C-----------------------------------------------------------------

      rtn = Update_LclInputGranID_MOD06(MET_Handles,
     1                                 NUM_Of_MODIS_InputFiles,
     2                                 LUN_MODIS_InputFiles,
     3                                 VRSN_MODIS_InputFiles)

      If (rtn .eq. FAIL) Then
         error_flag = .true.

         msgbuf = 'Update_LclInputGranID_MOD06 detected error incrementing ECS '
     1            // char(10) // 'LocalInputGranuleID metadata attribute. ' 
     2            // char(10) // 'Operator Action:  Refer to prior low level LogStatus ' 
     3            // char(10) // 'messages originating from call to routine ' 
     4            // 'Update_LclInputGranID_MOD06. ' 

      
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


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

         msgbuf = 
     1   'Number of archive PSAs (=' // msg25a(fbytea:lbytea) // ') exceeds Update_ArchMet_MOD06. ' 
     2   // char(10) // 'internal limit.  No PSAs set by Update_ArchMet_MOD06. '
     3   // char(10) // 'Operator Action:  Notify SDST.'

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
     1         // char(10) // 'PSA ' // msg25b(fbyteb:lbyteb) // ' not set.'
     2         // char(10) // 'Operator Action:  Notify SDST.'

               Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else
               write(msg8, '(F8.2)') PSA_V 
               AttrN = Name_ArchivePSA_SC(i)

c-----test for blank PSA name
               If (AttrN .eq. BLANK) Then
                  write(msg25,'(i25)') i
                  rtn = String_Loc(msg25,fbyte,lbyte)
                  error_flag = .true.

                  msgbuf = 'Archive PSA name is blank.  Update_ArchMet_MOD06 unable to set unknown PSA.'
     1            // char(10) // 'PSA number is ' // msg25(fbyte:lbyte)
     2            // char(10) // 'Operator Action:  Notify SDST.'
     
                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

C-----Set PSA value
               Else
                  rtn = pgs_met_setattr_s(MET_Handles(ARCHIVE), AttrN, msg8)

c-----set error flag if unable to set archive metadata
                  If (rtn.ne.PGS_S_SUCCESS) Then
                     error_flag = .true.
                     rtn = String_Loc(AttrN,fbyte,lbyte)

                     msgbuf = 
     1               'pgs_met_setattr_s unable to set archive PSA attribute ' // AttrN(fbyte:lbyte)
     2               // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3               // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4               // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                     Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
                  End If

               End If

            End If

 400     Continue

      End If

      If (.not.error_flag) Update_ArchMet_MOD06 = SUCCEED 

      Return

      End

C-----------------------------------------------------------------------

      Integer Function Update_LclInputGranID_MOD06(MET_Handles, 
     +                                             NUM_Of_MODIS_InputFiles,
     +                                             LUN_MODIS_InputFiles,
     +                                             VRSN_MODIS_InputFiles)

      implicit none
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Increment and set an updated list of MODIS input product 
C                file pointers to be stored as MOD06 product metadata (in 
C                Archive attribute LOCALINPUTGRANULEID) 
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
C  integer  NUM_Of_MODIS_InputFiles 
C                           Variable containing the number of MODIS input
C                           products used by process.
C
C  integer  LUN_MODIS_InputFiles(*)
C                           Array containing MODIS input product LUNs used by
C                           process (See also VRSN_MODIS_InputFiles below).
C
C  integer  VRSN_MODIS_InputFiles(*)
C                           Array containing file version numbers of MODIS input
C                           products (See LUN_MODIS_InputFiles above).
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
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Update_LclInputGranID_MOD06')

      integer ARCHIVEMETADATA
      Parameter (ARCHIVEMETADATA= 3)

      integer MAX_NO_PNTRS
      parameter (MAX_NO_PNTRS=50)

      integer LUN_MOD06
      Parameter (LUN_MOD06 = 412500)

      integer    FAIL,    SUCCEED
      parameter (FAIL=-1, SUCCEED=0)


C input argument declarations
      character*(*) MET_Handles(*)
      integer NUM_Of_MODIS_InputFiles
      integer LUN_MODIS_InputFiles(*)
      integer VRSN_MODIS_InputFiles(*)


C local variable declarations
      character*25     msg25_2,msg25_3,msg25_4,msg25_LUN,msg25_VRSN
      character*255    HDF_AttrName, ECS_AttrName
      character*(1024) msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) New_MODIS_Files(MAX_NO_PNTRS),
     1                                     MODIS_Files(MAX_NO_PNTRS) /MAX_NO_PNTRS*PGSd_MET_STR_END/ 

      integer rtn_flag, String_Loc, PGS_MET_GetPCAttr_s, PGS_MET_SetAttr_s
      integer i1, i2, rtn_loc, rtn, 
     1        fbyte,fbyte1,fbyte2,fbyte3,fbyte4,fbyte_LUN,fbyte_VRSN,
     2        lbyte,lbyte1,lbyte2,lbyte3,lbyte4,lbyte_LUN,lbyte_VRSN
      integer No_Duplicates,No_MODIS_Files,No_New_MODIS_Files,No_Pntr,No_Total_MODIS_Files,
     1        No_Unique,
     2        MODIS_Files_Array_Index,VRSN_MOD06

      logical error_flag, In_Existing_MODIS_File_Array 


c-----Initialization
      Update_LclInputGranID_MOD06 = SUCCEED
      error_flag = .false.
      No_Unique     = 0 
      No_Duplicates = 0
      VRSN_MOD06    = 1
      HDF_AttrName  = MECS_ARCHIVE
      ECS_AttrName  = MCORE_LOCALINPUTGRANULEID

c-----Set up variables for status messaging
      write(msg25_LUN,'(I25)') LUN_MOD06
      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

      write(msg25_VRSN,'(I25)') VRSN_MOD06 
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)
      
      rtn_loc = String_Loc(ECS_Attrname,fbyte1,lbyte1)


c-----read existing LocalInputGranID from MOD06 product file 
      rtn = PGS_MET_GetPCAttr_s(LUN_MOD06, VRSN_MOD06, HDF_AttrName, ECS_AttrName, MODIS_Files)  

      If (rtn .EQ. PGS_S_SUCCESS) Then
         No_MODIS_Files = 1

c--------search existing MODIS_Files array for PGSd_MET_STR_END 
         Do WHILE ( (MODIS_Files(No_MODIS_Files) .NE. PGSd_MET_STR_END) .AND.
     1              (No_MODIS_Files .LE. MAX_NO_PNTRS) ) 

            No_MODIS_Files = No_MODIS_Files + 1
         End Do

         No_MODIS_Files = No_MODIS_Files - 1
         No_Total_MODIS_Files = No_MODIS_Files

c--------retrieve LocalGranuleID set of MODIS files used in current process 
         No_Pntr = NUM_Of_MODIS_InputFiles
         If (NUM_Of_MODIS_InputFiles .GT. MAX_NO_PNTRS) No_Pntr = MAX_NO_PNTRS  

         Call Retr_LocalGranID( No_Pntr,
     1                          LUN_MODIS_InputFiles,
     2                          VRSN_MODIS_InputFiles,
     3                          No_New_MODIS_Files,
     4                          New_MODIS_Files,
     5                          rtn_flag )


         If (rtn_flag .EQ. SUCCEED) Then
 
            Do 100 i1 = 1, No_New_MODIS_Files
               In_Existing_MODIS_File_Array = .FALSE.

c--------------search against existing file name set for duplicate name  
               Do 200 i2 = 1, No_MODIS_Files

                  If (New_MODIS_Files(i1) .EQ. MODIS_Files(i2)) Then
                     In_Existing_MODIS_File_Array = .TRUE.
                     No_Duplicates = No_Duplicates + 1
                  End If

 200           Continue


c--------------check results of duplicate name search 
               If ( .NOT. In_Existing_MODIS_File_Array ) Then
                  No_Total_MODIS_Files = No_Total_MODIS_Files + 1
                  No_Unique   = No_Unique + 1

                  If (No_Total_MODIS_Files .LE. MAX_NO_PNTRS) Then 
                     MODIS_Files_Array_Index = No_Total_MODIS_Files
                     MODIS_Files(MODIS_Files_Array_Index) = New_MODIS_Files(i1)
                  Else
                     write(msg25_2,'(i25)') No_Total_MODIS_Files
                     rtn_loc = String_Loc(msg25_2,fbyte2,lbyte2) 

                     msgbuf = 'Total number of MODIS_Files ' // msg25_2(fbyte2:lbyte2)
     1               // ' exceeds '// FUNCNAME // ' internal buffer storage size.'  
     2               // char(10) // 'New MODIS_File ' //  New_MODIS_Files(i1)
     3               // ' not aggregated to ECS '// 'LocalInputGranuleID array.'

                     Call MODIS_SMF_SetDynamicMsg(MODIS_W_GENERIC,msgbuf,FUNCNAME)
                  End If  
               End If


100         Continue 

c-----------report the number of new and duplicate file references
            write(msg25_3,'(i25)') No_Unique
            rtn_loc = String_Loc(msg25_3,fbyte3,lbyte3)

            write(msg25_4,'(i25)') No_Duplicates
            rtn_loc = String_Loc(msg25_4,fbyte4,lbyte4)

            msgbuf = msg25_3(fbyte3:lbyte3) // ' unique file names added to '
     1      // 'LOCALINPUTGRANULEID array.  '
     2      // char(10) // msg25_4(fbyte4:lbyte4) // ' duplicate file names found.'
 
            Call MODIS_SMF_SetDynamicMsg(MODIS_A_GENERIC,msgbuf,FUNCNAME) 

         Else
            error_flag = .TRUE.

            msgbuf = 'Retr_LocalGranID detected error reading LocalGranuleIDs from MODIS input files.'
     1      // char(10) // 'ECS LocalInputGranuleID attribute not updated in MODIS Cloud Product. ' 
     2      // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3      // char(10) // 'messages originating from call to routine Retr_LocalGranID. '

 
            Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If
            

C--------If sufficient storage in MODIS_File buffer, add end-of-string marker to MODIS_File array 
         If ( No_Total_MODIS_Files .lt. MAX_NO_PNTRS) Then
            MODIS_Files_Array_Index = No_Total_MODIS_Files + 1
            MODIS_Files(MODIS_Files_Array_Index) = PGSd_MET_STR_END
         End If


         ECS_AttrName = MCORE_LOCALINPUTGRANULEID

         rtn = PGS_MET_SetAttr_s(MET_Handles(ARCHIVEMETADATA),ECS_AttrName,MODIS_Files)

         If (rtn .NE. PGS_S_SUCCESS) Then
             error_flag = .TRUE. 
             rtn_loc = String_Loc(ECS_AttrName,fbyte,lbyte)

             msgbuf = 'PGS_MET_SetAttr_s detected error setting attribute ' 
     1       // ECS_AttrName(fbyte:lbyte)
     2       // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3       // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4       // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 
         End If

      Else  
         error_flag = .TRUE.

         msgbuf = 'PGS_MET_GetPCAttr_s detected error reading ECS attribute '
     1   // ECS_Attrname(fbyte1:lbyte1) 
     2   // char(10) // 'from MODIS cloud product on LUN = ' // msg25_LUN(fbyte_LUN:lbyte_LUN) 
     3   // ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     4   // char(10) // 'LocalInputGranID not set to HDF product file. '
     5   // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     6   // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     7   // char(10) // 'fault is identified, stage correct PCF/input file and '
     8   // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

      End If   !  INPUTPOINTER read check

      If (error_flag) Update_LclInputGranID_MOD06 = FAIL 

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
C                Procedure return flag (SUCCEED = 0, FAIL=-1)  
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
C    Written by Rich Hucek,                     February 1998 
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (_1) if an error occurs
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C
C    Named Constant:
C        MCORE_*               (mapi.inc)
C        MECS_CORE             (mapi.inc)
C        MODIS_E_GENERIC       (PGS_MODIS_39500.f)
C        PGSd_PC_*             (PGS_PC.f)
C        PGSd_MET_*            (PGS_MET.f: included in mapi.inc)
C        PGS_S_SUCCESS         (PGS_SMF.f)
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
      parameter (FUNCNAME = 'Retr_LocalGranID')

      integer    FAIL,    SUCCEED
      parameter (FAIL=-1, SUCCEED=0)

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
      rtn_flag = SUCCEED 
      AttrN =  MCORE_LOCALGRANULEID 
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

  100 Continue

      Return

      End
