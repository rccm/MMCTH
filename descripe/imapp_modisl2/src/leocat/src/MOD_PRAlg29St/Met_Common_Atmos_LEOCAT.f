      SUBROUTINE CUT_NAME(string, fbyte, lbyte)
      IMPLICIT NONE

C----------------------------------------------------------------------
C !F77
C
C !Description:
C
C    Find the position (first and last bytes) of the file name within a
C    string buffer containing both path (file location) and file name.
C
C !INPUT PARAMETERS:
C
C    character string  A string variable of arbitrary length which
C                      shall consist of a unix path and file name.
C
C !OUTPUT PARAMETERS:
C
C    integer fbyte     The byte location of the first nonblank character
C                      of the input string.
C    integer lbyte     The byte location of the last nonblank character
C                      of the input string.
C
C !REVISION HISTORY:
c Revision 1.2  1997/02/25  18:55:14  vlin
c Changed from function to a subroutine
C
C !TEAM-UNIQUE HEADER:
C
C    This software was developed by the MODIS Science Data Support Team
C    (SDST) for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS
C
C    Written by Richard Hucek
C
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    There is no check to identify the unpredictable value of an
C    undefined string.  Consequently, users must take care to initialize
C    all string variables before passing them to CUT_NAME.
C
C !END
C-----------------------------------------------------------------------

      CHARACTER*(*) string
      INTEGER fbyte, lbyte, string_len

C Initialize variables
      string_len=len(string)
      lbyte = string_len

C Determine byte position of last non-blank, non-slash character.
C This is last character in file name
      DO WHILE ( (string(lbyte:lbyte).eq.' ') .and. (lbyte.ge.1) )
         lbyte=lbyte-1
      END DO

C Determine byte position of last slash (/) in character string
      fbyte = lbyte - 1

      DO WHILE ( (fbyte .ge. 1) .and. (string(fbyte:fbyte) .ne. '/') )
         fbyte=fbyte-1
      END DO

      fbyte = fbyte + 1

      RETURN
      END



c-----------------------------------------------------------------------
      subroutine Get_LUN_Of_LclVrsnID_Atmos( LUN_Of_NUM_Of_RP_Pairs,
     1                                       LUN_Of_LclVrsnID,
     2                                       rtn_flag)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Identify PCF Runtime Parameter (RP) LUN containing
C                the value of the LOCALVERSIONID (See design notes).
C
C
C !INPUT PARAMETERS:
C  integer LUN_Of_NUM_Of_RP_Pairs
C                     PCF RP LUN containing the number of ECS ODL
C                     objects implemented as successive name/value
C                     pairs in the RP section of the PCF.
C
C
C !OUTPUT PARAMETERS:
C  integer LUN_Of_LclVrsnID
C                     PCF RP LUN of LOCALVERSIONID.
C
C  integer rtn_flag   Subroutine return status: SUCCEED (0) or FAIL (-1).
C
C
C !REVISION HISTORY:
C 2001/03/30  rhucek
C Renamed subprogram "Get_LUN_Of_LocalVrsnID_Atmos" to "Get_LUN_Of_LclVrsnID_Atmos" 
C
C  11/27/00 rhucek: 
C  Corrected counting bug in loop index that was terminating Do While 
C  loop prematurely.  Loop index now starts at 1, is incremented by 2 on
C  each loop pass and is constrained to be < 2*NUM_Of_RP_Pairs.   
C
C
C !TEAM-UNIQUE HEADER:
C
C  This software is developed by the MODIS Science Data Support Team
C  for the National Aeronautics and Space Administration,
C  Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C
C !DESIGN NOTES:
C
C  ECS metadata attributes known prior to PGE startup are assumed to be 
C  implemented in the PCF as a series of ODL name/value pairs on 
C  successive RP LUNs.  At the top of the list is the number of 
C  name/value pairs in the list, an integer that is read and used to 
C  control looping over the set of LUNS.  Only RP names are actually 
C  read (on alternate LUNs), but each name is tested against the named 
C  constant MCORE_LOCALVERSIONID.  When a match is found, looping 
C  ends and the next LUN is assumed to contain the LOCALVERSIONID. 
C  
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        FAIL                      (Atmos_ECSMET.inc)
C        MCORE_LOCALVERSIONID      (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Get_LUN_Of_LclVrsnID_Atmos')

C Function argument declarations
      integer LUN_Of_NUM_Of_RP_Pairs
      integer LUN_Of_LclVrsnID
      integer rtn_flag


C other variable declarations

      character*25  msg25a
      character*512 msgbuf, StringValue_RP

      integer i, rtn, rtn_loc, rtn_status, NUM_Of_RP_Pairs
      integer fbyte_a, lbyte_a

      integer String_Loc

      logical error_flag, loopexit


C------------------------
C Initialization
C------------------------

      error_flag       = .FALSE.
      LUN_Of_LclVrsnID =  -1
      rtn_flag         =  FAIL


C-------------------------------------------------------------------------------
C Get LUN of RP LOCALVERSIONID
C-------------------------------------------------------------------------------

c-----get the number of name/value pairs in PCF RP ODL list

      Call Get_RP_Int( LUN_Of_NUM_Of_RP_Pairs,
     1                 NUM_Of_RP_Pairs,
     2                 rtn_status )

      If (rtn_status .ne. SUCCEED) Then
         error_flag  = .TRUE.

         write(msg25a,'(I25)') LUN_Of_NUM_Of_RP_Pairs
         rtn = string_loc(msg25a,fbyte_a,lbyte_a)

         msgbuf =
     1   'Get_RP_Int returned error reading RP on input argument '
     2   // '"LUN_Of_NUM_Of_RP_Pairs" (=' // msg25a(fbyte_a:lbyte_a) // ').'
     3   // char(10) // 'Operator Action:  Refer to prior low level error LogStatus '
     4   // char(10) // 'messages generated by routine Get_RP_Int.'

         Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c-----Find LUN_Of_LclVrsnID

      Else
         i = 1
         loopexit = .FALSE.

         Do While(i .lt. 2*NUM_Of_RP_Pairs .and. .NOT. loopexit)

            Call Get_RP_String( LUN_Of_NUM_Of_RP_Pairs+i,
     1                          StringValue_RP,
     2                          rtn_status )


c-----------error - get out of loop!
            If (rtn_status .eq. FAIL) Then
               error_flag = .TRUE.
               loopexit = .TRUE.

               write(msg25a,'(I25)') LUN_Of_NUM_Of_RP_Pairs+i
               rtn = string_loc(msg25a,fbyte_a,lbyte_a)

               msgbuf =
     1         'Get_RP_String returned error reading RP on LUN ' // msg25a(fbyte_a:lbyte_a)
     2         // char(10) // 'Operator Action:  Refer to prior low level error LogStatus '
     3         // char(10) // 'messages generated by routine Get_RP_String.'

               Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c-----------success reading RP
            Else

               If (StringValue_RP .eq. MCORE_LOCALVERSIONID) Then
                 LUN_Of_LclVrsnID = LUN_Of_NUM_Of_RP_Pairs + i + 1
                 loopexit = .TRUE.
               End If

            End If   ! check on Get_RP_String

            i=i+2
         End Do   ! do over number of RP name/value pairs

      End If   ! check on Get_RP_Int


      If (.not.error_flag .and. LUN_Of_LclVrsnID.ne.-1) rtn_flag = SUCCEED

      Return

      End  !Get_LUN_Of_LclVrsnID_Atmos



c---------------------------------------------------------------------------------
      Subroutine Get_RP_Int(LUN_Of_RP,IntValue_RP,rtn_flag)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Read and intepret a PCF Runtime Parameter (RP) string
C                as an integer value.
C
C
C !INPUT PARAMETERS:
C  Integer LUN_Of_RP     PCF LUN of RP value
C
C
C !OUTPUT PARAMETERS:
C  Integer IntValue_RP   Value of RP interpreted as an integer
C
C  Integer rtn_flag      Subroutine return flag  0 = SUCCESS
C                                               -1 = FAIL
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
C  Externals:
C
C    Functions and Subroutines:
C        Get_RP_String             (atmos shared code)
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        String_Loc                (atmos shared code)
C
C    Named Constant:
C        FAIL                      (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSd_PC_VALUE_LENGTH_MAX  (PGS_PC.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Get_RP_Int')

c-----Declaration of function arguments
      integer LUN_Of_RP,
     1        IntValue_RP,
     2        rtn_flag

c-----Declaration of local variables
      character*25   FORM, msg25, msg25a
      character*1024 msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) StringValue_RP

      integer string_loc

      integer fbyte,fbytea,
     1        lbyte,lbytea,
     2        IOS,lcl_rtn_flag,rtn

      logical error_flag


C------------------------
C Initialization
C------------------------

      rtn_flag     = FAIL
      lcl_rtn_flag = FAIL
      error_flag   = .FALSE.

      write(msg25a,'(I25)') LUN_Of_RP
      rtn = string_loc(msg25a,fbytea,lbytea)


C--------------------------------------------------------------
C Retrieve PCF RP string, then convert to integer value
C--------------------------------------------------------------

      Call Get_RP_String(LUN_Of_RP, StringValue_RP, lcl_rtn_flag)


c-----failed to successfully read RP string
      If (lcl_rtn_flag .eq. FAIL) Then
         error_flag = .TRUE.

         msgbuf =
     1   'Get_RP_String unable to read RP string on LUN ' // msg25a(fbytea:lbytea)
     2   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error'
     3   // char(10) // 'messages produced by routine Get_RP_String.'

         Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c-----construct integer FORMAT code to convert RP string to integer value
      Else
         write(msg25,'(I25)') PGSd_PC_VALUE_LENGTH_MAX
         rtn  = string_loc(msg25,fbyte,lbyte)
         FORM = '(I' // msg25(fbyte:lbyte) // ')'

         read(StringValue_RP,FORM,IOSTAT=IOS) IntValue_RP

         rtn = string_loc(StringValue_RP,fbyte,lbyte)


c--------test result of internal read for data conversion error
         If (IOS .ne. 0) then
            error_flag = .TRUE.

            msgbuf =
     1      'FORTRAN internal read error converting RP string '
     2      // char(10) // '(' // StringValue_RP(fbyte:lbyte) // ') to integer value.'
     3      // char(10) // 'Operator Action:  Notify SDST.'

            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If !-----data conversion error check

      End If !-----Get_RP_String


      If (.not.error_flag) rtn_flag = SUCCEED

      Return

      End  !Get_RP_Int



c---------------------------------------------------------------------------------
      subroutine Get_RP_String(LUN_RP, StringValue_RP, rtn_flag)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Read the string value associated with a PCF RP.
C
C
C !INPUT PARAMETERS:
C  Integer LUN_Of_RP     PCF LUN of RP value
C
C !OUTPUT PARAMETERS:
C  character*(*) StringValue_RP
C                        RP string value
C
C  Integer rtn_flag      Subroutine return flag (see Design Notes)
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
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_PC_GetConfigData      (libPGSTK.a)
C        String_Loc                (atmos shared code)
C
C    Named Constant:
C        BLANK                     (Atmos_ECSMET.inc)
C        FAIL                      (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Get_RP_String')


c-----Declaration of function arguments
      integer       LUN_RP, rtn_flag
      character*(*) StringValue_RP


c-----Declaration of local variables
      character*25   msg25
      character*1024 msgbuf

      integer pgs_pc_getconfigdata,string_loc

      integer fbyte,lbyte,rtn,rtn2

      logical error_flag


C------------------------
C Initialization
C------------------------

      rtn_flag = FAIL
      error_flag = .FALSE.

      write(msg25,'(I25)') LUN_RP
      rtn2 = String_Loc(msg25,fbyte,lbyte)


C----------------------------------------------------------------------
C Retrieve PCF RP string value
C----------------------------------------------------------------------

      rtn = pgs_pc_getconfigdata(LUN_RP, StringValue_RP)

c-----set error flag if unable to retrieve RP metadata value
      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf =
     1   'pgs_pc_getconfigdata unable to read RP string on LUN ' // msg25(fbyte:lbyte)
     2   // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     3   // char(10) // 'entry.  If LUN is nonexistent or value is blank, stage correct '
     4   // char(10) // 'PCF and rerun PGE.  Otherwise, notify SDST.'

         Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      Else If (StringValue_RP .eq. BLANK) Then
         error_flag = .true.

         msgbuf =
     1   'RP string on LUN ' // msg25(fbyte:lbyte) // ' is blank.  Cannot convert '
     2   // 'to valid integer value.'
     3   // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter '
     4   // char(10) // 'value.  If blank, stage correct PCF and rerun PGE.'
     5   // char(10) // 'Otherwise, notify SDST.'

         Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      End If !---pgs_pc_getconfigdata check



      If (.not.error_flag) rtn_flag = SUCCEED

      Return

      End  !Get_RP_String



c---------------------------------------------------------------------------------
      Subroutine Query_HDF_Attr( File_LUN,
     1                           File_VRSN,
     2                           HDF_FileAttrName,
     3                           rtn_status,
     4                           Attr_Status )

      Implicit None

c-----Insert Include files
      Include 'hdf.inc'
      Include 'PGS_MODIS_39500.f'
      Include 'PGS_PC.f'
      Include 'PGS_SMF.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Query an HDF file for the presence of global (file)
C                attribute identified by name.
C
C
C !INPUT PARAMETERS:
C
C    integer File_LUN:       PCF LUN number of file to be queried
C
C    integer File_VRSN:      PCF version number of file to be queried
C
C    character*(*) HDF_FileAttrName:
C                            Name of HDF file attribute search target
C
C
C !OUTPUT PARAMETERS:
C
C    integer rtn_status:     Subroutine return status (-1=fail; 0=success)
C
C    logical Attr_Status:    Attribute status (.true.=exists; .false.=does
C                            not exist)
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
c    HDF portions developed at the National Center for Supercomputing
c    Applications at the University of Illinois at Urbana-Champaign.
c
C !REFERENCES AND CREDITS:
C
C    Written by         Richard Hucek
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C !DESIGN NOTES:
C
C  Externals:
C
C    Functions and Subroutines:
C        modis_smf_setdynamicmsg   (atmos shared code)
C        pgs_pc_getreference       (sdptk library)
C        sffattr                   (hdf library)
C        sfend                     (hdf library)
C        sfstart                   (hdf library)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        DFACC_RDONLY            (hdf.inc)
C        FAIL                    (hdf.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        PGSd_PC_FILE_PATH_MAX   (PGS_PC.f)
C        PGS_S_SUCCESS           (PGS_SMF.f)
C        SUCCEED                 (hdf.inc)
C
C !END
C----------------------------------------------------------------------


c-----Parameter declarations
      Character*(*)  FUNCNAME
      Parameter     (FUNCNAME = 'Query_HDF_Attr')

      Character*(*)  BLANK
      Parameter     (BLANK = ' ')

      Logical        EXISTS
      Parameter     (EXISTS = .TRUE.)

      Logical        DOES_NOT_EXIST
      Parameter     (DOES_NOT_EXIST = .FALSE.)


c-----subroutine argument declarations
      Character*(*)  HDF_FileAttrName
      Integer        File_LUN, File_VRSN, rtn_status
      Logical        Attr_Status


c-----Local variable declarations
      Character*25   msg25_LUN, msg25_VRSN
      Character*1024 msgbuf
      Character*(PGSd_PC_FILE_PATH_MAX) HDF_FileName

      Integer        fbyte1,fbyte2,fbyte3,fbyte4,
     1               lbyte1,lbyte2,lbyte3,lbyte4,
     2               rtn,rtn2,local_vrsn
      Integer        hdf_fileID

c-----Functions declarations
      Integer        pgs_pc_getreference     ! SDPTK
      Integer        sfend,sffattr,sfstart   ! HDF
      Integer        string_loc              ! Other


C-----------------------------------------------------------------------
C Initialization
C-----------------------------------------------------------------------

       rtn_status  = SUCCEED
       Attr_Status = DOES_NOT_EXIST


C-----------------------------------------------------------------------
C Set status message variables
C-----------------------------------------------------------------------

c-----Set component status messages
      write(msg25_LUN, '(I25)') File_LUN
      write(msg25_VRSN,'(I25)') File_VRSN

      rtn = String_Loc(msg25_LUN,fbyte1,lbyte1)
      rtn = String_Loc(msg25_VRSN,fbyte2,lbyte2)
      rtn = String_Loc(HDF_FileAttrName,fbyte3,lbyte3)

C-----------------------------------------------------------------------
C Perform input argument checks and return if not valid
C-----------------------------------------------------------------------

c-----check for valid PCF LUN
      If (File_LUN .LE. 0) Then
         rtn_status = FAIL

         msgbuf =
     1   'PCF LUN number ' // msg25_LUN(fbyte1:lbyte1) // ' is out of bounds.'
     2   // char(10) // FUNCNAME // ' unable to identifiy input file'
     3   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     4   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


c-----check for valid PCF version number
      If (File_VRSN .LT. 1) Then
         rtn_status = FAIL

         msgbuf =
     1   'PCF version number (' // msg25_VRSN(fbyte2:lbyte2) // ') is out of bounds '
     2   // char(10) // 'on LUN number ' // msg25_LUN(fbyte1:lbyte1)
     3   // char(10) // FUNCNAME // ' unable to identify input file.'
     4   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     5   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


c-----check for blank file name
      If (HDF_FileAttrName .EQ. BLANK) Then
         rtn_status = FAIL

         msgbuf =
     1   'HDF file attribute name is blank.'
     2   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     3   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'

        Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


c-----------------------------------------------------------------------
c Return if invalid inputs found
c-----------------------------------------------------------------------

      If (rtn_status .EQ. FAIL) Return


c-----------------------------------------------------------------------
c Check if attribute is actually present in file as global attribute
c-----------------------------------------------------------------------

c-----get name of product file
      local_vrsn = File_VRSN  !WARNING, pgs_pc_getreference resets local_vrsn
      rtn  = pgs_pc_getreference(File_LUN,local_vrsn,HDF_FileName)
      rtn2 = string_loc(HDF_FileName,fbyte4,lbyte4)

      If (rtn .NE. PGS_S_SUCCESS) Then  ! unable to retrieve file name from PCF
          rtn_status = FAIL

          msgbuf =
     1    'pgs_pc_getreference() unable to retrieve filename '
     2    // char(10) // HDF_FileName(fbyte4:lbyte4)
     3    // char(10) // 'on PCF LUN ' // msg25_LUN(fbyte1:lbyte1) // ' and VRSN '
     4    // msg25_VRSN(fbyte2:lbyte2)
     5    // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     6    // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     7    // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

          Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      Else
         hdf_fileID = sfstart(HDF_FileName, DFACC_RDONLY)

c--------unable to open hdf product file
         If (hdf_fileID .EQ. FAIL) Then
            rtn_status = FAIL

            msgbuf =
     1      'HDF function sfstart unable to open file '
     2      // char(10) // HDF_FileName(fbyte4:lbyte4)
     3      // char(10) // 'on LUN ' // msg25_LUN(fbyte1:lbyte1) // ' and VRSN '
     4      // msg25_VRSN(fbyte2:lbyte2)
     5      // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     6      // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     7      // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'


            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------successfully opened HDF product file
         Else

c-----------NOTE: FAIL from sffattr is interpreted not as a S/W error, but
c-----------rather that the requested attribute is not in file.
            rtn = sffattr(hdf_fileID, HDF_FileAttrName)

c-----------attribute exists - change default status to EXISTS
            If (rtn .NE. FAIL) Attr_Status = EXISTS

            rtn = sfend(hdf_fileID)

c-----------unable to close hdf file
            If (rtn .EQ. FAIL) Then
               rtn_status = FAIL

               msgbuf = 'HDF function sfend unable to close HDF product file '
     1         // char(10) // HDF_FileName(fbyte4:lbyte4)
     2         // char(10) // 'on LUN ' // msg25_LUN(fbyte1:lbyte1) // ' and VRSN '
     3         // msg25_VRSN(fbyte2:lbyte2)

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            End If

         End If   ! open HDF file

      End If   ! read HDF file name from PCF

      Return

      End  !Query_HDF_Attr



c-----------------------------------------------------------------------
      Subroutine Retrieve_LclGranID_Set( NUM_LclInputPntr,
     1                                   LUN_LclInputPntr,
     2                                   VRSN_LclInputPntr,
     3                                   NUM_LocalGranID,
     4                                   LocalGranID,
     5                                   rtn_flag )

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC_9.f'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Retrieve a list of input file names from the PCF, each 
C                identified by LUN and version number.  MODIS L1 product 
C                file names are read from the LOCALGRANULEID attribute 
C                stored in the ECS Inventory metadata record of every 
C                MODIS standard data product.  All other file names, 
C                including MODIS and ancillary data products, are read
C                from field 2 of the PCF line entry matching the file's 
C                LUN and version number.  (See comments in code) 
C
C
C !INPUT PARAMETERS:
C
C  integer NUM_LclInputPntr
C                Variable containing the number of names to be retrieved 
C
C  integer LUN_LclInputPntr(*)
C                Array of input product LUNs 
C
C  integer VRSN_LclInputPntr(*)
C                Array of input product version numbers
C
C                A one-to-one correspondence between the elements of
C                the arrays LUN_LclInputPntr and
C                VRSN_LclInputPntr is assumed.
C
C !OUTPUT PARAMETERS:
C
C  character*(*) LocalGranID(*)
C                Array of file names read from the PCF or HDF file (for
C                L1 inputs) 
C
C  integer       NUM_LocalGranID 
C                Number of file names successfully retrieved from array
C                of input file references (LUN_LclInputPntr and
C                VRSN_LclInputPntr)
C
C  integer rtn_flag
C                Procedure return flag (SUCCEED=0, FAIL=-1)
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
C        MODIS_SMF_SETDYNAMICMSG    (atmos shared code)
C        PGS_MET_GetPCAttr_s        (libPGSTK.a)
C        PGS_PC_GetReference        (libPGSTK.a)
C        String_Loc                 (atmos shared code)
C
C    Named Constant:
C        BLANK                      (Atmos_ECSMET.inc)
C        FAIL                       (Atmos_ECSMET.inc)
C        LUN_IS_420000              (Atmos_ECSMET.inc)
C        LUN_IS_422500              (Atmos_ECSMET.inc)
C        LUN_IS_600000              (Atmos_ECSMET.inc)
C        LUN_IS_700000              (Atmos_ECSMET.inc)
C        LUN_IS_700001              (Atmos_ECSMET.inc)
C        LUN_IS_700002              (Atmos_ECSMET.inc)
C        MECS_CORE                  (mapi.inc)
C        MCORE_LOCALGRANULEID       (mapi.inc)
C        MODIS_E_GENERIC            (PGS_MODIS_39500.f)
C        MODIS_W_GENERIC            (PGS_MODIS_39500.f)
C        NO_SET_UR                  (Atmos_ECSMET.inc)
C        PGSPC_W_NO_REFERENCE_FOUND (PGS_PC_9.f)
C        PGSd_PC_VALUE_LENGTH_MAX   (PGS_PC.f)
C        PGS_S_SUCCESS              (PGS_SMF.f)
C        SUCCEED                    (Atmos_ECSMET.inc)
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*)   FUNCNAME
      parameter      (FUNCNAME = 'Retrieve_LclGranID_Set')

c-----Declaration of function arguments
      character*(*) LocalGranID(*)

      integer NUM_LclInputPntr
      integer LUN_LclInputPntr(*)
      integer VRSN_LclInputPntr(*)
      integer NUM_LocalGranID,rtn_flag

c-----Declaration of local variables
      character*25  msg25a, msg25b, msg25c
      character*60  AttrN
      character*512 msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) AttrV
      character*(PGSd_PC_FILE_PATH_MAX) FilePathName

      integer pgs_met_getpcattr_s,
     1        pgs_pc_getreference,
     2        string_loc

      integer fbyte,fbytea,fbyteb,fbytec,
     1        lbyte,lbytea,lbyteb,lbytec,
     2        icounter,localversion,LUN,VRSN,rtn


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

c-----loop over the number of MODIS file references
      Do 100 icounter = 1, NUM_LclInputPntr
         VRSN     = VRSN_LclInputPntr(icounter)
         LUN      = LUN_LclInputPntr(icounter)

c--------set message buffers for error reporting
         write(msg25a,'(i25)') LUN
         write(msg25b,'(i25)') VRSN
         write(msg25c,'(i25)') icounter

         rtn = String_Loc(msg25a,fbytea,lbytea)
         rtn = String_Loc(msg25b,fbyteb,lbyteb)
         rtn = String_Loc(msg25c,fbytec,lbytec)

c--------ignore InputPointers with LUN values equal to NO_SET_UR.
         If (LUN .EQ. NO_SET_UR) Then
            if (icounter.ne.6) then
               continue
            else
c --------- threshold file position
               NUM_LocalGranID = NUM_LocalGranID + 1
            end if

c--------check file version number
         Else If (VRSN .LT. 1)  Then
            rtn_flag = FAIL

            msgbuf =
     1      'MODIS file PCF VRSN number (' // msg25b(fbyteb:lbyteb) // ') on LUN '
     2      // msg25a(fbytea:lbytea) // ' is less than one.'
     3      // char(10) // 'Unable to retrieve LOCALINPUTGRANULEID array element '
     4      // msg25c(fbytec:lbytec)  // '. '
     5      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c--------check LUN constraints
         Else If ( (LUN .LE. 0) .OR. (LUN.GE.10000 .AND. LUN.LE.10999) ) Then
            rtn_flag = FAIL

            msgbuf =
     1      'MODIS input file LUN number (' // msg25a(fbytea:lbytea) // ') is out of bounds. '
     2      // char(10) // 'LOCALINPUTGRANULEID array element ' // msg25c(fbytec:lbytec)
     3      // ' not retrieved. '
     4      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c--------retrieve the LOCALGRANULEID of MODIS L1 products 
         Else If ( LUN .EQ. LUN_Is_420000 .OR. 
     1             LUN .EQ. LUN_Is_422500 .OR.
     2             LUN .EQ. LUN_Is_600000 .OR.
     3             LUN .EQ. LUN_Is_700000 .OR.
     4             LUN .EQ. LUN_Is_700001 .OR.
     5             LUN .EQ. LUN_Is_700002 ) Then

c           During PGE processing, output product files are represented in GDAAC PCF 
c           files with names very similar to their MODIS format, but with the 
c           distinction that production datetime field is set to the processing start,
c           not the processing end time used by MODIS.  
c
c           This is a concern for PGEs running more than one process where later 
c           processes accesss the output of former.  In this case, e.g. PGE03, 
c           reading the MOD35_L2 file name (as does MOD_PR07) from the GDAAC PCF 
c           would give different results than reading the MODIS file name directly 
c           from the LocalGranuleId stored in the MOD_PR35 HDF product metadata.  
c
c           To avoid ambiguity, the file names of MODIS products that could be run
c           at the GDAAC will always be read directly from the the HDF file 
c           LocalGranuleId metadata and not from the PCF.

            rtn = pgs_met_getpcattr_s(LUN,VRSN,MECS_CORE,AttrN,AttrV)

            If (rtn .ne. PGS_S_SUCCESS) Then
               rtn_flag = FAIL

               msgbuf =
     1         'pgs_met_getpcattr_s unable to retrieve LOCALGRANULEID from '
     2         // char(10) // 'MODIS product on PCF LUN ' // msg25a(fbytea:lbytea)
     3         // 'and  VRSN ' // msg25b(fbyteb:lbyteb) // '. '
     4         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5         // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     6         // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else If (AttrV .EQ. BLANK) Then

               msgbuf =
     1         'LOCALINPUTGRANULEID of MODIS file on PCF LUN ' // msg25a(fbytea:lbytea)
     2         // ' and VRSN ' // msg25b(fbyteb:lbyteb)
     3         // char(10) // 'is blank.  Blank LOCALINPUTGRANULEIDs not set.'
     4         // char(10) // 'Notify SDST. '

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)

            Else
               NUM_LocalGranID = NUM_LocalGranID + 1
               LocalGranID(NUM_LocalGranID) = AttrV 

            End If


c--------All other products 
         Else
            localversion = VRSN
            rtn = pgs_pc_getreference(LUN, localversion, FilePathName)

c-----------LUN not in PCF 
            If (rtn .eq. PGSPC_W_NO_REFERENCE_FOUND) Then
               msgbuf = 'No LocalInputGranuleID set to ECS metadata record for PCF# ' 
     1         // msg25a(fbytea:lbytea) // '.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)

c-----------error retrieving file name
            Else If (rtn .ne. PGS_S_SUCCESS) Then
               rtn_flag = FAIL

               msgbuf =
     1         'pgs_pc_getreference unable to retrieve file path name '
     2         // char(10) // 'from PCF on LUN ' // msg25a(fbytea:lbytea) // ' '
     3         // 'and VRSN ' // msg25b(fbyteb:lbyteb) // '. '
     4         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     5         // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     6         // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else 
               Call Cut_Name(FilePathName,fbyte,lbyte)

               If (fbyte .EQ. lbyte) Then 
                  msgbuf =
     1            'File name on PCF LUN ' // msg25a(fbytea:lbytea)
     2            // ' and VRSN ' // msg25b(fbyteb:lbyteb)
     3            // char(10) // 'is blank.  Blank LOCALINPUTGRANULEIDs not set.'
     4            // char(10) // 'Notify SDST. '

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)

               Else
                  NUM_LocalGranID = NUM_LocalGranID + 1
                  LocalGranID(NUM_LocalGranID) = FilePathName(fbyte:lbyte) 

               End If   ! test blank file name

            End If   ! test getreference return 

         End If  ! check on LUN and VRSN file identifier conditions

100   Continue

      Return

      End  !Retrieve_LclGranID_Set



c---------------------------------------------------------------------------------
      Subroutine Retrieve_UR_Set( NUM_InputPntr,
     1                            LUN_InputPntr,
     2                            VRSN_InputPntr,
     3                            NMBR_UR,
     4                            UR,
     5                            rtn_flag )

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'PGS_MODIS_39500.f'
      include 'PGS_PC_9.f'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Retrieve ECS Universal References (UR) for a specified
C                file set (identified by PCF LUN and version numbers)
C
C
C !INPUT PARAMETERS:
C  integer  NUM_InputPntr  Number of valid file references in LUN_InputPntr
C                          and VRSN_InputPntr array arguments.  Ancillary
C                          data files, look up tables, and MODIS product
C                          files are included.  System files, such as the
C                          MCF, are not included.
C
C  integer  LUN_InputPntr(*)
C                          Array containing a list of PCF file LUNs (See
C                          NUM_InputPntr above).
C
C                          A one-to-one correspondence between the elements
C                          of the arrays LUN_InputPntr and VRSN_InputPntr
C                          is assumed.
C
C  integer  VRSN_InputPntr(*)
C                          Array containing a list of PCF file version
C                          numbers  (See NUM_InputPntr and LUN_InputPntr
C                          above).
C
C
C !OUTPUT PARAMETERS:
C  character*(*) UR(*)     Array containing UR strings retrieved from PCF.
C                          Blanks (' '), but not missing URs, are returned.
C
C  integer NMBR_UR         Variable containing number of URs successfully
C                          retrieved from PCF.  This includes blank (' ')
C                          UR values.
C
C  integer rtn_flag        Procedure return status: SUCCEED (0) or FAIL (-1)
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
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_PC_GetUniversalRef    (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        BLANK                      (Atmos_ECSMET.inc)
C        FAIL                       (Atmos_ECSMET.inc)
C        IGNORE_ENTRY               (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC            (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS              (PGS_SMF.f)
C        PGSd_PC_UREF_LENGTH_MAX    (PGS_PC.f)
C        PGSPC_W_NO_REFERENCE_FOUND (PGS_PC_9.f)
C        SUCCEED                    (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Retrieve_UR_Set')

C passed-in argument declarations
      integer NUM_InputPntr, LUN_InputPntr(*), VRSN_InputPntr(*), NMBR_UR
      character*(*) UR(*)
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
      NMBR_UR  = 0


c-----loop over number of file references
      Do 100 icounter = 1, NUM_InputPntr
         Version_No = VRSN_InputPntr(icounter)
         LUN = LUN_InputPntr(icounter)

c-----set message buffers for error reporting
         write(msg25b,'(i25)') LUN
         write(msg25c,'(i25)') Version_No
         write(msg25d,'(i25)') icounter

         rtn = String_Loc(msg25b,fbyteb,lbyteb)
         rtn = String_Loc(msg25c,fbytec,lbytec)
         rtn = String_Loc(msg25d,fbyted,lbyted)


c--------ignore InputPointers with LUN entries equal to NO_SET_UR
         If (LUN .EQ. NO_SET_UR) Then
            continue

c--------set error flag if file version number is out of bounds
         Else If (Version_No .LT. 1) Then
            rtn_flag = FAIL

            msgbuf =
     1      'PCF file version number ' // msg25c(fbytec:lbytec) // ' on LUN '
     2      // msg25b(fbyteb:lbyteb) // ' is out of bounds.'
     3      // char(10) // 'UR of element ' // msg25d(fbyted:lbyted)
     4      // ' of input pointer array cannot be retrieved.'
     5      // char(10) // 'Operator Action:  Notify SDST.'


            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------set error flag if file LUN is out of bounds
         Else If (LUN.LT.1 .OR. (LUN.GE.10000 .AND. LUN.LE.10999)) Then
            rtn_flag = FAIL

            msgbuf = 'PCF file LUN number ' // msg25b(fbyteb:lbyteb)
     1      // ' is out of bounds.'
     2      // char(10) // 'UR of element ' // msg25d(fbyted:lbyted)
     3      // ' of input pointer array cannot be retrieved.'
     4      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------retrieve PCF UR
         Else
            rtn = PGS_PC_GetUniversalRef(LUN, Version_No, work_buf)

            If (rtn .eq. PGSPC_W_NO_REFERENCE_FOUND) Then
               msgbuf = 'No InputPointer set to ECS metadata record for PCF# ' 
     1         // msg25b(fbyteb:lbyteb) // '.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)

            Else If (rtn .ne. PGS_S_SUCCESS) Then
               rtn_flag = FAIL

               msgbuf =
     1         'PGS_PC_GetUniversalRef unable to retrieve UR on PCF LUN '
     2         // msg25b(fbyteb:lbyteb) // ' and VRSN ' // msg25c(fbytec:lbytec) // '. '
     3         // char(10) // 'Operator Action:  Check PCF Universal Reference Identifier (UR).'
     4         // char(10) // 'If invalid, stage correct PCF and rerun PGE.  Otherwise, notify SDST.'


               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----------if UR successfully read, insert it into output array
            Else
               NMBR_UR     = NMBR_UR + 1
               UR(NMBR_UR) = work_buf

               If (work_buf .EQ. BLANK) Then

                  msgbuf =
     1            'UR on LUN ' // msg25b(fbyteb:lbyteb) // ' and VRSN ' // msg25c(fbytec:lbytec)
     2            // ' is blank. '
     3            // char(10) // 'Operator Action:  Check PCF Universal Reference identifier. '
     4            // char(10) // 'If blank, stage correct PCF and rerun PGE.  Otherwise, '
     5            // char(10) // 'notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,FUNCNAME)
               End If   ! check on BLANK UR

            End If   ! check on PGS_PC_GetUniversalRef rtn

         End If !(Version_No .LT. 1)

100   Continue

      Return

      End  !Retrieve_UR_Set



c---------------------------------------------------------------------------------
      Integer Function Gen_Modis_FileName_V2(ESDT_Name,Local_VRSN_ID,
     1                                       FileNameSuffix, fname)

      implicit none

      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------------------------
C !F77
C
C !Description:  This function produces a file name in the MODIS V2 
C                naming convention.  It accepts an arbitrary 
C                3-character string (but no blanks) and appends it 
C                along with a dot "." as the file suffix. 
C
C
C !Input Parameters:
C       character*(*) ESDT_Name         ESDT shortname of the product 
C
C       character*(*) Local_VRSN_ID     MODIS 3-digit string 
C                                       representation of ECS metadata 
C                                       field LOCALVERSIONID 
C
C       character*(*) FileNameSuffix    string variable containing file 
C                                       suffix (e.g, 'dat', 'exe', 'hdf').
C                                       Suffix length up to 3 characters,
C                                       starting at element 1 and terminated
C                                       by first blank character. 
C                                        
C
C
C !Output Parameters: 
C
C       character*(*) fname             MODIS file name (LOCALGRANULEID)
C
C !Revision History:
c Revision 1.5  1998/10/28  13:08:49  rhucek
c Changed code line from "istart = istart + GMT_Delta*60" to
c                        "istart = istart - GMT_Delta*60".\
c
c Revision 1.4  1998/07/20  17:18:26  rhucek
c Revised code to extract granule start time from RangeBeginningDate/Time
c on MOD03 (Geolocation Product) rather than from CollectionStartTime in PCF.
c
c Revision 1.3  1998/06/18  15:53:23  fhliang
c Modified error messages and description of parameter MAX_MODIS_FNAME_LEN.
c
c Revision 1.2  1998/04/24  14:54:40  lma
c Updated SDPTK error messages with "Operator Action:" text.
c
c Revision 1.1  1997/11/13  20:17:48  lma
c Initial revision
c
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C !References and Credits:
C
C    Written by:
C
C    Liqun Ma         11/05/97
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    lma@ltpmail.gsfc.nasa.gov
C    rhucek@ltpmail.gsfc.nasa.gov
C
C!DESIGN NOTES:
C
C    Externals:
C
C       Functions and Subroutines:
C         pgs_td_asciitime_atob         (libPGSTK.a)
C         pgs_td_utctotai               (libPGSTK.a)
C         pgs_td_taitoutc               (libPGSTK.a)
C         pgs_met_getpcattr_s           (libPGSTK.a)
C         MODIS_SMF_SETDYNAMICMSG       (L2 atmos share code library) 
C         strlen                        (L2 atmos share code library) 
C         string_loc                    (L2 atmos share code library) 
C         Get_date_time                 (L2 atmos share code library)
C
C
C       Named Constant:
C         MCORE_RANGE_BEG_DATE          mapi.inc
C         MCORE_RANGE_BEG_TIME          mapi.inc
C         MECS_CORE                     mapi.inc
C         PGS_S_SUCCESS                 PGS_SMF.f
C
C    Internals:
C
C       Functions and Subroutines:
C
C
C !END
C-----------------------------------------------------------------------------------------


C PARAMETER declarations
      character*(*)  FUNCNAME
      parameter    ( FUNCNAME = 'Gen_Modis_FileName_V2' )

      character*1    BLANK
      parameter    ( BLANK = ' ' )

      integer        LRN_GEO
      parameter    ( LRN_GEO = 600000 )

      integer        VRSN_GEO
      parameter    ( VRSN_GEO = 1 )

      integer        FAIL,      SUCCEED
      parameter    ( FAIL = -1, SUCCEED = 0 )

      integer        MAX_ESDT_Name_LEN
      parameter    ( MAX_ESDT_Name_LEN = 8 )

      integer        MAX_MODIS_FNAME_LEN
      parameter    ( MAX_MODIS_FNAME_LEN = 44 )

C Declaration of function calling arguments
      character*(*) Local_VRSN_ID, ESDT_Name, FileNameSuffix, fname


C Declaration of local variables
      character date*(10),time*(10),work_buf*(28),collectn_b*(28),datetime_b*(27)
      character*3   version
      character*25  msg25,msg25a, msg25b, msg25c
      character*512 msgbuf
      character*128 AttrN, buf_char(2)

      integer buf_len,i,nchar,rtn,sl,string_len
      integer fbyte,fbytea,fbyteb,fbytec
      integer lbyte,lbytea,lbyteb,lbytec
      integer ESDT_Name_LEN,GMT_Delta

      double precision istart

      logical error_flag,VRSN_error_flag

C Function declarations
      integer pgs_td_asciitime_atob,pgs_td_utctotai,pgs_td_taitoutc,
     &        PGS_MET_GetPCAttr_s
      integer strlen,string_loc


*----------------------------------------------------------------------------*


      error_flag            = .FALSE.
      VRSN_error_flag       = .FALSE.
      Gen_Modis_FileName_V2 =  FAIL


C-----------------------------------------------------------------------------------------
C Perform initial input argument checks and return if not valid. 
C
C 1 - ESDT name must not be blank and may contain no more 
C     than 8 characters
C 2 - Local_VRSN_ID must be non-blank, and represent a positive integer  
C     from 1 to 999.
C 3 - String buffer fname must be at least 44 characters in size.
C-----------------------------------------------------------------------------------------

c.....Tests for blank ESDT_Name
      If (ESDT_Name .eq. BLANK) Then
         error_flag = .true.

         msgbuf = 'ESDT name is blank.'
     1   //  char(10) // 'Gen_Modis_FileName_V2 unable to construct MODIS file name ' 
     2   //  char(10) // 'without valid ESDT.'
     3   //  char(10)//'Operator Action: Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c.....ESDT name is non-blank; check name length
      Else 
         rtn = String_Loc(ESDT_Name,fbyte,lbyte)
         ESDT_Name_LEN = lbyte - fbyte + 1

         write(msg25a,'(I25)') ESDT_Name_LEN
         rtn = String_Loc(msg25a,fbytea,lbytea)
 
         If (ESDT_Name_LEN .gt. MAX_ESDT_Name_LEN) Then
            error_flag = .true.

            msgbuf = 'ESDT name contains ' // msg25a(fbytea:lbytea) // ' characters, ' 
     1      // char(10) // 'too long to comply with the MODIS file naming convention.'
     2      // char(10) // 'Operator Action: Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

         End If
      End If 


c.....Tests for blank Local_VRSN_ID
      If (Local_VRSN_ID .eq. BLANK) Then
         error_flag = .true.

         msgbuf = 'Local_VRSN_ID is blank.'
     1   // char(10) // 'Gen_Modis_FileName_V2 unable to construct MODIS file name ' 
     2   // char(10) // 'without valid Local version ID.'
     3   // char(10) // 'Operator Action: Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c.....Local_VRSN_ID is non-blank; check for valid range 
      Else 
         rtn = String_Loc(Local_VRSN_ID,fbyte,lbyte)
         string_len = lbyte - fbyte + 1
 
c........string length greater than 3 characters
         If (string_len .GT. 3) Then
            error_flag = .true.

            write(msg25a,'(I25)') string_len 
            rtn = String_Loc(msg25a,fbytea,lbytea)
 
            msgbuf = 'Local_VRSN_ID contains ' // msg25a(fbytea:lbytea)
     1      // ' characters, too long to ' 
     2      // char(10) // 'comply with 3-digit MODIS version number convention.'
     3      // char(10) // 'Operator Action: Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)
         Else 
            If (string_len .eq. 1) version = '00' // Local_VRSN_ID(fbyte:lbyte)    
            If (string_len .eq. 2) version = '0' // Local_VRSN_ID(fbyte:lbyte)    
            If (string_len .eq. 3) version = Local_VRSN_ID(fbyte:lbyte)

c...........check for non-digits in version number
            Do i=1,3
               If (version(i:i).lt.'0' .or. version(i:i).gt.'9') 
     1            VRSN_error_flag = .true.            
            End Do

            If (VRSN_error_flag) Then
               error_flag = .true.
         
               msgbuf = 'Local_VRSN_ID "' // version // '" does not '
     1         // 'represent a 3-digit integer.'
     2         // char(10) // 'Gen_Modis_FileName_V2 unable to construct MODIS file name.'
     3         // char(10) // 'Operator Action: Notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

             End If  ! check on number of valid digits 
         End If   ! check for 3 digits or less 
      End If   ! check for blank version id  



      buf_len = LEN(fname)

c.....Test for adequate fname buffer length
      If (buf_len .lt. MAX_MODIS_FNAME_LEN) Then
         error_flag = .true.

         write(msg25,'(I25)') buf_len
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'string buffer "fname" contains only ' // msg25(fbyte:lbyte) 
     1   // ' characters, ' 
     2   // char(10) // 'too small to contain MODIS file name.'
     3   // char(10) // 'Operator Action: Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


C-----------------------------------------------------------------------------------------
C Return if input arguments are invalid
C-----------------------------------------------------------------------------------------

      If (error_flag) Return


C-----------------------------------------------------------------------------------------
C Retrieve values of RangeBeginningDate and RangeBeginningTime ECS 
C metadata fields C from Geolocation Product.  Combine into CCSDS time 
C format A and convert time format A to B.
C-----------------------------------------------------------------------------------------

c.....setup for error reporting
      write(msg25b,'(I25)') LRN_GEO
      rtn = string_loc(msg25b,fbyteb,lbyteb)

      write(msg25c,'(I25)') VRSN_GEO
      rtn = string_loc(msg25c,fbytec,lbytec)

c.....Retrieve RangeBeginningDate
      AttrN = MCORE_RANGE_BEG_DATE
      rtn = PGS_MET_GetPCAttr_s(LRN_GEO,VRSN_GEO,MECS_CORE,AttrN,buf_char(1))

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.
         rtn = String_Loc(AttrN,fbyte,lbyte)

         msgbuf = 'PGS_MET_GetPCAttr_s unable to retrieve ECS attribute ' // AttrN(fbyte:lbyte) 
     1   // char(10) // 'from Geolocation product (LUN = ' // msg25b(fbyteb:lbyteb)
     2   // ', Version = ' // msg25c(fbytec:lbytec) // ').'
     3   // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     4   // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     5   // char(10) // 'fault is identified, stage correct PCF/input file and '
     6   // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

c.....Retrieve RangeBeginningTime
      Else  
         AttrN = MCORE_RANGE_BEG_TIME
         rtn = PGS_MET_GetPCAttr_s(LRN_GEO,VRSN_GEO,MECS_CORE,AttrN,buf_char(2))

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.
            rtn = String_Loc(AttrN,fbyte,lbyte)

            msgbuf = 'PGS_MET_GetPCAttr_s unable to retrieve ECS attribute ' // AttrN(fbyte:lbyte)
     1      // char(10) // 'from Geolocation product (LUN = ' // msg25b(fbyteb:lbyteb)
     2      // ', Version = ' // msg25c(fbytec:lbytec) // ').'
     3      // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     4      // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     5      // char(10) // 'fault is identified, stage correct PCF/input file and '
     6      // char(10) // 'rerun PGE.  Otherwise, notify SDST.'

            call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

c........Convert collection start time from CCSDS format A to format B.
         Else     
            rtn = string_loc(buf_char(1),fbyteb,lbyteb)
            rtn = string_loc(buf_char(2),fbytec,lbytec)

            work_buf = buf_char(1)(fbyteb:lbyteb) // 'T' // buf_char(2)(fbytec:lbytec)

            rtn=pgs_td_asciitime_atob(work_buf,collectn_b)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf = 'pgs_td_asciitime_atob unable to convert collection start time '
     1         // char(10) //'from CCSDS format A to format B.'
     2         // char(10) //'Operator Action: '
     3         // 'Check system and SDPTK configuration. If fault identified,'
     4         // char(10) //'correct and rerun PGE.  Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            End If   ! convert time format A to B
         End If   ! read RangeBeginningTime 
      End If   ! read RangeBeginningDate

C-----------------------------------------------------------------------
C Query system for current local time and transform it to current GMT 
C time in CCSDS B format.
C-----------------------------------------------------------------------

      If (.not. error_flag) Then

C ...... Get current local time 

         call Get_date_time(date,time,GMT_Delta)

C .......construct the current local time in CCSDS format A

         work_buf(1:4)=date(1:4)
         work_buf(5:5)='-'
         work_buf(6:7)=date(5:6)
         work_buf(8:8)='-'
         work_buf(9:10)=date(7:8)
         work_buf(11:11)='T'

         work_buf(12:13)=time(1:2)
         work_buf(14:14)=':'
         work_buf(15:16)=time(3:4)
         work_buf(17:17)=':'
         work_buf(18:23)=time(5:10)



C .......Convert current local time from CCSDS A to TAI format
         rtn=pgs_td_utctotai(work_buf,istart)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.

            msgbuf = 'pgs_td_utctotai unable to convert current'
     1      // char(10) // 'local time from CCSDS A to TAI format.'
     2      // char(10) // 'Operator Action:  Check system resources/environment, '
     3      // char(10) // 'PCF, and SDPTK configuration (including up-to-date Leap '
     4      // char(10) // 'Seconds and UT1 Pole files).  If a fault is identified, '
     5      // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

         Else

C ..........Convert current locat time to GMT time
            istart = istart - GMT_Delta*60

C ..........Convert current GMT time from tai to CCSDS format A
            rtn=pgs_td_taitoutc(istart,work_buf)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf = 'pgs_td_taitoutc unable to convert current GMT '
     1         // char(10) //'time from TAI to CCSDS A format.'
     2         // char(10) //'Operator Action: '
     3         // 'Check system and SDPTK configuration. If fault identified,'
     4         // char(10) //'correct and rerun PGE.  Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

            Else

C .............Convert current GMT time from CCSDS format A to format B
               rtn=pgs_td_asciitime_atob(work_buf,datetime_b)


               If (rtn.ne.PGS_S_SUCCESS) Then
                  error_flag = .true.

                  msgbuf = 'pgs_td_asciitime_atob unable to convert current '
     1            // char(10) // 'GMT time from CCSDS A to CCSDS B format.'
     2            // char(10) // 'Operator Action: '
     3            // 'Check system and SDPTK configuration. If fault identified,'
     4            // char(10) //'correct and rerun PGE.  Otherwise, notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

               End If   ! convert CCSDS from A to B format
            End If   ! Convert GMT from tai to CCSDS A
         End If   ! convert local time CCSDS A to tai
      End If   ! check error flag


C-----------------------------------------------------------------------------------------
C Construct text string containing file name in MODIS convention 
C-----------------------------------------------------------------------------------------

      If(.not.error_flag) Then

C .......Calculate the actual string length of short name.

         sl=strlen(ESDT_Name)

C........Identify consecutive non-blank characters from 1 to 3 in  
C........FileNameSuffix

         nchar=0
         if ( FileNameSuffix(1:1).ne.BLANK ) nchar=1 
         if ( FileNameSuffix(2:2).ne.BLANK .and. nchar.eq.1) nchar=2 
         if ( FileNameSuffix(3:3).ne.BLANK .and. nchar.eq.2) nchar=3 

C .......Construct the filename from EDST name, local version ID,
C .......collection start time format B and current GMT time format B.

         fname(1:sl)=ESDT_Name(1:sl)
         fname(sl+1:sl+1)='.'
         fname(sl+2:sl+2)='A'

         fname(sl+3:sl+6)=collectn_b(1:4)
         fname(sl+7:sl+9)=collectn_b(6:8)
         fname(sl+10:sl+10)='.'
         fname(sl+11:sl+12)=collectn_b(10:11)
         fname(sl+13:sl+14)=collectn_b(13:14)

         fname(sl+15:sl+15)='.'

         fname(sl+16:sl+18)=version(1:3)

         fname(sl+19:sl+19)='.'

         fname(sl+20:sl+23)=datetime_b(1:4)
         fname(sl+24:sl+26)=datetime_b(6:8)
         fname(sl+27:sl+28)=datetime_b(10:11)
         fname(sl+29:sl+30)=datetime_b(13:14)
         fname(sl+31:sl+32)=datetime_b(16:17)

C........Add suffix if 

         If (nchar.ge.1 .and. nchar.le.3) Then
            fname(sl+33:sl+33)='.'

            Do i = 1,nchar
               fname(sl+33+i:sl+33+i)=FileNameSuffix(i:i)
            End Do

         End If

      End If

      If (.not. error_flag)  Gen_Modis_FileName_V2 = SUCCEED

      Return

      End



c---------------------------------------------------------------------------------
      Integer Function Retrieve_OneMeasParm( FileLUN,
     1                                       FileVRSN,
     2                                       Class_MeasParm,
     3                                       Name_MeasParm,
     4                                       PctMissing_MeasParm,
     5                                       AutoFlag_MeasParm,
     6                                       AutoFlagExpl_MeasParm )

      implicit none

      include 'PGS_SMF.f'
      include 'PGS_MET_13.f'
      include 'Atmos_ECSMET.inc'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Retrieves a single MeasuredParameter from an HDF-EOS 
C               file given its class instance index in the file, and 
C               the PCF LUN and Version numbers of the file.
C               The four MeasuredParameter fields common to MODIS
C               atmosphere products are returned.  They are:
C                      PARAMETERNAME
C                      QAPERCENTMISSINGDATA
C                      AUTOMATICQUALITYFLAG
C                      AUTOMATICQUALITYFLAGEXPLANATION
C
C
C !INPUT PARAMETERS:
C  integer       FileLUN   Variable set to PCF LUN of file to be queried
C                          for named MeasuredParameter 
C
C  integer       FileVRSN  Variable set to PCF version number of the 
C                          file to be queried for named 
C                          MeasuredParameter
C                            
C  integer       Class_MeasParm
C                          Variable set the class instance index of the
C                          named MeasuredParameter in the file to be 
C                          queried. 
C
C !OUTPUT PARAMETERS:       
C  character*(*) Name_MeasParm
C                          Variable containing the name of the retrieved 
C                          MeasuredParameter.
C
C  character*(*) AutoFlag_MeasParm 
C                          Variable set to data quality flag of named 
C                          MeasuredParameter.  Valids are "Passed"/
C                          "Suspect"/"Failed".
C 
C  character*(*) AutoFlagExpl_MeasParm
C                          Variable containing the criteria used to set
C                          the data quality flag of the named 
C                          MeasuredParameter. 
C                 
C  integer PctMissing_MeasParm
C                          Percent missing data in field of named 
C                          MeasuredParameter. Note integer value. 
C
C !REVISION HISTORY:
C 09/12/2003 Richard Hucek:  initial version
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
C        PGS_MET_GetPCAttr_i       (libPGSTK.a)
C        PGS_MET_GetPCAttr_s       (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        string_loc                (science code)
C
C
C    Named Constant:
C        FAIL                      (Atmos_ECSMET.inc)
C        FUNCNAME_PGS_MET_GET_i    (Atmos_ECSMET.inc)
C        FUNCNAME_PGS_MET_GET_s    (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'Retrieve_OneMeasParm')


C function arguments declarations
      character*(*) Name_MeasParm,
     1              AutoFlag_MeasParm,
     2              AutoFlagExpl_MeasParm 

      integer       FileLUN, 
     1              FileVRSN, 
     2              Class_MeasParm,
     3              PctMissing_MeasParm


C other variable declarations
      character*25   msg25d,msg25_LUN,msg25_VRSN
      character*60   Name
      character*128  AttrN,AttrV_s,Field_Name(4),HDF_AttrName
      character*512  msgbuf_GET_err
      character*1024 msgbuf

      integer i,i1,rtn,rtn_loc
      integer AttrV_i
      integer fbyte,fbytea,fbyteb,lbyte,lbytea,lbyteb
      integer fbyte_LUN,fbyte_VRSN,lbyte_LUN,lbyte_VRSN

      logical error_flag

C function declarations
      integer String_Loc,
     1        PGS_MET_GetPCAttr_i,PGS_MET_GetPCAttr_s


C------------------------
C Initialization
C------------------------
      Retrieve_OneMeasParm = SUCCEED
      error_flag           = .FALSE.
      HDF_AttrName         = 'CoreMetadata.0' 

      Field_Name(1) = 'PARAMETERNAME' 
      Field_Name(2) = 'QAPERCENTMISSINGDATA'
      Field_Name(3) = 'AUTOMATICQUALITYFLAG'
      Field_Name(4) = 'AUTOMATICQUALITYFLAGEXPLANATION'

c-----------------------------------------------------------------------
c Set up status message variables
c-----------------------------------------------------------------------

      write(msg25_LUN, '(I25)') FileLUN 
      write(msg25_VRSN,'(I25)') FileVRSN 

      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)

c-----Retrieve MeasuredParameter in product file 
         Name = Name_MeasParm

         write(msg25d,'(i25)') Class_MeasParm 
         rtn_loc = String_Loc(msg25d,fbyte,lbyte)
         rtn_loc = String_Loc(Name,fbytea,lbytea)

c--------loop over the 4 ODL objects

         Do i=1,4
            rtn_loc = String_Loc(Field_Name(i),fbyteb,lbyteb)
            AttrN   = Field_Name(i)(fbyteb:lbyteb) // '.' // msg25d(fbyte:lbyte)
            rtn_loc = String_Loc(AttrN,fbyteb,lbyteb)

            msgbuf_GET_err =
     1      ' unable to retrieve value of MeasuredParameter ' 
     2      // char(10) // 'attribute '// AttrN(fbyteb:lbyteb) // ' from file on LUN '
     3      // msg25_LUN(fbyte_LUN:lbyte_LUN) // ' and VRSN ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.' 
     4      // char(10) // 'Operator Action: Check for flawed input product on specified PCF LUN'
     5      // char(10) // 'and VRSN numbers.  Otherwise, notify SDST.'


c-----------ODL object is type integer
            If (i.eq.2) Then
               rtn = PGS_MET_GetPCAttr_i(FileLUN, FileVRSN, HDF_AttrName, AttrN, AttrV_i)

c--------------set integer ODL object
               If (rtn .eq. PGS_S_SUCCESS) Then
                  PctMissing_MeasParm = AttrV_i                  
               Else
                  error_flag = .true.
                  msgbuf = FUNCNAME_PGS_MET_GET_i // msgbuf_GET_err
                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If

c-----------ODL object is type character
            Else
               rtn = PGS_MET_GetPCAttr_s(FileLUN, FileVRSN, HDF_AttrName, AttrN, AttrV_s)

c--------------set character ODL object
               If (rtn .eq. PGS_S_SUCCESS) Then
                  if (i.eq.1) Name_MeasParm         = AttrV_s
                  if (i.eq.3) AutoFlag_MeasParm     = AttrV_s 
                  if (i.eq.4) AutoFlagExpl_MeasParm = AttrV_s 
               Else
                  error_flag = .true.
                  msgbuf = FUNCNAME_PGS_MET_GET_s // msgbuf_GET_err
                  call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If

            End If   ! data type of ODL object (integer or string)

         END DO   ! number of ODL objects per MearuredParameter

      If (error_flag) Retrieve_OneMeasParm = FAIL

      Return

      END



c---------------------------------------------------------------------------------
      Integer Function Retrieve_PSAValues(FileLUN, FileVRSN,
     1                                    NUM_PSA, Name_PSA, Value_PSA)

      implicit  NONE

      include   'PGS_MODIS_39500.f'
      include   'Atmos_ECSMET.inc'

C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:   Retrieve a list of Product Specific Attribute (PSAs) 
C                 values given a list of PSA names and file identifiers
C                 (LUN and VRSN).  See DESIGN NOTES below. 
C
C
C !INPUT PARAMETERS:
C  integer FileLUN    PCF Logical Unit Number (LUN) of file to be 
C                     queried for PSAs. 
C
C  integer FileVRSN   PCF version number of file to be queried 
C
C  integer NUM_PSA    Variable containing the number of additional 
C                     Product-Specific Attributes (PSAs) to be 
C                     retrieved 
C 
C  character*(*) Name_PSA(*)  
C                     Array of PSA names to be retrieved 
C
C
C !OUTPUT PARAMETERS:
C
C  real Value_PSA(*)  Ordered list of PSA values corresponding to input
C                     PSA name list 
C
C
C !REVISION HISTORY:
C   Initial Development:  R. Hucek 4/17/03
C
C
C!TEAM-UNIQUE HEADER:
C
C    This software was developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C
C !REFERENCES AND CREDITS:
C
C    Written by         Richard Hucek 4/2003 
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    rhucek@ltpmail.gsfc.nasa.gov
C
C
C !DESIGN NOTES:
C  Blank PSA names within the first NUM_PSA elements is considered an
C  error.  Retrieve_PSAValues will return FAIL if even one of the 
C  names is blank.  It is also an error to pass a value of NUM_PSA that 
C  is less than or equal to zero.  No initialization of array PSA_Values
C  is performed prior to retrieval of PSA values. 
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_GetPCAttr_s       (libPGSTK.a)
C        MODIS_SMF_SETDYNAMICMSG   (science code)
C        string_loc                (science code)
C
C    Named Constants:
C        FAIL                      (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C-----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'Retrieve_PSAValues')

C Function argument declarations
      character*(*) Name_PSA(*)
      integer       FileLUN, FileVRSN, NUM_PSA 
      real          Value_PSA(*)

C other variable declarations
      character*8    char_buf
      character*25   msg25, msg25_LUN, msg25_VRSN
      character*100  HDF_FileAttrName
      character*1024 msgbuf

      integer i, rtn, rtn_loc
      integer fbyte,fbyte_LUN,fbyte_VRSN,
     1        lbyte,lbyte_LUN,lbyte_VRSN

      logical error_flag

C external function declarations
      Integer pgs_met_getpcattr_s,string_loc


c-----------------------------------------------------------------------
c Initializations
c-----------------------------------------------------------------------

      error_flag         = .false.
      Retrieve_PSAValues = SUCCEED     
      HDF_FileAttrName   = 'CoreMetadata.0'

c-----Set up variables for status messaging
      write(msg25_LUN,'(I25)') FileLUN
      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)

      write(msg25_VRSN,'(I25)') FileVRSN
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)

c-----------------------------------------------------------------------
c Perform input argument checks
c-----------------------------------------------------------------------

c --- check for a positive number of PSAs
      If (NUM_PSA .LE. 0) Then
         error_flag = .true.
         Retrieve_PSAValues = FAIL

         write(msg25,'(i25)') NUM_PSA
         rtn_loc = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'The number of PSAs passed (' // msg25(fbyte:lbyte) // ') is <= 0.'
     1   // char(10) // 'No PSAs retrieved.' 
     2   // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Return
      End If


c-----------------------------------------------------------------------
c  Retrieve PSA values
c-----------------------------------------------------------------------
      DO i = 1, NUM_PSA

         If (Name_PSA(i) .eq. BLANK) Then
            error_flag = .true.

            write(msg25,'(i25)')  i
            rtn = String_Loc(msg25,fbyte,lbyte)

            msgbuf = 'PSA name is blank.  Cannot retrieve PSA array element ' 
     1      //  msg25(fbyte:lbyte) // '.'
     2      // char(10) // 'Operator Action:  Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------PSA name not blank
         Else
            rtn = pgs_met_getpcattr_s( FileLUN,
     1                                 FileVRSN,
     2                                 HDF_FileAttrName,
     3                                 Name_PSA(i),
     4                                 char_buf )

c-----------errror returned by pgs_met_getpcattr_s
            If (rtn .ne. SUCCEED) Then
               error_flag = .true.
               rtn = String_Loc(Name_PSA(i),fbyte,lbyte)

               msgbuf = 'pgs_met_getpcattr_s unable to retrieve PSA ' 
     1         // char(10) // '"' // Name_PSA(i)(fbyte:lbyte) // '" from file on LUN '
     2         // msg25_LUN(fbyte_LUN:lbyte_LUN) // ' and VRSN ' 
     3         // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.'
     4         // char(10) // 'Operator Action: Check for flawed data product on specified PCF LUN'
     5         // char(10) // 'VRSN numbers.  Otherwise, notify SDST.'
   
               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            Else
               read(char_buf,'(f8.2)') Value_PSA(i)
            EndIf

         EndIf

      End Do

      If ( error_flag ) Retrieve_PSAValues = FAIL 

      Return

      End



c---------------------------------------------------------------------------------
      Integer Function Set_AncInputPntr(MET_Handles)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Set one instance of ECS AncillaryInputGranule metadata
C               group to internal memory.  The AncillaryInputGranule
C               group consists of the following 2 ODL objects:
C
C           AncillaryInputType     (set to "Geolocation")
C           AncillaryInputPointer  (contains UR of input Geolocation file)
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
C
C !OUTPUT PARAMETERS:       None
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
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        PGS_PC_GetUniversalRef    (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        FAIL                      (Atmos_ECSMET.inc)
C        INVENTORYMETADATA         (Atmos_ECSMET.inc)
C        LUN_GEO                   (Atmos_ECSMET.inc)
C        MCORE_*                   (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSd_PC_UREF_LENGTH_MAX   (PGS_PC.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C        VRSN_GEO                  (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Set_AncInputPntr')


C function argument declaration
      character*(*)  MET_Handles(*)


C other variable declarations
      character*25   msg25_LUN,msg25_VRSN
      character*128  AttrN,AttrV
      character*1024 msgbuf
      character*(PGSd_PC_UREF_LENGTH_MAX) UR

      integer rtn,rtn_loc,VRSN_local

      integer pgs_met_setattr_s, PGS_PC_GetUniversalRef, String_Loc

      integer fbyte,fbyte_LUN,fbyte_VRSN,
     1        lbyte,lbyte_LUN,lbyte_VRSN

      logical error_flag

C------------------------
C Initialization
C------------------------

      Set_AncInputPntr = FAIL
      error_flag       = .FALSE.

c-------------------------------------------------------------------------------
c set ODL object AncillaryInputType
c-------------------------------------------------------------------------------

      AttrN = MCORE_ANCIL_INPUT_TYPE // '.1'
      AttrV = 'Geolocation'

      rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV)

      If (rtn.ne.PGS_S_SUCCESS) Then
          error_flag = .TRUE.
          rtn = String_Loc(AttrN,fbyte,lbyte)

          msgbuf = 'PGS_MET_SetAttr_s unable to set '// AttrN(fbyte:lbyte)
     1     // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2     // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3     // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

          Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


c-------------------------------------------------------------------------------
c set ODL object AncillaryInputPointer
c 1 - retrieve Geolocation UR from PCF line entry
c 2 - set to memory as AncillaryInputPointer
c-------------------------------------------------------------------------------

c-----retrieve Geolocation UR
      VRSN_local = VRSN_GEO  !WARNING - VRSN_local is reset by PGS_PC_GetUniversalRef

      rtn = PGS_PC_GetUniversalRef(LUN_GEO,VRSN_local,UR)

      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .TRUE.

c--------Set up status message variables
         write(msg25_LUN, '(I25)') LUN_GEO
         write(msg25_VRSN,'(I25)') VRSN_GEO

         rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)
         rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)

         msgbuf = 'PGS_PC_GetUniversalRef unable to retrieve UR from Geolocation '
     1   // char(10) // 'product (LUN = ' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     2   // ' Version_No = '//msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // ').'
     3   // char(10) // 'Operator Action:  Check PCF Universal Reference Identifier.'
     4   // char(10) // 'If invalid, stage correct PCF and rerun PGE.  Otherwise, notify SDST.'


         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c-----successfully read Geolocation UR, now set it
      Else
         AttrN = MCORE_ANCIL_POINTER // '.1'

         rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,UR)

         If (rtn.ne.PGS_S_SUCCESS) Then
             error_flag = .TRUE.
             rtn        = String_Loc(AttrN,fbyte,lbyte)

             msgbuf = 'PGS_MET_SetAttr_s unable to set '// AttrN(fbyte:lbyte)
     1       // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2       // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3       // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

             Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If

      End If

      If (.NOT.error_flag) Set_AncInputPntr = SUCCEED

      Return

      End  ! Set_AncInputPntr



c---------------------------------------------------------------------
      Integer Function Set_ArchPSA( MET_Handles,
     1                              NUM_Of_ArchivePSA_SC,
     2                              Name_ArchivePSA_SC,
     3                              Value_ArchivePSA_SC )

      implicit none

      include 'Atmos_ECSMET.inc'
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
C  integer NUM_Of_ArchivePSA_SC
C                Variable containing the number of archive Product
C                Specific Attributes (PSA) with granule dependent
C                values determined by the science code (SC) to be set
C                by process.  If none, pass NUM_Of_ArchivePSA_SC=0.
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
C                 groups in the MCF.).
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
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        String_Loc                (atmos shared code)
C
C    Named Constant:
C        ARCHIVEMETADATA           (Atmos_ECSMET.inc)
C        FAIL                      (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C-----------------------------------------------------------------------

c-----PARAMETER declarations
      character*(*) FUNCNAME
      parameter (FUNCNAME = 'Set_ArchPSA')

      integer MAX_NUM_PSA
      parameter (MAX_NUM_PSA = 50)

c-----Declaration of function arguments
      character*(*) MET_Handles(*), Name_ArchivePSA_SC(*)
      integer NUM_Of_ArchivePSA_SC
      real Value_ArchivePSA_SC(*)

c-----other variable declarations
      character*60 AttrN
      real         PSA_V

      character*8 msg8
      character*25 msg25, msg25a, msg25b
      character*512 msgbuf

      integer fbyte,fbytea,fbyteb,
     2        lbyte,lbytea,lbyteb,
     3        i,rtn

      integer string_loc
      integer pgs_met_setattr_s

      logical error_flag


c------------------------
c Initialization
c------------------------

      Set_ArchPSA = FAIL
      error_flag = .FALSE.


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
                  rtn = pgs_met_setattr_s(MET_Handles(ARCHIVEMETADATA),AttrN, msg8)

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


      If (.not.error_flag) Set_ArchPSA = SUCCEED

      Return

      End  !Set_ArchPSA



c---------------------------------------------------------------------------------
      Integer Function Set_GRing(MET_Handles, MET_MasterGroup_ID)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Copy ECS GRING metadata attributes from the MODIS
C                Geolocation product to internal SDPTK memory.
C
C
C !INPUT PARAMETERS:
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
C  integer MET_MasterGroup_ID
C                 Array index referencing MCF MASTER group element to
C                 which GRING attributes are to be set.  2 -> Inventory;
C                 3 -> Archive. (See MET_Handles above.)
C
C
C !OUTPUT PARAMETERS: None.
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
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        String_Loc                (atmos shared code)
C
C    Named Constant:
C        FAIL                      (Atmos_ECSMET.inc)
C        LUN_GEO                   (Atmos_ECSMET.inc)
C        MCORE_*                   (mapi.inc)
C        MECS_CORE                 (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSd_PC_VALUE_LENGTH_MAX  (PGS_PC.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C        VRSN_GEO                  (Atmos_ECSMET.inc)
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Set_GRing')


c-----Declaration of function arguments
      character*(*) MET_Handles(*)
      integer       MET_MasterGroup_ID


c-----Declaration of local variables
      character*25   msg25a, msg25b
      character*60   AttrN, Field_Name(4)
      character*1012 msgbuf_GET_err,msgbuf_SET_err
      character*(PGSd_PC_VALUE_LENGTH_MAX) buf_s

      double precision buf_dbl(4)

      integer pgs_met_setattr_d, pgs_met_setattr_s, pgs_met_setattr_i,
     1        pgs_met_getpcattr_d, pgs_met_getpcattr_s,
     2        pgs_met_getpcattr_i, string_loc

      integer fbyte,fbytea,fbyteb,
     1        lbyte,lbytea,lbyteb,
     2        buf_i(4),i,rtn


      logical error_flag


C------------------------
C Initialization
C------------------------

      Set_GRing = FAIL
      error_flag = .false.


C----------------------------------------------------------------------------
C ECS GPolygon group attributes include the four metadata objects:
C
C 1 - ExclusionGRingFlag
C 2 - GRingPointLatitude
C 3 - GRingPointLongitude
C 4 - GRingPointSequenceNo
C----------------------------------------------------------------------------

      Field_Name(1) = MCORE_EXCLUS_GRING_FLG
      Field_Name(2) = MCORE_GRING_POINT_LAT
      Field_Name(3) = MCORE_GRING_POINT_LON
      Field_Name(4) = MCORE_GRING_POINT_NUM

c-----set message buffers for error reporting
      write(msg25a,'(I25)') LUN_GEO
      write(msg25b,'(I25)') VRSN_GEO
      rtn = String_Loc(msg25a,fbytea,lbytea)
      rtn = String_Loc(msg25b,fbyteb,lbyteb)

      Do 100 i = 1, 4
         AttrN = Field_Name(i)

         rtn = string_loc(AttrN,fbyte,lbyte)

         msgbuf_GET_err = 'pgs_met_getpcattr_s unable to retrieve ECS attribute ' // AttrN(fbyte:lbyte)
     1   // char(10) // 'from MODIS geolocation product on LUN = ' // msg25a(fbytea:lbytea)
     2   // ';  Version number = ' // msg25b(fbyteb:lbyteb)
     3   // char(10) // AttrN(fbyte:lbyte) // ' not set.'
     4   // char(10) // 'Operator Action:  Check for corrupted input file or PCF '
     5   // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     6   // char(10) // 'fault is identified, stage correct PCF/input file and '
     7   // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


         msgbuf_SET_err = 'pgs_met_setattr_s unable to set G-Ring attribute ' // AttrN(fbyte:lbyte)
     1   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'


c--------Retrieve and set ExclusionGRingFlag (Y=inner, N=outer)
         If (i .eq. 1) Then
            rtn = pgs_met_getpcattr_s(LUN_GEO, VRSN_GEO, MECS_CORE, AttrN, buf_s)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf_GET_err,FUNCNAME)
            Else
               rtn = pgs_met_setattr_s(MET_Handles(MET_MasterGroup_ID), AttrN, buf_s)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf_SET_err,FUNCNAME)
               End If
            End If


c--------Retrieve and set GRingPointLatitude and GRingPointLongitude
         Else If ( (i.eq.2) .or. (i.eq.3) ) Then
            rtn = pgs_met_getpcattr_d(LUN_GEO, VRSN_GEO, MECS_CORE, AttrN, buf_dbl)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf_GET_err,FUNCNAME)

            Else
               rtn = pgs_met_setattr_d(MET_Handles(MET_MasterGroup_ID), AttrN, buf_dbl)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf_SET_err,FUNCNAME)
               End If
            End If

c--------Retrieve and set GRingPointSequenceNo
         Else If (i.eq.4) Then

            rtn = pgs_met_getpcattr_i(LUN_GEO, VRSN_GEO, MECS_CORE, AttrN, buf_i)

c-----------set error flag if unable to retrieve metadata value
            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf_GET_err,FUNCNAME)

            Else
               rtn = pgs_met_setattr_i(MET_Handles(MET_MasterGroup_ID), AttrN, buf_i)

               If (rtn.ne.PGS_S_SUCCESS) Then
                   error_flag = .true.

                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf_SET_err,FUNCNAME)
               End If

            End If

         End If

100   Continue



      If (.not.error_flag) Set_GRing = SUCCEED

      Return

      End  !Set_GRing



c---------------------------------------------------------------------------------
      Integer Function Set_GeoData_Atmos(MET_Handles)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Set (to internal SDPTK memory) the MODIS atmosphere
C                product ECS inventory metadata that can be copied from
C                the MODIS Geolocation product.  (See design notes for
C                a list of these attributes.)
C
C
C !INPUT PARAMETERS:
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
C !OUTPUT PARAMETERS:       None
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
C  The following ECS metadata attributes are copied from the MODIS
C  Geolocation product:
C
C "DayNightFlag"
C "EastBoundingCoordinate"
C "NorthBoundingCoordinate"
C "SouthBoundingCoordinate"
C "WestBoundingCoordinate"
C "RangeBeginningTime"
C "RangeEndingTime"
C "RangeBeginningDate"
C "RangeEndingDate"
C "OrbitNumber"
C "EquatorCrossingDate"
C "EquatorCrossingTime".
C "EquatorCrossingLongitude"
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        FAIL                      (Atmos_ECSMET.inc)
C        INVENTORYMETADATA         (Atmos_ECSMET.inc)
C        LUN_GEO                   (Atmos_ECSMET.inc)
C        MCORE_*                   (mapi.inc)
C        MECS_CORE                 (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C        VRSN_GEO                  (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Set_GeoData_Atmos')

C Function argument declarations
      character*(*) MET_Handles(*)


C other variable declarations
      character*25   msg25_LUN,msg25_VRSN
      character*128  AttrN,Field_Name(8),buf_char
      character*1024 msgbuf
      character*512  msgbuf_GET_err, msgbuf_SET_err

      integer i,rtn,rtn_loc
      integer buf_int(4)

      integer pgs_met_setattr_d,
     1        pgs_met_setattr_s, pgs_met_setattr_i,
     2        PGS_MET_GetPCAttr_i,PGS_MET_GetPCAttr_d,
     3        PGS_MET_GetPCAttr_s, string_Loc

      integer fbyte,fbyte_LUN,fbyte_VRSN,
     1        lbyte,lbyte_LUN,lbyte_VRSN

      double precision buf_dbl(4)

      logical error_flag


C-------------------------------------------------------------------------------
C Initialization
C-------------------------------------------------------------------------------

      Set_GeoData_Atmos = FAIL
      error_flag        = .false.

      Do 10 i = 1, 4
         buf_int(i) = 0
         buf_dbl(i) = 0.0D0
   10 continue

c-----Set up variables for status messaging
      write(msg25_LUN,'(I25)')  LUN_GEO
      write(msg25_VRSN,'(I25)') VRSN_GEO

      rtn_loc = String_Loc(msg25_LUN,fbyte_LUN,lbyte_LUN)
      rtn_loc = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)

      msgbuf_GET_err =
     1            char(10) // 'from Geolocation product on LUN = '//msg25_LUN(fbyte_LUN:lbyte_LUN)
     2         //  ' and Version number = ' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN)
     3         // char(10) // 'Operator Action:  Check for corrupted input and/or PCF '
     4         // char(10) // 'file, or incorrect PCF reference to input file.  If a '
     5         // char(10) // 'fault is identified, stage correct PCF/input file and '
     6         // char(10) // 'rerun PGE.  Otherwise, notify SDST.'


      msgbuf_SET_err =
     1             char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2          // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3          // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'


C-------------------------------------------------------------------------------
C Set Bounding Coordinates
C-------------------------------------------------------------------------------

      Field_Name(1) = MCORE_EAST_BOUND
      Field_Name(2) = MCORE_NORTH_BOUND
      Field_Name(3) = MCORE_SOUTH_BOUND
      Field_Name(4) = MCORE_WEST_BOUND

      Do 20 i = 1, 4
         AttrN = Field_Name(i)
         rtn   = String_Loc(AttrN,fbyte,lbyte)

         rtn = PGS_MET_GetPCAttr_d(LUN_GEO,VRSN_GEO,MECS_ARCHIVE,AttrN,buf_dbl)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.

            msgbuf =
     1      'PGS_MET_GetPCAttr_d unable to retrieve ECS attribute '// AttrN(fbyte:lbyte)
     2      // msgbuf_GET_err

            call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else
            rtn = PGS_MET_SetAttr_d(MET_Handles(INVENTORYMETADATA),AttrN,buf_dbl)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf =
     1         'PGS_MET_SetAttr_d unable to set ECS attribute ' // AttrN(fbyte:lbyte)
     2         // ' retrieved ' // msgbuf_SET_err

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            End If

         End If

   20 Continue


C-------------------------------------------------------------------------------
C Set DAYNIGHTFLAG, RANGEBEGINNING/ENDING DATE & TIME, and EQUATORCROSSING
C DATE & TIME.
C-------------------------------------------------------------------------------

      Field_Name(1) = MCORE_DAYNIGHTFLAG
      Field_Name(2) = MCORE_RANGE_BEG_TIME
      Field_Name(3) = MCORE_RANGE_ENDING_TIME
      Field_Name(4) = MCORE_RANGE_BEG_DATE
      Field_Name(5) = MCORE_RANGE_ENDING_DATE
      Field_Name(6) = MCORE_EQUATCROSSINGDATE//'.1'
      Field_Name(7) = MCORE_EQUATCROSSINGTIME//'.1'

      Do 25 i = 1, 7
         buf_char = ' '
         AttrN    = Field_Name(i)
         rtn      = String_Loc(AttrN,fbyte,lbyte)

         rtn = PGS_MET_GetPCAttr_s(LUN_GEO,1,MECS_CORE,AttrN,buf_char)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.

            msgbuf =
     1      'PGS_MET_GetPCAttr_s unable to retrieve ECS attribute '// AttrN(fbyte:lbyte)
     2      // msgbuf_GET_err

            call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else
            rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,buf_char)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf =
     1         'PGS_MET_SetAttr_s unable to set ECS attribute ' // AttrN(fbyte:lbyte)
     2         // ' retrieved ' // msgbuf_SET_err

               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            End If

         End If

   25 Continue


C-------------------------------------------------------------------------------
C Set ORBITNUMBER
C-------------------------------------------------------------------------------

      AttrN = MCORE_ORBIT_NUM//'.1'
      rtn   = String_Loc(AttrN,fbyte,lbyte)

      rtn = PGS_MET_GetPCAttr_i(LUN_GEO,1,MECS_CORE,AttrN,buf_int)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf =
     1   'PGS_MET_GetPCAttr_i unable to retrieve ECS attribute '// AttrN(fbyte:lbyte)
     2   // msgbuf_GET_err


         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      Else
         rtn = PGS_MET_SetAttr_i(MET_Handles(INVENTORYMETADATA),AttrN,buf_int)

         If (rtn.ne.PGS_S_SUCCESS) Then
             error_flag = .true.

             msgbuf =
     1       'PGS_MET_SetAttr_i unable to set ECS attribute ' // AttrN(fbyte:lbyte)
     2       // ' retrieved ' // msgbuf_SET_err

             call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If

      End If


C-------------------------------------------------------------------------------
C Set ORBITCROSSINGLONGITUDE
C-------------------------------------------------------------------------------

      AttrN = MCORE_EQUATCROSSINGLONG//'.1'
      rtn = String_Loc(AttrN,fbyte,lbyte)

      rtn = PGS_MET_GetPCAttr_d(LUN_GEO,1,MECS_CORE,AttrN,buf_dbl)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf =
     1   'PGS_MET_GetPCAttr_d unable to retrieve ECS attribute '// AttrN(fbyte:lbyte)
     2    // msgbuf_GET_err

         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      Else
         rtn = PGS_MET_SetAttr_d(MET_Handles(INVENTORYMETADATA),AttrN,buf_dbl)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.
            rtn = String_Loc(AttrN,fbyte,lbyte)

            msgbuf =
     1      'PGS_MET_SetAttr_d unable to set ECS attribute ' // AttrN(fbyte:lbyte)
     2      // ' retrieved ' // msgbuf_SET_err

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         End If

      End If



      If (.not.error_flag) Set_GeoData_Atmos = SUCCEED

      Return

      End  !Set_GeoData_Atmos



c---------------------------------------------------------------------------------
      Integer Function Set_InputPntr_Atmos( NUM_InputPntr,
     1                                      LRN_InputPntr,
     2                                      VRSN_InputPntr,
     3                                      MET_Handles, NDVIfile,
     4                                      thresfile)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:   Set to internal memory the array of file identifiers
C                 (or Universal References) for inputs used by the
C                 current process.   This array is to be written to the
C                 product file as ECS Inventory metadata attribute
C                 INPUTPOINTER in a subsequent processing step.
C
C
C !INPUT PARAMETERS:
C
C  integer  NUM_InputPntr   Number of input datasets for granule.  This
C                           includes ancillary data files, look up
C                           tables, and other MODIS product files.
C                           System files, such as MCF, are not included.
C
C  integer  LRN_InputPntr(*)
C                           Array containing the PCF LRNs of all input
C                           datasets for granule (See NUM_InputPntr above).
C
C                           A one-to-one correspondence between the elements
C                           of the arrays LRN_InputPntr and VRSN_InputPntr
C                           is assumed.
C
C  integer  VRSN_InputPntr(*)
C                           Array containing the PCF version numbers of all
C                           input datasets for granule.  (See NUM_InputPntr
C                           and LRN_InputPntr above)
C
C
C  character MET_Handles(20)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           In FORTRAN, element 1 of the MET_Handles
C                           array is reserved for the MCF file.   ECS
C                           inventory metadata is referenced as element
C                           2 and archive metadata as element 3.
C
C
C !OUTPUT PARAMETERS:       None
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
C
C !DESIGN NOTES:
C
C  Returns:     SUCCEED if successful, FAIL if an error occurs
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        BLANK                     (Atmos_ECSMET.inc)
C        FAIL                      (Atmos_ECSMET.inc)
C        MAX_NUM_PNTRS             (Atmos_ECSMET.inc)
C        MCORE_INPUT_POINTER       (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSd_PC_UREF_LENGTH_MAX   (PGS_PC.f)
C        PGSd_MET_STR_END          (PGS_MET.f: included in mapi.inc)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)   FUNCNAME
      PARAMETER      (FUNCNAME = 'Set_InputPntr_Atmos')


C function argument declarations
      character*(*) MET_Handles(*)

      integer NUM_InputPntr, LRN_InputPntr(*),VRSN_InputPntr(*)


C other variable declarations
      character*25  msg25,msg25b
      character*512 msgbuf
      character*34 ndvifile
      character*23 thresfile
      character*(PGSd_PC_UREF_LENGTH_MAX) UR_Ref(MAX_NUM_PNTRS)

      integer pgs_met_setattr_s, String_Loc
      integer fbyte, fbyteb,
     1        lbyte, lbyteb
      integer i, NUM_Pointers, NUM_nonBlank_UR, NUM_UR, rtn, rtn_flag

      logical error_flag

C------------------------
C Initialization
C------------------------

      Set_InputPntr_Atmos =  FAIL
      error_flag          = .FALSE.


C-----------------------------------------------------------------
C Set Universal References (URs) or "InputPointers" to input data
C set.  URs are retrieved by reading field 5 of PCF file line
C entries.
C-----------------------------------------------------------------

      NUM_Pointers = NUM_InputPntr

c-----first check if number of input pointers exceed Set_CoreMetadata internal limit.
      If (NUM_Pointers .gt. MAX_NUM_PNTRS) Then
         error_flag   = .true.
         NUM_Pointers = MAX_NUM_PNTRS

         write(msg25,'(i25)') NUM_Pointers
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'Number of Universal References exceeds '
     1   // char(10) // 'number of elements that can be set by Set_InputPntr_Atmos.'
     2   // char(10) // 'Only ' // msg25(fbyte:lbyte) // ' URs set.'
     3   // char(10) // 'Operator Action: Notify SDST. '

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


c-----retrieve input URs and report any problem to log file.
      Call Retrieve_UR_Set( NUM_Pointers,
     1                      LRN_InputPntr,
     2                      VRSN_InputPntr,
     3                      NUM_UR,
     4                      UR_Ref,
     5                      rtn_flag)

      If (rtn_flag .ne. SUCCEED) Then
         error_flag = .true.

         msgbuf = 'Retrieve_UR_Set detected error acquiring process URs.'
     1   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     2   // char(10) // 'messages originating within subroutine Retrieve_UR_Set.'

         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End IF


c-----loop over retrieved URs and remove BLANK URs.
      NUM_nonBlank_UR = 0

      Do 100 i = 1, NUM_UR

         If (UR_Ref(i) .eq. BLANK) Then
            error_flag = .TRUE.

            msgbuf =
     1      'UR passed from Retrieve_UR_Set is blank - not aggregated to InputPointer array. '
     2      // char(10) // 'Operator Action:  Refer to prior "WARNING" messages from routine '
     3      // char(10) // 'Retrieve_UR_Set to identify PCF file line entry and required '
     4      // char(10) // 'actions. '

            Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         Else
            NUM_nonBlank_UR         = NUM_nonBlank_UR + 1
            if (i.eq.6) then
               UR_Ref(NUM_nonBlank_UR) = thresfile
            else
               UR_Ref(NUM_nonBlank_UR) = UR_Ref(i)
            end if
         End If

 100  Continue   ! loop over URs and remove BLANK ones

c----- add in NDVI file name from LEOCAT processing 
      NUM_nonBlank_UR         = NUM_nonBlank_UR + 1
      UR_Ref(NUM_nonBlank_UR) = ndvifile

c-----add end-of-data marker to UR array, space permitting
      If (NUM_nonBlank_UR .lt. MAX_NUM_PNTRS) UR_Ref(NUM_nonBlank_UR+1) = PGSd_MET_STR_END

      rtn   = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),MCORE_INPUT_POINTER,UR_Ref)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag =  .TRUE.

         msgbuf = 'PGS_MET_SetAttr_s detected error setting ' // MCORE_INPUT_POINTER // '.'
     1   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c rhucek 09/02/03: commented out success message block
c     Else
c        write(msg25b,'(i25)') NUM_nonBlank_UR
c        rtn = String_Loc(msg25b,fbyteb,lbyteb)
c
c        msgbuf = msg25b(fbyteb:lbyteb) // ' URs set to ' // MCORE_INPUT_POINTER // '. '
c
c        Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,msgbuf,FUNCNAME)
      End If


      If (.NOT.error_flag) Set_InputPntr_Atmos = SUCCEED

      Return

      End  !Set_InputPntr_Atmos



c---------------------------------------------------------------------------------
      Integer Function Set_InvPSA_Atmos( Num_PSA,
     1                                   Name_PSA,
     2                                   Value_PSA,
     3                                   MET_Handles )

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Set to internal memory the Inventory PSAs (Product
C                Specific Attributes) defined by the current process.
C                These data are to be written to the product file as
C                ECS metadata fields ADDITIONALATTRIBUTENAME and
C                PARAMETERVALUE in a subsequent processing step.
C
C !INPUT PARAMETERS:
C
C  integer Num_PSA          Variable containing the number of
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
C  character MET_Handles(20)
C                           Array containing the names of the "MASTER"
C                           groups as defined in the MCF.  There may be
C                           up to 20 "MASTER" groups in the MCF.
C
C                           In FORTRAN, element 1 of the MET_Handles
C                           array is reserved for the MCF file.   ECS
C                           inventory metadata is referenced as element
C                           2 and archive metadata as element 3.
C
C
C !OUTPUT PARAMETERS:       None
C
C
C !REVISION HISTORY:
C rhucek 07/18/03:  Commented out Success messages previously written to LogStatus
C   file when setting individual PSAs
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
C  Returns:     SUCCEED if successful, FAIL if an error occurs
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        Set_PSA                   (atmos shared code)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        FAIL                      (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        MODIS_W_GENERIC           (PGS_MODIS_39500.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Set_InvPSA_Atmos')

      real   NoSet_Flag, Rel_equality_EPS
      parameter ( Rel_equality_EPS=0.000001, NoSet_Flag = -99999.0 )

C other variable declarations
      integer       Num_PSA
      character*(*) Name_PSA(*)
      real          Value_PSA(*)
      character*(*) MET_Handles(*)

      character*8   msg8
      character*25  msg25,msg25b
      character*512 msgbuf

      integer  String_Loc, Set_PSA

      integer  fbyte,fbyteb,fbytec,
     1         lbyte,lbyteb,lbytec,
     2         icounter,rtn

      logical error_flag


C------------------------
C Initialization
C------------------------

      Set_InvPSA_Atmos = FAIL
      error_flag = .false.


C-----------------------------------------------------------------------
C  Set Product Specific Attributes (PSAs).  There are two fields:
C  "AdditionalAttributeName"  & "ParameterValue"
C-----------------------------------------------------------------------

C --- check for positive number of PSAs
      If (Num_PSA .LE. 0) Then
         write(msg25,'(i25)') Num_PSA
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'The number of PSAs passed in this call is ' // msg25(fbyte:lbyte) // '.'
     1   // char(10) // 'No PSAs set by Set_InvPSA_Atmos.'


         Call MODIS_SMF_SETDYNAMICMSG(MODIS_N_GENERIC,msgbuf,FUNCNAME)

c-----number of PSAs is positive
      Else

         Do 50 icounter = 1, Num_PSA
            write(msg25,'(i25)') icounter
            rtn = String_Loc(msg25,fbyte,lbyte)

            write(msg25b,'(E25.3)') Value_PSA(icounter)
            rtn = String_Loc(msg25b,fbyteb,lbyteb)
            rtn = String_Loc(Name_PSA(icounter),fbytec,lbytec)


            If ( abs( (Value_PSA(icounter) - NoSet_Flag)/NoSet_Flag) .gt. Rel_equality_EPS) Then

c--------------check for blank PSA name
               If (Name_PSA(icounter) .eq. ' ') Then
                  error_flag = .true.

                  msgbuf = 'PSA name is blank.' // char(10) // 'Set_InvPSA_Atmos unable to set '
     1            // msg25(fbyte:lbyte) // 'th element of PSA array.'
     2            // char(10) // 'Operator Action:  Notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------------PSA value out of range
               Else If ( Abs(Value_PSA(icounter)) .gt. MAX_ABS_VALUE_PSA) Then
                  error_flag = .true.

                  msgbuf = 'PSA value ' // msg25b(fbyteb:lbyteb) // ' is out of bounds.' // char(10)
     1            // 'PSA ' // Name_PSA(icounter)(fbytec:lbytec) // ' not set.'
     2            // char(10) // 'Operator Action:  Notify SDST.'


                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------------set PSA name and value
               Else
                 write(msg8, '(f8.2)') Value_PSA(icounter)

                 rtn = Set_PSA(MET_Handles, Name_PSA(icounter), icounter, msg8)

                 If (rtn .EQ. FAIL) Then
                    error_flag = .TRUE.

                    msgbuf = 'Set_PSA detected error setting PSA '
     1              // Name_PSA(icounter)(fbytec:lbytec) // ' Name/Value pair. '
     2              // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     3              // char(10) // 'messages originating within routine Set_PSA. '

                    Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c rhucek 07/18/03: commented out PSA success message block 
c                Else
c                   msgbuf = 'PSA ' // Name_PSA(icounter)(fbytec:lbytec) // ' Name/Value pair successfully set.'
c
c                   Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
                 End If

               End If        ! check for blank PSA name

            End If       ! check for PSA NoSet_Flag

50       continue    ! Loop on PSAs

      End If  ! check for PSA < 0



      If (.not.error_flag) Set_InvPSA_Atmos = SUCCEED

      Return

      End  !Set_InvPSA_Atmos



c---------------------------------------------------------------------------------
      Integer Function Set_LclGranID(MET_Handles, LUN_Of_LclVrsnID)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Construct and set ECS LOCALGRANULEID inventory metadata
C                attribute to internal memory.  The LOCALGRANULEID is the
C                granule name in MODIS file naming convention.  See also
C                design notes.
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
C  integer LUN_Of_LclVrsnID PCF RP LUN containing the value of LOCALVERSIONID.
C
C
C !OUTPUT PARAMETERS:  None
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
C  The code design assumes that ECS attribute LOCALVERSIONID (used
C  to construct the LOCALGRANULEID) is implemented as a PCF Runtime
C  Parameter (RP) - thus the need to pass the RP LUN of the
C  LOCALVERSIONID.
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        Gen_Modis_FileName        (atmos shared code)
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        PGS_MET_GetSetAttr_s      (libPGSTK.a)
C        PGS_PC_GetConfigData      (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        BLANK                     (Atmos_ECSMET.inc)
C        FAIL                      (Atmos_ECSMET.inc)
C        INVENTORYMETADATA         (Atmos_ECSMET.inc)
C        MCORE_SHORT_NAME          (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        PGSd_PC_VALUE_LENGTH_MAX  (PGS_PC.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Set_LclGranID')


C Function argument declarations
      character*(*) MET_Handles(*)
      integer       LUN_Of_LclVrsnID


C other variable declarations
      character*8   ESDT_Name

      character*25  msg25
      character*128 FileName,ECS_AttrN
      character*512 msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) LclVrsnID

      integer    rtn, rtn_loc

      integer    pgs_met_setattr_s,pgs_met_getsetattr_s,
     1           PGS_PC_GetConfigData,
     2           Gen_Modis_FileName,String_Loc

      integer    fbyte,fbyte_name,
     1           lbyte,lbyte_name

      logical    error_flag


C------------------------
C Initialization
C------------------------

      Set_LclGranID =  FAIL
      FileName      =  BLANK
      error_flag    = .FALSE.

c-----set message buffers for error reporting
      write(msg25,'(I25)') LUN_Of_LclVrsnID
      rtn = string_loc(msg25,fbyte,lbyte)


C-----------------------------------------------------------------------
C Set "LocalGranuleID" in 4 steps.  The LocalGranuleID is the product
C file name in MODIS file naming convention.
C
C 1 - retrieve ESDT (Earth Science Data Type) SHORTNAME from the MCF.
C 2 - retrieve LOCALVERSIONID from PCF
C 3 - construct MODIS file name
C 4 - set LocalGranuleID
C-----------------------------------------------------------------------

      ECS_AttrN = MCORE_SHORT_NAME

      rtn = pgs_met_getsetattr_s(MET_Handles(INVENTORYMETADATA),ECS_AttrN,ESDT_Name)

c.....failed to read ESDT_Name
      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .TRUE.

         msgbuf =
     1   'pgs_met_getsetattr_s unable to read the ESDT shortname from product MCF.'
     2   // char(10) // 'Operator Action:  Check for correct MCF.  If wrong, stage '
     3   // char(10) // 'correct MCF and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c.....Oh dear, name is blank
      Else If (ESDT_Name .eq. BLANK) Then
         error_flag = .TRUE.

         msgbuf =
     1   'ESDT_Name retrieve from product MCF is blank.'
     2   // char(10) // 'Operator Action:  Check for correct MCF.  If wrong, stage '
     3   // char(10) // 'correct MCF and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c.....Retrieve value of LOCALVERSIONID from PCF runtime parameter.
      Else

         rtn = pgs_pc_getconfigdata(LUN_Of_LclVrsnID,LclVrsnID)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.

            msgbuf =
     1      'pgs_pc_getconfigdata unable to read LOCALVERSIONID on input argument '
     2      // char(10) // '"LUN_Of_LclVrsnID" (=' // msg25(fbyte:lbyte) // ').'
     3      // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     4      // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     5      // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else

c...........Construct LOCALGRANULEID
c...........note: function Gen_Modis_FileName does its own input argument checking
            rtn = Gen_Modis_FileName(ESDT_Name,LclVrsnID,FileName)

            If (rtn.ne.SUCCEED) Then
               error_flag = .TRUE.
               rtn_loc    =  string_loc(ESDT_Name,fbyte_name,lbyte_name)

               msgbuf =
     1         'Gen_Modis_FileName detected error constructing MODIS file name '
     2         // char(10) // 'for ESDT ' // ESDT_Name(fbyte_name:lbyte_name) // '.'
     3         // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4         // char(10) // 'messages originating in function Gen_Modis_FileName.'


               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else

c .............Set LOCALGRANULEID
               ECS_AttrN = MCORE_LOCALGRANULEID

               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),ECS_AttrN,FileName)

               If (rtn .ne. PGS_S_SUCCESS) Then
                  error_flag = .TRUE.
                  rtn = String_Loc(ECS_AttrN,fbyte,lbyte)

                  msgbuf = 'PGS_MET_SetAttr_s unable to set '// ECS_AttrN(fbyte:lbyte)
     1             // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2             // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3             // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'



                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               End If   ! set LOCALGRANULEID

            End If   ! retrieve MODIS file name

         End If   ! retrieve LOCALVERSIONID from PCF

      End If   ! retrieve SHORTNAME from MCF


      If (.not.error_flag) Set_LclGranID = SUCCEED

      Return

      End  !Set_LclGranID



c---------------------------------------------------------------------------------
      Integer Function Set_LclGranID_MCF( MET_Handles )

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  
C   This second version of Set_LclGranID, like the original, generates 
C   a file name in MODIS convention (referred to as the LOCALGRANULEID) 
C   and sets the value to shared memory.  It differs from the original 
C   by using the MODIS data collection ID retrieved from the Metadata 
C   Configuration File (MCF) as the file version identifier.  The 
C   original code used Runtime Parameter value read in the Process 
C   Control File (PCF).  See also design notes.
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
C
C !OUTPUT PARAMETERS:  None
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
C  This code used the MODIS data collection version as the file
C  version identifier (a 3 character string) in the MODIS LOCALGRANULEID. 
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        Gen_Modis_FileName        (atmos shared code)
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        PGS_MET_GetSetAttr_s      (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        BLANK                     (Atmos_ECSMET.inc)
C        FAIL                      (Atmos_ECSMET.inc)
C        INVENTORYMETADATA         (Atmos_ECSMET.inc)
C        MCORE_SHORT_NAME          (mapi.inc)
C        MCORE_VERSIONID           (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        PGSd_PC_VALUE_LENGTH_MAX  (PGS_PC.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C-----------------------------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Set_LclGranID_MCF')


C Function argument declarations
      character*(*) MET_Handles(*)


C other variable declarations
      character*8   ESDT_Name

      character*3   My_VersionID 
      character*25  msg25, msg25a
      character*128 FileName,ECS_AttrN
      character*512 msgbuf

      integer    VersionID, string_len, rtn, rtn_loc

      integer    pgs_met_setattr_s,pgs_met_getsetattr_i,
     1           pgs_met_getsetattr_s,
     2           Gen_Modis_FileName,String_Loc

      integer    fbyte,fbyte_name,fbytea,
     1           lbyte,lbyte_name,lbytea

      logical    error_flag


C------------------------
C Initialization
C------------------------

      Set_LclGranID_MCF =  FAIL
      FileName          =  BLANK
      error_flag        = .FALSE.


C-----------------------------------------------------------------------------------------
C Set "LocalGranuleID" in 4 steps.  The LocalGranuleID is the product
C file name in MODIS convention.
C
C 1 - retrieve ESDT (Earth Science Data Type) SHORTNAME from the MCF.
C 2 - retrieve VERSIONID from the MCF 
C 3 - construct MODIS file name
C 4 - set LocalGranuleID
C-----------------------------------------------------------------------------------------

      ECS_AttrN = MCORE_SHORT_NAME

      rtn = pgs_met_getsetattr_s(MET_Handles(INVENTORYMETADATA),ECS_AttrN,ESDT_Name)

c.....failed to read ESDT_Name
      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .TRUE.

         msgbuf =
     1   'pgs_met_getsetattr_s unable to read the ESDT shortname from product MCF.'
     2   // char(10) // 'Operator Action:  Check for correct MCF.  If wrong, stage '
     3   // char(10) // 'correct MCF and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c.....Oh dear, name is blank
      Else If (ESDT_Name .eq. BLANK) Then
         error_flag = .TRUE.

         msgbuf =
     1   'ESDT_Name retrieve from product MCF is blank.'
     2   // char(10) // 'Operator Action:  Check for correct MCF.  If wrong, stage '
     3   // char(10) // 'correct MCF and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


C-----------------------------------------------------------------------------------------
c.....Retrieve MODIS Collection VERSIONID from the MCF.
C-----------------------------------------------------------------------------------------
      Else

         ECS_AttrN = MCORE_VERSIONID
         rtn = pgs_met_getsetattr_i(MET_Handles(INVENTORYMETADATA),ECS_AttrN,VersionID)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.

            msgbuf =
     1      'pgs_met_getsetattr_i unable to read VERSIONID from product MCF.'
     2      // char(10) // 'Operator Action:  Check for correct MCF.  If wrong, stage '
     3      // char(10) // 'correct MCF and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else
            write(msg25,'(I25)') VersionID
            rtn = string_loc(msg25,fbyte,lbyte)
            
            string_len = lbyte - fbyte + 1

c...........string length greater than 3 characters
            If (string_len .GT. 3) Then
               error_flag = .true.
            
               write(msg25a,'(I25)') string_len
               rtn = String_Loc(msg25a,fbytea,lbytea)
 
               msgbuf = 'VersionID contains ' // msg25a(fbytea:lbytea)
     1         // ' characters, too long to '
     2         // char(10) // 'comply with 3-digit MODIS version number convention.'
     3         // char(10) // 'Operator Action:  Check for valid MCF file. If wrong or ' 
     4         // char(10) // 'corrupted, stage correct MCF/PCF and rerun PGE. ' 
     5         // 'Otherwise, notify SDST.'

 
               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else
               If (string_len .eq. 1) My_VersionID = '00' // msg25(fbyte:lbyte)
               If (string_len .eq. 2) My_VersionID =  '0' // msg25(fbyte:lbyte)
               If (string_len .eq. 3) My_VersionID = msg25(fbyte:lbyte)

            End If   ! check for 3 digits or less


C-----------------------------------------------------------------------------------------
c...........Generate LOCALGRANULEID
c           note: function Gen_Modis_FileName does its own input argument checking
C-----------------------------------------------------------------------------------------
            rtn = Gen_Modis_FileName(ESDT_Name,My_VersionID,FileName)

            If (rtn.ne.SUCCEED) Then
               error_flag = .TRUE.
               rtn_loc    =  string_loc(ESDT_Name,fbyte_name,lbyte_name)

               msgbuf =
     1         'Gen_Modis_FileName detected error constructing MODIS file name '
     2         // char(10) // 'for ESDT ' // ESDT_Name(fbyte_name:lbyte_name) // '.'
     3         // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4         // char(10) // 'messages originating in function Gen_Modis_FileName.'


               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else

C-----------------------------------------------------------------------------------------
c .............Set LOCALGRANULEID
C-----------------------------------------------------------------------------------------
               ECS_AttrN = MCORE_LOCALGRANULEID

               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),ECS_AttrN,FileName)

               If (rtn .ne. PGS_S_SUCCESS) Then
                  error_flag = .TRUE.
                  rtn = String_Loc(ECS_AttrN,fbyte,lbyte)

                  msgbuf = 'PGS_MET_SetAttr_s unable to set '// ECS_AttrN(fbyte:lbyte)
     1             // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2             // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3             // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               End If   ! set LOCALGRANULEID

            End If   ! retrieve MODIS file name

         End If   ! retrieve LOCALVERSIONID from PCF

      End If   ! retrieve SHORTNAME from MCF


      If (.not.error_flag) Set_LclGranID_MCF = SUCCEED

      Return

      End  !Set_LclGranID_MCF



c---------------------------------------------------------------------------------
      Integer Function Set_LclGranID_Sfx( MET_Handles, 
     1                                    LUN_Of_LclVrsnID, 
     2                                    FileNameSuffix )

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Sets ECS LOCALGRANULEID inventory metadata attribute to 
C                internal memory.  The LOCALGRANULEID is the file name 
C                in MODIS naming convention. Set_LclGranID_Sfx accepts 
C                an arbitrary file name suffix to characterize file 
C                type. See also design notes.
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
C  integer LUN_Of_LclVrsnID PCF RP LUN containing the value of LOCALVERSIONID.
C
C  character*(*) FileNameSuffix    
C                           String variable containing file suffix (e.g, 
C                           'dat', 'exe', 'hdf').  Suffix length up to 3 
C                           characters, starting at element 1 and terminated
C                           by first blank character. 
C
C
C !OUTPUT PARAMETERS:  None
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
C  The code design assumes that ECS attribute LOCALVERSIONID (used
C  to construct the LOCALGRANULEID) is implemented as a PCF Runtime
C  Parameter (RP) - thus the need to pass the RP LUN of the
C  LOCALVERSIONID.
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        Gen_Modis_FileName_V2     (atmos shared code)
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        PGS_MET_GetSetAttr_s      (libPGSTK.a)
C        PGS_PC_GetConfigData      (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        BLANK                     (Atmos_ECSMET.inc)
C        FAIL                      (Atmos_ECSMET.inc)
C        INVENTORYMETADATA         (Atmos_ECSMET.inc)
C        MCORE_SHORT_NAME          (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        PGSd_PC_VALUE_LENGTH_MAX  (PGS_PC.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Set_LclGranID_Sfx')


C Function argument declarations
      character*(*) MET_Handles(*), FileNameSuffix
      integer       LUN_Of_LclVrsnID


C other variable declarations
      character*8   ESDT_Name

      character*25  msg25
      character*128 FileName,ECS_AttrN
      character*512 msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) LclVrsnID

      integer    rtn, rtn_loc

      integer    pgs_met_setattr_s,pgs_met_getsetattr_s,
     1           PGS_PC_GetConfigData,
     2           Gen_Modis_FileName_V2,String_Loc

      integer    fbyte,fbyte_name,
     1           lbyte,lbyte_name

      logical    error_flag


C------------------------
C Initialization
C------------------------

      Set_LclGranID_Sfx =  FAIL
      FileName          =  BLANK
      error_flag        = .FALSE.

c-----set message buffers for error reporting
      write(msg25,'(I25)') LUN_Of_LclVrsnID
      rtn = string_loc(msg25,fbyte,lbyte)


C-----------------------------------------------------------------------
C Set "LocalGranuleID" in 4 steps.  The LocalGranuleID is the product
C file name in MODIS file naming convention.
C
C 1 - retrieve ESDT (Earth Science Data Type) SHORTNAME from the MCF.
C 2 - retrieve LOCALVERSIONID from PCF
C 3 - construct MODIS file name
C 4 - set LocalGranuleID
C-----------------------------------------------------------------------

      ECS_AttrN = MCORE_SHORT_NAME

      rtn = pgs_met_getsetattr_s(MET_Handles(INVENTORYMETADATA),ECS_AttrN,ESDT_Name)

c.....failed to read ESDT_Name
      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .TRUE.

         msgbuf =
     1   'pgs_met_getsetattr_s unable to read the ESDT shortname from product MCF.'
     2   // char(10) // 'Operator Action:  Check for correct MCF.  If wrong, stage '
     3   // char(10) // 'correct MCF and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c.....Oh dear, name is blank
      Else If (ESDT_Name .eq. BLANK) Then
         error_flag = .TRUE.

         msgbuf =
     1   'ESDT_Name retrieve from product MCF is blank.'
     2   // char(10) // 'Operator Action:  Check for correct MCF.  If wrong, stage '
     3   // char(10) // 'correct MCF and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c.....Retrieve value of LOCALVERSIONID from PCF runtime parameter.
      Else

         rtn = pgs_pc_getconfigdata(LUN_Of_LclVrsnID,LclVrsnID)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.

            msgbuf =
     1      'pgs_pc_getconfigdata unable to read LOCALVERSIONID on input argument '
     2      // char(10) // '"LUN_Of_LclVrsnID" (=' // msg25(fbyte:lbyte) // ').'
     3      // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     4      // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     5      // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Else

c...........Construct LOCALGRANULEID
c...........note: function Gen_Modis_FileName_V2 does its own input argument checking
            rtn = Gen_Modis_FileName_V2( ESDT_Name,     
     1                                   LclVrsnID, 
     2                                   FileNameSuffix, 
     3                                   FileName )

            If (rtn.ne.SUCCEED) Then
               error_flag = .TRUE.
               rtn_loc    =  string_loc(ESDT_Name,fbyte_name,lbyte_name)

               msgbuf =
     1         'Gen_Modis_FileName_V2 detected error constructing MODIS file name '
     2         // char(10) // 'for ESDT ' // ESDT_Name(fbyte_name:lbyte_name) // '.'
     3         // char(10) // 'Operator Action:  Refer to prior low level LogStatus error '
     4         // char(10) // 'messages originating in function Gen_Modis_FileName_V2.'


               call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else

c .............Set LOCALGRANULEID
               ECS_AttrN = MCORE_LOCALGRANULEID

               rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),ECS_AttrN,FileName)

               If (rtn .ne. PGS_S_SUCCESS) Then
                  error_flag = .TRUE.
                  rtn = String_Loc(ECS_AttrN,fbyte,lbyte)

                  msgbuf = 'PGS_MET_SetAttr_s unable to set '// ECS_AttrN(fbyte:lbyte)
     1             // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     2             // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     3             // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'



                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               End If   ! set LOCALGRANULEID

            End If   ! retrieve MODIS file name

         End If   ! retrieve LOCALVERSIONID from PCF

      End If   ! retrieve SHORTNAME from MCF


      If (.not.error_flag) Set_LclGranID_Sfx = SUCCEED

      Return

      End  !Set_LclGranID_Sfx



c---------------------------------------------------------------------------------
      Integer Function Set_LclInputGranID( MET_Handles,
     1                                     NUM_LclInputPntr,
     2                                     LUN_LclInputPntr,
     3                                     VRSN_LclInputPntr, NDVIfile,
     4                                     thresfile)

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Identify a list of MODIS and ancillary input files by 
C               their standard MODIS names and for ancillary products by 
C               the local names assigned at the processing center.  The 
C               file name list is to be stored in ECS Archive metadata 
C               record written to the HDF product under the attribute 
C               name LOCALINPUTGRANULEID.
C
C
C !INPUT PARAMETERS:
C  character MET_Handles(*)
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
C  integer NUM_LclInputPntr
C                Variable containing the number of input file names 
C                to be set.  If NUM_LclInputPntr <= 0, Set_LclInputGranID
C                returns immediately and does not set the LOCALINPUTGRANULEID.
C
C  integer LUN_LclInputPntr(*)
C                Array of input product LUNs  
C
C  integer VRSN_LclInputPntr(*)
C                Array of input product file version numbers 
C
C                A one-to-one correspondence between the elements of
C                the arrays LUN_LclInputPntr and
C                VRSN_LclInputPntr is assumed.
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
C !DESIGN NOTES:  Except for L1 MODIS products, MODIS and ancillary 
C                 product file names are retrieved from the PCF by 
C                 reading field 2 of the line entry matching the LUN and 
C                 Version identifiers of the file.  MODIS L1 product 
C                 file names are read from the LOCALGRANULEID attribute
C                 stored in the ECS Inventory metadata record of every 
C                 MODIS standard data product 
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        Retrieve_LclGranID_Set    (atmos shared code)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        ARCHIVEMETADATA           (Atmos_ECSMET.inc)
C        FAIL                      (Atmos_ECSMET.inc)
C        MAX_NUM_PNTRS             (Atmos_ECSMET.inc)
C        MCORE_LOCALINPUTGRANULEID (mapi.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSd_MET_STR_END          (PGS_MET.f via mapi.inc)
C        PGSd_PC_FILE_NAME_MAX     (PGS_PC.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

c-----PARAMETER declarations
      character*(*)   FUNCNAME
      parameter      (FUNCNAME = 'Set_LclInputGranID')

c-----Declaration of function arguments
      character*(*) MET_Handles(*)

      integer NUM_LclInputPntr
      integer LUN_LclInputPntr(*)
      integer VRSN_LclInputPntr(*)
      character*34 ndvifile
      character*23 thresfile


c-----Declaration of local variables
      character*25   msg25
      character*60   AttrN
      character*512  msgbuf
      character*(PGSd_PC_FILE_NAME_MAX) MODIS_FileNames(MAX_NUM_PNTRS)

      integer fbyte, lbyte, rtn, rtn_flag
      integer i
      integer NUM_MODIS_FileNames,
     1        NUM_MODIS_InputPntrs

      integer pgs_met_setattr_s, string_loc

      logical error_flag


C------------------------
C Initialization
C------------------------

      error_flag           = .false.
      NUM_MODIS_InputPntrs = NUM_LclInputPntr
      Set_LclInputGranID   = FAIL

C-----------------------------------------------------------------------
C Perform input argument checks 
C-----------------------------------------------------------------------

c-----if the number of input files is =< 0, return immediately 
      If (NUM_LclInputPntr .le. 0) Then
         Set_LclInputGranID = SUCCEED
         Return
      End If
 
c-----check if number of MODIS input files exceeds Set_ArchiveMetadata internal limit.
      If (NUM_LclInputPntr .gt. MAX_NUM_PNTRS) Then
         NUM_MODIS_InputPntrs = MAX_NUM_PNTRS
         error_flag           = .true.

         write(msg25,'(i25)') MAX_NUM_PNTRS
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf =
     1   'Number of MODIS input files exceeds internal code limit. Only '
     2   // char(10) // msg25(fbyte:lbyte) // ' elements of LOCALINPUTGRANULEID array to be set. '
     3   // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


c-----retrieve LocalInputGranID set
      Call Retrieve_LclGranID_Set( NUM_MODIS_InputPntrs,
     1                             LUN_LclInputPntr,
     2                             VRSN_LclInputPntr,
     3                             NUM_MODIS_FileNames,
     4                             MODIS_FileNames,
     5                             rtn_flag )

      If (rtn_flag .ne. SUCCEED) Then
         error_flag = .TRUE.

         msgbuf =
     1   'Retrieve_LclGranID_Set detected error reading local name of input product.' 
     2   // char(10) // 'Operator Action:  Refer to prior low level LogStatus messages'
     3   // char(10) // 'originating from call to routine Retrieve_LclGranID_Set.'

         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      End If

c ----  put in thresholds file name for position 6
      MODIS_FileNames(6) = thresfile
c ---- add in NDVI filename as used from LEOCAT (if sufficient space
      If (NUM_MODIS_FileNames .lt. MAX_NUM_PNTRS) then
         NUM_MODIS_FileNames = NUM_MODIS_FileNames + 1
         MODIS_FileNames(NUM_MODIS_FileNames) = ndvifile
      end if


c-----If sufficient space in MODIS_FileNames work buffer, place end-of-data
c-----marker into element NUM_MODIS_FileNames+1.

      If (NUM_MODIS_FileNames .lt. MAX_NUM_PNTRS)
     1   MODIS_FileNames(NUM_MODIS_FileNames+1) = PGSd_MET_STR_END

      AttrN = MCORE_LOCALINPUTGRANULEID
      rtn   = PGS_MET_SetAttr_s(MET_Handles(ARCHIVEMETADATA),AttrN,MODIS_FileNames)

      If (rtn.ne.PGS_S_SUCCESS) Then
          error_flag = .true.
          rtn = String_Loc(AttrN,fbyte,lbyte)

          msgbuf =
     1    'PGS_MET_SetAttr_s unable to set ' // AttrN(fbyte:lbyte) // '.'
     2    // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3    // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4    // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

          Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


      If (.not.error_flag) Set_LclInputGranID = SUCCEED

      Return

      End  !Set_LclInputGranID



c---------------------------------------------------------------------------------
      Integer Function Set_MeasParm_Atmos( MET_Handles,
     1                                     NUM_NewMeasParm,
     2                                     Name_NewMeasParm,
     3                                     PctMissing_NewMeasParm,
     4                                     AutoFlag_NewMeasParm,
     5                                     AutoFlagExp_NewMeasParm )


      implicit none

      include 'Atmos_ECSMET.inc'
      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Set ECS MeasuredParameter groups computed and passed
C                by current process to internal memory.
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
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_GetPCAttr_*       (libPGSTK.a)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        BLANK                   (Atmos_ECSMET.inc)
C        FAIL                    (Atmos_ECSMET.inc)
C        FUNCNAME_PGS_MET_SET_i  (Atmos_ECSMET.inc)
C        FUNCNAME_PGS_MET_SET_s  (Atmos_ECSMET.inc)
C        INVENTORYMETADATA       (Atmos_ECSMET.inc)
C        MCORE_*                 (mapi.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS           (PGS_SMF.f)
C        SUCCEED                 (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------

c PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'Set_MeasParm_Atmos')

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

      integer i1,i2,rtn,rtn_loc
      integer AttrV_i
      integer Class_Index
      integer pgs_met_setattr_s, pgs_met_setattr_i,String_Loc
      integer fbytea,fbyteb,fbytec,fbyted,fbytee,
     1        lbytea,lbyteb,lbytec,lbyted,lbytee

      logical error_flag, error_flag_group


C------------------------
C Initialization
C------------------------

      Set_MeasParm_Atmos =  FAIL
      error_flag        = .FALSE.


c-------------------------------------------------------------------------------
c Loop on list of new Measured Parameters
c-------------------------------------------------------------------------------

      Do 100 i1 = 1, NUM_NewMeasParm
         Name    = Name_NewMeasParm(i1)
         rtn_loc = string_loc(Name,fbytea,lbytea)
         Class_Index = i1


c--------invalid Measure Parameter Name
         If (Name .eq. BLANK) Then
            error_flag = .TRUE.

            msgbuf =
     1      'MeasuredParameter group name is blank - cannot be set. '
     2      // char(10) // 'Operator Action:  Notify SDST.'

           Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------valid Measure Parameter Name
         Else
            write(Class, '(I25)') Class_Index
            rtn_loc = string_loc(Class, fbyteb, lbyteb)


c-----------loop over 4 MeasuredParameter ODL objects used by MODIS atmosphere
            error_flag_group = .FALSE.

            Do 200 i2 = 1, 4

               If (i2 .eq. 1) Then   ! is PARAMETERNAME
                  AttrN = MCORE_PARAMETERNAME // '.' // Class(fbyteb:lbyteb)
                  AttrV_s = Name_NewMeasParm(i1)
                  rtn = PGS_MET_SetAttr_s(MET_Handles(INVENTORYMETADATA),AttrN,AttrV_s)


               Else If (i2 .eq. 2) Then   ! is QAPERCENTMISSINGDATA
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


               Else If (i2 .eq. 3) Then   ! is AUTOMATICQUALITYFLAG
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
     2                  // char(10) // 'is blank - not one of valid values:  Passed/Failed/Suspect. '
     3                  // 'To be set anyway. '
     4                  // char(10) // 'Operator Action:  Notify SDST.'
                     Else

                        msgbuf =
     1                  AttrN(fbytee:lbytee) // ' in MeasuredParameter group ' // Name(fbytea:lbytea)
     2                  // char(10) // 'is:  ' // AttrV_s(fbyted:lbyted)
     3                  // char(10) // '- not one of valid values:  Passed/Failed/Suspect. '
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
                     error_flag = .TRUE.

                     msgbuf =
     1               AttrN(fbytee:lbytee) // ' in MeasuredParameter group ' // Name(fbytea:lbytea)
     2               // char(10) // 'is blank - not a valid value.  To be set anyway. '
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
                  If (i2 .eq. 2) msgbuf = FUNCNAME_PGS_MET_SET_i // msgbuf_SET_err

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               End If

200         Continue

c rhucek 09/02/03: commented success message block
c-----------report success message to LogStatus
c           If (.not. error_flag_group) Then
c               msgbuf =
c    1          'MeasuredParameter group ' // Name(fbytea:lbytea) // ' successfully set. '
c
c               Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
c           End If

         End If

100   Continue


c-----return SUCCEED = 0 if no errors found
      If (.not.error_flag) Set_MeasParm_Atmos = SUCCEED

      Return

      END  !Set_MeasParm_Atmos



c---------------------------------------------------------------------------------
      integer function Set_PSA(MET_Handles,PSA_Name,PSA_Class,PSA_Value)

      Implicit None

C-----Insert Include files
      Include 'Atmos_ECSMET.inc'
      Include 'mapi.inc'
      Include 'PGS_SMF.f'
      Include 'PGS_MODIS_39500.f'

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
C  Externals:
C
C    Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        string_loc                (atmos shared code)
C
C    Named Constant:
C        FAIL                    (Atmos_ECSMET.inc)
C        INVENTORYMETADATA       (Atmos_ECSMET.inc)
C        MCORE_*                 (mapi.inc)
C        MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C        PGS_S_SUCCESS           (PGS_SMF.f)
C        SUCCEED                 (Atmos_ECSMET.inc)
C
C !END
C----------------------------------------------------------------------


C-----Parameter declarations
      Character*(*)  FUNCNAME
      Parameter     (FUNCNAME = 'Set_PSA')


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

      END  !Set_PSA



c---------------------------------------------------------------------------------
      Integer Function Set_RP_Data_Atmos( MET_Handles,
     1                                    LUN_Of_NUM_Of_RP_Pairs,
     2                                    MET_MasterGroup_ID )

      implicit none

      include 'Atmos_ECSMET.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Read and set a list of ECS metadata name/value pairs
C                tabulated as Runtime Parameters (RPs) in the PCF.  (See
C                design notes.)
C
C
C !INPUT PARAMETERS:
C
C  character*(*) MET_Handles(*)
C                    Array containing the names of the "MASTER"
C                    groups as defined in the MCF.  There may be
C                    up to 20 "MASTER" groups in the MCF.
C
C                    MET_Handles is defined during initialization
C                    of the ECS Metadata Configuration File (MCF).
C                    In FORTRAN, element 1 of the MET_Handles
C                    array refers to the MCF file.  ECS inventory
C                    metadata are referenced as element 2 and
C                    archive metadata as element 3.
C
C  integer LUN_Of_Num_Of_RP
C                    PCF LUN that identifies the head of the RP
C                    metadata list to be set.  The value on
C                    LUN_Of_Num_Of_RP is the number of Name/Value
C                    pairs.  The metadata names and values follow
C                    immediately on successive LUNs incrementing by 1.
C
C  integer MET_MasterGroup_ID
C                    Array index referencing MCF MASTER group element
C                    to which RP metadata attributes are to be set.
C                    2 -> Inventory; 3 -> Archive. (See MET_Handles
C                    above.)
C
C
C !OUTPUT PARAMETERS:  None
C
C
C !REVISION HISTORY:
C  rhucek 11/01/04
C  Added logic to return immediately in the case when input argument
C  NUM_Of_RP_Pairs is 0.  A notice message highlighting this incident
C  is written to the LogStatus file. 
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
C    A configured list of alternating metadata names and values on successive
C    RP LUNs is used.  The list is headed by the number of valid name/value
C    pairs in the sequence.
C
C  Returns:     SUCCEED (0) if successful, FAIL (-1) if error
C
C  Externals:
C
C    Functions and Subroutines:
C        Get_RP_Int                (atmos shared code)
C        Get_RP_String             (atmos shared code)
C        MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C        PGS_MET_SetAttr_s         (libPGSTK.a)
C        String_Loc                (atmos shared code)
C
C    Named Constant:
C        FAIL                      (Atmos_ECSMET.inc)
C        MAX_NUM_RP_PAIRS          (Atmos_ECSMET.inc)
C        MODIS_E_GENERIC           (PGS_MODIS_39500.f)
C        PGSd_PC_VALUE_LENGTH_MAX  (PGS_PC.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C        SUCCEED                   (Atmos_ECSMET.inc)
C
C !END
C-----------------------------------------------------------------------


c-----PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Set_RP_Data_Atmos')


c-----Declaration of function arguments
      character*(*) MET_Handles(*)
      integer       LUN_Of_NUM_Of_RP_Pairs, MET_MasterGroup_ID


c-----Declaration of local variables
      character*25   msg25, msg25a
      character*1012 msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) RP_Name, RP_Value

      integer LUN_N,LUN_V,NUM_Of_RP_Pairs
      integer pgs_met_setattr_s, string_loc

      integer fbyte,fbytea,fbyteb,fbytec,
     1        lbyte,lbytea,lbyteb,lbytec,
     2        i,rtn,rtn_flag,rtn_flag_RP_Name, rtn_flag_RP_Value

      logical error_flag


C------------------------
C Initialization
C------------------------

      Set_RP_Data_Atmos = FAIL
      error_flag = .false.

      write(msg25,'(I25)') LUN_Of_NUM_Of_RP_Pairs
      rtn = string_loc(msg25,fbyte,lbyte)


c-----check for valid inputs
      If (LUN_Of_NUM_Of_RP_Pairs .LE. 0) Then
         error_flag = .true.

         msgbuf =
     1   'Input PCF LUN ' // msg25(fbyte:lbyte) // ' is out of bounds.'
     2   // char(10) // 'Operator Action:  Check for valid PCF file.  If wrong or '
     3   // char(10) // 'corrupted, stage correct PCF and rerun PGE.  Otherwise, '
     4   // char(10) // 'notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


      If (MET_MasterGroup_ID .LE. 0) Then
         write(msg25a,'(I25)') MET_MasterGroup_ID
         rtn = string_loc(msg25a,fbytea,lbytea)

         error_flag = .true.

         msgbuf =
     1   'Input MET_MasterGroup_ID ' // msg25a(fbytea:lbytea) // ' is out of bounds.'
     2   // char(10) // 'Operator Action: notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf

c-----Return if invalid inputs found
      If (error_flag) Return


      Call Get_RP_Int(LUN_Of_NUM_Of_RP_Pairs, NUM_Of_RP_Pairs, rtn_flag)

      If (rtn_flag .eq. FAIL) Then
         error_flag = .true.

         msgbuf =
     1   'Get_RP_Int returned error reading integer RP on LUN ' // msg25(fbyte:lbyte)
     2   // char(10) // 'Operator Action:  Refer to prior low level error LogStatus '
     3   // char(10) // 'messages generated by routine Get_RP_Int.'

         Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      Else
         write(msg25a,'(I25)') NUM_Of_RP_Pairs
         rtn = string_loc(msg25a,fbytea,lbytea)

c--------return if no RP listed 
         If (NUM_Of_RP_Pairs .eq. 0) Then

            msgbuf = 
     1      'The number of RP pairs on LUN ' // msg25(fbyte:lbyte) // ' is 0.' 
     2      // char(10) // 'No RP pairs set in this call.'  

            Call modis_smf_setdynamicmsg(MODIS_N_GENERIC,msgbuf,FUNCNAME)

c--------test for out of bounds NUM_Of_RP_Pairs
         Else If (NUM_Of_RP_Pairs .lt. 0) Then 
            error_flag = .true.

            msgbuf =
     1      'NUM_Of_RP_Pairs on LUN = ' // msg25(fbyte:lbyte) // ' ('
     2      // msg25a(fbytea:lbytea) // ') is less than or equal to 0.'
     3      // char(10) // 'No metadata RP name/value pairs set by ' // FUNCNAME
     4      // char(10) // 'Operator Action:  Check for valid PCF file.  If wrong or '
     5      // char(10) // 'corrupted, stage correct PCF and rerun PGE.  Otherwise, '
     6      // char(10) // 'notify SDST.'

            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

c--------test for too many RPs
         Else If (NUM_Of_RP_Pairs .gt. MAX_NUM_RP_PAIRS) Then
            error_flag = .true.

            msgbuf =
     1      'Number of metadata RP name/value pairs '
     2      // 'exceeds ' // FUNCNAME // ' internal limit.'
     3      // char(10) // 'RP string on LUN ' // msg25(fbyte:lbyte) // ' has value of '
     4      // msg25a(fbytea:lbytea)
     5      // char(10) // 'No metadata RP name/value pairs set by ' // FUNCNAME
     6      // char(10) // 'Operator Action:  Check for valid PCF file.  If wrong or '
     7      // char(10) // 'corrupted, stage correct PCF and rerun PGE.  Otherwise, '
     8      // char(10) // 'notify SDST.'

            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c--------Valid number of RPs...Set RP pairs
         Else

c-----------loop over and set RP name/value pairs.
            DO 100 i = 2, NUM_Of_RP_Pairs * 2, 2
               LUN_N = LUN_Of_NUM_Of_RP_Pairs + (i - 1)
               LUN_V = LUN_Of_NUM_Of_RP_Pairs + i

               Call Get_RP_String(LUN_N, RP_Name, rtn_flag_RP_Name)

               If (rtn_flag_RP_Name .ne. FAIL) Then
                  Call Get_RP_String(LUN_V, RP_Value, rtn_flag_RP_Value)

c-----------------Set ECS metadata objects that are specified as name/value pairs in
c-----------------the USER DEFINED RUNTIME PARAMETERS section of the PCF.
                  If (rtn_flag_RP_Value .ne. FAIL) Then
                     rtn = pgs_met_setattr_s(MET_Handles(MET_MasterGroup_ID),RP_Name,RP_Value)

c--------------------set error flag if unable to set metadata name/value pair
                     If (rtn .ne. PGS_S_SUCCESS) Then
                        error_flag = .true.

                        rtn = String_Loc(RP_Name,fbyteb,lbyteb)
                        rtn = String_Loc(RP_Value,fbytec,lbytec)

                        msgbuf =
     1                  'pgs_met_setattr_s unable to set metadata RP name/value '
     2                  // char(10) // 'RP name is '// RP_Name(fbyteb:lbyteb)
     3                  // char(10) // 'RP value is ' // RP_Value(fbytec:lbytec)
     4                  // char(10) // 'Operator Action:  Check for valid MCF file and correct PCF '
     5                  // char(10) // 'reference to MCF file.  If MCF/PCF files are wrong or '
     6                  // char(10) // 'corrupted, stage correct files and rerun PGE.  Otherwise, '
     7                  // char(10) // 'notify SDST.'


                        Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

                     End If !---pgs_met_setattr_s check

                  End If !---return of Get_RP_String: RP_Value

               End If !---return of Get_RP_String: RP_Name

100         Continue

         End If !---test for 0 Runtime Parameters 

      End If   !---check on read RP on LUN_Of_NUM_RP_Pairs


      If (.not.error_flag) Set_RP_Data_Atmos = SUCCEED

      Return

      End  !Set_RP_Data_Atmos
