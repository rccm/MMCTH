      Subroutine Parse_ECS_DateTime( File_LUN, File_VRSN, 
     1                               HDF_AttributeName,
     2                               ECS_DateTimeGroup,
     3                               year, month, day, 
     4                               hour, minute, second, fsecond,
     5                               rtn_code )

C---------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: Parse an ECS "DateTime" metadata group into separate 
C              components rendering the year, month, day, hour, minute, 
C              second, and fractional second.
C
C
C!INPUT PARAMETERS:
C  integer File_LUN       PCF LUN number of file to be queried 
C  integer File_VRSN      PCF Version number of file to be queried 
C
C  character*(*) HDF_AttributeName
C                         HDF file attribute name containing ODL object
C                         ECS_DateTimeGroup (see below)
C
C  character*(*) ECS_DateTimeGroup 
C                         ECS DateTime metadata group to be parsed.  Case 
C                         sensitive valid values are:  'SingleDateTime', 
C                         'BeginningDateTime' and 'EndingDateTime'.
C
C
C!OUTPUT PARAMETERS:
C  integer year           year parsed from ECS DateTime group 
C  integer month          month parsed from ECS DateTime group 
C  integer day            day parsed from ECS DateTime group 
C  integer hour           hour parsed from ECS DateTime group 
C  integer minute         minute parsed from ECS DateTime group 
C  integer second         second parsed from ECS DateTime group 
C  integer fsecond        fractional seconds in units of microsecond parsed 
C                         from ECS DateTime group 
C  integer rtn_code       subroutine return value: SUCCEED (0) or FAIL (-1) 
C
C
C!REVISION HISTORY:
C
C
C!TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C!REFERENCES and CREDITS:
C    Written by
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C DESIGN NOTES:  Time data for all ECS DateTime groups is assumed to be 
C                stored in the format HH:MM:SS.dddddd.  In particular, 6 
C                digits of fractional seconds are assumed present. 
C
C
C  Returns:     SUCCEED( 0 )if successful, FAIL( -1) if an error occurs
C               or the granule time is not in the process interval.
C
C    Externals:
C       Functions and Subroutines:
C          pgs_met_getpcattr_s   (libPGSTK.a)
C
C       Named Constants:
C
C          MODIS_E_GENERIC          (PGS_MODIS_39500.f)
C          PGS_S_SUCCESS            (PGS_SMF.f)
C          MCORE_RANGE_BEG_DATE     (mapi.inc)
C          MCORE_RANGE_BEG_TIME     (mapi.inc)
C          MCORE_RANGE_ENDING_TIME  (mapi.inc)
C          MCORE_RANGE_ENDING_DATE  (mapi.inc)
C
C
C    Internals:
C       Functions and Subroutines:
C          CCSDS_FullTimeFormat_A    (atmos shared code)
C          MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C          string_loc                (atmos shared code)
C
C
C!END
C----------------------------------------------------------------------

      IMPLICIT NONE

      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'


C PARAMETER declarations
      character*(*)     BLANK
      PARAMETER        (BLANK = ' ')

      character*(*)     FUNCNAME
      PARAMETER        (FUNCNAME = 'Parse_ECS_DateTime')

      integer           FAIL,      SUCCEED
      PARAMETER        (FAIL = -1, SUCCEED = 0)


C Declaration of function arguments
      character*(*)  ECS_DateTimeGroup, HDF_AttributeName

      integer File_LUN,File_VRSN,
     1        year,month,day,hour,minute,second,fsecond,
     2        rtn_code


C Declaration of local variables
      character*15   time,date
      character*25   msg25_LUN,msg25_VRSN
      character*27   FullDateTime
      character*256  ECS_DateParm,ECS_TimeParm
      character*2048 msgbuf,hdfattrname

      integer ios,rtn
      integer fbyte,fbyte1,fbyte2,fbyte_a,fbyte_d,fbyte_t,fbyte_LUN,fbyte_VRSN,
     1        lbyte,lbyte1,lbyte2,lbyte_a,lbyte_d,lbyte_t,lbyte_LUN,lbyte_VRSN

      logical error_flag

C Function declarations 
      character*27 CCSDS_FullTimeFormat_A

      integer pgs_met_getpcattr_s, string_loc


C-----------------------------------------------------------------------
C Initialization
C-----------------------------------------------------------------------
      error_flag  = .false.
      hdfattrname = HDF_AttributeName
      rtn_code    = FAIL


C-----------------------------------------------------------------------
C Identify desired ECS DateTime group and select correct time and data 
C attributes (note use of MAPI parameter names.  Return fail if input 
C argument "ECS_DateTimeGroup" does not match one of valid values. 
C-----------------------------------------------------------------------
      If (ECS_DateTimeGroup .eq. 'SingleDateTime') Then

c........MAPI parameter not yet available for SingleDateTime attributes 
         ECS_TimeParm = 'TIMEOFDAY'
         ECS_DateParm = 'CALENDARDATE'

      ElseIf (ECS_DateTimeGroup .eq. 'BeginningDateTime') Then
         ECS_TimeParm = MCORE_RANGE_BEG_TIME
         ECS_DateParm = MCORE_RANGE_BEG_DATE 

      ElseIf (ECS_DateTimeGroup .eq. 'EndingDateTime') Then 
         ECS_TimeParm = MCORE_RANGE_ENDING_TIME
         ECS_DateParm = MCORE_RANGE_ENDING_DATE 

      Else
         rtn = String_Loc(ECS_DateTimeGroup,fbyte,lbyte)

         msgbuf =
     1   'Invalid input argument: ECS_DateTimeGroup = "' // ECS_DateTimeGroup(fbyte:lbyte) // '".'
     2   // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME) 

         Return
      EndIf


C-----------------------------------------------------------------------
C Define status message component terms
C-----------------------------------------------------------------------
      write(msg25_LUN, '(I25)') File_LUN
      write(msg25_VRSN,'(I25)') File_VRSN

      rtn = String_Loc(msg25_LUN, fbyte_LUN, lbyte_LUN)
      rtn = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)
      rtn = String_Loc(hdfattrname,fbyte_a,lbyte_a)
      rtn = String_Loc(ECS_TimeParm,fbyte_t,lbyte_t)
      rtn = String_Loc(ECS_DateParm,fbyte_d,lbyte_d)
          

C-----------------------------------------------------------------------
C Retrieve time and date values from stored ECS metadata.  Separate into 
C year, month, day, hour, minute, second and fractional second.  If 
C error, set error flag. 
C-----------------------------------------------------------------------

c-----Get ECS time attribute value
      rtn = pgs_met_getpcattr_s( File_LUN,
     1                           File_VRSN,
     2                           HDF_AttributeName,
     3                           ECS_TimeParm, 
     4                           time )

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf =
     1   'pgs_met_getpcattr_s unable to read ODL object ' 
     2   // char(10) // '"' // ECS_TimeParm(fbyte_t:lbyte_t) // '" '
     3   // 'in HDF attribute "' // hdfattrname(fbyte_a:lbyte_a) // '" '  
     4   // char(10) // 'in file on PCF LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     5   // ', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.'
     6   // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     7   // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     8   // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c-----Get granule end date
      Else
         rtn = pgs_met_getpcattr_s( File_LUN,
     1                              File_VRSN,
     2                              HDF_AttributeName,
     3                              ECS_DateParm,
     4                              date )

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .true.

            msgbuf =
     1      'pgs_met_getpcattr_s unable to read ODL object ' 
     2      // char(10) // '"' // ECS_DateParm(fbyte_d:lbyte_d) // '" '
     3      // 'in HDF attribute "' // hdfattrname(fbyte_a:lbyte_a) // '" ' 
     4      // char(10) // 'in file on PCF LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     5      // ', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.'
     6      // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     7      // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     8      // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c--------Extract time/date components
         Else
            rtn = String_Loc(date,fbyte1,lbyte1)
            rtn = String_Loc(time,fbyte2,lbyte2)
                      
            FullDateTime = CCSDS_FullTimeFormat_A(date(fbyte1:lbyte1)//'T'//time) 

            If (FullDateTime .eq. BLANK) Then
               error_flag = .true.

               msgbuf = 
     1         'CCSDS_FullTimeFormat_A unable to convert retrieved date '
     2         // char(10) // '("' // date(fbyte1:lbyte1) // '") and time ' 
     3         // '("' // time(fbyte2:lbyte2) // '") to full CCSDS code A format. '
     4         // char(10) // 'Operator Action:  Notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            Else
               read(FullDateTime(1:4),   '(I4)', IOSTAT=ios)  year
               If (ios .ne. 0) error_flag = .true.
   
               read(FullDateTime(6:7),   '(I2)', IOSTAT=ios)  month
               If (ios .ne. 0) error_flag = .true.
   
               read(FullDateTime(9:10),  '(I2)', IOSTAT=ios)  day
               If (ios .ne. 0) error_flag = .true.
   
               read(FullDateTime(12:13),   '(I2)', IOSTAT=ios)  hour
               If (ios .ne. 0) error_flag = .true.
   
               read(FullDateTime(15:16),   '(I2)', IOSTAT=ios)  minute
               If (ios .ne. 0) error_flag = .true.
   
               read(FullDateTime(18:19),   '(I2)', IOSTAT=ios)  second
               If (ios .ne. 0) error_flag = .true.
   
               read(FullDateTime(21:26), '(I6)', IOSTAT=ios)  fsecond
               If (ios .ne. 0) error_flag = .true.
   
               If (error_flag) Then
                  msgbuf =
     1            'One or more internal read errors converting CCSDS '
     2            // char(10) // 'full code A date/time string to numeric integers.'
     3            // char(10) // 'Operator Action:  Notify SDST.'
   
                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
               EndIf !If (error_flag)

            EndIf  !Test for BLANK datetime

         EndIf !If (rtn .ne. PGS_S_SUCCESS)

      EndIf !If (rtn.ne.PGS_S_SUCCESS)


      If (.not.error_flag) rtn_code = SUCCEED

      Return

      End
