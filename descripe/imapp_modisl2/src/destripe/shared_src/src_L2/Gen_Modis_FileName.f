      integer function Gen_Modis_FileName(ESDT_Name,Local_VRSN_ID,fname)

      implicit none

      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !Description:  This function generates a file name in the MODIS V2
C                naming convention.
C
C
C !Input Parameters:
C       character*(*) ESDT_Name         ESDT shortname of the product
C
C       character*(*) Local_VRSN_ID     MODIS 3-digit string
C                                       representation of ECS metadata
C                                       field LOCALVERSIONID
C
C
C !Output Parameters:
C
C       character*(*) fname             version 2 product filename
C                                       (LOCALGRANULEID)
C
C !Revision History:
c Revision 1.6  2009/06/12 gfireman
c Add logic to insert "NRT" in fname if LUN 800504 starts with "N".
c
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
C         pgs_pc_getconfigdata          (libPGSTK.a)
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
C----------------------------------------------------------------------------


C PARAMETER declarations
      character*(*)  FUNCNAME
      parameter    ( FUNCNAME = 'Gen_Modis_FileName' )

      character*1    BLANK
      parameter    ( BLANK = ' ' )

      integer        LUN_GEO
      parameter    ( LUN_GEO = 600000 )

      integer        VRSN_GEO
      parameter    ( VRSN_GEO = 1 )

      integer        FAIL,      SUCCEED
      parameter    ( FAIL = -1, SUCCEED = 0 )

      integer        MAX_ESDT_Name_LEN
      parameter    ( MAX_ESDT_Name_LEN = 8 )

      integer        MAX_MODIS_FNAME_LEN
      parameter    ( MAX_MODIS_FNAME_LEN = 44 )

      integer        LUN_PROCTYPE
      parameter    ( LUN_PROCTYPE = 800504 )

C Declaration of function calling arguments
      character*(*) Local_VRSN_ID
      character*(*) ESDT_Name
      character*(*) fname


C Declaration of local variables
      character date*(10),time*(10),work_buf*(28),collectn_b*(28),datetime_b*(27)
      character*3   version
      character*25  msg25,msg25a, msg25b, msg25c
      character*512 msgbuf
      character*128 AttrN, buf_char(2), proctype

      integer buf_len,i,rtn,sl,string_len
      integer fbyte,fbytea,fbyteb,fbytec
      integer lbyte,lbytea,lbyteb,lbytec
      integer ESDT_Name_LEN,GMT_Delta

      double precision istart

      logical error_flag,VRSN_error_flag

C Function declarations
      integer pgs_td_asciitime_atob,pgs_td_utctotai,pgs_td_taitoutc,
     &        PGS_MET_GetPCAttr_s
      integer pgs_pc_getconfigdata
      integer strlen,string_loc


*----------------------------------------------------------------------------*


      error_flag=.FALSE.
      VRSN_error_flag=.FALSE.
      Gen_MODIS_FileName=FAIL


C-----------------------------------------------------------------------
C Perform initial input argument checks and return if not valid.
C
C 1 - ESDT name must not be blank and may contain no more
C     than 8 characters
C 2 - Local_VRSN_ID must be non-blank, and represent a positive integer
C     from 1 to 999.
C 3 - String buffer fname must be at least 44 characters in size.
C-----------------------------------------------------------------------

c.....Tests for blank ESDT_Name
      If (ESDT_Name .eq. BLANK) Then
         error_flag = .true.

         msgbuf = 'ESDT name is blank.'
     1   //  char(10) // 'Gen_MODIS_FileName unable to construct MODIS file name '
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
     1   // char(10) // 'Gen_MODIS_FileName unable to construct MODIS file name '
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
     2         // char(10) // 'Gen_Modis_FileName unable to construct MODIS file name.'
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


C-----------------------------------------------------------------------
C Return if input arguments are invalid
C-----------------------------------------------------------------------

      If (error_flag) Return


C-----------------------------------------------------------------------
C Retrieve values of RangeBeginningDate and RangeBeginningTime ECS
C metadata fields C from Geolocation Product.  Combine into CCSDS time
C format A and convert time format A to B.
C-----------------------------------------------------------------------

c.....setup for error reporting
      write(msg25b,'(I25)') LUN_GEO
      rtn = string_loc(msg25b,fbyteb,lbyteb)

      write(msg25c,'(I25)') VRSN_GEO
      rtn = string_loc(msg25c,fbytec,lbytec)

c.....Retrieve RangeBeginningDate
      AttrN = MCORE_RANGE_BEG_DATE
      rtn = PGS_MET_GetPCAttr_s(LUN_GEO,VRSN_GEO,MECS_CORE,AttrN,buf_char(1))

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
         rtn = PGS_MET_GetPCAttr_s(LUN_GEO,VRSN_GEO,MECS_CORE,AttrN,buf_char(2))

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


C-----------------------------------------------------------------------
C Construct text string containing file name in MODIS convention
C-----------------------------------------------------------------------

      if(.not.error_flag) then

C .......Calculate the actual string length of short name.

         sl=strlen(ESDT_Name)


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

C     Add suffix.  If ReprocessingActual begins with "N", add "NRT" (Near Real Time).
         rtn = pgs_pc_getconfigdata(LUN_PROCTYPE,proctype)
         if ( (rtn .eq. 0) .and. (proctype(1:1) .eq. 'N')) then
            fname(sl+33:) = '.NRT.hdf' // char(0)
         else
            fname(sl+33:) = '.hdf' // char(0)
         end if

      End If

      If (.not. error_flag)  Gen_Modis_FileName = SUCCEED

      Return

      End
