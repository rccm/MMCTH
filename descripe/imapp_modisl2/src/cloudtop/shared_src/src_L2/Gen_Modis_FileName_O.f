      integer function Gen_Modis_FileName_O(ESDT_Name,Local_VRSN_ID,fname)

      implicit none

      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !Description:  This function generates a file name in the MODIS V2 
C                naming convention.  It relies on the PCF 
C                CollectionStartTime to set the start time of the granule. 
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
c Revision 1.3  1998/10/28  13:00:33  rhucek
c Changed line from "istart = istart + GMT_Delta*60" to "istart = istart - GMT_Delta*60"
c
c Revision 1.2  1998/07/22  19:08:52  rhucek
c Renamed function Gen_Modis_FileName_O (was Gen_Modis_FileName)
c
c Revision 1.1  1998/07/20  21:19:00  rhucek
c Initial revision
c
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
C        pgs_td_asciitime_atob         (libPGSTK.a)
C        pgs_pc_getconfigdata          (libPGSTK.a)
C        pgs_td_utctotai               (libPGSTK.a)
C        pgs_td_taitoutc               (libPGSTK.a)
C
C
C       Named Constant:
C
C        MODIS_E_GENERIC               PGS_MODIS_39500.f 
C        PGS_S_SUCCESS                 PGS_SMF.f
C
C    Internals:
C
C       Functions and Subroutines:
C        MODIS_SMF_SETDYNAMICMSG
C        strlen
C        string_loc
C        Get_date_time
C
C       Variables:
C
C !END
C----------------------------------------------------------------------------


C PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Gen_Modis_FileName_O')

      character*1    BLANK
      parameter     (BLANK = ' ')

      integer        LUN_COLLECTN_START
      parameter     (LUN_COLLECTN_START=10258)

      integer        FAIL,      SUCCEED
      parameter     (FAIL = -1, SUCCEED = 0)

      integer        MAX_ESDT_NAME_LEN
      parameter     (MAX_ESDT_NAME_LEN = 8)

      integer        MAX_MODIS_FNAME_LEN
      parameter     (MAX_MODIS_FNAME_LEN = 44)

C Declaration of function arguments
      character*(*) Local_VRSN_ID
      character*(*) ESDT_Name
      character*(*) fname


C Declaration of local variables
      character date*(10),time*(10),work_buf*(28),
     &          collectn_b*(28),datetime_b*(27)
      character*3 version
      character*25 msg25,msg25a
      character*(512) msgbuf

      integer pgs_pc_getconfigdata,pgs_td_asciitime_atob,
     &        pgs_td_utctotai,pgs_td_taitoutc,
     &        strlen,string_loc
      integer buf_len,i,rtn,sl,string_len
      integer fbyte,fbytea,lbyte,lbytea
      integer ESDT_Name_LEN,GMT_Delta

      double precision istart

      logical errflag,VRSN_errflag


*----------------------------------------------------------------------------*


      errflag=.FALSE.
      VRSN_errflag=.FALSE.
      Gen_Modis_FileName_O=FAIL


C-----------------------------------------------------------------------
C Perform initial input argument checks and return if not valid. 
C
C 1 - ESDT name must not be blank and may contain no more 
C     than 8 characters
C 2 - Local_VRSN_ID must be non-blank, and represent a positive integer  
C     from 1 to 999.
C 3 - String buffer fname must be at least 44 characters in size.
C-----------------------------------------------------------------------

c ... Tests on ESDT_Name
C-----------------------

      If (ESDT_Name .eq. BLANK) Then
         errflag = .true.

         msgbuf = 
     1   'ESDT name is blank.'
     2   //  char(10)
     3   // 'Gen_Modis_FileName_O unable to construct MODIS file name ' 
     4   //  char(10) // 'without valid ESDT.'
     5   //  char(10)//'Operator Action: Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


      Else 
c ... ESDT name is non-blank; check name length
         rtn = String_Loc(ESDT_Name,fbyte,lbyte)
         ESDT_Name_LEN = lbyte - fbyte + 1

         write(msg25a,'(I25)') ESDT_Name_LEN
         rtn = String_Loc(msg25a,fbytea,lbytea)
 
         If (ESDT_Name_LEN .gt. MAX_ESDT_NAME_LEN) Then
            errflag = .true.

            msgbuf = 
     1      'ESDT name contains ' // msg25a(fbytea:lbytea) // ' ' 
     2      // 'characters, ' 
     3      //  char(10)
     4      // 'too long to comply with the MODIS file '
     5      // 'naming convention.'
     6      //  char(10)//'Operator Action: Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

         End If
      End If 


C ... Tests on Local_VRSN_ID
C---------------------------

      If (Local_VRSN_ID .eq. BLANK) Then
         errflag = .true.

         msgbuf = 
     1   'Local_VRSN_ID is blank.'
     2   //  char(10)
     3   // 'Gen_Modis_FileName_O unable to construct MODIS file name ' 
     4   //  char(10) // 'without valid Local version ID.'
     5   //  char(10)//'Operator Action: Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


      Else 
c ... Local_VRSN_ID is non-blank; check for valid range 
         rtn = String_Loc(Local_VRSN_ID,fbyte,lbyte)
         string_len = lbyte - fbyte + 1
 
         If (string_len .GT. 3) Then
            errflag = .true.

            write(msg25a,'(I25)') string_len 
            rtn = String_Loc(msg25a,fbytea,lbytea)
 
            msgbuf = 
     1      'Local_VRSN_ID contains ' // msg25a(fbytea:lbytea) // ' ' 
     2      // 'characters, too long to ' 
     3      //  char(10) // 'comply with 3-digit MODIS version number convention.'
     4      //  char(10) //'Operator Action: Notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         Else 
            If (string_len .eq. 1) 
     1          version = '00' // Local_VRSN_ID(fbyte:lbyte)    
            If (string_len .eq. 2) 
     1          version = '0' // Local_VRSN_ID(fbyte:lbyte)    
            If (string_len .eq. 3) 
     1          version = Local_VRSN_ID(fbyte:lbyte)

            Do i=1,3
               If (version(i:i).lt.'0' .or. version(i:i).gt.'9') 
     1            VRSN_errflag = .true.            
            End Do

            If (VRSN_errflag) Then
               errflag = .true.
         
               msgbuf = 
     1         'Local_VRSN_ID "' // version // '" does not '
     2         // 'represent a 3-digit integer.'
     3         // char(10) // 'Gen_Modis_FileName_O unable to construct '
     4         // 'MODIS file name.'
     5         //  char(10)//'Operator Action: Notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

             End If  ! check on number of valid digits 
         End If   ! check for 3 digits or less 
      End If   ! check for blank version id  


c ... Test for adequate fname buffer length
C------------------------------------------

      buf_len = LEN(fname)

      If (buf_len .lt. MAX_MODIS_FNAME_LEN) Then
         errflag = .true.

         write(msg25,'(I25)') buf_len
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf =
     1   'string buffer "fname" contains only ' // msg25(fbyte:lbyte) 
     2   // ' characters, ' // char(10) // 'too small to contain '
     3   // 'MODIS file name.'
     4   //  char(10)//'Operator Action: Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      End If


C-----------------------------------------------------------------------
C Return if input arguments are invalid
C-----------------------------------------------------------------------

      If (errflag) Return


C-----------------------------------------------------------------------
C Read value of ECS collection start time from PCF and convert from
C CCSDS format A to B.
C-----------------------------------------------------------------------

C ... read value of ECS collection start time from PCF

      rtn = PGS_PC_GetConfigData(LUN_COLLECTN_START,work_buf)

      If (rtn .ne. PGS_S_SUCCESS) Then
         errflag = .true.
         write(msg25,'(I25)') LUN_COLLECTN_START
         rtn = string_loc(msg25,fbyte,lbyte)

         msgbuf =
     2   'PGS_PC_GetConfigData unable to read collection start'
     3   //char(10)//'time on LUN '
     4   //msg25(fbyte:lbyte)
     5   //char(10)//'Operator Action: Check for valid "Collection Start Time"'
     6   //' Runtime Parameter'//char(10)//'entry in PCF.  If incorrect, '
     7   //'stage correct PCF and rerun'//char(10)//'PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      Else

C .......Convert collection start time from CCSDS format A to format B.

         rtn=pgs_td_asciitime_atob(work_buf,collectn_b)

         If (rtn.ne.PGS_S_SUCCESS) Then
            errflag = .true.

            msgbuf = 'pgs_td_asciitime_atob unable to convert'
     2               //' collection start time '
     3               //char(10)
     4               //'from CCSDS format A to format B.'
     5               //  char(10)//'Operator Action: '
     6               //'Check system and SDPTK configuration. If fault identified,'
     7               //char(10)
     8               //'correct and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         End If
      End If


C-----------------------------------------------------------------------
C Query system for current local time and transform it to current GMT 
C time in CCSDS B format.
C-----------------------------------------------------------------------

      If (.not. errflag) Then

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
            errflag = .true.

            msgbuf = 'pgs_td_utctotai unable to convert current'
     2      // char(10) // 'local time from CCSDS A to TAI format.'
     3      // char(10) // 'Operator Action:  Check system resources/environment, '
     4      // char(10) // 'PCF, and SDPTK configuration (including up-to-date Leap '
     5      // char(10) // 'Seconds and UT1 Pole files).  If a fault is identified, '
     6      // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

         Else


C ..........Convert current locat time to GMT time
            istart = istart - GMT_Delta*60

C ..........Convert current GMT time from tai to CCSDS format A
            rtn=pgs_td_taitoutc(istart,work_buf)

            If (rtn.ne.PGS_S_SUCCESS) Then
               errflag = .true.

               msgbuf = 
     1         'pgs_td_taitoutc unable to convert current GMT '
     2         //char(10)
     3         //'time from TAI to CCSDS A format.'
     4         //  char(10)//'Operator Action: '
     5         //'Check system and SDPTK configuration. If fault identified,'
     6         //char(10)
     7         //'correct and rerun PGE.  Otherwise, notify SDST.'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

            Else

C .............Convert current GMT time from CCSDS format A to format B
               rtn=pgs_td_asciitime_atob(work_buf,datetime_b)


               If (rtn.ne.PGS_S_SUCCESS) Then
                  errflag = .true.

                  msgbuf = 
     1            'pgs_td_asciitime_atob unable to convert current '
     2            // char(10) // 'GMT time from CCSDS A to CCSDS '
     3            //'B format.'
     4            //  char(10)//'Operator Action: '
     5            //'Check system and SDPTK configuration. If fault identified,'
     6            //char(10)
     7            //'correct and rerun PGE.  Otherwise, notify SDST.'

                  Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf, FUNCNAME)

               End If
            End If
         End If
      End If


C-----------------------------------------------------------------------
C Construct text string containing file name in MODIS convention 
C-----------------------------------------------------------------------

      if(.not.errflag) then

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

         fname(sl+33:sl+36)='.hdf'
         fname(sl+37:sl+37)=char(0)
      End If

      If (.not. errflag)  Gen_Modis_FileName_O = SUCCEED

      Return
      End
