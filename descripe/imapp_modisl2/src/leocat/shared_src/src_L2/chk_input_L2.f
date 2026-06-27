      integer function chk_input_L2 ( LUN_Input_Granule,
     1                                VRSN_Input_Granule )

      IMPLICIT NONE
      include 'PGS_SMF.f'
      include 'PGS_PC.f'
      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'

C---------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C
C    This function reads the temporal coverage start and stop time 
C    metadata of an ECS input granule, and verifies that (to within a 
C    delta epsilon) these times fall within the processing interval 
C    of the current output collection granule.
C
C
C!INPUT PARAMETERS:
C
C   integer            LUN_Input_Granule   PCF LUN number of input
C                                          granule
C   integer            VRSN_Input_Granule  PCF Version number of 
C                                          input granule
C
C
C!OUTPUT PARAMETERS: N/A
C
C
C!REVISION HISTORY:
c Revision 1.4  1998/12/02  22:55:55  rhucek
c Relaxed granule level temporal "match" criteria.
c
c Revision 1.3  1998/04/28  20:48:22  rhucek
c Updated error messages with "Operator Action:' strings
c
c Revision 1.2  1997/11/03  18:03:53  rhucek
c Updated descriptions of module "returns"
c
c Revision 1.1  1997/11/01  13:08:31  rhucek
c Initial revision
c
C
C
C!TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C!REFERENCES and CREDITS:
C
C
C    Written by Liqun Ma       October 1997
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C    lma@ltpmail.gsfc.nasa.gov
C
C
C DESIGN NOTES:
C
C    The code reads the input granule's ECS metadata to extract the
C    temporal coverage start and stop datetimes.  These are converted 
C    to seconds from Jan 1, 1993 and compared against the current 
C    granule processing interval, which is also expressed in seconds 
C    from Jan 1, 1993.  The processing interval is supplied by ECS 
C    through the PCF Runtime Parameter LUNs 10258 and 10259 which are 
C    read by the code.  For an input granule to be accepted for 
C    processing, its start and stop times must fit (after allowing for
C    a small error epsilon) within the data processing time window.  
C
C  Returns:     SUCCEED( 0 )if successful, FAIL( -1) if an error occurs
C               or the granule time is not in the process interval.
C
C    Externals:
C       Functions and Subroutines:
C          pgs_met_getpcattr_s   (libPGSTK.a)
C          pgs_td_utctotai       (libPGSTK.a)
C          PGS_PC_GetConfigData  (libPGSTK.a)
C
C       Named Constants:
C
C          MODIS_E_GENERIC          (PGS_MODIS_39500.f)
C          PGS_S_SUCCESS            (PGS_SMF.f)
C          PGSd_PC_*                (PGS_PC.f)
C          MECS_CORE                (mapi.inc)
C          MCORE_RANGE_BEG_TIME     (mapi.inc)
C          MCORE_RANGE_BEG_DATE     (mapi.inc)
C          MCORE_RANGE_ENDING_TIME  (mapi.inc)
C          MCORE_RANGE_ENDING_DATE  (mapi.inc)
C          
C
C    Internals:
C       Functions and Subroutines:
C          MODIS_SMF_SETDYNAMICMSG
C          String_Loc
C
C       Variables:
C
C
C!END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*)     FUNCNAME
      PARAMETER        (FUNCNAME = 'chk_input_L2')
  
c ... ECS Runtime Parameter LUNs for Collection Start & Stop date/times 
      integer           LUN_IS_10258,         LUN_IS_10259
      PARAMETER        (LUN_IS_10258 = 10258, LUN_IS_10259 = 10259)

      integer           FAIL,      SUCCEED 
      PARAMETER        (FAIL = -1, SUCCEED = 0)

      double precision  DELTA 
      PARAMETER        (DELTA = 1.5d0)


C Declaration of function arguments
      integer LUN_Input_Granule,VRSN_Input_Granule


C Declaration of local variables
      character*15   stime,sdate,etime,edate
      character*25   msg25, msg25_LUN, msg25_VRSN
      character*256  hdfattrname,parmname
      character*2048 msgbuf
      character*(PGSd_PC_VALUE_LENGTH_MAX) work_buf1,work_buf2
      character*(PGSd_PC_VALUE_LENGTH_MAX) utcatm

      integer pgs_td_utctotai,pgs_met_getpcattr_s
      integer PGS_PC_GetConfigData,String_Loc
      integer fbyte,fbytd,fbytt,fbyte_LUN,fbyte_VRSN,
     1        lbyte,lbytd,lbytt,lbyte_LUN,lbyte_VRSN,
     2        rtn

      double precision gra_start,gra_stop
      double precision collectn_start,collectn_stop

      logical error_flag



C Initialization
      chk_input_L2 = FAIL
      error_flag = .false.

c-----Identify hdf attribute containing ECS inventory metadata.
      hdfattrname = MECS_CORE


C-----------------------------------------------------------------------
C Perform input argument checks and return if not valid 
C-----------------------------------------------------------------------

c-----Set status message component terms 
      write(msg25_LUN, '(I25)') LUN_Input_Granule
      write(msg25_VRSN,'(I25)') VRSN_Input_Granule

      rtn = String_Loc(msg25_LUN, fbyte_LUN, lbyte_LUN )
      rtn = String_Loc(msg25_VRSN,fbyte_VRSN,lbyte_VRSN)


c-----Check for valid input granule PCF LUN
      If (LUN_Input_Granule .LE. 0) Then
         error_flag = .true.

         msgbuf = 
     1   'Function input argument "LUN_Input_Granule" (=' // msg25_LUN(fbyte_LUN:lbyte_LUN) 
     2   // ') is out of bounds.' 
     3   // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


c-----Check for valid input granule PCF version number 
      If (VRSN_Input_Granule .LE. 0) Then
         error_flag = .true.

         msgbuf = 
     1   'Function input argument "VRSN_Input_Granule" (=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) 
     2   // ') is out of bounds. ' 
     3   // char(10) // 'Operator Action:  Notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf

c-----Return if invalid inputs found
      If (error_flag) Return
         
 
C-----------------------------------------------------------------------
C Read the collection start and stop datetimes from PCF.
C Convert PCF ascii datetimes to seconds from Jan 1, 1993.
C Return if error incurred.
C-----------------------------------------------------------------------
 
c-----Read PCF collection start datetime 
      rtn = PGS_PC_GetConfigData(LUN_IS_10258,work_buf1)

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.

         write(msg25,'(I25)') LUN_IS_10258
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 
     1   'PGS_PC_GetConfigData unable to read Collection Start Time on '
     2   // 'PCF RP LUN='// msg25(fbyte:lbyte) // '.'
     3   // char(10) // 'Operator Action:  Check for valid "Collection Start Time" value.  If '
     4   // char(10) // 'incorrect, stage correct PCF and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c-----Read PCF collection stop datetime 
      Else
         rtn = PGS_PC_GetConfigData(LUN_IS_10259,work_buf2)

         If (rtn.ne.PGS_S_SUCCESS) Then
            error_flag = .true.

            write(msg25,'(I25)') LUN_IS_10259
            rtn = String_Loc(msg25,fbyte,lbyte)

            msgbuf = 
     1      'PGS_PC_GetConfigData unable to read "Collection Stop Time" on '
     2      // 'PCF RP LUN='// msg25(fbyte:lbyte) // '.'
     3      // char(10) // 'Operator Action:  Check for valid "Collection Stop Time" value.  If '
     4      // char(10) // 'incorrect, stage correct PCF and rerun.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c--------Convert ASCII start/stop datetimes to TAI seconds 
         Else
            rtn = pgs_td_utctotai(work_buf1,collectn_start)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               write(msg25,'(I25)') LUN_IS_10258
               rtn = String_Loc(msg25,fbyte,lbyte)

               msgbuf = 
     1         'pgs_td_utctotai failed to convert Collection Start Time (read from '
     2         // char(10) // 'PCF RP LUN=' // msg25(fbyte:lbyte) // ') to TAI. '
     3         // char(10) // 'Operator Action:  Check for valid PCF (including up-to-date '
     4         // char(10) // 'Leap Seconds and UT1 Pole files) and input granule.  If a fault ' 
     5         // char(10) // 'is identified, correct and rerun PGE.  Otherwise, notify SDST.' 

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            End If

            rtn = pgs_td_utctotai(work_buf2,collectn_stop)

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               write(msg25,'(I25)') LUN_IS_10259
               rtn = String_Loc(msg25,fbyte,lbyte)

               msgbuf = 
     1         'pgs_td_utctotai failed to convert Collection Stop Time (read from '
     2         // char(10) // 'PCF RP LUN=' // msg25(fbyte:lbyte) // ') to TAI. '
     3         // char(10) // 'Operator Action:  Check for valid PCF (including up-to-date '
     4         // char(10) // 'Leap Seconds and UT1 Pole files) and input granule.  If a fault ' 
     5         // char(10) // 'is identified, correct and rerun PGE.  Otherwise, notify SDST.' 

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            EndIf

         EndIf

      EndIf


c-----If error reading and converting collection start/stop date&times, return 
      If (error_flag) Return
 

C-----------------------------------------------------------------------
C Get input granule beginning time and date from stored ECS metadata.
C Convert ascii datetimes to seconds from Jan 1, 1993.
C Return if error incurred
C-----------------------------------------------------------------------

c-----Get granule begin time; use MAPI parameter name
      parmname = MCORE_RANGE_BEG_TIME

      rtn = pgs_met_getpcattr_s( LUN_Input_Granule,
     1                           VRSN_Input_Granule,
     2                           hdfattrname,
     3                           parmname,
     4                           stime )

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf = 
     1   'pgs_met_getpcattr_s unable to read ECS ' // MCORE_RANGE_BEG_TIME 
     2   // char(10) // 'from input product on PCF LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN) 
     3   // ', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.'
     4   // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     5   // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     6   // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c-----Get granule begin date.  Use MAPI parameter name
      Else
         parmname = MCORE_RANGE_BEG_DATE
   
         rtn = pgs_met_getpcattr_s( LUN_Input_Granule,
     1                              VRSN_Input_Granule,
     2                              hdfattrname,
     3                              parmname,
     4                              sdate )
   
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .true.
   
            msgbuf = 
     1      'pgs_met_getpcattr_s unable to read ECS ' // MCORE_RANGE_BEG_DATE
     2      // char(10) //'from input product on PCF LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN) 
     3      //', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.'
     4      // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     5      // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     6      // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c--------Convert ASCII start and stop datetimes to TAI seconds 
         Else 
            rtn = String_Loc(stime,fbytt,lbytt)
            rtn = String_Loc(sdate,fbytd,lbytd)

            utcatm(1:10)  = sdate(fbytd:lbytd)
            utcatm(11:11) = 'T'
            utcatm(12:26) = stime(fbytt:lbytt)
            utcatm(27:27) = 'Z'

c-----convert ASCII time to TAI seconds 
            rtn = pgs_td_utctotai( utcatm,
     1                             gra_start )

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf = 
     1         'pgs_td_utctotai failed to convert granule StartDateTime (read from '
     2         // char(10) // 'input file on PCF LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     3         // ', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // ') to TAI. '
     4         // char(10) // 'Operator Action:  Check for valid PCF (including up-to-date '
     5         // char(10) // 'Leap Seconds and UT1 Pole files) and input granule.  If a fault ' 
     6         // char(10) // 'is identified, correct and rerun PGE.  Otherwise, notify SDST.' 

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            EndIf 

         EndIf   ! check return read MCORE_RANGE_BEG_DATE

      EndIf   ! check return read MCORE_RANGE_BEG_TIME 


c-----If error reading and converting granule start datetime, return 
      If (error_flag) Return


C-----------------------------------------------------------------------
C Get input granule stop time and date from stored ECS metadata.
C Convert ascii datetimes to seconds from Jan 1, 1993.
C Return if error incurred
C-----------------------------------------------------------------------

c-----Get granule ending time; use MAPI parameter name
      parmname = MCORE_RANGE_ENDING_TIME

      rtn = pgs_met_getpcattr_s( LUN_Input_Granule,
     1                           VRSN_Input_Granule,
     2                           hdfattrname,
     3                           parmname,
     4                           etime )

      If (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf = 
     1   'pgs_met_getpcattr_s unable to read ECS ' // MCORE_RANGE_ENDING_TIME 
     2   // char(10) // 'from input product on PCF LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN) 
     3   // ', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.'
     4   // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     5   // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     6   // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c-----Get granule end date.  Use MAPI parameter name
      Else
         parmname = MCORE_RANGE_ENDING_DATE
   
         rtn = pgs_met_getpcattr_s( LUN_Input_Granule,
     1                              VRSN_Input_Granule,
     2                              hdfattrname,
     3                              parmname,
     4                              edate )
   
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .true.
   
            msgbuf = 
     1      'pgs_met_getpcattr_s unable to read ECS ' // MCORE_RANGE_ENDING_DATE 
     2      // char(10) //'from input product on PCF LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN) 
     3      //', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // '.'
     4      // char(10) // 'Operator Action:  Check for corrupted PCF and/or input '
     5      // char(10) // 'file.  If a fault is identified, stage correct PCF/input '
     6      // char(10) // 'file and rerun PGE.  Otherwise, notify SDST.'
   
            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c--------Successfully read granule begin datetime.  Now convert ASCII start and stop  
c--------datetimes to TAI seconds 
         Else 
            rtn = String_Loc(etime,fbytt,lbytt)
            rtn = String_Loc(edate,fbytd,lbytd)

            utcatm(1:10)  = edate(fbytd:lbytd)
            utcatm(11:11) = 'T'
            utcatm(12:26) = etime(fbytt:lbytt)
            utcatm(27:27) = 'Z'

c-----------convert ASCII time to TAI seconds 
            rtn = pgs_td_utctotai( utcatm,
     1                             gra_stop )

            If (rtn.ne.PGS_S_SUCCESS) Then
               error_flag = .true.

               msgbuf = 
     1         'pgs_td_utctotai failed to convert granule StopDateTime (read from '
     2         // char(10) // 'input file on PCF LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN)
     3         // ', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // ') to TAI. '
     3         // char(10) // 'Operator Action:  Check for valid PCF (including up-to-date Leap Seconds '
     4         // char(10) // 'and UT1 Pole files) and input granule RangeDateTime ECS metadata.  If a ' 
     5         // char(10) // 'fault is identified, correct and rerun PGE.  Otherwise, notify SDST.' 

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            EndIf 

         EndIf   ! check return read MCORE_RANGE_ENDING_DATE

      EndIf   ! check return read MCORE_RANGE_ENDING_TIME


c-----check consistency of granule start and stop times.
      If (gra_start .ge. gra_stop) Then
         error_flag = .true.

         msgbuf = 
     1   'Input granule (LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN) // ', VRSN=' 
     2   // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // ') StartTime greater or equal to StopTime. '
     3   // char(10) // 'Operator Action:  Check for valid PCF (including up-to-date Leap Seconds '
     4   // char(10) // 'and UT1 Pole files) and input granule RangeDateTime ECS metadata.  If a ' 
     5   // char(10) // 'fault is identified, correct and rerun PGE.  Otherwise, notify SDST.' 

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


c-----If error reading and converting granule start datetime, return 
      If (error_flag) Return
 

C-----------------------------------------------------------------------
c Compare the input granule start and stop times with the processing
c interval (collection) start and stop times.  The granule times 
C must fit (to within a "delta") within the processing time window 
C-----------------------------------------------------------------------

      If ( (collectn_start-gra_start .gt. DELTA)   .OR.
     1     (gra_stop-collectn_stop   .gt. DELTA) ) Then

         error_flag = .true.

         msgbuf = 
     1   'Input granule (LUN=' // msg25_LUN(fbyte_LUN:lbyte_LUN) 
     2   // ', VRSN=' // msg25_VRSN(fbyte_VRSN:lbyte_VRSN) // ') start and stop times do not ' 
     3   // char(10) // 'fit within processing interval defined on PCF LUNs 10258 and 10259.'
     4   // char(10) // 'Operator Action:  Check for wrong or corrupted PCF and/or input '
     5   // char(10) // 'file.  If a fault is identified, stage correct PCF/input file '
     6   // char(10) // 'and rerun PGE.  Otherwise, notify SDST.'


         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


      If (.not.error_flag) chk_input_L2 = SUCCEED

      Return

      End
