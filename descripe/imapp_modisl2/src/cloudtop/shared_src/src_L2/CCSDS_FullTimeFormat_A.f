      Character*27 Function CCSDS_FullTimeFormat_A( DateTimeString ) 

C---------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: Converts ECS ASCII "DateTime" string to full CCSDS time 
C              code A format. 
C
C
C!INPUT PARAMETERS:
C  character*(*) DateTimeString: 
C                  Any ECS ASCII datetime string acceptable to the SDPTK 
C                  Date/Time conversions routines. (See SDPTK Users Guide 
C                  for description of acceptable full and truncated CCSDS 
C                  time code formats.)
C
C
C!OUTPUT PARAMETES: None (see Function Return Value)
C
C
C
C!REVISION HISTORY:
C
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
C    Richard Hucek 1/2000
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C DESIGN NOTES: The input ECS ASCII datetime string is converted to a double 
C               precision numeric value and then transformed back to ASCII 
C               standard CCSDS time code format A.  Both transformations  
C               are handled using SDPTK function calls.
C
C
C  Returns:     27-character "datetime" string in full CCSDS time code A 
C               format if successful.  A blank string ' ' if error.
C
C    Externals:
C       Functions and Subroutines:
C          pgs_td_utctotai           (libPGSTK.a)
C          pgs_td_taitoutc           (libPGSTK.a)
C          MODIS_SMF_SETDYNAMICMSG   (atmos shared code)
C          string_loc                (atmos shared code)
C
C       Named Constants:
C          MODIS_E_GENERIC          (PGS_MODIS_39500.f)
C          PGS_S_SUCCESS            (PGS_SMF.f)
C
C
C!END
C----------------------------------------------------------------------

      IMPLICIT NONE

      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

C PARAMETER declarations
      character*(*)     FUNCNAME
      PARAMETER        (FUNCNAME = 'CCSDS_FullTimeFormat_A')

      character*(*)     BLANK 
      PARAMETER        (BLANK = ' ')

C Function argument declarations
      character*(*)  DateTimeString

C Local variable declarations
      character*27   FullTimeFormat
      character*2048 msgbuf

      double precision tai_time
      integer    fbyte,lbyte,rtn
      logical    error_flag

C Function declarations
      integer String_Loc, pgs_td_utctotai, pgs_td_taitoutc

C-----------------------------------------------------------------------
C Initialization
C-----------------------------------------------------------------------
      CCSDS_FullTimeFormat_A = BLANK 
      error_flag             = .false.

      rtn = String_Loc(DateTimeString,fbyte,lbyte)

      rtn = pgs_td_utctotai(DateTimeString,tai_time)

      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .true.
   
         msgbuf = 
     1   'Conversion of input datetime string ("' // DateTimeString(fbyte:lbyte) // '") '  
     2   // char(10) // 'to CCSDS full time code A format failed on call to pgs_td_utctotai. ' 
     3   // char(10) // 'Operator Action:  Check for valid, up-to-date leapsec.dat file '
     4   // char(10) // 'on PCF file LUN 10301.  If leapsec.dat file is missing or out-of-date, '
     5   // char(10) // 'stage correct PCF and rerun PGE.  Otherwise, notify SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


      Else
         rtn = pgs_td_taitoutc(tai_time,FullTimeFormat)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .true.
   
            msgbuf = 
     1      'Conversion of input datetime string ("' // DateTimeString(fbyte:lbyte) // '") '  
     2      // char(10) // 'to full CCSDS time code A format failed on call to pgs_td_taitoutc. ' 
     3      // char(10) // 'Operator Action:  Check for valid, up-to-date leapsec.dat file '
     4      // char(10) // 'on PCF file LUN 10301.  If leapsec.dat file is missing or out-of-date, '
     5      // char(10) // 'stage correct PCF and rerun PGE.  Otherwise, notify SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)


         Endif 

      Endif

      If (.not. error_flag) CCSDS_FullTimeFormat_A = FullTimeFormat 

      Return

      End
