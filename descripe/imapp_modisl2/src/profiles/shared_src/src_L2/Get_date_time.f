      SUBROUTINE GET_DATE_TIME(DATE,LOCAL_TIME,GMT_DELTA)

C------------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C     This subroutine get the current time.
C
C !INPUT PARAMETERS: N/A
C
C
C !OUTPUT PARAMETERS:
C
C     character*10     date                 current date
C     character*10     local_time           current local time
C     character*10     UTC_Local_Time_Diff  The time different with respect
C                                           to UTC time(GMT)
C
C !REVISION HISTORY:
C
C
C !TEAM-UNIQUE HEADER:
C
C      THIS SOFTWARE WAS DEVELOPED BY THE MODIS SCIENCE DATA SUPPORT
C      TEAM FOR THE NATIONAL AERONAUTICS AND SPACE ADMINISTRATION,
C      GODDARD SPACE FLIGHT CENTER, UNDER CONTRACT NAS5-32373.
C
C
C !REFERENCES AND CREDITS:
C
C      WRITTEN BY Liqun Ma  11/05/97
C      RESEARCH AND DATA SYSTEMS CORPORATION
C      SAIC/GSC MODIS SCIENCE DATA SUPPORT OFFICE
C      7501 FORBES BLVD, SEABROOK MD 20706
C
C      lma@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C
C    Externals:     
C
C       Functions and Subroutines: N/A
C
C       Named Constant: N/A
C
C    Internals:
C
C       Functions and Subroutines:
C
C       Variables: N/A
C
C
C Internals:
C
C
C !END
C-----------------------------------------------------------------------

      implicit none
      
      character*10 date
      character*10 local_time
      character*10 UTC_Local_Time_Diff

      integer  date_time(8)
      integer  GMT_delta

      call date_and_time(date,local_time,  
     &     UTC_Local_Time_Diff,date_time)

      GMT_delta=date_time(4)

      END
