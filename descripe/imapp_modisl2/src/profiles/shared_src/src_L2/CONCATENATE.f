      SUBROUTINE CONCATENATE(instring1, instring2, outstring)
      IMPLICIT NONE

C------------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: 
C
C     Concatenate locates the message text in 2 separate input strings,
C     removes blank characters before and after each message, and 
C     joins the messages with a blank between them.  
C
C !INPUT PARAMETERS:  
C
C     instring1 - input character string 1
C     instring2 - input character string 2
C
C !OUTPUT PARAMETERS:  
C
C     outstring - string consisting of concatenated strings 1 and 2
C
C !REVISION HISTORY:
c Revision 2.2  1997/11/13  20:20:54  lma
c move out the initialization of fbytes.
c
c Revision 2.1  1996/04/24  18:01:00  vlin
c put in "!" sign
c
c Revision 1.6  1995/08/24  16:50:50  rhucek
c Revised call to function STRING_LOC which previously was implemented
c as a subroutine.  CONCATENATE now contains additional logic to check
c and act upon the function return status from STRING_LOC.
c
C !TEAM-UNIQUE HEADER: 
C
C      THIS SOFTWARE WAS DEVELOPED BY THE MODIS SCIENCE DATA SUPPORT
C      TEAM FOR THE NATIONAL AERONAUTICS AND SPACE ADMINISTRATION,
C      GODDARD SPACE FLIGHT CENTER, UNDER CONTRACT NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C      WRITTEN BY Richard Hucek 
C      RESEARCH AND DATA SYSTEMS CORPORATION
C      SAIC/GSC MODIS SCIENCE DATA SUPPORT OFFICE
C      7501 FORBES BLVD, SEABROOK MD 20706
C
C      rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C      If one of the input strings is blank, then the output contains 
C      only the text from the second string.  If both strings are blanks,
C      a blank output string is returned.
C 
C Externals:     None
C
C Internals:
C
C      Functions and Subroutines:
C               String_Loc     define byte position of nonblank text
C                              within string buffer
C
C      Variables:
C               fbyte1         beginning byte of nonblank text of 
C                              string1
C               fbyte2         beginning byte of nonblank text of 
C                              string2
C               lbyte1         ending byte of nonblank text of string1
C               lbyte2         ending byte of nonblank text of string2
C               rtn1           String_Loc return value for string1
C               rtn2           String_Loc return value for string2
C
C !END
C-----------------------------------------------------------------------

C  Declaration of variables

      CHARACTER*(*) instring1, instring2, outstring 
      INTEGER fbyte1, fbyte2, lbyte1, lbyte2, rtn1, rtn2, String_Loc

C  Initialize variables

      outstring = ' '
      rtn1 = String_Loc(instring1,fbyte1,lbyte1)
      rtn2 = String_Loc(instring2,fbyte2,lbyte2)

C  Concatenate input strings where possible

      IF ((rtn1 .EQ. 0) .AND. (rtn2 .EQ. 0)) THEN
         outstring = instring1(fbyte1:lbyte1) // ' ' // 
     &               instring2(fbyte2:lbyte2)

      ELSE IF ((rtn1 .EQ. 0) .AND. (rtn2 .EQ. -1)) THEN
         outstring = instring1(fbyte1:lbyte1)

      ELSE if ((rtn1 .EQ. -1) .AND. (rtn2 .EQ. 0)) THEN
         outstring = instring2(fbyte2:lbyte2)

      END IF

      RETURN
      END
