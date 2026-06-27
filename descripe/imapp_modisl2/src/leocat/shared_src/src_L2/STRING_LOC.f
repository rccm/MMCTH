      INTEGER FUNCTION STRING_LOC(string, fbyte, lbyte)
      IMPLICIT NONE

C-----------------------------------------------------------------------
C !F77
C
C !Description:  
C
C     STRING_LOC finds the position (first and last bytes) of the 
C     nonblank characters within a string buffer.  If the string 
C     consists of all blank characters, both fbyte (first byte) and 
C     lbyte (last byte) are set to 1, and the return value is -1.  
C     Otherwise, the return value is 0. 
C
C !INPUT PARAMETERS:
C
C     string - a character variable of arbitrary length
C
C !OUTPUT PARAMETERS:
C
C     fbyte - byte location of first nonblank character in input string
C
C     lbyte - byte location of last nonblank character in input string.
C
C !REVISION HISTORY:
c Revision 2.2  1998/12/09  14:08:47  rhucek
c 'For BLANK input string, fbyte and lbyte are set to 1, rather than -999.
c
c Revision 2.1  1996/04/24  17:50:17  vlin
c put in "!" sign
c
c Revision 1.3  1995/08/24  16:49:50  rhucek
c Converted subroutine to function.  Function returns 0
c for nonblank string and -1 for an all blank string.
c
c Revision 1.2  1995/07/05  16:43:27  vlin
c substitute tabs with spaces
c
c Revision 1.1  1995/06/21  17:21:07  rhucek
c Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C     This software was developed by the MODIS Science Data Support
C     Team for the National Aeronautics and Space Administration,
C     Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS
C
C     Written by Richard Hucek
C     Research and Data systems Corporation
C     SAIC/GSC MODIS Science Data Support Office
C     7501 Forbes Blvd, Seabrook MD 20706
C
C     rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C     There is no check to identify the unpredictable value of
C     an undefined string.  Consequently, users must take care to 
C     initialize all string variables before passing them to 
C     STRING_LOC.
C
C Externals:    None
C
C Internals:
C
C     Variables:
C          string_len       size of input string buffer
C     
C     STRING_LOC Return Values:
C              0            nonblank text found
C             -1            all blank string buffer
C
C !END
C-----------------------------------------------------------------------

      CHARACTER*(*) string
      INTEGER fbyte, lbyte, string_len

C Initialize variables
      string_len=len(string)
      fbyte=1
      lbyte = string_len
      STRING_LOC = 0

C Determine byte position of first non-blank character.
      DO WHILE ( (string(fbyte:fbyte).eq.' ') .and. (fbyte.le.lbyte) )
         fbyte=fbyte+1
      END DO

C For an all blank string, return STRING_LOC=-1 
      IF (fbyte .gt. string_len) THEN
         fbyte = 1 
         lbyte = 1 
         STRING_LOC = -1
         RETURN
      END IF
 
C Determine byte position of last non-blank character.
      DO WHILE ( (string(lbyte:lbyte).eq.' ').and.(lbyte.ge.1) )
         lbyte=lbyte-1
      END DO

      RETURN
      END
