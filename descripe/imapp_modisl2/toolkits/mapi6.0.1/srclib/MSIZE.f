      FUNCTION MSIZE(dtype)
C
C-----------------------------------------------------------------------
C!F77
C 
C !Description:Function MSIZE is part of a larger software system 
C               called the MODIS Applications Programming Interface 
C               (API) Utility, abbreviated M-API.  The M-API Utility con-
C               sists of subroutines which allow MODIS Science Team-sup-
C               plied software to read in Level 1B radiance bands and 
C               write out output products and metadata to HDF files. 
C 
C               The MODIS API uses a set of standard strings to describe
C               the data types stored in arrays and table structures.   
C               These strings are returned, for example, by the routine
C               getMODISardims (GMARDM) to describe the data type of the
C               target array structure.  MSIZE (MSIZE) returns the number 
C               of bytes required to store a data type given this data
C               type string.  The input string may be a series of comma
C               delimited data type strings, in which case the total num-
C               ber of bytes to store the record described by the string   
C               is returned.
C
C !Input Parameters:
C            dtype      IN:    String of comma-delimited data types.
C
C                       Permitted FORTRAN data types:
C                              'CHARACTER*(*)'
C                              'INTEGER*1'
C                              'UINTEGER*1'
C                              'INTEGER*2'
C                              'UINTEGER*2'
C                              'INTEGER*4'
C                              'UINTEGER*4'
C                              'INTEGER*8'
C                              'UINTEGER*8'
C                              'REAL*4'
C                              'REAL*8'
C
C !Output Parameters:NONE
C
C RETURNS:    Number of bytes required for a string of data types, or
C             MFAIL if the data type was not identified.
C
C EXTERNAL REFERENCES:
C
C               I8                (mapi.inc) 
C               I16               (mapi.inc) 
C               I32               (mapi.inc) 
C               I64               (mapi.inc) 
C               UI8               (mapi.inc) 
C               UI16              (mapi.inc)
C               UI32              (mapi.inc)
C               UI64              (mapi.inc)
C               R32               (mapi.inc)
C               R64               (mapi.inc)
C               TXT               (mapi.inc)
C               MFAIL             (mapi.inc)
C               MAPIOK            (mapi.inc)
C		DATATYPELENMAX	  (mapic.inc)
C
C !Revision History:
C
C             Revision 01.00  1995/9/22
C             Qi Huang/RDC    (qhuang@modis-xl.gsfc.nasa.gov)
C             Original Development/testing
C
C             Liping Di, Hughes STX, Oct 17, 1995
C             Set dstring to ' ' instead of '' for the portability
C
C	      Jeffrey J. Blanchette, Oct 31, 1995
C	      Set dstring character array with parameter
C	      Installed array dimension check for dstring
C
C $Id: MSIZE.f,v 1.1 1999/04/15 19:26:51 jayshree Exp $
C
C $Log: MSIZE.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.1  1995/10/12  15:53:55  qhuang
c Initial revision
c
C
C !Team-unique Header:
C
C               This software is developed by the MODIS Science Data Support
C               Team for the National Aeronautics and Space Administration,
C               Goddard Space Flight Center, under contract NAS5-32373.
C
C               Portions developed at the National Center for Supercomputing
C               Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !References and Credits:
C
C !Design Notes:
C
C-----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE 'mapic.inc'

      CHARACTER*(*) dtype
      CHARACTER*(DATATYPELENMAX) dstring
      INTEGER       bytsum
      INTEGER       start, eend
      INTEGER       count, i

C     Set MSIZE to MFAIL
      MSIZE = MFAIL

C     Sets byte sum to 0
      bytsum = 0

      CALL STRIPLEN(dtype,start,eend)
    
      IF (eend .GT. start) THEN
C       Empty the local declared string dstring
        dstring = ' '
 
C       Set count to one
        count = 1

        DO i = start, eend
          IF ( (dtype(i:i) .EQ. ',') .OR. (i .EQ. eend) ) THEN
            IF ( (i .EQ. eend) .AND. (dtype(i:i) .NE. ',') ) THEN
               dstring(count:count) = dtype(i:i)
            END IF

            IF (dstring .EQ. I8) THEN
              bytsum = bytsum + 1

            ELSE IF (dstring .EQ. UI8) THEN
              bytsum = bytsum + 1

            ELSE IF (dstring .EQ. I16) THEN
              bytsum = bytsum + 2
   
            ELSE IF (dstring .EQ. UI16) THEN
              bytsum = bytsum + 2

            ELSE IF (dstring .EQ. I32) THEN
              bytsum = bytsum + 4

            ELSE IF (dstring .EQ. UI32) THEN
              bytsum = bytsum + 4

C  Comment out 8 byte integer conversions.  DFNT_INT64 and DFNT_UINT64
C  operations are not supported by HDF 3.3r4.
C           ELSE IF (dstring .EQ. I64) THEN
C             bytsum = bytsum + 8
C
C           ELSE IF (dstring .EQ. UI64) THEN
C             bytsum = bytsum + 8

            ELSE IF (dstring .EQ. R32) THEN
              bytsum = bytsum + 4

            ELSE IF (dstring .EQ. R64) THEN
              bytsum = bytsum + 8

            ELSE IF (dstring .EQ. TXT) THEN
              bytsum = bytsum + 1

            ELSE
              RETURN

            ENDIF
             
C           Empty dstring
            dstring = ' '
 
C           Set count to one
            count = 1

          ELSE
            dstring(count:count) = dtype(i:i)
            count = count + 1
C	    Return with error if data type is too long
            IF (count .GT. DATATYPELENMAX) RETURN

          END IF

        END DO

      END IF

C   Set MSIZE to byte sum.
      MSIZE = bytsum

      RETURN
      END
