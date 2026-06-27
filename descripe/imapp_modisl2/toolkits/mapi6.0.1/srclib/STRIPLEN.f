      SUBROUTINE STRIPLEN(string,beg,eend)
C
C-----------------------------------------------------------------------
C!F77
C
C!Description: 	Subroutine STRIPLEN is part of a larger software system called
C               the MODIS Applications Programming Interface (API) Utility,
C               abbreviated M-API.  The M-API Utility consists of subroutines 
C		which allow MODIS Science Team-supplied software to read in 
C		Level 1B radiance bands and write out output products and
C               metadata to HDF files.  The functionality of the M-API is 
C		defined in the MODIS API User's Guide, Version 1, dated 4/3/95.
C
C               STRIPLEN finds the positions of the first and last nonblank 
C		characters of a character string.
C
C !Input Parameters:
C
C		CHARACTER*(*) string: A string variable of arbitrary length. 
C
C !Output Parameters:
C
C               INTEGER beg:          The byte location of the first nonblank
C                                     character of the input string.
C               INTEGER eend:          The byte location of the last nonblank
C                                     character of the input string.
C
C !Revision History:
C $Log: STRIPLEN.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
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
C               WRITTEN BY:  
C        
C		Dave Lorenzi		  05/16/95
C		Vicky Lin 		  06/12/95
C               SAIC/GSC MODIS Science Data Support Office
C               7501 Forbes Blvd
C               Seabrood, MD 20706
C
C               vlin@modis-xl.gsfc.nasa.gov
C
C !Design Notes:
C
C	Externals:
C
C	   Function:
C		LEN	Returns the length of an arbitrary character string.
C
C-----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
        CHARACTER*(*) string
        INTEGER beg, eend

C  Initialization

        beg=1
        eend=len(string)

C  Look for blank characters

        DO WHILE ((string(beg:beg).eq.' ').and.(beg.le.eend))
           beg=beg+1
        END DO

        eend=len(string)

        DO WHILE ((string(eend:eend).eq.' ').and.(eend.ge.1))
           eend=eend-1
        END DO

        RETURN
        END
