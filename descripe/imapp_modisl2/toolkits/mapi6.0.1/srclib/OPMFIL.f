      FUNCTION OPMFIL(fname,aaccess,modfil)
C
C-----------------------------------------------------------------------
C!F77
C
C!Purpose:	Opens a file for access by the MODIS API routines.
C 
C!Description:	Function OPMFIL is part of a larger software system called the 
C               MODIS Applications Programming Interface (API) Utility, 
C               abbreviated M-API. The M-API Utility consists of subroutines 
C               which allow MODIS Science Team-supplied software to read  and 
C               write data and metadata from/to HDF files. The functionality of 
C               the M-API is defined in the MODIS Application Program Interface 
C               (API) Specification.
C 
C               OPMFIL opens a file and creates the HDF structures to support the 
C               MODIS API routines access to it.  OPMFIL must be called to produce 
C               the modfil array before any of these routine can access it.  Note 
C               that setting the file access to ‘w’ creates a file and will overwrite 
C               a pre-existing one.  OPMFIL will close the file and return a zeroed 
C               array if an error occurs.
C 
C FORTRAN:	INTEGER FUNCTION OPMFIL (fname, acces, modfil)
C		      CHARACTER*(*)	(fname,aaccesss,modfil)
C		      INTEGER		modfil(aMODFIs
C 
C!Input Parameters:
C		 fname	IN:  Complete path and filename for the file to be opened.
C		 aaccess	IN:  Access mode.
C 
C		 	One of:
C 				"r"	n for read only.
C 				"w"   reate for read/write.
C 				"a"    Open for read/write.
C 
C!Output Parameters:
C		 modfil	OUT:	Array that is used to reference the file in all 
C		 other MODIS API routines.  The array will return all 
C		 zeroes if an error occurs.
C 
C Returns:	MAPIOK if successful, MFAIL if an error occurs. 
C
C External references:
C			copmfil		(mapic.h)
C 
C !Revision History:
C $Log: OPMFIL.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.2  1996/06/06  20:14:43  qhuang
c ..
c
C
C !Team-unique Header:TEAM-UNIQUE HEADER:

C             This software is developed by the MODIS Science Data Support
C             Team for the National Aeronautics and Space Administration,
C             Goddard Space Flight Center, under contract NAS5-32373.
C
C             Portions developed at the National Center for Supercomputing
C             Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !References and Credits:
C
C !Design Notes:
C
C----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      CHARACTER*(*)    fname,aaccess
      INTEGER          modfil(MODFILLEN)
      
      INTEGER  ret

      CALL copmfil(fname,len(fname),aaccess,len(aaccess),modfil,ret)
      OPMFIL = ret
      RETURN
      END
