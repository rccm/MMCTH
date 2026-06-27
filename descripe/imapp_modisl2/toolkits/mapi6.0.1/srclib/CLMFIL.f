      FUNCTION CLMFIL(modfil)
C
C-----------------------------------------------------------------------
C!F77
C
C!Purpose:	Closes a file accessed by the MODIS API routines.
C 
C!Description:	Function CLMFIL is part of a larger software system 
C called the MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API. The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read  and 
C write data and metadata from/to HDF files. The functionality of 
C the M-API is defined in the MODIS Application Program Interface 
C (API) Specification.
C 
C CLMFIL terminates the access of MODIS API routines to a MODIS-
C HDF file opened using OPMFIL.  CLMFIL may fail to close the file 
C if an error occurs.
C 
C!Input Parameters:
C 	modfil	IN:	Array that is used to reference a MODIS-HDF 
C 			file in all other MODIS API routines.
C 
C!Output Parameters:	None.
C 
C Returns:	MAPIOK if successful, MFAIL if an error occurs. 
C 
C External references:
C			MODFILLEN		(mapi.inc)
C			cclmfil			(mapic.h)
C
C!Revision History:
C $Log: CLMFIL.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.2  1996/06/06  20:14:43  qhuang
!C*** empty log message ***
c
c Revision 1.1  1996/01/02  18:53:14  qhuang
c Initial revision
c
C 		
C!Team-unique Header:
C This software is developed by the MODIS Science Data Support Team for 
C the National Aeronautics and Space Administration, Goddard Space 
C Flight Center, under contract NAS5-32373.
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
      INTEGER          modfil(MODFILLEN)
      
      INTEGER  ret

      CALL cclmfil(modfil,ret)
      CLMFIL = ret
      RETURN
      END
