       FUNCTION PMFIN(modfil,attr,dtype,nms,value)
C-----------------------------------------------------------------------
C!F77
C!Purpose:	Writes a MODIS-HDF file attribute/value pair to a file.
C 
C!Description: Function PMFIN is part of a larger software system called the 
C MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API.  The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read  and 
C write data and metadata from/to MODIS HDF files. The 
C functionality of the M-API is defined in the MODIS Application 
C Program Interface (API) Specification.
C 
C PMFIN stores an attribute = value(s) attribute pair to the indicated 
C MODIS-HDF file.  If the attribute already exists, the value(s) will 
C be updated.
C 
C!Input parameters:
C	 modfil	IN:	Array that is used to reference the MODIS-HDF 
C		file.
C	 attr	IN:	Name to assign the attribute.  Provided macros 
C		for accepted MODIS file metadata names are listed in 
C		Appendix A, MODIS API Supplied Constants.
C	 dtype	IN:	Data Type of the value.
C 
C		Permitted FORTRAN data types:
C			CHARACTER*(*)
C			INTEGER*1
C			UINTEGER*1
C			INTEGER*2
C			UINTEGER*2
C			INTEGER*4
C			UINTEGER*4
C			INTEGER*8
C			REAL*4
C			REAL*8
C 
C	 nms	IN:	Number of metadata values to extract from 
C		value.
C	 value	IN:	The data to store in the attribute.  If the 
C		attribute already exists, the value will be updated.  
C		Values should conform to the data types, formats 
C		and/or those values enumerated for the attribute in 
C		Appendix A, MODIS_API - Supplied Constants.
C 
C!Output parameters:		NONE
C 
C Returns:	MAPIOK if successful, MFAIL if an error occurs.
C
C External references:
C 			MODFILLEN	(mapi.inc)
C			cpmfin		(mapic.h)
C
C!Revision History:
C		Qi Huang	1996/03/18
C		Version 2.0
C		Original development and testing
C 
C $Log: PMFIN.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.4  1996/04/26  20:13:24  qhuang
c Added 'UINTEGER*1 UINTEGER*2 UINTEGER*4' in prolog.
c
c Revision 1.3  1996/04/26  18:43:55  qhuang
c Added bangs, etc. in prolog, removed 'NTEGERR*1' in prolog.
c
c Revision 1.2  1996/03/19  15:06:57  qhuang
c Version 2.0 original development.
c
C
C!Team-unique header:
C	 This software is developed by the MODIS Science Data SupportTeam for 
C	 the National Aeronautics and Space Administration, Goddard Space 
C	 Flight Center, under contract NAS5-32373.
C 
C!References and Credits:
C
C!Design Notes:
C
C-----------------------------------------------------------------------
C!END
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER  modfil(MODFILLEN), nms
      CHARACTER*(*) attr, dtype
      BYTE          value(*)

      INTEGER  ret
      CALL cpmfin(modfil,attr,len(attr),dtype,len(dtype),nms,value,ret)
      PMFIN = ret
      RETURN
      END
