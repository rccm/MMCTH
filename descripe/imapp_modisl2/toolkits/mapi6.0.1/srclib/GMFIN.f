       FUNCTION GMFIN(modfil,attr,dtype,nms,value)
C------------------------------------------------------------------------
C!F77 
C!Purpose:  	Reads the value(s) of a global attribute from a MODIS HDF file.
C 
C!Description:	Function GMFIN is part of a larger software system called the 
C MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API.  The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read  and 
C write data and metadata from/to HDF files. The functionality of 
C the M-API is defined in the MODIS Application Program Interface 
C (API) Specification.
C 
C GMFIN retrieves the value(s) associated with an attribute = value 
C pair attached to a MODIS HDF file by giving the attribute name.  If 
C the attribute cannot be found, the routine will return MFAIL and 
C the passed variable unchanged.
C 
C The routine will also fail if the provided dtype is found to be 
C different than the attribute’s data type or the number of 
C elements (nms) is found to be too small to contain the attribute’s 
C value.  GMFIN replaces this input information with the actual 
C data type and number of elements contained in the attribute 
C value (in the case of character data, it is the length of the string).  
C These attribute's attributes may be used to properly retrieve the 
C attribute value(s) with a second call to the routine.
C 
C A variable of the proper data type should be passed for the value 
C parameter.  The data type information required to properly use 
C either routine may be found in Appendix A, M-API-Supplied 
C Constants, and Appendix C, MODIS Data Product File Definitions, of 
C M-API User's Guide.  Appendix A has a listing for each M-API 
C provided attribute’s attributes that include the data type, the 
C format, and/or specific values associated with it.
C 
C!Input parameters:
C	 modfil	IN: 	array that is used to reference a MODIS HDF 
C		file containing the attribute.
C	 attr	IN:	ASCII string name of the attribute.  Provided 
C		macros for accepted MODIS HDF file attr names are 
C		listed in Appendix A, M-API-Supplied Constants.
C	 dtype	IN/OUT: Data type of the value output.  Output 
C		replaces with the data type of the retrieved 
C		attribute. This memory size of this character should 
C		be at least 13 characters long.
C 
C		 Permitted FORTRAN data types:
C			CHARACTER*(*)
C			INTEGER*1
C			INTEGER*2
C			UINTEGER*2
C			INTEGER*4
C			UINTEGER*4
C			INTEGER*8
C			REAL*4
C			REAL*8
C 
C	 nms	IN/OUT: Number of elements the value buffer can 
C		contain.  Output replaces with the number of 
C		elements required to contain the attribute. If a 
C		function error is occur, the value is set to zero. If 
C		the requested attribute does not exist, nms will be 
C		unchanged.
C 
C!Output parameters:
C	value	OUT:	Values associated with the attribute.
C 
C Returns:	MAPIOK if successful, MFAIL if value cannot contain the 
C		retrieved attribute value, the attr cannot be found, or an error 
C		occurs.
C 
C External references:
C		MODFILLEN		(mapi.inc)
C		getMODISfileinfo	(mapi.h)	
C
C!Revision History:
C 		Qi Huang	1996/03/19
C		Version 2.0
C		Original development and testing
C
C $Log: GMFIN.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.4  1996/04/26  20:15:21  qhuang
c Added 'UINTEGER*1 UINTEGER*2 UINTEGER*4' in prolog.
c
c Revision 1.3  1996/04/26  18:56:28  qhuang
c Added bangs,etc. in prolog, removed 'NTEGERR' in prolog.
c
c Revision 1.2  1996/03/19  22:20:54  qhuang
c Version 2.0.  Original development and testing.
c
c Revision 1.1  1996/03/19  22:20:28  qhuang
c Initial revision
c
C
C!Team-unique header:
C	This software is developed by the MODIS Science Data SupportTeam for 
C	the National Aeronautics and Space Administration, Goddard Space 
C	Flight Center, under contract NAS5-32373.
C 
C!REFERENCES AND CREDITS
C
C!DESIGN NOTES:
C
C-----------------------------------------------------------------------
C!END
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER  modfil(MODFILLEN), nms
      CHARACTER*(*) attr, dtype
      BYTE          value(*)

      INTEGER  ret

      CALL cgmfin(modfil,attr,len(attr),dtype,len(dtype),nms,value,ret)
      GMFIN = ret
      RETURN
      END
