       FUNCTION GMARIN(modfil,arrnm,grpnm,attr,dtype,nms,value)
C
C-----------------------------------------------------------------------
C!F77
C!Purpose: 	Reads the value(s) of a local attribute attached to an SDS (array) 
C 		from a MODIS HDF file.
C 
C!Description:	Function GMARIN is part of a larger software system called the 
C		MODIS Applications Programming Interface (API) Utility, 
C		abbreviated M-API.  The M-API Utility consists of subroutines 
C		which allow MODIS Science Team-supplied software to read  and 
C		write data and metadata from/to HDF files. The functionality of 
C		the M-API is defined in the MODIS Application Program Interface 
C		(API) Specification.
C		
C		GMARIN retrieves the value(s) associated with an attribute = 
C		values pair attached to an SDS (array) by giving the array name 
C		and attribute name.  If the attribute cannot be found, the routine 
C		will return MFAIL and the passed variable unchanged.
C		
C		The routine will also fail if the provided dtype is found to be 
C		different from the attribute data type or the number of 
C		elements (nms) is found to be too small to contain the attribu
C		value.  GMARIN replaces this input information with the actual 
C		data type and number of elements contained in the attribute 
C		value (in the case of character data, it is the length of the string).  
C		These attr may be may be used to properly retrieve the 
C		attribute value(s) with a second call to the routine. If a function 
C		failure occurs,  nms will be set to zero.
C		
C		A variable of the proper data type should be passed for the value 
C		parameter.  The data type information required to properly use 
C		either routine may be found in Appendix A, M-API-Supplied 
C		Constants, and Appendix C, MODIS Data Product File Definitions of 
C		M-API User's Guide  Appendix A has a listing for each M-API 
C		provide attributes that include the data type, the 
C		format, and/or specific values associated with it.
C 
C!Input parameters:
C     modfile	IN: 	array  that is used to reference a MODIS HDF 
C		file containing the attribute.
C 	arrnm	IN:	ASCII string name of the array.  Provided 
C		macros for accepted MODIS HDF file array names are 
C 		listed in Appendix A, M-API-Supplied Constants.
C 	grpnm	IN:	ASCII string name of the data group 
C		containing the array structure to which the 
C		attribute is attached.  If set to NULL the entire file 
C		will be searched for the array structure named 
C		arrnm.
C 	attr	IN:	ASCII string name of the attribute.  Provided 
C 		macros for accepted MODIS HDF file attribute names 
C 		are listed in Appendix A, M-API-Supplied Constants.
C 	dtype	IN/OUT: ASCII string of data type of the value output.  
C 		Output replaces with the data type of the retrieved 
C 		attribute. The memory size of dtype should be at least 
C 		13 characters long. 
C 
C		 Permitted FORTRAN data types:
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
C	nms	IN/OUT: The number of elements the value buffer 
C		can contain.  Output replaces with the number of 
C		elements required to contain the attribute. If a 
C		function failure occurs, the value will be set to zero.
C 
C!Output parameters:
C	 value	OUT:	Values associated with the attribute.
C 
C Returns:	MAPIOK if successful, MFAIL if value cannot contain the 
C 		retrieved attribute value, the data type is different, the attr 
C		 cannot be found, or an error occurs.
C
C External references:
C 			MODFILLEN	(mapi.inc)
C			cgmarinn	(mapic.h)
C
C!Revision History:
C		Qi Huang	1996/03/14
C		Version 2.0
C		Original development and testing
C $Log: GMARIN.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.1  1996/04/26  19:36:44  qhuang
c Initial revision
c
C 
C!Team-unique header:
C 	This software is developed by the MODIS Science Data Support Team for 
C 	the National Aeronautics and Space Administration, Goddard Space 
C 	Flight Center, under contract NAS5-32373.
C 
C!REFERENCES AND CREDITS
C
C!DESIGN NOTES:
C
C-----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER  modfil(MODFILLEN), nms
      CHARACTER*(*) arrnm, grpnm, attr, dtype
      BYTE          value(*)

      INTEGER  ret

      CALL cgmarin(modfil,arrnm,len(arrnm),grpnm,len(grpnm),attr,
     1             len(attr),dtype,len(dtype),nms,value,ret)
      GMARIN = ret
      RETURN
      END
