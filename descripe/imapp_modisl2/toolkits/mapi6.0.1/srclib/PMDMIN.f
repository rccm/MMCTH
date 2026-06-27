      FUNCTION PMDMIN(modfil, arrnm, grpnm, moddim, 
     *                        attr, dtype, nms, value) 
C---------------------------------------------------------------------------
C!F77
C
C!Purpose:	Attach a local attribute/value pair to a specific dimension of a 
C MODIS array.
C 
C!Description: 	Function PMDMIN is part of a larger software system called the 
C MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API.  The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read  and 
C write data and metadata from/to HDF files. The functionality of 
C the M-API is defined in the MODIS Application Program Interface 
C (API) Specification.
C 
C PMDMIN stores an attribute = value(s) pair to the indicated 
C dimension of an array.  If the attribute already exists, the 
C value(s) will be updated.
C 
C!Input Parameters:
C modfil	IN:	Array that is used to reference the MODIS-HDF 
C 		file containing the target array structure.
C arrnm	IN:	ASCII string name of the array.  Provided 
C		macros for accepted MODIS HDF file array names are 
C		listed in Appendix A of the User’s Guide, M-API-
C		Supplied Constants.
C grpnm	IN:	ASCII string name of the data group 
C		containing the array structure to which the 
C		attribute is attached.  If set grpnm=‘’, the entire file 
C		will be searched for the array structure named 
C		arrnm.
C moddim	IN:	The dimension to which the attribute will be 
C 		attached (0-based).
C attr	IN:	Name to assign the attribute.  Provided macros 
C		for accepted MODIS file attribute names are listed in 
C		Appendix A of the User’s Guide, MODIS API Supplied 
C		Constants.
C dtype	IN:	Data Type of the value.
C 
C              Permitted FORTRAN data types:
C                'CHARACTER*(*)'
C                'INTEGER*1'
C                'UINTEGER*1'
C                'INTEGER*2'
C                'UINTEGER*2'
C                'INTEGER*4'
C                'UINTEGER*4'
C                'INTEGER*8'
C                'REAL*4'
C                'REAL*8'
C 
C nms	IN:	Number of attribute values in value.
C		value	IN:	The data to store in the attribute.  If the 
C		attribute already exists, the value will be updated.  
C		Values should conform to the data types, formats 
C		and/or those values enumerated for the attribute in 
C		Appendix A of the User’s Guide, MODIS_API - 
C		Supplied Constants.
C 
C!Output Parameters:		NONE
C 
C Returns:	MAPIOK if successful, MFAIL if an error occurs.
C 
C External references:
C		MODFILLEN(mapi.inc)
C		cpmdmin(mapic.h)
C
C!Revision History:
C		Qi Huang	1996/08/22
C		Version 2.1
C		Original development and testing
C
C $Log: PMDMIN.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.6  1996/08/30  19:39:45  qhuang
c Version 2.1
c
c Revision 1.5  1996/08/30  19:39:25  qhuang
c *** empty log message ***
c
C
C!Team-unique header:
C               This software is developed by the MODIS Science Data Support
C               Team for the National Aeronautics and Space Administration,
C               Goddard Space Flight Center, under contract NAS5-32373.
C
C Refence and credits:
C
C!Design Notes:
C
C---------------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE 'mapic.inc'
      INTEGER       modfil(MODFILLEN), moddim, nms
      CHARACTER*(*) grpnm, arrnm, attr, dtype
      BYTE          value(*)

      INTEGER       ret

      CALL cpmdmin(modfil,arrnm,len(arrnm),grpnm,len(grpnm),moddim,
     +             attr,len(attr),dtype,len(dtype),nms,value,ret)
      PMDMIN = ret

      RETURN 
      END
