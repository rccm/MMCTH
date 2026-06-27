C----------------------------------------------------------------------
C!F77
C!Purpose:	Attach a local metadata attribute/value pair to a MODIS array.
C 
C!Description: Function PMARIN is part of a larger software system called the 
C MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API.  The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read  and 
C write data and metadata from/to MODIS HDF files. The 
C functionality of the M-API is defined in the MODIS Application 
C Program Interface (API) Specification.
C 
C 
C PMARIN (PMARIN) stores an attribute = value(s) metadata pair to 
C the indicated array.  If the attribute already exists, the value(s) 
C will be updated.
C 
C!Input Parameters:
C	modfil	IN:	Array that is used to reference the MODIS-HDF 
C 		file containing the target array structure.
C 	arrnm	IN:	ASCII string name of the array.  Provided 
C 		macros for accepted MODIS HDF file array names are 
C		listed in Appendix A of the Use's Guide, M-API-
C		Supplied Constants.
C	grpnm	IN:	ASCII string name of the data group 
C		containing the array structure to which the 
C		metadata is attached.  If set grpnm = '',  the entire file 
C		will be searched for the array structure named 
C		arrnm.
C	attr	IN:	Name to assign the attribute.  Provided macros 
C		for accepted MODIS file metadata names are listed in 
C		Appendix A of the Usr's Guide, MODIS API Supplied 
C		Constants.
C	dtype	IN:	Data Type of the value.
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
C	nms	IN:	Number of metadata values in value.
C	value	IN:	The data to store in the attribute.  If the 
C		attribute already exists, the value will be updated.  
C		Values should conform to the data types, formats 
C		and/or those values enumerated for the attribute in 
C		Appendix A of ths GuideGuide, MODIS_API - 
C		Supplied Constants.
C 
C!Output Parameters:		NONE
C 
C Returns:	MAPIOK if successful, MFAIL if an error occurs.
C
C External references:
C		 MODFILLEN       (mapi.inc)
C		 cpmarin	 (mapic.h)
C
C!Revision History:
C		Qi Huang	1996/03/14
C		Version 2.0
C		Original development and testing
C 
C $Log: PMARIN.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.3  1996/04/26  20:09:00  qhuang
c Added 'UINTEGER*1 UINTEGER*2 UINTEGER*4' in prolog.
c
c Revision 1.2  1996/04/26  19:07:28  qhuang
c In prolog: added bangs, etc. , removed 'NTEGERR', changed grpn'' to
c grpnm = ''.
c
c Revision 1.1  1996/04/26  19:06:39  qhuang
c Initial revision
c
C
C!Team-unique header:
C	This software is developed by the MODIS Science Data Support Team for 
C	the National Aeronautics and Space Administration, Goddard Space 
C	Flight Center, under contract NAS5-32373.
C
C!References and Credits:
C
C!Design Notes:
C
C-----------------------------------------------------------------------
C!END
C
      FUNCTION PMARIN(modfil,arrnm,grpnm,attr,dtype,nms,value)
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER  modfil(MODFILLEN), nms
      CHARACTER*(*) arrnm, grpnm, attr, dtype
      BYTE          value(*)

      INTEGER  ret

      CALL cpmarin(modfil,arrnm,len(arrnm),grpnm,len(grpnm),attr,
     1             len(attr),dtype,len(dtype),nms,value,ret)
      PMARIN = ret
      RETURN
      END
