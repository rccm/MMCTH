      FUNCTION GMHOID(modfil,objname,group,type,aaccess) 
C
C--------------------------------------------------------------------------
C!F77
C
C!Purpose:	 get the id of an existed HDF object (either array or 
C table) in a MODIS HDF file so that a M-API application can directly
C call the HDF library functions to access the data in MODIS HDF file in
C case that the application program needs to do so.
C 
C!Description: 	Function GMHOID is part of a larger software system 
C called the MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API.  The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read  and 
C write data and metadata from/to HDF files. The functionality of M-
C API is defined in the MODIS Application Program Interface (API) 
C Specification.
C 
C GMHOID returns the HDF object id. The application program can 
C directly use this id in the native HDF library function calls 
C because the HDF object is already attached. If the type is specified 
C as MODIS_ARRAY, the return value is an array id (sds id). If the 
C type is MODIS_TABLE, the return value is the table id (vdata id). 
C The intended use of this function is for M-API applications which 
C want  to call HDF functions directly. The call sequence in an M-
C API applications should be:
C 1). Use OPMFIL to open a MODIS HDF file.
C 2). Call GMHOID to obtain an object's HDF id.
C 3). Call HDF library functions by using the obtained id.
C 4). Call CLMFIL or CPMFIL to end the access.
C The calls to HDF library functions can be mixed with calls to M-
C API functions.
C NOTE: Please do not use the HDF vsfdtch (for table) or sfendacc 
C (for array) to end the access to the object. The CLMFIL/CPMFIL 
C will take care of all opened objects. In case  of application 
C programs need to end the access to an object while keeping the 
C MODIS file open,  use the M-API function EMOBJ.
C 
C!Input Parameters:
C 
C modfil	IN:  Array that is used to reference the MODIS-HDF file.
C objname		IN:  The name of the object.
C group		IN:  The name of the group to which the object 
C 		belongs. If set to ' ' the entire file will be searched 
C 		for the object named objname.
C type		IN:  The type of object: MODIS_ARRAY (for SDS, 
C		the numerical value is 720, the same value as 
C		DFTAG_NDG),  MODIS_TABLE(for Vdata, the 
C		numerical value is 1962, the same value as 
C		DFTAG_VH).
C aaccess	IN:  'r'	The object is attached for ready only.
C		'w'	The object is attached for read/write.
C		If the MODIS file is opened for read only, and the 
C		application program calls this function to obtain an 
C		object id with write operation, the function will 
C		return an error.
C 
C!Output parameters:		None
C 
C 
C Returns:	The object id or MFAIL.
C 
C
C External reference:
C		MODFILLEN	(mapi.inc)
C     		cgmhoid		(mapic.h)
C
C!Revision History:
C   		Qi Huang	1996/11/26
C		Version 2.2
C		Ring super structure and other changes make this version
C		much faster.
C
C $Log: GMHOID.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
C
C !Team-unique Header:
C             This software is developed by the MODIS Science Data Support
C             Team for the National Aeronautics and Space Administration,
C             Goddard Space Flight Center, under contract NAS5-32373.
C
C !References and Credits:
C             Portions developed at the National Center for Supercomputing
C             Applications at the Univ. of Illinois at Urbana-Champaign.
C
C----------------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapi.inc'
      INTEGER        modfil(MODFILLEN), type
      CHARACTER*(*)  objname, group, aaccess

      INTEGER    ret
 
      call cgmhoid(modfil,objname,len(objname),group,len(group),
     +             type,aaccess,len(aaccess),ret)
      GMHOID = ret
 
      RETURN
      END
