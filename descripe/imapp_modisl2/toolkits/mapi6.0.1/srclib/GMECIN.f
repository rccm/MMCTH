       FUNCTION GMECIN(modfil,pvlname,pname,nms,dtype,pvalue)
C-----------------------------------------------------------------------
C!F77
C!Purpose:	To retrieve  the value of a parameter from a HDF PVL text block.
C 
C!Description: Function GMECIN is part of a larger software system called the 
C MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API. The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read and 
C write data and metadata from/to HDF files. The functionality of 
C the M-API is defined in the MODIS Application Program Interface 
C (API) Specification.
C 	
C In HDF-EOS, parameters are collected together to form a text block 
C using PVL. Then the text block is stored in HDF as a single 
C attribute. GMECIN  retrieve the value of a parameter from the PVL 
C text block.
C 
C In order to obtain value of a parameter inside a PVL text block, 
C the function reads the PVL text block specified by pvlname from 
C the MODIS file, creates the internal ODL tree structure from the 
C PVL text block, and search the tree structure to retrieve the value 
C of a parameter. The tree structure is then saved internally for 
C consecutive searches in the same PVL text block for code 
C efficiency. If multiple parameters will be retrieved from the 
C same PVL block, just set pvlname to the HDF PVL attribute name 
C in the first call and set to ' '  in the consecutive calls. If the next 
C call is to retrieve the value of a parameter in a different PVL text 
C block, set the pvlname to the new PVL attribute name. The saved 
C old tree structure will be deleted automatically and a new ODL tree 
C will be created and saved. If you will no longer call GMECIN in 
C your program and want to release the memory occupied by the 
C saved tree, just set both pvlname and pname to ' ' .
C 
C!Input parameters:
C modfil	IN: 	MODFILE file array that is used to reference 
C 		the MODIS-HDF file containing the target PVL 
C 		attribute.
C pvlname	IN:	ASCII string name of the HDF attribute which 
C		contains the PVL text block. Set pvlname to ' ' while 
C		pname is  not equal to ' ' will result in searching the 
C		last PVL text block for the value of pname 
C		parameter.
C pname	IN:	ASCII string name of a parameter whose 
C		value will be retrieved.  Set both pvlname and 
C		pname to ' ' will release the memory occupied by the 
C		internal ODL tree. The pname could parameter name 
C		only or combination of name and class represented 
C		as "name.class".
C nms	IN/OUT: The number of memory elements as dtype 
C		available in the value array. The attribute's value 
C		will not be retrieved unless nms indicates that there 
C		is sufficient space available in value. 
C		getMODISECSinfo replaces this input with the 
C		number of elements required to contain the 
C		metadata. If the parameter cannot be found, *nms 
C		will be left unchanged, or set to 0 if a function error 
C		occurs. This argument must be a variable.
C 
C		SPECIAL CASE for multiple strings:
C		If there are multiple character strings for the 
C		parameters, strings will be packed together and 
C		returned in value . The separator between strings is 
C		'\0' (numerical value 0). The low 16 bit of nms will 
C		return the total bytes in value, including the '\0's. 
C		The part above the low 16 bits will return ( number 
C		of strings packed - 1). To obtain how many string 
C		retrieved, do the calculation:
C		n_strings  =  nms/65536  + 1
C		n_bytes  =  MOD(nms, 65536)
C		Therefore, if nms is less than 65536, there is only 
C		one strings in value and nms is the number of bytes 
C		(characters) in the string.
C		dtype	IN/OUT: Data type of value. Output replaces with the 
C		data type of the retrieved  metadata. There are only 3 
C		data types in PVL:
C				'CHARACTER*(*) '
C				'INTEGER*4'
C				'REAL*8'
C		But you might use "REAL*4" in input. If the value of 
C		the parameter is in 'REAL*8' type, this function will 
C		return in 'REAL*4' if users set the input value of 
C		dtype to 'REAL*4'. This argument must be a 
C		variable.The memory size for dtype should be at 
C		least DATATYPELENMAX bytes long
C 
C !Output Parameters:
C pvalue	OUT: 	buffer for the value. User should allocate 
C		enough memory for this buffer. If there are 
C		multiple data values, the value will be placed 
C		consecutively. If the data value type is 
C		"CHARACTER*(*)", strings will be separated by value 
C		0.
C 
C Returns:	MAPIOK if successful, MFAIL if an error occurs.
C
C External References:
C		MODFILLEN	(mapi.inc)
C		cgmecin		(mapic.h)
C
C!REVISION HISTORY:
C		Qi Huang	1996/05/07
C		Version 2.0
C		Original development and testing
C $Log: GMECIN.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.1  1996/05/07  17:57:01  qhuang
c Initial revision
c
C
C!Team-unique header:
C	This software is developed by the MODIS Science Data Support Team for 
C	the National Aeronautics and Space Administration, Goddard Space 
C	Flight Center, under contract NAS5-32373.
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
      CHARACTER*(*) pvlname, pname, dtype
      BYTE          pvalue(*)

      INTEGER  ret

      CALL cgmecin(modfil,pvlname,len(pvlname),pname,len(pname),dtype,
     +             len(dtype),nms,pvalue,ret)
      GMECIN = ret
      RETURN
      END
