#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mapic.h"

void ncpmfin(intf modfil[MODFILLEN], _fcd attr, intf *lat, _fcd dtype, 
		intf *ldt, intf *nms, VOIDP value, intf *ret)
/*
!C***********************************************************************
* 
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* putMODISfileinfo. This C function is called only by FORTRAN function 
* PMFIN. cpmfin is a M-API low-level internal routine.
* 
*!Description: Function cpmfin is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read  and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cpmfin is a C function callable from FORTRAN. This function will 
* call putMODISfileinfo to put a global (file) attribute.
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cpmfin, this function is named ncpmfin in 
* the actual C code. ncpmfin is redefined in mapic.h according to 
* FORTRAN compilers' naming conventions/conversion of each 
* platform, so that the object name of C ncpmfin will always be the 
* same as the object name of a FORTRAN function named cpmfin.
* 
*!Input parameters:
*	 modfil	IN: 	FORTRAN integer array used to reference the 
*		MODIS-HDF file.
*	 attr	IN:	FORTRAN character string  of the attribute 
*		name.
*	 lat	IN:	FORTRAN integer address of the memory size 
*		of attr.
*	 dtype	IN:	FORTRAN character string of the memory size 
*		of the attribute data type.
*	 ldt	IN:	FORTRAN integer address of the memory size 
*		of dtype.
*	 nms	IN:	FORTRAN integer address of the number of 
*		elements in value.
*	 value	IN:	value buffer.
* 
*!Output parameters:
*	 ret	OUT:	FORTRAN integer address of the return 
*		status(MFAIL, 				MAPIOK) 
* 
* Returns:	none
* 
* External references:
*			   MODFILLEN			 (mapic.h)
*			   MODFILE			 (mapi.h)
*			   HDf2cstring			 (hproto.h)
*			   DATATYPELENMAX		 (mapic.h)
*			   P_ADDR			 (mapic.h)
*			   MTYPEf2c			 (mapic.h)
*			   HDfreespace			 (hproto.h)
*			   VOIDP			 (hdfi.h)	
*			   putMODISfileinfo		 (mapi.h)
*!Revision history:
*		Qi Huang	1996/08/27
*		Version 2.1
*
*		Ring super structure and other changes make this
*		version much faster.
*
*		Qi Huang	1996/03/18
*		Version 2.0
*		Original development and testing
* $Log: cpmfin.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.5  1996/09/09  20:27:24  qhuang
 * Version 2.1
 *
 * Revision 1.4  1996/05/07  19:41:12  qhuang
 * Initialized variable c_length as DATATYPELENMAX.
 *
 * Revision 1.3  1996/04/26  18:49:12  qhuang
 * Added '!'s in prolog.
 *
 * Revision 1.2  1996/03/19  15:00:29  qhuang
 * Version 2.0 original development.
 *
 * Revision 1.1  1996/03/19  15:00:03  qhuang
 * Initial revision
 *
* 
*!Team-unique header:
*	This software is developed by the MODIS Science Data SupportTeam for 
*	the National Aeronautics and Space Administration, Goddard Space 
*	Flight Center, under contract NAS5-32373.
* 
*!References and Credits
*    Portions developed at the National Center for Supercomputing
*    Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes
*
!END********************************************************************
*/
{
  MODFILE	*file;
  char		cdtype[DATATYPELENMAX];
  char		*attribute;
  char		*cfdtype;
  long int	c_length = DATATYPELENMAX;

  /* Set file by memcpy */
  memcpy(&file,&modfil[P_ADDR],sizeof(MODFILE *));

  /* Convert FORTRAN strings to C strings using HDf2cstring. */ 
  attribute = HDf2cstring(attr, (intn)*lat);
  cfdtype = HDf2cstring(dtype, (intn)*ldt);

  cdtype[0] = '\0';

  MTYPEf2c(cfdtype,cdtype,&c_length); 

  *ret = putMODISfileinfo(file,attribute,cdtype,*nms,value);

  /* Free the C strings allocated by HDf2cstring if they are not null by
     using HDfreespace. */

  
  if (attribute) HDfreespace((VOIDP)attribute);
  if (cfdtype) HDfreespace((VOIDP)cfdtype);

  return;
}
