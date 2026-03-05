#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mapic.h"

void ncpmarin(intf modfil[MODFILLEN], _fcd arrnm, intf *lar,
		_fcd grpnm, intf *lgr, _fcd attr, intf *lat, _fcd dtype, 
		intf *ldt, intf *nms, VOIDP value, intf *ret)
/*
!C**********************************************************************
* 
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* putMODISarinfo. This C function is called only by FORTRAN function 
* PMARIN. cpmarin is a M-API low-level internal routine.
* 
*!Description: Function cpmarin is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read  and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cpmarin is a C function callable from FORTRAN. This function 
* will call putMODISarinfo to put an array attribute.
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cpmarin, this function is named ncpmarin 
* in the actual C code. ncpmarin is redefined in mapic.h according 
* to FORTRAN compilers' naming conventions/conversion of each 
* platform, so that the object name of C ncpmarin will always be 
* the same as the object name of a FORTRAN function named 
* cpmarin.
* 
*!Input parameters:
*	 modfil	IN: 	FORTRAN integer array used to reference the 
*		MODIS-HDF file.
*	 arrnm	IN:	FORTRAN character string of  the name of the 
*		array.
*	 lar	IN:	FORTRAN integer address of the memory size 
*		of arrnm.
*	 grpnm	IN:	FORTRAN string name of the data group the 
*		array belongs to.
*	 lgr	IN:	FORTRAN integer address of the memory size 
*		of grpnm.
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
*			   P_ADDR			 (mapic.h)
*			   DATATYPELENMAX		 (mapic.h)
*			   MTYPEf2c			 (mapic.h)
*			   HDfreespace			 (hproto.h)
*			   VOIDP			 (hdfi.h)	
*			   putMODISarinfo		 (mapi.h)
*
*!Revision History:
*		Qi Huang	1996/08/08
*		Version 2.1
*		Ring super structure and other changes make this
*		version much faster.
*
*		Qi Huang	1996/03/15
*		Version 2.0
*		Original development and testing
* $Log: cpmarin.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.4  1996/08/15  18:41:55  qhuang
 * Added MODFILLEN and deleted P_HDFID, P_SDID and P_ACCESS in External
 * reference in prolog.
 *
 * Revision 1.3  1996/05/07  19:39:41  qhuang
 * Initialized variable c_length as DATATYPELENMAX.
 *
 * Revision 1.2  1996/04/26  19:04:43  qhuang
 * Added '!'s in prolog.
 *
 * Revision 1.2  1996/04/26  19:04:43  qhuang
 * Added '!'s in prolog.
 *
 * Revision 1.1  1996/03/18  18:14:56  qhuang
 * Initial revision
 *
* 
*!Team-unique header:
*	 This software is developed by the MODIS Science Data SupportTeam for 
*	 the National Aeronautics and Space Administration, Goddard Space 
*	 Flight Center, under contract NAS5-32373.
* 
*!References and Credits
*
*!Design Notes
*
!END********************************************************************
*/
{
  MODFILE	*file;
  char		cdtype[DATATYPELENMAX];
  char		*cfdtype;
  long int	c_length = DATATYPELENMAX;
  char		*cgroupname;
  char		*carrname;
  char		*cattrname;

  
  /* Set file by memcpy. */
  memcpy(&file,&modfil[P_ADDR],sizeof(MODFILE *));

  /* Convert FORTRAN strings to C strings using HDf2cstring. */ 
  cgroupname = HDf2cstring(grpnm, (intn)*lgr);
  carrname = HDf2cstring(arrnm, (intn)*lar);
  cattrname = HDf2cstring(attr, (intn)*lat);
  cfdtype = HDf2cstring(dtype, (intn)*ldt);

  *cdtype = '\0';

  MTYPEf2c(cfdtype,cdtype,&c_length); 

  *ret = putMODISarinfo(file,carrname,cgroupname,cattrname,cdtype,*nms,
                   	value);

  /* Free the C strings allocated by HDf2cstring if they are not null by
     using HDfreespace. */

  
  if (cgroupname) HDfreespace((VOIDP)cgroupname);
  if (carrname) HDfreespace((VOIDP)carrname);
  if (cattrname) HDfreespace((VOIDP)cattrname);
  if (cfdtype) HDfreespace((VOIDP)cfdtype);

  return;
}
