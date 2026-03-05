#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mapic.h"

void ncgmfin(intf modfil[MODFILLEN], _fcd attr, intf *lat, _fcd dtype,
		 intf *ldt, intf *nms, VOIDP value, intf *ret)
/*
!C***********************************************************************
* 
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* getMODISfileinfo. This C function is only called by FORTRAN function 
* GMFIN. This function is a M-API internal routine.
* 
*!Description:	Function cgmfin is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API.  The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read  and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cgmfin is a C function which is callable from FORTRAN. This 
* function will call getMODISfileinfo to get a file (global) attribute. 
* In M-API, cgmfin is low-level routine called only by GMFIN. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cgmfin, this function is called ncgmfin in 
* the actual C code. ncgmfin is redefined in mapic.h according to 
* compiler's FORTRAN naming conventions/conversion of each 
* platform, so that the object name of ncgmfin will always be the 
* object name of a FORTRAN function named cgmfin.
* 
*!Input parameters:
*	 modfil	IN:	FORTRAN array of MODIS file structure.
*	 attr	IN:	FORTRAN ASCII string of  name of the 
*		attribute to be retrieved.
*	 lat	IN:	Address for the number of bytes in attr.
*	 dtype	IN/OUT: FORTRAN ASCII string of data type.
*	 ldt	IN:	Address for the number of bytes in dtype.
*	 nms	IN/OUT:	Address for the number of elements can be 
*		stored in value. Output replaces with the actual 
*		number of elements in value.
* 
*!Output parameters:
* 	value		OUT: pointer to the data buffer.
* 	ret		OUT: MAPIOK if successful, otherwise MFAIL
* 
* Returns:	None.
*
* External references:
*			   MODFILLEN			 (mapic.h)
*			   HDf2cstring			 (hproto.h)
*			   MODFILE			 (mapi.h)
*			   DATATYPELENMAX		 (mapic.h)
*			   FDATATYPELENMAX		 (mapic.h)
*			   P_ADDR			 (mapic.h)
*			   MTYPEf2c			 (mapic.h)
*			   MTYPEc2f			 (mapic.h)
*			   HDfreespace			 (hproto.h)
*			   VOIDP			 (hdfi.h)	
*			   TXT				 (mapi.h)
*			   getMODISfileinfo		 (mapi.h)
* 
*!Revision History:
*		Qi Huang	1996/08/27
*		Version 2.1
*
*		Ring super structure and other changes make this
*		version much faster.
*
*		Qi Huang	1996/03/19
*		Version 2.0
*		Original development and testing
*
* $Log: cgmfin.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.5  1996/09/09  20:14:25  qhuang
 * Version 2.1
 *
 * Revision 1.4  1996/05/07  19:36:43  qhuang
 * Initialized variables c_length and f_length as DATATYPELENMAX and
 * FDATATYPELENMAX.
 *
 * Revision 1.3  1996/04/26  18:57:49  qhuang
 * Added '!'s and '!END' in prolog, removed 'NTEGERR*1' in prolog.
 *
 * Revision 1.2  1996/03/19  22:17:12  qhuang
 * Version 2.0.  Original development and testing.
 *
 * Revision 1.1  1996/03/19  22:16:30  qhuang
 * Initial revision
 *
* 	
*!Team-unique header:
*	This software is developed by the MODIS Science Data Support Team for 
*	the National Aeronautics and Space Administration, Goddard Space 
*	Flight Center, under contract NAS5-32373.
* 
*!References and Credits:
*
*!Design Notes:
*
!END********************************************************************
*/
{
  MODFILE	*file;
  char		*cattrname, *fdtype, *fcdtype;
  long int	n_elements; 
  char		cdtype[DATATYPELENMAX]; 
  char		cfdtype[FDATATYPELENMAX]; 
  long int	c_length = DATATYPELENMAX; 
  long int	f_length = FDATATYPELENMAX; 
  long int	i; 

  /* Set file by memcpy */ 
  memcpy(&file,&modfil[P_ADDR],sizeof(MODFILE *)); 
  n_elements = *nms; 

  /* Convert FORTRAN strings to C strings by using HDf2cstring */ 
  fdtype = HDf2cstring(dtype, (intn) *ldt); 
  cattrname = HDf2cstring(attr, (intn) *lat); 

  cdtype[0] = '\0'; 
  MTYPEf2c(fdtype,cdtype,&c_length);

  /* Convert dtype to c pointer fcdtype by _fcdtocp */
  fcdtype = _fcdtocp(dtype);

  *ret = getMODISfileinfo(file,cattrname,cdtype,&n_elements,value);

  if ( n_elements != 0 )
  {
    /* (Note: make FORTRAN dtype) */
    if ( *cdtype != '\0' )
    {
      MTYPEc2f(cdtype,cfdtype,&f_length);
      memset(fcdtype,' ',(size_t)*ldt);
      memcpy(fcdtype,cfdtype,strlen(cfdtype));
    }

    if ( strcmp(cdtype,TXT) == 0 )
    {
      n_elements = n_elements - 1;
      if ( n_elements < *nms )
	for ( i=n_elements; i< *nms; i++)
	  *( (char *)value + i ) = ' ';
    }
  }

  *nms = (intf)n_elements;

  if ( cattrname ) HDfreespace((VOIDP)cattrname);
  if ( fdtype ) HDfreespace((VOIDP)fdtype);

  return;
}
