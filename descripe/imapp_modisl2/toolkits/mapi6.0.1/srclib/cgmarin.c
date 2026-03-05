#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mapic.h"

void ncgmarin(intf modfil[MODFILLEN], _fcd arrnm, intf *lar, _fcd grpnm,
	     intf *lgr, _fcd attr, intf *lat, _fcd dtype, intf *ldt,
	     intf *nms, VOIDP value, intf *ret)
/*
!C**********************************************************************
* 
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
*	 getMODISarinfo. This C function is only called by FORTRAN function 
* 	GMARIN. This function is a M-API internal routine.
* 
*!Description:	Function cgmarin is part of a larger software system called the 
*		MODIS Applications Programming Interface (API) Utility, 
*		abbreviated M-API.  The M-API Utility consists of subroutines 
*		which allow MODIS Science Team-supplied software to read  and 
*		write data and metadata from/to HDF files. The functionality of 
*		the M-API is defined in the MODIS Application Program Interface 
*		(API) Specification.
*		
*		cgmarin is a C function which is callable from FORTRAN. This 
*		function will call getMODISarinfo to create a group structure. In 
*		M-API, cgmarin is low-level routine called only by GMARIN. 
*		
*		In order to be callable from the FORTRAN in different platforms 
*		using function name cgmarin, this function is called ncgmarin 
*		in the actual C code. ncgmarin is redefined in mapic.h according 
*		to compiler's FORTRAN naming conventions/conversion of each 
*		platform, so that the object name of ncgmarin will always be the 
*		object name of a FORTRAN function named cgmarin.
* 
*!Input parameters:
*	 modfil	IN:	FORTRAN array of MODIS file structure.
*	 arrnm	IN:	FORTRAN ASCII string of array name.
*	 lar	IN:	Address for the number of bytes in arrnm.
*	 grpnm	IN:	FORTRAN ASCII string of group name.
*	 lgr	IN:	Address for the number of bytes in 
*	 groupname.
*	 attr	IN:	FORTRAN ASCII string of  name of the 
*			 attribute to be retrieved.
*	 lat	IN:	Address for the number of bytes in attr.
*	 dtype	IN/OUT: FORTRAN ASCII string of data type.
*	 ldt	IN:	Address for the number of bytes in dtype.
*	 nms	IN/OUT:	Address for the number of elements can be 
*	 		stored in value. Output replaces with the actual 
*	 		number of elements in value.
* 
*!Output parameters:
* 	value		OUT: pointer to the data buffer.
* 	ret		OUT: MAPIOK if successful, otherwise MFAIL
* 
* Returns:	None.
* 
* External references:
*			   MODFILLEN			 (mapic.h)
*			   MODFILE			 (mapi.h)
*			   HDf2cstring			 (hproto.h)
*			   P_ADDR 			 (mapic.h)
*			   DATATYPELENMAX		 (mapic.h)
*			   FDATATYPELENMAX		 (mapic.h)
*			   _fcdtocp			 (hdfi.h)
*			   MTYPEf2c			 (mapic.h)
*			   MTYPEc2f			 (mapic.h)
*			   HDfreespace			 (hproto.h)
*			   VOIDP			 (hdfi.h)	
*			   TXT				 (mapi.h)
*			   getMODISarinfo		 (mapi.h)
*!Revision History:
*		Qi Huang	1996/08/28
*		Version 2.1
*
*		Ring super structure and other changes make this
*		version much faster.
*
*		Qi Huang	1996/03/14
*		Version 2.0
*		Original development and testing
*
* $Log: cgmarin.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.3  1996/09/09  20:09:04  qhuang
 * Version 2.1
 *
 * Revision 1.2  1996/03/19  17:06:53  qhuang
 * *** empty log message ***
 *
 * Revision 1.1  1996/03/15  15:36:55  qhuang
 * Initial revision
 *
*
*!Team-unique header:
* 	This software is developed by the MODIS Science Data Support Team for 
*	 the National Aeronautics and Space Administration, Goddard Space 
* 	 Flight Center, under contract NAS5-32373.
* 
*!References and Credits
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes
*
!END********************************************************************
*/
{
  MODFILE	*file;
  char		*cgroupname, *carrayname, *cattrname, *fdtype, *fcdtype;
  long int n_elements;
  char		cdtype[DATATYPELENMAX+1];
  char		cfdtype[FDATATYPELENMAX+1];
  long int	c_length = DATATYPELENMAX;
  long int	f_length = FDATATYPELENMAX;
  int		i;

  /* Set file by memcpy */
  memcpy(&file,&modfil[P_ADDR],sizeof(MODFILE *));

  n_elements = *nms;

  /* Convert dtype to fdtype using HDf2cstring */
  fdtype = HDf2cstring(dtype, (intn)*ldt);

  cdtype[0] = '\0';
  MTYPEf2c(fdtype,cdtype,&c_length);

  /* Convert dtype to c pointer fcdtype by call to _fcdtocp */
  fcdtype = _fcdtocp(dtype);

  /* Convert other FORTRAN strings to C strings */ 
  cgroupname = HDf2cstring(grpnm, (intn)*lgr);
  carrayname = HDf2cstring(arrnm, (intn)*lar);
  cattrname = HDf2cstring(attr, (intn)*lat);

  *ret = getMODISarinfo(file,carrayname,cgroupname,cattrname,cdtype,
			&n_elements,value);

  if ( n_elements != 0 )
  {
    /* Note: make FORTRAN dtype) */
    if ( cdtype[0] != '\0' )
    {
      /* Convert C type cdtype to FORTRAN type cfdtype */
      MTYPEc2f(cdtype,cfdtype,&f_length);
      memset(fcdtype,' ',(size_t)*ldt);
      memcpy(fcdtype,cfdtype,strlen(cfdtype));
    }

    if ( strcmp(cdtype,TXT) == 0 )
    {
      n_elements = n_elements - 1;
      if ( n_elements < *nms )
        for ( i=(int)n_elements; i< *nms; i++)
	  *( (char *)value + i) = ' ';
    }
  }

  *nms = (intf)n_elements;

  if (cgroupname) HDfreespace((VOIDP)cgroupname);
  if (carrayname) HDfreespace((VOIDP)carrayname);
  if (cattrname) HDfreespace((VOIDP)cattrname);
  if (fdtype) HDfreespace((VOIDP)fdtype);

  return;
}

