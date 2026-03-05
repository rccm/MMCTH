#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include "mapic.h"

void ncpmdmin(intf modfil[], _fcd arrnm, intf *lar, _fcd grpnm, intf *lgr,
	      intf *dim, _fcd attrnm, intf *lat, _fcd dtype, 
	      intf *ldt, intf *nms, void* value, intf *ret)
/*
!C**********************************************************************
* 
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* putMODISdiminfo. This C function is only called by FORTRAN function 
* PMDMIN. This function is a M-API internal routine.
* 
*!Description: Function cpmdmin is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cpmdmin is a wrapper which is callable from FORTRAN. This 
* function will call putMODISdiminfo to write a dimensional 
* attribute. In M-API, cpmdmin is low-level routine which is called 
* only by PMAR. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cpmdmin, this function is called ncpmdmin 
* in the actual C code. ncpmdmin is redefined in mapic.h according 
* to compiler's FORTRAN naming conventions/conversion of each 
* platform, so that the object name of ncpmdmin will always be the 
* object name of a FORTRAN function named cpmdmin.
* 
* !Input Parameters:
* modfil	IN: 	FORTRAN integer array that is used to 
* 		reference the MODIS-HDF file.
* arrnm		IN:	FORTRAN character string that is the name of 
* 		the array.
* lar		IN:	FORTRAN integer address of the memory size 
* 		of arrnm.
* grpnm		IN:	FORTRAN character string which is the  name 
*		of the data group to which the array (SDS) 
*		belongs.
* lgr		IN:	FORTRAN integer address of the memory size 
* 		of grpnm.
* dim		IN:	The dimension to which the attribute will be 
* 		attached (0-based).
* attrnm	IN:	Name to assign the attribute.  Provided macros 
*		for accepted MODIS file attribute names are 
*		listed in Appendix A of the UserÕs Guide, 
*		MODIS API Supplied Constants.
* lat		IN:	FORTRAN address of the string length of 
* 		attrnm
* dtype		IN:	Data Type of the value.
* ldt		IN:	FORTRAN address of the string length of dtype
* nms		IN:	Number of attribute values in value.
* value		IN:	The data to store in the attribute.  If the 
*		attribute already exists, the value will be 
*		updated.
* 
* !Output Parameters:
* 
* ret		OUT:	FORTRAN integer address of the status(MFAIL, 
* 		MAPIOK) 
* 
* Returns:	none
* 
* External references:
*			   MODFILE			 (mapi.h)
*			   DATAID			 (mapi.h)
*			   P_ADDR			 (mapic.h)
*			   HDf2cstring			 (hproto.h)
*			   DATATYPELENMAX		 (mapic.h)
*			   getMODISarrayid		 (mapic.h)
*			   MFAIL			 (mapi.h)
*			   SDSINFO			 (mapic.h)
*			   MTYPEf2c			 (mapic.h)
*			   HDfreespace			 (hproto.h)
*			   putMODISdiminfo		 (mapi.h)
*
*!Revision History:
*		Qi Huang	1996/08/22
*		Version 2.1
*		Original development and testing
*
*		Ring super structure and other changes make this
*		version much faster.
*
* $Log: cpmdmin.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/08/30  19:32:53  qhuang
 * Version 2.1
 *
 * Revision 1.1  1996/08/30  19:28:49  qhuang
 * Initial revision
 *
* 
*!Team-unique header:
*	 This software is developed by the MODIS Science Data SupportTeam for 
*	 the National Aeronautics and Space Administration, Goddard Space 
*	 Flight Center, under contract NAS5-32373.
* 
*!References and Credits
*		Qi Huang	1996/08/22
*		Version 2.1
*		Original development and testing
*
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes
*
!END********************************************************************
*/
{
  MODFILE	*mfile;
  char		*carrnm, *cgrpnm, *cattrnm, *cfdtype;
  DATAID	*did;
  char		cdtype[DATATYPELENMAX];
  long int	rank, cdim;
  long int	c_length = DATATYPELENMAX;

  /* Convert FORTRAN strings to C strings using HDf2cstring. */ 
  cgrpnm = HDf2cstring(grpnm, (intn)*lgr);
  carrnm = HDf2cstring(arrnm, (intn)*lar);
  cattrnm = HDf2cstring(attrnm, (intn)*lat);
  cfdtype = HDf2cstring(dtype, (intn)*ldt);

  /* Set file by memcpy. */
  memcpy(&mfile,&modfil[P_ADDR],sizeof(MODFILE *));

  did = getMODISarrayid(mfile,carrnm,cgrpnm);
  if ( did == NULL )
  {
    if (carrnm) HDfreespace(carrnm);
    if (cgrpnm) HDfreespace(cgrpnm);
    if (cattrnm) HDfreespace(cattrnm);
    if (cfdtype) HDfreespace(cfdtype);
    *ret = MFAIL;
    return;
  }

  /* Convert FORTRAN dimension to C dimension */
  rank = ((SDSINFO *)did->info)->rank;
  cdim = rank - *dim -1;

  /* Convert FORTRAN data type to C data type */
  cdtype[0] = '\0';
  MTYPEf2c(cfdtype,cdtype,&c_length); 

  *ret = putMODISdiminfo(mfile,carrnm,cgrpnm,cdim,cattrnm,cdtype,*nms,
                   	value);

  /* Free the C strings allocated by HDf2cstring if they are not null by
     using HDfreespace. */
  if (carrnm) HDfreespace(carrnm);
  if (cgrpnm) HDfreespace(cgrpnm);
  if (cattrnm) HDfreespace(cattrnm);
  if (cfdtype) HDfreespace(cfdtype);

  return;
}
