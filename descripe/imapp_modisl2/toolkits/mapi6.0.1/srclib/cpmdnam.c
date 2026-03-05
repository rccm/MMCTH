#include <stdio.h> 
#include <stdlib.h>
#include "mapic.h"

void ncpmdnam(intf *modfil, _fcd arrnm, intf *lar, _fcd grpnm, intf *lgr,
		intf *dim, _fcd dname, intf *ldn, intf *ret)
/*
!C**********************************************************************
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* putMODISdimname. This C function is only called by FORTRAN function 
* PMDNAM. This function is a M-API internal routine.
* 
*!Description: Function cpmdnam is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cpmdnam is a C function which is callable from FORTRAN. This 
* function will call putMODISdimname to read a dimension name. 
* In M-API, cpmdnam is low-level routine which is called only by 
* PMDNAM. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cpmdnam, this function is called ncpmdnam 
* in the actual C code. ncpmdnam is redefined in mapic.h according 
* to compiler's FORTRAN naming conventions/conversion of each 
* platform, so that the object name of ncpmdnam will always be the 
* object name of a FORTRAN function named cpmdnam.
* 
*!Input Parameters:
* modfil	IN: 	FORTRAN integer array that is used to 
* 		reference the MODIS-HDF file.
* arrnm	IN:	FORTRAN character string that is the name of 
* 		the array.
* lar	IN:	FORTRAN integer address of the memory size 
* 		of arrnm.
* grpnm	IN:	FORTRAN character string which is the  name 
* 		of the data group to which the array (SDS) 
* 		belongs.
* lgr	IN:	FORTRAN integer address of the memory size 
* 		of group.
* 	IN:	FORTRAN integer array to indicate the start 
* 		location to retrieve.
* dim	IN:	FORTRAN integer address specifying the 
* 		dimension.
* dname	IN: 	FORTRAN string containing the dimension 
* 		name.
* 
*!Output Parameters:
* 
* ret		OUT:	FORTRAN integer address of the status(MFAIL, 
* 		MAPIOK) 
* 
* Returns:	none
* External Reference:
*	     DATAID			(mapi.h)
*	     MODFILE			(mapi.h)
*	     MFAIL			(mapi.h)
*            HDf2cstring                (hproto.h)
*	     getMODISarrayid		(mapic.h)
*            HDfreespace                (hproto.h)
*	     SDSINFO			(mapic.h)
*            putMODISdimname           	(mapi.h)
*	     VOIDP			(hdfi.h)
*
*!Revision History:
*		Qi Huang	1996/08/19
*		Version 2.1
*		Original development and testing
*
*		Ring super structure and other changes make
*		this version much faster.
*
* $Log: cpmdnam.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/08/26  21:05:59  qhuang
 * Version 2.1
 *
 * Revision 1.1  1996/08/26  21:05:41  qhuang
 * Initial revision
 *
*
*!Team-unique Header:
*    This software is developed by the MODIS Science Data Support Team for 
*    the National Aeronautics and Space Administration, Goddard Space 
*    Flight Center, under contract NAS5-32373.
*
*!References and Credits:
*		Qi Huang	1996/08/19
*		Version 2.1
*		Original development and testing
*
*    HDF portions developed at the National Center for Supercomputing
*    Applications at the University of Illinois at Urbana-Champaign.
*
*!Design Notes:
*
!END********************************************************************
*/
{
  MODFILE  *mfile;            /* pointer to MODFILE structure */
  char	   *carrnm, *cgrpnm, *cdname;
  DATAID   *did;
  long int cdim;

  /* Convert the FORTRAN character strings to C characater
     strings using HDf2cstring. */
  carrnm = HDf2cstring(arrnm, (intn)*lar);
  cgrpnm = HDf2cstring(grpnm, (intn)*lgr);
  cdname = HDf2cstring(dname, (intn)*ldn);
  
  /* Set mfile by memcpy */
  memcpy(&mfile,&modfil[P_ADDR],sizeof(MODFILE *));

  did = getMODISarrayid(mfile,carrnm,cgrpnm);
  if ( did == NULL )
  {
    *ret = MFAIL;
    if (carrnm) HDfreespace((VOIDP)carrnm);
    if (cgrpnm) HDfreespace((VOIDP)cgrpnm);
    if (cdname) HDfreespace((VOIDP)cdname);
    return;
  }

  cdim = ((SDSINFO *)did->info)->rank - *dim - 1;

  /* Call putMODISdimname to write the dimensional name. */
  *ret = putMODISdimname(mfile,carrnm,cgrpnm,cdim,cdname);

  if (carrnm) HDfreespace((VOIDP)carrnm);
  if (cgrpnm) HDfreespace((VOIDP)cgrpnm);
  if (cdname) HDfreespace((VOIDP)cdname);
  
  return;
}
