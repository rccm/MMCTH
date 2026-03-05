#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"

void ncadmgrp(intf modfil[MODFILLEN], _fcd groupname, intf *lgr, _fcd classname,
             intf *lcl, intf *tag, intf *ref, intf *ret)
/*
!C**********************************************************************
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* addMODISgroup. This C function is only called by FORTRAN function 
* ADMGRP. This function is a M-API internal routine.
* 
*!Description:	Function cadmgrp is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API.  The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read  and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cadmgrp is a C function which is callable from FORTRAN. This 
* function will call addMODISgroup to insert an object into a group 
* structure. In M-API, cadmgrp is low-level routine called only by 
* SRMGRP. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cadmgrp, this function is called ncadmgrp 
* in the actual C code. ncadmgrp is redefined in mapic.h according 
* to compiler's FORTRAN naming conventions/conversion of each 
* platform, so that the object name of ncadmgrp will always be the 
* object name of a FORTRAN function named cadmgrp.
* 
* !Input Parameters:
*	modfil		IN:	FORTRAN array of MODIS file structure.
*	groupname	IN:	FORTRAN ASCII string of group name.
*	lgr		IN:	Address for the number of bytes in 
*				groupname.
*	classname	IN:	FORTRAN ASCII string of class name.
*	lcl		IN:	Address for the number of bytes in 
*				classname.
*	tag		IN:	The HDF tag of the object.
*	ref		IN:	The HDF reference number of the object.
* 
* !Output Parameters:
* 	ret		OUT: 	MAPIOK if successful, otherwise MFAIL
* 
* External references:
*	     MODFILLEN			(mapi.inc)
*	     MODFILE			(mapi.h)
*            HDf2cstring                (hproto.h)
*            HDfreespace                (hproto.h)
*	     VOIDP			(hdfi.h)
*	     addMODISgroup		(mapic.h)
*
* Returns:	None.
*  
* !Revision History:
*		Qi Huang	1996/10/10
*		Version 2.2
*		Ring super structure and other changes make this
*		version much faster.
*
*		Qi Huang	01/24/96
*		Version 2.0
*		Original development and testing
*
* $Log: cadmgrp.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.3  1996/09/26  17:02:25  qhuang
 * Insert bangs in prolog.
 *
 * Revision 1.2  1996/05/31  14:07:20  qhuang
 * Added external references, added prolog, removed print diagnostics, passed
 * strings converted to C format to addMODISgroup.
 *
 * Revision 1.1  1996/01/29  19:07:54  qhuang
 * Initial revision
 *
*
* !Team-unique Header:
* This software is developed by the MODIS Science Data Support Team for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
* !References and Credits:
*      Portions developed at the National Center for Supercomputing
*      Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END********************************************************************
*/
{
  MODFILE   *file;
  char      *cgroupname;
  char      *cclassname;

  /* Set file by memcpy */
  memcpy(&file,&modfil[3],sizeof(MODFILE *));

  /* Convert groupname to cgroupname and classname to cclassname by 
     using HDf2cstring. */
  cgroupname = HDf2cstring(groupname, (intn)*lgr);
  cclassname = HDf2cstring(classname, (intn)*lcl);

  if (strlen(cclassname) == 0)
  {
    HDfreespace((VOIDP)cclassname);
    cclassname = NULL;
  }

  /* Call addMODISgroup to add the object to the group. */
  *ret = addMODISgroup(file,cgroupname,cclassname,*tag,*ref);

  /* Free cgroupname and cclassname by using HDfreespace if they are not 
     NULL */
  if (cgroupname) HDfreespace((VOIDP)cgroupname);
  if (cclassname) HDfreespace((VOIDP)cclassname);

  return;
}
