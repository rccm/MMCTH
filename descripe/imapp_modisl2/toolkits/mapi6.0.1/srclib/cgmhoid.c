#include <stdlib.h>
#include "mapic.h"

void ncgmhoid(intf *modfil, _fcd name, intf *lna, _fcd group, intf *lgr,
                intf *type, _fcd access, intf *lac, intf *ret)

/***********************************************************************
!C
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* getMODIShobjid. This C function is only called by FORTRAN function 
* GMHOID. This function is a M-API internal routine.
* 
*!Description: Function cgmhoid is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read  and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cgmhoid is a C function which is callable from FORTRAN. This 
* function will call getMODIShobjid to obtain the HDF id for a data 
* object in the MODIS file. In M-API, cgmhoid is low-level routine 
* which is called only by GMHOID. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cgmhoid, this function is called ncgmhoid 
* in the actual C code. ncgmhoid is redefined in mapic.h according 
* to compiler's FORTRAN naming conventions/conversion of each 
* platform, so that the object name of ncgmhoid will always the 
* object name of a FORTRAN function named cgmhoid.
* 
*!Input parameters:
* modfil	IN: 	Array  of FORTRAN integer array that is used 
* 		to reference the MODIS-HDF file.
* name		IN:	FORTRAN character string of the object name.
* lna		IN:	FORTRAN integer address of the memory size 
* 		of name.
* group		IN:	ASCII string name of the data group to which 
* 		the object belongs.
* lgr		IN:	FORTRAN integer address of the memory size 
* 		of group.
* type		IN:	FORTRAN integer address of the object type.
* access	IN:	FORTRAN character string of the access type.
* lac		IN:	FORTRAN inteher address of the memory size 
* 		of access
* 
*!Output parameters:
* ret		OUT:	FORTRAN integer address of the return value 
* 		(id or MFAIL)
* 
* Returns:	none
* 
* Externals:
*	     MODFILE			(mapi.h)
*            HDf2cstring                (hproto.h)
*	     getMODIShobjid		(mapi.h)
*            HDfreespace                (hproto.h)
*	     VOIDP			(hdfi.h)
*	     P_ADDR			(mapic.h)
*
*!Revision History:
*		Qi Huang	1996/11/25
*		Version 2.2
*		Ring super structure and other changes make this
*		version much faster.
*
* $Log: cgmhoid.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
*
*!Team-unique Header:
*
*    This software is developed by the MODIS Science Data Support Team 
*    for the National Aeronautics and Space Administration, 
*    Goddard Space Flight Center, under contract NAS5-32373.
*
*!References and Credits:
*
*    Portions developed at the National Center for Supercomputing
*    Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes:
*
*!END*****************************************************************
*/

 { 
   /* declare local variables */
   MODFILE *file;
   char *cname, *cgroup, *caccess;
  
   /* Convert the FORTRAN character strings to C strings */
   cname = HDf2cstring(name,(intn)*lna);
   cgroup = HDf2cstring(group,(intn)*lgr);
   caccess = HDf2cstring(access,(intn)*lac);
  
   /* Set file by memcpy */
   memcpy(&file,&modfil[P_ADDR],sizeof(MODFILE *));

   *ret = getMODIShobjid(file,cname,cgroup,*type,caccess);
   
   /* Free all allocated memory */ 
   if ( cname ) HDfreespace((VOIDP)cname);
   if ( cgroup ) HDfreespace((VOIDP)cgroup);
   if ( caccess ) HDfreespace((VOIDP)caccess);
    
   return;
 }
