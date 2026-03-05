#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include "mapic.h"

void ncpmar(intf modfil[], _fcd arrnm, intf *lar, _fcd grpnm, intf *lgr,
     	    intf start[], intf dims[], void *data, intf *ret)

/***********************************************************************
!C
*
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* putMODISarray. This C function is only called by FORTRAN function 
* PMAR. This function is a M-API internal routine.
* 
*!Description: Function cpmar is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cpmar is a wrapper which is callable from FORTRAN. This 
* function will call putMODISarray to write array data. In M-API, 
* cpmar is low-level routine which is called only by PMAR. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cpmar, this function is called ncpmar in the 
* actual C code. ncpmar is redefined in mapic.h according to 
* compiler's FORTRAN naming conventions/conversion of each 
* platform, so that the object name of ncpmar will always be the 
* object name of a FORTRAN function named cpmar.
* 
*!Input Parameters:
* modfil	IN: 	FORTRAN integer array that is used to 
*		reference the MODIS-HDF file.
* arrnm	IN:	FORTRAN character string that is the name of 
* 		the array.
* lar	IN:	FORTRAN integer address of the memory size 
* 		of arrnm.
* grpnm	IN:	FORTRAN character string which is the  name 
* 		of the data group to which the array (SDS) 
* 		belongs.
* lgr	IN:	FORTRAN integer address of the memory size 
* 		of grpnm.
* start	IN: 	FORTRAN integer array containing the array 
*		structure location to begin writing the data to 
*		the array structure.  start  must have the 
*		same number of elements as the target array 
*		has dimensions.
* dims	IN:	FORTRAN integer array describing the size of 
*		the array being written to the array 
*		structure. dims must have the same number 
*		of elements as the target array structure has 
*		dimensions and the product of the array 
*		dimensions must equal the number of 
*		elements in data.
* data	IN:	The data array.
* 
*!Output Parameters:
* 
* ret		OUT:	FORTRAN integer address of the status(MFAIL, 
* 		MAPIOK) 
* 
* Returns:	none
* 
* Externals:
*	     DATAID			(mapi.h)
*	     MODFILE			(mapi.h)
*            MAX_VAR_DIMS               (netcdf.h)
*	     MFAIL			(mapi.h)
*            HDf2cstring                (hproto.h)
*	     getMODISarrayid		(mapic.h)
*            HDfreespace                (hproto.h)
*	     SDSINFO			(mapic.h)
*            putMODISarray           	(mapi.h)
*	     VOIDP			(hdfi.h)
*
*!Revision History:
*		Qi Huang	1996/08/16
*		Version 2.1
*		Original development and testing
*		Ring super structure and other changes make this
*		version much faster.
*
* $Log: cpmar.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/08/26  21:19:03  qhuang
 * Version 2.1
 *
 * Revision 1.1  1996/08/26  21:18:23  qhuang
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
*		Qi Huang	1996/08/16
*		Version 2.1
*		Original development and testing
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
  
   MODFILE *mfile;
   char *carrnm, *cgrpnm;
   DATAID  *did;
   long int cstart[MAX_VAR_DIMS], cdims[MAX_VAR_DIMS];
   long int rank, i;
  
   /* Convert the FORTRAN character strings to C strings */
  
   carrnm = HDf2cstring(arrnm,(intn)*lar);
   cgrpnm = HDf2cstring(grpnm,(intn)*lgr);
  
   /* Set mfile by memcpy */
   memcpy(&mfile,&modfil[P_ADDR],sizeof(MODFILE *));

   did = getMODISarrayid(mfile,carrnm,cgrpnm);
   if ( did == NULL )
   {
     if ( carrnm ) HDfreespace((VOIDP)carrnm);
     if ( cgrpnm ) HDfreespace((VOIDP)cgrpnm);
     *ret = MFAIL;
     return;
   }

   rank = ((SDSINFO *)did->info)->rank;

   /* first dimension in FORTRAN becomes last dimension in C */
   for (i=0; i<rank; i++)
   {
     cstart[i] = start[rank-i-1];
     cdims[i] = dims[rank-i-1];
   }

   *ret = putMODISarray(mfile,carrnm,cgrpnm,cstart,cdims,data);

   if (carrnm) HDfreespace((VOIDP)carrnm);
   if (cgrpnm) HDfreespace((VOIDP)cgrpnm);

   return;
 }
