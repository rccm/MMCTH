#include <stdio.h> 
#include <stdlib.h>
#include "mapic.h"

void ncclmfil(intf modfil[MODFILLEN], intf *ret)
/*
!C**********************************************************************
* 
*!Purpose:	A wrapping function interfacing between C closeMODISfile
* and FORTRAN CLMFIL. This C function is only called by FORTRAN function 
* CLMFIL. This function is a M-API internal routine.
* 
*!Description:	Function cclmfil is part of a larger software system called 
* the MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API.  The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read  and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cclmfil is a C function which is callable from FORTRAN. This 
* function will call closeMODISfile to terminate the access of a 
* MODIS-HDF file. In M-API, cclmfil is a low-level routine called 
* only by CLMFIL. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cclmfil, this function is called ncclmfil in 
* the actual C code. ncclmfil is redefined in mapic.h according to 
* compiler's FORTRAN naming conventions/conversion of each 
* platform, so that the object name of ncclmfil will always the 
* object name of a FORTRAN function named cclmfil.
* 
* !Input Parameters:
* 	modfil	IN/OUT: FORTRAN MODIS file array. If all successful, 
* 			the elements of the array will be filled with 
* 			zero after closing the MODIS file.
* !Output Parameters:
* 	ret	OUT: 	if function is successful, ret is set to 0, 
* 			otherwise -1
* 
* Returns:	None.
* 
* External references:
*			MFAIL		(mapi.h)
*			MAPIOK		(mapi.h)
*			P_HDFID		(mapic.h)
*			P_SDID		(mapic.h)
*			P_ACCESS	(mapic.h)
*			P_ADDR		(mapic.h)
*			closeMODISfile  (mapi.h)
* !Revision History:
*		Qi Huang	1996/07/24
*		Version 2.1
*		Ring super structure and other changes make this
*		version much faster.
*
*		Qi Huang	01/02/96
*		Version 2.0
*		Original development and testing
*
* $Log: cclmfil.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/08/08  16:05:35  qhuang
 * Changed modfil[3] to modfil[P_ADDR].
 *
 * Revision 1.1  1996/08/08  16:04:57  qhuang
 * Initial revision
 *
 *
 * Revision 1.2  1996/05/31  14:36:58  qhuang
 * Put procedure header above prolog.
 *
 * Revision 1.1  1996/01/02  18:53:14  qhuang
 * Initial revision
 *
* 	
* !Team-unique Header:
* This software is developed by the MODIS Science Data Support Team for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
* !References and Credits:
*
* !Design Notes:
*
!END********************************************************************
*/
{
  MODFILE *mfile;           /* pointer to MODFILE */

  *ret = MFAIL;

  /* Set mfile by memcpy */
  memcpy(&mfile,&modfil[P_ADDR],sizeof(MODFILE *));

  /* Call closeMODISfile to close the file. */
  *ret = closeMODISfile(&mfile);

  if (*ret == MAPIOK)
  {
    modfil[P_HDFID] = 0;
    modfil[P_SDID]  = 0;
    modfil[P_ACCESS] = 0;
    modfil[P_ADDR] = 0;
  }

  return;
}
