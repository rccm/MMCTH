#include <stdio.h> 
#include <stdlib.h>
#include "mapic.h"

void ncopmfil(_fcd fname, intf *lfn, _fcd access, intf *lac,
	      intf modfil[MODFILLEN], intf *ret)
/*
!C**********************************************************************
*!Purpose: A wrapping function interfacing between C openMODISfile and 
*          FORTRAN OPMFIL. This C function is only called by FORTRAN function 
*          OPMFIL. This function is a M-API internal routine.
* 
*!Description: Function copmfil is part of a larger software system called 
*              the MODIS Applications Programming Interface (API) Utility, 
*              abbreviated M-API. The M-API Utility consists of subroutines 
*              which allow MODIS Science Team-supplied software to read  and 
*              write data and metadata from/to HDF files. The functionality of 
*              the M-API is defined in the MODIS Application Program Interface 
*              (API) Specification.
* 
*              copmfil is a C function which is callable from FORTRAN. This 
*              function will call openMODISfile to open/create a MODIS-HDF file. 
*              In M-API, copmfil is low-level routine which is called only by 
*              CRMTBL. 
* 
*              In order to be callable from the FORTRAN in different platforms 
*              using function name copmfil, this function is called ncopmfil in 
*              the actual C code. ncopmfil is redefined in mapic.h according to 
*              compiler's FORTRAN naming conventions/conversion of each 
*              platform, so that the object name of ncopmfil will always the 
*              object name of a FORTRAN function named copmfil.
* 
* !Input Parameters:
*              fname	IN:	FORTRAN complete path and filename for the 
*		                file to be opened.
*              lfn	IN:	Number of bytes in fname.
*              access	IN:	FORTRAN character string of access mode.
* 
*                               Permitted access mode:
*                              	"r"	Open for read only.
*                         	"w"     Create for read/write.
*                         	"a"     Open for read/write.
* 
*	       lac	IN:	Number of bytes in access.
* 
* !Output Parameters:
*              modfil	OUT:	Array that is used to reference the file in all 
*            		        other MODIS API routines in FORTRAN.  The 
*              			array will return all zeroes if an error occurs.
*              ret	OUT: 	if function is successful, ret is set to 0, otherwise -1
* 
* Returns:	None.
* 
* External references:
*		  	   MODFILLEN			 (mapic.h)	
*			   MODFILE			 (mapi.h)
*                          MFAIL              		 (mapi.h)
*                          MAPIOK             		 (mapi.h)
*			   P_SDID			 (mapic.h)
*			   P_HDFID			 (mapic.h)
*			   P_ACCESS			 (mapic.h)
*			   P_ADDR			 (mapic.h)
*			   HDf2cstring			 (hproto.h)
*			   openMODISfile		 (mapi.h)
*			   HDfreespace			 (hproto.h)
*			   VOIDP			 (hdfi.h)	
*
* !Revision History:
*		Arvind Solanki	1999/02/04
*		Fixed typos in prolog.
*
*		Qi Huang	1996/07/24
*		Version 2.1
*		Ring super structure and other changes make
*		this version much faster.
*
*		Qi Huang	04/10/96
*		Version 2.0
*		Original development and testing
*
* $Log: copmfil.c,v $
* Revision 1.2  1999/04/26 16:56:29  jayshree
* reorganizing mapi RCS
*
 * Revision 1.4  1996/08/09  20:30:27  qhuang
 * Cast some variables in assignment statements.
 *
 * Revision 1.3  1996/08/08  16:04:09  qhuang
 * Changed modfil[3] to modfil[P_ADDR].
 *
 * Revision 1.2  1996/05/31  15:48:56  qhuang
 * Added NCSA credit in prolog.
 *
 * Revision 1.1  1996/05/31  15:26:13  qhuang
 * Initial revision
 *
 * Revision 1.1  1996/05/31  15:26:13  qhuang
 * Initial revision
 *
* 
* !Team-unique Header:
*    This software is developed by the MODIS Science Data Support Team for 
*    the National Aeronautics and Space Administration, Goddard Space 
*    Flight Center, under contract NAS5-32373.
*
* !References and Credits:
*    HDF portions developed at the National Center for Supercomputing
*    Applications at the University of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END********************************************************************
*/
{
  char     *cfname;           /* pointer to C MODIS HDF file name */
  char     *caccess;          /* pointer to C access mode */
  MODFILE  *mfile;            /* pointer to MODFILE structure */

  /* Zero modfil. */
  modfil[P_SDID] = 0;
  modfil[P_HDFID] = 0;
  modfil[P_ACCESS] = 0;
  modfil[P_ADDR] = 0;
 
  /* Convert fname to cfname and access to caccess using HDf2cstring. */
  cfname = HDf2cstring(fname, (intn)*lfn);
  caccess = HDf2cstring(access, (intn)*lac);
  
  /* Call openMODISfile to obtain the MODFILE structure. */
  mfile = openMODISfile(cfname,caccess);

  /* Free caccess and cfname using HDfreespace if they are not NULL */
  if (caccess)
    HDfreespace((VOIDP)caccess);

  if (cfname)
    HDfreespace((VOIDP)cfname);

  if ( mfile == NULL )
    *ret = MFAIL;  
  else
  {
    *ret = MAPIOK; 
    /* Set the elements in modfil according to values in mfile. */
    modfil[P_SDID] = (intf)mfile->sd_id;
    modfil[P_HDFID] = (intf)mfile->hdf_id;
    modfil[P_ACCESS] = (intf)mfile->access;
    memcpy(&modfil[P_ADDR],&mfile,sizeof(MODFILE *));
  }
  
  return;
}

