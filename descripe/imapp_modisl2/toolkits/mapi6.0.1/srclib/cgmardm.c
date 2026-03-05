#include <stdio.h> 
#include <stdlib.h>
#include "mapic.h"

void ncgmardm(intf *modfil,_fcd arrnm,intf *lar,_fcd grpnm,intf *lgr,
              _fcd dtype,intf *ldt,intf *rank,intf dims[],intf *ret)
/*
!C**********************************************************************
*
*!Description:
*
*   Function cgmardm is part of a larger software system called the 
*   MODIS Applications Programming Interface (API) Utility, abbreviated 
*   M-API.  The M-API Utility consists of subroutines which allow MODIS 
*   Science Team-supplied software to read and write data and metadata 
*   from/to HDF files.  The functionality of the M-API is defined in 
*   the MODIS Application Program Interface Specification.
*
*   cgmardm is a C function which is callable from FORTRAN. This 
*   function will call getMODISardims to read data from an array. In 
*   M-API, cgmardm is low-level routine which is called only by GMARDM.  
*   This function is a M-API internal routine.
*
*   In order to be callable from the FORTRAN in different platforms 
*   using function name cgmardm, this function is called ncgmardm in the
*   actual C code.  ncgmardm is redefined in mapic.h according to 
*   compiler's FORTRAN naming conventions of each platform, so that the 
*   object name of ncgmardm will always be the object name of a FORTRAN 
*   function named cgmardm.
*
* !Input Parameters:
*
*   modfil IN:  FORTRAN integer array used to reference the MODIS-HDF file.
*   arrnm  IN:  FORTRAN character string, the name of the array.
*   lar    IN:  FORTRAN integer address of the memory size of arrnm.
*   grpnm  IN:  FORTRAN character string, the name of the data group to 
*          which the array (SDS) belongs.
*   lgr    IN:  FORTRAN integer address of the memory size of group.
*   dims   IN:  FORTRAN integer array to indicate the size of each dimension.
*
* !Output Parameters:
*
*   dtype  OUT:  FORTRAN string, the data type of the array structure. 
*   rank   IN/OUT:  The address of a FORTRAN integer. 
*          Input with the number of elements in the dims array, and 
*          output replaced with actually number of dimension in the array.
*          If a function error, the value is set to zero.  If the dims is 
*          not large enough, this value will be set the actually rank in 
*          the array in the HDF while the function will return MFAIL.
*   dims   IN:  FORTRAN integer array to place the size for each dimension. 
*   ret    OUT:  FORTRAN integer address of the status(MFAIL, MAPIOK) 
*
* Returns:	none
*
* External:
*          MODFILE                 (mapi.h)
*          MFAIL                   (mapi.h)
*          FDATATYPELENMAX	   (mapic.h)
*	   MAX_VAR_DIMS		   (netcdf.h)
*	   HDf2cstring		   (hproto.h)
*	   P_ADDR		   (mapic.h)
*	   _fcdtocp		   (hdfi.h)
*	   getMODISardims	   (mapi.h)
*	   MTYPEc2f		   (mapic.h)
*	   HDfreespace		   (hproto.h)
*	   VOIDP		   (hdfi.h)
*
* !Revision History:
*		Qi Huang	1996/09/03
*		Version 2.1
*
*		Ring super structure and other changes make
*		this version much faster.
*
* $Log: cgmardm.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/09/18  16:45:05  qhuang
 * Version 2.1
 *
 * Revision 1.1  1996/09/18  16:44:52  qhuang
 * Initial revision
 *
*
* !Team-unique Header:
*
*   This software is developed by the MODIS Science Data Support Team 
*   for the National Aeronautics and Space Administration, 
*   Goddard Space Flight Center, under contract NAS5-32373.
*
* !References and Credits:
*
*   Written by   Vicky Lin        02/06/96
*   Research and Data Systems Corporation
*   SAIC/GSC MODIS SCIENCE DATA SUPPORT OFFICE
*   7501 FORBES BLVD, SEABROOK MD 20706
*
*   vlin@ltpmail.gsfc.nasa.gov
*
*   Portions developed at the National Center for Supercomputing
*   Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END*****************************************************************
*/

{
  /* declare local variables */
  MODFILE	*mfile;
  char *carrnm, *cgrpnm, *fcdtype;
  char fdtype[FDATATYPELENMAX + 1];
  long int dimsizes[MAX_VAR_DIMS];
  int i;
  long int crank;
  long int f_length = FDATATYPELENMAX + 1;
  
  /* convert FORTRAN input strings to C strings */
  carrnm = HDf2cstring(arrnm,(intn)*lar);
  cgrpnm = HDf2cstring(grpnm,(intn)*lgr);

  /* Set mfile by memcpy */
  memcpy(&mfile,&modfil[P_ADDR],sizeof(MODFILE *));
       
  /* convert FORTRAN string to C string */
  fcdtype = _fcdtocp(dtype);

  crank = *rank;

  *ret = getMODISardims(mfile,carrnm,cgrpnm,fcdtype,&crank,dimsizes);
  *rank = (intf)crank;

  if (crank != 0L)	/* crank will not be zero if call to 
			getMODISardims succeed or partially failed. */
  {
    /* Convert C type fcdtype to FORTRAN type fdtype */
    MTYPEc2f(fcdtype,fdtype,&f_length);    

    memset(fcdtype,' ',(size_t)*ldt);     /* set the memory of 
                                           fcdtype to blank */
    memcpy(fcdtype,fdtype,strlen(fdtype)); 
                                          /* memcpy the content of 
                                             fdtype to fcdtype */
    if (*ret != MFAIL) 
    {
      for (i=0; i<crank; i++)
        dims[i] = (intf)dimsizes[crank-i-1];
    }

  } /* end of if */

  if (carrnm) HDfreespace((VOIDP)carrnm);
  if (cgrpnm) HDfreespace((VOIDP)cgrpnm);
  return;

 }

