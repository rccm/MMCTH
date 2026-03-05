#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"
#include <ctype.h>

void  ncgmdnam(intf *modfil, _fcd arrnm, intf *lar, _fcd grpnm, intf *lgr,  intf *dim, _fcd dname,  intf *ldn, intf *ret)
/*
!C**********************************************************************
* 
* !Purpose: A wrapping function interfacing between C and FORTRAN for
*           getMODISdimname. This C function is only called by FORTRAN
*           function GMDNAM. This function is a M-API internal routine.
*
*!Description: Function cgmdnam is part of a larger software system called
*              the MODIS Applications Programming Interface (API) Utility,
*              abbreviated M-API. The M-API Utility consists of subroutines
*              which allow MODIS Science Team-supplied software to read and
*              write data and metadata from/to HDF files. The functionality
*              of the M-API is defined in the MODIS Application Program
*              Interface (API) Specification.
*
*              cgmdnam is a wrapper which is callable from FORTRAN. This
*              function will call getMODISdimname to read a dimension name.
*              In M-API, cgmdnam is low-level routine which is called only
*              by GMDNAM.
*
*              In order to be callable from the FORTRAN in different
*              platforms using function name cgmdnam, this function is called 
*              ncgmdnam in the actual C code. ncgmdnam is redefined in mapic.h
*              according to compiler's FORTRAN naming conventions/conversion
*              of each platform, so that the object name of ncgmdnam will
*              always be the object name of a FORTRAN function named cgmdnam.
*
*!Input Parameters:
*
*       modfil  IN: FORTRAN integer array that is used to reference the
*                   MODIS-HDF file.
*       arrnm   IN: FORTRAN character string that is the name of the array.
*       lar     IN: FORTRAN integer address of the memory size of arrnm.
*       grpnm   IN: FORTRAN character string which is the  name of the data
*                   group to which the array (SDS) belongs.
*       lgr     IN: FORTRAN integer address of the memory size of grpnm.
*       dim     IN: FORTRAN integer address specifying the dimension.
*       ldn     IN: FORTRAN integer address of the memory size of dname.
*
*!Output Parameters:
        name   OUT: FORTRAN string containing the dimension name.
        ret    OUT: FORTRAN integer address of the status(MFAIL, MAPIOK)
*
* Returns: none
*
* External references:
*
*    SDSINFO                                  (mapic.h)
*    HDFreespace                              (hdf.h)
*    MFAIL                                    (mapi.h)
*    MAPIOK                                   (mapi.h)
*    MODFIL                                   (mapi.h)
*    DATAID                                   (mapic.h)
*    getMODISarrayid                          (mapic.h)
*    getMODISdimname                          (mapic.h)
*    fcdtocp                                  (hdf.h)
*    HDf2cstring                              (hdf.h)
*    HDc2fstring                              (hdf.h)
*    P_ADDR                                   (mpaic.h)
*
*!Revision History:
*   $LOG$
* 
*!Team-unique Header:
* This software is developed by the MODIS Science Data SupportTeam for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
*!References and Credits:
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes:
*
!END********************************************************************
*/
{
    MODFILE  *mfile;
    char *carrnm, *cgrpnm, *cdname;
    DATAID *did;
    long int rank, cdim;
/*
Convert the FORTRAN character strings arrnm and  grpnm to C character strings carrnm and cgrpnm by using HDf2cstring. 
*/
    carrnm = HDf2cstring(arrnm,(intn)*lar);
    cgrpnm = HDf2cstring(grpnm,(intn)*lgr);

/* Convert dname to c pointer cdname by _fcdtocp */
    cdname  = (char *)_fcdtocp(dname);

/* Set mfile by memcpy (memcpy(&mfile, &modfil[3], sizeof(MODFILE *)) ) */
    memcpy(&mfile, &modfil[P_ADDR], sizeof(MODFILE *));
	
/* Set did  to the return value of getMODISarrayid  */
   did = getMODISarrayid( mfile, carrnm, cgrpnm);

   if( did == NULL ) {
       if(carrnm) HDfreespace( (VOIDP)carrnm);
       if(cgrpnm) HDfreespace( (VOIDP)cgrpnm);
       *ret = MFAIL;
       return;
   }
   
   rank = ((SDSINFO *)did->info)->rank;
   cdim = (rank - *dim) - 1;

   /* Call getMODISdimname to retrieve the dimensional name.
    The return status of getMODISdimname should be assigned to the
     return variable *ret
     (Note: using cdim, carrnm, cgrpnm, mfile, and cdname in the call).
   */
   *ret = getMODISdimname(mfile, carrnm, cgrpnm, cdim, cdname);

   /* HDfreespace carrnm and cgrpnm                          */

   if (carrnm) HDfreespace((VOIDP)carrnm);
   if (cgrpnm) HDfreespace((VOIDP)cgrpnm);
 
   if (*ret != MFAIL) {
       HDc2fstr(cdname, (intn) *ldn);    
   }
   
  return;
}
