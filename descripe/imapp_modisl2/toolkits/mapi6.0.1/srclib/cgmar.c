#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"

void  ncgmar(intf modfil[], _fcd arrnm, intf *lar, _fcd grpnm, intf *lgr,  intf start[],  intf dims[],  void *data, intf *ret)
/**************************************************************************************
!C
*
* !Purpose: A wrapping function interfacing between C and FORTRAN for getMODISarray.
*           This C function is only called by FORTRAN function GMAR. This function
*           is a M-API internal routine.
*
*!Description: Function cgmar is part of a larger software system called the MODIS
*               Applications Programming Interface (API) Utility, abbreviated M-API.
*               The M-API Utility consists of subroutines which allow MODIS Science
*               Team-supplied software to read and write data and metadata from/to
*               HDF files. The functionality of the M-API is defined in the MODIS
*               Application Program Interface (API) Specification.
*
*               cgmar is a wrapper which is callable from FORTRAN. This function will
*               call getMODISarray to read array data. In M-API, cgmar is low-level
*               routine which is called only by GMAR. 
*
*               In order to be callable from the FORTRAN in different platforms using
*               function name cgmar, this function is called ncgmar in the actual C
*               code. ncgmar is redefined in mapic.h according to compiler's FORTRAN
*               naming conventions/conversion of each platform, so that the object
*               name of ncgmar will always be the object name of a FORTRAN function
*               named cgmar.
*
*!Input Parameters:
*
*               modfil	IN: FORTRAN integer array that is used to reference the
*                           MODIS-HDF file.
*               arrnm	IN: FORTRAN character string that is the name of the array.
*               lar	IN: FORTRAN integer address of the memory size of arrnm.
*               grpnm	IN: FORTRAN character string which is the  name of the data
*                           group to which the array (SDS) belongs.
*               lgr	IN: FORTRAN integer address of the memory size of grpnm.
*               start	IN: FORTRAN integer array containing the array structure
*                       location to begin reading the data from the array structure.
*                       start must have the same number of elements as the target
*                       array has dimensions.
*               dims	IN: FORTRAN integer array describing the size of the array
*                       being retrieved from ther array structure. dims must have
*                       the same number of elements as the target array structure
*                       has dimensions and the product of the array dimensions must
*                       equal the number of elements in data.
*
*!Output Parameters:
*               data	OUT: Multi-dimensional data array.
*
*               ret	OUT: FORTRAN integer address of the status(MFAIL, MAPIOK) 
*
* !Returns:	none
* 
* !External references:
*                      SDSINFO                  (mapic.h)
*                      MODFILE                  (mapi.h)
*                      DATAID                   (mapic.h)
*                      MAX_VAR_DIMS             (netcdf.h)
*                      P_ADDR                   (mapic.h)
*                      HDf2string               (hproto.h)
*                      HDfreespace              (hproto.h)
*                      getMODISarrayid          (mapic.h)
*                      VOIDP                    (hdfi.h)
*                      getMODISarray            (mapic.h)
*                      MFAIL                    (mapic.h)
*                      
*
*!Revision History:
*               Frederick J. Shaw               Sept 9, 1996
*               Version 2.1
*               Orginal Development and Testing
*  $LOG$
*
*!Team-unique Header:
*  This software is developed by the MODIS Science Data Support Team for the
*  National Aeronautics and Space Administration, Goddard Space Flight Center,
*  under contract NAS5-32373.
*
*  Portions developed at the National Center for Supercomputing Applications
*  at the Univ. of Illinois at Urbana-Champaign.
*
*!References and Credits:
*   
*
*!Design Notes:
*
* !END*************************************************************
*/

{
    MODFILE  *mfile;
    char *carrnm, *cgrpnm;
    DATAID *did;
    long int cstart[MAX_VAR_DIMS], cdims[MAX_VAR_DIMS];
    long int rank, i;

/*
  Convert the FORTRAN character strings arrnm and grpnm to C character
  strings carrnm and cgrpnm by using HDf2cstring. 
*/

    carrnm = HDf2cstring(arrnm,(intn)*lar);
    cgrpnm = HDf2cstring(grpnm,(intn)*lgr);

    memcpy(&mfile, &modfil[P_ADDR], sizeof(MODFILE *));
	
    did = getMODISarrayid( mfile, carrnm, cgrpnm);

    if( did == NULL ) {
        if (carrnm) HDfreespace((VOIDP)carrnm);
        if (cgrpnm) HDfreespace((VOIDP)cgrpnm);
        *ret = MFAIL;
        return;
    }

   rank = ((SDSINFO *)did->info)->rank;

/* first dimension in FORTRAN becomes last dimension in C */
   for( i = 0; i <= rank - 1; i++ ) { 
       cstart[i] = start[(rank-i)-1];
       cdims[i] = dims[(rank-i)-1];
   }

   *ret = getMODISarray(mfile, carrnm, cgrpnm, cstart, cdims, data);

   if (carrnm) HDfreespace((VOIDP)carrnm);
   if (cgrpnm) HDfreespace((VOIDP)cgrpnm);

   return; 
}

