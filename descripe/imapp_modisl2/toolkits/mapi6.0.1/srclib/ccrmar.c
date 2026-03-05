#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mapic.h"

void nccrmar(intf *modfil, _fcd arrnm, intf *lar, _fcd grpnm, 
     intf *lgr, _fcd dtype, intf *ldt, intf *rank, intf *dims, intf *ret)

/***********************************************************************
!C
*
*!Description: 
*
*    ccrmar is a wrapping function interfacing between C and FORTRAN for 
*    createMODISarray.  This C function is only called by FORTRAN 
*    function CRMAR.  It is a M-API internal routine.
*
*    Function ccrmar is part of a larger software system called the
*    MODIS Application Program Interface (API) Utility, abbreviated M-API.
*    The M-API Utility consists of subroutines which allow MODIS Science
*    Team supplied software to read and write data and metadata from/to
*    HDF files.  The functionality of the M-API is defined in the MODIS
*    API Specification.
*
*    ccrmar is a C function callable from FORTRAN.  This function will 
*    call createMODISarray to form an array structure.  In M-API, 
*    ccrmar is a low-level routine which is called only by CRMAR. 
*
*    In order to be callable from the FORTRAN in different platforms 
*    using function name ccrmar, this function is called nccrmar in the
*    actual C code.  nccrmar is redefined in mapic.h according to 
*    compiler's FORTRAN naming conventions of each platform, so that the
*    object name of nccrmar will always the object name of a FORTRAN 
*    function named ccrmar.
*
* !Input Parameters:
*    modfil IN:  Array of FORTRAN integer array that is used to reference 
*           the MODIS-HDF file receiving the new array.
*    arrnm  IN:  FORTRAN character string that will be the name of the array.
*    lar    IN:  FORTRAN integer address of the memory size of arrnm.
*    grpnm  IN:  ASCII string name of the data group to place the new array in.
*    lgr    IN:  FORTRAN integer address of the memory size of group.
*    dtype  IN:  FORTRAN character string of data types for the array.
*    ldt    IN:  FORTRAN integer address of the memory size of dtype.
*    rank   IN:  FORTRAN integer address of the number of dimensions in the 
*           array
*    dims   IN:  The array containing the size of each dimension.
*
* !Output Parameters:
*    ret    OUT:  FORTRAN integer address of the status(MFAIL,MAPIOK) 
*
* Returns:	none
*
* Externals:
*	     MODFILE			(mapi.h)
*	     MTYPEf2c			(mapic.h)
*	     MFAIL			(mapi.h)
*            createMODISarray           (mapi.h)
*            HDf2cstring                (hproto.h)
*            HDfreespace                (hproto.h)
*            MAX_VAR_DIMS               (netcdf.h)
*            DATATYPELENMAX             (mapic.h)
* !Revision History:
*		Qi Huang	1996/07/26
*		Version 2.1
*		Ring super structure and other changes make this
*		version much faster.
*
* $Log: ccrmar.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.3  1996/08/08  16:09:34  qhuang
 * Removed redundant statement in prolog.
 *
*
* !Team-unique Header:
*
*    This software is developed by the MODIS Science Data Support Team 
*    for the National Aeronautics and Space Administration, 
*    Goddard Space Flight Center, under contract NAS5-32373.
*
*    Portions developed at the National Center for Supercomputing
*    Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !References and Credits:
*
*    Written by   Vicky Lin        01/31/96
*    Research and Data Systems Corporation
*    SAIC/GSC MODIS SCIENCE DATA SUPPORT OFFICE
*    7501 FORBES BLVD, SEABROOK MD 20706
*
*    vlin@ltpmail.gsfc.nasa.gov
*
*
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes
*
* END *****************************************************************
*/

 { 
/*
 *  declare local variables 
 */
    MODFILE *mfile;
    int  i;
    long int c_length, cdims[MAX_VAR_DIMS];
    char cdtype[DATATYPELENMAX];
    char *carrnm, *cgrpnm, *fdtype;

/* 
 *  convert FORTRAN input strings to C strings 
 */
    carrnm = HDf2cstring(arrnm,(intn)*lar);
    cgrpnm = HDf2cstring(grpnm,(intn)*lgr);
    fdtype = HDf2cstring(dtype,(intn)*ldt);

/* 
 *  Set mfile by memcpy
 */
    memcpy(&mfile,&modfil[3],sizeof(MODFILE *));

    /* Copy dims to cdims(Note: total rank elements. Copy inverselyl. */
    for (i=0; i<*rank; i++) 
      cdims[i] = (long int) dims[*rank-i-1];

    c_length = DATATYPELENMAX;
    if (MTYPEf2c(fdtype, cdtype, &c_length) == MFAIL) 
    {
      *ret = MFAIL ;
      return;                                       
    }

    *ret = createMODISarray(mfile,carrnm,cgrpnm,cdtype,(int32)*rank,cdims);

    if (carrnm) HDfreespace((VOIDP)carrnm);
    if (cgrpnm) HDfreespace((VOIDP)cgrpnm);
    if (fdtype) HDfreespace((VOIDP)fdtype);

    return;
 }
