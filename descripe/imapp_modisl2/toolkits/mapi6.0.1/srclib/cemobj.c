#include <stdlib.h>
#include <stdio.h>
#include "mapic.h"

void  ncemobj( intf *modfil, _fcd name, intf *lna, _fcd group, intf *lgr, intf *type, intf *ret)

/*
!C
*******************************************************************************************************
*
*
*!Purpose:   A wrapping function interfacing between C and FORTRAN for endMODISobjaccess.
*            This C function is only called by FORTRAN function EMOBJ. This function is a M-API internal
*            routine.
*
*!Description: Function cemobj is part of a larger software system called the MODIS Applications
*              Programming Interface (API) Utility, abbreviated M-API. The M-API Utility consists of
*              subroutines which allow MODIS Science Team-supplied software to read and write data and
*              metadata from/to HDF files. The functionality of the M-API is defined in the MODIS
*              Application Program Interface (API) Specification.
*
*            cemobj is a C function which is callable from FORTRAN. This function will 
*            call endMODISobjaccess to end the access to objects in a MODIS HDF file. 
*            In M-API, cemobj is low-level routine which is called only by EMOBJ. 
*
*            In order to be callable from the FORTRAN in different platforms using function name
*            cemobj, this function is called ncemobj in the actual C code. ncemobj is redefined in
*            mapic.h according to compiler's FORTRAN naming conventions/conversion of each platform,
*            so that the object name of ncemobj will always be the object name of a FORTRAN function
*            named cemobj.
*
*!Input Parameters:
*            modfil   IN:  FORTRAN integer array that is used to reference the MODIS-HDF file.
*
*            name     IN:  FORTRAN character string that is the name of the object.
*
*            lna      IN:  FORTRAN integer address of the memory size of name.
*
*            group    IN:  FORTRAN character string which is the name of the data group to which
*                           the object (eithe an SDS or a Vdata) belongs.
*
*            lgr      IN:  FORTRAN integer address of the memory size of group.
*
*            type     IN:  FORTRAN integer address specifying the type of objects to be closed.
*
*!Output Parameters:
*
*             ret    OUT: FORTRAN integer address of the status( number of objects closed or MFAIL
*                         if errors) 
*
* Returns:     none
*
*
* Externally defined:
*                       MODFILE                         (mapi.h)
*                       NULL                            (stdio.h)
*                       MFAIL                           (mapi.h)
*                       MAPIOK                          (mapi.h)
*                       endMODISobjaccess.c             (mapi.h)
*                       HDfreespace                     (hproto.h)                    
*                       HDf2cstring                     (hproto.h)
*
*!Revision History:
*             1996/08/06
*             Frederick J. Shaw
*             Original development for M-API version 2.1
*
*!Team-unique Header:
*      This software is developed by the MODIS Science Data Support Team for the National
*      Aeronautics and Space Administration, Goddard Space Flight Center, under contract NAS5-32373.
*
*!References and Credits:
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
*
*!Design Notes:
*
*!END
******************************************************************************************
*/

{

  MODFILE *mfile;
  char *cname, *cgroup;
 /*
  * Convert the FORTRAN character strings name and  group to C character strings cname
  * and cgroup by using HDf2cstring.     
  */
  cname = HDf2cstring(name,(intn)*lna);
  cgroup = HDf2cstring(group,(intn)*lgr);

 /*     Set mfile by memcpy                */
  memcpy(&mfile, &modfil[P_ADDR], sizeof(MODFILE *));

  *ret = endMODISobjaccess(mfile, cname, cgroup, *type);

  if (cname) HDfreespace((VOIDP)cname);
  if (cgroup) HDfreespace((VOIDP)cgroup);

  return;
}

