#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"

void ncpmtbl(intf *modfil, _fcd tbname, intf *ltb, _fcd group, intf *lgr,
	     intf *start, intf *recno, VOIDP data, intf *ret)

/*
!C************************************************************************
*
*!Purpose:    A wrapping function interfacing between C and FORTRAN for
*             putMODIStable.  This C function is called only by FORTRAN
*             function PMTBL.  cpmtbl is a M-API low level internal routine.
*
*!Description:Function cpmtbl is part of a larger software system called
*             the MODIS Application Programming Interface (API) Utility,
*             abbreviated M-API.  The M-API Utility consists of subroutines 
*             which allow MODIS Science Team-supplied software to read in 
*             Level 1B radiance bands and write out output products and 
*             metadata to HDF files.  The functionality of the M-API is 
*             defined in the MODIS API specification.
*
*             cpmtbl is a C function callable from FORTRAN. This function
*             will call putMODIStable to write data records in a specified
*             HDF Vdata.
*
*             In order to be callable from FORTRAN in different platforms
*             using function name cpmtbl, this function is named ncpmtbl
*             in the actual C code.  ncpmtbl is redefined in mapic.h 
*             according to FORTRAN compilers' naming conventions/conversion
*             of each platform, so that the object name of C ncpmtbl will
*             always be the same as the object name of a FORTRAN function 
*             named cpmtbl.
*
* !Input Parameters:
* modfil      IN:  FORTRAN integer array that is used to reference the MODIS-HDF
*             file.
* tbname      IN:  FORTRAN character string of the name of the table.
* ltb         IN:  FORTRAN integer address of the memory size of tbname.
* group       IN:  ASCII string name of the data group the table will be put.
* lgr         IN:  FORTRAN integer address of the memory size of group.
* start       IN:  FORTRAN integer address of the start record, 0-based.
* recno       IN:  FORTRAN integer address of the number of records to be put.
* data        IN:  Data buffer of the table.
*
*!Output Parameters:
* ret         OUT:  FORTRAN integer address of the return status(MFAIL, MAPIOK)
*
* Returns:    None.
*
* External reference:
*	      MODFILE		  (mapi.h)
*             HDf2cstring         (hproto.h)
*             HDfreespace         (hproto.h)
*             putMODIStable       (mapi.h)
*
* !Revision History:
*		Qi Huang	1996/10/02
*		Version 2.2
*		Ring super structure and other changes make this version
*		much faster.
* $Log: cpmtbl.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
*
*              Xia W. Chang, Dec 5 1995, xia@ltpmail.gsfc.nasa.gov
*              Version 1.4 original development from prototype by Liping Di.
*
* !Team-unique Header:
*              This software is developed by the MODIS Science Data Sup-
*              port Team for the National Aeronautics and Space Admini-
*              stration, Goddard Space Flight Center, under contract 
*              NAS5-32373.
*
* !References and Credits:
*              Portions developed at the National Center for Supercompu-
*              ting Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END**************************************************************************
*/
{
  MODFILE *file;
  char *tablename, *groupname;
  
  /* Convert group and tbname to C string  */
  groupname = HDf2cstring(group, (intn)*lgr);
  tablename = HDf2cstring(tbname,(intn)*ltb);

  /* Set file by memcpy */
  memcpy(&file,&modfil[3],sizeof(MODFILE *));
  *ret = putMODIStable(file, tablename, groupname, (long int)*start,
		      (long int)*recno, (unsigned char *)data); 
 
  /* free tablename and groupname if they are not NULL   */
  if (tablename) HDfreespace((VOIDP)tablename);
  if (groupname) HDfreespace((VOIDP)groupname);

  return;
}
