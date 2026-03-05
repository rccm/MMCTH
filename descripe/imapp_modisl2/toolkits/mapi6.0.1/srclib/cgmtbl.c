#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "mapic.h"

void ncgmtbl(intf *modfil, _fcd tbname, intf *ltb, _fcd group, intf *lgr,
             _fcd fldnm, intf *lfl, intf *start, intf *recno, intf *bsize, 
             VOIDP data, intf *ret)

/*
!C************************************************************************
* 
*!Purpose:    A wrapping function interfacing between C and FORTRAN for
*             getMODIStable.  This C function is called only by FORTRAN
*             function GMTBL.  cgmtbl is a M-API low level internal routine.
*
*!Description:Function cgmtbl is part of a larger software system called
*             the MODIS Application Programming Interface (API) Utility,
*             abbreviated M-API.  The M-API Utility consists of subroutines 
*             which allow MODIS Science Team-supplied software to read in 
*             Level 1B radiance bands and write out output products and 
*             metadata to HDF files.  The functionality of the M-API is 
*             defined in the MODIS API specification.
*
*             cgmtbl is a C function callable from FORTRAN. This function
*             will call getMODIStable to get data from a Vdata.
*
*             In order to be callable from FORTRAN in different platforms
*             using function name cgmtbl, this function is named ncgmtbl
*             in the actual C code.  ncgmtbl is redefined in mapic.h 
*             according to FORTRAN compilers' naming conventions/conversion
*             of each platform, so that ncgmtbl compiled by C compilers
*             will always have the object name the same as that of
*             a FORTRAN function named cgmtbl.
*
* !Input Parameters:
* mfile		IN: 	FORTRAN integer array that is used to 
* 		reference the MODIS-HDF file.
* tbname	IN:	FORTRAN character string of the name of the 
* 		table.
* ltb		IN:	FORTRAN integer address of the memory size 
* 		of tbname.
* group		IN:	ASCII string name of the data group the table 
* 		will be retrieved.
* lgr		IN:	FORTRAN integer address of the memory size 
* 		of group.
* fldnm		IN	FORTRAN character string of comma delimited 
* 		name of fields to be retrieved.
* lfl		IN:	FORTRAN integer address of the length of 
* 		fldnm.
* start		IN:	FORTRAN integer address of the start record, 
* 		0-based.
* recno		IN:	FORTRAN integer address of the number of 
* 		records to be retieved.
* bsize		IN/OUT: FORTRAN integer address containing the 
*		memory size of the data buffer. Output will be 
*		replaced with the actual bytes required to 
*		hold the data.
* 
* !Output Parameters:
* data		OUT:	The data buffer to store the retrieved table 
* 		data. 
* ret		OUT:	FORTRAN integer address of the return 
* 		status(MFAIL,MAPIOK) 
* 
* Returns:	none
*
* External reference:
*	      VOIDP		  (hdfi.h)
*	      MODFILE		  (mapi.h)
*             HDf2cstring         (hproto.h)
*             getMODIStable       (mapi.h)
*             HDfreespace         (hproto.h)
*
* !Revision History:
*		Qi Huang	1996/09/30
*		Version 2.2
*		Ring super structure and other changes make this version 
*		much faster.
*
* $Log: cgmtbl.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/10/28  19:01:40  qhuang
 * Version 2.2
 *
 * Revision 1.1  1996/10/28  19:00:53  qhuang
 * Initial revision
 *
*
* Revision 01.00 12/01/95
* Xia W. Chang, xia@ltpmail.gsfc.nasa.gov
* Version 1.4 original development from Prototype by Liping Di.
* 
* !Team-unique Header:
*              This software is developed by the MODIS Science Data Sup-
*              port Team for the National Aeronautics and Space Admini-
*              stration, Goddard Space Flight Center, under contract 
*              NAS5-32373.
* !References and Credits:
*              Portions developed at the National Center for Supercompu-
*              ting Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !Disign Notes:
*
!END**************************************************************************
*/
    
{
 /* declare local variables   */
  MODFILE *file;
  char *tablename, *groupname, *fieldname;
  long int buffersize;

  /* Convert FORTRAN strings group, tbname, fldnm to C string groupname,
     tablename, and fieldname using HDf2cstring                         */
  groupname = HDf2cstring(group, (intn)*lgr);
  tablename = HDf2cstring(tbname, (intn)*ltb);
  fieldname = HDf2cstring(fldnm, (intn)*lfl);
  
  /* Set file by memcpy */
  memcpy(&file,&modfil[3],sizeof(MODFILE *));

  /* set buffersize to *bsize       */
  buffersize = *bsize;

  /* call getMODIStable to get the table  */
  *ret = getMODIStable(file, tablename, groupname, fieldname, 
                       (long int)*start, (long int)*recno, &buffersize,
                       (unsigned char *)data );
  
  /* free tablename, grupname, and fieldname if they are not equal to 
     null by calling HDfreesapce */                                     
  if (tablename)  HDfreespace((VOIDP)tablename);
  if (groupname)  HDfreespace((VOIDP)groupname);
  if (fieldname)  HDfreespace((VOIDP)fieldname);

  *bsize = (intf)buffersize;

  return;

}




