#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"

void nccrmgrp(intf *modfil, _fcd groupname, intf *lgr, _fcd classname, 
             intf *lcl, intf *ret)

/*
!C************************************************************************
* 
*!Purpose:    A wrapping function interfacing between C and FORTRAN for
*             createMODISgroup.  This C function is called only by FORTRAN
*             function GMTBL.  cgmtbl is a M-API low level internal routine.
*
*!Description:Function ccrmgrp is part of a larger software system called
*             the MODIS Application Programming Interface (API) Utility,
*             abbreviated M-API.  The M-API Utility consists of subroutines 
*             which allow MODIS Science Team-supplied software to read in 
*             Level 1B radiance bands and write out output products and 
*             metadata to HDF files.  The functionality of the M-API is 
*             defined in the MODIS API specification.
*
*             ccrmgrp is a C function callable from FORTRAN. This function
*             will call createMODISgroup to create a group structure.  In
*             M-API, ccrmgrp is low-level routine called only by CRMGRP.
*
*             In order to be callable from FORTRAN in different platforms
*             using function name ccrmgrp, this function is named nccrmgrp
*             in the actual C code.  nccrmgrp is redefined in mapic.h 
*             according to FORTRAN compilers' naming conventions/conversion
*             of each platform, so that nccrmgrp compiled by C compilers
*             will always have the object name the same as that of
*             a FORTRAN function named ccrmgrp.
*
* !Input Parameters:
* modfil      IN:  FORTRAN array of MODIS file structure.
* groupname   IN:  FORTRAN ASCII string of group name.
* lgr         IN:  Address for the numbewr of bytes in groupname.
* classname   IN:  FORTRAN ASCII string of class name.
* lcl         IN:  Address for the number of bytes in classname.
*
* !Output Parameters:
* ret         OUT:  FORTRAN integer address of the return status(MFAIL,MAPIOK)
*
* Returns:    None.
*
* External reference:
*             P_HDFID                (mapic.h)
*             P_SDID                 (mapic.h)
*             P_ACCESS               (mapic.h)
*             HDf2cstring            (hproto.h)
*             HDfreespace            (hproto.h)
*             createMODISgroup       (mapi.h)
*           
* !Revision History:
*		Qi Huang	1996/11/18
*		Version 2.2
*		Ring super structure and other changes make this
*		version much faster.
* $Log: ccrmgrp.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
*
* Xia W. Chang, Dec 26, 1995 , xia@ltpmail.gsfc.nasa.gov
* Version 2.0 original development from Prototype by Liping Di.
* 
* !Team-unique Header:
*              This software is developed by the MODIS Science Data Sup-
*              port Team for the National Aeronautics and Space Admini-
*              stration, Goddard Space Flight Center, under contract 
*              NAS5-32373.
* Refereces and Credits:
*              Portions developed at the National Center for Supercompu-
*              ting Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END**************************************************************************
*/
{
  /* declare local variables   */
  MODFILE *file;
  char *cgroupname, *cclassname;

  /* Set file by memcpy. */
  memcpy(&file,&modfil[3],sizeof(MODFILE *));

  /* Convert FORTRAN strings groupname and classname to C string. */
  cgroupname = HDf2cstring(groupname, (intn)*lgr);
  cclassname = HDf2cstring(classname, (intn)*lcl);

  if (strlen(cclassname) == 0) 
    {
      HDfreespace(cclassname);
      cclassname = NULL;
    }

  *ret = createMODISgroup(file, cgroupname, cclassname);
  
  /* free the cgroupname, and cclassname  */
  if (cgroupname) HDfreespace(cgroupname);
  if (cclassname) HDfreespace(cclassname);

  return;
}
