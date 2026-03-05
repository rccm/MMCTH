#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "mapic.h"
#include "local_nc.h"

int createMODISgroup(MODFILE *file, char const *groupname, char const *classname)

/*
!C*****************************************************************************
* 
*!Purpose:    creates a group structure in a MODIS-HDF file to sotre
*             data objects in.
*
*!Description:Function createMODISgroup is part of a larger software system 
*             called the MODIS Applications Programming Interface (API) 
*             utility, abbreviated M-API.  The M-API Utility consists of 
*             subroutines which allow MODIS Science Team-supplied software to 
*             read and write data and metadata from/to HDF files. The 
*             functionality of the M-API is defined in the MODIS API 
*             specification.
*
*             createMODISgroup(CRMGRP) creates an HDF Vgroup structure in a
*             MODIS HDF file to store table and array structures.  It must
*             be called before any of the data objects to be aggregated in it
*             are created.  The use of data groups is optional, but data
*             objects stored in them are documented in the MODIS Product File
*             Definitions in Appendix C.  A data group with the name groupname
*             must be unique in a file. This prevents confusion that is caused
*             by multiple data groups with the same name.
*
* !Input Parameters:
* file        Address of MODFILE structure that is used to reference the 
*             MODIS-HDF file receiving the new group.
* groupname   ASCII string that will be the name of the data group.
* classname   (Optional) additional ASCII string that will be the class name
*             of the group.  Set to NULL for default.
*
* !Output Parameters:None.
*
* Return      MAPIOK if successful, otherwise MFAIL.
*
* External:   PGS_SMF_MAX_MSGBUF_SIZE    (mapic.h)    
*             VGNAMELENMAX               (hdf.h)
*             NULLstr                    (mapic.h)
*             NULLMODFIL                 (mapic.h)
*             MAPIERR                    (mapic.h)
*             MFAIL                      (mapi.h)
*             MAPIOK                     (mapi.h)
*             DFACC_READ                 (hdf.h)
*             Vgetid                     (vproto.h) 
*             Vattach                    (vproto.h)
*             Vgetclass                  (vproto.h)
*             vgetname                   (vproto.h)
*             vdetach                    (vproto.h)
*             Vsetname                   (vproto.h)
*             VARIABLE                   (local_nc.h)
*             DIMENSION                  (local_nc.h)
*             CDF                        (locla_nc.h)
*
* !Revision History:
* $Log: createMODISgroup.c,v $
* Revision 5.1  2005/04/04 18:17:57  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.3  1996/05/08  17:29:11  qhuang
 * Changed Vgroup class check logic from ORs to ANDs.
 *
 * Revision 1.2  1996/04/02  15:17:44  qhuang
 * Added #include "local_nc.h"
 *
 * Revision 1.2  1996/04/02  15:17:44  qhuang
 * Added #include "local_nc.h"
 *
 * Revision 1.1  1996/04/02  15:17:17  qhuang
 * Initial revision
 *
*
*             Xia W. Chang, xia@ltpmail.gsfc.nasa.gov  , 12/12/95
*             Version 2.0 original development.
*
* !Team-unique Header:
*             This software is developed by the MODIS Science Data Support
*             Team for the National Aeronautics and Space Administration,
*             Goddard Space Flight Center, under contract NAS5-32373.
*             Portions developed at the National Center for Supercomputing
*             Applications at the Univ. of Illinois at Urbana-Champaign.
!END***********************************************************************
*/

{
  char buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char *funcname = "createMODISgroup";
  char vgroup_name[VGNAMELENMAX+1];
  char vgroup_class[VSNAMELENMAX+1]; 
  int32 vgroup_ref = -1;
  int32 vgroup_id = -1;

  if (NULLstr(groupname))
    {
      sprintf(buff, "ERROR: createMODISgroup unable to continue with"
	      " an empty group name.\n");
      MAPIERR(buff, funcname);
      return(MFAIL);
    }

  if (NULLMODFIL(file))
    {
      sprintf(buff, "ERROR: createMODISgroup unable to continue to make a"
	      " new data group %.*s with an empty file pointer.\n",
	      VGNAMELENMAX, groupname);
      MAPIERR(buff, funcname);
      return(MFAIL);
    }

  /* check whether the file is opened for read only  */
  if (file->access == DFACC_READ)
    {
      sprintf(buff, "ERROR: createMODISgroup unable to make new data group"
	      " %.*s in file opened for read only. \n", 
	      VGNAMELENMAX, groupname);
      MAPIERR(buff, funcname);
      return(MFAIL);
    }

  /* Loop calling Vgetid until the return value is -1  */
  while((vgroup_ref=Vgetid((int32)file->hdf_id, vgroup_ref)) != FAIL)
    {
      if ((vgroup_id = Vattach((int32)file->hdf_id, vgroup_ref, "r")) == FAIL )
	{
	  sprintf(buff, "ERROR: createMODISgroup detected FAIL from HDF"
		  " procedure Vattach.\n");
	  MAPIERR(buff, funcname);
	  return(MFAIL);
	}
      /* check whether the vgroup is already exists   */
      Vgetclass(vgroup_id, vgroup_class);
      if ((strcmp(vgroup_class, VARIABLE) != 0) && 
	  (strcmp(vgroup_class, DIMENSION) != 0) && 
	  (strcmp(vgroup_class, CDF) != 0) )
	{
	  Vgetname(vgroup_id, vgroup_name);  /* get vgroup_name  */
	  if (strcmp(vgroup_name, groupname) == 0)
	    {
	      Vdetach(vgroup_id);
	      sprintf(buff, "ERROR: createMODISgroup found group"
		  " %.*s already exists.\n", VGNAMELENMAX, groupname);
	      MAPIERR(buff, funcname);
	      return(MFAIL);
	    }
	}
      Vdetach(vgroup_id);
    }     /* end of while .....   */

  /* create a new vgroup using Vattach   */
  if ((vgroup_id = Vattach((int32)file->hdf_id, -1, "w")) == FAIL )
    {
      sprintf(buff, "ERROR: createMODISgroup detected FAIL from HDF"
	      " procedure Vattach while create group %.*s .\n",
	      VGNAMELENMAX, groupname);
      MAPIERR(buff, funcname);
      return(MFAIL);
    }

  Vsetname(vgroup_id, groupname);          /* set the Vgroup name  */

  if (classname != NULL) Vsetclass(vgroup_id, classname); /* set class name */
  
  Vdetach(vgroup_id);
  return(MAPIOK);
}
	  

