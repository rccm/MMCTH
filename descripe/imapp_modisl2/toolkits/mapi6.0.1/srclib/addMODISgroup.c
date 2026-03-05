#include <stdio.h>
#include <stdlib.h>
#include "mapic.h" 
#include "local_nc.h"

int addMODISgroup(MODFILE *file, char const *groupname, char const *classname,
		  int32 tag, int32 ref)
/*
!C**********************************************************************
* 
*!Purpose:	To add an object to an HDF Vgroup.
* 
*!Description:	Function addMODISgroup is part of a larger software system 
* called the MODIS Applications Programming Interface (API) 
* Utility, abbreviated M-API. The M-API Utility consists of 
* subroutines which allow MODIS Science Team-supplied software 
* to read  and write data and metadata from/to HDF files. The 
* functionality of the M-API is defined in the MODIS Application 
* Program Interface (API) Specification.
* 
* addMODISgroup is an MAPI internal function which adds an HDF 
* object into a Vgroup. The group is specified by its name and class 
* name. However, the classname is an optional feature. If class name is 
* set to NULL, the function will find the Vgroup with the same name as 
* groupname without considering the classname of the group to insert 
* the object into the group. The object is specified by its tag and 
* reference number. The tag for SDS is DFTAG_NDG and for Vdata is 
* DFTAG_VH. To obtain the reference, using SDidtoref for SDS and 
* VSQueryref for Vdata. No check for the tag or reference is 
* performed.
* 
* !Input Parameters:
* 
* file		IN: 	Address of MODFILE structure that is used to 
*	                reference the MODIS-HDF file.
* groupname	IN:	ASCII string of the name of the data group.
* classname	IN:	(Optional) ASCII string of the class name of 
* 			the data group. Set to NULL for not comparing the 
* 			group classname.
* tag		IN:	The tag of the object. Valid values include:
* 			DFTAG_NDG, DFTAG_VH, and DFTAG_VG. 
* ref		IN:	The reference number of the object.
* 
* !Output Parameters:.:		none
* 
* Returns:	MAPIOK if successful, otherwise, MFAIL.
* 
* External references:
*    		PGS_SMF_MAX_MSGBUF_SIZE		PGS_SMF.h
*		VGNAMELENMAX			hdf.h
*		MFAIL				mapi.h
*		NULLMODFIL			mapic.h
*		MAPIERR				mapic.h
*		NULLstr				mapic.h
*		vgetid				hproto.h
*		Vattach				hproto.h
*		Vgetclass			hproto.h
*		Vgetname			hproto.h
*		Vdetach				hproto.h
*		Vaddtagref			proto.h
*		MAPIOK				mapi.h
*		CDF				local_nc.h
*		DIMENSION			local_nc.h
*		VARIABLE			local_nc.h
*
* !Revision History:
* $Log: addMODISgroup.c,v $
* Revision 5.1  2005/04/04 17:58:32  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.4  1996/04/02  15:43:13  qhuang
 * Added #include "local_nc.h"
 *
 * Revision 1.3  1996/01/29  17:44:58  qhuang
 * Casted hdf_id in calls Vgetid and Vattach.
 * Used HDF macros for specific group classes.
 *
 * Revision 1.2  1996/01/24  15:25:04  qhuang
 * Original development and testing for mapi version 2.0.
 *
 * Revision 1.1  1996/01/24  15:24:02  qhuang
 * Initial revision
 * 
* 
* !Team-unique Header:
* This software is developed by the MODIS Science Data Support Team for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
* 
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !References and Credits:
*
* !Design Notes:
*
!END********************************************************************
*/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning
                                           message */
  char *funcname="addMODISgroup";       /* name of this module */
  char  vgname[VGNAMELENMAX + 1];
  char  vgclass[VGNAMELENMAX + 1];
  int32 vgroup_ref;                     /* Vgroup referece number */
  int32 vg_id = 0;                          /* Vgroup identifier */
  int status_code;

  status_code = MFAIL;

  if ( NULLMODFIL(file) )
  { 
    sprintf(buff,"ERROR: addMODISgroup unable to continue with an empty\n"
			"\t file pointer.\n");
    MAPIERR(buff,funcname);
    return(status_code);
  }

  if ( NULLstr(groupname) )
  {
    
    sprintf(buff,"ERROR: addMODISgroup unable to continue with an empty\
n"
                         "\t group name.\n");
    MAPIERR(buff,funcname);
    return(status_code);	
  }

  vgroup_ref = -1;
 
  /* Loop calling Vgetid until the return value is FAIL. */
  while ( (vgroup_ref = Vgetid( (int32)file->hdf_id,vgroup_ref)) != -1 )
  {
    if ( (vg_id = Vattach( (int32)file->hdf_id,vgroup_ref,"w")) == FAIL )
    {
      sprintf(buff,"ERROR: addMODISgroup detected failures in HDF\n"
			  "\t function Vattach\n");
      MAPIERR(buff,funcname);
      return(status_code);
    }

    /* Get Vgroup classname vgclass by using Vgetclass. */
    Vgetclass(vg_id,vgclass);
    if ( (classname != NULL) || !( (strcmp(vgclass,VARIABLE) == 0)
         || (strcmp(vgclass,DIMENSION) == 0)
         || (strcmp(vgclass,CDF) == 0) ) )
    {
      /* Get vgname using Vgetname. */
      Vgetname(vg_id,vgname);
  
      if (strcmp(vgname,groupname) == 0)
        if ( (classname == NULL) || (strcmp(vgclass,classname) == 0 ) )
          break;
    }
    Vdetach(vg_id);      /* Vdetach to detach the unmatched Vgroup. */
  } /* end of while */

  if (vgroup_ref == -1)
  {
    sprintf(buff,"ERROR: addMODISgroup unable to find the specified\n"
			"\t Vgroup group %s\n",groupname);
    MAPIERR(buff,funcname);
    return(status_code);
  }

  if ( Vaddtagref(vg_id,tag,ref) == FAIL)
  {
    sprintf(buff,"ERROR: addMODISgroup detected failures in HDF function\n"
			"\t Vaddtagref\n");
    MAPIERR(buff, funcname);
  }

  else
    status_code = MAPIOK;

  /* Vdetach the matched Vgroup */
  Vdetach(vg_id);

  return(status_code);
}
