#include <stdlib.h>  
#include <stdio.h>
#include <string.h>
#include "mapic.h"

DATAID *searchid(MODFILE *file, char const *name, char const *group, 
                 int32 type, char const *access)
/*
!C**********************************************************************
* 
*!Purpose:	 Searches the ring super structure to find if an object is opened. 
*		If it is opened, return the address of the DATAID structure.
*		
*!Description: 	Function searchid is part of a larger software system called the 
*		MODIS Applications Programming Interface (API) Utility, 
*		abbreviated M-API.  The M-API Utility consists of subroutines 
*		which allow MODIS Science Team-supplied software to read  and 
*		write data and metadata from/to HDF files. The functionality of M-
*		API is defined in the MODIS Application Program Interface (API) 
*		Specification.
*		
*		searchid searches the ring super structure to find if an object is 
*		opened. If it is opened, return the address of the DATAID 
*		structure. The function also maintains the ring so that the first 
*		object in the ring is always the current accessed object. Please 
*		note any routines using contents in the DATAID structure should 
*		neither change the contents nor release the DATAID structure.   
*		
*		searchid is an internal M-API routine which is used only by M-
*		API routines. Application programs should not call this routine 
*		directly. Since this function will be only used by M-API functions 
*		and the speed of this function is important for overall 
*		performance of M-API, minimum error checking is performed in 
*		this function because all possible errors are checked in the 
*		calling routines.
* 
* !Input Parameters:
* 	file	IN: 	Address of MODFILE structure that is used to 
* 		reference a MODIS HDF file containing opened data 
* 		objects (Vdata or SDS).
* 	name	IN:	The name of the object.
* 	group	IN:	The name of the group to which the object 
*		belongs. If set to NULL the object is a lone object.
*		type	IN:	The object type: either DFTAG_NDG (for SDS) 
*		or DFTAG_VH (for Vdata).
* 	access	IN:	The access type, either "w" or "r". This 
* 		argument is ignored when type is DFTAG_NDG.
* 
*!Output Parameters:		None
* 
* 
* Returns:	NULL if the object is not found in the ring super structure. 
* 		The address of the DATAID structure if the object is found.
* External reference:
*		MODFILE				(mapi.h)
*		DATAINFO			(mapi.h)
*		DATAID				(mapi.h)
*		DFTAG_NDG			(hdf.h)
*		VDINFO				(mapic.h)
*
* !Revision History:
* 		Qi Huang	1996/07/22
*		Version 2.1
*		Original development and testing
*		Ring super structure and other changes make this
*		this version much faster.
*
*		Qi Huang	1996/09/18
*		Version 2.2
*		Implementation of Vdata part of ring super structure.
* 
* $Log: searchid.c,v $
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.3  1996/08/09  16:20:43  qhuang
 * Small changes in prolog.
 *
 * Revision 1.2  1996/07/30  21:27:23  qhuang
 * Removed redundant variable declarations, added NCSA acknowledgement and
 * missing ! in prolog.
 *
 * Revision 1.1  1996/07/30  21:26:42  qhuang
 * Initial revision
 *
*
* !Team-unique Header:
* This software is developed by the MODIS Science Data Support Team for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
* 
*!References and Credits
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes
*
!END********************************************************************
*/
{ 
  DATAINFO	*dinfo;
  DATAID	*did, *didroot;
  int		ndata, i;

  dinfo = file->dinfo;
  if ( dinfo == NULL )
    return(NULL);

  if ( type == DFTAG_NDG )
  {
    did = dinfo->sds;
    ndata = dinfo->nsds;
  }
  else
  {
    did = dinfo->vd;
    ndata = dinfo->nvd;
  }
  didroot = did;

  for ( i=0; i<ndata; i++)	/* to find it the object is already opened */
  {
    if ( ((group == NULL) && (did->group == NULL)) ||
	 ((group != NULL) && (did->group != NULL) && 
	  (strcmp(group,did->group) == 0)) )
      if ( strcmp(name,did->name) == 0 )
	if ( (type == DFTAG_NDG) || (access[0] == 'r')
	     || (((VDINFO *)did->info)->access == access[0]) )
	  break;
    did = did->prev;
  }

  if ( i == 0 )		/* current object is the neede one */
    return(did);

  if ( i == ndata )	/* no object is found */
    return(NULL);

  /* There are more than 1 object between the last accessed object and
     current object in the ring super structure.  Therefore, reconnect the
     ring so that the last accessed object and current object become neighbor */
  if ( i > 1 )
  {
    did->prev->next = did->next;
    did->next->prev = did->prev;
    didroot->prev->next = did;
    did->prev = didroot->prev;
    did->next = didroot;
    didroot->prev = did;
  }

  /* point to the current object */
  if ( type == DFTAG_NDG )
    dinfo->sds = did;
  else
    dinfo->vd = did;

  return(did);
}
