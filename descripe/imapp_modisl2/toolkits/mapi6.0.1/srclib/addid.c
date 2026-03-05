#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"

int addid(MODFILE *file, char const *name, char const *group, 
          int32 id, int32 type, char const *access)
/*
!C**********************************************************************
* 
*!Purpose:	Add newly opened data object (vdata or SDS) id into the DATAID 
* 		structure in the MODFILE structure.
* 
*!Description: 	Function addid is part of a larger software system called the 
*		MODIS Applications Programming Interface (API) Utility, 
*		abbreviated M-API.  The M-API Utility consists of subroutines 
*		which allow MODIS Science Team-supplied software to read  and 
*		write data and metadata from/to HDF files. The functionality of M-
*		API is defined in the MODIS Application Program Interface (API) 
*		Specification.
*		
*		addid creates a DATAID structure which contains the id of the 
*		opened object and inserts the structure to the ring super 
*		structure. In SDS case, the DATAID structure also includes a void 
*		pointer which points to an SDSINFO structure created and filled 
*		by addid. The SDSINFO structure contains the rank, dimension size 
*		and data type information for the opened SDS.   
*		
*		addid is an internal M-API routine which is used only by M-API 
*		routines. Application programs should not call this routine 
*		directly. Minimum input argument checking is performed since  
*		those arguments are checked in the calling routine.
* 
* !Input Parameters:
*	 file	IN: 	Address of MODFILE structure that is used to 
* 		reference a MODIS HDF file containing the newly 
* 		opened data object (Vdata or SDS).
* 	 name	IN:	The name of the object.
* 	 group	IN:	The name of the group to which the object 
* 		belongs. If set to NULL the object is a lone object.
*	 id	IN:	The id for the object returned from HDF object 
* 		open routines such as VSattach or SDselect..
*  	 type	IN:	The object type: either DFTAG_NDG (for SDS) 
* 		or DFTAG_VH (for Vdata).
* 	 access	IN:	The access type of the object: "r" or "w". This 
* 		parameter is ignored when type is  DFTAG_NDG.
* 
* !Output Parameters:		NONE
* 
* Returns:	MAPIOK if successful, MFAIL if an error occurs.
* 
* External references:
*		MODFILE				(mapi.h)
*		PGS_SMF_MAX_MSGBUF_SIZE		(mapic.h)
*		DATAINFO			(mapi.h)
*		DATAID				(mapi.h)
*		SDSINFO				(mapic.h)
*		VDINFO				(mapic.h)
*		MAX_VAR_DIMS			(netcdf.h)
*		DFTAG_VH			(hdf.h)
*		MAPIERR				(mapic.h)
*		MFAIL				(mapi.h)
*		MAX_NC_NAME			(netcdf.h)
*		DFTAG_NDG			(hdf.h)
*		MAPIOK				(mapi.h)
*		SDgetinfo			(mfhdf.h)
*		VDINFO				(mapic.h)
*
* !Revision History:
*		Qi Huang	1996/09/16
*		Version 2.2
*		Implementing Vdata part.  Initialiaing the VDINFO structure.
*
* 		Qi Huang	1996/07/10
*		Version 2.1 
*		Original development and testing
*
* $Log: addid.c,v $
* Revision 5.1  2005/04/04 18:03:52  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/08/09  17:06:07  qhuang
 * Changed coupled assignments to seperated, individual assignments, cast
 * some argument in malloc call.
 *
 * Revision 1.1  1996/08/09  17:05:36  qhuang
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
{ char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning
                                           message */
  char *funcname="addid";       /* name of this module */

  DATAINFO *dinfo;
  DATAID *did, *didroot;
  SDSINFO *sinfo;
  VDINFO *vinfo;
  int ndata;
  int32 dimsizes[MAX_VAR_DIMS];
  int32 nattrs;			/* number of "netCDF-style attributes
				   for the data set */
  int i;			/* loop control variable */

  /* Input check: */
  if ( (type == DFTAG_VH) && ((access[0] != 'r') && (access[0] != 'w')) )
  {
    sprintf(buff, "ERROR: addid found the access type %s is\n"
			"\t unrecognizable.\n",access);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  dinfo = file->dinfo; /* dinfo will be released when the file is closed */

  if ( dinfo == NULL ) /* allocate the DATAINFO structure */
  {
    if ( (dinfo = (DATAINFO *)malloc(sizeof(DATAINFO))) == NULL )
    {
      sprintf(buff, "ERROR: addid unable to allocate a DATAINFO\n"
  		"\t structure while attempting to insert data\n"
  		"\t id of %.*s into MODFILE structure.\n",
  		MAX_NC_NAME,name);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }
  
    dinfo->nsds = 0;
    dinfo->nvd = 0;
    dinfo->vd = 0;
    dinfo->sds = NULL;
    file->dinfo = dinfo;
  }

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

  for (i=0; i<ndata; i++) /* to find if the object is already opened */
  {
     if ( did->id == id )
       if ( ((group == NULL) && (did->group == NULL))
	    || ( (group != NULL) && (did->group != NULL) &&
		 (strcmp(group,did->group) == 0) ) )
	 if ( strcmp(name,did->name) == 0 )
	   return(MAPIOK);
     did = did->next;
  }

  /* prepare for the new DATAID structure */
  if ( (did = (DATAID *)malloc(sizeof(DATAID))) == NULL )
  {
    sprintf(buff,"ERROR: addid unable to allocate a DATAID structure\n"
			"\t while attempting to insert data id of %.*s into\n"
			"\t MODFILE structure.\n",MAX_NC_NAME,name);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( (did->name = (char *)malloc(strlen(name) + 1)) == NULL )
  {
    free(did);
    sprintf(buff,"ERROR: addid unable to allocate space for storing\n"
			"\t object's name while attempting to insert data id\n"
			"\t of %.*s into MODFILE structure.\n",MAX_NC_NAME,name);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  strcpy(did->name,name);
  did->id = id;
  did->info = NULL;
  
  if ( group != NULL )
  {
    if ( (did->group = (char *)malloc(strlen(group) + 1)) == NULL )
    {
      free(did->name);
      free(did);
      sprintf(buff,"ERROR: addid unable to allocate space for\n"
			"\t storing object's group name while\n"
			"\t attempting to insert data id of %.*s into\n"
			"\t MODFILE structure.\n",MAX_NC_NAME,name);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }
    strcpy(did->group,group);
  }

  else
    did->group = NULL;

  /* Set the SDSINFO structure if the object is an SDS */
  if ( type == DFTAG_NDG )
  {
    if ( (sinfo = (SDSINFO *)malloc(sizeof(SDSINFO))) == NULL )
    {
      free(did->name);
      if ( did->group != NULL )
	free(did->group);
      free(did);
      sprintf(buff,"ERROR: addid unable to allocate an SDSINFO\n"
			"\t structure while attempting to insert data\n"
			"\t id of %.*s into MODFILE structure.\n",MAX_NC_NAME,name);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }

    if ( SDgetinfo(id,did->name,&(sinfo->rank),dimsizes,
		  &(sinfo->ntype),&nattrs) == FAIL )
    {
      free(sinfo);
      free(did->name);
      if (did->group != NULL )
	free(did->group);
      free(did);
      sprintf(buff,"ERROR: addid detected FAIL from HDF procedure\n"
			"\t SDgetinfo while attempting to  insert data\n"
			"\t id of %.*s into MODFILE structure.\n",
			MAX_NC_NAME,name);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }

    if ( (sinfo->dimsizes = (int32 *)malloc((size_t)sinfo->rank * sizeof(int32))) == NULL )
    {
      free(sinfo);
      free(did->name);
      if ( did->group != NULL )
	free(did->group);
      free(did);
      sprintf(buff,"ERROR: addid unable to allocate space for\n"
			"\t storing dimsizes while attempting to \n"
			"\t insert data id of %.*s into MODFILE\n"
			"\t structure.\n",MAX_NC_NAME,name);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }

    for (i=0;i<sinfo->rank;i++)
      sinfo->dimsizes[i] = dimsizes[i];
    did->info = sinfo;
  }

  else
  {
    if ( (vinfo = (VDINFO *)malloc(sizeof(VDINFO))) == NULL )
    {
      free(did->name);
      if (did->group != NULL)
	free(did->group);
      free(did);
      sprintf(buff,"ERROR: addid unable to allocate VDINFO structure\n"
			"\t while attempting to  insert data id of %.*s\n"
			"\t into MODFILE structure.\n",MAX_NC_NAME,name);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }

    else
    {
      vinfo->access = access[0];
      vinfo->read_fields = NULL;
      vinfo->size = 0;
      did->info = vinfo;
    }
  }

  /* Insert did into the ring super structure */
  if ( didroot == NULL )
  {
    did->next = did;
    did->prev = did;
  }

  else
  {
    did->next = didroot;
    did->prev = didroot->prev;
    didroot->prev->next = did;
    didroot->prev = did;
  }

  /* Update DATAINFO structure */
  if ( type == DFTAG_NDG )
  {
    dinfo->nsds ++;
    dinfo->sds = did;
  }

  else
  {
    dinfo->nvd ++;
    dinfo->vd = did;
  }

  return(MAPIOK);
}


