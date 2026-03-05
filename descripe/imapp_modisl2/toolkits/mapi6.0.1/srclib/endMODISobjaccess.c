#include <stdlib.h>   
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include "mapic.h"

int endMODISobjaccess(MODFILE *file, char const *name, char const *group, 
                      long int type)
/*
!C**********************************************************************
* 
* !Purpose:	ends the access to an individual or a group of opened objects 
* 		(close the object and release resources).
* 
* !Description:  Function endMODISobjaccess is part of a larger software system 
*		 called the MODIS Applications Programming Interface (API) 
*		 Utility, abbreviated M-API.  The M-API Utility consists of 
*		 subroutines which allow MODIS Science Team-supplied software 
*		 to read  and write data and metadata from/to HDF files. The 
*		 functionality of M-API is defined in the MODIS Application 
*		 Program Interface (API) Specification.
*		 
*		 A Vdata or SDS is opened and remains opened when application 
*		 programs call any of M-API routines to access the object. The M-
*		 API access routines automatically call M-API internal routine 
*		 addid to insert the HDF object ID and related information into the 
*		 ring super structure in MODFILE. In the subsequent calls to 
*		 access the same object, M-API routines will use the id stored in 
*		 the ring super structure for fast data access.   
*		 
*		 endMODISobjaccess ends the access to an individual or a group of 
*		 opened objects by deleting objects' DATAID structure from the 
*		 ring super structure, releasing memory, and detaching (for 
*		 Vdata) or end-accessing (for SDS) the objects. This routine is 
*		 called by closeMODISfile or completeMODISfile, which calls 
*		 closeMODISfile, so that all opened objects will be closed 
*		 automatically before the MODIS HDF file is closed. As long as an 
*		 application program calls closeMODISfile or completeMODISfile, 
*		 the application does not need to worry  when or how to use this 
*		 routine to close an object or a group of object. However, if an 
*		 application program determines an object will no longer be 
*		 accessed and wish to end  the access to the object for releasing 
*		 computer resource, the application program can call this 
*		 routines any time before the application program closes the 
*		 opened HDF file.
* 
* !Input parameters:
* file	IN: 	Address of MODFILE structure that is used to 
* 		reference a MODIS HDF file containing objects to be 
* 		closed (Vdata or SDS).
* name	IN:	The name of the object. If the name is set to 
* 		NULL, all objects matched with group and object 
* 		type will be closed.
* group	IN:	The name of the group to which objects 
* 		belongs. If group is set to NULL, all lone objects 
*		(objects belonging to no group) matched with  name 
* 		and type will be closed. If both name and group are 
* 		NULL, all objects matched with type will be closed.
* type	IN:	The object type: MODIS_ARRAY (for SDS, the 
*		numerical value is 720, the same value as 
*		DFTAG_NDG),  MODIS_TABLE(for Vdata, the 
*		numerical value is 1962, the same value as 
*		DFTAG_VH), or MODIS_ALL_TYPES (the numerical 
*		values is 0). If type is 720, only SDS objects will be 
*		closed.  If type is 1962, only Vdata will be closed. If 
*		type is MODIS_ALL_TYPES, both Vdata and SDS 
*		objects specified by name and group will be closed. 
*		Therefore, to close all opened objects, set both name 
*		and group to NULL and set type to 0.
* 
* !Output Parameters:		None
* 
* 
* Returns:	Number of objects closed, or MFAIL if errors.
* 
* External references:
*		MODFILE                         (mapi.h)
*               PGS_SMF_MAX_MSGBUF_SIZE         (mapic.h)
*               DATAINFO                        (mapi.h)
*               DATAID                          (mapi.h)
*		SDSINFO				(mapic.h)
*               MODIS_ARRAY                     (mapi.h)
*		MODIS_TABLE			(mapi.h)
*		MODIS_ALL_TYPES			(mapi.h)
*               MAPIERR                         (mapic.h)
*               MFAIL                           (mapi.h)
*               MAX_NC_NAME                     (netcdf.h)
*               MAPIOK                          (mapi.h)
*		SDendaccess			(mfhdf.h)
*		VSdetach			(vproto.h)
*		VSNAMELENMAX			(hdf.h)
*		
* !Revision History:
*		Qi Huang	1996/11/14
*		Version 2.2
*		Implementing the Vdata part of the ring superstructure.
*
* 		Qi Huang	1996/07/09
*		Version 2.1
*		Original development and testing
* $Log: endMODISobjaccess.c,v $
* Revision 5.1  2005/04/04 18:21:25  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/07/16  18:15:31  qhuang
 * Added SDSINFO in external reference in prolog; changed did->group != NULL to
 * did->group == NULL in an if statement; added some commnets.
 *
 * Revision 1.1  1996/07/16  18:14:40  qhuang
 * Initial revision
 *
*
* !Team-unique header:
*	This software is developed by the MODIS Science Data Support Team for 
* 	the National Aeronautics and Space Administration, Goddard Space 
* 	Flight Center, under contract NAS5-32373.
* 
*!References and Credits
*
*!Design Notes
*
!END********************************************************************
*/
{
   char  buff[PGS_SMF_MAX_MSGBUF_SIZE];         /* buffer to hold the error */
                                                /* warning message      */
   char *funcname="endMODISobjaccess";              /* name of this routine */
   DATAINFO *dinfo;
   DATAID   *did, *didroot, *didprev;
   int	    ndata;	/* number of opened objects before closing */
   int	    ncurrent;   /* number of currently opened objects */
   int	    i,j;
   int	    ret;	/* reutrn value */
   int	    start, end; /* loop start and end variables for closing 
			   SDS, Vdata or both */

   /* Input check: */
   if ( NULLMODFIL(file) )
   {
     sprintf(buff,"ERROR: endMODISobjaccess unalbe to close object %.*s\n"
				"\t with an invalid MODIS file structure input.\n",
				MAX_NC_NAME,name);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( (type != MODIS_ARRAY) && (type != MODIS_TABLE) 
	&& (type != MODIS_ALL_TYPES) )
   {
     sprintf(buff,"ERROR: endMODISobjaccess unable to close objects with\n"
			"\t an invalid MODIS object type: %ld.\n",type);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   ret = 0;			/* 0 object closed */
   dinfo = file->dinfo;
   if ( dinfo == NULL )
     return(ret);

   if ( (name != NULL) && (*name == '\0') )
     name = NULL;
   if ( (group != NULL) && (*group == '\0') )
     group = NULL;

   /* Set the loop variables */
   start = 0;
   end = 2;
   if ( type == MODIS_ARRAY )
     end = 1;
   if ( type == MODIS_TABLE )
     start = 1;

   /* Delete and close the opened objects */
   /* Loop i from start to end - 1 */
   for ( i=start; i<end; i++)
   {
     if ( i == 0 )
     {
       ncurrent = ndata = dinfo->nsds;
       didroot = did = dinfo->sds;
     }
     else
     {
       ncurrent = ndata = dinfo->nvd;
       didroot = did = dinfo->vd;
     }
 
     for ( j=0; j<ndata; j++)
     /* Loop j from 0 to ndata - 1 */
     {
       if ( ((name == NULL) && (group == NULL)) ||
	    ((name == NULL) && (group != NULL) && 
	     (did->group != NULL) && (strcmp(group,did->group) == 0)) ||
	    ((name != NULL) && (group == NULL) &&
	     (did->group == NULL) && (strcmp(name,did->name) == 0))   ||
	    ((name != NULL) && (group != NULL) &&
	     (did->name != NULL) && (did->group != NULL) &&
	     (strcmp(name,did->name) == 0) && (strcmp(group,did->group) == 0)) )
       {
	 ncurrent --;
	 ret ++;

	 /* reconnet the ring super structure */
	 if ( ncurrent > 0 )
	 {
	   did->prev->next = did->next;
	   did->next->prev = did->prev;
	   if ( didroot == did )
	     didroot = did->prev;
	 }

	 /* close the object */
	 if ( i == 0 )		/* SDS */
	 {
	   if ( SDendaccess((int32)did->id) == FAIL )
	   {
	     sprintf(buff,"ERROR: endMODISobjaccess  detected FAIL from HDF"
				"procedure SDendaccess\n"
				"\t while closing access to array %.*s.\n",
				MAX_NC_NAME,did->name);			
	     MAPIERR(buff,funcname);
	     ret = MFAIL;
	   }
	 }
	 else		/* Vdata */
	 {
	   if ( VSdetach((int32)did->id) == FAIL )
	   {
	     sprintf(buff,"ERROR: endMODISobjaccess  detected FAIL from HDF"
				"procedure VSdetach\n"
				"\t while closing access to table %.*s\n",
				VSNAMELENMAX,did->name);
	     MAPIERR(buff,funcname);
	     ret = MFAIL;
	   }
	 }

	 /* free the memory allocated to did */
	 didprev = did->prev;
	 free(did->name);
	 if ( did->group )
	   free(did->group);
	 if ( i == 0 )		/* free SDSINFO */
	 {
	   free(((SDSINFO *)did->info)->dimsizes);
	   free(did->info);
	 }
	 else
	 {
	   if ( ((VDINFO *)did->info)->read_fields )
	     free( ((VDINFO *)did->info)->read_fields );
	   free(did->info);	/* free VDINFO */
	 }
	 free(did);

	 if ( (i == 0) && (name) )	/* For SDS, break */
	   break;

	 if ( ret == MFAIL )
	   break;
	 did = didprev;
       }
       else
	 did = did->prev;
     } /* end of loop, for (j=0... */

     /* update dinfo */
     if (i == 0)	/* SDS */
     {
       dinfo->nsds = ncurrent;
       if ( ncurrent == 0 )
	 dinfo->sds = NULL;
       else
	 dinfo->sds = didroot;
     }
     else
     {
       dinfo->nvd = ncurrent;
       if ( ncurrent == 0 )
	 dinfo->vd = NULL;
       else
	 dinfo->vd = didroot;
     }

     if (ret == MFAIL)
       break;
   }/* end of loop, for (i=start... */

   if ( (dinfo->nsds == 0) && (dinfo->nvd == 0) )
   {
     free(dinfo);
     file->dinfo = NULL;
   }

   return(ret);
 } 
