#include <stdlib.h>
#include <stdio.h>
#include "mapic.h"

int getMODIShobjid(MODFILE *file, char const *name, char const *group, 
                   long int type, char const *access)
/*
!C****************************************************************************
*
*!Purpose:	 get the id of an existed HDF object (either array or table) in a 
* MODIS HDF file so that a M-API application can directly call the HDF 
* library functions to access the data in MODIS HDF file in case that the 
* application program needs to do so.
* 
*!Description: 	Function getMODIShobjid is part of a larger software 
*	system called the MODIS Applications Programming Interface (API) 
*	Utility, abbreviated M-API.  The M-API Utility consists of 
*	subroutines which allow MODIS Science Team-supplied software 
*	to read  and write data and metadata from/to HDF files. The 
*	functionality of M-API is defined in the MODIS Application 
*	Program Interface (API) Specification.
*	
*	getMODIShobjid returns the HDF object id. The application 
*	program can directly use this id in the native HDF library 
*	function calls because the HDF object is already attached. If the 
*	type is specified as MODIS_ARRAY, the return value is an array id 
*	(sds id). If the type is MODIS_TABLE, the return value is the table 
*	id (vdata id). The intended use of this function is for M-API 
*	applications which want  to call HDF functions directly. The call 
*	sequence in an M-API applications should be:
*	1). Use openMODISfile to open a MODIS HDF file.
*	2). Call getMODIShobjid to obtain an object's HDF id.
*	3). Call HDF library functions by using the obtained id.
*	4). Call closeMODISfile or completeMODISfile to end the 
*	access.
*	The calls to HDF library functions can be mixed with calls to M-
*	API functions.
*	NOTE: Please do not use the HDF VSdetach (for table) or 
*	SDendaccess (for array) to end the access to the object. The 
*	closeMODISfile/completeMODISfile will take care of all opened 
*	objects. In case  of application programs need to end the access
*	to an object while keeping the MODIS file open,  use the M-API 
*	function endMODISobjaccess.
*
* !Input Parameters:
* file		IN: 	Address of MODFILE structure that is used to 
* 		reference a MODIS HDF file containing objects to be 
* 		accessed (Vdata or SDS).
* name		IN:	The name of the object.
* group		IN:	The name of the group to which objects 
* 		belongs. If set to NULL the entire file will be 
* 		searched for the object named name.
* type		IN:	The type of object: MODIS_ARRAY (for SDS, 
*		the numerical value is 720, the same value as 
*		DFTAG_NDG),  MODIS_TABLE(for Vdata, the 
*		numerical value is 1962, the same value as 
*		DFTAG_VH).
* access	IN:	"r"	The object is attached for ready only.
*		"w"	The object is attached for read/write.
*		If the MODIS file is opened for read only, and the 
*		application program calls this function to obtain an 
*		object id with write operation, the function will 
*		return an error.
* 
*!Output Parameters:		None
* 
* 
* Returns:	The object id or MFAIL.
* 
*
*Externally defined:
*			MODFILE                (mapi.h)
*                       PGS_SMF_MSGBUF_SIZE    (mapic.h)
*			DATAID		       (mapi.h)
*                       MAPIERR                (mapic.h)
*                       MFAIL                  (mapi.h)
*			MODIS_ARRAY	       (mapi.h)
*			MODIS_TABLE	       (mapi.h)
*			NULLMODFIL	       (mapic.h)
*			DFACC_READ	       (hdf.h)
*			getMODISarrayid	       (mapic.h)
*			getMODIStableid	       (mapic.h)
*			MAX_NC_NAME	       (netcdf.h)
*
* !Revision History:
*		Qi Huang	1996/11/20
*		Version 2.2
*		Ring super structure and other changes make this 
*		version much faster.
*
* $Log: getMODIShobjid.c,v $
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
*
* !Team-unique Header:
*	      This software is developed by the MODIS Science Data Support
*	      Team for the National Aeronautics and Space Administration,
*	      Goddard Space Flight Center, under contract NAS5-32373.
*
* !References and Credits:
*	      Portions developed at the National Center for Supercomputing
*	      Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END*********************************************************************
*/
{ 
   char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
   char *funcname="closeMODISfile"; 
   DATAID *did;

   /* Input check */
   if ( (name == NULL) || (*name == '\0'))
   {
     sprintf(buff,"ERROR: getMODIShobjid unble to obtain the id \n"
			"\t with an empty name input.\n");
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( NULLMODFIL(file) )	      
   {
     sprintf(buff,"ERROR: getMODIShobjid unable to obtain the id for\n"
			"\t object %.*s with an invalid MODIS file structure input.\n",
			MAX_NC_NAME,name);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( !( (type == MODIS_ARRAY) || (type == MODIS_TABLE) ) )
   {
     sprintf(buff, "ERROR: getMODIShobjid unable to obtain the id for\n"
			"\t object %.*s with an invalid MODIS object type: %ld\n",
			MAX_NC_NAME,name,type);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( !( (strcmp(access,"w") == 0) || (strcmp(access,"r") == 0) ) )
   {
     sprintf(buff, "ERROR: getMODIShobjid unable to obtain the id for\n"
			"\t object %.*s with an invalid access mode:%s\n",
			MAX_NC_NAME,name,access);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( (file->access == DFACC_READ) && (strcmp(access,"w") ==  0) )
   {
     sprintf(buff, "ERROR: getMODIShobjid unable to obtain the id of object\n"
			"\t %.*s for write while the file is opened for read only.\n",
			MAX_NC_NAME,name);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( type == MODIS_ARRAY )
     did = getMODISarrayid(file,name,group);
   else
     did = getMODIStableid(file,name,group,access);

   if ( did == NULL )
   {
     sprintf(buff, "ERROR: getMODIShobjid detected errors while attempting\n"
			"\t to obtain the id for object %.*s (Most likely\n"
			"\t the object does not exist).\n",
			MAX_NC_NAME,name);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   return((int)did->id);
}
