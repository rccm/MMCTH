#include <stdlib.h>  
#include <stdio.h>
#include <string.h>
#include "mapic.h"

DATAID *getMODISarrayid(MODFILE *file, char const *arrayname, 
                        char const *groupname)
/*
!C**********************************************************************
* 
*!Purpose:	get the address of the DATAID structure for a MODIS array (SDS) 
* 
*!Description: 	Function getMODISarrayid is part of a larger software system 
*		called the MODIS Applications Programming Interface (API) 
*		Utility, abbreviated M-API.  The M-API Utility consists of 
*		subroutines which allow MODIS Science Team-supplied software 
*		to read  and write data and metadata from/to HDF files. The 
*		functionality of the M-API is defined in the MODIS Application 
*		Program Interface (API) Specification.
*		
*		getMODISarrayid returns the address of DATAID structure for the 
*		array. The DATAID structure contains the HDF access ID as well as 
*		rank and dimension information. This is an internal M-API 
*		routine used by M-API array access routines. It is suggested that 
*		getMODISarray and putMODISarray should not call but in-line the 
*		code of this routine since the performance of getMODISarray and 
*		putMODISarray is critical for the overall performace of M-API.
*		
*		Because this routine is only used by M-API internally, and all 
*		inputs of this routine have been checked in the calling routines, 
*		this routine does not performace the input parameter checking.
* 
* !Input Parameters:
* file	IN: 	Address of MODFILE structure that is used to 
* 		reference a MODIS HDF file receiving the dimension 
* 		name.
* arrayname	IN:	ASCII string name of the target array 
* 		structure.
* groupname	IN:	ASCII string name of the data group 
* 		containing the target
* 
* !Output Parameters:		NONE
* 
* Returns:	(DATAID *)  if successful, NULL if an error occurs.
*
* External reference:
* 
*		MODFILE				(mapi.h)
*		PGS_SMF_MAX_MSGBUF_SIZE		(mapic.h)
*		DATAID				(mapi.h)
*		NULLstr				(mapic.h)
*		searchid			(mapic.h)
*		searchMODISgroup		(mapic.h)
*		NO_OBJECT			(mapic.h)
*		MAPIERR				(mapic.h)
*		MFAIL				(mapi.h)
*		MAX_NC_NAME			(netcdf.h)
*		VSNAMELENMAX			(hdf.h)
*		DFTAG_NDG			(hdf.h)
*		MAPIOK				(mapi.h)
*		MAPIERR				(mapic.h)
*		SDnametoindex			(mfhdf.h)
*		SDselect			(mfhdf.h)
*		addid				(mapic.h)
*		SDendaccess			(mfhdf.h)
*		searchid			(mapic.h)
*
* !Revision History:
* 		Qi Huang	1996/07/22
*		Version 2.1
*		Original development and testing
*		Ring super structure and other changes make
*		this version much faster.
* $Log: getMODISarrayid.c,v $
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.3  1996/08/09  16:38:25  qhuang
 * Cast some arguments in HDF routine calls.
 *
 * Revision 1.2  1996/07/30  21:29:17  qhuang
 * Added NCSA acknowledgement in prolog.
 *
 * Revision 1.1  1996/07/30  21:28:48  qhuang
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
  char *funcname="getMODISarrayid";       /* name of this module */

  DATAID	*did;
  char		*access="";
  int32		sds_index;
  int32		sds_id;

  if ( NULLstr(groupname) )
    groupname = NULL;

  did = searchid(file,arrayname,groupname,DFTAG_NDG,access);

  if ( did == NULL )
  {
    if ( groupname )
    {
      if ( (sds_index = searchMODISgroup(file,groupname,NULL,arrayname,
					 NULL,DFTAG_NDG)) == NO_OBJECT )
      {
	sprintf(buff,"ERROR: getMODISarrayid cannot find the %.*s array\n"
				"\t in the %.*s data group.\n",
				MAX_NC_NAME,arrayname,VSNAMELENMAX,groupname);
	MAPIERR(buff,funcname);
	return(NULL);
      }
      else if ( sds_index == MFAIL )
      {
	sprintf(buff,"ERROR: getMODISarrayid unable to find the %.*s\n"
				"\t data group containing the %.*s array.\n",
				VSNAMELENMAX,groupname,MAX_NC_NAME,arrayname);
	MAPIERR(buff,funcname);
	return(NULL);
      }
    }
    else if ( (sds_index = SDnametoindex((int32)file->sd_id,arrayname)) == FAIL )
    {
      sprintf(buff,"ERROR: getMODISarrayid cannot find array %.*s\n",
				MAX_NC_NAME,arrayname);
      MAPIERR(buff,funcname);
      return(NULL);
    }

    if ( (sds_id = SDselect((int32)file->sd_id,sds_index)) == FAIL )
    {
      sprintf(buff,"ERROR: getMODISarrayid detected FAIL from HDF\n"
			"\t procedure SDselect attempting to access\n"
			"\t the %.*s array.\n",MAX_NC_NAME,arrayname);
      MAPIERR(buff,funcname);
      return(NULL);
    }

    if ( addid(file,arrayname,groupname,sds_id,DFTAG_NDG,access) == MFAIL )
    {
      SDendaccess(sds_id);
      sprintf(buff,"ERROR: getMODISarrayid detected FAIL from M-API\n"
			"\t internal routine addid while attempting to\n"
			"\t access the %.*s array.\n",MAX_NC_NAME,arrayname);
      MAPIERR(buff,funcname);
      return(NULL);
    }

    did = searchid(file,arrayname,groupname,DFTAG_NDG,access);
  }/* end of if; if ( did == NULL) */

  return(did);
}
