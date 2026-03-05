#include <stdlib.h>  
#include <stdio.h>
#include <string.h>
#include "mapic.h"

DATAID *getMODIStableid(MODFILE *file, char const *tablename, 
                        char const *groupname, char const *access)
/*
!C**********************************************************************
* 
*!Purpose:	get the address of the DATAID structure for a MODIS table (Vdata)
* 
*!Description: 	Function getMODIStableid is part of a larger software system 
*		called the MODIS Applications Programming Interface (API) 
*		Utility, abbreviated M-API.  The M-API Utility consists of 
*		subroutines which allow MODIS Science Team-supplied software 
*		to read  and write data and metadata from/to HDF files. The 
*		functionality of the M-API is defined in the MODIS Application 
*		Program Interface (API) Specification.
*		
*		getMODIStableid returns the address of DATAID structure for  
*		a MODIS table. The DATAID structure contains the HDF access ID as well as 
*		vdata specific information. This is an internal M-API 
*		routine used by M-API table access routines.
*
*		
*		Because this routine is only used by M-API internally, and all 
*		inputs of this routine have been checked in the calling routines, 
*		this routine does not performace the input parameter checking.
* 
* !Input Parameters:
* file		IN: 	Address of MODFILE structure that is used to 
* 		reference a MODIS HDF file.
* tablename	IN:	ASCII string name of the table structure
* 		to be searched.
* groupname	IN:	ASCII string name of the data group 
* 		containing the table
* access	IN:	ASCII string of the table access type: "r"
*		for read only, "w" for write.
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
*		VSNAMELENMAX			(hdf.h)
*		VGNAMELENMAX			(hdf.h)
*		DFTAG_VH 			(hdf.h)
*		MVSfind				(mapic.h)
*		VSattach			(vproto.h)
*		VSdetach			(vproto.h)
*		addid				(mapic.h)
*
* !Revision History:
* 		Qi Huang	1996/09/17
*		Version 2.2
*		Original development and testing
*		Ring super structure and other changes make
*		this version much faster.
* $Log: getMODIStableid.c,v $
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
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
  char *funcname="getMODIStableid";       /* name of this module */

  DATAID	*did;
  int32		ref_id;
  int32		vdata_id;

  if ( NULLstr(groupname) )
    groupname = NULL;

  did = searchid(file,tablename,groupname,DFTAG_VH,access);

  if ( did == NULL )
  {
    if ( groupname )
    {
      if ( (ref_id = searchMODISgroup(file,groupname,NULL,tablename,
					 NULL,DFTAG_VH)) == MFAIL )
      {
	sprintf(buff,"ERROR: getMODIStableid unable to find the %.*s\n"
				"\t data group containing the %.*s table.\n",
				VGNAMELENMAX,groupname,VSNAMELENMAX,tablename);
	MAPIERR(buff,funcname);
	return(NULL);
      }
      else if ( ref_id == NO_OBJECT )
      {
	sprintf(buff,"ERROR: getMODIStableid cannot find the %.*s table\n"
				"\t in the %.*s data group.\n",
				VSNAMELENMAX,tablename,VGNAMELENMAX,groupname);
	MAPIERR(buff,funcname);
	return(NULL);
      }
    }
    else if ( (ref_id = MVSfind((int32)file->hdf_id,tablename)) == FAIL )
    {
      sprintf(buff,"ERROR: getMODIStableid cannot find %.*s table.\n",
				VSNAMELENMAX,tablename);
      MAPIERR(buff,funcname);
      return(NULL);
    }

    if ( (vdata_id = VSattach((int32)file->hdf_id,ref_id,access)) == FAIL )
    {
      sprintf(buff,"ERROR: getMODIStableid detected FAIL from HDF\n"
			"\t procedure VSattach attempting to access\n"
			"\t the %.*s table.\n",VSNAMELENMAX,tablename);
      MAPIERR(buff,funcname);
      return(NULL);
    }

    if ( addid(file,tablename,groupname,vdata_id,DFTAG_VH,access) == MFAIL )
    {
      VSdetach(vdata_id);
      sprintf(buff,"ERROR: getMODIStableid detected FAIL from M-API\n"
			"\t internal routine addid while attempting to\n"
			"\t access the %.*s table.\n",VSNAMELENMAX,tablename);
      MAPIERR(buff,funcname);
      return(NULL);
    }

    did = searchid(file,tablename,groupname,DFTAG_VH,access);
  }/* end of if; if ( did == NULL) */

  return(did);
}
