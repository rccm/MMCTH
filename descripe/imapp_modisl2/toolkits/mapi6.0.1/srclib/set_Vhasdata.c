#include <stdlib.h>
#include <stdio.h>
#include "mapic.h"


short int set_Vhasdata (MODFILE *file, int32 vdata_ref) 

/*
!C****************************************************************************
*
*Routine:     short int set_Vhasdata(MODFILE *file, int32 vdata_ref)
*
*!Purpose:     Creates a semaphore attribute to indicate that the first record
*             of particular Vdata is not a dummy record.
*					
*!Description: Function set_Vhasdata is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
* 	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      functions which allow MODIS Science Team-supplied software
*	      to read in Level 1B radiance bands and write out output 
*	      products and metadata to HDF files.  The functionality of the
*	      M-API is defined in the MODIS API User's Guide, Version 1, 
*	      dated 4/3/95.
*
*             set_Vhasdata (HASDAT) sets a semaphore HDF global attribute
*             associated with a Vdata to indicate whether the Vdata is not
*             'empty'. An empty Vdata may not exist, so a dummy record is
*             inserted into it when it is first created. This dummy record
*             should be overwritten with the first call from putMODIStable.
*             When the dummy record is overwritten by a single record the
*             content of the Vdata is ambiguous, so this routine creates
*             a metadata ('NOT EMPTY') associated with the Vdata to indicate
*             that the Vdata is no longer empty. This routine is called ONLY
*             if the Vdata contains one non-dummy record; multiple records
*             in a Vdata implies that the dummy record has been overwritten
*             and no semaphore metadata needs to be checked (or written).
*
* !Input Parameters:
*             file	IN: Address of MODFILE structure that is used to
*			reference the MODIS-HDF file receiving the new 
*			table.
*             vdata_ref	IN: Reference number of the Vdata to interrogate.
*
* !Output Parameters:	NONE.
*
*Returns:     MAPIOK if the semaphor attribute is written successfully,
*             MFAIL if an error occurs.
*
*Externally defined:   
*			Vdatattr_name	(mapic.h)
*			SDsetattr	(mfhdf.h)
*			NULL		(stdio.h)
*			MFAIL		(mapi.h)
*			MAPIOK		(mapi.h)
*			DATA_LABEL	(mapic.h)
*			VSNAMELENMAX	(hdf.h)
*
* !Revision History:
* $Log: set_Vhasdata.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
*
*	1.    Baochun Ge/GSC		8/15/95
*	     
*	      Original Development / Testing
*
* !Team-unique Header:
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !References and Credits:
*              This software is developed by the MODIS Science Data Support
*	       Team for the National Aeronautics and Space Administration,
*	       Goddard Space Flight Center, under contract NAS5-32373.
* !Design Notes:
*
!END********************************************************************
*/
{
   char attribute_name[VSNAMELENMAX];	/* global attribute name	*/
   int status;				/* status code			*/

/* Use Vdatattr_name to retrieve the global attribute name of the metadata DATA_LABEL
   associated with the Vdata having the vdata_ref reference id. The output string must
   be at least VSNAMELENMAX elements long.
*/
   status = Vdatattr_name(DATA_LABEL, vdata_ref, attribute_name);
   if(status == FAIL)
      return(MFAIL);

/* Create a global attribute with SDsetattr with the global attribute name and the
   string value DATA_LABEL */
   status = SDsetattr((int32)file->sd_id, attribute_name, (int32)DFNT_CHAR8,
			 (int32)strlen(DATA_LABEL), DATA_LABEL);

/* IF SDsetattr FAILS to create the global attribute, RETURN MFAIL */
   if(status == FAIL)
      return(MFAIL);
   else
      return(MAPIOK);
}

/* End set_Vhasdata.c */

