#include "mapic.h" 

int getMODISardims(MODFILE *file, char const *arrayname, char const *groupname,
                   char *data_type, long int *rank, long int dimsizes[])
/*
!C**********************************************************************
* 
*!Purpose: Retrieves the rank, dimensions and data type of a specified
*          MODIS-HDF file array structure.
* 
*!Description: Subroutine getMODISardims (GMARDM) is part of a larger software 
*              system called the MODIS Applications Programming Interface (API)
*              Utility, abbreviated M-API.  The M-API Utility consists of
*              subroutines which allow MODIS Science Team-supplied software
*              to read in Level 1B radiance bands and write out output
*              products and metadata to HDF files.  The functionality of the
*              M-API is defined in the MODIS API User's Guide, Version 1,
*              dated 4/3/95.
*
*              getMODISardims retrieves the essential characteristics
*              of an HDF SDS array structure contained in a MODIS-HDF
*              file.  This provides the information needed for properly
*              reading data from the array structure using getMODISarray
*              (GMAR). 
*
*              The groupname string provides the facility to select an array
*              structure placed in a particular HDF 'Vgroup' data group.
*              Alternatively, the entire file will be searched for an array
*              structure named arrayname if groupname=NULL.
*
*              Proper dimensioning of dimsizes to provide sufficient elements
*              for the dimensions of the array structure may at first appear
*              to require precognition.  The easiest solution is to provide a 
*              generous dimsizes array. Another approach is to use the rank
*              variable as an input containing the number of elements in 
*              dimsizes.  If dimsizes is inadequate for the multi-dimensional 
*              array structure in question, getMODISardims will fail gracefully
*              but will return the rank of the array structure, allowing for the
*              dimension information to be retrieved on the second call.  
* 
* 
* !Input Parameters:
*
*
*     file	       IN:  Address of MODFILE structure that is used to reference 
*                      the MODIS-HDF file containing the target array structure
* 
*     arrayname	       IN:  ASCII string name of the target array structure
*     groupname	       IN:  ASCII string name of the data group containing the target
*                      array structure. If set to NULL the entire file will be 
*                      searched for the array structure name arrayname.
*
* !Output Parameters:
*
*     data_type        OUT:  String describing the data type of the array.
*                        permitted types: "int8"
*                                       "int16"
*                                       "int32"
*                                       "int64"
*                                       "float32"
*                                       "float64"
*                                       "char *"
*
*     rank	       IN/OUT:   (Optional)The number of elements in the array 
		       dimsizes on input.  This is replaced with the number of 
		       dimensions in the array structure for output.  It is
		       set to 0 if a functional error occurs.  No dimensional 
		       information will be provided if rank = NULL.
*     dimsizes	       OUT:  Array describing the size of each dimension of the
*                      array structure. The dimensions will not be provided
		       unless dimsizes contains sufficient elements for the
		       rank of the array.  If the dimsizes is too small,
		       the array rank information will be returned in the 
		       rank argument and the function will return MFAIL.
*
* Externals:
*     MODFILE				(mapi.h)
*     PGS_SMF_MAX_MSGBUF_SIZE		(PGS_SMF.h)
*     DATAID				(mapic.h)
*     SDSINFO				(mapic.h)
*     MAX_NC_NAME			(netcdf.h)
*     NULLstr				(mapic.h)
*     MAPIERR				(mapic.h)
*     MAPIOK				(mapi.h)
*     MFAIL				(mapi.h)
*     NULLMODFIL			(mapic.h)
*     getMODISarrayid			(mapic.h)
*     DFNT_to_datatype			(mapic.h)
*     
* !Revision History:
*		Qi Huang	1996/07/23
*		Version 2.1
*		Ring super structure and other changes make this
*		version much faster.
*
*		Qi Huang  1996/02/08
*		Version 2.0
*
* $Log: getMODISardims.c,v $
* Revision 5.1  2005/04/04 18:47:54  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.6  1996/08/15  18:31:38  qhuang
 * Version 2.1
 *
 * Revision 1.5  1996/03/14  14:57:56  qhuang
 * removed redundant expression (ret == MFAIL) in if condition.
 *
 * Revision 1.4  1996/03/13  21:59:14  qhuang
 * removed array brackets from sds_name initialization,
 * cast setting of in_rank.
 *
 * Revision 1.3  1996/02/29  22:40:37  qhuang
 * Check the return value of DFNT_to_datatype.
 *
 * Revision 1.2  1996/02/09  18:04:49  qhuang
 * Version2.0  Added group handling, changed messages, separated the input checks and
 * changed the return of rank.
 *
 * Revision 1.1  1996/02/09  18:01:49  qhuang
 * Initial revision
 *
 * Revision 1.1  1996/02/09  18:01:49  qhuang
 * Initial revision
 *
 * Revision 1.2  1995/10/30  21:33:07  qhuang
 * Added capability to pass status messages to log files.
 *
 * Revision 1.1  1995/10/30  21:30:54  qhuang
 * Initial revision
 *
 * Revision 1.1  1995/10/30  21:30:54  qhuang
 * Initial revision
 *
* 
*      Revision 01.00 1995/04/19
*      Angela Sigmund 
*      Original development. 
*  
*      Revision 01.01 1995/04/25
*      Angela Sigmund 
*      Added minor comments, deleted ret stmt at SDselect call.
*
*      Revision 01.02 1995/05/05
*      Angela Sigmund
*      Added description of externals, removed check for
*      groupname!=NULL and set 'long int' to int32, moved
*      location of SDendaccess call. 
*
*      Revision 01.03 1995/10/25
*      Qi Huang
*      Added error message in case when users request an array
*      with an incorrect "rank" specification.
* 
*      Revision 01.04 1995/10/30
*      SMF(Status Message Files) integration.
*
* !Team-unique Header:
* 
* !References and Credits:
*      This software is developed by the MODIS Science Data Support
*      Team for the National Aeronautics and Space Administration,
*      Goddard Space Flight Center, under contract NAS5-32373.
*
*      HDF portions developed at the National Center for Supercomputing
*      Applications at the University of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END********************************************************************
*/

{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char *funcname="getMODISardims"; 
  
  int32 in_rank = 0;
  DATAID	*did;
  SDSINFO	*sinfo;
   
  int ret;                  /* status code to be returned by routine */
  int  i;                  /* loop counter */
  
  if ( rank != NULL )
  {
    in_rank = (int32)*rank;
    *rank = 0;
  }
  
  /* Input checks: */
  if ( NULLstr(arrayname) )
  {
    sprintf(buff,"ERROR: getMODISardims unable to access an array\n"
			"\t without an array name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) )
  {
    sprintf(buff,"ERROR: getMODISardims unable to access the %.*s\n"
			"\t array with a NULL file MODFILE structure.\n",
		 	MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( ( rank != NULL ) && ( dimsizes == NULL ) )
  {
    sprintf(buff,"ERROR: getMODISardims unable to return the %.*s\n"
			"\t array's dimensions without a dimsizes array.\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  ret=MFAIL;	 /* error */

  did = getMODISarrayid(file,arrayname,groupname);
  if ( did == NULL )
  {
    sprintf(buff,"ERROR: getMODISardims detected errors from\n"
			"\t getMODISarrayid while attempting to access the\n"
			"\t %.*s array.\n",MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
  }
  else
  {
    sinfo = (SDSINFO *)did->info;
    /* (Note: sinfo contains SDS's rank, numerical type, and dimension sizes) */
    ret = MAPIOK;

    if ( data_type != NULL )
      if ( DFNT_to_datatype(sinfo->ntype,data_type) == MFAIL )
      {
	sprintf(buff,"ERROR: getMODISardims does not recognize\n"
			"\t HDF number type %ld which\n"
			"\t is not implemented in M-API while\n"
			"\t attempting to access the %.*s\n"
			"\t array.\n",(long)sinfo->ntype,MAX_NC_NAME,arrayname);
	MAPIERR(buff,funcname);
	ret = MFAIL;
      }

    if ( (rank != NULL) && (ret == MAPIOK) )
    {
      if ( in_rank >= sinfo->rank )
	for (i=0; i < sinfo->rank; i++)
	  dimsizes[i] = sinfo->dimsizes[i];
     
      else
      {
	sprintf(buff,"ERROR: getMODISardims unable to return\n"
			"\t %.*s array's %ld dimension\n"
			"\t sizes in a %ld element dimsizes\n"
			"\t array.\n",MAX_NC_NAME,arrayname,
			(long)sinfo->rank, (long)in_rank);
	MAPIERR(buff,funcname);
	ret = MFAIL;
      }

      *rank = sinfo->rank;

    } /* end of if */
  }/* end of else */

  return(ret);
}
