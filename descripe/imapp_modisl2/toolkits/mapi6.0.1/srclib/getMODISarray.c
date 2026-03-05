#include <stdlib.h> 
#include <stdio.h>
#include "mapic.h"

 int getMODISarray (MODFILE *file, char const *arrayname, char const *groupname,
                    long int start[], long int dimsizes[], void *data)

/*
!C****************************************************************************
*
*!Purpose:    Retrieves an array or subarray of data from a MODIS-HDF file
*	      array structure.
*					
*!Description: Function getMODISarray is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
* 	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      functions which allow MODIS Science Team-supplied software
*	      to read in Level 1B radiance bands and write out output 
*	      products and metadata to HDF files.  The functionality of the
*	      M-API is defined in the MODIS API User's Guide, Version 1, 
*	      dated 4/3/95.
*
*	      getMODISarray retrieves an array or sub-array of data from
*	      a MODIS-HDF file array structure.  The data array must be of 
*	      the same data type as the data in the target array structure.
*	      In addition, the dimensions and array region requested from
*	      the array structure must be consistent with the structure's 
*	      rank and dimensions.  (The array structure's data type, rank,
*	      and dimensions may be retrieved using getMODISarinfo).  If a
*  	      getMODISarray error occurs, the data retreival will not be
*	      performed.  See Chapter 6, "Accessing the Arrays" for more
*	      information.
*
*	      The groupname string provides the facility to select an array
*	      structure placed in a particular HDF Vgroup data group.  
*	      Alternatively, the entire file will be searched for an array
*	      structure named arrayname if groupname==NULL in C and a blank
*	      string, " " in FORTRAN.
*
*!Input Parameters:MODFILE *file	: MODISfile structure that references
*					  the MODIS-HDF file containing the
*					  target array structure
*		     char *arrayname	: target array name
*		     char *groupname	: data group name. If NULL, then the
*					  entire file is searched for arrayname
*		     long int start[]   : array structure locations to begin
*					  reading data from the array struct.
*					  Must have the same number of elements
*					  as the target array has dimensions
*	             long int dimsizes[]: array describing the size of the 
*					  array being received from the array
*					  structure.  dimsize must have the
*					  same number of elements as the 
*					  target array has dimensions and the
*					  product of the array dimensions must 
*					  equal the number of elements of data.
*						
*				
*!Output Parameters:void *data		: Output data
*
*		     Returns: MAPIOK if successful, MFAIL on error
*
*Externally defined:	
*			MODFILE				(mapi.h)
*			PGS_SMF_MAX_MSGBUF_SIZE 	(mapic.h)
*			NULLstr				(mapic.h)
*			MAPIERR				(mapic.h)
*			MAX_VAR_DIMS			(netcdf.h)
*			MAX_NC_NAME			(netcdf.h)
*			SDS_footprintOK 		(mapic.h)
*			SDreaddata			(mfhdf.h)
*			NULL				(stdio.h)
*			MFAIL				(mapi.h)
*			MAPIOK				(mapi.h)
*                       DATAID                          (mapi.h)
*                       SDSINFO                         (mapic.h)
*                       NULLMOFIL                       (mapic.h)
*                       getMODISarrayid.c               (mapi.h)
*
*
*!Revision History:
*		1996/02/01
*		Qi Huang
*		Original development for M-API version 2.0
*$Log: getMODISarray.c,v $
*Revision 5.1  2005/04/04 18:49:11  vlin
*constant safe for pointer arguments.
*
*Revision 1.1  1998/02/06 22:26:06  fshaw
*Initial revision
*
 * Revision 1.6  1996/08/06  12:28:30  fshaw
 * corrected E & W from pr:qa
 *
 * Revision 1.5  1996/08/06  12:19:06  fshaw
 * version 2.1
 *
 * Revision 1.4  1996/08/01  18:26:17  fshaw
 * version 2.0
 *
 * Revision 1.3  1996/05/08  17:45:48  qhuang
 * (This log message replaces "minor changes" below.) Inserted arrayname in
 * to SDS_footprintOK ERROR message.
 *
 * Revision 1.2  1996/02/07  22:40:52  qhuang
 * minor changes.
 *
 * Revision 1.1  1996/02/01  21:33:26  qhuang
 * Initial revision
 *
*
*!Team-unique Header:
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!References and Credits:
*              This software is developed by the MODIS Science Data Support
*	       Team for the National Aeronautics and Space Administration,
*	       Goddard Space Flight Center, under contract NAS5-32373.
*!Design Notes:
*
!END*********************************************************************
*/

   {
   char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
   char *funcname="getMODISarray"; 
   int32 sds_dimsizes[MAX_VAR_DIMS];	/* int32 HDF SD array		*/
   int32 sds_start[MAX_VAR_DIMS];	/* local array specifying the   */
                                        /*  array starting point        */
   int status;				/* function return status	*/
   int i;				/* loop counting variable	*/
   DATAID *did;
   SDSINFO *sinfo;
   
  /* Input checks: */
  if ( NULLstr(arrayname) ) {
    sprintf(buff,"ERROR: getMODISarray unable to read from an array\n"
			"\t without an array name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) ) {
    sprintf(buff,"ERROR: getMODISarray unable to read from the\n"
		"\t %.*s array with a NULL file MODFILE structure.\n",
			 MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( dimsizes == NULL ) {
    sprintf(buff,"ERROR: getMODISarray unable to read from the \n"
		"\t %.*s array without array dimension input.\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( start == NULL ) {
    sprintf(buff,"ERROR: getMODISarray unable to read from the \n"
		"\t %.*s array without array start input.\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( data == NULL ) {
    sprintf(buff,"ERROR: getMODISarray unable to read from the \n"
		"\t %.*s array without a data buffer.\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  status = MFAIL;

  did = getMODISarrayid( file, arrayname, groupname);

  if ( did  == NULL) {
     sprintf(buff, "ERROR: getMODISarray detected errors from"
                   "\t getMODSIarrayid while attempting to read from" 
                   "\t the %.*s array.\n", MAX_NC_NAME,arrayname);
     MAPIERR(buff,funcname);
  }
      
  else {
     sinfo = (SDSINFO *)did->info;
     if ( !SDS_footprintOK(sinfo->dimsizes,sinfo->rank,start,dimsizes)) {
        sprintf(buff, "ERROR: getMODISarray unable to read data from\n"
		    "\t invalid array structure locations in the\n"
                    "\t %.*s array.\n",MAX_NC_NAME,arrayname);
        MAPIERR(buff,funcname);
     }

     else {
        for(i=0; i < sinfo->rank; i++) {
          sds_start[i] = (int32) start[i];
          sds_dimsizes[i] = (int32) dimsizes[i];
        }

        if ( SDreaddata((int32)did->id, sds_start,NULL, sds_dimsizes, data) == FAIL){
          sprintf(buff, "ERROR: getMODISarray detected FAIL\n"
                        "\t from HDF procedure SDreaddata\n"
                        "\t while attempting to read from\n"
                        "\t the %.*s array.\n", MAX_NC_NAME,arrayname);
          MAPIERR(buff,funcname);
        }

	else
	  status = MAPIOK;	/* successful operation */
      } /* end of else;  for(i=0 ... */
  }

  return(status);
}
