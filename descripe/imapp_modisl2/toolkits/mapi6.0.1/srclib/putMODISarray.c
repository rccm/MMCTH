#include <stdio.h> 
#include <stdlib.h>
#include "mapic.h"

int putMODISarray(MODFILE *file, char const *arrayname, 
                  char const *groupname, long int const start[], 
                  long int const dimsizes[], void const *data)
/*
!C****************************************************************************
*
*!Purpose:	Inserts an array or subarray of data into a MODIS-HDF file array 
* 		structure.
* 
*!Description:  Function putMODISarray is part of a larger software system called 
*		the MODIS Applications Programming Interface (API) Utility, 
*		abbreviated M-API. The M-API Utility consists of subroutines 
*		which allow MODIS Science Team-supplied software to read  and 
*		write data and metadata from/to HDF files. The functionality of 
*		the M-API is defined in the MODIS Application Program Interface 
*		(API) Specification.
*		
*		putMODISarray (PMAR) places a multi-dimensional array of data 
*		into an HDF SDS array structure previously created using 
*		createMODISarray (CRMAR).  The data in the array must be of the 
*		data type the target array structure was created for.  In addition, 
*		the dimensions and placement of the input array in the array 
*		structure must be consistent with the structureÕs rank and 
*		dimensions.  If a putMODISarray error message occurs, the data 
*		insertion will not be performed.  See Chapter 6, ÒAccessing 
*		ArraysÓ for additional information.  This routine may be called 
*		multiple times to fill the array structure.  Data previously in the 
*		array structure may be overwritten.
*		
*		The groupname  string provides the facility to select an array 
*		structure placed in a particular HDF ÔVgroupÕ data group.  The 
*		entire file will be searched for an array structure named 
*		arrayname  if groupname = NULL in C and a blank string (Ò Ò) in 
*		FORTRAN.
* 
* !Input Parameters:
* file	IN: 	Address of MODFILE structure that is used to 
* 		reference the MODIS-HDF file containing the target 
* 		array structure.
* arrayname	IN:	ASCII string name of the target array 
* 		structure.
* groupname	IN:	ASCII string name of the data group 
*		containing the target array structure.  If set to 
*		NULL the entire file will be searched for the array 
*		structure named arrayname.
* start	IN: 	Array containing the array structure 
*		location to begin placing the data into the array 
*		structure.  start  must have the same number of 
*		elements as the target array has dimensions.
* dimsizes	IN:	Array describing the size of the array being 
*		inserted into the array structure.  dimsizes  must 
*		have the same number of elements as the target 
*		array structure has dimensions and the product of 
*		the array dimensions must equal the number of 
*		elements in data.
* data	IN:	Address of the data buffer.
* 
* !Output Parameters:		NONE.
* 
* Returns:	MAPIOK if successful, MFAIL if an error occurs.
*
* External references: 
*			MODFILE				(mapic.h)
*			PGS_SMF_MAX_MSGBUF_SIZE		(mapic.h)
*			NULLMODFIL			(mapic.h)
*		        NULLstr          		(mapic.h)
*			MAX_VAR_DIMS	 		(netcdf.h)
*			MAX_NC_NAME	 		(netcdf.h)
*			getMODISarrayid			(mapic.h)
*			SDSINFO				(mapic.h)
*                       DFACC_READ       		(hdf.h)
*                       SDS_footprintOK  		(mapic.h)
*                       SDwritedata      		(mfhdf.h)
*                	(VOIDP)		                (hdfi.h)
*                       FAIL             		(hdf.h)
*			NULL		 		(stdio.h)
*			MFAIL		 		(mapi.h)
*			MAPIOK		 		(mapi.h)
* 
* Returns:     MAPIOK if successful
*              MFAIL  if an error occurs
*
* !Revision History:
* 		Qi Huang	1996/07/22
*		Version 2.1
*		Ring super structure and other changes make
*		this version much faster.
*
*		1996/02/01
* 		Qi Huang
*		Original development for M-API version 2.0
* 
* $Log: putMODISarray.c,v $
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.4  1996/08/09  15:22:11  qhuang
 * Casted some arguments in HDF routine calls, casted some variables in some
 * assignments, removed some unused variables in variable declaration.
 *
 * Revision 1.3  1996/07/30  21:38:40  qhuang
 * Removed redundant statement in prolog.
 *
 *
 * Revision 1.2  1996/05/23  20:25:27  qhuang
 * Limited string variable output inserted into messages.
 *
 * Revision 1.1  1996/02/01  20:13:49  qhuang
 * Initial revision
 *
*
* !Team-unique Header:
* This software is developed by the MODIS Science Data SupportTeam for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
* 
*
* !References and Credits:
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END*********************************************************************
*/

{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning
                                           message */
  char *funcname="putMODISarray";       /* name of this module */
  int32   sds_dimsizes[MAX_VAR_DIMS];
  int32   sds_start[MAX_VAR_DIMS];
  int     i,status;
  DATAID  *did;
  SDSINFO *sinfo;
  
  /* Input checks: */
  if ( NULLstr(arrayname) )
  {
    sprintf(buff,"ERROR: putMODISarray unable to write to an array\n"
			"\t without an array name input\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) )
  { 
    sprintf(buff,"ERROR: putMODISarray unable to write to the %.*s\n"
			"\t array with a NULL file MODFILE structure\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( dimsizes == NULL )
  {
    sprintf(buff,"ERROR: putMODISarray unable to write to the %.*s\n"
			"\t array without array dimension input.\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( start == NULL )
  {
    sprintf(buff,"ERROR: putMODISarray unable to write to the %.*s\n"
			"\t array without array start input.\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( data == NULL )
  {
    sprintf(buff,"ERROR: putMODISarray unable to write to the %.*s\n"
			"\t array without a data buffer\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  status = MFAIL;		/* error */

  if ( file->access == DFACC_READ )
  {
    sprintf(buff,"ERROR: putMODISarray unable to write to the %.*s\n"
			"\t array in file opened for read only.\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( (did = getMODISarrayid(file,arrayname,groupname)) == NULL )
  {
    sprintf(buff,"ERROR: putMODISarray detected errors from\n"
			"\t getMODISarrayid while attempting to write\n"
			"\t data to %.*s array.\n",MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
  }
  else
  {
    sinfo = (SDSINFO *)did->info;

    /* (Note: did->id is the sds id, and sinfo contains rank, dimsizes,
       and num_type for the SDS) */
    if ( !SDS_footprintOK(sinfo->dimsizes,sinfo->rank,start,dimsizes) )
    {
      sprintf(buff,"ERROR: putMODISarray unable to write data to\n"
			"\t invalid array structure locations in the\n"
			"\t %.*s array.\n",MAX_NC_NAME,arrayname);
      MAPIERR(buff,funcname);
    }
    else
    {
      for (i=0; i<sinfo->rank; i++)
      {
	sds_start[i] = (int32)start[i];
	sds_dimsizes[i] = (int32)dimsizes[i];
      }

      if ( SDwritedata((int32)did->id,sds_start,NULL,sds_dimsizes,
           (char *)data) == FAIL )
      {
	sprintf(buff,"ERROR: putMODISarray detected FAIL from\n"
			"\t HDF procedure SDwritedata while\n"
			"\t attempting to write data to %.*s array.\n",
			MAX_NC_NAME,arrayname);
	MAPIERR(buff,funcname);
      }
      else
	status = MAPIOK;	/* successful operation */
    }
  }/* end of else, else sinfo = (SDSINFO *)did->info; */

  return(status);
}
