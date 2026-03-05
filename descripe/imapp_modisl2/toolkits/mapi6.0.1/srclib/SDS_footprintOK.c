#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "mapic.h"

short int SDS_footprintOK(int32 const sds_dimsizes[], int32 sds_rank, 
                          long int const start[], long int const dimsizes[])
/*
!C**********************************************************************
* 
*!Purpose:	Determines whether an array access will fall within the 
* boundaries of the target SDS.
* 
*!Description: Function SDS_footprintOK is part of a larger software system 
* called the MODIS Applications Programming Interface (API) 
* Utility, abbreviated M-API. The M-API Utility consists of 
* subroutines which allow MODIS Science Team-supplied software 
* to read  and write data and metadata from/to HDF files. The 
* functionality of the M-API is defined in the MODIS Application 
* Program Interface (API) Specification.
* 
* SDS_footprintOK checks the start position and dimensions of a 
* hyperslab to be accessed from an SDS against the dimensions of 
* that SDS.
* 
* !Input Parameters:
* 
* sds_dimsizes	IN:	Array containing the size of the SDS in each 
*		dimension.
* sds_rank	IN:	number of dimensions in the SDS. 
* 		sds_dimsizes, start, and dimsizes all have sds_rank 
* 		elements
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
* !Output Parameters:
*		None
* 
* Returns:	1 if proposed hyperslab consistent with SDS dimensions, 0 if not 
* or an error occurs.
*
* External refereces:
*		PGS_SMF_MAX_MSGBUF_SIZE		(mapic.h)
*		MAPIERR				(mapic.h)
* 
* !Revision History:
*		Qi Huang	1996/07/22
*		Version 2.1
*		Original development and testing
* 		Ring super structure and other changes make
*		this version much faster.
*
* $Log: SDS_footprintOK.c,v $
* Revision 5.1  2005/04/04 16:50:23  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.1  1996/07/31  17:31:06  qhuang
 * Initial revision
 *
*
* !Team-unique Header:
* This software is developed by the MODIS Science Data SupportTeam for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
*!References and Credits
*
*!Design Notes
*
!END********************************************************************
*/
{ char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning
                                           message */
  char *funcname="SDS_footprintOK";       /* name of this module */
  short int status;
  int		i;

  /* Input checks: */
  if ( (sds_dimsizes == NULL) || (dimsizes == NULL) || (start == NULL) )
    return(0);

  status = 1;		/* proposed access to SDS OK */

  for (i=0; (i<sds_rank) && (status == 1); i++)
  {
    if ( (start[i] < 0) || (start[i] > (sds_dimsizes[i]-1)) ) /* inclusive */
    {
      sprintf(buff,"ERROR: SDS_footprintOK found unable to access\n"
			"\t data to array structure location \"%ld\n"
			"\t ...%ld\".\n",start[0],start[sds_rank-1]);
      MAPIERR(buff,funcname);
      status = 0;	/* access to SDS NOT advised */
    }
    else if ( ((start[i] + dimsizes[i]) < 1) || 
	      ((start[i] + dimsizes[i]) > sds_dimsizes[i]) ) /* inclusive */
      status = 0;	/* access to SDS NOT advised */
  }
  
  return(status);
}



  
