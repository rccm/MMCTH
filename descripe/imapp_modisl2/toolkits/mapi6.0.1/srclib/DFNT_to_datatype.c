#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"

int DFNT_to_datatype (int32 dfnt, char *datatype)

/*
!C****************************************************************************
*
*Routine:     DFNT_to_datatype
*
*!Description: Subroutine DFNT_to_datatype is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
*	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      subroutines which allow MODIS Science Team-supplied software
*	      to read from/write to HDF format files. 
*
*	      DFNT_to_datatype receives the HDF number type constant
*             from the calling functions, determines the corresponding
*	      M-API data type string, and returns it to the calling function.
*	      If there no M-API data type string corresponds to the input HDF
*	      number type constant, then this routine returns with datatype
*	      pointing to NULL.  Note that the input datatype must be a string
*	      of at least 8 characters.
*
*	      This routine returns MAPIOK if successful, MFAIL on failure.
*
* !Input Parameters:	int dfnt	: HDF number type constant
*
* !Output Parameters:	char *datatype  : M-API data type string
*
*Externals:     	DFNT_INT8       (hdf.h)
*                       DFNT_UINT8      (hdf.h)
*			DFNT_INT16	(hdf.h)
*                       DFNT_UINT16	(hdf.h)
*			DFNT_INT32	(hdf.h)
*                       DFNT_UINT32     (hdf.h)
*               	DFNT_INT64      (hdf.h)
*                       DFNT_UINT64     (hdf.h)
*               	DFNT_FLOAT32    (hdf.h)
*               	DFNT_FLOAT64    (hdf.h)
*               	DFNT_CHAR8      (hdf.h)
*			I8		(mapi.h)
*                       UI8             (mapi.h)
*			I16             (mapi.h)
*                       UI16            (mapi.h)
*                       I32             (mapi.h)
*                       UI32            (mapi.h)
*                       I64             (mapi.h)
*                       UI64            (mapi.h)
*                       R32             (mapi.h)
*                       R64             (mapi.h)
*                       TXT             (mapi.h)
*                       strcpy          <string.h>
*                       MFAIL           (mapi.h)
*                       MAPIOK          (mapi.h)
*
* !Revision History:
* $Log: DFNT_to_datatype.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
*
*       2.    Paul Fisher               04 July 95
*             Research and Data Systems Corporation
*             SAIC/GSC MODIS Support Office
*             Revised to reflect updated PDL
*
*	1.    Dave Lorenzi/GSC		28 April 95
*	     
*	      Original Development / Testing
*
* !Team-unique Header:
*             Portions developed at the National Center for Supercomputing
*             Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !References and Credits:
*             This software is developed by the MODIS Science Data Support
*	      Team for the National Aeronautics and Space Administration,
*	      Goddard Space Flight Center, under contract NAS5-32373.
*
* !Design Notes:
*
!END*********************************************************************
*/
{
  int output_status=MAPIOK;                  /* return status */

  if (datatype==NULL)return(MFAIL);	     /* check input datatype */

   switch (dfnt)				/* Find the correct string*/
      {
      case DFNT_INT8    :  
	strcpy (datatype, I8);		        /* signed 8 bit integer	*/
	 break;
      case DFNT_UINT8   :  
         strcpy (datatype, UI8);	      /* unsigned 8 bit integer	*/
	 break;
      case DFNT_INT16   :  
         strcpy (datatype, I16);	       /* signed 16 bit integer	*/
         break;				 
      case DFNT_UINT16  :  
         strcpy (datatype, UI16);	     /* unsigned 16 bit integer	*/
         break;				 
      case DFNT_INT32   :  
	 strcpy (datatype, I32);	       /* signed 32 bit integer	*/
    	 break;	
      case DFNT_UINT32   :  
	 strcpy (datatype, UI32);	     /* unsigned 32 bit integer	*/
    	 break;	
      case DFNT_INT64   :  
	 strcpy (datatype, I64);	      /* signed 64 bit integer */
      	 break;
      case DFNT_UINT64   :  
	 strcpy (datatype, UI64);	    /* unsigned 64 bit integer */
      	 break;
      case DFNT_FLOAT32 :  
	 strcpy (datatype, R32);		/* 32 bit float		*/
	 break;
      case DFNT_FLOAT64 :  
	 strcpy (datatype, R64);		/* 64 bit float		*/
	 break;
      case DFNT_CHAR8   :  
	 strcpy (datatype, TXT);		/* character string	*/
	 break;
      default		:  
	 *datatype='\0';		       	/* No match, error 	*/ 
	 output_status=MFAIL; 
	 } 

   return(output_status);			/* Normal return	*/ 
   } 						/* End DFNT_to_dataype  */
