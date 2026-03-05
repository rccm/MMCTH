#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"


/*
!C****************************************************************************
*
*Routine:     short int emptyVdata (MODFILE *file, int32 vdata_ref) 
*					
*!Description: Function emptyVdata is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
* 	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      subroutines which allow MODIS Science Team-supplied software
*	      to read in Level 1B radiance bands and write out output 
*	      products and metadata to HDF files.
*
*             Function emptyVdata (EMPTV) indicates whether a Vdata is 
*             'empty' by searching for an HDF global metadata associated 
*             with the Vdata. An empty Vdata may not exist, so a dummy
*             record is inserted into it when it is first created. This 
*             dummy record should be overwritten with the first call from
*             putMODIStable. When the dummy record is overwritten by a 
*             single record the content of the Vdata is ambiguous, so a 
*             metadata ('NOT EMPTY') is associated with the Vdata to 
*             indicate that the Vdata is no longer empty. Function
*             emptyVdata searches for the existence of this label and reads
*             it as evidence that the dummy record has been oveerwritten.
*             This check is performed ONLY if the Vdata contains one 
*             record; multiple records in a Vdata implies that the dummy
*             record has been overwritten and no semaphore metadata
*             needs to be checked (or written).
*
*
*	      NOTE :  An assumption is made that int32 arrays are equivalent 
*		      to long int arrays.
* 
* !Input Parameters:MODFILE *file	: Address of MODISfile structure
*					  that references the MODIS-HDF file
*                                         receiving the new table.
*		      vdata_ref 	: Reference number of the Vdata to
*                                         interrogate.
*						
* !Output Parameters:NONE
*
*Returns: 1 if the Vdata is empty, 0 if it contains data,
*                             MFAIL (-1), if an error occurs.
*
*Externally defined:	SDfindattr      (mfhdf.h) 
*                       Vdatattr_name   (mapic.h)
*                       VSNAMELENMAX    (hdf.h)
*                       MFAIL           (mapi.h)
*                       DATA_LABEL      (mapic.h)
*
* !Revision History:
*
*	1.    Mitchell Weiss/RDC		07/28/95
*	     
*	      Original Development / Testing
*
*       2.   Xiao-Yang Ding/RDC                 07/28/95
*
*            Revision after walk through.
* !Team-unique Header:
*
* !References and Credits:
*              This software is developed by the MODIS Science Data Support
*	       Team for the National Aeronautics and Space Administration,
*	       Goddard Space Flight Center, under contract NAS5-32373.
*
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
!END*********************************************************************
*/

   short int  emptyVdata (MODFILE *file, int32 vdata_ref)
   {
   int   status=MFAIL;                          /* status set to MFAIL */
   char attribute_name[VSNAMELENMAX] = "";    /* attribute name assigned to
                                                   a particular Vdata       */
						/* check input for errors   */
   if (Vdatattr_name(DATA_LABEL,vdata_ref,attribute_name) != MAPIOK) 
      return status;
                                                /* check n_elements         */
   if (SDfindattr((int32)file->sd_id,attribute_name) == FAIL)
      status=1;                                 /* (the Vdata is empty)     */
   else
      status=0;                                 /* attribute exists and the */
   return status;
   }                                           /* end emptyVdata            */
