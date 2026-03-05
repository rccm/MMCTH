#include <stdlib.h>
#include <stdio.h>
#include "mapic.h"

int closeMODISfile (MODFILE **file)
/*
!C****************************************************************************
*
*!Description: Subroutine closeMODISfile is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
*	     Utility, abbreviated M-API.  The M-API Utility consists of 
*	     subroutines which allow MODIS Science Team-supplied software 
*	     to read in Level 1B radiance bands and write out output 
*	     products and metadata to HDF files.  The functionality of the
*	     M-API is defined in the MODIS API User's Guide, Version 1, 
*	     dated 4/3/95.
*
*	     closeMODISfile terminates the access of MODIS API routines
*	     to a MODIS HDF file opened using openMODISfile.  Only
*	     pre-existing files should be closed by closeMODISfile.
*	     completeMODISfile should be used to end access to a new
*	     MODIS HDF file so that the file's header information can 
*	     be completed.  closeMODISfile may fail to close the file 
*	     if an error occurs.
*
*
* !Input Parameters:	MODFILE **file      : MODIS file structure
*
* !Output Parameters:	NONE
* 
*Returns:               MAPIOK if successful, MFAIL if an error
*                       occurs.
*
*Externally defined:	MODFILE		       (mapi.h)
*			MODIS_ALL_TYPES	       (mapi.h)
*			endMODISobjaccess      (mapi.h)
*                       MAPIERR                (mapic.h)
*                       MFAIL                  (mapi.h)
*			MAPIOK		       (mapi.h)
*			NULLMODFIL	       (mapic.h)
*			SDend		       (libdf.a)       
*                       DFACC_CREATE	       (hdf.h) 	
*                       FAIL                   (hdf.h)
*                       PGS_SMF_MSGBUF_SIZE    (mapic.h)
* !Revision History:
*		Qi Huang	07/09/1996
*		Version 2.1. Removed the special case for FORTRAN.
*
* $Log: closeMODISfile.c,v $
* Revision 6.2  2010/07/12 16:23:17  kuyper
* Made argument conversion explicit.
*
* Revision 6.1  2010/04/20 19:16:45  kuyper
* Changed to invert the changes made to openMODISfile() for compatibility with
*   versions of HDF greater than 4.2r1.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.7  1997/05/20  21:22:15  qhuang
 * Added 'file' in if check for MODFILE structure.
 *
 * Revision 1.6  1996/10/15  19:52:55  qhuang
 * In SMF message, changed MAPIERR to MAPIWARN.
 *
 * Revision 1.5  1996/08/09  18:42:20  qhuang
 * Removed unused variables in variable declaration.
 *
 * Revision 1.4  1996/07/16  18:52:07  qhuang
 * Ring super structure and other changes make this version much faster.
 *
 * Xia W. Chang  1995/12/27  xia@modis-xl.gsfc.nasa.gov
 * version 2.0 design. (Using hdf_id associated with sd_id instead
 * of reopening the hdf file)  
 *
 * Revision 1.1  1995/11/29  17:22:56  qhuang
 * Initial revision
 *
 * Revision 1.3  1995/11/01  21:53:51  qhuang
 * minor changes
 *
 * Revision 1.2  1995/10/30  22:39:30  qhuang
 * Added capability to pass status messages to log files.
 *
*
*       1.    Revision 01.00  1995/04/12   0800 EDT
*	      Jenny Glenn/GSC  (jglenn@modis-xl.gsfc.nasa.gov)
*
*	      Create the initial template for MODIS API utility subroutine
*	      closeMODISfile.
*
*	2.    Dave Lorenzi/GSC		28 Apr 95
*
*	      Original Development / Testing
*
*       3.    Vicky Lin / RDC
*             Check the return status of function "Hfidinquire" 
*             using "FAIL" instead of "NULL".
*
* !Team-unique Header:
*	      Portions developed at the National Center for Supercomputing
*	      Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !References and Credits:
*	      This software is developed by the MODIS Science Data Support
*	      Team for the National Aeronautics and Space Administration,
*	      Goddard Space Flight Center, under contract NAS5-32373.
*
* !Design Notes:
*
!END*********************************************************************
*/
{ 
   char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
   char *funcname="closeMODISfile"; 
   int ret= MFAIL;  		/* return code, 0=ok, -1=error	   */
   int MAXLEN = 100;

   if ( file && (!NULLMODFIL(*file)) )	      
   {
     if ( endMODISobjaccess(*file,NULL,NULL,MODIS_ALL_TYPES) == MFAIL )
     {
       sprintf(buff,"ERROR: closeMODISfile detects errors in\n"
			"\t endMODISobjaccess while trying to close\n"
			"\t all opend objects for file %.*s.\n",
			MAXLEN, (*file)->filename);
       MAPIERR(buff,funcname);
       return(MFAIL);
     }

     if ( SDend((int32)(*file)->sd_id) == FAIL ) 
     {
       sprintf(buff,"ERROR: closeMODISfile detected FAIL from HDF\n"
           		"\t procedure SDend, unable to close file %.*s.\n",
                	MAXLEN, (*file)->filename);
       MAPIERR(buff,funcname);
     }
    else if(Vend((HFILEID)(*file)->hdf_id) != SUCCEED)
    {
	sprintf(buff, "ERROR: Vend(\"%.*s\") failed.\n",
	    MAXLEN, (*file)->filename);
	MAPIERR(buff,funcname);
    }
    else if(Hclose((int32)(*file)->hdf_id) != SUCCEED)
    {
	sprintf(buff,"ERROR: Hopen(\"%.*s\") failed.\n",
	    MAXLEN, (*file)->filename);
	MAPIERR(buff,funcname);
    }
     else
     {
       if ((*file)->access == DFACC_CREATE )     /* Close new file  */
       {
 	 sprintf(buff,"WARNING: closeMODISfile closed new file\n"

			"\t %.*s without complete header information.\n",
			MAXLEN, (*file)->filename);
 	 MAPIWARN(buff,funcname);
       }

       if ((*file)->filename) free((*file)->filename); 
       free(*file);           /* file data structure           */
       *file = NULL;
      
       ret = MAPIOK;
     }
   }/* end of if; if ( !NULLMODFIL(*file) ) */

   else
   {
     sprintf(buff, "ERROR: closeMODISfile can not close a NULL file\n");
     MAPIERR(buff, funcname);
   }

   return(ret);

}
