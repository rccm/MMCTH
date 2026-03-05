#include <stdlib.h>  
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include "mapic.h"
#include "local_nc.h"

MODFILE *openMODISfile(
	char const	*filename,
	char const	*access
)

/*
!C****************************************************************************
*
*!Purpose:     Opens a file for access by the MODIS API routines.
*
*!Description: Function openMODISfile is part of a larger software system 
*	      called the MODIS Applications Programming Interface (API)
*	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      subroutines which allow MODIS Science Team-supplied software
*	      to read and write data and metadata from/to HDF files.
*	      The functionality of the M-API is defined in the MODIS 
*             Application Program Interface (API) specification.
*
*             openMODISfile opens a file and creates the HDF structures to
*	      support the MODIS API routines' access to it.  openMODISfile
*	      must be called to produce the MODFILE file structure before any
*	      of these C routines can access it.  Note that setting the
*	      file access to "w" creates a file and will overwrite a
*	      pre-existing one.  openMODISfile will close the file and 
*	      return null outputs if an error occurs.
*
* !Input Parameters:
*  filename       Complete path and filename for the file to be opened.
*  access         Standard C access mode.
*
*                 Permitted access modes:
*                   "r"  Open for read only.
*                   "w"  Create for read/write.
*                   "a"  Open for read/write.
*
* !Output Parameters:	NONE		
*
*Returns:     MODFILE structure that is used to reference the file in all
*             other MODIS API C routines, or NULL if an error occurs. 
*
*Externally defined:	PGS_SMF_MAX_MSGBUF_SIZE (mapic.h)
*                       MODFILE		        (mapi.h)
*                       NC                      (local_nc.h)
*			NULLstr		        (mapic.h)
*                       MAPIERR                 (mapic.h)
*			DFACC_CREATE	        (hdf.h)
*                       DFACC_RDWR              (hdf.h)
*                       DFACC_READ              (hdf.h)
*			SDstart		        (libdf.a)
*                       NC_check_id             (local_nc.h)
*			Hishdf			(hproto.h)
*                       
* !Revision History:
*		Qi Huang	1996/07/09
*		Version 2.1.
*		Ring super structure and other changes make this 
*		version	much faster.
*
* $Log: openMODISfile.c,v $
* Revision 6.2  2010/07/13 14:10:58  kuyper
* Made argument conversions explicit.
* Corrected type and value of size arguments to sprintf() calls.
*
* Revision 6.1  2010/04/20 20:19:42  kuyper
* Changed to be compatible with versions of HDF > 4.2r1: replaced extraction
*   of H-interface file handle from SDS interface ID, with a seperate call
*   to Hopen. This necessitates a call to Vstart(), as well.
* Rearranged logic to ensure more reliable error detection.
*
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.3  1996/08/09  19:30:38  qhuang
 * Cast some arguments in HDF routine calls and some bitwise operation.
 *
 * Revision 1.2  1996/07/16  18:38:35  qhuang
 * Added some comments in Revision history in prolog.
 *
 * Revision 1.1  1996/07/16  18:38:13  qhuang
 * Initial revision
 *
 * Revision 1.7  1996/06/03  19:36:05  qhuang
 * Added Hishdf in external reference section in prolog.
 *
 * Revision 1.6  1996/05/23  17:40:23  qhuang
 * (Note: This information replaces 'minor changes' in revision 1.3.) Removed
 * one of the 'Revision 1.1' sections in the Revision History section, added
 * comments for the first two character declarations, stated which parameter
 * was empty in the first SMF message, in the code initializing internal HDF
 * structures, correct SMF message to state '...unable to allocate memory for
 * a MODIS file', deleted the word 'name'.
 *
 * Revision 1.5  1996/05/17  16:00:28  qhuang
 * Used Hishdf function call to check file's existance instead of using
 * fopen and fclose.
 *
 * Revision 1.5  1996/05/17  16:00:28  qhuang
 * Used Hishdf function call to check file's existance instead of using
 * fopen and fclose.
 *
 * Revision 1.4  1996/04/02  15:09:50  qhuang
 * Added #include "local_nc.h"
 *
 * Revision      1995/12/19  Xia W. Chang
 * Version 2 design. Use the hdf_id associated with sd_id instead of Hopen.
 * Modification of error message.
 *
 * Revision 1.3  1995/11/07  20:14:37  qhuang
 * minor changes
 *
 * Revision 1.2  1995/10/31  19:53:41  qhuang
 * Added capability to pass status messages to log files.
 *
 * Revision 1.1  1995/10/31  19:51:30  qhuang
 * Initial revision
 *
*
*       1.    Revision 01.00  1995/04/13   2200 EDT
*	      Jenny Glenn/GSC  (jglenn@modis-xl.gsfc.nasa.gov)
*
*	      Create the initial template for MODIS API utility subroutine
*	      openMODISfile.
*
*	2.    Dave Lorenzi/GSC 		28 Apr 95
*	      Original Development / Testing
*
*       3.    Gary Fu/GSC               5 Oct 95
*             Replace FEXIST macro with fopen call.
*
* !Team-unique Header:
*              This software is developed by the MODIS Science Data Support
*              Team for the National Aeronautics and Space Administration,
*              Goddard Space Flight Center, under contract NAS5-32373.
*
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
*
*
!END*********************************************************************
*/
{
   char  buff[PGS_SMF_MAX_MSGBUF_SIZE];         /* buffer to hold the error */
                                                /* warning message      */
   char *funcname="openMODISfile";              /* name of this routine */
   MODFILE *mfile = NULL;			/* MODISfile structure	*/
   MODFILE *retval = NULL;
   int     MAXILEN=100;

   if ( NULLstr(filename))               	/* Check input parameters*/
      {
      sprintf(buff, "ERROR: openMODISfile unable to access a file without"
	            " a filename input.\n");
      MAPIERR(buff,funcname);
      return(NULL);
      }
   if ( NULLstr(access))               
      {
      sprintf(buff, "ERROR: openMODISfile unable to open file \"%.*s\" " 
	            " without access mode input\n", MAXILEN, filename);
      MAPIERR(buff,funcname);
      return(NULL);
      }
							
   /* allocate MODFILE structure to output pointer  */
   if ((mfile= (MODFILE *) malloc(sizeof(MODFILE))) == NULL) 
      {
      sprintf(buff, "ERROR: openMODISfile unable to allocate memory for a"
	  " MODIS file structure for file \"%.*s\" .\n", MAXILEN, filename);
      MAPIERR(buff,funcname);
      return(NULL);
      }
    retval = mfile;
 
    /* Make assignment to structure's elements:
      Set dinfo component in the MODFILE structure to NULL */
    mfile->dinfo = NULL;

    /*
    **  Check to see if the file already exists, then set the file access
    **  mode according to the user input.  We must use two different file 
    **  access modes, one each for Hopen and SDstart.
     */

    /* Call Hishdf to check file's existance. */
    if ( (intn)Hishdf(filename) == TRUE )	/* file exists */
    {
	switch(*access)
	{
	    case 'r':
		mfile->access= DFACC_READ;		/* Read only old file */
		break;
	    case 'w':
		mfile->access= DFACC_CREATE;	/* Truncates old file */
		break;
	    case 'a':
		mfile->access= DFACC_RDWR;		/* Append to old file */
		break;
	    default:				/* Unkown access type */
		retval = NULL;
		sprintf(buff,"ERROR: openMODISfile unable to recognize access"
		    " type \"%.*s\" to open file \"%.*s\" .\n", 
		    MAXILEN, access, MAXILEN, filename);
		MAPIERR(buff,funcname);
	 }
    }
    else 				/* Get access mode for new file */
    {
	switch(*access)
	{
	    case 'w':
		mfile->access= DFACC_CREATE;	/* Open new file "r+" */
		break;
	    case 'a':
		mfile->access= DFACC_CREATE;	/* Open new file "r+" */
		break;
	    default:
		retval = NULL;
		if ( strcmp(access,"r") == 0 )
		{
		    sprintf(buff, "ERROR: openMODISfile unable to find "
			"file \"%.*s\" .\n", MAXILEN, filename);
		    MAPIERR(buff,funcname);
		}
		else
		{
		    sprintf(buff,"ERROR: openMODISfile unable to recognize "
			    "access type \"%.*s\" to open file \"%.*s\" .\n",
			    MAXILEN, access, MAXILEN, filename);
		    MAPIERR(buff,funcname);
		}
		break;
	}
    } 	/* Is an HDF file. */

    if (retval != NULL )
    {
	mfile->hdf_id = Hopen(filename, (int32)mfile->access, 0);
	if(mfile->hdf_id == FAIL)
	{
	    retval = NULL;
	    sprintf(buff, "ERROR: openMODISfile fails to get hdf_id while "
		"opening file \"%.*s\". ", MAXILEN, filename);
	    MAPIERR(buff, funcname);
	}
	else
	{
	    if(Vstart((HFILEID)mfile->hdf_id) != SUCCEED)
	    {

		retval = NULL;
		sprintf(buff, "ERROR: Vstart(\"%.*s\") failed",
		    MAXILEN, filename);
		MAPIERR(buff, funcname);
	    }
	    else
	    {
		int32 access = (int32)mfile->access;
		if(access == DFACC_CREATE)
		    access = DFACC_WRITE;

		mfile->sd_id= (long int)SDstart(filename, access);
		if (mfile->sd_id==FAIL)
		{
		    retval = NULL;
		    sprintf(buff, "ERROR: SDstart(\"%.*s\") failed.\n",
			MAXILEN, filename);
		    MAPIERR(buff,funcname);
		}
		else
		{
		    mfile->filename=(char *)malloc(strlen(filename)+1);
		    if (mfile->filename == NULL)
		    {
			 retval = NULL;
			 sprintf(buff, "ERROR: openMODISfile unable to "
			     "allocate memory for a MODIS file \"%.*s\".\n",
			     MAXILEN, filename);
			 MAPIERR(buff,funcname);
    
		    }
		    else
		       strcpy(mfile->filename,filename);

		    if(retval==NULL)
		    {
			free(mfile->filename);
			if(SDend((int32)mfile->sd_id) != SUCCEED)
			{
			    sprintf(buff,
				"ERROR: SDend(\"%.*s\") failed.\n",
				(int)(sizeof buff -
				sizeof "ERROR: SDend(\"\")"), filename);
			    MAPIERR(buff,funcname);
			}
		    }
		}

		if(retval==NULL && Vend((HFILEID)mfile->hdf_id) != SUCCEED)
		{
		    sprintf(buff, "ERROR: Vend(\"%.*s\") failed.\n",
			(int)(sizeof buff - sizeof "ERROR: Vend(\"\")"),
			filename);
		    MAPIERR(buff,funcname);
		}
	    }	/* Vstart() successful */
	    if(retval==NULL && Hclose((int32)mfile->hdf_id) != SUCCEED)
	    {
		sprintf(buff, "ERROR: Hclose(\"%.*s\") failed.\n",
		    (int)(sizeof buff - sizeof "ERROR: Hclose(\"\")"),
		    filename);
		MAPIERR(buff,funcname);
	    }
	}	/* Hopen() successful */
    }	/* valid access mode */
    if(retval==NULL)
	free(mfile);
    return retval ;
}
