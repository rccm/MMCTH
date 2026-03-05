#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"

int completeMODISfile (MODFILE **file,
                       PGSt_MET_all_handles mdHandles,
                       ECSattr_names_for_all_handles HDFattrNames,
                       long int NumHandles)
/*
!C**********************************************************************
*
*!Description: Subroutine completeMODISfile is part of a larger software
*              system called the MODIS Applications Programming Interface
*              (API) Utility, abbreviated M-API.  The M-API Utility
*              consists of subroutines which allow MODIS Science
*              Team-supplied software to read in Level 1B radiance
*              bands and write out output products and metadata to
*              HDF files.  The functionality of the M-API is defined
*              in the MODIS API User's Guide, Version 1, dated 4/3/95.
*
*	           completeMODISfile terminates the access of MODIS API
*              routines to a MODIS HDF file opened using openMODISfile
*              (OPMFIL).  In addition to closing the file, the MODIS
*              file's standard header information is inserted. A
*              pre-existing MODIS-HDF file should be closed by
*              closeMODISfile(CLMFIL) or some of its header information
*              will be over-written.  completeMODISfile may fail to
*              close the file if an error occurs.
*
*              See Chapter 8, Accessing Metadata, for a complete list of
*              metadata completeMODISfile writes to the MODIS-HDF file
*              before closing it.
*
*!Input parameters:
*  MODFILE: MODIS file structure
*           Pointer to MODFILE structure's address
*           used to reference a file in all MODIS
*           API routines.
*
*  PGSt_MET_all_handles: A character array with size of
*                        [PGSd_MET_NUM_OF_GROUPS][PGSd_MET_GROUP_NAME_L],
*                        where PGSd_MET_NUM_OF_GROUPS is 20 and
*                        PGSd_MET_GROUP_NAME_L is 50.  This array
*                        is typedef-ed as PGSt_MET_all_handles.
*                        Each row in the array stores a handles to an
*                        internal ODL tree structure which will be written
*                        out a an ECS PVL attribute.  Each handles,
*                        which is a string, should be less than
*                        50 characters and occupy, one row in the
*                        array.   Therefore, the maxmum number of
*                        handles should be 20.
*
*  ECSattr_names_for_all_handles: A character array with size of
*                                [PGSd_MET_NUM_OF_GROUPS][MAX_ECS_NAME_L],
*                                where PGSd_MET_NUM_OF_GROUPS is 20 and
*                                MAX_ECS_NAME_L is 50.  This array is
*                                typedef-ed as ECSattr_names_for all_handles.
*                                Each row in this array is a character string
*                                used as a global attribute name for storing
*                                an ECS PVL text block which has a handle in
*                                the corresponding row in mdHandles array.
*                                Each name, which is a string, should be less
*                                that MAX_ECS_NAME_L characters and occupies
*                                one row in the array.
*
*  long int: Specifie the number of actual handles contained
*            in mdHandles.  This may be set from 0 to
*            PGSd_MET_NUM_OF_GROUPS.
*
*!Output parameters:	NONE
* 
* Returns: MAPIOK if successful, MFAIL if an error occurs.
*
* Externally defined:	MODFILE		          (mapi.h)
*                       closeMODISfile()          (mapi.h)              
*                       DFACC_READ                (hdf.h)
*                       DFACC_RDWR                (hdf.h)
*                       PGS_MET_write             (mapic.h)
*                       PGSd_MET_NUM_OF_GROUPS    (PGS_MET.h) 
*                       PGS_SMF_MAX_MSGBUF_SIZE   (PGS_SMF.h) 
*                       PGS_S_SUCCESS             (PGS_SMF.h)
*                       MAPIOK                    (mapi.h)
*                       MFAIL                     (mapi.h)
*                       MAPIERR                   (mapic.h)
*                       MAPIWARN                  (mapic.h)
*
*!Revision History:
* $Log: completeMODISfile.c,v $
* Revision 6.1  2010/07/13 14:03:08  kuyper
* Removed const from multi-dimensional function parameters, which cannot work
*   the way it was intended to.
*
* Revision 5.1  2005/04/04 18:17:57  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.9  1997/03/13  18:50:29  fshaw
 * *** empty log message ***
 *
 * Revision 1.8  1997/03/12  21:10:26  fshaw
 * change error to warning
 *
 * Revision 1.7  1996/10/11  20:04:09  fshaw
 * *** empty log message ***
 *
 * Revision 1.1  1996/09/26  17:31:58  qhuang
 * Initial revision
 *
 * Revision 1.2  1996/06/06  20:14:43  qhuang
 * ..
 *
 * Revision 1.4  1996/05/22  14:14:42  fshaw
 * corrected module to match the PDL
 *
 * Revision 1.3  1996/04/24  19:37:41  fshaw
 * no change
 *
 * Revision 1.2  1996/04/23  17:29:16  fshaw
 * compiled successfully with make 4/23/96 FJS
 *
 * Revision 1.1  1996/04/18  21:04:28  fshaw
 * Initial revision
 *
*
*
*!Team-unique Header:
*	      Portions developed at the National Center for Supercomputing
*	      Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!References and Credits:
*	      This software is developed by the MODIS Science Data Support
*	      Team for the National Aeronautics and Space Administration,
*	      Goddard Space Flight Center, under contract NAS5-32373.
*
*!Design Notes:
*
*!END*********************************************************************
*/

{
    char buff[PGS_SMF_MAX_MSGBUF_SIZE];
    char *funcname = "completeMODISfile";
    int status = MAPIOK;		/* subroutine return status  */
    int i;
    int MAXILEN = 100;

/* Check for NULL input file or invalid *file          */


    if ((!file) || !(*file) || !((*file)->access)){                        
	    sprintf(buff,"ERROR: completeMODISfile cannot close a null file.");
	    MAPIERR(buff, funcname);
	    return(MFAIL);
    }
    if ((*file)->access == DFACC_READ){
	    sprintf(buff,"Warning: completeMODISfile unable to revise header"
                "of %.*s file open for read-only.\n", MAXILEN, (*file)->filename);
	    MAPIWARN(buff, funcname);
    }
    else{
	    if ((*file)->access == DFACC_RDWR){
	        sprintf(buff,"Warning: completeMODISfile revised header"
                    "data of pre-existing %.*s file.\n", MAXILEN, (*file)->filename);
	        MAPIWARN(buff,funcname);
	    }
	    if ((NumHandles < 0) || (NumHandles >= PGSd_MET_NUM_OF_GROUPS)){
	        sprintf(buff,"ERROR: completeMODISfile found the number of"
                    "handles in the input argument out of range.\n");
	        MAPIERR(buff,funcname);
	        status = MFAIL;
	    }
	    if ( status == MAPIOK){
	        for ( i = 1; i <= NumHandles; i++){
	    	    if (PGS_MET_Write((char *)mdHandles[i], 
                        (char *)HDFattrNames[i], (PGSt_integer)(*file)->sd_id)!= PGS_S_SUCCESS){
		           sprintf(buff,"ERROR: completeMODISfile detected failure"
                   "from SDP function PGS_MET_Write to store the %.*s files's"
                   "ECS PVL metadata %.*s.", MAXILEN, (*file)->filename, MAXILEN, 
                    HDFattrNames[i]);
		           MAPIERR(buff,funcname);
		           status = MFAIL;
		        }
            }
        }
    (*file)->access = DFACC_RDWR;
    }
    if (closeMODISfile(file) == MFAIL){
	sprintf(buff,"ERROR: completeMODISfile detected FAIL from "
               "MAPI function closeMODISfile while trying to"
               " completeMODISfile %.*s\n", MAXILEN, (*file)->filename);
        MAPIERR(buff,funcname);
        status = MFAIL;
    }
    return(status);
}
