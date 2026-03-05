#include	"mapic.h" 
#include	"HdfEosDef.h"

MODFILE	*createMAPIfilehandle(int32 fid)
/* 
!C**********************************************************************
* 
*!Description:	Function createMAPIfilehandle provide a MAPI file handle for a 
*		previously opened HDF-EOS file (swath,grid or point). A routine
*	           called releaseMAPIfilehandle must be called for each call to this
*	           routine or memory leak is likely happen. NOTE: This file handle
*		should not be passed to MAPI routine closeMODISfile or 
*		completeMODISfile.  The file must be opened and closed using 
*		one API before it can be opened using a different API.
*
*		Users may open an HDF-EOS file using HDF-EOS opening routine;
*		call this routine to create the MAPI file handle; then use M-API
*		routines to access data object(s) in the opened file; after finishing
*		object access call releaseMAPIfilehandle to release this file handle
*		created by calling this routine (releaseMAPIfilehandle also 
*		automatically closes all the objects opened with M-API routines);
*		finally call HDF-EOS closing routine to close the file.
*
*
*!Input parameters:
*
*    fid		IN:	HDF-EOS file ID.
*
*!Output parameters:		None
*
*Returns:	Address of MODFILE structure that is used to reference the file in
*		all other MODIS API C routines, or NULL if an error occurs.
*
* External references:
*			MODFILE                         (mapic.h)
*                       PGS_SMF_MAX_MSGBUF_SIZE         (mapic.h)
*			DATAINFO			(mapi.h)
*			EHidinfo			(HdfEosDef.h)
*			Hfidinquire			(hproto.h)
*
*!Revision History:
*$Log: createMAPIfilehandle.c,v $
*Revision 1.1  1998/02/06 22:26:06  fshaw
*Initial revision
*
 * Revision 1.1  1997/09/08  14:34:03  fshaw
 * Initial revision
 *
 * Revision 1.3  1997/08/08  15:01:43  qhuang
 * Corrected a typo.
 *
 * Revision 1.2  1997/05/28  20:08:29  qhuang
 * Version 2.3.
 *
 * Revision 1.1  1997/05/28  20:07:51  qhuang
 * Initial revision
 *
*	Qi Huang	1997/05/28
*	version 2.3
*	Original development and testing
* 
*!Team-unique header:
* This software is developed by the MODIS Science Data SupportTeam for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
*!References and Credits
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes
*
!END********************************************************************
*/
{
  char		buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char		*funcname = "createMAPIfilehandle";
  MODFILE	*mfile = NULL;
  int32		sd_id;
  int32		hdf_id;
  char		*filename;
  intn		access, attach;


  if (EHidinfo(fid,&hdf_id,&sd_id) == FAIL)
  {
    sprintf(buff,"ERROR: createMAPIfilehandle detected FAIL from \n"
			"\t HDF-EOS procedure EHidinfo.\n");
    MAPIERR(buff,funcname);
    return(NULL);
  }

  if (Hfidinquire(hdf_id,&filename,&access,&attach) == FAIL)
  {
    sprintf(buff,"ERROR: createMAPIfilehandle detected FAIL from\n"
			"\t HDF procedure Hfidinquire while retrieving\n"
			"\t filename and access mode.\n");
    MAPIERR(buff,funcname);
    return(NULL);
  }

  /* Allocate MODFILE structure to pointer mfile */
  if ( (mfile = (MODFILE *)malloc(sizeof(MODFILE))) == NULL)
  {
    sprintf(buff,"ERROR: createMAPIfilehandle is unable to allocate\n"
			"\t memory for a MODIS file structure for file %s.\n",
			filename);
    MAPIERR(buff,funcname);
    return(NULL);
  }

  /* Make assignments to structure's elements: */
  mfile->dinfo = NULL;
  mfile->access = access;

  /* Allocate memory for mfile->filename */
  if ( (mfile->filename = (char *)malloc(strlen(filename)+1)) == NULL)
  {
    free(mfile);
    mfile = NULL;
    sprintf(buff,"ERROR: createMAPIfilehandle is unable to allocate \n"
			"\t memory for the MODIS-HDF file structure's \n"
			"\t filename %s.\n",filename);
    MAPIERR(buff,funcname);
  }

  else
    strcpy(mfile->filename,filename);

  if (mfile)
  {
    mfile->sd_id = sd_id;
    mfile->hdf_id = hdf_id;
  }

  return(mfile);
}

