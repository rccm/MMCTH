#include "mapic.h"  

int releaseMAPIfilehandle(MODFILE **file)
/*
!C**********************************************************************
* 
*!Description:	Function releaseMAPIfilehandle releases MODFILE
*		structure created by routine createMAPIfilehandle. 
*		This routine also automatically closes all the 
*		objects opend with M-API routines.  The MAPI file
*		handle created by calling createMODISfilehandle
*		should not be passed to M-API routine closeMODISfile
*		or completeMODISfile.
*
*               Users may open an HDF-EOS file using HDF-EOS opening routine;
*               call this routine to create the MAPI file handle; then use M-API
*               routines to access data object(s) in the opened file; after fini
shing
*               object access call releaseMAPIfilehandle to release this file ha
ndle
*               created by calling this routine (releaseMAPIfilehandle also 
*               automatically closes all the objects opened with M-API routines)
* 
*!Input parameters:
*  **file	IN:	Address of pointer to MODFILE structure.
*
*!Output parameters:	None
*
* Returns:		MAPIOK if successful, or MFAIL if an error occurs.
*
* External references:
*                       MODFILE                         (mapic.h)
*                       PGS_SMF_MAX_MSGBUF_SIZE         (mapic.h)
*			MAPIOK				(mapi.h)
*			endMODISobjaccess		(mapi.h)
*			MFAIL				(mapi.h)
*			MODIS_ALL_TYPES			(mapi.h)
*
*!Revision history:
*$Log: releaseMAPIfilehandle.c,v $
*Revision 1.1  1998/02/06 22:26:06  fshaw
*Initial revision
*
 * Revision 1.2  1997/05/28  20:11:56  qhuang
 * Version 2.3.
 *
 * Revision 1.1  1997/05/28  20:11:32  qhuang
 * Initial revision
 *
*		Qi Huang	1997/05/28
*		Version 2.3
*		Original development and testing
* 
*!Team-unique header:
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
{
  char          buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char          *funcname = "releaseMAPIfilehandle";
  int		status;

  status = MAPIOK;

  if (file && *file)
  {
    if (endMODISobjaccess(*file,NULL,NULL,MODIS_ALL_TYPES) == MFAIL)
    {
      sprintf(buff,"ERROR: releaseMAPIfilehandle detected MFAIL from\n"
			"\t endMODISobjaccess while trying to close\n"
			"\t all opened objects for file.\n");
      MAPIERR(buff,funcname);
      return(MFAIL);
    }

    if ((*file)->filename)
      free((*file)->filename);

    free(*file);
    (*file) = NULL;
  }

  else
  {
    sprintf(buff,"ERROR: releaseMAPIfilehandle can not release a NULL\n"
			"\t file handle.\n");
    MAPIERR(buff,funcname);
    status = MFAIL;
  }

  return(status);
}

