#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include "mapic.h"
#define  MAXDIMNAMELEN	60

int  putMODISdimname(MODFILE *file, char const *arrayname, 
                     char const *groupname, long int dimension, 
                     char const *dimname)
/*
!C****************************************************************************
*
*!Purpose:	Assigns a name to a specific dimension of an array structure..
* 
*!Description: 	Function putMODISdimname is part of a larger software system 
*		called the MODIS Applications Programming Interface (API) 
*		Utility, abbreviated M-API.  The M-API Utility consists of 
*		subroutines which allow MODIS Science Team-supplied software 
*		to read  and write data and metadata from/to HDF files. The 
*		functionality of the M-API is defined in the MODIS Application 
*		Program Interface (API) Specification.
*		
*		putMODISdimname (PMDNAM) associates an HDF dimension name 
*		with a specified SDS array structure dimension. The SDS array 
*		must be created (using createMODISarray) before it is possible to 
*		name any of its dimensions. This routine does not create a 
*		"long_name" dimension attribute. putMODISdimname can produce 
*		such a dimension label, however.
*		
*		putMODISdimname does more than apply an appellation to a 
*		dimension. An HDF dimension name is an independent data 
*		object. It may be shared by several array structure dimensions, 
*		but they all must be of the same size. Any dimension attribute 
*		that is associated with any one of these dimensions is immediately 
*		associated with all the dimensions having that name. Likewise, 
*		updating a dimension attribute for one dimension updates it for 
*		all dimensions having the same name (they could only have one 
*		"long name" dimension shared between them).
*		
*		Naming an SDS dimension will also cause any dimension 
*		attributes currently associated with that dimension to be lost. 
*		Therefore it is most practical to name an array's dimensions, if 
*		necessary, immediately after the array structure's creation and 
*		before creating dimension attributes for it.
*		
*		The groupname string provides the facility to select an array 
*		structure placed in a particular HDF 'Vgroup' data group. 
*		Alternatively, the entire file will be searched for an array 
*		structure named arrayname if the argument groupname = NULL 
*		in C or grpnm is a blank string ( " ") in FORTRAN.
* 
* !Input Parameters:
* file		IN: 	Address of MODFILE structure that is used to 
* 		reference a MODIS HDF file receiving the dimension 
* 		name.
* arrayname	IN:	ASCII string name of the target array 
* 		structure.
* groupname	IN:	ASCII string name of the data group 
*		containing the target array structure. If set to NULL 
*		the entire file will be searched for the array 
*		structure named arrayname.
* dimension	IN:	The dimension number which is to be named 
*		(0-based). The 0 dimension of an HDF SDS array 
*		structure is associated with the least rapidly varying 
*		array index.
* dimname	IN:	ASCII string name to give to the dimension.
* 
* !Output Parameters:		NONE
* 
* Returns:	MAPIOK if successful, MFAIL if an error occurs.
* 
* External references:
*              MODFILE(mapi.h)
*	       PGS_SMF_MAX_MSGBUF_SIZE(mapic.h)
*              MAPIOK(mapi.h)
*              MFAIL(mapi.h)
*	       MAPIERR(mapic.h)
*              NULLstr(mapic.h)
*	       NULLMODFIL(mapic.h)
*	       MAX_NC_NAME(netcdf.h)
*              DFACC_READ(hdf.h)
*	       getMODISarrayid(mapic.h)
*	       SDdiminfo(mfhdf.h)
*	       MAPIWARN(mapic.h)
*	       SDSINFO(mapic.h)
*              SDgetdimid(mfhdf.h)
*              SDsetdimname(mfhdf.h)
*
* !Revision History:
* 		Qi Huang	1996/07/23
*		Version 2.1
*		Original development and testing
*		Ring super structure and other changes make this
*		version much faster.
*
* $Log: putMODISdimname.c,v $
* Revision 6.1  2010/07/13 19:27:32  kuyper
* Corrected format specifier and argument type to match.
*
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.5  1997/12/02  21:10:35  fshaw
 * removed extraneous commnets
 *
 * Revision 1.4  1997/12/01  21:51:53  fshaw
 * corrected defect in code replace l_dimname with NULL
 *
 * Revision 1.3  1996/08/09  14:57:22  qhuang
 * Casted some arguments in HDF routine calls, changed .s*s to .*s in a error
 * message.
 *
 * Revision 1.2  1996/08/08  16:45:02  qhuang
 * Added developer's name, etc. to Reference and Credits section in prolog,
 * added comments to body of code near major blocks, removed some comment in
 * variable dimname declaration in prolog.
 *
 * Revision 1.1  1996/08/08  16:44:42  qhuang
 * Initial revision
 *
*
* !Team-unique Header:
*              This software is developed by the MODIS Science Data Support
*              Team for the National Aeronautics and Space Administration,
*              Goddard Space Flight Center, under contract NAS5-32373.
*
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!References and Credits
*
* 		Qi Huang	1996/07/23
*		Version 2.1
*		Original development and testing
*		Ring super structure and other changes make this
*		version much faster.
*
*!Design Notes
*
!END***************************************************************************
*/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning message */
  char *funcname="putMODISdimname";     /* name of this routine */
  int status_code;              /* return value for routine.MAPIOK = successful
                                                            MFAIL = fail */
  int32 sds_id;                 /* SDS id */
  int32 dim_id;                 /* SDS dimension id */
  DATAID *did;
  SDSINFO *sinfo;
  int32 count;
  int32 number_type;
  int32 nattrs;

  status_code = MFAIL;		/* error */

  /* Input checks: */
  if ( NULLstr(dimname) )
  {
    sprintf(buff,"ERROR: putMODISdimname unable to name a dimension\n"
			"\t without a dimension name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLstr(arrayname) )
  {
    sprintf(buff,"ERROR: putMODISdimname unable to name a dimension\n"
			"\t %.*s without the name of array the dimension\n"
			"\t is associated with.\n",MAX_NC_NAME,dimname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) )
  {
    sprintf(buff,"ERROR: putMODISdimname unable to name a dimension\n"
			"\t %.*s in the %.*s array with an invalid\n"
			"\t MODIS file structure input.\n",
			MAXDIMNAMELEN,dimname,MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( file->access == DFACC_READ )	/* read only */
  {
    sprintf(buff,"ERROR: putMODISdimname unable to name a dimension\n"
			"\t %.*s in a file opened for read only.\n",
			MAX_NC_NAME,dimname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* Get SDS array id */
  did = getMODISarrayid(file,arrayname,groupname);

  if ( did == NULL )
  {
    sprintf(buff,"ERROR: putMODISdimname detected MFAIL from M-API\n"
			"\t internal function getMODISarrayid while\n"
			"\t attempting to name a dimension %.*s in the\n"
			"\t %.*s array.\n",MAXDIMNAMELEN,dimname,
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }
  
  sinfo = (SDSINFO *)did->info;
  sds_id = (int32)did->id;

  /* Check dimension */
  if ( (dimension < 0) || (dimension > (sinfo->rank - 1)) )
  {
    sprintf(buff,"ERROR: putMODISdimname unable to name a non-existing\n"
			"\t dimension %.*s dimname.\n",MAX_NC_NAME,dimname);
    MAPIERR(buff,funcname);
  }
  /* Get dimension id */
  else if ( (dim_id = SDgetdimid(sds_id,(intn)dimension)) == FAIL) 
  {
    sprintf(buff, "ERROR: putMODISdimname detected FAIL from HDF procedure\n"
			"\t SDgetdimid attempting to name a dimension %.*s.\n",
                        MAX_NC_NAME,dimname);
    MAPIERR(buff,funcname);
  }
  /* Retrives number of attributes */
  else if (SDdiminfo(dim_id, NULL, &count, &number_type, &nattrs) == FAIL)
  {
    sprintf(buff,"ERROR: putMODISdimname detected FAIL from HDF\n"
			"\t procedure SDdiminfo attempting to name a\n"
			"\t dimension %.*s.\n",MAX_NC_NAME,dimname);
    MAPIERR(buff,funcname);
  }
  else
  /* If attribute(s) already exsit(s), give warning */
  {
    if ( nattrs > 0 )
    {
      sprintf(buff,"WARNING: putMODISdimname detected %ld\n"
			"\t attributes currently attached to the\n"
			"\t dimension. Naming the %ld dimension\n"
			"\t of the %.*s array %.*s will lose\n"
			"\t those attributes.\n", (long)nattrs,dimension,
			MAX_NC_NAME,arrayname,MAXDIMNAMELEN,dimname);
      MAPIWARN(buff,funcname);
    }

    /* Set dimension name */
    if (SDsetdimname(dim_id,dimname) == FAIL) 
    {
      sprintf(buff, "ERROR: putMODISdimname detected FAIL from HDF procedure\n"
			"\t SDsetdimname attempting to name a dimension %.*s.\n",
                        MAX_NC_NAME,dimname);
      MAPIERR(buff,funcname);
    }
    else
      status_code = MAPIOK;

  }/* end of else */

  return(status_code);
}
