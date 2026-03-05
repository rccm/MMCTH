#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mapic.h"
#define  MAXATTRLEN	120

int  putMODISdiminfo(MODFILE *file, char const *arrayname, 
                     char const *groupname, long int dimension, 
                     char const *attribute, char const *data_type,
                     long int n_elements, void const *value)

/*
!C****************************************************************************
*
*!Description: Function putMODISdiminfo is part of a larger software system
*              called the MODIS Applications Programming Interface(API)
*              Utility, abbreviated M-API.  The M-API Utility consists of
*              subroutines which allow MODIS Science Team-supplied software
*              to read and write data and metadata from/to HDF files.
*              The functionality of the M-API is defined in the MODIS
*              Application Program Interface(API) Specification.
*
*              putMODISdiminfo stores an attribute = value(s) metadata pair
*              to the indicated dimension of a MODIS array.  If the attribute
*              already exists, the value(s) will be updated.
*
*!Input Parameters:
* file         Address of MODFILE structurethat is used to reference the 
*              MODIS-HDF file receiving the metadata.
* arrayname    ASCII string name of the array.  Provided macros for accepted
*              MODIS HDF file array names are listed in Appendix A, M-API-
*              Supplied Constants.
* groupname    ASCII string name of the data group containing the array
*              structure to which the metadata is attached.  If set to NULL,
*              the entire file will be searched for the array structure named
*              arrayname.
* attribute    Name to assign the attribute.  Provided macros for accepted
*              MODIS file metadata names are listed in Appendix A, MODIS 
*              API Supplied Constants.
* dimension    The dimension to which the metadata is attached (0-based).
* data_type    Data Type of the value.
*
*              Permitted C data types:
*                "int8"
*                "uint8"
*                "int16"
*                "uint16"
*                "int32"
*                "uint32"
*                "int64"
*                "float32"
*                "float64"
*                "char *"
*
* n_elements   Number of metadata values in value.
* value        Address of the data to store in the attribute.  If the attribute
*              already exists, the value will be updated.  Values should 
*              conform to the data types, formats and/or those values
*              enumerated for the attribute in Appendix A, MODIS API - Supplied
*              Constants.
*
*!Output Parameters:NONE
*
* External references:
*              MODFILE(mapi.h)
*	       PGS_SMF_MAX_MSGBUF_SIZE(PGS_SMF.h)
*	       DATAID(mapi.h)
*              MAPIOK(mapi.h)
*              MFAIL(mapi.h)
*              MAX_REC(mapic.h)
*              NULLstr(mapic.h)
*              datatype_to_DFNT(mapic.h)
*              DFACC_READ(hdf.h)
*              DFKNTsize(hproto.h)
*	       getMODISarrayid(mapic.h)
*	       SDSINFO(mapic.h)
*              SDgetdimid(mfhdf.h)
*              SDsetattr(mfhdf.h)
*	       MAX_NC_NAME(netcdf.h)
*	       NULLMODFIL(mapic.h)
*	       DATATYPELENMAX(mapic.h)
*	       VOIDP(hdfi.h)
*	       MAPIERR(mapic.h)
*
* Returns:     MAPIOK if successful
*              MFAIL if an error occurs.
*
*!Revision History:
*		Qi Huang	1996/08/21
*		Version 2.1
*		Ring super structure and other changes make this
*		version much faster.
*
* $Log: putMODISdiminfo.c,v $
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.5  1996/08/30  19:40:39  qhuang
 * Version 2.1
 *
 * Revision 1.4  1996/03/07  21:43:21  qhuang
 * Version 2 design, group handling, separated input checking, new messages,
 * removed the specifical case for TXT type which uses strlen to obtain the
 * length, removed the warning message about overwriting the existed attribute.
 *
 * Revision 1.3  1995/11/07  20:42:00  qhuang
 * minor changes
 *
 * Revision 1.2  1995/10/31  22:30:22  qhuang
 * Added capability to pass status messages to log files.
 *
 * Revision 1.1  1995/10/31  22:29:17  qhuang
 * Initial revision
 *
*
*              Frank Chen frank@modis-xl.gsfc.nasa.gov  Sept. 23, 1995
*              initial implementation
*
*!Team-unique Header:
*              This software is developed by the MODIS Science Data Support
*              Team for the National Aeronautics and Space Administration,
*              Goddard Space Flight Center, under contract NAS5-32373.
*
*!References and Credits:
*
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes:
*
!END***************************************************************************
*/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning message */
  char *funcname="putMODISdiminfo";     /* name of this routine */

  DATAID	*did;
  int	status_code;              /* return value for routine.MAPIOK = successful
                                                            MFAIL = fail */
  long int size;          /* total number of byte to write */
  int32 dtype;                  /* HDF attribute data type */
  int32 dim_id;                 /* SDS dimension id */

  /* initialize local variables */
  size = 0;
  dtype = 0;
  dim_id = 0;

  /* initialize return status code */
  status_code = MFAIL;

  /* Input checks: */
  if ( NULLstr(attribute) )
  {
    sprintf(buff,"ERROR: putMODISdiminfo unable to write an dimension\n"
			"\t attribute without an attribute name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLstr(data_type) )
  {
    sprintf(buff,"ERROR: putMODISdiminfo unable to write the %.*s\n"
			"\t dimension attribute without data type\n"
			"\t information\n",MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( value == NULL )
  {
    sprintf(buff,"ERROR: putMODISdiminfo unable to write the %.*s\n"
			"\t dimension attribute without the value buffer.\n",
			MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLstr(arrayname) )
  {
    sprintf(buff,"ERROR: putMODISdiminfo unable to write the %.*s\n"
			"\t dimension attribute without the name of the\n"
			"\t array it is associated with.\n",
			MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) )
  {
    sprintf(buff,"ERROR: putMODISdiminfo unable to write the %.*s\n"
			"\t dimension attribute with an invalid MODIS file\n"
			"\t structure input.\n",MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( n_elements < 1 )
  {
    sprintf(buff,"ERROR: putMODISdiminfo unable to write %ld\n"
			"\t %.*s dimension attribute values.\n",
		 	n_elements,MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);	
    return(MFAIL);
  }

  if (file->access == DFACC_READ)
  {
    sprintf(buff,"ERROR: putMODISdiminfo unable to write the %.*s\n"
			"\t dimension attribute in a file opened for read\n"
			"\t only.\n",MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ((dtype = datatype_to_DFNT(data_type)) == MFAIL) {
    sprintf(buff, "ERROR: putMODISdiminfo unable to write the %.*s\n"
			"\t dimension attribute of data type %.*s.\n",
                    	MAX_NC_NAME,attribute,DATATYPELENMAX, data_type);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( (size =  n_elements * DFKNTsize(dtype) ) > MAX_REC )
  {
    sprintf(buff,"ERROR: putMODISdiminfo unable to write the %.*s\n"
			"\t dimension attribute with a %ld byte value.\n",
			MAX_NC_NAME,attribute,size);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  did = getMODISarrayid(file,arrayname,groupname);

  if ( did == NULL )
  {
    sprintf(buff,"ERROR: putMODISdiminfo detected errors from\n"
			"\t getMODISarrayid attempting to write the\n"
			"\t %.*s dimension attribute.\n",
			MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
  }
  /* else check value of dimension is valid */
  else if ((dimension < 0) || ( dimension > (((SDSINFO *)did->info)->rank - 1))) 
  {
    sprintf(buff, "ERROR: putMODISdiminfo unable to write the %.*s\n"
			"\t attribute to non-existing dimension %ld of\n"
                    	"the %.*s array.\n",MAXATTRLEN,attribute,
			dimension,MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
  }

  /* else get the dimension id */
  else if ((dim_id = SDgetdimid((int32)did->id,(int32)dimension)) == FAIL) 
  {
    sprintf(buff, "ERROR: putMODISdiminfo detected FAIL from HDF\n"
			"\t procedure SDgetdimid attempting to write the\n"
			"\t %.*s dimension attribute.\n",
                    	MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
  }

  else
  {
    if (SDsetattr(dim_id,attribute,dtype,(int32)n_elements,(VOIDP)value) == FAIL) 
    {
      sprintf(buff,"ERROR: putMODISdiminfo detected FAIL from HDF\n"
			"\t procedure SDsetattr attempting to write\n"
			"\t the %.*s dimension attribute.\n",
                    	MAX_NC_NAME,attribute);
      MAPIERR(buff,funcname);
    }

    else 
      status_code = MAPIOK;

  } /* end of else */

  return(status_code);
}
