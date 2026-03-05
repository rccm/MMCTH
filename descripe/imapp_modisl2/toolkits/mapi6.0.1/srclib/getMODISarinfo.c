#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mapic.h"
#define  MAXATTRNAMELEN  80

int  getMODISarinfo(MODFILE *file, char const *arrayname, char const *groupname, 
                   char const *attribute, char *data_type, long int *n_elements,
                   void *value)

/*
!C****************************************************************************
*
*!Purpose:     Reads the value(s) of a local attribute attached to an 
*              SDS (array) from a MODIS HDF file.
*
*!Description: Function getMODISarinfo is part of a larger software system
*              called the MODIS Applications Programming Interface(API)
*              Utility, abbreviated M-API.  The M-API Utility consists of
*              subroutines which allow MODIS Science Team-supplied software
*              to read and write data and metadata from/to HDF files.
*              The functionality of the M-API is defined in the MODIS
*              Application Program Interface(API) Specification.
*
*              getMODISarinfo retrieves the value(s) associated with an
*              attribute = value(s) pair attatched to an SDS (array) by 
*              giving the array name and attribute name.  If the attribute 
*              cannot be found, the routine will return MFAIL, and the 
*              passed variable unchanged.
*
*              The routine will also fail if the provided data_type is found
*              to be different from the attribute's data type or the n_elements
*              is found to be too small to contain the attribute's value(s).
*              getMODISarinfo replaces this input information with the actual
*              data type and number of elements contained in the attribute 
*              value(s) (in the case of character data, it is the length of 
*              the string, including the '\0' character).  These attribute's 
*              attribute may be used to properly retrieve the metadata value 
*              with a second call to the routine.  If a function failure 
*              occurs, n_elements will be set to zero.
*
*              A variable of the proper data type should be passed for the
*              value parameter.  The data type information required to properly
*              use either routine may be found in Appendix A, M-API-Supplied 
*              Constants, and Appendix C, MODIS Data Product File Definitions
*              of M-API User's Guide.  Appendix A has a listing for each M-API
*              provided metadata attribute that includes the data type, the
*              format, and/or specific values associated with it.
*
*!Input Parameters:
* file         IN:  Address of MODFILE structurethat is used to reference the 
*              MODIS HDF file containing the attribute.
* arrayname    IN:  ASCII string name of the array.  Provided macros for accepted
*              MODIS HDF file array names are listed in Appendix A, M-API-
*              Supplied Constants.
* groupname    IN:  ASCII string name of the data group containing the array
*              structure to which the attribute is attached.  If set to NULL,
*              the entire file will be searched for the array structure named
*              arrayname.
* attribute    IN:  ASCII string name of the attribute.  Provided macros for accepted
*              MODIS HDF file attribute names are listed in Appendix A, MODIS 
*              API Supplied Constants. (Note: Limit the length of attribute
*	       to 80 bytes)
* data_type    IN/OUT:  Address of data type of the value output.  Output replaces with
*              the data type of the retrieved attribute.  The memory size of this
*              arguments should be at least 8 bytes long.
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
* n_elements   IN/OUT:  Address of the number of elements the  value buffer can contain.
*              Output replaces with the number of elements required to
*              contain the attribute.  If a function failure occurs, the value
*              will be set to zero.
*
*!Output Parameters:
* value        OUT:  Address of value(s) associated with the attribute.
*
* Returns:     MAPIOK if successful, MFAIL if value cannot contain the 
*              retrieved attribute value, the data type is different, 
*              the attribute cannot be found, or an error occurs.
*
* External references:
*	       DATAID(mapi.h)
*              NULL(stdlib.h)
*              MODFILE(mapi.h)
*              MAPIOK(mapi.h)
*              MFAIL(mapi.h)
*              NULLstr(mapic.h)
*              DATATYPELENMAX (mapic.h)
*              DFNT_to_datatype(mapic.h)
*              DFNT_CHAR8(hdf.h)
*              MAX_NC_NAME(netcdf.h)
*	       getMODISarrayid(mapic.h)
*              SDfindattr(mfhdf.h)
*              SDattrinfo(mfhdf.h)
*              SDreadattr(mfhdf.h)
*	       PGS_SMF_MAX_MSGBUF_SIZE(PGS_SMF.h)
*	       NULLMODFIL(mapic.h)
*	       VOIDP(hdfi.h)
*	       MAPIERR(mapic.h)
*
*!Revision History:
*		Qi Huang	1996/07/26
*		Version 2.1
*		Ring super structure and other changes make
*		this version much faster.
*
* $Log: getMODISarinfo.c,v $
* Revision 5.1  2005/04/04 18:47:54  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.4  1996/08/12  19:47:50  qhuang
 * Version 2.1.
 *
 * Revision 1.3  1996/03/25  14:44:40  qhuang
 * cast long int(l_count) to assignment of *n_elements, added L to long int type
 * variables' assignment.
 *
 * Revision 1.2  1996/03/08  23:10:11  qhuang
 * Version 2, group handling, new messages, input checking, added the error
 * checking after DFNT_to_datatype.
 *
 * Revision 1.1  1996/03/08  23:07:45  qhuang
 * Initial revision
 *
 * Revision 1.1  1996/03/08  23:07:45  qhuang
 * Initial revision
 *
*
*              Frank Chen frank@modis-xl.gsfc.nasa.gov  Sept. 17, 1995
*              initial implementation
*
*!Team-unique Header:
*              This software is developed by the MODIS Science Data Support
*              Team for the National Aeronautics and Space Administration,
*              Goddard Space Flight Center, under contract NAS5-32373.
*
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
*
!END**************************************************************************
*/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char *funcname = "getMODISarinfo";  /* name of this module */

  DATAID *did;
  int status_code;                /* return value for routine.
                                     MAPIOK = successful
                                     MFAIL = fail */
  long int nelements_in;          /* initial input value of *n_elements */
  int32 l_count;                  /* number of data count from SDattrinfo() */
/*  int32 sds_index; */             /* SDS index number */
/*  int32 sds_id;    */             /* SDS id */
  int32 attr_index;               /* Attribute index */
  int32 l_n_dtype;                /* HDF number type from SDattrinfo() */
  char dtype[DATATYPELENMAX];     /* data type string */
  char l_attribute[MAX_NC_NAME+1];/* local buffer to hold attribute name */


  /* initialize local variables */
  nelements_in = 0;
/*  sds_index = 0;
  sds_id = 0;*/
  attr_index = 0;
  l_n_dtype = 0;
  l_count = 0;
  memset(dtype, 0, DATATYPELENMAX);
  memset(l_attribute, 0, MAX_NC_NAME + 1);

  /* initialize return status code */
  status_code = MAPIOK;

  /* Input checks: */
  if ( n_elements == NULL )
  {
    sprintf(buff,"ERROR: getMODDISarinfo unable continue with empty\n"
			"\t n_elements.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  nelements_in = *n_elements;
  *n_elements = 0L;

  if ( NULLstr(attribute) )
  {
    sprintf(buff,"ERROR: getMODISarinfo unable to access an array\n"
			"\t attribute without an attribute name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLstr(arrayname) )
  {
    sprintf(buff,"ERROR: getMODISarinfo unable to access the %.*s\n"
			"\t attribute without the name of the array it is\n"
			"\t associated with.\n",MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) )
  {
    sprintf(buff,"ERROR: getMODISarinfo unable to access the %.*s\n"
			"\t attribute in the %.*s array with a NULL\n"
			"\t file MODFILE structure.\n",MAX_NC_NAME,
			attribute,MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( data_type == NULL )
  {
    sprintf(buff,"ERROR: getMODISarinfo unable to access the %.*s\n"
			"\t attribute in the %.*s array without data\n"
			"\t type input.\n",MAX_NC_NAME,attribute,
			MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  did = getMODISarrayid(file,arrayname,groupname);

  if ( did == NULL )
  {
    sprintf(buff,"ERROR: getMODISarinfo detected errors from\n"
			"\t getMODISarrayid attempting to read the %.*s\n"
			"\t attribute from the %.*s array.\n",
			MAXATTRNAMELEN,attribute,MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* find attribute */
  if ( ( attr_index = SDfindattr((int32)did->id, attribute)) == FAIL )
  {
    sprintf(buff, "ERROR: getMODISarinfo cannot find local attribute\n"
			"\t %.*s attached to the %.*s array.\n",
            		MAXATTRNAMELEN,attribute,MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    *n_elements = nelements_in;
    status_code = MFAIL;
  }

  else if ( SDattrinfo((int32)did->id,attr_index,l_attribute,&l_n_dtype,& l_count) == FAIL )
  {
    sprintf(buff, "ERROR: getMODISarinfo detected FAIL from HDF procedure\n"
			"\t SDattrinfo attempting to read the %.*s\n"
			"\t attribute from the %.*s array.\n",
                        MAXATTRNAMELEN,attribute,MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    status_code = MFAIL;
  }

  else if ( DFNT_to_datatype(l_n_dtype,dtype) == MFAIL )
  {
    sprintf(buff,"ERROR: getMODISarinfo unable to recognize the HDF\n"
			"\t numerical data type %ld while\n"
			"\t attempting to read the %.*s attribute from\n"
			"\t the %.*s array.\n",(long)l_n_dtype,
			MAXATTRNAMELEN,attribute,MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    status_code = MFAIL;
  }

  else
  {
    if (strcmp(data_type, dtype) != 0) 
    /* (Note: Unsigned attribute stored in the HDF become singed attribute. ) */
      if (strstr(data_type,dtype) == NULL) 
      {
        status_code = MFAIL;
        strcpy(data_type,dtype);
      }

    /* if char string, increment count to include terminating null char */
    if (l_n_dtype == DFNT_CHAR8) 
      l_count++;

    if (l_count > nelements_in) 
      status_code = MFAIL;

    *n_elements = (long int)l_count;
  }/* end of else */

  /* enough buffer to hold value */
  if (status_code != MFAIL) 
  {
    if (value == NULL) 
    {
      sprintf(buff, "ERROR: getMODISarinfo unable to read local array\n"
			"\t attribute without output buffer for\n"
			"\t %.*s.\n",MAXATTRNAMELEN,attribute);
      MAPIERR(buff,funcname);
      status_code = MFAIL;
      *n_elements = 0L;
    }

    else if (SDreadattr((int32)did->id,attr_index,(VOIDP)value) == FAIL) 
    {
      sprintf(buff, "ERROR: getMODISarinfo detected FAIL from HDF\n"
			"\t procedure SDreadattr attempting to\n"
			"\t retrieve the %.*s attribute from the\n"
			"\t %.*s array.\n",MAXATTRNAMELEN,attribute,
			MAX_NC_NAME,arrayname);
      MAPIERR(buff,funcname);
      status_code = MFAIL;
      *n_elements = 0L;
    }

    else if (l_n_dtype == DFNT_CHAR8) 
      ((char *)value)[l_count-1] = '\0';
   
  } /* end of if */

  return(status_code);
}

