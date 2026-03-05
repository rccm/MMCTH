#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"
#define MAXDIMNAMELEN 100

int getMODISdiminfo(MODFILE *file, char const *arrayname, char const *groupname,
                    long int dimension, char const *attribute, char *data_type,
                    long int *n_elements, void *value)
/*
!C************************************************************************
*!Purpose:     Reads the value(s) of a local attribute attached to a spe- 
*              cific dimension of an SDS(array) from a MODIS HDF file. 
*
*!Description: Function getMODISdiminfo is part of a larger software sys-
*              tem called the MODIS Applications Programming Interface
*              (API) Utility, abbreviated M-API.  The M-API Utility con_
*              sists of subroutines which allow MODIS Science Team-sup-
*              plied software to read and write data and metadata from/to
*              HDF files.  The functionality of the M-API is defined in 
*              the MODIS Application Program Interface(API) Specification.
*
*              getMODISdiminfo (GMDMIN) retrieves the value(s) associated
*              with an attribute = value pair for a dimension attribute by
*              giving the array's name and dimension number and the attri-
*              bute name.  If the attribute cannot be found, the routine 
*              will return MFAIL and the passed variable unchanged.
*
*              The routine will also fail if the provided data_type is 
*              found to be different from the attribute's data type or the
*              n_elements is found to be too small to contain the attri-
*              bute's value.  getMODISdiminfo replaces this input informa-
*              tion with the actual data type and number of elements con-
*              tained in the attribute value(in the case of character data
*              , it is the length of the string, including the '\0' cha-
*              racter).  These attribute's attribute may be used to pro-
*              perly retrieve the attribute value with a second call to 
*              the routine.  If a function failure occurs or specified 
*              array or dimension does not exist, n_elements will be set
*              to zero.
*
*              A variable of the proper data type should be passed for the
*              value parameter.  The data type information required to pro
*              perly use either routine may be found in Appendix A, M-API-
*              Supplied Constants, and Appendix C, MODIS Data Product File
*              Definitions of M-API Use's Guide.  Appendix A has a listing
*              for each M-API provided attribute that includes the data
*              type, the format, and/or specific values associated with it.
*
*!Input Parameters: 
*       file         IN: Address of MODFILE structure that is used to 
*                    reference a MODIS HDF file containing the attribute.
*       arrayname    IN: ASCII string name of the array.  Provided macros 
*                    for accepted MODIS HDF file array names are listed in
*                    Appendix A, M-API-Supplied Constants.
*       groupname    IN: ASCII string name of the data group containing 
*                    the array structure to which the attribute is atta-
*                    ched.  If set to NULL, the entire file will be sear-
*                    ched for the array structure named arrayname.
*       dimension    IN: The dimension number from which the attribute is
*                    retrieved(0-based).
*       attribute    IN: ASCII string name of the attribute.  Provided 
*                    macros for accepted MODIS HDF file attribute names 
*                    are listed in Appendix A, M-API-Supplied Constants.
*       data_type    IN/OUT: Address of data type of the value output.  
*                    Output replaces with the data type of the retrieved 
*                    attribute.  The memory size of this arguments should
*                    be at least 8 bytes long.
*
*                    Permitted C data types:
*                             "int8"
*                             "uint8"
*                             "int16"
*                             "uint16"
*                             "int32"
*                             "uint32"
*                             "int64"
*                             "float32"
*                             "float64"
*                             "char *"
*
*      n_elements   IN/OUT: Address of the number of elements the  value
*                   buffer can contain.  Output replaces with the number
*                   of elements required to contain the attribute.  If a
*                   function failure occurs, the value will be set to zero.
*
*!Output Parameters: 
*      value        OUT: Address of value associated with the attribute.
*
* Returns:     MAPIOK if successful, MFAIL if value cannot contain the 
*              retrieved attribute value, the data type is different, the
*              attribute cannot be found, or an error occurs.
*
* External references:
*                      MAPIOK                   (mapi.h)
*                      MFAIL                    (mapi.h)
*                      NULLstr                  (mapic.h)
*                      DFNT_to_datatype         (mapic.h)
*                      DFNT_CHAR8               (hdf.h)
*                      DATATYPELENMAX           (mapic.h)
*                      MAX_NC_NAME              (netcdf.h)
*                      SDfindattr               (mfhdf.h)
*                      SDSINFO                  (mapic.h)
*                      SDattrinfo               (mfhdf.h)
*                      SDreadattr               (mfhdf.h)
*                      SDgetdimid               (mfhdf.h)
*		       MODFILE		        (mapi.h)
*		       PGS_SMF_MAX_MSGBUF_SIZE  (PGS_SMF.h)
*		       NULLMODFIL               (mapic.h)
*		       MAPIERR                  (mapic.h)
*                      DATAID                   (mapic.h) 
*                      getMODISarrayid          (mapic.h)
*
*!Revision History: 
* $Log: getMODISdiminfo.c,v $
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.2  1996/08/16  11:46:15  fshaw
 * version 2.1
 * updated to use the ring super structure.
 *
 * Revision 1.2  1996/08/16  11:46:15  fshaw
 * version 2.1
 * updated to use the ring super structure.
 *
 * Revision 1.1  1996/08/12  17:24:39  fshaw
 * Initial revision
 *
 * Revision 1.5  1996/03/25  14:26:59  qhuang
 * removed redundant 'arrayname' in error message, cast variable 'count',
 * added 'L' to assignment of '*n_elements'.
 *
 * Revision 1.4  1996/03/11  22:23:01  qhuang
 * Version 2 design, group handling, new messages, input checking, added error
 * checking after DFNT_to_datatype.
 *
 * Revision 1.3  1995/10/30  20:43:47  qhuang
 * Added capability to pass status message to log files.
 *
*
*                      Revision 01.00 1995/10/5
*                      Qi Huang
*                      Version 1.3 original development/testing.
*
*!Team-unique Header: 
*              This software is developed by the MODIS Science Data Sup-
*              port Team for the National Aeronautics and Space Admini-
*              stration, Goddard Space Flight Center, under contract 
*              NAS5-32373.
*
*              Portions developed at the National Center for Supercompu-
*              ting Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!References and Credits:
*
*!Design Notes:
*
!END**************************************************************************
*/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char *funcname="getMODISdiminfo"; 
  int  status;             /* variable to hold the returned value */
  long int  nelements_in;  /* variable to keep the value of n_elements
                              passed in  */
  int32  dim_id = 0L;
  int32  attr_index = 0L;
  char   attr_name[MAX_NC_NAME];
  int32  num_type = 0L; /* number type */
  int32  count = 0L;
  char   dtype[DATATYPELENMAX];
  DATAID *did;


  status = MFAIL;        /* error */

  /* Input checks: */
  if ( n_elements == NULL ) {
    sprintf(buff,"ERROR: getMODISdiminfo unable continue with empty\n"
			"\t n_elements.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  nelements_in = *n_elements;
  *n_elements = 0L;

  if ( NULLstr(attribute) ) {
    sprintf(buff,"ERROR: getMODISdiminfo unable to access an array\n"
			"\t dimension without an attribute name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLstr(arrayname) ) {
    sprintf(buff,"ERROR: getMODISdiminfo unable to access dimension\n"
			"\t attribute %.*s attribute without the name of the\n"
			"\t array it is associated with.\n",
			MAXDIMNAMELEN,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) ) {
    sprintf(buff,"ERROR: getMODISdiminfo unable to access dimension\n"
			"\t attribute %.*s in the %.*s array with\n"
			"\t a NULL file MODFILE structure.\n",
			MAXDIMNAMELEN,attribute,MAXDIMNAMELEN,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( data_type == NULL ) {
    sprintf(buff,"ERROR: getMODISdiminfo unable to access dimension\n"
			"\t attribute %.*s in the %.*s array\n"
                        "\t without data type input.\n",
                         MAXDIMNAMELEN,attribute,MAXDIMNAMELEN,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  did = getMODISarrayid( file, arrayname, groupname);

  if ( did == NULL) {
    sprintf(buff,"ERROR: getMODISdiminfo detected FAIL from M-API\n"
			"\t internal function getMODISarrayid attempting to\n"
			"\t read the %.*s attribute.\n",
			     MAXDIMNAMELEN,attribute);
    MAPIERR(buff,funcname);

  } else if ( dimension < 0 || dimension > ((SDSINFO *)did->info)->rank - 1) {
    sprintf(buff,"ERROR: getMODISdiminfo unable to retrieve an %.*s\n"
			"\t attribute for %ld dimension The %.*s \n"
			"\t array has %ld dimensions.\n",
		        MAXDIMNAMELEN,attribute, dimension, MAXDIMNAMELEN,
                        arrayname,(long)((SDSINFO *)did->info)->rank);
    MAPIERR(buff,funcname);
    
  }else if ( (dim_id = SDgetdimid((int32)did->id, (int32)dimension)) == FAIL ){
         sprintf(buff,"ERROR: getMODISdiminfo detected FAIL from HDF \n"
                      "\t procedure SDgetdimid attempting to read the \n"
                      "\t %.*s attribute.\n",
                      MAXDIMNAMELEN,attribute);
         MAPIERR(buff,funcname);

  }else if ( (attr_index = SDfindattr( dim_id, attribute)) == FAIL ){
         /* can't find the attribute */
         sprintf(buff,"ERROR: getMODISdiminfo cannot find local array\n"
                      "\t dimension attribute %.*s\n",
                      MAXDIMNAMELEN, attribute);
         MAPIERR(buff,funcname);

         *n_elements = nelements_in;

  }else if ( SDattrinfo(dim_id, attr_index, attr_name, &num_type, &count) == FAIL ) {
         sprintf(buff,"ERROR: getMODISdiminfo detected FAIL from HDF\n"
                      "\t procedure SDattrinfo attempting to read the\n"
                      "\t %.*s attribute.\n",
                      MAXDIMNAMELEN, attribute);
         MAPIERR(buff,funcname);

  }else  if ( DFNT_to_datatype(num_type,dtype) == MFAIL ) {
         sprintf(buff,"ERROR: getMODISdiminfo unable to recognize the HDF\n"
                      "\t numerical data type %ld while attempting to read \n"
                      "\t the %.*s dimension attribute from the %.*s array.\n",
                 (long)num_type, MAXDIMNAMELEN, attribute, MAXDIMNAMELEN, arrayname);
         MAPIERR(buff,funcname);
  }

  else {
    status = MAPIOK;

    if ( strcmp(data_type,dtype) != 0 )
    /* (NOTE: Unsigned attributes stored in the HDF become signed attribute.)*/

      if ( strstr(data_type, dtype) == NULL )
      {
        strcpy(data_type, dtype);
        status = MFAIL;
      }

    if ( num_type == DFNT_CHAR8 )
      count ++;

    if ( (long int)count > nelements_in )
      status = MFAIL;            /* attribute retrieval failed. */
 
    /* Set the value of n_elements to the number of values. */
    *n_elements = (long int)count;
  }

  if ( status != MFAIL ) {
  /* the value buffer should be capable of holding
      the attribute's values(s)). */

    if ( value == NULL ) {
      sprintf(buff,"ERROR: getMODISdiminfo unable to read local\n"
                   "\t array dimension attribute without output\n"
                   "\t buffer for %.*s\n", MAXDIMNAMELEN, attribute);
      MAPIERR(buff,funcname);
      status = MFAIL;
      *n_elements = 0L;

    } else if ( SDreadattr(dim_id, attr_index, value) == FAIL ) {
      sprintf(buff,"ERROR: getMODISdiminfo detected FAIL from HDF\n"
                   "\t procedure SDreadattr attempting to read\n"
                   "\t the %.*s attribute.\n", MAXDIMNAMELEN, attribute);
      MAPIERR(buff,funcname);
      status = MFAIL;          /* attribute retrieval failed. */
      *n_elements = 0L; 

    } else if ( num_type == DFNT_CHAR8 ) 
      /* insert string terminator */
      *((char *)value + count - 1) = '\0'; 
  }

  return(status);
}
