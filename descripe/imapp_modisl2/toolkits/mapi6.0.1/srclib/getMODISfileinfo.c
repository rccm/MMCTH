/*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"

/*****************************************************************************/

int getMODISfileinfo (MODFILE *file, char const *attribute, char *data_type,
                    long int *n_elements, void *value)

/*
!C****************************************************************************
*
*!Purpose:     Reads the value(s) of a global attribute from a MODIS HDF file.
*
*!Description: Function getMODISfileinfo is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
* 	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      subroutines which allow MODIS Science Team-supplied software
*	      to read and write data and metadata from/to HDF files.  The 
*             functionality of the M-API is defined in the MODIS Application 
*             Program Interface (API) specification. 
*
*	      getMODISfileinfo retrieves the value(s) associated with an
*	      attribute = value(s) pair given the attribute name. If the
*	      attribute cannot be found, the routine will return MFAIL and the
*	      passed variable unchanged. 
*
*	      The routine will also fail if the provided data_type is found
*	      to be different with the attribute's data type or the n_elements
*  	      is found to be too small to contain the attribute's value(s).
*	      getMODISfileinfo replaces this input information with the actual
*	      data type and number of elements contained in the attribute value(s)
*	      (in the case of character data, it is the length of the string,
*             including the '\0' character).  These attributes's attributes
*             may be used to properly retrieve the attribute value with a 
*             second call to the routine.  If a function failure occurs, 
*             n_elements will be set to zero.
*
*	      A variable of the proper data type should be passed for the value
*             parameter. The data type information required to properly use
*             either routine may be found in Appendix A, M-API-Supplied
*             Constants, and Appendix C, MODIS Data Product File Definitions
*             of M-API User's Guide.  Appendix A has a listing for each M-API
*             provided attribute's attributes that includes the data type, 
*             the format, and/or specific values associated with it.
*
* !Input Parameters:
* file         IN:  Address of MODFILE structure that is used to reference the
*              MODIS HDF file containing the attribute.
* attribute    IN:  ASCII string name of the attribute.  Provided macros for accepted
*              MODIS HDF file attribute names are listed in Appendix A, MODIS
*              API Supplied Constants.
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
* n_elements   IN/OUT:  Address of the number of elements available in the value array.
*              Output replaces it with the number of elements required to
*              contain the attribute.  If a function failure occurs, the value
*              will be set to zero.  If the requested attribute does not
*	       exist, n_elements will be unchanged.
*
* !Output Parameters:
* data_type    Address of data type of the value output.  Output replaces with
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
* n_elements   Address of the number of elements available in the value array.
*              Output replaces with the number of elements required to
*              contain the attribute.  If a function failure occurs, the value
*              will be set to zero.
*
* value        Address of value(s) associated with the attribute.
*
*Externally defined:	NULLstr		(mapic.h)
*                       DATATYPELENMAX  (mapic.h)
*                       DFNT_to_datatype(mapic.h)
*                       MAX_NC_NAME     (netcdf.h)
*			SDattrinfo	(mfhdf.h)
*			SDfindattr	(mfhdf.h)
*                       SDreadattr      (mfhdf.h)
*			DFNT_CHAR8	(hdf.h)
*                       MFAIL           (mapi.h)
*                       MAPIOK          (mapi.h)
*			MODFILE		(mapi.h)
* 	 	        PGS_SMF_MAX_MSGBUF_SIZE(PGS_SMF.h)
* 	 	        NULLMODFIL(mapic.h)
* 	 	        MAPIERR(mapic.h)
*
* Returns:     MAPIOK if successful, MFAIL if value cannot contain the 
*              retrieved metadata value, the attribute cannot be found, or
*              an error occurs.
*
* !Revision History:
*		Qi Huang	1996/03/20
*		Version 2.0 
*		Original development and testing
*
*$Log: getMODISfileinfo.c,v $
*Revision 6.1  2010/07/13 13:57:33  kuyper
*Removed inappropriate use of NULL in an arithmetic context.
*
*Revision 5.1  2005/04/04 18:49:11  vlin
*constant safe for pointer arguments.
*
*Revision 1.1  1998/02/06 22:26:06  fshaw
*Initial revision
*
 * Revision 1.6  1996/05/08  15:35:55  qhuang
 * Set *n_elements to nelements_in when attribute is not found.
 *
 * Revision 1.5  1996/03/25  14:35:24  qhuang
 * *** empty log message ***
 *
 * Revision 1.4  1996/03/08  19:25:20  qhuang
 * Version 2, set n_elements unchnged when the requested attribute does not
 * exist, separated input checking, new messages, check the return of DFNT_
 * to_datatype.
 *
 * Revision 1.3  1995/11/07  18:39:47  qhuang
 * minor changes
 *
 * Revision 1.2  1995/10/31  16:08:47  qhuang
 * Added capability to pass status messages to log files.
 *
* Revision 1.1  1995/10/31  16:04:40  qhuang
* Initial revision
*
*	1.    Mitchell Weiss/RDC		Jul 95
*	      Original Development / Testing
*     
*       2.    Gary Fu/GSC                       5 OCT 95
*             update based on the updated PDL
*
* !Team-unique Header:
*              This software is developed by the MODIS Science Data Support
*              Team for the National Aeronautics and Space Administration,
*              Goddard Space Flight Center, under contract NAS5-32373.
*
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
*
!END*********************************************************************
*/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning message */
  char *funcname="getMODISfileinfo";    /* name of this routine */
  char attr_name[MAX_NC_NAME] = "";      /* local buf for attribute name */
  int  status;                         /* function return code   */
  long int nelements_in;             /* initial input value of *n_elements */
  int32 number_type = 0;                       /* data type(type of data)  */ 
  int32 count = 0;                             /* number of values         */
  int32 attr_index = 0;                        /* attribute index          */
  char  dtype[DATATYPELENMAX] = "";     /* character number_type    */

  /* initialize return status code */
  status = MAPIOK;

  /* Input checks: */
  if ( n_elements == NULL )
  {
    sprintf(buff,"ERROR: getMODISfileinfo unable to continue with empty\n"
			"\t n_elements.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  nelements_in = *n_elements;
  *n_elements = 0;

  if ( NULLstr(attribute) )
  {
    sprintf(buff,"ERROR: getMODISfileinfo unable to access a file\n"
			"\t attribute without an attribute name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( data_type == NULL )
  {
    sprintf(buff,"ERROR: getMODISfileinfo unable to access the %.*s\n"
			"\t file attribute without data type input.\n",
			MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) )
  {
    sprintf(buff,"ERROR: getMODISfileinfo unable to continue with an\n"
			"\t invalid MODIS file structure input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* find attribute */
  if ( (attr_index = SDfindattr((int32)file->sd_id, attribute) ) == FAIL ) 
  {
    sprintf(buff, "ERROR: %s cannot find file attribute %.*s\n",
			funcname,MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    *n_elements = nelements_in;
    return(MFAIL);
  }

  /* After getting attr_index retrieve the attribute's attribute.   */
  if ( (SDattrinfo((int32)file->sd_id,attr_index,attr_name,&number_type,
                  &count) ) == FAIL )
  {
    sprintf(buff,"ERROR: %s detected FAIL from HDF\n"
			"\t procedure SDattrinfo attempting to read the\n"
			"\t %.*s file attribute.\n",
             		funcname,MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* test attribute components for input and output consistency */
  /* First test data_type     */

  if ( DFNT_to_datatype(number_type, dtype) == MFAIL )
  {
    sprintf(buff,"ERROR: getMODISfileinfo unable to recognize the HDF\n"
			"\t numerical data type %ld while\n"
			"\t attempting to read the %.*s file attribute.\n",
			(long)number_type,MAX_NC_NAME,attribute);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if (strcmp(data_type,dtype) != 0) 
    if ( strstr(data_type,dtype) == NULL )
    {
      status = MFAIL;
      strcpy(data_type,dtype);            /* copy actual data type    */
    }
 

  if (number_type == DFNT_CHAR8)          /* number_type is character */
    count++;
 
  if (count > nelements_in) 
    status = MFAIL;
  
  *n_elements = count;

  if (status != MFAIL) 
  {
    if ( value == NULL )
    {
      sprintf(buff, "ERROR: %s unable to read file\n"
			"\t attribute without output buffer for\n"
			"\t %.*s.\n", funcname,MAX_NC_NAME,attribute);
      MAPIERR(buff,funcname);
      status = MFAIL;
      *n_elements = 0;
    }
    else if ( SDreadattr((int32)file->sd_id,attr_index,value) == FAIL )
    {
      sprintf(buff, "ERROR: %s detected FAIL from HDF\n"
			"\t procedure SDreadattr attempting to read\n"
			"\t the %.*s file attribute.\n",
                        funcname,MAX_NC_NAME,attribute);
      MAPIERR(buff,funcname);
      status = MFAIL;
      *n_elements = 0;
    }
    else if ( number_type == DFNT_CHAR8 )
      ((char *)value)[count-1] = '\0';  
  } /* end of if */

  return status;
}                                       /* end getMODISfileinfo      */

