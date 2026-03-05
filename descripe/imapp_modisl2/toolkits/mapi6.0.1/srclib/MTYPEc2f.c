#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"

int MTYPEc2f(char *c_data_type, char *f_data_type, long int *f_length)
/*
!C**********************************************************************
* 
*!Description: Convert comma delimited M-API C data type string to M-API
*              FORTRAN data type string.
*              FunctionMTYPEc2f is part of a larger software system 
*              called the MODIS Applications Programming Interface
*              (API) Utility, abbreviated M-API.  The M-API Utility con-
*              sists of subroutines which allow MODIS Science Team-sup-
*              plied software to read and write data and metadata from/to
*              HDF files.  The functionality of the M-API is defined in 
*              the MODIS API Specification.
*
*              The MODIS API uses a set of standard strings to describe
*              the data types stored in array and table structures.  For
*	       the same data type, the name is different between M-API C
*              and FORTRAN interface.  Because part of the FORTRAN inter
*              face calls the C interface by using the wrapping method, 
*              the data type string returned from the C interface is in 
*              C notation.  This routine is a low level C function which
*              will convert the C data type notation to FORTRAN data type 
*	       notation.  The function will return MFAIL if c_data_type
*	       contains an unrecognized data type or the memory of f_data_
*	       type is not large enough to hold the FORTRAN data string.
*
*!Input Parameters:
*            c_data_type IN:   String of comma-delimited data types.
*
*			 Permitted C data types are:
*				 "int8"
*				 "unit8"
*				 "int16"
*				 "uint16"
*				 "int32"
*				 "uint32"
*				 "int64"
*				 "float32"
*				 "float64"
*				 "char *"
*	     f_length    IN/OUT: Address of the memory size(bytes) of the
*			 f_data_type string.  Output will replaced with
*			 the actual length of f_data_type + 1.  If an error
*			 occurs f_length will not be changed.  If the f_
*			 data_type is not large enough, f_length will be
*			 returned with the actual size needed to hold the 
*			 FORTRAN data type string.  If the memory size of
*			 f_data_type is three times as large as the string
*			 length of f_data_type, the memory will always be
*			 large enough.  To get the actual size needed to 
*			 hold the FORTRAN data type string, user may set
*			 the f_data_type to NULL and f_length to 0.
*
*!Output Parameters:
*	     f_data_type OUT:  Memory to hold comma-delimited M-API
*			       FORTRAN data types.  If f_length is set
*			       to zero, this parameter can be set to NULL.	
*         
* Return values:    
*             MFAIL     the f_data_type is not large enough or  
*                       an unrecognizable data type in c_data_type string
*             MAPIOK    successful             
*
* Externally defined:
*                          MFAIL               (mapi.h)
*                          MAPIOK              (mapi.h)
*                          FAIL                (hdf.h)
*                          NULLstr             (mapic.h)
*			   FDATATYPELENMAX     (mapic.h)
*			   DATATYPELENMAX      (mapic.h)
*			   I8		       (mapi.h)
*			   UI8		       (mapi.h)
*			   I16		       (mapi.h)
*			   UI16		       (mapi.h)
*			   I32		       (mapi.h)
*			   UI32		       (mapi.h)
*			   R32		       (mapi.h)
*			   R64		       (mapi.h)
*			   TXT		       (mapi.h)
*			   FI8		       (mapic.h)
*			   FUI8		       (mapic.h)
*			   FI16		       (mapic.h)
*			   FUI16	       (mapic.h)
*			   FI32		       (mapic.h)
*			   FUI32	       (mapic.h)
*			   FR32		       (mapic.h)
*			   FR64		       (mapic.h)
*			   FTXT		       (mapic.h)
*!Revision History:
* $Log: MTYPEc2f.c,v $
* Revision 5.1  2005/04/04 16:42:59  vlin
* remove use of function strtok().
*
* Revision 4.1  2003/01/27 19:17:18  kuyper
* Corrected size of allocation for C data type.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
*!Team-unique Header:
*
*      This software is developed by the MODIS Science Data Support
*      Team for the National Aeronautics and Space Administration,
*      Goddard Space Flight Center, under contract NAS5-32373.
*
*!END
**************************************************************************/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];   /* buffer to hold error/warning 
                                            message */
  char *funcname="MTYPEc2f";             /* name of this module */
  int   f_len;
  char  f_type[FDATATYPELENMAX + 1];     /* string to hold a single 
                                            FORTRAN data type string
                                            sometimes with a comma */
  char  *strptr;               /* copy of c_data_type */
  char  *data_type;         /* string to hold a single c data type string */

  f_len = 1;
  *f_type = '\0';

  if ( NULLstr(c_data_type) || (f_length == NULL) )
  { 
    sprintf(buff, "ERROR: MTYPEc2f unable to use empty input parameter."
                          "\n");
    MAPIERR(buff,funcname);
    return(MFAIL);    /* error condition */
  }

  if ( (*f_length > 0) && (f_data_type == NULL) )
  {
    sprintf(buff, "ERROR: MTYPEc2f unable to continue with NULL\n"
			  "\t f_data_type and larger than 0 f_length\n");
    MAPIERR(buff,funcname);
    return(MFAIL);    /*error condition */
  }

  if ( *f_length > 0 ) {
     *f_data_type = '\0';
     *f_length = 0L;
  }

  /* Copy c_data_type string into dynamically allocated memory, perform
     parsing on this copy. */
  if ( (strptr = (char *)malloc((size_t)(strlen(c_data_type)+1))) )
    strcpy(strptr,c_data_type);
  else
  {
    sprintf(buff, "ERROR: MTYPEc2f unable to continue with failure\n"
				"\t of dynamic memory allocation\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* Use strchr to point to the first comma-delimited data type string. */
  while (strptr != NULL) {
    data_type = strchr(strptr,',');
    if (data_type)
       *data_type++ = '\0';
  
  /* if this string is not NULL or empty */
    if (strcmp(strptr,I8) == 0 )
      strcat(f_type,FI8);
    else if (strcmp(strptr,UI8) == 0 )
      strcat(f_type,FUI8);
    else if (strcmp(strptr,I16) == 0 )
      strcat(f_type,FI16);
    else if (strcmp(strptr,UI16) == 0 )
      strcat(f_type,FUI16);
    else if (strcmp(strptr,I32) == 0 )
      strcat(f_type,FI32);
    else if (strcmp(strptr,UI32) == 0 )
      strcat(f_type,FUI32);
    else if (strcmp(strptr,R32) == 0 )
      strcat(f_type,FR32);
    else if (strcmp(strptr,R64) == 0 )
      strcat(f_type,FR64);
    else if (strcmp(strptr,TXT) == 0 )
      strcat(f_type,FTXT);
    else
    {
       sprintf(buff, "ERROR: MTYPEc2f finds unknown M-API C data\n"
               "\t type %.*s\n",DATATYPELENMAX,strptr);
       MAPIERR(buff,funcname);
       free(strptr);
       return(MFAIL);    /*error condition */
    }
 
    /* Increase f_len. */
    f_len += (int)strlen(f_type);
    strcat(f_data_type,f_type);
    strcpy(f_type,",");

    strptr = data_type;
  }

  /* release string copy */
  free(strptr);

    *f_length = f_len;
    return(MAPIOK);
}
