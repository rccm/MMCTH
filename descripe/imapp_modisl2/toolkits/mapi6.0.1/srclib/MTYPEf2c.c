#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h" 

int MTYPEf2c(char *f_data_type, char *c_data_type, long int *c_length)
/*
**********************************************************************
*!C
* 
*!Description: Convert comma delimited M-API FORTRAN data type string 
*	       to M-API C data type string.
*              Function MTYPEf2c is part of a larger software system 
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
*              and FORTRAN interface.  Because part of the FORTRAN interface
*              calls the C interface by using the wrapping method, the data
*              type string returned from the C interface is in C nota-
*	       tion.  This routine is a low level C function which will
*              convert the FORTRAN data type notation to data type 
*	       notation.  The function will return MFAIL if c_data_type
*	       contains an unrecognized data type or the memory of c_data_
*	       type is not large enough to hold the C data string.
*
*!Input parameters:
*            f_data_type IN:   String of comma-delimited data types.
*
*			 Permitted FORTRAN data types are:
*				  'INTEGER*1'
*				  'UINTEGER*1'
*				  'INTEGER*2'
*				  'UINTEGER*2'
*				  'INTEGER*4'
*				  'UINTEGER*4'
*				  'REAL*4'
*				  'REAL*8'
*				  'CHARACTER*(*)'
*	     c_length    IN/OUT: Address of the memory size(bytes) of the
*			 c_data_type string.  Output will replaced with
*			 the actual length of c_data_type + 1.  If an error
*			 occurs c_length will not be changed.  If the c_
*			 data_type is not large enough, c_length will be
*			 returned with the actual size needed to hold the 
*			 C data type string.  If the memory size of
*			 c_data_type is twice as large as the string
*			 length of c_data_type, the memory will always be
*			 large enough.  To get the actual size needed to 
*			 hold the C data type string, user may set
*			 the c_data_type to NULL and c_length to 0.
*
*!Output parameters:
*	     c_data_type OUT:  Memory to hold comma-delimited M-API
*			       C data types.  If c_length is set
*			       to zero, this parameter can be set to NULL.	
*         
*Return values:    
*          MFAIL     the c_data_type is not large enough or  
*                    an unrecognizable data type in f_data_type string
*          MAPIOK    successful             
*
*Externally defined:
*                          MFAIL               (mapi.h)
*                          MAPIOK              (mapi.h)
*                          FAIL                (hdf.h)
*                          NULLstr             (mapic.h)
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
* !Revision history:
* $Log: MTYPEf2c.c,v $
* Revision 6.1  2010/07/13 19:33:24  kuyper
* Clarified code by moving assignment expression out of if() condition.
*
* Revision 5.1  2005/04/04 16:42:59  vlin
* remove use of function strtok().
*
*
* !Team-unique header:
*
*       This software is developed by the MODIS Science Data Support
*       Team for the National Aeronautics and Space Administration,
*       Goddard Space Flight Center, under contract NAS5-32373.
* 
*!END
**************************************************************************/

{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];   /* buffer to hold error/warning 
                                            message */
  char *funcname="MTYPEf2c";             /* name of this module */
  int   c_len;
  char  c_type[DATATYPELENMAX + 1];     /* string to hold a single 
                                            c data type string
                                            sometimes with a comma */
  char  *strptr;               /* copy of f_data_type */
  char  *data_type;         /* string to hold a single f data type string */

  c_len = 1;
  *c_type = '\0';

  if ( NULLstr(f_data_type) || (c_length == NULL) )
  { 
    sprintf(buff, "ERROR: MTYPEf2c unable to use empty input parameter."
                          "\n");
    MAPIERR(buff,funcname);
    return(MFAIL);    /* error condition */
  }

  if ( (*c_length > 0) && (c_data_type == NULL) )
  {
    sprintf(buff, "ERROR: MTYPEf2c unable to continue with NULL\n"
			  "\t c_data_type and larger than 0 c_length\n");
    MAPIERR(buff,funcname);
    return(MFAIL);    /*error condition */
  }

  if ( *c_length > 0 )
    *c_data_type = '\0';

  /* Copy f_data_type string into dynamically allocated memory, perform
     parsing on this copy. */
  strptr = (char*)malloc((size_t)(strlen(f_data_type)+1)); 
  if ( strptr )
    strcpy(strptr,f_data_type);
  else
  {
    sprintf(buff, "ERROR: MTYPEf2c unable to continue with failure\n"
			"\t of dynamic memory allocation\n");
    MAPIERR(buff,funcname);
    return MFAIL;    /*error condition */
  }

  /* Use strchr to point to the first comma-delimited data type string. */
  while (strptr != NULL) {
    data_type = strchr(strptr, ',');
    if (data_type)
       *data_type++ = '\0';

  /* if this string is not NULL or empty */
    if (strcmp(strptr,FI8) == 0 )
      strcat(c_type,I8);
    else if (strcmp(strptr,FUI8) == 0 )
      strcat(c_type,UI8);
    else if (strcmp(strptr,FI16) == 0 )
      strcat(c_type,I16);
    else if (strcmp(strptr,FUI16) == 0 )
      strcat(c_type,UI16);
    else if (strcmp(strptr,FI32) == 0 )
      strcat(c_type,I32);
    else if (strcmp(strptr,FUI32) == 0 )
      strcat(c_type,UI32);
    else if (strcmp(strptr,FR32) == 0 )
      strcat(c_type,R32);
    else if (strcmp(strptr,FR64) == 0 )
      strcat(c_type,R64);
    else if (strcmp(strptr,FTXT) == 0 )
      strcat(c_type,TXT);
    else
    {
       sprintf(buff, "ERROR: MTYPEf2c finds unknown M-API\n"
               "\t FORTRAN data type %.*s\n",FDATATYPELENMAX,strptr);
       MAPIERR(buff,funcname);
       free(strptr);
       return(MFAIL);    /*error condition */
    }
 
    /* Increase c_len. */
    c_len += strlen(c_type);
    strcat(c_data_type,c_type);
    strcpy(c_type,",");

    strptr = data_type;
  }

  /* release string copy */
  free(strptr);

  *c_length = c_len;
  return(MAPIOK);
}
