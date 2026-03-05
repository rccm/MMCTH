#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mapic.h"

long int MODISsizeof(char const *data_type)

/*****************************************************************************
*!C
*
*!Description: Subroutine MODISsizeof is part of a larger software system
*            called the MODIS Applications Programming Interface (API)
*	     Utility, abbreviated M-API.  The M-API Utility consists of
*	     subroutines which allow MODIS Science Team-supplied software
*	     to read in/write to the NCSA HDF format.
*
*	     The MODIS API uses a set of standard strings to describe the
*	     data types stored in array and table structures.  These 
*	     strings are returned, for example, by the routine getMODISardims
*	     to describe the data type of the target array structure.
*	     MODISsizeof returns the number of bytes required to store a 
*	     data type given this data string.  The input string may be a
*	     series of comma delimited data type strings, in which case
*	     the total number of bytes to store the record described by the
*	     string is returned.
*
*!Input Parameters:	
*     data_type        String of comma-delimited data types.
*                      No space is allowed.  permitted types: 
*                      "int8"    "uint8"
*                      "int16"   "uint16"
*                      "int32"   "uint32"
*                      "int64"   "uint64"
*                      "float32" "float64"
*                      "char *"
*
*!Output Parameters:	NONE
*
*Return values:
*     byte_sum	       number of bytes needed to store data_type
*     MFAIL	       error
*
*Externally defined:	
*                       I8               (mapi.h)
*                       UI8              (mapi.h)
*                       I16              (mapi.h)
*                       UI16             (mapi.h)
*                       I32              (mapi.h)
*                       UI32             (mapi.h)
*                       R32              (mapi.h)
*                       R64              (mapi.h)
*                       TXT              (mapi.h)
*                       NULLstr          (mapic.h)
*                       strcpy           <string.h>
*                       malloc           <stdlib.h>
*                       strlen           <string.h>
*                       MFAIL            (mapi.h)
*                       datatype_to_DFNT (mapic.h)
*                       DFKNTsize        (hproto.h)
*                       free             <stdlib.h>
*
*!Revision History:
*
*       Revision 1.1  1995/7/5
*       Paul Fisher -- Research and Data Systems Corporation
*       SAIC/GSC MODIS Support Office
*       Revised MODISsizeof based on updated PDL.
*
*	Revision 01.00 1995/04/25 1639EST
*	Joan R. Baden/RDC
*
*	Create initial template for MODIS API utility subroutine, MODISsizeof.
*
*!Team-unique Header:
*
*	This software is developed by the MODIS Science Data Support Team for
*	the National Aeronautics and Space Administration, Goddard Space Flight
*	Center, under contract NAS5-32373.
*
*Design Notes:
*       Function doesn't support 64-bit integers, due to the current
*       limitation of DFKNTsize (7/6/95).
        
*
*****************************************************************************/

{
  char *local_string, *next_string;  
  int32 HDF_num_type=0L;  /* HDF number type returned from datatype_to_DFNT */
  long int byte_sum=0L;                           /* output number of bytes */ 

  if(!NULLstr(data_type)){
    local_string = (char *)malloc(strlen(data_type)+1); 
    if(!local_string) return(MFAIL);  

    strcpy(local_string, data_type);
    while (local_string != NULL && byte_sum != MFAIL) {
        next_string = strchr(local_string, ',');
        if (next_string)
           *next_string++ = '\0';
        HDF_num_type = datatype_to_DFNT(local_string);
        if (HDF_num_type==MFAIL)
           byte_sum = MFAIL;
        else
           byte_sum += (long int)DFKNTsize(HDF_num_type);
        local_string = next_string;
    }
    free(local_string);
  }
  return(byte_sum);
}
