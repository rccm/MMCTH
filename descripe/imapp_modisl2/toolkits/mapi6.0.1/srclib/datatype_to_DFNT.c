#include "mapic.h"
#include <string.h>

int32 datatype_to_DFNT(char const *datatype)
/*
!C**********************************************************************
* Routine: datatype_to_DFNT
*
*!Purpose: Returns the HDF number type associated with a M-API data
*          type string.
* 
*!Description: Subroutine datatype_to_DFNT is part of a larger
*              software system called the MODIS Applications
*              Programming Interface (API) Utility, abbreviated M-API.
*              The M-API Utility consists of subroutines which allow
*              MODIS Science Team supplied software to read and write
*              data from/to the HDF format.
*
*              datatype_to_DFNT takes a M-API data type string and
*              returns the HDF number type integer used by HDF routines 
*              to identify number types.
* 
* !Input Parameters:
*      char *datatype      pointer to M-API data type string 
*                          permitted types: "int8"
*                                           "uint8"
*                                           "int16"
*                                           "uint16"
*                                           "int32"
*                                           "uint32"
*                                           "int64"
*                                           "uint64"
*                                           "float32"
*                                           "float64"
*                                           "char *"
* 
* !Output Parameters:None
*
* Returns:    HDF number type or MFAIL if an error occurs
*      
* !Revision History:
* 
*      Revision 1.1   1995/06/04
*      Paul Fisher -- Research and Data Systems Corporation
*      SAIC/GSC MODIS support office
*      Revised to reflect updated PDL, includes support for unsigned 
*      integers.
*
*      Revision 01.01 1995/04/25
*      Angela Sigmund 
*      Added minor comments. 
*
*      Revision 01.00 1995/04/19
*      Angela Sigmund 
*      Initial creation. 
* 
* !Team-unique Header:
* 
* !References and Credits:
*      This software is developed by the MODIS Science Data Support
*      Team for the National Aeronautics and Space Administration,
*      Goddard Space Flight Center, under contract NAS5-32373.
* 
*      HDF portions developed at the National Center for Supercomputing
*      Applications at the University of Illinois at Urbana-Champaign.
*
* Externals:     	DFNT_INT8       (hdf.h)
*                       DFNT_UINT8      (hdf.h)
*			DFNT_INT16	(hdf.h)
*                       DFNT_UINT16	(hdf.h)
*			DFNT_INT32	(hdf.h)
*                       DFNT_UINT32     (hdf.h)
*               	DFNT_INT64      (hdf.h)
*                       DFNT_UINT64     (hdf.h)
*               	DFNT_FLOAT32    (hdf.h)
*               	DFNT_FLOAT64    (hdf.h)
*               	DFNT_CHAR8      (hdf.h)
*			I8		(mapi.h)
*                       UI8             (mapi.h)
*			I16             (mapi.h)
*                       UI16            (mapi.h)
*                       I32             (mapi.h)
*                       UI32            (mapi.h)
*                       I64             (mapi.h)
*                       UI64            (mapi.h)
*                       R32             (mapi.h)
*                       R64             (mapi.h)
*                       TXT             (mapi.h)
*                       strcmp          <string.h>
*                       MFAIL           (mapi.h)
*
* !Design Notes:
*
!END********************************************************************
*/
{
  int32 output=0;   /* used to return HDF datatype */

  /* if datatype points to a nullstring return MFAIL */  
  if (NULLstr(datatype)) return(MFAIL);

  /* select correct HDF datatype */
  if ((strcmp(datatype,I8)==0)) output=DFNT_INT8;
  else if ((strcmp(datatype,UI8)==0)) output=DFNT_UINT8;
  else if ((strcmp(datatype,I16)==0)) output=DFNT_INT16;
  else if ((strcmp(datatype,UI16)==0)) output=DFNT_UINT16;
  else if ((strcmp(datatype,I32)==0)) output=DFNT_INT32;
  else if ((strcmp(datatype,UI32)==0)) output=DFNT_UINT32;
/************************************************************************ 
   64 bit integer support disabled until HDF decides to support it.
   To re-enable, remove comments around the following block of code.

  else if ((strcmp(datatype,I64)==0)) output=DFNT_INT64;
  else if ((strcmp(datatype,UI64)==0)) output=DFNT_UINT64;
************************************************************************/
  else if ((strcmp(datatype,R32)==0)) output=DFNT_FLOAT32;
  else if ((strcmp(datatype,R64)==0)) output=DFNT_FLOAT64;
  else if ((strcmp(datatype,TXT)==0)) output=DFNT_CHAR8;
  else output=MFAIL;

  return(output);

} /* end datatype_to_DFNT */

