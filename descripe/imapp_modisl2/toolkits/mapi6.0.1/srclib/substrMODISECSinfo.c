#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"
int substrMODISECSinfo(char const *char_value, long int n_elements,
                       long int *n_strings, char *substr[])
/*
!C**********************************************************************
* 
*!Purpose:  Decomposes a multiple character string ECS metadata (retrieved by 
* getMODISECSinfo) into individual strings.
* 
*!Description:  ECS metadata values may be integer, floating point, or character
*               string values or arrays of values.  Some may be multiple strings.
*               The routine getMODISECSinfo retrieves such strings into a one-
*               dimension character array with the individual strings separated 
*               by nulls('\0').  substrMODISECSinfo breaks this ‘packed’ 
*               character array into its constituent substrings.  
*               substrMODISECSinfo sets the pointers in a provided output array 
*               to the beginning of each substring in the char_value array.
*               
*!Input parameters:
* char_value	IN: 	Character string containing t'packed' muultiple 
*		 substrings of ECS metadata retrieved with getMODISECSinfo.
* 		Do not deallocate char_value until substr array gets correct 
*		values.
* n_elements	IN: The composite output dimensions, from getMODISECSinfo, 
* 		containing (in the case of character string metadata) the total 
* 		length (in bytes) of the string in char_value in its lower two bytes 
* 		and the number of substrings packed into char_value  - 1 in the 
* 		upper two bytes.  The calculations
* 		int n_strings  =  n_elements/65536  + 1
* 		n_bytes  =  n_elements%65536
* 		provide the number of substrings and the total length, 
* 		respectively, of the data in char_value.  When there is only one 
* 		string in char_value, n_elements will be less than 65536 and 
* 		there is no need to use substrMODISECSinfo.
* n_strings	IN/OUT: Address of the number of pointers available  in the 
* 		substr  array.  The substr pointers will not be set to the substrings 
* 		in char_value unless there are sufficient pointers available in the 
* 		pointer array.  substrMODISECSinfo  replaces this input with the 
* 		number of substrings pointers have been set to in the 
* 		char_value array.  *n_strings will be set to 0 if a function error 
* 		occurs.  This argument must not be the address of a  constant.
*!Output Parameters:
* substr	OUT:	Array of poiners to the constituent substrings contained 
*		in the char_value array.
* 
* Returns:	MAPIOK (0) if the address of each substring in char_value is in 
*		the substr pointer array,  MFAIL (-1) if substr  contains insufficient 
*		pointers for the addresses of all of the substrings in char_value or 
*		an error occurs. 
* 
* External references:
*		PGS_SMF_MAX_MSGBUF_SIZE		(mapic.h)
*		NULLstr				(mapic.h)
*		MAPIERR				(mapic.h)
*		MAPIOK				(mapi.h)
*		MFAIL				(mapi.h)
*		MAX_NC_NAME			(netcdf.h)
*
*!Revision History:
*	 Qi Huang, RDC		May 16, 1996
* 	 Version 2.0
*	 Original development and testing
* $Log: substrMODISECSinfo.c,v $
* Revision 6.1  2010/07/12 16:15:46  kuyper
* Removed useless 'const' on multidimensional array parameter.
*
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.2  1999/04/26 17:01:49  jayshree
* reorganizing mapi RCS
*
 * Revision 1.4  1998/01/13  17:56:26  fshaw
 * corrected defect for reading NULL ECS-metadata strings
 *
 * Revision 1.4  1998/01/13  17:56:26  fshaw
 * corrected defect for reading NULL ECS-metadata strings
 *
 * Revision 1.3  1996/06/04  13:30:00  qhuang
 * Added comments to local variables.
 *
 * Revision 1.2  1996/06/03  21:38:46  qhuang
 * Added int before n_strings = ... in prolog, removed commented block.
 *
 * Revision 1.1  1996/05/17  13:44:41  qhuang
 * Initial revision
 *
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
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning message */
  char  *funcname="substrMODISECSinfo";    /* name of this routine */
  int32 loc_n_strings;			/* local variable for number of substrings */
  char const *p;				  /* pointer used to parse strings */
  int   status;				/* return status */
  int   i;				  /* for loop counter */

  status = MFAIL;	/* error */

  /* Input checks: */
  if (char_value == NULL) 
  {
    sprintf(buff,"ERROR: substrMODISECSinfo unable to continue without\n"
			"\t char_value input\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( (n_strings == NULL) || (*n_strings < 1) )
  {
    sprintf(buff,"ERROR: substrMODISECSinfo unable to continue without\n"
			"\t n_strings input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if (substr == NULL)
  {
    sprintf(buff,"ERROR: substrMODISECSinfo unable to continue without\n"
			"\t substr pointer array.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if (n_elements < 0)
  {
    sprintf(buff,"ERROR: substrMODISECSinfo unable to continue with\n"
			"\t invalid %ld n_elements\n",n_elements);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  loc_n_strings = (int32)(n_elements/65536 + 1);

  if ( *n_strings < loc_n_strings )
  {
    sprintf(buff,"ERROR: substrMODISECSinfo unable to fit %ld\n"
			"\t substrings into %ld pointers %.*s array.\n",
		(long)loc_n_strings,*n_strings,MAX_NC_NAME,(char *)substr);
    MAPIERR(buff,funcname);
    *n_strings = loc_n_strings;
    return(MFAIL);
  }

  /* Parse char_value string */
  p = char_value;
  substr[0] = (char*)p;
  for (i=1;i<loc_n_strings;i++)
  {
    while (*p != '\0')
      p++;
    p++;
    substr[i] = (char*)p;
  }

  *n_strings = loc_n_strings;
  status = MAPIOK;
  return(status);
}


