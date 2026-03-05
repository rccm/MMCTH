#include <stdio.h>
#include <string.h>
#include "mapic.h"

int parse_string(char target_string[], char const delimiter[], 
                 int n_tokens, char *tokens[])
/****************************************************************************
*!C
*
*!Description: Parse_string sets the tokens pointers to the first        
* n_tokens delimiter delimited string token in target_string.  Each 
* terminating delimiter character in target_string is replaced by a 
* null character so that each token pointed to in tokens may be handled
* as an independent string.  The tokens array should have at least  
* n_tokens elements.  If there are fewer tokens in target_string than
* n_tokens the unused pointers are set to NULL.
*
*!Input Parameters:
*	delimiter	string of acceptable delimiter characters
*	n_tokens	maximum number of pointers to string tokens
*			found in target_string.  Minimum number of 
*			available pointers in tokens array.
* 
*!Output Parameters:
*	tokens		array of pointers to the parsed string tokens from 
*			target_string.		
*
*!Input/Output Parameters
*	target_string	string to be parsed
*
*Return values:	
*	positive number	Number of tokens found in target_string
*	MFAIL			if an input error occurs
*
*Externally defined:
*		   NULL                         <stdio.h>
*                  NULLstr                      "mapic.h"
*		   MFAIL                        "mapi.h"
*
*!Revision History:
* $Log: parse_string.c,v $
* Revision 6.1  2010/07/13 18:03:12  kuyper
* Corrected to reflect fact that the function writes to target_string.
*
* Revision 5.1  2005/04/04 18:49:11  vlin
* avoid use of C function strtok().
*
* Revision 1.1 6/26/95
* Paul S. Fisher
* fisher@modis-xl.gsfc.nasa.gov
* Complete rewrite of parse_string
* 
*!Team-unique Header:
*
*          This software is developed by the MODIS Science Data Support
*          Team for the National Aeronautics and Space Administration,
*          Goddard Space Flight Center, under contract NAS5-32373.
*
*!END
*****************************************************************************/

{
  int  found_tokens = 0;       /* numbers of tokens found*/
  register int i = 0;          /* generic counter */
  char *temp = NULL;

  /* If tokens = NULL or n_tokens < 1 or delimeter = NULL or delimeter
     is an empty string, return (MFAIL) */
  if (tokens == NULL || n_tokens < 1 || NULLstr(delimiter))
     return(MFAIL);

  temp = tokens[0] = ( char *)target_string;
  if (temp) found_tokens++;

  for (i=1; i<n_tokens; i++) {
      tokens[i] = NULL;
      temp = strchr(temp, (int)*delimiter);
      if (temp){
         *temp++ = '\0';
         tokens[i] = temp;
         found_tokens++;
      }
      else
         break;
  }
  return(found_tokens);
}
