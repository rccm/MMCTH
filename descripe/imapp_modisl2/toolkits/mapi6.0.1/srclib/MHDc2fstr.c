#include "mapic.h"

intn MHDc2fstr(char *str, intn len)
/*
!C**********************************************************************
* 
*!Purpose:   This routine is temporarily used to replace a same named
* routine in HDF4.0r1p1 because of a bug in the routine.
* 
*!Description:   Function HDc2fstr convert a C string into a Fortran 
* string IN PLACE, the NULL of the C string is ripped off and is padded
* with spaces.
* 
*!Input parameters:
*	str	IN:  C string to be converted.
*	len	IN:  Memory size of FORTRAN string.
*
*!Output parameters:	None
*
* Returns:		MAPIOK
*
* External references:
*
*!Revision history:
* 		Qi Huang	1996/09/09
*		Version 2.1
*		Original development and testing
*
* $Log: MHDc2fstr.c,v $
* Revision 5.1  2005/04/04 16:37:23  vlin
* cast to int before assign to i.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
*
*!Team-unique header:
* This software is developed by the MODIS Science Data SupportTeam for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
*!References and Credits
* 		Qi Huang	1996/09/09
*		Version 2.1
*		Original development and testing
*
*!Design Notes
*
!END********************************************************************
*/
{
    int         i;

    i=(int)strlen(str);
    for (; i < len; i++)
        str[i] = ' ';

    return MAPIOK;
}
