#include <stdio.h>
#include <string.h>
#include "mapic.h"

short int Vdatattr_name(char const *attribute, int32 vdata_ref, 
		      char *attribute_name)

/*
!C*************************************************************************
*
*Routine:     Vdatattr_name
*
*!Description: Subroutine Vdatattr_name is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
*	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      subroutines which allow MODIS Science Team-supplied software
*	      to read from/write to HDF format files. 
*
*             Vdattattr_name creates an attribute name that is
*             associated with a particular Vdata.  This allows the
*             global metadata facility to store Vdata metadata.  The
*             Vdata tag and the Vdata's ref-id are appended to the
*             attribute name to make a Vdata-specific metadata name.
*             If the attribute name is very long, these two integers
*             will overwrite the last few characters of the attribute
*             so that the langth of attribute_name remains within the
*             VSNAMELENMAX character limit for HDF global attribute
*             names.
*
*
*
* !Input Parameters:
*            attribute: Standard attribute name of a Vdata metadata.
*
*            vdata_ref: Reference number of the Vdata to associate
*                       the metadata with.
*
* !Output Parameters:
*            attribute_name: Name to assign the attribute in the HDF
*                            file so that it may be associated with a 
*                            particular Vdata.
*
*Externals:
*             MFAIL              "mapi.h"
*             NULLstr            "mapic.h"
*             sprintf            <stdio.h>
*             DFTAG_VH           "hdf.h"
*             strcpy             <string.h>
*             strlen             <string.h>
*             VSNAMELENMAX       "hdf.h"
*             strcat             <string.h>
*             MAPIOK             "mapi.h"
*
* !Revision History:
* $Log: Vdatattr_name.c,v $
* Revision 5.1  2005/04/04 17:58:32  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.1  1996/09/16  18:53:27  qhuang
 * Initial revision
 *
*
*       1.    Paul Fisher               20 July 95
*             Research and Data Systems Corporation
*             SAIC/GSC MODIS Support Office
*
*	      Original Development / Testing
*
* !Team-unique Header:
*             This software is developed by the MODIS Science Data Support
*	      Team for the National Aeronautics and Space Administration,
*	      Goddard Space Flight Center, under contract NAS5-32373.
*
* !References and Credits:
*             HDF Portions developed at the National Center for Supercomputing
*             Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !Design Notes:
*
!END**************************************************************************
*/
{
  short int output_status=MFAIL;           /* return status, set to MFAIL */
  char vdref[VSNAMELENMAX]={0};            /* array to hold ascii
					      representation of vdata_ref */
  char dftag[VSNAMELENMAX]={0};            /* array to hold ascii
					      representation of DFTAG_VH */
  register int i=0, k=0;                   /* generic counters: 
					      k for placed string
					      i for the receiving string*/


  /* input checking */
  if (NULLstr(attribute)||attribute_name==NULL) return (output_status);

  /* copy variables into strings */
  sprintf(dftag, "%d", (int)DFTAG_VH);
  sprintf(vdref, "%ld", (long int)vdata_ref);

  /* copy the attribute string into attribute_name */
  strcpy(attribute_name, attribute);

  if((strlen(attribute)+strlen(vdref)+strlen(dftag))>(VSNAMELENMAX-1)){
    /* place dftag and vdrref at the end of attribute_name */
    i=(VSNAMELENMAX-(int)(strlen(dftag)+strlen(vdref))-1);
    for(k=0;k<(int)strlen(dftag);k++){
      attribute_name[i]=dftag[k];
      i++;
    }
    for(k=0;k<(int)strlen(vdref);k++){
      attribute_name[i]=vdref[k];
      i++;
    }
    /* Place null string */
    attribute_name[VSNAMELENMAX]='\0';
  }

  else{
    strcat(attribute_name, dftag);
    strcat(attribute_name, vdref);
  }

  return (MAPIOK);
}

