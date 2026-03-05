#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "mapic.h"
#define ATTRIBUTE "Attr0.0"
#include "local_nc.h"

int32 MVSfind(int32 fid,  char const *vsname)

/*****************************************************************************
*!C
*
*!Description: 
*      To replace the HDF VSfind for obtaining the reference id of a Vdata.
*            Function MVSfind is part of a larger software system called the 
*            MODIS Applications Programming Interface (API) Utility,abbreviated
*            M-API. The M-API Utility consists of subroutines which allow MODIS
*            Science Team-supplied software to read and write data and metadata
*            from/to HDF files. The functionality of the M-API is defined in 
*            the MODIS Application Program Interface (API) Specification.
*
*            MVSfind is an MAPI internal function for obtaining the reference 
*            id of a Vdata set. MVSfind is used to replace HDF VSfind which can
*            not distinguish the user-defined Vdata set from attributes.
*
*!Input Parameters:
*   fid	     HDF file id returned from Hopen
*   vsname   Name of the Vdata
*
*!Output Parameters:   none
*
* Return values:
*            Vdata reference number   if successful
*               FAIL                  otherwise
*
* Externally defined:
*            VSNAMELENMAX       (hdf.h)
*            VGNAMELENMAX       (hdf.h)
*            VSgetid            (vproto.h)
*            VSattach           (vproto.h)
*            VSgetclass         (vproto.h)
*            VSgetname          (vproto.h)
*            VSdetach           (vproto.h)
*            VSQueryref         (vg.h)
*            DIM_VALS01         (local_nc.h)
*
*!Revision History:
* $Log: MVSfind.c,v $
* Revision 5.1  2005/04/04 16:42:59  vlin
* Replaced "DIM_VALS" with "DIM_VALS01".
* constant safe for pointer arguments.
*
*            Xia W. Chang,  12/29/95  , xia@modis-xl.gsfc.nasa.gov
*            Version 2 original development.
*	
*!Team-unique Header:
*
*           This software is developed by the MODIS Science Data Support 
*           Team for the National Aeronautics and Space Administration, 
*           Goddard Space Flight Center, under contract NAS5-32373.
*
*           HDF portions developed at the National Center for Supercomputing
*           Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!END
**************************************************************************/

{
  int32 vsref = -1;
  int32 vsid;
  int32 ret_ref = FAIL;
  char  classname[VSNAMELENMAX+1];
  char  vdataname[VGNAMELENMAX+1];

  /* Loop calling VSgetid until the return value is not FAIL  */
  while((vsref = VSgetid(fid, vsref)) != FAIL) {
     /* get the Vdata access id   */
     if ((vsid = VSattach(fid,vsref,"r")) == FAIL) 
	return(FAIL);

     /* get the classname of the attached Vdata  */
     VSgetclass(vsid, classname);
     if ((strcmp(classname,ATTRIBUTE)!=0) && (strcmp(classname,DIM_VALS01)!=0))
     {
        VSgetname(vsid, vdataname);
        /* if vdataname is same as vsname, get the reference id  */
        if (strcmp(vdataname, vsname) == 0 )
        {
           ret_ref = VSQueryref(vsid);
           VSdetach(vsid);
           return(ret_ref);
        }
        VSdetach(vsid);
     }
    }      /* end of while  */
  return(FAIL);
}
