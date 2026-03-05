#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "mapic.h"

int getMODISdimname(MODFILE *file, char const *arrayname, char const *groupname,
                  long int dimension, char *dimname)
/*
*************************************************************************
*!C 
*
* Purpose:     Reads the name (string) associated with a spe- 
*              cific dimension of an SDS(array) from a MODIS HDF file. 
*
*!Description: Function getMODISdimname is part of a larger software sys-
*              tem called the MODIS Applications Programming Interface
*              (API) Utility, abbreviated M-API.  The M-API Utility con_
*              sists of subroutines which allow MODIS Science Team-sup-
*              plied software to read and write data and metadata from/to
*              HDF files.  The functionality of the M-API is defined in 
*              the MODIS Application Program Interface (API) Specification.
*
*              getMODISdimname (GMDMIN) retrieves the name associated with
*              an array dimension by giving the array's name and dimension
*              number. If the name cannot be found, the routine 
*              will return MFAIL and the passed variable unchanged.

*             The groupname string provides the facility to select an array
*             structure placed in a particular HDF Vgroup data group.
*             Alternatively, the entire file will be searched for an array
*             structure named arrayname if groupname = NULL in C and a blank
*             string, " " in FORTRAN.
*
*
*!Input Parameters:
*       file         IN: Address of MODFILE structure that is used to 
*                    reference a MODIS HDF file containing the attribute.
*       arrayname    IN: ASCII string name of the array.  Provided macros 
*                    for accepted MODIS HDF file array names are listed in
*                    Appendix A, M-API-Supplied Constants.
*       groupname    IN: ASCII string name of the data group containing 
*                    the array structure to which the attribute is atta-
*                    ched.  If set to NULL, the entire file will be sear-
*                    ched for the array structure named arrayname.
*       dimension    IN: The dimension number from which the attribute is
*                    retrieved(0-based).
*
*!Output Parameters: 
*       dimname      OUT: ASCII string name of the dimension.  Provided 
*                    macros for accepted MODIS HDF file attribute names 
*                    are listed in Appendix A, M-API-Supplied Constants.
*                   function failure occurs, the value will be set to zero.
*
* Returns:     MAPIOK if successful, MFAIL if an error occurs.
*
* External references:
*                      MAPIOK             (mapi.h)
*                      MFAIL              (mapi.h)
*                      NULLstr            (mapic.h)
*                      MAX_NC_NAME        (netcdf.h)
*                      SDgetdimid         (mfhdf.h)
*                      MODFILE            (mapi.h)
*                      PGS_MAX_MSGBUF_SIZE (mapic.h)
*                      DATAID              (mapi.h)
*                      SDSINFO             (mapic.h)
*                      NULLMODFIL          (mapic.h)
*                      getMODISarrayid.c   (mapi.h)
*                      SDgetdimid          (mfhdf.h)
*                      FAIL                (hdf.h)
*                      SDdiminfo           (mfhdf.h)
*
*!Revision History:
* Jeffrey J. Blanchette 1996/1/1 jjb@modis-xl.gsfc.nasa.gov
*               PROTOTYPE
*$Log: getMODISdimname.c,v $
*Revision 5.1  2005/04/04 18:49:11  vlin
*constant safe for pointer arguments.
*
*Revision 1.1  1998/02/06 22:26:06  fshaw
*Initial revision
*
 * Revision 1.3  1996/08/06  12:16:06  fshaw
 * corrected E & W from pr:qa
 *
 * Revision 1.2  1996/08/06  12:08:32  fshaw
 * version 2.1
 *
*
*
*!Team-unique Header:
*
*              Portions developed at the National Center for Supercompu-
*              ting Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!References and Credits:
*
*              This software is developed by the MODIS Science Data Sup-
*              port Team for the National Aeronautics and Space Admini-
*              stration, Goddard Space Flight Center, under contract 
*              NAS5-32373.
*
!END**************************************************************************
*/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char *funcname="getMODISdimname"; 
  int  status;             /* variable to hold the returned value */
  int32  sds_id;
  int32  nattrs; 
  int32  dim_id;
  int32  nnumber_type;
  int32  count;
  DATAID *did;
  SDSINFO *sinfo;

  status = MFAIL;        /* error */

  /* input checks */
  if ( NULLstr(arrayname) ){
    sprintf(buff,"ERROR: getMODISdimname unable to read the name of a \n"
                 "\t dimension without the name of array it is \n"
                 "\t associated with. \n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }


  if ( NULLMODFIL(file) ) {
    sprintf(buff,"ERROR: getMODISdimname unable to read the name of a \n"
                 "\t dimension in the %.*s array with an invalid \n"
                 "\t MODIS file structure input.\n",
                  MAX_NC_NAME,arrayname);
    return(MFAIL);
  }

  if ( dimname == NULL ) {
    sprintf(buff,"ERROR: getMODISdimname unable to read the name of a \n"
                 "\t %.*s array's dimension name without an \n"
                 "\t output character string. \n",
                  MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  did = getMODISarrayid( file, arrayname, groupname);

  if ( did  == NULL) {
    sprintf(buff,"ERROR: getMODISdimname detected MFAIL from M-API \n"
                 "\t internal function getMODISarrayid while \n"
                 "\t attempting to obtain the name of dimension \n"
                 "\t %ld in the %.*s array.\n",
                  dimension, MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
    return (MFAIL);
  }

  sinfo = (SDSINFO *)did->info;
  sds_id = (int32)did->id;

  if ( (dimension < 0) || (dimension > sinfo->rank-1 ) ) {
         sprintf(buff,"ERROR: getMODISdimname unable to read the dimension\n"
                      "\t name of the non-existent dimension %ld of\n"
                      "\t %.*s array.\n",
                      dimension, MAX_NC_NAME,arrayname);
         MAPIERR(buff,funcname);
  }else if ( (dim_id = SDgetdimid(sds_id,(int32)dimension)) == FAIL ){
         sprintf(buff,"ERROR: getMODISdimname detected FAIL from HDF \n"
                      "\t procedure SDgetdimid attempting to read the name \n"
                      "\t of an %.*s array's dimension.\n",
                      MAX_NC_NAME,arrayname);
         MAPIERR(buff,funcname);
  }else if ( SDdiminfo(dim_id,dimname,&count,&nnumber_type,&nattrs) == FAIL ){
         sprintf(buff,"ERROR: getMODISdimname detected FAIL from HDF\n"
                      "\t procedure SDdiminfo attempting to read the \n"
                      "\t name of an %.*s array's dimension.\n",
                      MAX_NC_NAME,arrayname);
         MAPIERR(buff,funcname);
  }else
    status = MAPIOK;

  return(status);
}
