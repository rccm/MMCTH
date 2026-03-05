#include "mapic.h"

void nccmfh(intf *SWfid, intf modfil[MODFILLEN], intf *ret)

/****************************************************************
 * !C
 *
 * !Purpose:	Provides a function interfacing between C and FORTRAN.
 *              This C function is only called by the FORTRAN function
 *              CMFH. This function is a M-API internal routine.
 *
 * !Description: Function ccmfh is part of a larger software system
 *               called the MODIS Applications Programming Interface (MAPI)
 *               utility. The M-API utility consists of
 *               subroutines which allow MODIS Science Team-supplied software
 *               to read  and write data and metadata from/to HDF files.  The
 *               functionality of the M-API is defined in the MODIS Application
 *               Program Interface (API) Specification.
 *
 *               ccmfh is a C function which is callable from FORTRAN. This
 *               function will call createMAPIfilehandle to open/create a MODIS-HDF
 *               file. In M-API, ccmfh is low-level routine which is called only by CMFH. 
 *
 *               In order to be callable from the FORTRAN in different platforms using
 *               function name ccmfh, this function is called nccmfh in the actual C code.
 *               nccmfh is redefined in mapic.h according to the compiler's FORTRAN naming
 *               conventions/conversion for each platform, so that the object name of
 *               nccmfh will always be the object name of a FORTRAN function named ccmfh.
 *
 * !Input parameters:
 *      SWfid    IN:   FORTRAN integer for the HDF-EOS file id handle. 
 *
 * !Output parameters:
 *
 *      modfil   OUT:  Array that is used to reference the file in all other MODIS API
 *                     routines in FORTRAN.  The array will return all zeroes if an error occurs.
 *
 *      ret      OUT:  if function is successful, ret is set to MAPIOK,	otherwise MFAIL.
 *
 * !Returns:   none.
 *
 * !Revision History:
 *      $Log: ccmfh.c,v $
 *      Revision 1.1  1998/02/06 22:26:06  fshaw
 *      Initial revision
 *
 * Revision 1.5  1997/12/23  14:53:27  fshaw
 * corrected typo's found in 12/17 W/T
 *
 * Revision 1.4  1997/12/22  21:47:03  fshaw
 * corrected defects found in 12/17 W/T.
 *
 * Revision 1.3  1997/12/09  19:57:53  fshaw
 * fixed defect in memcpy and removed magic numbers
 *
 * Revision 1.1  1997/11/24  15:59:32  fshaw
 * Initial revision
 *
 *      Frederick J. Shaw  (fshaw@ltpmail.gsfc.nasa.gov)
 *      General Sciences Corp.
 *
 * !Team-unique header:
 *      This software is developed by the MODIS Science Data Support Team for the National
 *      Aeronautics and Space Administration, Goddard Space Flight Center, under contract
 *      NAS5-32373.
 *
 * !References and Credits:
 *    HDF portions developed at the National Center for Supercomputing
 *    Applications at the University of Illinois at Urbana-Champaign.
 *
 * !Design Notes:
 *    none.
 * 
 * !END
 ****************************************************************************************/
 {
  MODFILE  *mfile;
  int       i;
  int32     swfid;

  for ( i = 0; i < MODFILLEN; i++) 
              modfil[i] = 0L;

  swfid = (int32)(*SWfid);
  mfile = createMAPIfilehandle(swfid);

  if( mfile == NULL)
     *ret = MFAIL;
  else{
     *ret = MAPIOK;
     /* Set the elements in modfil according to values in mfile. */
     modfil[P_SDID] = (intf)mfile->sd_id;
     modfil[P_HDFID] = (intf)mfile->hdf_id;
     modfil[P_ACCESS] = (intf)mfile->access;
     memcpy( &modfil[P_ADDR], &mfile, sizeof(MODFILE *));
  };

  return;
}

