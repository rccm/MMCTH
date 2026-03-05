#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"

void  ncgmdmin(intf modfil[], _fcd arrnm, intf *lar, _fcd grpnm, intf *lgr,  intf *dim,  _fcd attrnm, intf *lat, _fcd dtype, intf *ldt, intf *nms,  void *value, intf *ret)

/****************************************************************************************
!C 
* 
* !Purpose: A wrapping function interfacing between C and FORTRAN for getMODISdiminfo.
*          This C function is only called by FORTRAN function GMDMIN. This function is
*          a M-API internal routine.
*
*!Description: Function cgmdmin is part of a larger software system called the MODIS
*              Applications Programming Interface (API) Utility, abbreviated M-API. The
*              M-API Utility consists of subroutines which allow MODIS Science
*              Team-supplied software to read and write data and metadata from/to HDF
*              files. The functionality of the M-API is defined in the MODIS Application
*              Program Interface (API) Specification.
*
*              cgmdmin is a wrapper which is callable from FORTRAN. This function will
*              call getMODISdiminfo to read a dimensional attribute. In M-API, cgmdmin
*              is low-level routine which is called only by PMAR. 
*
*              In order to be callable from the FORTRAN in different platforms using
*              function name cgmdmin, this function is called ncgmdmin in the actual
*              C code. ncgmdmin is redefined in mapic.h according to compiler's FORTRAN
*              naming conventions/conversion of each platform, so that the object name
*              of ncgmdmin will always be the object name of a FORTRAN function named
*              cgmdmin.
*
*!Input Parameters:
*
*              modfil IN: FORTRAN integer array that is used to reference the MODIS-HDF
*                         file.
*              arrnm  IN: FORTRAN character string that is the name of the array.
*              lar    IN: FORTRAN integer address of the memory size of arrnm.
*              grpnm  IN: FORTRAN character string which is the  name of the data group
*                         to which the array (SDS) belongs.
*              lgr    IN: FORTRAN integer address of the memory size of grpnm.
*              dim    IN: The dimension from which the attribute will be read (0-based).
*              attrnm IN: Name of the attribute.
*              lat    IN: FORTRAN address of the string length of attrnm
*              dtype  IN/OUT:  Data Type of the value.
*              ldt    IN: FORTRAN address of the string length of dtype
*              nms    IN/OUT: Number of attribute values in value.
*
*!Output Parameters:
*              value  IN: The data to store in the attribute.
*              ret    OUT: FORTRAN integer address of the status(MFAIL, MAPIOK) 
*
*!Returns: none
*
*!External References:
*		MODFILLEN                         (mapic.h)
*               MODFILE                           (mapi.h)
*               HDf2cstring                       (hproto.h)
*               P_ADDR                            (mapic.h)
*               DATATYPELENMAX                    (mapic.h)
*               FDATATYPELENMAX                   (mapic.h)
*               _fcdtocp                          (hdfi.h)
*               MTYPEf2c                          (mapic.h)
*               MTYPEc2f                          (mapic.h)
*               HDfreespace                       (hproto.h)
*               VOIDP                             (hdfi.h)
*               TXT                               (mapi.h) 
*               getmodisdiminfo                   (mapic.h)
*
*!Revision History:
*     $LOG:$
*
*!Team-unique Header:
*      This software is developed by the MODIS Science Data Support Team for the
*      National Aeronautics and Space Administration, Goddard Space Flight Center,
*      under contract NAS5-32373.
*
*!References and Credits:
*
*!Desing Notes
*
*
!END
*/
{

   MODFILE  *mfile;
   char *carrnm, *cgrpnm, *cattrnm, *fcdtype, *cfdtype;
   DATAID *did;
   char cdtype[DATATYPELENMAX], fdtype[FDATATYPELENMAX + 1];
   long int i;
   long int rank, cdim, n_elements;
   long int c_length = DATATYPELENMAX, f_length = FDATATYPELENMAX;

   n_elements = *nms;

   /*
   Convert the FORTRAN character strings arrnm, grpnm, attrnm, dtype to C character
   strings carrnm, cgrpnm, cattrnm, cfdtype by using HDf2cstring. 
   */

   /* Convert FORTRAN strings to C strings using HDf2cstring. */
   cgrpnm = HDf2cstring(grpnm, (intn)*lgr);
   carrnm = HDf2cstring(arrnm, (intn)*lar);
   cattrnm = HDf2cstring(attrnm, (intn)*lat);
   cfdtype = HDf2cstring(dtype, (intn)*ldt);

   /* Set file by memcpy. */
   memcpy(&mfile, &modfil[P_ADDR], sizeof(MODFILE *));

   did = getMODISarrayid( mfile, carrnm, cgrpnm);

   if( did == NULL ) {
      if (carrnm) HDfreespace((VOIDP)carrnm);
      if (cgrpnm) HDfreespace((VOIDP)cgrpnm);
      if (carrnm) HDfreespace((VOIDP)carrnm);
      if (cfdtype) HDfreespace((VOIDP)cfdtype);
      *ret = MFAIL;
      *nms = 0;
      return;
   }

   /* convert FORTRAN dimension to C dimension */

   rank = (long int)((SDSINFO *)did->info)->rank;
   cdim = rank - (*dim) - 1;

   /* convert FORTRAN data type to C data type */
   cdtype[0] = '\0';

   /* Call to MTYPEf2c to convert cfdtype to cdtype */
   MTYPEf2c(cfdtype, cdtype, &c_length);

   *ret = getMODISdiminfo( mfile, carrnm, cgrpnm, cdim, cattrnm, cdtype, &n_elements, value);

   /* make the output for nms and attribute value */
   if (n_elements != 0) { 

        /*   (Note: make FORTRAN dtype)
              Convert dtype to c pointer fcdtype by _fcdtocp ) */
      fcdtype = (char *)_fcdtocp(dtype);

      if (cdtype[0] != '\0'){
          
        /* Convert C type cdtype to FORTRAN type fdtype by using MTYPEc2f.  */
          MTYPEc2f(cdtype, fdtype, &f_length);        
          
        /* Set the memory of fcdtype to blanks by memset.                     */
          memset(fcdtype, ' ', (size_t)*ldt);

       /* memcpy the content of cfdtype to fcdtype */
          memcpy(fcdtype, fdtype, strlen(fdtype));
      }
      if ( strcmp(cdtype, TXT) == 0){
          n_elements = n_elements - 1;
          if ( n_elements < *nms){
            for ( i =  n_elements; i< *nms; i++){
                *((char *)value + i) = ' ';     
            }
          }
      }
   }
   *nms = (intf)n_elements;

   /* HDfreespace carrnm, cgrpnm, cattrnm, cfdtype */

   if (cattrnm) HDfreespace((VOIDP)cattrnm);
   if (cfdtype) HDfreespace((VOIDP)cfdtype);
   if (carrnm) HDfreespace((VOIDP)carrnm);
   if (cgrpnm) HDfreespace((VOIDP)cgrpnm);

   return;
}

