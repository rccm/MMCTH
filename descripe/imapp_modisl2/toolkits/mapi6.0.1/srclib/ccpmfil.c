#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"
#include <ctype.h>

void nccpmfil(intf modfil[MODFILLEN], _fcd mdHandles, intf *lhandles,
         _fcd HDFattrnms, intf *lattrs, intf *NumHandles, intf *ret)
/*
!C
  **********************************************************************
* 
* Purpose: A wrapping function interfacing between C and FORTRAN for 
*          completeMODISfile. This C function is only called by FORTRAN
*          function CPMFIL. This function is a M-API internal routine.
*
* 
*!Description:
*              Function ccpmfil is part of a larger software system
*              called the MODIS Applications Programming Interface
*              (API) Utility, abbreviated M-API. The M-API Utility
*              consists of subroutines which allow MODIS Science
*              Team-supplied software to read and write data and
*              metadata from/to HDF files. The functionality of the
*              M-API is defined in the MODIS Application Program
*              Interface (API) Specification.
*
*              ccpmfil is a C function which is callable from FORTRAN.
*              This function will call completeMODISfile to form an array
*              structure. In M-API, ccpmfil is low-level routine which is
*              called only by CPMFIL. 
*
*              In order to be callable from the FORTRAN in different
*              platforms using function name ccpmfil, this function is
*              called nccpmfil in the actual C code. nccpmfil is redefined
*              in mapic.h according to compiler's FORTRAN naming
*              conventions/conversion of each platform, so that the object
*              name of nccpmfil will always the object name of a FORTRAN
*              function named ccpmfil.
* 
*!Input Parameters:
*  modfil: FORTRAN MODIS file array. If all successful, 
* 	   the elements of the array will be filled with 
* 	   zero after closing the MODIS file.
*
*  mdHandles: A character array with size of
*             [PGSd_MET_NUM_OF_GROUPS][PGSd_MET_GROUP_NAME_L-1], where
*             PGSd_MET_NUM_OF_GROUPS is 20 and
*             PGSd_MET_GROUP_NAME_L is 50.  This array is typedef-ed
*             as PGSt_MET_all_handles. Each row in the array stores
*             a handles to an internal ODL tree structure which will
*             be written out a an ECS PVL attribute.  Each handle,
*             which is a string, should be less than PGSd_MET_GROUP_NAME_L
*             characters and occupy, one row in the array.  Therefore, the
*             maxmum number of  handles should be 20.
*
*  HDFattrNames: A character array with size of 
*                [PGSd_MET_NUM_OF_GROUPS][MAX_ECS_NAME_L-1],
*                where PGSd_MET_NUM_OF_GROUPS is 20 and MAX_ECS_NAME_L
*                is 256.  This array is typedef-ed as ECSattr_names_for
*                all_handles.  Each row in this array is a character
*                string used as a global attribute name for storing
*                an ECS PVL text block which has a handle in the
*                corresponding row in mdHandles array. Each name, which
*                is a string, should be less that MAX_ECS_NAME_L
*                characters and occupies one row in the array.
*
*  NumHandles:  Specifie the number of actual handles contained in
*              mdHandles.  This may be set from 0 to PGSd_MET_NUM_OF_GROUPS.
*
*!Output Parameters:None.
* 
* Returns:	None.
* 
* External References:
*           ECSattr_names_for_all handles   (mapi.h)
*           PGSt_MET_all_handles            (mapi.h)
*           PGSt_MET_GROUP_NAME_L           (mapi.h)
*           MAX_ECS_NAME_L                  (mapi.h)
*           MAPIOK	                    (mapi.h)
*           completeMODISfile               (mapi.h)
*           PGSd_MET_NUM_OF_GROUPS          (PGS_MET.h)
*           _fcdtocp                        (hdfi.h)
*           memcpy                          (ctype.h)
*           MODFILLEN                       (mapic.h)
*
*!Revision History:
* $Log: ccpmfil.c,v $
* Revision 6.1  2010/07/13 19:57:24  kuyper
* Removed call to function that is not a part of the C standard library;
*   replaced with equivalent code that's slightly more portable.
*
* Revision 1.2  1999/04/26 15:38:41  jayshree
* reorganising mapi RCS
*
 * Revision 1.11  1999/04/09  20:15:12  jayshree
 * added lhandles, lattrs
 * modified for mapi2.3.2 release
 *
 * Revision 1.9  1996/08/16  12:56:20  fshaw
 * ver 2.1
 * modified to use ring super structure
 *
 * Revision 1.8  1996/08/16  12:09:49  fshaw
 * version 2.0
 *
 * Revision 1.7  1996/08/13  19:17:35  fshaw
 * changed mfile to modfile
 *
 * Revision 1.6  1996/08/06  15:15:26  fshaw
 * vers 2.0
 *
* Revision 1.6  1996/08/06 15:15:26  fshaw
* vers 2.0
*
 * Revision 1.5  1996/08/01  15:20:02  fshaw
 * updated prologue
 *
 * Revision 1.4  1996/07/29  13:53:09  fshaw
 * updated to correct problems
 * with character array sizes
 *
 * Revision 2.0  1996/07/22  20:51:45  fshaw
 * changed while test to be less than
 *
 * Revision 1.3  1996/07/16  20:32:32  fshaw
 * added validation check for cNumhandles
 *
 * Revision 1.2  1996/07/16  18:41:05  fshaw
 * working version
 *
 * Revision 1.1  1996/06/14  19:21:57  fshaw
 * Updated prolog. 
 *
 * Revision 1.0  1996/04/29  18:10:34  fshaw
 * Initial revision
 *
* 	
*!Team-unique Header:
*
*!References and Credits:
* This software is developed by the MODIS Science Data Support Team for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
*!Design Notes:
*
*!END********************************************************************
*/

{
    char *cmdHandles, *cHDFattrnms, *cp;
    long int cNumHandles;
    int      i, j;
    PGSt_MET_all_handles Handles;
    ECSattr_names_for_all_handles  attrNames;
    MODFILE  *mfile;
        
    cNumHandles = *NumHandles;
        
    cmdHandles  = (char *)_fcdtocp(mdHandles);
    cHDFattrnms = (char *)_fcdtocp(HDFattrnms);
  
    memcpy( &mfile, &modfil[3], sizeof(MODFILE *));    

    if ((cNumHandles > 0) && (cNumHandles < PGSd_MET_NUM_OF_GROUPS )){

	for ( i=0; i<=cNumHandles; i++){

           cp  = cmdHandles + i  * (*lhandles);
           
           for ( j = *lhandles - 1;
               j>=0 && ((unsigned)cp[j] > 0x7FU || !isgraph(cp[j])); j--)
                  /*EMPTY */;
           if(j >=0) {
	       if(j > PGSd_MET_GROUP_NAME_L - 2)
	       {
	       j = PGSd_MET_GROUP_NAME_L - 2; }
	       Handles[i][0]='\0';
               memcpy(Handles[i], cp, (size_t)(j + 1));
               Handles[i][j + 1] = '\0';
           }
           else
              Handles[i][0] = '\0';
            
           cp  = cHDFattrnms  + i  * (*lattrs);
           
           for ( j = *lattrs - 1;
                 j>=0 && (!isascii(cp[j]) || !isgraph(cp[j])); j--)
                  /*EMPTY */;
           if(j >=0) {
	       if (j > MAX_ECS_NAME_L - 2)
	       {
	       j = MAX_ECS_NAME_L - 2; }
	       attrNames[i][0] = '\0';
               memcpy(attrNames[i], cp, (size_t)(j + 1));
               attrNames[i][j + 1] = '\0';
           }
           else
              attrNames[i][0] = '\0';
        }
    }

      *ret = completeMODISfile (&mfile, Handles, attrNames, cNumHandles);
      if ( *ret == MAPIOK) {
          for ( i = 0; i < MODFILLEN; i++){
              modfil[i] = 0L;
          }
      } 
      return;
   }
