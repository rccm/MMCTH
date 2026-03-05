#include "mapic.h"

void ncrmfh ( intf  modfil[MODFILLEN], intf *ret)

/*   !C
 *
 *   !Purpose:	Provides a wrapping  function interfacing between C and FORTRAN.
 *              This C function is only called by the FORTRAN function RMFH.
 *              This function is a M-API internal routine.
 *
 *   !Description:  Function RMFH is part of a larger software system called
 *              the MODIS Applications Programming Interface (M-API) utility.
 *              The M-API utility consists of subroutines
 *              which allow MODIS Science Team-supplied software to read and
 *              write data, and metadata from/to HDF files. The functionality
 *              of the M-API is defined in the MODIS Application Program Interface
 *              (API) Specification
 * 
 *              Once  a  HDF-EOS file has been opened using an HDF-EOS open file 
 *              routine.  The user must call CMFH to create the M-API file id handle.
 *              The M-API routines can then be used to access the data object(s) in the
 *              opened  file.  Once the user is finished accessing the objects they must
 *              call RMFH to release this file id handle  created  by calling CMFH.  RMFH
 *              writes the objects accessed using M-API routines to a file, ends access to
 *              them, and frees the memory allocated.  Finally before exiting the user must
 *              call an HDF-EOS close file routine to close the HDF-EOS file.
 *
 *              crmfh is a C function which is callable from FORTRAN. This function will call
 *              releaseMAPIfilehandle to release the modfil array created using
 *              createMAPIfilehandle.
 *              
 *              In M-API, crmfh is a low-level routine called only by RMFH. 
 *              In order to be callable from the FORTRAN on different platforms using
 *              function name crmfh, this function is called ncrmfh in the actual C code.
 *              ncrmfh is redefined in mapic.h according to the compiler's FORTRAN naming
 *              conventions conversion for each platform, so that the object name of ncrmfh
 *              will always be the object name of a FORTRAN function named crmfh.
 *
 *   !Input  parameters:
 *        modfil IN/OUT: FORTRAN MODIS file array. If all successful, the elements of the
 *              array will be filled with zero after closing the MODIS file.
 *
 *   !Output  parameters:	
 *
 *        ret OUT: 	if function is successful, ret is set to MAPIOK, otherwise MFAIL.
 *
 *   !Returns:	None.
 *
 *   !Revision History:
 *     $Log: crmfh.c,v $
 *     Revision 1.2  1998/11/12 20:14:36  solanki
 *     Added '!' where needed.
 *
 * Revision 1.4  1997/12/23  22:59:01  fshaw
 * correctd typos
 *
 * Revision 1.3  1997/12/22  21:35:47  fshaw
 * made 12/17 W/T changes and corrections.
 *	
 *
 *   !Team_unique header:
 *     This software is developed by the MODIS Science Data Support Team for the
 *     National Aeronautics and Space Administration, Goddard Space Flight Center,
 *     under contract NAS5-32373.
 *
 *   !References and Credits:
 *     HDF portions developed at the National Center for Supercomputing
 *     Applications at the University of Illinois at Urbana-Champaign.
 *
 *   !Design Notes:
 *     None.
 *     
 *   !END
 ****************************************************************************************/
 {
	MODFILE *mfile;
        int i;

        *ret =  MFAIL;

	memcpy( &mfile, &modfil[P_ADDR], sizeof(MODFILE *));
		
        *ret = releaseMAPIfilehandle(&mfile);

        if ( *ret == MAPIOK){
            for ( i = 0; i < MODFILLEN; i++)
              modfil[i] = 0L;
        }

        return;
}
