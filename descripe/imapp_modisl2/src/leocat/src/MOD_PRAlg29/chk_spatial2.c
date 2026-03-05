/*
************************************************************************
!F77 

 !Description:

   Routine for checking spatial variability of ocean scenes over
   a 3X3 pixel region.

!Input parameters:
   none

   Inputs stored in structure variable 'pxin' of type "pixel_in" defined
   in pixel.h. 

!Output Parameters:
   none

   Outputs stored in structure variable 'rg_var' of type "regional_var"
   in pixel.h.

!Revision History:

   06/04 Collection 5  R. Frey:
   11/07 Converted to C  R. Frey

!Team-unique Header:

!References and Credits:
   See Cloud Mask ATBD-MOD-35.

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include "pixel.h"


void chk_spatial2()


{

/*  Declarations */

    extern void get_regdif(int);
    extern int spatial_var();

    int band;
    int var;

/**********************************************************************/

//  Get brightness temperature differences between pixel of interest
//  and the 8 surrounding it.

    band = 31;
    (void) get_regdif(band);

//  Check variation in the region.
    var = spatial_var();
        
    return;

}
