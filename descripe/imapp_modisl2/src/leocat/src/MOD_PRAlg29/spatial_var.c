/*
***********************************************************************
!F77

!DESCRIPTION:
   Routine for finding the number of pixels surrounding the center
   pixel in a context (nlcntx*necntx pixels) that satisfy a difference
   criterion.

!Input Parameters:
   None

   Inputs stored in structure variable 'pxin' of type "pixel_in" defined
   in pixel.h.

!Output Parameters:
   None

   Outputs stored in structure variable 'rg_var' of type "regional_var"
   defined in pixel.h

!Revision History:
   06/04 Collection 5  R. Frey
   Updated argument list.
   11/07 Converted to C   R. Frey

!Team-unique Header:

!References and Credits:
 See Cloud Mask ATBD  
 
!END
***********************************************************************/


/* Includes */

#include <stdio.h>
#include <math.h>
#include "pixel.h"
#include "thresholds.h"
#include "mask_processing_constants.h"


int spatial_var()


{

/*  Declarations */

	int i;
    int ipt;
    int result;

/*  Initializations */

//  printf("Executing spatial_var \n");
    
    ipt = 0;

/**********************************************************************/


//  Compare surrounding differences to threshold.
    for ( i=0; i<num_diffs_surr; ++i) {

      if ( rintf(rg_var.diff[i]) != rintf(bad_data) ) {
        if ( rg_var.diff[i] <= dovar11[0] ) {
          ipt++;
        }
      }

    }
    rg_var.num_small_diffs = ipt;

//  If all surrounding pixel differences were less than the threshold
//  value, the scene is declared to be uniform.
    if (ipt == num_diffs_surr) {
      result = 1;
    }
    else {
      result = 0;
    }

//  printf("ipt, result %d %d \n", ipt, result);
    
    return result;

}
