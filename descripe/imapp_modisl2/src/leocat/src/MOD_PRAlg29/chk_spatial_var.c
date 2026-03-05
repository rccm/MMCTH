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
   conf    Modified confidence of clear sky value 

   Returned through function call

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


float chk_spatial_var()


{

/*  Declarations */

    extern void set_bit(int, unsigned char[]);
    extern void get_regdif(int);
    extern int spatial_var();

    int band;
    int var;
    float conf;

/**********************************************************************/

//  Initializations.
    
//  printf("Executing chk_spatial_var at %d %d \n", pxin.scan, pxin.elem);
    conf = pxout.intermediate_conf;
    
/**********************************************************************/   
    
//  Get brightness temperature differences between pixel of interest
//  and the 8 surrounding it.

    band = 31;
    (void) get_regdif(band);

//  Check variation in the region.
    var = spatial_var();
        
//  Set bit that indicates this test was applied.
    (void) set_bit(25, pxout.qabits);

//  Increase confidence of clear sky if spatial variability test showed
//  uniform conditions. Set bit #25 if confidence is increased.
    
    if ( (var == 1) && (pxout.init_conf > 0.66) ) {

      conf = 0.96;
      (void) set_bit(25, pxout.testbits);
    
    } 
    else if ( (var == 1) && (pxout.init_conf <= 0.66) ) {

      conf = 0.67;
      (void) set_bit(25, pxout.testbits);

    }

//  printf("var, conf: %d %f %f \n", var, pxout.init_conf, conf);
    
    return conf;

}
