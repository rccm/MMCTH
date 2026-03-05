/*
************************************************************************
!F77 

!Description:

   Routine which checks for strong temperature inversions. When the 6.7
   micron BT is greater than the 11 micron BT by a predetermined threshold,
   then there is very likely a strong temperature inversion and therefore
   clear skies.  The same applies to the 13.3-11 micron and 7.2-11 micron
   BTDs.

!Input parameters:
   none

   Inputs stored in structure variable 'pxin' of type "pixel_in" defined
   in pixel.h. 

!Output Parameters:
   conf        modified clear sky confidence for pixel

   Returned through function call

!Revision History:
    Converted to C         11/07      R. Frey

!Team-unique Header:

!References and Credits:
   See Cloud Mask ATBD.

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>
#include "pixel.h"
#include "thresholds.h"
#include "mask_processing_constants.h"


float NightSnow_Inv_check()


{

/***********************************************************************/
	
//  Declarations
	
    extern void set_bit(int, unsigned char[]);
    
    float m27, m28, m31, m33;
    float conf;
    
/**********************************************************************/

//  Initializations

//  printf("Executing NightSnow_Inv_check \n");
    conf = pxout.intermediate_conf;
    
    m27 = pxin.rad[26];
    m28 = pxin.rad[27];
    m31 = pxin.rad[30];
    m33 = pxin.rad[32];
    
/**********************************************************************/    
    
    if(pxin.polar) {

//    Tests for polar regions.    	
    	
      if (rintf(m31) != rintf(bad_data) && rintf(m27) != rintf(bad_data)) {
    	(void) set_bit(26, pxout.qabits);
        if ( (m27 - m31) > pn65_11[0]) { 
          conf = 1.0;
          (void) set_bit(26, pxout.testbits);
        }
      }

      if(conf < 1.0) {
        if (rintf(m31) != rintf(bad_data) && rintf(m33) != rintf(bad_data)) {
          (void) set_bit(26, pxout.qabits);
          if ( (m33 - m31) > pn13_11[0]) { 
            conf = 1.0;
            (void) set_bit(26, pxout.testbits);
          }
        }
      }

      if(conf < 1.0) {
        if (rintf(m31) != rintf(bad_data) && rintf(m28) != rintf(bad_data)) {
          (void) set_bit(26, pxout.qabits);
          if ( (m28 - m31) > pn7_11[0]) { 
            conf = 1.0;
            (void) set_bit(26, pxout.testbits);
          }
        }
      }

    } 
    else {
    	
//    Tests for non-polar areas.    	
    	
      if (rintf(m31) != rintf(bad_data) && rintf(m27) != rintf(bad_data)) {
    	(void) set_bit(26, pxout.qabits);
    	if ( (m27 - m31) > n65_11[0]) { 
    	  conf = 1.0;
    	  (void) set_bit(26, pxout.testbits);
    	}
      }	
    	
    }
    
//  printf("Inversion Tests: %d %f %f %f %f %f %f %f %f %f \n", pxin.polar,
//    	m31, m27, m28, m33, pn65_11[0], pn13_11[0], pn7_11[0], n65_11[0], conf);
    
/**********************************************************************/    
    
    return conf;
    
/**********************************************************************/    
    
}
