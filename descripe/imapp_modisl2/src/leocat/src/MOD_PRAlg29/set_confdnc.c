/*
************************************************************************
!F77 
 
!Description:
  Routine for setting output "bit" flags according
  to final confidence of clear sky.

!Input Parameters:
  None
  pxin.confdnc (final confidence of clear sky) is stored in pixel.h.

!Output Parameters:
  None
  pxout.testbits (6 byte array containing cloud mask bit flags) is
  stored in pixel.h.

!Revision History:
  R. Frey    10-25-07  Converted to C

!Team-Unique Header:

!References and Credits:

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include "pixel.h"
#include "mask_processing_constants.h"
      
void set_confdnc()
      
{      

   extern void set_bit(int, unsigned char[]);
   
	
   if(pxout.confidence > conf_clr_thresh) {
     (void) set_bit(1, pxout.testbits);
	 (void) set_bit(2, pxout.testbits);
   }
   else if(pxout.confidence > prob_clr_thresh) {
	 (void) set_bit(2, pxout.testbits);
   }	 
   else if(pxout.confidence > prob_cld_thresh) {
	 (void) set_bit(1, pxout.testbits);
   }


}
