/*
************************************************************************
!F77 

!Description:
  Sets bit flags indicating processing path.

!Input integer parameters found in pixel.h (1=true, 0=false):
  pxin.water         water background
  pxin.land          land background
  pxin.day           sza < 85
  pxin.ice           ice background
  pxin.snow          snow background
  pxin.snglnt        sun glint
  pxin.coast         coast background
  pxin.desert        desert background

!Output Parameters:
 testbits      6 byte array containing cloud mask bit flags

!Revision History:
  R. Frey      10-25-07        Converted to C

!Team-Unique Header:

!References and Credits:

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include "pixel.h"
#include "mask_processing_constants.h"
      
void set_proc_path()
      
{      

   extern void set_bit(int, unsigned char[]);


// Set snow/ice bit.
   if( !pxin.snow && !pxin.ice) {
     (void) set_bit(5, pxout.testbits);
   }

// Set day/night flag.
   if (pxin.day) {
     (void) set_bit(3, pxout.testbits);	   
   }

// Set sunglint flag.
   if(!pxin.sunglint || !pxin.water) {
     (void) set_bit(4, pxout.testbits);
   }

// Set coast, desert, or land processing path flags.  Default is water (00).
   if (pxin.coast) {
	 (void) set_bit(6, pxout.testbits);
   }
   else if (pxin.desert) {
	 (void) set_bit(7, pxout.testbits);
   }
   else if (pxin.land) {
     (void) set_bit(6, pxout.testbits);
	 (void) set_bit(7, pxout.testbits);
   }

// Set NISE snow/ice flag.  Set corresponding QA flag for display purposes.
   (void) set_bit(10, pxout.qabits);
   if(!pxin.map_snow) 
     (void) set_bit(10, pxout.testbits);

   
}
