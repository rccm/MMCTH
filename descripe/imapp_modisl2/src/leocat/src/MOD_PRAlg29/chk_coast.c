/*
************************************************************************
!F77 

!Description:

   Routine which checks for extreme values of NDVI in coastal situations.
   Sets confidence to "confident clear" if very low or
   very high values are found (clear water or clear land, respectively).

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


float chk_coast()


{

/***********************************************************************/
	
//  Declarations
	
    extern int check_bits(int, unsigned char[]);
    extern void set_bit(int, unsigned char[]);
    
    float m01, m02;
    float ndvi;
    float conf;
    
    int bit_test;
    int qa_test;
    int irclr;
    
/**********************************************************************/

//  Initializations

//  printf("Executing chk_coast \n");
    conf = pxout.intermediate_conf;
    
    ndvi = 0.0;
    
    m01 = pxin.rad[0];
    m02 = pxin.rad[1];
    
/**********************************************************************/    
    
//  Determine the logical flag 'irclr' - true if ir cloud tests below
//  have all been passed. IR thin cirrus test makes final decision for
//  bit 18.
    
    irclr = 1;
    bit_test = check_bits(14, pxout.testbits);
    qa_test = check_bits(14, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    bit_test = check_bits(15, pxout.testbits);
    qa_test = check_bits(15, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    bit_test = check_bits(18, pxout.testbits);
    qa_test = check_bits(18, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;

    if(irclr) {
      if ( rintf(m02) != rintf(bad_data) &&
              rintf(m01) != rintf(bad_data) ) {

//      Check for very low or very high ndvi values.
        (void) set_bit(22, pxout.qabits);
        ndvi = (m02 - m01) / (m02 + m01);
        if(ndvi <= swc_ndvi[0] || ndvi >= swc_ndvi[1]) {
          conf = 1.0;
          (void) set_bit(22, pxout.testbits);
        }

      }
    }

/**********************************************************************/

//  printf("CC NDVI: %d %f %f %f %f %f %f \n", irclr, m01, m02, ndvi,
//  		swc_ndvi[0], swc_ndvi[1], conf);

    return conf;
    
/**********************************************************************/    
    
}
