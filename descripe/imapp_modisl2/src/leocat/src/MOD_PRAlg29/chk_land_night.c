/*
************************************************************************
!F77

!Description:
   Perform final clear-sky confidence check on nighttime land pixels.
   If confidence of clear sky is low but temp is warm and if the
   IR clear sky tests all passed, then use the 11 um brightness
   temperature to assign a final confidence.

!Input Parameters:
   none

!Output Parameters:
   conf

!Revision History:
   Converted to C      11/07        R. Frey

!Team-Unique Header:

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


float chk_land_night()


{
	

/**********************************************************************/
	
//  Declarations
	
    extern int check_bits(int, unsigned char[]);
    extern void set_bit(int, unsigned char[]);

    int bit_test;
    int qa_test;
    int irclr;
    
    float m31;
    float hds11[3];
    float conf;
	
/**********************************************************************/	

//	Initializations
    
//  printf("Executing chk_land_night \n");
    conf = pxout.intermediate_conf;
	
    m31 = pxin.rad[30];
    
    hds11[0] = 0.0;
    hds11[1] = 0.0;
    hds11[2] = 0.0;
    
/**********************************************************************/	    

//  Check IR clear sky tests.
    irclr = 1;
    bit_test = check_bits(14, pxout.testbits);
    qa_test = check_bits(14, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    bit_test = check_bits(15, pxout.testbits);
    qa_test = check_bits(15, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    bit_test = check_bits(17, pxout.testbits);
    qa_test = check_bits(17, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    bit_test = check_bits(23, pxout.testbits);
    qa_test = check_bits(23, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;

/**********************************************************************/    
    
    if (irclr) {

      if (rintf(m31) != rintf(bad_data)) {

//      Get elevation-adjusted 11 micron brightness temperature threshold.

        (void) set_bit(26, pxout.qabits);

        hds11[0] = lnbt11[0] - pxin.tbadj;
        hds11[1] = lnbt11[1] - pxin.tbadj;
        hds11[2] = lnbt11[2] - pxin.tbadj;

//      Check for hot scene.
        if(m31 > hds11[0]) {

//        Assign confidence level based on 11 micron Tbb.
          if(m31 > hds11[2]) {
//          Assign pixel to confident clear, set bit #26.
            conf = 1.0;
            (void) set_bit(26, pxout.testbits);
          }
          else if(m31 > hds11[1]) {
//          Assign pixel to probably clear, set bit #26.
            conf = 0.96;
            (void) set_bit(26, pxout.testbits);
          }
          else {
//          Assign pixel to uncertain, do not set bit #26.
            conf = 0.95;
          }
          
        }

      }
      
//    printf("chk_land_night m31: %f %f %f %f %f %f \n", pxin.tbadj,
//               hds11[0], hds11[1], hds11[2], m31, conf);

    }  

    return conf;

}
