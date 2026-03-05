/*
************************************************************************
!F77 

 !Description:

   Performs clear sky restoral tests in sun-glint conditions.

!Input parameters:
   none

   Inputs stored in structure variable 'pxin' of type "pixel_in" defined
   in pixel.h. 

!Output Parameters:
   conf        modified clear sky confidence for pixel

   Returned through function call

!Revision History:

   06/04 Collection 5  R. Frey:
   11/07 Converted to C  R. Frey

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


float chk_sunglint()


{

/***********************************************************************/
	
/*  Declarations */
	
    extern void get_regdif(int);
    extern float get_regstd(int);
    extern int spatial_var();
    extern int check_bits(int, unsigned char[]);
    extern void set_bit(int, unsigned char[]);

    int band;
    int var;
    int irclr;
    int bit_test;
    int qa_test;

    float m20, m31, m17, m18, m09;
    float d37_11;
    float rat;
    float stdev;
    float conf;

/**********************************************************************/
    
/*  Initializations */

//  printf("Executing chk_sunglint \n");
    conf = pxout.intermediate_conf;
    d37_11 = -99.9;
    rat = -99.9;
    stdev = 0.0;
    rg_var.mean = 0.0;
    
/**********************************************************************/    

//  Get brightness temperature differences between pixel of interest
//  and the 8 surrounding it.

    band = 31;
    (void) get_regdif(band);

//  Check variation in the region.
    var = spatial_var();
        
    if(var == 1) {

//    Region is uniform in the 11 micron IR window.
//    Determine the logical flag 'irclr' - true if ir cloud tests below
//    have all been passed. 11-12 thin cirrus test makes final decision
//    for bit 18.

      irclr = 1;
      bit_test = check_bits(13, pxout.testbits);
      qa_test = check_bits(13, pxout.qabits);
      if( qa_test && !bit_test ) irclr = 0;
      bit_test = check_bits(14, pxout.testbits);
      qa_test = check_bits(14, pxout.qabits);
      if( qa_test && !bit_test ) irclr = 0;
      bit_test = check_bits(15, pxout.testbits);
      qa_test = check_bits(15, pxout.qabits);
      if( qa_test && !bit_test ) irclr = 0;
      bit_test = check_bits(18, pxout.testbits);
      qa_test = check_bits(18, pxout.qabits);
      if( qa_test && !bit_test ) irclr = 0;
      bit_test = check_bits(27, pxout.testbits);
      qa_test = check_bits(27, pxout.qabits);
      if( qa_test && !bit_test ) irclr = 0;

      if( irclr ) {

        m20 = pxin.rad[19];
        m31 = pxin.rad[30];
        m17 = pxin.rad[16];
        m18 = pxin.rad[17];
        m09 = pxin.rad[8];

        if(rintf(m20) != rintf(bad_data) && rintf(m31) != rintf(bad_data) &&
           rintf(m17) != rintf(bad_data) && rintf(m18) != rintf(bad_data) ) {

//        Set bit to indicate test was performed.
          (void) set_bit(26, pxout.qabits);

          d37_11 = m20 - m31;

//        Set bit if tests passed.
          if(d37_11 >= sg_tbdfl[0]) {

            conf = 0.67;
            rat = m17 / m18;

//          if(rat > snglrat[0] && (rintf(m09) != rintf(bad_data)) ) {
//          Reflectance value of 0.568 simulates saturation (bad data).            
            if(rat > snglrat[0] && m09 < 0.568) {
              
              (void) set_bit(26, pxout.testbits);
              conf = 0.96;
              
            }
            else {
            	
              band = 2;
              stdev = get_regstd(band);
            
              if( rintf(rg_var.mean) != rintf(bad_data) ) {
            	if( (stdev * rg_var.mean) < 0.001) {  
                	  
                  (void) set_bit(26, pxout.testbits);
                  conf = 0.96;
            	
            	}
              }
              
            }
            
          }

        }

      }
      
//    printf("Sun-glint CSR: %d %d %f %f %f %f %f %f %f %f %f \n", var,irclr,d37_11,sg_tbdfl[0],rat,snglrat[0],
//      		              m09,stdev,rg_var.mean,stdev*rg_var.mean,conf);  

    }

    
    
    return conf;

}
