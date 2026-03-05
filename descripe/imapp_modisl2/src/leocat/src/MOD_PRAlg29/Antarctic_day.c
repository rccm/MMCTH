/*
!C ********************************************************************
!Description:

  Integer function Antarctic_day.c
  Computes initial clear sky confidence for daytime Antarctic surfaces.
  Called from polar_day.c

!Input arguments:
  none

  Inputs stored in structure variable 'pxin' of type "pixel_in" defined
  in pixel.h.
  Test thresholds accessed through thresholds.h.

!Output arguments:
  none

  int return_code   returns number of cloud tests performed
                    (returned through function call)

  Output confidence and associated products stored in structure
  variable 'pxout' of type "pixel_out" defined in pixel.h

!Revision History:
  Original version taken from MODIS Cloud Mask FORTRAN code developed
  by CIMSS, UW-Madison.
  Converted to C            R. Frey           10/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>
#include "pixel.h"
#include "thresholds.h"
#include "mask_processing_constants.h"

/*********************************************************************/


int Antarctic_day()


/*********************************************************************/

{

/*  Declarations */

    extern void get_pn_thresholds(float, float*, float*, float*, float*,
                      float*, float*, float*, float*, float*, float*);
    extern float conf_test(float, float, float, float, float);
    extern void set_bit(int, unsigned char[]);
    
    int groups;
    int k;
    int ngtests[num_test_groups];
    int num_tests;
    int return_code;

    float cmin1, cmin2;
    float conf;
    float fac;
    float m22, m27, m31; 
    float m22_m31;
    float pre_confidence;
    float midpt;
    float locut;
    float hicut;
    float power;
    
/*********************************************************************/

/*  Initializations */   
    return_code = 0;

/*  Counter for total number tests performed. */
    num_tests = 0;

/*  Test group confidences. */
    cmin1 = 1.0;
    cmin2 = 1.0;

/*  Number of test groups used. */
    groups = 0;

/*  Inverse of 'groups' used as power in calculating confidence. */
    fac = 1.0;

/*  Counter for number tests performed in each test group. */
    k = 0;
    while(k < num_test_groups) {
      ngtests[k] = 0;
      k++;
    }

/*  Input radiance data. */
    m22 = pxin.rad[21];
    m27 = pxin.rad[26];
    m31 = pxin.rad[30];

/*********************************************************************/

/*  Perform cloud tests. */
//  printf("\nPerforming cloud tests in Antarctic_day.c\n");

/*********************************************************************/

/*  Group 1 tests. */

/*********************************************************************/

/*  Water vapor high cloud test. */

    if(rintf(m27) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(15, pxout.qabits);
      if(m27 > anth20[1]) (void) set_bit(15, pxout.testbits);
      conf = conf_test(m27, anth20[0], anth20[2], anth20[3],
                                  anth20[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("h20 %f %f %f %d %d\n", m27,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

/*  Group 2 tests. */

/*********************************************************************/

/*  11 minus 4 micron BTDIF fog and low cloud test. */

    if (pxin.visusd) {
      if (rintf(m31) != rintf(bad_data) && rintf(m22) != rintf(bad_data)) {

    	if(m31 > antbt1[0]) {
    		
          num_tests++;
          (void) set_bit(19, pxout.qabits);
          m22_m31 = m22 - m31;

//        Get thresholds for land (snow) surface.
          (void) get_pn_thresholds(m31,bt_11_bnds4, ant4_11l, ant4_11m1,
                                   ant4_11m2, ant4_11m3, ant4_11h, &locut,
                                   &hicut, &midpt, &power);

          if (m22_m31 <= midpt) (void) set_bit(19, pxout.testbits);
          conf = conf_test(m22_m31, locut, hicut, 1.0, midpt);
          cmin2 = min(cmin2, conf);
          ngtests[1]++;
//        printf("4-11: %f %f %f %f %f %f %f %d %d\n", antbt1[0], m22_m31, locut,
//                   midpt, hicut, conf, cmin2, num_tests, ngtests[1]);
          
    	}

      }
    }

/*********************************************************************/
    
/*
    printf("\nvalue of qa byte 0 %d\n", pxout.qabits[0]);
    printf("value of qa byte 1 %d\n", pxout.qabits[1]);
    printf("value of qa byte 2 %d\n", pxout.qabits[2]);
    printf("value of qa byte 3 %d\n", pxout.qabits[3]);
    printf("value of qa byte 4 %d\n", pxout.qabits[4]);
    printf("value of qa byte 5 %d\n", pxout.qabits[5]);

    printf("\nvalue of cm byte 0 %d\n", pxout.testbits[0]);
    printf("value of cm byte 1 %d\n", pxout.testbits[1]);
    printf("value of cm byte 2 %d\n", pxout.testbits[2]);
    printf("value of cm byte 3 %d\n", pxout.testbits[3]);
    printf("value of cm byte 4 %d\n", pxout.testbits[4]);
    printf("value of cm byte 5 %d\n", pxout.testbits[5]);
*/
/*********************************************************************/

/*  Determine initial confidence based on group values. */

    pre_confidence = cmin1 * cmin2;

/*  Find the number of test groups used. */
    k = 0;
    while(k < num_test_groups) {
      if(ngtests[k] > 0) groups++;
      k++;
    }

/*  Get power for confidence calculation. */
    if(groups > 0) fac = 1.0 / groups;

//  printf("\npre_conf, groups, fac: %f %d %f\n", pre_confidence, groups, fac);

/*  Get initial confidence of clear sky. */
    pxout.init_conf = powf(pre_confidence, fac);

    return_code = num_tests;
    return return_code;

}
