/*
!C ********************************************************************
!Description:

  Integer function polarday_snow.c
  Computes initial clear sky confidence for polar daytime snow/ice surfaces.
  Called from polar_day.c.

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


int polarday_snow()


/*********************************************************************/

{

/*  Declarations */

	extern void get_pn_thresholds(float, float*, float*, float*, float*,
	            float*, float*, float*, float*, float*, float*);
    extern float cithr(int, float, float);
    extern float conf_test(float, float, float, float, float);
    extern void set_bit(int, unsigned char[]);

    int groups;
    int k;
    int ngtests[num_test_groups];
    int num_tests;
    int return_code;

    float btd_thr;
    float cmin1, cmin2, cmin4;
    float conf;
    float cosvza;
    float fac;
    float m22, m26, m27, m31, m32;
    float m31_m32, m22_m31;
    float locut, hicut, midpt, power;
    float pre_confidence;
    float schi;
    float thr;

/*********************************************************************/

/*  Initializations */
    return_code = 0;

/*  Counter for total number tests performed. */
    num_tests = 0;

/*  Test group confidences. */
    cmin1 = 1.0;
    cmin2 = 1.0;
    cmin4 = 1.0;

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
    m26 = pxin.rad[25];
    m27 = pxin.rad[26];
    m31 = pxin.rad[30];
    m32 = pxin.rad[31];

/*********************************************************************/

/*  Perform cloud tests. */
//  printf("\nPerforming cloud tests in polarday_snow.c\n");

/*********************************************************************/

/*  Group 1 tests. */

/*********************************************************************/

/*  Water vapor high cloud test. */

    if(rintf(m27) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(15, pxout.qabits);
      if(m27 >= dpsh20[1]) (void) set_bit(15, pxout.testbits);
      conf = conf_test(m27, dpsh20[0], dpsh20[2], dpsh20[3],
                                  dpsh20[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("h20 %f %f %f %d %d\n", m27,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

/*  Group 2 tests. */

/*********************************************************************/

/*  11-12 micron BTD test for transmissive cirrus clouds. */

    if(rint(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data) &&
                    pxin.vza > 0.0) {

      m31_m32 = m31 - m32;

      num_tests++;
      (void) set_bit(18, pxout.qabits);

/*    Get secant of viewing zenith angle. */
      cosvza = cosf(pxin.vza * dtr);
      if(fabsf(cosvza) > 0.0) {
        schi = 1.0 / cosvza; }
      else {
        schi = 99.0;
      }

/*    Interpolate look-up table values of 11-12 micron BTD thresholds
      (functions of viewing zenith angle and 11 micron BT).
*/
      thr = cithr(1, schi, m31);
      if(thr < 0.1 || fabsf(schi-99.0) < 0.0001) {
        btd_thr = ds11_12hi[0]; }
      else {
/*      Add adjustment for snow cover. */
        btd_thr = thr + ds11_12adj[0];
      }

      if (m31_m32 <= btd_thr) (void) set_bit(18, pxout.testbits);

      locut = btd_thr + (0.3 * btd_thr);
      hicut = btd_thr - (0.3 * btd_thr);
      conf = conf_test(m31_m32,locut,hicut,1.0,btd_thr);
      cmin2 = min(cmin2, conf);
      ngtests[1]++;

//    printf("m31_32: %f %f %f %f %f %f %d %d\n", m31_m32, locut,
//              btd_thr, hicut, conf, cmin2, num_tests, ngtests[1]);

    }

/**********************************************************************/

/*  11 minus 4 micron BTDIF fog and low cloud test. */

    if (pxin.visusd) {
      if (rintf(m31) != rintf(bad_data) && rintf(m22) != rintf(bad_data)) {

    	if(m31 > dpsbt1[0]) {

          num_tests++;
          (void) set_bit(19, pxout.qabits);
          m22_m31 = m22 - m31;

//        Get thresholds for land (snow) surface.
          (void) get_pn_thresholds(m31, bt_11_bnds3, dps4_11l, dps4_11m1,
                                   dps4_11m2, dps4_11m3, dps4_11h, &locut,
                                   &hicut, &midpt, &power);

          if (m22_m31 <= midpt) (void) set_bit(19, pxout.testbits);
          conf = conf_test(m22_m31, locut, hicut, 1.0, midpt);
          cmin2 = min(cmin2, conf);
          ngtests[1]++;
//        printf("snow 4-11: %f %f %f %f %f %f %f %d %d\n", dpsbt1[0], m22_m31, locut,
//                   midpt, hicut, conf, cmin2, num_tests, ngtests[1]);

    	}

      }
    }

/**********************************************************************/

/*  Group 4 tests. */

/*********************************************************************/

/*  Near-infrared (1.38 micron) high cloud test. */

    if (pxin.visusd && pxin.hi_elev == 0 && pxin.tpw > dps_ref3_tpw[0]) {
      if (rintf(m26) != rintf(bad_data)) {

        num_tests++;
        (void) set_bit(16, pxout.qabits);

        if (m26 <= dpsref3[1]) (void) set_bit(16, pxout.testbits);
        conf = conf_test(m26, dpsref3[0], dpsref3[2], 1.0, dpsref3[1]);
        cmin4 = min(cmin4, conf);
        ngtests[3]++;
//      printf("1.38: %f %f %f %f %f %f %d %d\n", m26,dpsref3[0],
//           dpsref3[1],dpsref3[2],conf,cmin4,num_tests,ngtests[3]);

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

    pre_confidence = cmin1 * cmin2 * cmin4;

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
