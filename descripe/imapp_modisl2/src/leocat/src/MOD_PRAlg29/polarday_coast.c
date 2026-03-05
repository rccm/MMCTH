/*
!C ********************************************************************
!Description:

  Integer function polarday_coast.c
  Computes initial clear sky confidence for polar coastal, non-arid
    daytime surfaces.
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
  Original version taken from MODIS cloud mask FORTRAN code developed
  at CIMSS, UW-Madison.
  Convertedd to C             R. Frey           11/2007


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


int polarday_coast()


/*********************************************************************/

{

/*  Declarations */

    extern int get_b1_thresholds(float*, float*, float*, float*);
    extern float cithr(int, float, float);
    extern float conf_test(float, float, float, float, float);
    extern void set_bit(int, unsigned char[]);

    int irtn;
    int groups;
    int k;
    int ngtests[num_test_groups];
    int num_tests;
    int return_code;

    float b1_locut, b1_midpt, b1_hicut, b1_power;
    float btd_thr;
    float cmin1, cmin2, cmin3, cmin4;
    float conf;
    float cosvza;
    float fac;
    float m01, m02, m22, m26, m27, m31, m32;
    float m31_m32, m31_m22;
    float locut, hicut;
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
    cmin3 = 1.0;
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
    m01 = pxin.rad[0];
    m02 = pxin.rad[1];
    m22 = pxin.rad[21];
    m26 = pxin.rad[25];
    m27 = pxin.rad[26];
    m31 = pxin.rad[30];
    m32 = pxin.rad[31];

/*********************************************************************/

/*  Perform cloud tests. */
//  printf("\nPerforming cloud tests in polarday_coast.c\n");

/*********************************************************************/

/*  Group 1 tests. */

/*********************************************************************/

/*  Water vapor high cloud test. */

    if(rintf(m27) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(15, pxout.qabits);
      if(m27 >= pdlh20_t2[1]) (void) set_bit(15, pxout.testbits);
      conf = conf_test(m27, pdlh20_t2[0], pdlh20_t2[2], pdlh20_t2[3],
                                  pdlh20_t2[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("h20 %f %f %f %d %d\n", m27,conf,cmin1,num_tests,
//            ngtests[0]);
    }

/*********************************************************************/

/*  Group 2 tests. */

/*********************************************************************/

/*  11-12 micron BTD test for transmissive cirrus clouds. */

    if(rintf(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data) &&
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
        btd_thr = dl11_12hi_t2[0]; }
      else {
        btd_thr = thr;
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

/*    11 minus 4 micron BTDIF fog and low cloud test. */

    if (pxin.visusd) {
      if (rintf(m31) != rintf(bad_data) && rintf(m22) != rintf(bad_data)) {

        num_tests++;
        (void) set_bit(19, pxout.qabits);
        m31_m22 = m31 - m22;

        if(m31_m22 >= pdl11_4lo_t2[1]) (void) set_bit(19, pxout.testbits);
        conf = conf_test(m31_m22, pdl11_4lo_t2[0], pdl11_4lo_t2[2], 1.0,
                         pdl11_4lo_t2[1]);
        cmin2 = min(cmin2, conf);
        ngtests[1]++;
//      printf("11-4: %f %f %f %f %f %f %d %d\n", m31_m22, pdl11_4lo_t2[0],
//       pdl11_4lo_t2[1], pdl11_4lo_t2[2], conf, cmin2, num_tests,
//       ngtests[1]);

      }
    }

/**********************************************************************/

/*  Group 3 tests. */

/**********************************************************************/

/*  Visible reflectance threshold test. */

    if (pxin.visusd) {
      if (rintf(m01) != rintf(bad_data)) {

        num_tests++;
        (void) set_bit(20, pxout.qabits);

        irtn = get_b1_thresholds(&b1_locut, &b1_hicut, &b1_midpt, &b1_power);

        if(irtn != 0) {
          b1_locut = pdlref1[0];
          b1_midpt = pdlref1[1];
          b1_hicut = pdlref1[2];
          b1_power = pdlref1[3];
        }

        if (m01 <= b1_midpt) (void) set_bit(20, pxout.testbits);
        conf = conf_test(m01, b1_locut, b1_hicut, b1_power, b1_midpt);
        cmin3 = min(cmin3, conf);
        ngtests[2]++;

//      printf("B1: %d %f %f %f %f %f %f %f\n", irtn, m01, b1_locut,
//        		    b1_midpt, b1_hicut, b1_power, conf, cmin3);
//      printf("Vis: %f %f %f %f %f %f %d %d\n", m01, pdlref1_t2[0],
//        pdlref1_t2[1], pdlref1_t2[2], conf, cmin3, num_tests,
//        ngtests[2]);

      }
    }

/*********************************************************************/

/*  Group 4 tests. */

/*********************************************************************/

/*  Near-infrared (1.38 micron) high cloud test. */

    if (pxin.visusd && pxin.hi_elev == 0 && pxin.tpw > pdl_ref3_tpw_t2[0]) {
      if (rintf(m26) != rintf(bad_data)) {

        num_tests++;
        (void) set_bit(16, pxout.qabits);

        if (m26 <= pdlref3_t2[1]) (void) set_bit(16, pxout.testbits);
        conf = conf_test(m26, pdlref3_t2[0], pdlref3_t2[2], 1.0,
                    pdlref3_t2[1]);
        cmin4 = min(cmin4, conf);
        ngtests[3]++;
//      printf("1.38: %f %f %f %f %f %f %d %d\n", m26,pdlref3_t2[0],
//        pdlref3_t2[1],pdlref3_t2[2],conf,cmin4,num_tests,ngtests[3]);

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

    pre_confidence = cmin1 * cmin2 * cmin3 * cmin4;

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
