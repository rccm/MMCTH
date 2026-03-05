/*
!C ********************************************************************
!Description:

  Integer function polarnight_snow.c
  Computes initial clear sky confidence for polar nighttime snow/ice 
    surfaces.
  Called from polar_night.c.

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
  Converted to C            R. Frey           11/2007


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


int polarnight_snow()


/*********************************************************************/

{

/*  Declarations */

    extern float cithr(int, float, float);
    extern float conf_test(float, float, float, float, float);
    extern void set_bit(int, unsigned char[]);
    extern void get_pn_thresholds(float, float*, float*, float*, float*,
                        float*, float*, float*, float*, float*, float*);
    
    int groups;
    int k;
    int ngtests[num_test_groups];
    int num_tests;
    int return_code;

    float btd_thr;
    float cmin1, cmin2, cmin5;
    float conf;
    float cosvza;
    float fac;
    float m22, m27, m28, m31, m32; 
    float m31_m32, m31_m22, m22_m32, m28_m31;
    float midpt, locut, hicut;
    float power;
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
    cmin5 = 1.0;

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
    m28 = pxin.rad[27];
    m31 = pxin.rad[30];
    m32 = pxin.rad[31];

/*********************************************************************/

/*  Perform cloud tests. */
//  printf("\nPerforming cloud tests in polarnight_snow.c\n");

/*********************************************************************/

/*  Group 1 tests. */

/*********************************************************************/

/*  Water vapor high cloud test. */

    if( !(pxin.Antarctic && pxin.land) ) {
      if(rintf(m27) != rintf(bad_data)) {
    	  
        num_tests++;
        (void) set_bit(15, pxout.qabits);
        if(m27 >= pnsh20[1]) (void) set_bit(15, pxout.testbits);
        conf = conf_test(m27, pnsh20[0], pnsh20[2], pnsh20[3],
                                  pnsh20[1]);
        cmin1 = min(cmin1, conf);
        ngtests[0]++;
//      printf("h20 %f %f %f %d %d\n", m27,conf,cmin1,ngtests[0],num_tests);
        
      }
    }

/*********************************************************************/

/*  Group 2 tests. */

/*********************************************************************/

/*  11-12 micron BTD test for transmissive cirrus clouds. */

    if( !(pxin.Antarctic && pxin.land) ) {
      if(rintf(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data)
                 && pxin.vza > 0.0) {

        m31_m32 = m31 - m32;

        num_tests++;
        (void) set_bit(18, pxout.qabits);

/*      Get secant of viewing zenith angle. */
        cosvza = cosf(pxin.vza * dtr);
        if(fabsf(cosvza) > 0.0) {
          schi = 1.0 / cosvza; }
        else {
          schi = 99.0;
        }

/*      Interpolate look-up table values of 11-12 micron BTD thresholds
        (functions of viewing zenith angle and 11 micron BT).
*/
        thr = cithr(1, schi, m31);
        if(thr < 0.1 || fabsf(schi-99.0) < 0.0001) {
          btd_thr = pns11_12hi[0]; }
        else {
/*        Add adjustment for snow cover. */
          btd_thr = thr + pn11_12adj[0];
        }
      
        locut = btd_thr;
        midpt = btd_thr - (0.3 * btd_thr);
        if(m31 < pnbt1[0]) {
          hicut = midpt - (0.2 * btd_thr); }
        else {
          hicut = midpt - 1.25;
        }
        
        if (m31_m32 <= midpt) (void) set_bit(18, pxout.testbits);
        
        conf = conf_test(m31_m32,locut,hicut,1.0,midpt);
        cmin2 = min(cmin2, conf);
        ngtests[1]++;

//      printf("m31_32: %f %f %f %f %f %f %f %f %d %d\n", pnbt1[0], m31_m32,
//          btd_thr, locut, midpt, hicut, conf, cmin2, num_tests, ngtests[1]);
        
      }
    }
    
/**********************************************************************/

/*  11 minus 4 micron BTDIF fog and low cloud test. */

    if (rintf(m31) != rintf(bad_data) && rintf(m22) != rintf(bad_data)) {

      if(m31 > pnbt3[0]) {

        num_tests++;
        (void) set_bit(19, pxout.qabits);
        m31_m22 = m31 - m22;

//      Get thresholds for polar snow/ice surface.
        (void) get_pn_thresholds(m31,bt_11_bounds,pn_11_4l,pn_11_4m1,
                                 pn_11_4m2,pn_11_4m3,pn_11_4h,&locut,&hicut,
                                 &midpt,&power);     
      
        if (m31_m22 <= midpt) (void) set_bit(19, pxout.testbits);
        conf = conf_test(m31_m22, locut, hicut, 1.0, midpt);
        cmin2 = min(cmin2, conf);
        ngtests[1]++;

//      printf("11-4: %f %f %f %f %f %f %d %d\n", m31_m22, locut, midpt,
//            hicut, conf, cmin2, num_tests, ngtests[1]);

      }
    }

/**********************************************************************/

/*  7.3 minus 11 micron cloud test. */
    if (rintf(m31) != rintf(bad_data) && rintf(m28) != rintf(bad_data)) {
      if(m31 <= pnbt2[0]) {

        num_tests++;
        (void) set_bit(23, pxout.qabits);
        m28_m31 = m28 - m31;

        if ( pxin.land ) {
/*        Get thresholds for land (snow) surface. */
          (void) get_pn_thresholds(m31,bt_11_bnds2,pn_7_11l,pn_7_11m1,
                                  pn_7_11m2,pn_7_11m3,pn_7_11h,&locut,&hicut,
                                  &midpt,&power);
        }
        else {
/*        Get thresholds for water (ice) surface. */
          (void) get_pn_thresholds(m31,bt_11_bnds2,pn_7_11lw,pn_7_11m1w,
                                  pn_7_11m2w,pn_7_11m3w,pn_7_11hw,&locut,&hicut,
                                  &midpt,&power);
        }

        if (m28_m31 > midpt) (void) set_bit(23, pxout.testbits);
        conf = conf_test(m28_m31, locut, hicut, 1.0, midpt);
        cmin2 = min(cmin2, conf);
        ngtests[1]++;

//      printf("7-11: %f %f %f %f %f %f %f %d %d\n", pnbt2[0], m28_m31,
//         locut, midpt, hicut, conf, cmin2, num_tests, ngtests[1]);

      }
    }

/**********************************************************************/

/*  Group 5 tests. */

/*********************************************************************/

/*  3.7-12um brightness temperature difference test for thin cirrus. */

    if (rintf(m32) != rintf(bad_data) && rintf(m22) != rintf(bad_data)
    		       && rintf(m31) != rintf(bad_data)) {
      if( !pxin.hi_elev) {
    	  
        num_tests++;
        (void) set_bit(17, pxout.qabits);
        m22_m32 = m22 - m32;

//      Get thresholds for polar snow/ice surface.
        (void) get_pn_thresholds(m31,bt_11_bounds,pn_4_12l,pn_4_12m1,
                                      pn_4_12m2,pn_4_12m3,pn_4_12h,&locut,
                                      &hicut,&midpt,&power);
      
        if (m22_m32 <= midpt) (void) set_bit(17, pxout.testbits);
        conf = conf_test(m22_m32, locut, hicut, 1.0, midpt);
        cmin5 = min(cmin5, conf);
        ngtests[4]++;

//      printf("3.7-12: %f %f %f %f %f %f %d %d\n", m22_m32, locut,
//         midpt, hicut, conf, cmin5, num_tests, ngtests[4]);
        
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

    pre_confidence = cmin1 * cmin2 * cmin5;

/*  Find the number of test groups used. */
    k = 0;
    while(k < num_test_groups) {
      if(ngtests[k] > 0) groups++;
      k++;
    }

/*  Get power for confidence calculation. */
    if(groups > 0) fac = 1.0 / groups;

//  printf("\npre_conf, groups, fac: %f %d %f\n", pre_confidence, groups,
//            fac);

/*  Get initial confidence of clear sky. */
    pxout.init_conf = powf(pre_confidence, fac);

    return_code = num_tests;
    return return_code;

}
