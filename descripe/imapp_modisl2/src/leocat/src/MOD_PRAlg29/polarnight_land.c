/*
!C ********************************************************************
!Description:

  Integer function polarnight_land.c
  Computes initial clear sky confidence for polar nighttime land surfaces.
  Called from polar_night.c

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
  Converted to C             R. Frey           11/2007


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


int polarnight_land()


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

    float a;
    float btd_thr;
    float cmin1, cmin2, cmin5;
    float conf;
    float corr;
    float cosvza;
    float df1;
    float df2;
    float fac;
    float hicut;
    float locut;
    float lst_thresh;
    float m22, m27, m28, m31, m32; 
    float m31_m32, m31_m22, m28_m31, m22_m32;
    float midpt;
    float power;
    float pre_confidence;
    float schi;
    float sfcdif;
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
//  printf("\nPerforming cloud tests in landnight.c\n");

/*********************************************************************/    

//  Group 1 Tests
    
/*********************************************************************/     

/*  Water vapor high cloud test. */

    if(rintf(m27) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(15, pxout.qabits);
      if(m27 >= pnlh20[1]) (void) set_bit(15, pxout.testbits);
      conf = conf_test(m27, pnlh20[0], pnlh20[2], pnlh20[3],
                       pnlh20[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("h20 %f %f %f %d %d\n", m27,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

//  Surface Temperature Test.

    if ( rintf(m31) != rintf(bad_data) && !pxin.hi_elev &&
         rintf(m32) != rintf(bad_data) && pxin.eco_type != 8 ) {

      if (pxin.sfctmp > 0.0 && pxin.sfctmp < 350.0) {

        df1 = m31 - m32;
        df2 = m31 - m22;

        num_tests++;
        (void) set_bit(27, pxout.qabits);

        if(pxin.desert) {
          lst_thresh = 20.0;
        }
        else if(df1 >= 0.0 || (df1 < 0.0 && (df2 <= -0.5 || df2 >= 1.0)) ) {
          lst_thresh = 12.0;
        }
        else {
          lst_thresh = 20.0;
        }

        if(df1 >= 1.0) {
          midpt = lst_thresh + (2.0 * rintf(df1));
        }
        else {
          midpt = lst_thresh;
        }
        
        a = pxin.vza / max_vza;
        corr = powf(a, 4) * 3.0;
        midpt = midpt + corr;
        locut = midpt + 2.0;
        hicut = midpt - 2.0;

        sfcdif = pxin.sfctmp - m31;

        if( sfcdif < midpt ) (void) set_bit(27, pxout.testbits);

        conf = conf_test(sfcdif, locut, hicut, 1.0, midpt);
        cmin1 = min(cmin1, conf);
        ngtests[0]++;
//      printf("SFCT: %f %f %f %f %f %f %f %d %d \n", m31, pxin.sfctmp, sfcdif,
//        		locut, midpt, hicut, conf, ngtests[0], num_tests);

      }

    }

/*********************************************************************/    

//  Group 2 Tests
        
/*********************************************************************/   

//  11-12um BTD test for thin cirrus.
    
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
        btd_thr = pnl11_12hi[0]; }
      else {
//      Add 0.2 for likely snow cover.    	  
        btd_thr = thr + 0.2;
      }
      
      locut = btd_thr;
      midpt = btd_thr - (0.3 * btd_thr);
      if(m31 < pnlbt1[0]) {
        hicut = midpt - (0.2 * btd_thr);
      }
      else {
    	hicut = midpt - 1.25;
      }
      
      if (m31_m32 <= midpt) (void) set_bit(18, pxout.testbits);
      
      conf = conf_test(m31_m32, locut, hicut, 1.0, midpt);
      cmin2 = min(cmin2, conf);
      ngtests[1]++;

//    printf("m31_32: %f %f %f %f %f %f %f %f %d %d\n", pnlbt1[0], pnl11_12hi[0], m31_m32,
//              locut, midpt, hicut, conf, cmin2, num_tests, ngtests[1]);

    }
    
/**********************************************************************/

/*  11 minus 4 micron BTDIF fog and low cloud test. */
      
/**********************************************************************/

    if (rintf(m31) != rintf(bad_data) && rintf(m22) != rintf(bad_data)) {

      num_tests++;
      (void) set_bit(19, pxout.qabits);
      
      m31_m22 = m31 - m22;
      df1 = m31 - m32;
      (void) get_pn_thresholds(m31,bt_11_bounds,pn_11_4l,pn_11_4m1,
                               pn_11_4m2,pn_11_4m3,pn_11_4h,&locut,
                               &hicut,&midpt,&power);     
      
      if (m31_m22 <= midpt) (void) set_bit(19, pxout.testbits);
      conf = conf_test(m31_m22, locut, hicut, 1.0, midpt);
      cmin2 = min(cmin2, conf);
      ngtests[1]++;
//    printf("11-4: %f %f %f %f %f %f %d %d\n", m31_m22, locut, midpt,
//           hicut, conf, cmin2, num_tests, ngtests[1]);

    }

/**********************************************************************/
    
//  7.3-11um BTD test for thick, mid-level clouds.
    
/**********************************************************************/

    if (rintf(m31) != rintf(bad_data) && rintf(m28) != rintf(bad_data)) {

      if(m31 < pnlbt2[0]) {

        num_tests++;
        (void) set_bit(23, pxout.qabits);
        m28_m31 = m28 - m31;

//      Use polar night snow thresholds here. Logic is that land surfaces
//      poleward of 60 latitude are snow covered most of the year and
//      many times ancillary snow map is incomplete.        
        (void) get_pn_thresholds(m31,bt_11_bnds2,pn_7_11l,pn_7_11m1,
                                 pn_7_11m2,pn_7_11m3,pn_7_11h,&locut,
                                 &hicut,&midpt,&power);  
        
        if ( m28_m31 > midpt ) (void) set_bit(23, pxout.testbits);
        conf = conf_test(m28_m31, locut, hicut, 1.0, midpt);
        cmin2 = min(cmin2, conf);
        ngtests[1]++;
//      printf("7.3-11: %f %f %f %f %f %f %f %f %f %d %d \n", pnlbt2[0], m28, m31, m28_m31,
//        		nl7_11s[0], nl7_11s[1], nl7_11s[2], conf, cmin2, num_tests,
//        		ngtests[1]);

      }

    }

/*********************************************************************/    

//  Group 5 Tests
            
/*********************************************************************/   

//  4-12 um BTD thin cirrus test.
    
/*********************************************************************/

    if (rintf(m32) != rintf(bad_data) && rintf(m22) != rintf(bad_data)) {

      m22_m32 = m22 - m32;
      num_tests++;
      (void) set_bit(17, pxout.qabits);
      
      (void) get_pn_thresholds(m31,bt_11_bounds,pn_4_12l,pn_4_12m1,
                               pn_4_12m2,pn_4_12m3,pn_4_12h,&locut,
                               &hicut,&midpt,&power);  
      
      if (m22_m32 <= midpt) (void) set_bit(17, pxout.testbits);
      conf = conf_test(m22_m32, locut, hicut, 1.0, midpt);
      cmin5 = min(cmin5, conf);
      ngtests[2]++;
//    printf("4-12: %f %f %f %f %f %f %d %d \n", m22_m32, nl4_12hi[0],
//   		  nl4_12hi[1], nl4_12hi[2], cmin5, conf, num_tests, ngtests[2]);
      
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

//  printf("\npre_conf, groups, fac: %f %d %f\n", pre_confidence, groups, fac);

/*  Get initial confidence of clear sky. */
    pxout.init_conf = powf(pre_confidence, fac);

    return_code = num_tests;
    return return_code;

}
