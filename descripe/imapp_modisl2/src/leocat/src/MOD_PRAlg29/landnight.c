/*
!C ********************************************************************
!Description:

  Integer function landnight.c
  Computes initial clear sky confidence for nighttime land surfaces.
  Called from land_night.c

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
  Converted to C             R. Frey           10/2007


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


int landnight()


/*********************************************************************/

{

/*  Declarations */

    extern float cithr(int, float, float);
    extern float conf_test(float, float, float, float, float);
    extern float conf_test_dble(float, float[], float[], float, float[]);
    extern void set_bit(int, unsigned char[]);
    extern void get_nl_thresholds(float, float *, float *, float *, float *);
    
    int Australia;
    int groups;
    int k;
    int ngtests[num_test_groups];
    int num_tests;
    int return_code;

    float locutr[2], midptr[2], hicutr[2];

    float a;
    float b1=-0.0077, b2=1.1234, b3=-0.3403;
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
    float m22, m27, m28, m31, m32, m35; 
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
    m35 = pxin.rad[34];

/*********************************************************************/

/*  Perform cloud tests. */
//  printf("\nPerforming cloud tests in landnight.c\n");

/*********************************************************************/    

//  Group 1 Tests
    
/*********************************************************************/     

/*  Co2 high cloud test. */

    if(rintf(m35) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(14, pxout.qabits);
      if(m35 >= nlco2[1]) (void) set_bit(14, pxout.testbits);
      conf = conf_test(m35, dlco2[0], nlco2[2], nlco2[3],
                                      nlco2[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("co2 %f %f %f %d %d\n", m35,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

/*  Water vapor high cloud test. */

    if(rintf(m27) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(15, pxout.qabits);
      if(m27 >= nlh20[1]) (void) set_bit(15, pxout.testbits);
      conf = conf_test(m27, nlh20[0], nlh20[2], nlh20[3],
                                      nlh20[1]);
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
          lst_thresh = nl_sfct1[0];
        }
        else if(df1 >= nl_df1[0] || (df1 < nl_df1[0] && (df2 <= nl_df2[0] || df2 >= nl_df2[1])) ) {
          lst_thresh = nl_sfct2[0];
        }
        else {
          lst_thresh = nl_sfct1[0];
        }

        if(df1 >= nl_df1[1]) {
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
//      printf("SFCT: %f %f %f %f %f %f \n", nl_sfct1[0], nl_sfct2[0], 
//        nl_df1[0], nl_df1[1], nl_df2[0], nl_df2[1]);
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
        btd_thr = nl11_12hi[0]; }
      else {
        btd_thr = thr;
      }
      
      locut = btd_thr;
      midpt = btd_thr - (0.3 * btd_thr);
      if(m31 < nlbt1[0]) {
        if(fabs(pxin.lat) <= nl_lat[0]) {
          hicut = midpt - 1.25;
        }
        else {
          a = (90.0 - fabs(pxin.lat)) / 60.0;
          hicut = -0.1 - (powf(a, 4) * 1.15);
        }
      }
      else {
    	hicut = midpt - 1.25;
      }
      
      if (m31_m32 <= midpt) (void) set_bit(18, pxout.testbits);
      
      conf = conf_test(m31_m32, locut, hicut, 1.0, midpt);
      cmin2 = min(cmin2, conf);
      ngtests[1]++;

//    printf("m31_32: %f %f %f %f %f %f %f %f %f %d %d\n", nlbt1[0], nl_lat[0],
//       m31_m32, btd_thr, locut, midpt, hicut, conf, cmin2, num_tests, ngtests[1]);

    }
    
/**********************************************************************/

/*  11 minus 4 micron BTDIF fog and low cloud test. */
      
/**********************************************************************/

    if (rintf(m31) != rintf(bad_data) && rintf(m22) != rintf(bad_data)
             && rintf(m32) != rintf(bad_data) ) {

      num_tests++;
      (void) set_bit(19, pxout.qabits);
      
      m31_m22 = m31 - m22;

//    Find if pixel is within Australian region.
      Australia = 0;
      if(pxin.lat <= -11.0 && pxin.lat > -40.0 && pxin.lon <= 155.0 && pxin.lon > 110.0) {
        Australia = 1;
      }

      if( (Australia == 0 && pxin.ndvibk < nl_ndvi[0] && pxin.coast == 0) ||
          (Australia == 1 && pxin.ndvibk < nl_ndvi_Aust[0] && pxin.coast == 0) ) {

        df1 = m31 - m32;
        (void) get_nl_thresholds(df1, &locut, &hicut, &midpt, &power);

        if(pxin.sh_lake) {
          locut = locut + 2.0;
          midpt = midpt + 2.0;
          hicut = hicut + 2.0;
        }

        if (m31_m22 <= midpt) (void) set_bit(19, pxout.testbits);
        conf = conf_test(m31_m22, locut, hicut, 1.0, midpt);

//      printf("11-4a: %d %f %f %f %f %f %f %f %f %d %d\n", Australia, pxin.ndvibk, nl_ndvi_Aust[0],
//             m31_m22, locut, midpt, hicut, conf, cmin2, num_tests, ngtests[1]);

      }
      else {

        thr = b1 + nl_intadj[0] + (b2 * pxin.tpw) + (b3 * pxin.tpw * pxin.tpw);
        hicutr[0] = thr - 1.0;
        hicutr[1] = thr + 1.0;
        midptr[0] = hicutr[0] - 0.50;
        midptr[1] = hicutr[1] + 0.50;
        locutr[0] = hicutr[0] - 1.0;
        locutr[1] = hicutr[1] + 1.0;

        if(pxin.sh_lake) {
          locutr[1] = locutr[1] + 1.0;
          locutr[0] = locutr[0] - 1.0;
          midptr[1] = midptr[1] + 1.0;
          midptr[0] = midptr[0] - 1.0;
          hicutr[1] = hicutr[1] + 1.0;
          hicutr[0] = hicutr[0] - 1.0;
        }
      
        if (m31_m22 >= midptr[0] && m31_m22 <= midptr[1]) (void) set_bit(19, pxout.testbits);
        conf = conf_test_dble(m31_m22, locutr, hicutr, 1.0, midptr);

//      printf("11-4b: %d %d %f %f %f %f %f %f %f %f %f %f %f %d %d\n", Australia, pxin.coast,
//            pxin.ndvibk, nl_ndvi[0], m31_m22, locutr[0], locutr[1], midptr[0], midptr[1],
//            hicutr[0], hicutr[1], conf, cmin2, num_tests, ngtests[1]);

      }

      cmin2 = min(cmin2, conf);
      ngtests[1]++;
    }

/**********************************************************************/
    
//  7.3-11um BTD test for thick, mid-level clouds.
    
/**********************************************************************/

    if (rintf(m31) != rintf(bad_data) && rintf(m28) != rintf(bad_data) &&
        rintf(m22) != rintf(bad_data) ) {

      m31_m22 = m31 - m22;
      if(m31_m22 <= nl_btd1[0]) {

        num_tests++;
        (void) set_bit(23, pxout.qabits);
        m28_m31 = m28 - m31;

        if ( m28_m31 <= nl7_11s[1] ) (void) set_bit(23, pxout.testbits);
        conf = conf_test(m28_m31, nl7_11s[0], nl7_11s[2], nl7_11s[3],
                         nl7_11s[1]);
        cmin2 = min(cmin2, conf);
        ngtests[1]++;
//      printf("7.3-11: %f %f %f %f %f %f %f %f %f %d %d \n", nl_btd1[0],
//        m28, m31, m28_m31,
//        nl7_11s[0], nl7_11s[1], nl7_11s[2], conf, cmin2, num_tests,
//        ngtests[1]);

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
      if (m22_m32 <= nl4_12hi[1]) (void) set_bit(17, pxout.testbits);
      conf = conf_test(m22_m32, nl4_12hi[0], nl4_12hi[2], nl4_12hi[3],
                      nl4_12hi[1]);
      cmin5 = min(cmin5, conf);
      ngtests[2]++;
//    printf("4-12: %f %f %f %f %f %f %d %d \n", m22_m32, nl4_12hi[0],
//    		  nl4_12hi[1], nl4_12hi[2], cmin5, conf, num_tests, ngtests[2]);
      
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
