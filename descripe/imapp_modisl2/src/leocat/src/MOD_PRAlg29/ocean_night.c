/*
!C ********************************************************************
!Description:

  Integer function ocean_night.c
  Computes initial clear sky confidence for nighttime water surfaces.
  Called from water_night.c

!Input arguments:
  none

  Inputs stored in structure variable 'pxin' of type "pixel_in" defined
  in pixel.h.
  Test thresholds accessed through thresholds.h.

!Output arguments:
  none

  int return_code   successful completion is zero, otherwise non-zero
                    (returned through function call)

  Output confidence and associated products stored in structure
  variable 'pxout' of type "pixel_out" defined in pixel.h

!Revision History:
  R. Frey           07/2007


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


int ocean_night()


/*********************************************************************/

{

/*  Declarations */

    extern int check_bits(int, unsigned char[]);
    extern void chk_spatial2();
    extern float cithr(int, float, float);
    extern void clear_bit(int, unsigned char[]);
    extern float conf_test(float, float, float, float, float);
    extern float conf_test_dble(float, float[], float[], float, float[]);
    extern void set_bit(int, unsigned char[]);
    extern float trispc(float);
    
    int groups;
    int k;
    int ngtests[num_test_groups];
    int num_tests;
    int return_code;

    float a1=1.8860 , a2=0.9380, a3=0.1280, a4=1.094;
    float b1=0.5972, b2=-0.2460, b3=-0.1501;
    float btd_thr;
    float cmin1, cmin2;
    float conf;
    float cosvza;
    float fac;
    float m22, m27, m28, m29, m31, m32, m35; 
    float m29_m28, m29_m31, m31_m32, m31_m22;
    float m31c, m32c, m31c_m32c;
    float locut, hicut;
    float midptr[2], locutr[2], hicutr[2];
    float modsst;
    float mu;    
    float np;
    float pre_confidence;
    float schi;
    float sfcdif;
    float sfct;
    float sstc;
    float thr;
    
/*********************************************************************/

/*  Initializations */
    return_code = 0;

/*  Counter for total number tests performed. */
    num_tests = 0;

/*  Test group confidences. */
    cmin1 = 1.0;
    cmin2 = 1.0;
//  cmin3 = 1.0;
//  cmin4 = 1.0;

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
    m29 = pxin.rad[28];
    m31 = pxin.rad[30];
    m32 = pxin.rad[31];
    m35 = pxin.rad[34];

/*********************************************************************/

/*  Perform cloud tests. */
//  printf("\nPerforming cloud tests in ocean_night.c\n");

/*********************************************************************/

/*  Group 1 tests. */

/*  11 micron brightness temperature threshold test. */

    if(rintf(m31) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(13, pxout.qabits);
      if(m31 >= nobt11[1]) (void) set_bit(13, pxout.testbits);
      conf = conf_test(m31, nobt11[0], nobt11[2], nobt11[3],
                                  nobt11[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("bt11 %f %f %f %d %d\n", m31,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

/*  Co2 high cloud test. */

    if(rintf(m35) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(14, pxout.qabits);
      if(m35 >= noco2[1]) (void) set_bit(14, pxout.testbits);
      conf = conf_test(m35, noco2[0], noco2[2], noco2[3],
                                  noco2[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("co2 %f %f %f %d %d\n", m35,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

/*  Water vapor high cloud test. */

    if(rintf(m27) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(15, pxout.qabits);
      if(m27 >= noh20[1]) (void) set_bit(15, pxout.testbits);
      conf = conf_test(m27, noh20[0], noh20[2], noh20[3],
                                  noh20[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("h20 %f %f %f %d %d\n", m27,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

/*  SST Test. */

    if(rint(m31) != rint(bad_data) && rint(m32) != rint(bad_data)) {
      
//    if(pxin.sfctmp > 0.0 && pxin.sfctmp < 350.0 && pxin.sfct_flag == 0) {
      if(pxin.sfctmp > 0.0 && pxin.sfctmp < 350.0) {

        num_tests++;
        (void) set_bit(27, pxout.qabits);

//      Limit ancillary surface temperature to 295K.
        sfct = pxin.sfctmp;
//      if(sfct > 295.0) sfct = 295.0;
        m31c = m31 - 273.16;
        m32c = m32 - 273.16;      
        m31c_m32c = m31c - m32c;
        sstc = sfct - 273.16;
        mu = cos(pxin.vza * 3.14159 / 180.0);        

        modsst = 273.16 + a1 + a2*m31c + a3*m31c_m32c*sstc + a4*m31c_m32c*((1.0/mu)-1.0);      
        
        sfcdif = sfct - modsst;        
        
        if(sfcdif < nosst[1]) (void) set_bit(27, pxout.testbits);
        conf = conf_test(sfcdif, nosst[0], nosst[2], nosst[3], nosst[1]);
        
        cmin1 = min(cmin1, conf);
        
        ngtests[0]++;
        
//      printf("SST %f %f %f %f %f\n", m31c_m32c,sstc,mu,modsst,sfct);
//      printf("SST %f %f %f %f %f %f %f %d %d\n", nosst[0], nosst[1],
//        nosst[2], nosst[3], sfcdif, conf, cmin1, ngtests[0], num_tests);
        
      }

    }

/*********************************************************************/

/*  Group 2 tests. */

/*  11-8.6 micron BTD test. */

    if(rint(m31) != rint(bad_data) && rint(m32) != rint(bad_data) &&
       rint(m29) != rint(bad_data)) {

      num_tests++;
      (void) set_bit(24, pxout.qabits);

      m29_m31 = m29 - m31;

      if(m29_m31 < no11_86[1]) (void) set_bit(24, pxout.testbits);
      conf = conf_test(m29_m31, no11_86[0], no11_86[2], no11_86[3], no11_86[1]);
      
      cmin2 = min(cmin2, conf);
      
      ngtests[1]++;
      
//    printf("11-8.6 %f %f %f %f %f %f %d %d\n", m29_m31,no11_86[0],no11_86[1],
//              no11_86[2],conf,cmin2,num_tests,ngtests[1]);
    
    }

/*********************************************************************/

/*  11-12 micron BTD test for transmissive cirrus clouds. */

    if(rint(m31) != rint(bad_data) && rint(m32) != rint(bad_data) &&
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
      if(thr < 0.1 || fabsf(schi - 99.0) < 0.0001) {
        btd_thr = no11_12hi[0];
      }
      else {
        btd_thr = thr;
      }
      
      if (m31_m32 < btd_thr) {
    	(void) set_bit(18, pxout.testbits);
      }
    
      locut = btd_thr + (0.3 * btd_thr);
      hicut = btd_thr - 1.25;
      conf = conf_test(m31_m32,locut,hicut,1.0,btd_thr);
      
      cmin2 = min(cmin2, conf);
      
      ngtests[1]++;

//    printf("m31_32: %f %f %f %f %f %f %d %d\n", m31_m32, locut,
//              btd_thr, hicut, conf, cmin2, num_tests, ngtests[1]);

    }
    
/**********************************************************************/

/*  11 minus 4 micron BTDIF cloud test. */

    if (rint(m31) != rint(bad_data) && rint(m22) != rint(bad_data)) {

      num_tests++;
      (void) set_bit(19, pxout.qabits);
      m31_m22 = m31 - m22;

      thr = b1 + no_intadj[0] + (b2 * pxin.tpw) + (b3 * pxin.tpw * pxin.tpw);
      hicutr[0] = thr - 1.5;
      hicutr[1] = thr + 1.5;
      midptr[0] = hicutr[0] - 2.5;
      midptr[1] = hicutr[1] + 2.5;
      locutr[0] = hicutr[0] - 4.0;
      locutr[1] = hicutr[1] + 4.0;
      
      if (m31_m22 >= midptr[0] && m31_m22 <= midptr[1]) (void) set_bit(19, pxout.testbits);
      conf = conf_test_dble(m31_m22, locutr, hicutr, 1.0, midptr);
      
      cmin2 = min(cmin2, conf);
      
      ngtests[1]++;
      
//    printf("11-4: %f %f %f %f %f %f %f %f %f %f %f %f %d %d\n", m31_m22, pxin.tpw,
//           no_intadj[0], thr, locutr[0], midptr[0], hicutr[0], locutr[1],
//           midptr[1], hicutr[1], conf, cmin2, num_tests, ngtests[1]);
      
    }

/**********************************************************************/

/*  Water vapor cloud test. */

    if (rint(m28) != rint(bad_data) && rint(m29) != rint(bad_data)) {

      num_tests++;
      (void) set_bit(29, pxout.qabits);
      m29_m28 = m29 - m28;

      if (m29_m28 > no86_73[1]) (void) set_bit(29, pxout.testbits);
      conf = conf_test(m29_m28, no86_73[0], no86_73[2], 1.0,
                       no86_73[1]);
      cmin2 = min(cmin2, conf);
      
      ngtests[1]++;
      
//    printf("8.6-7.3: %f %f %f %f %f %f %d %d\n", m29_m28, no86_73[0],
//       no86_73[1], no86_73[2], conf, cmin2, num_tests, ngtests[1]);

    }

/*********************************************************************/

/*  11 micron BT variability test. */

    if (pxin.uniform == 1) {

      num_tests++;
      (void) set_bit(30, pxout.qabits);

//    Find number of surrounding pixels that satify 11 um variability
//    (BTD) constraint for clear sky.
      (void) chk_spatial2();
      np = rg_var.num_small_diffs * 1.0;

      if (np > no_11var[1]) (void) set_bit(30, pxout.testbits);
      conf = conf_test(np, no_11var[0], no_11var[2], 1.0,
                       no_11var[1]);
      
      cmin2 = min(cmin2, conf);
      
      ngtests[1]++;
      
//    printf("11 um var: %f %f %f %f %f %f %d %d\n", np, no_11var[0],
//       no_11var[1], no_11var[2], conf, cmin2, num_tests, ngtests[2]);

    }

/*********************************************************************/
    
/*  11 minus 4 micron BTDIF oceanic stratus (low emissivity water cloud) test. */

    if (rint(m31) != rint(bad_data) && rint(m22) != rint(bad_data)) {

      num_tests++;
      (void) set_bit(31, pxout.qabits);
      m31_m22 = m31 - m22;
          
      if (m31_m22 <= no11_4lo[1]) (void) set_bit(31, pxout.testbits);
      conf = conf_test(m31_m22, no11_4lo[0], no11_4lo[2], 1.0,
                       no11_4lo[1]);
      
      cmin2 = min(cmin2, conf);
      
      ngtests[1]++;
              
//    printf("11-4: %f %f %f %f %f %f %d %d\n", m31_m22, no11_4lo[0],
//            no11_4lo[1], no11_4lo[2], conf, cmin2, num_tests, ngtests[2]);

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
//  pre_confidence = cmin1 * cmin3 * cmin4;

/*  Find the number of test groups used. */
    k = 0;
    while(k < num_test_groups) {
      if(ngtests[k] > 0) groups++;
      k++;
    }

/*  Get power for confidence calculation. */
    if(groups > 0) fac = 1.0 / groups;

//  printf("pre_conf, groups, fac: %f %d %f\n", pre_confidence, groups,
//            fac);

/*  Get initial confidence of clear sky. */
    pxout.init_conf = powf(pre_confidence, fac);

    return_code = num_tests;
    return return_code;

}
