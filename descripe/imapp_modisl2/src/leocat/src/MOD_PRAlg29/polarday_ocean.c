/*
!C ********************************************************************
!Description:

  Integer function polarday_ocean.c
  Computes initial clear sky confidence for polar daytime water surfaces.
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
  R. Frey           05/2007


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


int polarday_ocean()


/*********************************************************************/

{

/*  Declarations */

    extern int check_bits(int, unsigned char[]);
    extern float cithr(int, float, float);
    extern void clear_bit(int, unsigned char[]);
    extern float conf_test(float, float, float, float, float);
    extern float conf_test_dble(float, float[], float[], float, float[]);
    extern void get_sg_thresholds(float, float*);
    extern void set_bit(int, unsigned char[]);
    
    int groups;
    int k;
    int ngtests[num_test_groups];
    int num_tests;
    int return_code;

    float a1=1.8860 , a2=0.9380, a3=0.1280, a4=1.094;
    float btd_thr;
    float cmin1, cmin2, cmin3, cmin4;
    float conf;
    float cosvza;
    float dosgref2[4];
    float fac;
    float m01, m02, m22, m26, m27, m29, m31, m32;
    float m29_m31, m31_m32, m31_m22;
    float m31c, m32c, m31c_m32c;
    float midpt, locut, hicut;
    float midptr[2], locutr[2], hicutr[2];
    float modsst;
    float mu;
    float pre_confidence;
    float schi;
    float sfcdif;
    float sfct;
    float sstc;
    float thr;
    float vrat;
    
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
    m29 = pxin.rad[28];
    m31 = pxin.rad[30];
    m32 = pxin.rad[31];

/*********************************************************************/

/*  Perform cloud tests. */
//  printf("\nPerforming cloud tests in polarday_ocean.c\n");

/*********************************************************************/

/*  Group 1 tests. */

/*  11 micron brightness temperature threshold test. */

    if(rintf(m31) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(13, pxout.qabits);
      if(m31 >= pdobt11[1]) (void) set_bit(13, pxout.testbits);
      conf = conf_test(m31, pdobt11[0], pdobt11[2], pdobt11[3],
                                  pdobt11[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("bt11 %f %f %f %d %d\n", m31,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

/*  Water vapor high cloud test. */

    if(rintf(m27) != rintf(bad_data)) {
      num_tests++;
      (void) set_bit(15, pxout.qabits);
      if(m27 >= pdoh20[1]) (void) set_bit(15, pxout.testbits);
      conf = conf_test(m27, pdoh20[0], pdoh20[2], pdoh20[3],
                                  pdoh20[1]);
      cmin1 = min(cmin1, conf);
      ngtests[0]++;
//    printf("h20 %f %f %f %d %d\n", m27,conf,cmin1,ngtests[0],num_tests);
    }

/*********************************************************************/

/*  SST Test. */

    if(rintf(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data)) {
      
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
        
        if(sfcdif < pdosst[1]) (void) set_bit(27, pxout.testbits);
        conf = conf_test(sfcdif, pdosst[0], pdosst[2], pdosst[3], pdosst[1]);
        
        cmin1 = min(cmin1, conf);
        ngtests[0]++;
        
//      printf("SST %f %f %f %f %f\n", m31c_m32c,sstc,mu,modsst,sfct);
//      printf("SST %f %f %f %f %f %f %f %d %d\n", sfcdif, pdosst[0], pdosst[1],
//           pdosst[2], pdosst[3], conf, cmin1, ngtests[0], num_tests);
        
      }

    }

/*********************************************************************/

/*  Group 2 tests. */
    
/*********************************************************************/

/*  8.6-11 micron BTD test. */

    if(rintf(m31) != rintf(bad_data) && rintf(m29) != rintf(bad_data)) {

      num_tests++;
      (void) set_bit(24, pxout.qabits);

      m29_m31 = m29 - m31;

      if(m29_m31 < pdo11_86[1]) (void) set_bit(24, pxout.testbits);

      conf = conf_test(m29_m31,pdo11_86[0], pdo11_86[2], pdo11_86[3], pdo11_86[1]);
      cmin2 = min(cmin2, conf);
      ngtests[1]++;
//    printf("11-8.6 %f %f %f %f %f %f %f %d %d\n", m29_m31,pdo11_86[0],pdo11_86[1],
//                pdo11_86[2],pdo11_86[3],conf,cmin2,num_tests,ngtests[1]);
    
    }

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
        btd_thr = do11_12hi[0]; }
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
//                btd_thr, hicut, conf, cmin2, num_tests, ngtests[1]);

    }
    
/**********************************************************************/

/*    11 minus 4 micron BTDIF fog and low cloud test. */

    if (pxin.visusd && (pxin.sunglint == 0)) {
      if (rintf(m31) != rintf(bad_data) && rintf(m22) != rintf(bad_data)) {

        num_tests++;
        (void) set_bit(19, pxout.qabits);
        m31_m22 = m31 - m22;

        if (m31_m22 >= pdo11_4lo[1]) (void) set_bit(19, pxout.testbits);
        conf = conf_test(m31_m22, pdo11_4lo[0], pdo11_4lo[2], 1.0,
                         pdo11_4lo[1]);
        cmin2 = min(cmin2, conf);
        ngtests[1]++;
//      printf("11-4: %f %f %f %f %f %f %d %d\n", m31_m22, pdo11_4lo[0],
//         pdo11_4lo[1], pdo11_4lo[2], conf, cmin2, num_tests, ngtests[1]);

      }
    }

/**********************************************************************/

/*  Group 3 tests. */

/**********************************************************************/
    
/*  NIR reflectance threshold test. */

    if (pxin.visusd) {
      if (rintf(m02) != rintf(bad_data)) {

        num_tests++;
        (void) set_bit(20, pxout.qabits);

/*      Different thresholds apply in sun-glint conditions. */
        if(pxin.sunglint) {
          (void) get_sg_thresholds(pxin.refang, dosgref2);
          locut = dosgref2[0];
          midpt = dosgref2[1];
          hicut = dosgref2[2]; }
        else {
          locut = pdoref2[0];
          midpt = pdoref2[1];
          hicut = pdoref2[2];
        }

        if (m02 <= midpt) (void) set_bit(20, pxout.testbits);
        conf = conf_test(m02, locut, hicut, 1.0, midpt);
        cmin3 = min(cmin3, conf);
        ngtests[2]++;
//      printf("NIR: %f %f %f %f %f %f %d %d\n", m02, locut, midpt,
//                hicut, conf, cmin3, num_tests, ngtests[2]);

      }
    }

/**********************************************************************/

/*  Vis/NIR ratio test. */

    if( pxin.visusd && rintf(m01) != rintf(bad_data) &&
                rintf(m02) != rintf(bad_data)  ) {

      num_tests++;
      (void) set_bit(21, pxout.qabits);

      if(pxin.sunglint) {
        locutr[0] = snglntvcl[0];
        locutr[1] = snglntvcl[1];
        hicutr[0] = snglntvch[0];
        hicutr[1] = snglntvch[1];
        midptr[0] = snglntv[0];
        midptr[1] = snglntv[1]; }
      else {
        locutr[0] = pdovratlo[0];
        locutr[1] = pdovrathi[0];
        hicutr[0] = pdovratlo[2];
        hicutr[1] = pdovrathi[2];
        midptr[0] = pdovratlo[1];
        midptr[1] = pdovrathi[1];
      }

      vrat = m02 / m01;
      if(vrat < midptr[0] || vrat > midptr[1])
                  (void) set_bit(21, pxout.testbits);
      conf = conf_test_dble(vrat, locutr, hicutr, 1.0, midptr);
      cmin3 = min(cmin3, conf);
      ngtests[2]++;
//    printf("vrat %f %f %f %f %f %d %d\n", m01,m02,vrat,conf,cmin3,
//                  num_tests,ngtests[2]);
    }

/*********************************************************************/

/*  Group 4 tests. */

/*  Near-infrared (1.38 micron) high cloud test. */

    if (pxin.visusd) {
      if (rintf(m26) != rintf(bad_data)) {

        num_tests++;
        (void) set_bit(16, pxout.qabits);

        if (m26 <= pdoref3[1]) (void) set_bit(16, pxout.testbits);
        conf = conf_test(m26, pdoref3[0], pdoref3[2], 1.0, pdoref3[1]);
        cmin4 = min(cmin4, conf);
        ngtests[3]++;
//      printf("1.38: %f %f %f %f %f %f %d %d\n", m26,pdoref3[0],
//             pdoref3[1],pdoref3[2],conf,cmin4,num_tests,ngtests[3]);

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
