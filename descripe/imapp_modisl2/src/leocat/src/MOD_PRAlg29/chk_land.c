/*
************************************************************************
!F77

!Description:
   Perform final clear-sky confidence check on daytime land pixels.
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


float chk_land()


{


/**********************************************************************/

//  Declarations

    extern int check_bits(int, unsigned char[]);
    extern void set_bit(int, unsigned char[]);

    int bit_test;
    int bit_test2;
    int qa_test;
    int irclr;

    float m04, m05, m20, m22, m31;
    float m5_4_thr;
    float m5_4;
    float md1, md2;
    float hds11[3];
    float conf;

/**********************************************************************/

//	Initializations

//  printf("Executing chk_land \n");
    conf = pxout.intermediate_conf;

    m05 = pxin.rad[4];
    m04 = pxin.rad[3];
    m31 = pxin.rad[30];
    m22 = pxin.rad[21];
    m20 = pxin.rad[19];

    hds11[0] = 0.0;
    hds11[1] = 0.0;
    hds11[2] = 0.0;

    m5_4_thr = 0.0;
    m5_4 = 0.0;
    md1 = 0.0;
    md2 = 0.0;

/**********************************************************************/

//  Check IR clear sky tests.
    irclr = 1;
    bit_test = check_bits(14, pxout.testbits);
    qa_test = check_bits(14, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    bit_test = check_bits(15, pxout.testbits);
    qa_test = check_bits(15, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    bit_test = check_bits(16, pxout.testbits);
    qa_test = check_bits(16, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    bit_test = check_bits(18, pxout.testbits);
    qa_test = check_bits(18, pxout.qabits);
    if( qa_test && !bit_test ) irclr = 0;
    if(pxin.desert == 0) {
      bit_test2 = 1;
      bit_test = check_bits(20, pxout.testbits);
      qa_test = check_bits(20, pxout.qabits);
      if( qa_test && !bit_test ) bit_test2 = 0;
      bit_test = check_bits(21, pxout.testbits);
      qa_test = check_bits(21, pxout.qabits);
      if( qa_test && !bit_test ) bit_test2 = 0;
      bit_test = check_bits(19, pxout.testbits);
      qa_test = check_bits(19, pxout.qabits);
      if( qa_test && !bit_test && !bit_test2) irclr = 0;
    }

/**********************************************************************/

    if(irclr) {

      if (rintf(m31) != rintf(bad_data)) {

//      Get elevation-adjusted 11 micron brightness temperature threshold.

        (void) set_bit(26, pxout.qabits);

        if(pxin.eco_type == 8) {
          hds11[0] = ldsbt11bd[0] - pxin.tbadj;
          hds11[1] = ldsbt11bd[1] - pxin.tbadj;
          hds11[2] = ldsbt11bd[2] - pxin.tbadj;
        }
        else {
          hds11[0] = ldsbt11[0] - pxin.tbadj;
          hds11[1] = ldsbt11[1] - pxin.tbadj;
          hds11[2] = ldsbt11[2] - pxin.tbadj;
        }

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

//    printf("chk_land m31: %d %f %f %f %f %f %f \n", pxin.eco_type, pxin.tbadj,
//    		    hds11[0], hds11[1], hds11[2], m31, conf);

/**********************************************************************/

//    Perform additional clear-sky spectral tests.

      if (conf <= 0.95) {

        if (rintf(m20) != rintf(bad_data) &&
            rintf(m22) != rintf(bad_data) &&
            rintf(m31) != rintf(bad_data) &&
            rintf(m05) != rintf(bad_data) &&
            rintf(m04) != rintf(bad_data)) {

          if(pxin.desert) {
            m5_4_thr = ldsr5_4_thr[0];
          }
          else {
            m5_4_thr = ldr5_4_thr[0];
          }

          m5_4 = m05 / m04;
          md1 = m20 - m22;
          md2 = m22 - m31;

          if (md1 < ld20m22[0] && md2 < ld22m31[0] && m5_4 > m5_4_thr) {
            conf = 0.96;
            (void) set_bit(26, pxout.testbits);
          }

        }

      }

//    printf("chk_land spec: %d %f %f %f %f %f %f %f \n", pxin.desert,
// 		    m5_4_thr, m5_4, md1, md2, ld20m22[0], ld22m31[0], conf);

    }

/**********************************************************************/

    return conf;

}
