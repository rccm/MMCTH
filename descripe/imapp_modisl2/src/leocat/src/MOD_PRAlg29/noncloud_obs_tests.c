/*
************************************************************************
!F77

!Description:

   Performs non-cloud obstruction tests.  These results of these tests
   do not contribute to the final confidence of clear sky but are reported
   in bits 8 (smoke) and 28 (suspended dust).

!Input parameters:
   none

   Inputs stored in structure variable 'pxin' of type "pixel_in" defined
   in pixel.h.

!Output Parameters:
   none

   Output values stored in structure variable 'pxout' of type "pixel.out".

!Revision History:
    Converted to C              11/07      R. Frey
    Added "GOES-R" dust tests   12/08      R. Frey

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


void noncloud_obs_tests()


{

/***********************************************************************/

//  Declarations

    extern void set_bit(int, unsigned char[]);
    extern void clear_bit(int, unsigned char[]);
    extern int check_bits(int, unsigned char[]);
    extern float get_regstd(int);

    int bit_test;
    int qa_test;
    int irclr;
    int fire;
    int band;
    int thk_smoke;
    int thk_smk_water;
    int dust;

    float a;
    float corr;
    float coef;
    float m01, m02, m03, m04, m05, m07, m20, m21, m26, m29, m31, m32;
    float midpt;
    float ndvi;
    float newrat1, newrat2;
    float rat1, rat2, rat3, rat4, rat35, rat75;
    float sfcdif;
    float smkrat;
    float sst_thrsh;
    float stdev;
    float tdif, tdiff1, tdiff2, tdiff3;
    float vrat;
    float mx_vza = 65.49;

/**********************************************************************/

//  Initializations

//  printf("Executing noncloud_obs_tests \n");

    m01 = pxin.rad[0];
    m02 = pxin.rad[1];
    m03 = pxin.rad[2];
    m04 = pxin.rad[3];
    m05 = pxin.rad[4];
    m07 = pxin.rad[6];
    m20 = pxin.rad[19];
    m21 = pxin.rad[20];
    m26 = pxin.rad[25];
    m29 = pxin.rad[28];
    m31 = pxin.rad[30];
    m32 = pxin.rad[31];

    dust = 0;
    fire = 0;
    thk_smoke = 0;
    thk_smk_water = 0;

/**********************************************************************/

//  Thick smoke test (land).

    if(pxin.land && pxin.day) {

//    Check IR clear sky tests.
      irclr = 1;
      bit_test = check_bits(15, pxout.testbits);
      qa_test = check_bits(15, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
      bit_test = check_bits(16, pxout.testbits);
      qa_test = check_bits(16, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
//    if( qa_test && !bit_test ) irclr = 0;
      bit_test = check_bits(18, pxout.testbits);
      qa_test = check_bits(18, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
      bit_test = check_bits(19, pxout.testbits);
      qa_test = check_bits(19, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;

//    Do not perform fire or smoke tests if tests reported in bits
//    15, 16, 18, or 19 reported cloud.

      if(irclr) {

//      Set bit indicating smoke test attempt.
        (void) set_bit(8, pxout.qabits);

//      Check for fires (hot spots).

        if(m21 != rintf(bad_data) && m31 != rintf(bad_data)) {
          tdif = m21 - m31;
    	  if (m21 > nc_bt37[0] && tdif > nc37_11[0]) {
    	    fire = 1;
    	  }
        }

//      Test for thick smoke.

        if (rintf(m01) != rintf(bad_data) && rintf(m07) != rintf(bad_data) &&
     	    rintf(m03) != rintf(bad_data) && rintf(m02) != rintf(bad_data)) {

          coef = 0.0;
          smkrat = 0.0;
          vrat = 0.0;
          stdev = 0.0;

    	  if( (m07 * 100.0) < nc21[0]) {

    	    coef = 6.0 + (m07 * 100.0);
    	    if( (m01 * 100.0) > coef) {

    	      smkrat = m03 / m01;
    	      if(smkrat >= ncrat[0]) {

    	        vrat = m02 / m01;
    	        if(vrat >= ncvrat[0]) {

    	          band = 1;
    	          stdev = get_regstd(band);
    	          if(stdev <= ncsig[0]) {
                    thk_smoke = 1;
    	          }
    	        }
    	      }
    	    }
    	  }
        }
      }

      if(!thk_smoke && !fire) (void) set_bit(8, pxout.testbits);

//    printf("Smoke test: %d %f %d %f %f %f %f %f %f %d \n",
//    		  irclr, tdif, fire, m07, coef, m01, smkrat, vrat, stdev, thk_smoke);

    }

//  Thick smoke/aerosol test (water).

    if(pxin.water && pxin.day) {

      if(rintf(m03) != rintf(bad_data) && rintf(m05) != rintf(bad_data) &&
         rintf(m02) != rintf(bad_data) ) {

//      Check IR clear sky tests.
        irclr = 1;
        bit_test = check_bits(15, pxout.testbits);
        qa_test = check_bits(15, pxout.qabits);
        if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
        bit_test = check_bits(9, pxout.testbits);
        qa_test = check_bits(9, pxout.qabits);
        if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
        bit_test = check_bits(16, pxout.testbits);
        qa_test = check_bits(16, pxout.qabits);
        if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
        bit_test = check_bits(18, pxout.testbits);
        qa_test = check_bits(18, pxout.qabits);
        if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
        bit_test = check_bits(27, pxout.testbits);
        qa_test = check_bits(27, pxout.qabits);
        if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
//      This check (11-3.9 micron) also screens out sun-glint.
        bit_test = check_bits(19, pxout.testbits);
        qa_test = check_bits(19, pxout.qabits);
        if( (qa_test && !bit_test) || !qa_test ) irclr = 0;

//      Do not perform tests if cloud tests reported in bits
//      9, 15, 16, 18, 19, or 27 indicated cloud.

        if(irclr == 1) {

//        Set bit indicating smoke/pollution test attempt.
          (void) set_bit(8, pxout.qabits);

//        Initialize to "no smoke/aerosols". Bit is cleared if smoke is found.
          (void) set_bit(8, pxout.testbits);

//        Compute band 2 standard deviation over 3x3 pixel area.
          band = 2;
          stdev = get_regstd(band);

          rat35 = m03/m05;
          rat75 = m07/m05;

          if(stdev < 0.003) {

            if(rat35 >= 5.0 && m03 >= 0.12 && m05 >= 0.022 && m05 < 0.050 && rat75 < 0.5) {
              (void) clear_bit(8, pxout.testbits);
              thk_smk_water = 1;
            }

//          Check for positive NIR reflectance cloud test (0.86 um).
            bit_test = check_bits(20, pxout.testbits);
            qa_test = check_bits(20, pxout.qabits);
            if( qa_test && !bit_test) {
              if(rat35 >= 2.5) {
                (void) clear_bit(8, pxout.testbits);
                thk_smk_water = 1;
              }
            }


          }

          bit_test = check_bits(20, pxout.testbits);
          qa_test = check_bits(20, pxout.qabits);
          if( qa_test && !bit_test) {
            if(rat35 >= 2.5 && rat75 < 0.3) {
              (void) clear_bit(8, pxout.testbits);
              thk_smk_water = 1;
            }
          }

//        printf("smoke test: %f %f %f %f %d %d\n", stdev, rat35, m03, m05, qa_test, bit_test);

        }


      }
    }

/**********************************************************************/

//  Perform dust tests.

/**********************************************************************/

    (void) set_bit(28, pxout.testbits);

    if (pxin.land && pxin.day) {

//    Check clear sky tests.
      irclr = 1;
      bit_test = check_bits(15, pxout.testbits);
      qa_test = check_bits(15, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
      bit_test = check_bits(18, pxout.testbits);
      qa_test = check_bits(18, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;

//    Do not perform dust tests if tests reported in bits
//    15 or 18 reported cloud.

      if(irclr) {

        if (rintf(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data) &&
	        rintf(m20) != rintf(bad_data) && rintf(m01) != rintf(bad_data) &&
	        rintf(m03) != rintf(bad_data) && rintf(m02) != rintf(bad_data) &&
	        rintf(m26) != rintf(bad_data) && rintf(m05) != rintf(bad_data) ) {

//        Set bit to indicate that dust test was attempted.
          (void) set_bit(28, pxout.qabits);

          tdiff1 = m31 - m32;
          tdiff2 = m20 - m31;
          rat1 = (m05 - m02) / (m05 + m02);
          rat2 = (m01 - m03) / (m01 + m03);
          newrat2 = powf(rat2,2) / powf(m01,2);
          ndvi = (m02 - m01) / (m02 + m01);
          newrat1 = powf(ndvi,2) / powf(m01,2);

          if(tdiff1 <= ncd_lddf1[0] && tdiff2 >= ncd_lddf2[0] && m26 < ncd_ldm26[0]) {

            if(newrat2 > ncd_ldrat2[0] && newrat1 < ncd_ldrat1[0]) {
              (void) clear_bit(28, pxout.testbits);
              dust = 1;
            }
            else if(tdiff2 >= ncd_lddf2[1]) {
              (void) clear_bit(28, pxout.testbits);
              dust = 1;
            }

          }
//        Test for very thick dust.
          else if(tdiff1 < ncd_lddf1[1] && tdiff2 >= ncd_lddf2[2] && m26 < ncd_ldm26[1] && newrat1 < ncd_ldrat1[1]) {

            (void) clear_bit(28, pxout.testbits);
            dust = 1;

          }

        }
      }


//    printf("Dust test for land, day: %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n",
//    		irclr, dust, m31, m32, m20, m26, m01, m02, m03, m05, tdiff1, tdiff2, rat1, rat2,
//    		newrat1, newrat2);

    }

/**********************************************************************/

    else if(pxin.water && pxin.day) {

//    Day water dust tests.

//    Check clear sky tests.
      irclr = 1;
      bit_test = check_bits(15, pxout.testbits);
      qa_test = check_bits(15, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
      bit_test = check_bits(18, pxout.testbits);
      qa_test = check_bits(18, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;

//    Do not perform dust tests if tests reported in bits
//    15 or 18 reported cloud (6.7 um high cloud and 11-12um BTD thin cirrus tests).

      if(irclr) {

        if (rintf(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data) &&
  	        rintf(m20) != rintf(bad_data) && rintf(m01) != rintf(bad_data) &&
  	        rintf(m03) != rintf(bad_data) && rintf(m02) != rintf(bad_data) &&
  	        pxin.sfctmp > 0.0 && pxin.sfctmp < 350.0) {

//        Set bit to indicate that dust test was attempted.
          (void) set_bit(28, pxout.qabits);

//        Perform dust test.

//        SST test to screen out low water clouds.
          if(pxin.sh_ocean) {
    	    sst_thrsh = ncd_wdsst[0];
          }
    	  else {
    	    sst_thrsh = ncd_wdsst[1];
    	  }

    	  tdiff1 = m31 - m32;
    	  if(tdiff1 >= ncd_wddf1[0]) {
    	    midpt = sst_thrsh + (2.0 * rintf(tdiff1));
    	  }
    	  else {
    	    midpt = sst_thrsh;
    	  }

    	  a = pxin.vza / max_vza;
    	  corr = ( powf(a,4) ) * 3.0;
    	  midpt = midpt + corr;

    	  sfcdif = pxin.sfctmp - m31;

    	  if( sfcdif < midpt ) {

//          Check band 2 (0.86 um) reflectance spatial variability.
    	    band = 2;
    	    stdev = get_regstd(band);

    	    if( rintf(rg_var.mean) != rint(bad_data) && stdev <= ncd_wdsig2[0]) {

//            Gross blue channel (0.47 um) test.
    	      if(m03 <= ncd_wdm03[0]) {

//              Gross blue/vis ratio (0.47/0.65 um) test.
    	        rat4 = m03 / m01;
    	        if(rat4 < ncd_wdrat4[0]) {

    	          tdiff2 = m20 - m31;
    	          rat3 = (m02 - m01) / (m02 + m01);

//                Gross 3.9-11 um BTD test.
    	          if(tdiff2 > ncd_wddf2[0] && tdiff2 <= ncd_wddf2[1]) {

//                  11-12 um and NDVI tests.
    	            if(tdiff1 < ncd_wddf1[1] && rat3 <= ncd_wdrat3[0] && rat3 >= ncd_wdrat3[1]) {
    	              (void)clear_bit(28, pxout.testbits);
    	              dust = 1;
    	            }

//                  Blue/vis ratio (0.47/0.65 um) test.
    	            if(rat4 < ncd_wdrat4[1] ) {
    	              (void)clear_bit(28, pxout.testbits);
    	              dust = 1;
    	            }

//                  3.9-11 um BTD test.
    	            if(tdiff2 > ncd_wddf2[2] && tdiff1 < ncd_wddf1[2]) {
    	              (void)clear_bit(28, pxout.testbits);
    	              dust = 1;
    	            }

//    	            3.9-11 um BTD test for very thick (Asian?) dust.
    	            else if(tdiff2 > ncd_wddf2[3]) {

//                    11-12 um and NDVI tests.
    	              if(tdiff1 < ncd_wddf1[3] && rat3 >= ncd_wdrat3[2] && rat3 <= ncd_wdrat3[3]) {
    	                (void)clear_bit(28, pxout.testbits);
    	                dust = 1;
    	              }
    	            }
    	          }

    	        }
    	      }
//            printf("Dust tests for water day: %d %f %f %f %f \n", dust, tdiff2, rat3, tdiff1, rat4);
    	    }

    	  }
          if( dust == 0 ) {

            float m14;
            float m43;
            m14 = m01 / m04;
            m43 = m04 / m03;
            if(m14 >= 1.15 && m43 >= 1.15) {
              dust = 1;
              (void)clear_bit(28, pxout.testbits);
            }
//          printf("new dust %f %f %d\n", m14, m43, dust);
          }

    	}
      }

//    printf("Dust tests (setup) for water day: %d %d %f %f %f %f %f %f %f %f %f %f \n",
//     		irclr, pxin.sh_ocean, sst_thrsh, a, corr, midpt, sfcdif, pxin.sfctmp, m31,
//     		stdev, m03, rat4);

    }

/**********************************************************************/

    else if(pxin.land && pxin.night) {


//    Check clear sky tests.
      irclr = 1;
      bit_test = check_bits(15, pxout.testbits);
      qa_test = check_bits(15, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
      bit_test = check_bits(17, pxout.testbits);
      qa_test = check_bits(17, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
      bit_test = check_bits(18, pxout.testbits);
      qa_test = check_bits(18, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;

//    Do not perform dust tests if tests reported in bits
//    15, 17, or 18 reported cloud.

      if(irclr) {

//      Perform dust test.
    	if ( rintf(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data) &&
    		 rintf(m32) != rintf(bad_data) && rintf(m29) != rintf(bad_data) ) {

//        Set bit to indicate that dust test was attempted.
    	  (void)set_bit(28, pxout.qabits);

    	  tdiff1 = m31 - m32;
    	  tdiff2 = m29 - m31;
    	  tdiff3 = m20 - m31;

    	  if(m31 > ncd_lnm31[0]) {

    	    band = 31;
    	    stdev = get_regstd(band);

    	    if(stdev < ncd_lnsig31[0]) {

    	      if(tdiff1 < ncd_lndf1[0] && tdiff2 > ncd_lndf2[0]) {
    	        if(tdiff3 <= ncd_lndf3[0]) {
    	          (void)clear_bit(28, pxout.testbits);
    	          dust = 1;
    	        }
    	        else if(tdiff3 > ncd_lndf3[1]) {
    	          (void)clear_bit(28, pxout.testbits);
    	          dust = 1;
    	        }
    	      }
    	    }
    	  }

    	}
      }

//    printf("Dust tests for land night: %d %d %f %f %f %f %f %f %f %f \n",
//     		irclr, dust, m31, m32, m29, m20, tdiff1, tdiff2, tdiff3, stdev);

    }

/**********************************************************************/

    else if(pxin.water && pxin.night) {

//    Night water dust tests.

//    Check clear sky tests.
      irclr = 1;
      bit_test = check_bits(15, pxout.testbits);
      qa_test = check_bits(15, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;
      bit_test = check_bits(18, pxout.testbits);
      qa_test = check_bits(18, pxout.qabits);
      if( (qa_test && !bit_test) || !qa_test ) irclr = 0;

//    Do not perform dust tests if tests reported in bits
//    15 or 18 reported cloud.

      if(irclr) {

    	if (rintf(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data) &&
    	    rintf(m29) != rintf(bad_data) && rintf(m20) != rintf(bad_data) &&
    	    pxin.sfctmp > 0.0 && pxin.sfctmp < 350.0) {

//        Set bit to indicate that dust test was attempted.
    	  (void) set_bit(28, pxout.qabits);

//        Perform dust tests.

          if(pxin.sh_ocean) {
            sst_thrsh = ncd_wnsst[0];
          }
          else {
            sst_thrsh = ncd_wnsst[1];
          }

          tdiff1 = m31 - m32;
          if(tdiff1 >= ncd_wndf1[0]) {
            midpt = sst_thrsh + (2.0 * rintf(tdiff1) );
          }
          else {
            midpt = sst_thrsh;
          }

          a = pxin.vza / mx_vza;
          corr = ( powf(a,4) ) * 3.0;
          midpt = midpt + corr;

          sfcdif = pxin.sfctmp - m31;

          if( sfcdif < midpt ) {

            band = 31;
            stdev = get_regstd(band);

            if(rintf(rg_var.mean) != rintf(bad_data) && stdev <= ncd_wnsig31[0]) {

//            Set bit to indicate that dust test was attempted.
              (void) set_bit(28, pxout.qabits);

              tdiff2 = m29 - m31;
              tdiff3 = m20 - m31;

              if(tdiff1 < ncd_wndf1[1] && tdiff3 < ncd_wndf3[0]) {
                (void)clear_bit(28, pxout.testbits);
                dust = 1;
              }
              else if(tdiff1 < ncd_wndf1[2] && tdiff3 < ncd_wndf3[1] && tdiff3 > ncd_wndf3[2]) {
                (void)clear_bit(28, pxout.testbits);
                dust = 1;
              }
            }
          }

    	}
      }

//    printf("Dust tests for water night: %d %d %d %f %f %f %f %f %f %f %f %f %f %f \n",
//     		irclr, dust, pxin.sh_ocean, sst_thrsh, a, corr, midpt, sfcdif, pxin.sfctmp, m31,
//     		stdev, tdiff1, tdiff2, tdiff3);

    }

    if (thk_smoke || fire || dust || thk_smk_water) pxin.nonCloudObstruction = 1;

/**********************************************************************/

    return;

/**********************************************************************/

}
