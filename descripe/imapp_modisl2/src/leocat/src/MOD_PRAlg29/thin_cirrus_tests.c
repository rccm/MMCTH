/*
************************************************************************
!F77

!Description:

   Performs thin cirrus tests.  These results of these tests do not
   contribute to the final confidence of clear sky but are reported in
   bits 9 (1.38 test, day only) and 11 (11-12 um test).

!Input parameters:
   none

   Inputs stored in structure variable 'pxin' of type "pixel_in" defined
   in pixel.h.

!Output Parameters:
   none

   Output values stored in structure variable 'pxout' of type "pixel.out".

!Revision History:
    Converted to C         11/07      R. Frey

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


void thin_cirrus_tests()


{

/***********************************************************************/

//  Declarations

    extern void set_bit(int, unsigned char[]);
    extern float cithr(int, float, float);

    int code;

    float diff;
    float thr;
    float btd_thr;
    float schi;
    float ci1, ci2;
    float cosvza;
    float m26, m31, m32;
    float lo_thresh, hi_thresh;

/**********************************************************************/

//  Initializations

//  printf("Executing thin_cirrus_tests \n");

    m26 = pxin.rad[25];
    m31 = pxin.rad[30];
    m32 = pxin.rad[31];

    code = 1;

/**********************************************************************/

    if( (!pxin.hi_elev) && (!(pxin.Antarctic && pxin.land)) && pxin.day && pxin.visusd) {

//    1.38 Um test for daytime scenes.

      if (rintf(m26) != rintf(bad_data)) {

//      Obtain proper thresholds.

//    	printf("Getting thin cirrus test thresholds\n");
        if( pxin.snow || pxin.ice) {
//        Snow/ice.
          hi_thresh = dstci[0];
          lo_thresh = dstci[1];
        }
    	else if(pxin.water) {
//        Ocean.
    	  hi_thresh = dotci[0];
    	  lo_thresh = dotci[1];
    	}
    	else if(pxin.land) {
//        Land.
    	  hi_thresh = dltci[0];
          lo_thresh = dltci[1];
    	}

//      printf("Data for thin cirrus test %d %d %d %f %f %f %f\n", pxin.snow, pxin.land, pxin.water, hi_thresh,
//        		                 lo_thresh, tci_ref3_tpw[0], pxin.tpw);
        if ( (pxin.land && !pxin.snow && !pxin.ice) ||
             ((pxin.snow || pxin.ice) && pxin.tpw > tci_ref3_tpw[0]) ||
             (pxin.water && !pxin.snow && !pxin.ice) ){

//        printf("Performing thin cirrus test\n");
          (void) set_bit(9, pxout.qabits);
          if (m26 < lo_thresh || m26 >= hi_thresh) {
            (void) set_bit(9, pxout.testbits);
          }
          else {
            pxin.thinCirrusSolar = 1;
          }

        }

      }

    }

    if(!pxin.snow && !pxin.ice) {

//    11-12 um test.

      if (rintf(m31) != rintf(bad_data) && rintf(m32) != rintf(bad_data)) {

        diff = m31 - m32;
        cosvza = cos(pxin.vza * dtr);
        if (fabsf(cosvza) > 0.0) {
          schi = 1.0 / cosvza;
        }
        else {
          schi = 99.0;
        }

/*      Interpolate look-up table values of 11 - 12 micron BT
        difference thresholds (function of viewing zenith and
        11 micron BT).
*/
        thr = cithr(1, schi, m31);
        if(thr < 0.1 || fabsf(schi-99.0) < 0.0001) {
          code = 0;
        }
        else {
          btd_thr = thr;
        }

//      Use a threshold range indicating very thin cirrus.
        if (code) {
          (void) set_bit(11, pxout.qabits);
          ci1 = btd_thr;
          ci2 = btd_thr + (0.3 * btd_thr);

          if (diff <= ci1 || diff > ci2) {
            (void) set_bit(11, pxout.testbits);
          }
          else {
            pxin.thinCirrusIR = 1;
          }

        }
      }

    }


/**********************************************************************/

    return;

/**********************************************************************/

}
