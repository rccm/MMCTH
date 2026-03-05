/*$Id: modis_ir_cloud_typing.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <hdf.h>
#include "common_leocat.h"
#include "nwp_leocat.h"
#include "rtm_leocat.h"
#include "radutils_leocat.h"
#include "imagerL1_leocat.h"
#include "imagerL2_leocat.h"
#include "imagerL3_leocat.h"
#include "sounderL1_leocat.h"
#include "sounderL2_leocat.h"
#include "sounderL3_leocat.h"


void modis_ir_cloud_typing (unsigned char verbose,
                            imagerL1* imgr1, imagerL2* imgr2, 
			    sounderL1* sndr1, sounderL2* sndr2,
			    nwp_params nwp, rtm_profiles** rtm, 
			    rtm_toa rclr, rad_utils rutil)

{
  long int
    n;
    
  float
    btd8511;
    
  printf("loop through all pixels\n");
  (void)fflush(stdout);

  for (n=0; n < imgr1->npts; n++) {
  
    btd8511 = imgr1->bt29[n] - imgr1->bt31[n];

//  printf("print1: %ld %f %d %f %f\n", n,btd8511,imgr1->cldmask[n],imgr2->cld_emiss11[n],
//          imgr2->cld_emiss12[n]);
//  (void)fflush(stdout);
  
    if (imgr1->cldmask[n] == PROB_CLEAR || imgr1->cldmask[n] == CLEAR) {
      
//    printf("print5: %ld %d \n", n,imgr1->cldmask[n]);
//    (void)fflush(stdout);
  
      imgr2->cldtype[n] = CLEAR_TYPE;
      imgr2->cldphase[n] = CLEAR_PHASE;
      
//    printf("print2: %ld %f %d \n", n,btd8511,imgr2->cldphase[n]);
//    (void)fflush(stdout);
  
    }
    else {
     
      if (imgr1->bt31[n] <= 238.0 || btd8511 >= 0.5) {
        imgr2->cldtype[n] = TICE_TYPE;
	imgr2->cldphase[n] = ICE_PHASE;
      }
      else if (imgr1->bt31[n] > 238.0 && imgr1->bt31[n] < 268.0 && 
               btd8511 >= -0.25 && btd8511 < 0.5) {
        imgr2->cldtype[n] = MIXED_TYPE;
	imgr2->cldphase[n] = MIXED_PHASE;
      }
      else if (imgr1->bt31[n] > 238.0 && imgr1->bt31[n] < 268.0 && 
               btd8511 < -0.25 && btd8511 >= -1.0) {	   
        imgr2->cldtype[n] = UNKNOWN_TYPE;
	imgr2->cldphase[n] = UNKNOWN_PHASE;
      }
      else if (imgr1->bt31[n] > 238.0 && btd8511 <= -1.0) {
        imgr2->cldtype[n] = WATER_TYPE;
	imgr2->cldphase[n] = WATER_PHASE;
      }
      else if (imgr1->bt31[n] >= 285.0 && btd8511 <= -0.5) {
        imgr2->cldtype[n] = WATER_TYPE;
	imgr2->cldphase[n] = WATER_PHASE;
      }
      else {
        imgr2->cldtype[n] = UNKNOWN_TYPE;
	imgr2->cldphase[n] = UNKNOWN_PHASE;
      }
      
//    printf("print3: %ld %f %d \n", n,btd8511,imgr2->cldphase[n]);
//    (void)fflush(stdout);
  
    }

//  printf("print4: %ld %f %d \n", n,btd8511,imgr2->cldphase[n]);
//  (void)fflush(stdout);

  }
      
}
