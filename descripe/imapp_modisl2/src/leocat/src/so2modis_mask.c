/*$Id: so2modis_mask.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

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


void so2modis_mask (unsigned char verbose,
                    imagerL1* imgr1, imagerL2* imgr2, 
	            sounderL1* sndr1, sounderL2* sndr2,
	            nwp_params nwp, rtm_profiles** rtm, 
	            rtm_toa rclr, rad_utils rutil)
		    
{
  size_t i;
  int i_nwp, ivza;
  float bias29, bias31, bias27, bias28, diff1, diff2;
  
  printf("foo\n");
  for (i=0; i<imgr1->npts; i++) {
    i_nwp = imgr1->index_nwp[i];
    ivza = imgr1->index_vza[i];
    bias29 = rtm[i_nwp][ivza].bt29_clr - imgr1->bt29[i];
    bias31 = rtm[i_nwp][ivza].bt31_clr - imgr1->bt31[i];
    diff1 = bias29 - bias31;
    
    bias27 = rtm[i_nwp][ivza].bt27_clr - imgr1->bt27[i];
    bias28 = rtm[i_nwp][ivza].bt28_clr - imgr1->bt28[i];
    diff2 = bias27 - bias28;
    diff1 = bias27 - bias28;
    diff1 = imgr1->bt27[i] - imgr1->bt28[i];
    imgr2->so2mask[i] = 0;
    if (diff1 < -11.0 && diff2 < -10.0) imgr2->so2mask[i] = 1;
    imgr2->col_aer[i] = diff1;
  }

}
