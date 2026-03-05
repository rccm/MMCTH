/*$Id: clavrx_volash_mask.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

/****************************************************************************/
/* Ash Detection Routine (01/06/2005)                                       */
/* Author:                                                                  */
/* Michael J. Pavolonis, CIMSS/SSEC/UW (mpav@ssec.wisc.edu)                 */
/****************************************************************************/	   

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

#include "clavrx_volash_mask.h"

void clavrx_volash_mask (unsigned char verbose,
                         imagerL1* imgr1, imagerL2* imgr2, 
			 sounderL1* sndr1, sounderL2* sndr2,
			 nwp_params nwp, rtm_profiles** rtm, 
			 rtm_toa rclr, rad_utils rutil)

{
    
  extern double 
    A_rat3865_water[],
    B_rat3865_water[],
    C_rat3865_water[],
    D_rat3865_water[],
    E_rat3865_water[];   
  
  unsigned char
    gflag,
    *hot_flg,
    *certain_ash_flg,

    YES_ASH;
    
  int
    DAY=0,
    n,i,
    index1,
    index2,
    index3,
    x,y,
    icol,irow,
    start_col,end_col,
    start_row,end_row,

    current_line,current_pix,
    ash_check,
    warm_check,
    splitwin_check,
    avg_count;
    
  float
    rat3865,
    RAT3865_THRES,
    MAX_BTD1112,
    MAX_BTD1112_RAT3865_QUALITY,
    RAT3865_DELTA,
    foo,
    *ash_frac,
    *btd1112,
    avg_bt31,
    avg_btd1112;
    

    
 
    
  char *rout = {"clavrx_volash_mask.c"};
    
  if ((hot_flg = (unsigned char *) malloc(imgr1->npts*sizeof(unsigned char))) == NULL)
    error_allo(rout,"hot_flg");
  if ((certain_ash_flg = (unsigned char *) malloc(imgr1->npts*sizeof(unsigned char))) == NULL)
    error_allo(rout,"certain_ash_flg");
  /*if ((ash_quality = (unsigned char *) malloc(imgr1->npts*sizeof(unsigned char))) == NULL)
    error_allo(rout,"ash_quality");*/
  if ((ash_frac = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"ash_frac");  
  if ((btd1112 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"btd1112");
    
  for (n=0; n < imgr1->npts; n++) {
    hot_flg[n] = NO;
    ash_frac[n] = 0.0;
    btd1112[n] = imgr1->bt31[n] - imgr1->bt32[n];
    rat3865 = imgr1->ref20[n]/imgr1->ref1[n];
    YES_ASH = NO;
    if (fabs(imgr1->lat[n]) < 30.0) {
      if ((imgr1->bt31[n] < 280.0 && rat3865 > 1.0 && btd1112[n] < 0.0) ||
          (imgr1->bt31[n] < 285.0 && rat3865 > 1.0 && btd1112[n] < -1.0) ||
          (imgr1->bt31[n] < 277.0 && btd1112[n] < -2.0 && rat3865 > 0.7) ||
	  (imgr1->bt31[n] < 233.0 && imgr1->ref20[n] > 0.20 && imgr1->ref1[n] < 0.60 && imgr1->sfc_type[n] != DESERT_IGBP)) YES_ASH = YES;
    } else if (fabs(imgr1->lat[n]) >= 30.0 && fabs(imgr1->lat[n]) < 60.0) {
      if ((imgr1->bt31[n] < 270.0 && rat3865 > 1.0 && btd1112[n] < -0.5 && imgr1->sfc_type[n] != DESERT_IGBP) || 
          (imgr1->bt31[n] < 270.0 && rat3865 > 0.7 && btd1112[n] < -1.0 && imgr1->sfc_type[n] != DESERT_IGBP) ||
          (imgr1->bt31[n] < 277.0 && btd1112[n] < -2.0 && rat3865 > 0.7) ||
	  (imgr1->bt31[n] < 235.0 && imgr1->ref20[n] > 0.20 && imgr1->ref1[n] < 0.60 && imgr1->sfc_type[n] != DESERT_IGBP)) YES_ASH = YES;
    } else if (fabs(imgr1->lat[n]) >= 60.0) {
      if ((imgr1->bt31[n] < 270.0 && rat3865 > 1.1 && btd1112[n] < -0.5 && imgr1->solzen[n] < 70.0) ||
          (imgr1->bt31[n] < 240.0 && imgr1->ref20[n] > 0.20 && imgr1->solzen[n] < 70.0 && imgr1->ref1[n] < 0.80) ||
	  (imgr1->bt31[n] < 277.0 && btd1112[n] < -3.0) ||
	  (imgr1->bt31[n] < 245.0 && btd1112[n] < -0.5 && imgr1->ref20[n] > 0.10 && imgr1->solzen[n] < 70.0)) YES_ASH = YES;
    }
    if (YES_ASH == YES) {
      printf("Found good ash pixel - lat=%f  lon=%f\n",imgr1->lat[n],imgr1->lon[n]);
      current_line = n/imgr1->ncol;
      current_pix = n % imgr1->ncol;
      start_row = max(0,current_line-200);
      end_row = min(imgr1->nrow-1,current_line+200);
      start_col = max(0,current_pix-200);
      end_col = min(imgr1->ncol-1,current_pix+200);
      for (y=start_row; y<=end_row; y++) {
        for (x=start_col; x<=end_col; x++) {
          i = (y*imgr1->ncol) + x;
	  certain_ash_flg[i] = YES;
        }
      }
    }
  }
  
  /*--------------------------------------------------------------------------*/
  /* Loop over all pixels stored in the pointers.                             */
  /*--------------------------------------------------------------------------*/
  
  for (n=0; n < imgr1->npts; n++) {
    
    /*--------------------------------------------------------------------------*/
    /* Check solar zenith angle to determine if it is day or night.             */
    /*--------------------------------------------------------------------------*/
      
    DAY = YES;
    if (imgr1->solzen[n] >= 70.0)
    DAY = NO;  
    
    /*--------------------------------------------------------------------------*/
    /* Determine which viewing zenith, solar zenith, and scattering angle bin   */
    /* the current pixel resides in.                                            */
    /*--------------------------------------------------------------------------*/
      
    index1 = min(NVZA-1,max(0,(int)(imgr1->satzen[n]/10.0)));
    index2 = min(NSZA-1,max(0,(int)(imgr1->solzen[n]/10.0)));
    index3 = min(NSCT-1,max(0,(int)(imgr1->scatzen[n]/10.0)));
    
    if (fabs(imgr1->lat[n]) < 20.0) {
      MAX_BTD1112 = MAX_BTD1112_TROPICS;
      MAX_BTD1112_RAT3865_QUALITY = MAX_BTD1112_RAT3865_QUALITY_TROPICS;
      RAT3865_DELTA = RAT3865_DELTA_TROPICS;
      if (imgr1->sfc_type[n] == WATER_IGBP && imgr1->glintzen[n] >= 150.0 && certain_ash_flg[n] == NO) MAX_BTD1112 = 0.7;
    }
    else if (fabs(imgr1->lat[n]) >= 20.0 && fabs(imgr1->lat[n]) < 45.0) {
      MAX_BTD1112 = MAX_BTD1112_MIDLAT;
      MAX_BTD1112_RAT3865_QUALITY = MAX_BTD1112_RAT3865_QUALITY_MIDLAT;
      RAT3865_DELTA = RAT3865_DELTA_MIDLAT;
      /*if (DAY == YES && imgr1->ref1[n] > 0.10) {
        MAX_BTD1112 -= 0.5;
	MAX_BTD1112_RAT3865_QUALITY -= 0.5;
      }*/
      if (imgr1->sfc_type[n] == WATER_IGBP && imgr1->glintzen[n] >= 150.0 && certain_ash_flg[n] == NO) MAX_BTD1112 = 0.0;
    }
    else {
      MAX_BTD1112 = MAX_BTD1112_HIGHLAT;
      MAX_BTD1112_RAT3865_QUALITY = MAX_BTD1112_RAT3865_QUALITY_HIGHLAT;
      RAT3865_DELTA = RAT3865_DELTA_HIGHLAT;
    }
  
    /*ash_quality[n] = 0;*/
    imgr2->aeromask[n] = NO_AEROSOL;
    
    /*--------------------------------------------------------------------------*/
    /* If daytime, use the daytime ash detection algorithm.                     */
    /*--------------------------------------------------------------------------*/
      
    if (DAY == YES) {
    
      rat3865 = imgr1->ref20[n]/imgr1->ref1[n];
    
      /*--------------------------------------------------------------------------*/
      /* Get the threshold RAT(3.75,0.65) value for daytime ash detection.        */
      /*--------------------------------------------------------------------------*/ 
      
      if (imgr1->sfc_type[n] == WATER_IGBP) {
        RAT3865_THRES = poly4(A_rat3865_water[index3],B_rat3865_water[index3],
	                      C_rat3865_water[index3],D_rat3865_water[index3],
	         	      E_rat3865_water[index3],imgr1->ref1[n]) + 0.025;
      }
      else if (imgr1->sfc_type[n] == DESERT_IGBP) {
        /*RAT3865_THRES = 999.0;*/
	RAT3865_THRES = poly4(A_rat3865_water[index3],B_rat3865_water[index3],
	                      C_rat3865_water[index3],D_rat3865_water[index3],
	         	      E_rat3865_water[index3],imgr1->ref1[n]) + 0.025;
      }
      else {
	RAT3865_THRES = poly4(A_rat3865_water[index3],B_rat3865_water[index3],
	                      C_rat3865_water[index3],D_rat3865_water[index3],
	         	      E_rat3865_water[index3],imgr1->ref1[n]) + 0.025;
      }
      
      foo = 285.0;
      if (imgr1->satzen[n] > 45.0) foo = 283.0;
      if (imgr1->satzen[n] > 58.0) foo = 282.0;
      
      gflag = NO;
      if (imgr1->glintzen[n] >= 150.0) gflag = YES;
      if (gflag == YES && fabs(imgr1->lat[n]) < 30.0 && imgr1->bt31[n] < 280.0) gflag = NO;
      
      /****************************GROUP I TESTS****************************/
      /*RAT[3.75,0.65] dominated tests*/
      
      /*WATER SURFACE*/
      if (certain_ash_flg[n] == YES) {
        if (rat3865 > RAT3865_THRES-0.1 && btd1112[n] < MAX_BTD1112 &&
            imgr1->ref1[n] > 0.04 && imgr1->ref1[n] < 0.3 && imgr1->bt31[n] < 295.0 && 
	    imgr1->satzen[n] < 55.0 && imgr1->sfc_type[n] == WATER_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
	    
	if (imgr1->glintzen[n] > 160.0 && imgr1->bt31[n] > 293.0 && imgr1->sfc_type[n] == WATER_IGBP) imgr2->aeromask[n] = NO_AEROSOL;
	  
	if (rat3865 > 1.20 && fabs(imgr1->lat[n]) < 20.0 && imgr1->ref1[n] > 0.1 &&
	    imgr1->ref1[n] < 0.2 && imgr1->bt31[n] < 283.0 && imgr1->satzen[n] < 55.0 &&
	    imgr1->sfc_type[n] == WATER_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
      }
      else {
        if (rat3865 > (RAT3865_THRES + 0.1) &&
	    btd1112[n] < MAX_BTD1112 && imgr1->sfc_type[n] == WATER_IGBP &&
	    imgr1->bt31[n] < 290.0 && imgr1->ref1[n] < 0.20 &&
	    imgr1->ref1[n] > 0.06 && gflag == NO) imgr2->aeromask[n] = MOSTLY_ASH;
      }
      
      /*LAND SURFACE*/
      if (certain_ash_flg[n] == YES) {
        if (fabs(imgr1->lat[n]) > 20.0) {
	  if (rat3865 > RAT3865_THRES-0.025 && btd1112[n] < MAX_BTD1112+0.0 &&
              imgr1->ref1[n] > 0.04 && imgr1->ref1[n] < 0.4 && imgr1->bt31[n] < 295.0 && 
	      imgr1->satzen[n] < 55.0 && imgr1->sfc_type[n] != WATER_IGBP &&
	      imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
	}
	else {
	  if (rat3865 > RAT3865_THRES-0.025 && btd1112[n] < MAX_BTD1112-0.5 &&
              imgr1->ref1[n] > 0.04 && imgr1->ref1[n] < 0.4 && imgr1->bt31[n] < 295.0 && 
	      imgr1->satzen[n] < 55.0 && imgr1->sfc_type[n] != WATER_IGBP &&
	      imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
	}
	  
	if (rat3865 > 1.20 && fabs(imgr1->lat[n]) < 20.0 && imgr1->ref1[n] > 0.1 &&
	    imgr1->ref1[n] < 0.2 && imgr1->bt31[n] < 283.0 && imgr1->satzen[n] < 55.0 &&
	    imgr1->sfc_type[n] != WATER_IGBP &&
	    imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
      }
      else {
        if (rat3865 > RAT3865_THRES && btd1112[n] < MAX_BTD1112 &&
            imgr1->ref1[n] > 0.06 && imgr1->ref1[n] < 0.4 && imgr1->bt31[n] < 290.0 && 
	    imgr1->satzen[n] < 55.0 && imgr1->sfc_type[n] != WATER_IGBP &&
	    imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
      }
      
      /****************************GROUP II TESTS****************************/
      /*BTD[11,12] dominated tests*/
      
      if (btd1112[n] < -3.0 && imgr1->bt31[n] < 270.0 && 
          imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
      
      if (btd1112[n] < -2.0 && rat3865 > 0.95 && imgr1->ref1[n] < 0.20) imgr2->aeromask[n] = MOSTLY_ASH;
      
      if (btd1112[n] < -0.5 && rat3865 > 0.95 && imgr1->ref1[n] < 0.10) imgr2->aeromask[n] = MOSTLY_ASH;
      
      /*High Latitude Focus*/
      if (btd1112[n] < 0.0 && imgr1->bt31[n] < 277.0 && 
          imgr1->sfc_type[n] != DESERT_IGBP &&  rat3865 > 0.6) imgr2->aeromask[n] = MOSTLY_ASH;
      
      /*if (fabs(imgr1->lat[n]) < 20.0 && btd1112[n] < -1.0  && 
          imgr1->sfc_type[n] != DESERT_IGBP && imgr1->bt31[n] < foo) imgr2->aeromask[n] = MOSTLY_ASH;*/
	  
      /*Added vis ref threshold to avoid flagging stratocumulus decks*/
      /*Using 3.75/0.65 ratio as so to avoid finding dust*/
      /*if (fabs(imgr1->lat[n]) < 60.0 && btd1112[n] < 0.0 &&
          imgr1->bt31[n] > 280.0 && imgr1->ref1[n] < 0.3 && 
	  imgr1->sfc_type[n] == WATER_IGBP && rat3865 > 1.0) imgr2->aeromask[n] = MOSTLY_ASH;*/
	  
      if (fabs(imgr1->lat[n]) < 20.0 && btd1112[n] < -0.5 &&
          rat3865 > 0.55 && imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
	  
      /*if (imgr1->bt31[n] < 290.0 && imgr1->ref1[n] > 0.04 && imgr1->sfc_type[n] == WATER_IGBP &&
          imgr1->glintzen[n] < 150.0 && rat3865 > 1.0 && imgr1->ref1[n] < 0.06 &&
	  btd1112[n] < 0.5) imgr2->aeromask[n] = MOSTLY_ASH;*/
      
      
      if (certain_ash_flg[n] == YES) {
        if (btd1112[n] < 0.0 && rat3865 > 0.5 && imgr1->bt31[n] < 290.0 && imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
	/*if (fabs(imgr1->lat[n]) > 50.0 && btd1112[n] < -0.5 && imgr1->sfc_type[n] != DESERT_IGBP && imgr1->satzen[n] < 50.0 &&
	    imgr1->ref20[n] > 0.03 && imgr1->ref20[n] < 0.15 && imgr1->bt31[n] < 270.0) imgr2->aeromask[n] = MOSTLY_ASH;
	if (fabs(imgr1->lat[n]) > 50.0 && btd1112[n] < -0.2 && imgr1->sfc_type[n] != DESERT_IGBP && imgr1->satzen[n] < 50.0 &&
	    imgr1->ref20[n] > 0.05 && imgr1->ref20[n] < 0.15 && rat3865 < 0.1) imgr2->aeromask[n] = MOSTLY_ASH;*/
	if (fabs(imgr1->lat[n]) > 50.0 && btd1112[n] < -0.2 && imgr1->sfc_type[n] != DESERT_IGBP && rat3865 < 0.2 &&
	    imgr1->satzen[n] < 50.0 && imgr1->ref20[n] > 0.03) imgr2->aeromask[n] = MOSTLY_ASH;
        /*if (fabs(imgr1->lat[n]) > 50.0 && btd1112[n] < -0.5 && imgr1->bt31[n] < 270.0 && 
            imgr1->ref20[n] > 0.05 && imgr1->ref20[n] < 0.18 && imgr1->ref1[n] < 0.75 && rat3865 > 0.10 &&
	    imgr1->satzen[n] < 50.0 && imgr1->sfc_type[n] == LAND_SFC) imgr2->aeromask[n] = MOSTLY_ASH;*/
	if (btd1112[n] < 0.5 && rat3865 > 0.7 && imgr1->bt31[n] < 290.0 && imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = MOSTLY_ASH;
      }
      
      
      /****************************GROUP III TESTS****************************/
      /*REF[3.75] dominated tests*/
      
      /*if (imgr1->ref20[n] > 0.18 && imgr1->bt31[n] < 235.0 &&
          rat3865 > 0.3 && fabs(imgr1->lat[n]) < 50.0) imgr2->aeromask[n] = ASH_ICE;
      
      if (imgr1->ref20[n] > 0.20 && imgr1->bt31[n] < 235.0 &&
          rat3865 > 0.4 && fabs(imgr1->lat[n] >= 50.0)) imgr2->aeromask[n] = ASH_ICE;*/
	  
      if (imgr1->ref20[n] > 0.18 && imgr1->bt31[n] < 235.0) imgr2->aeromask[n] = ASH_ICE;
      
      if (imgr1->ref1[n] < 0.40 && imgr1->bt31[n] < 210.0 && imgr1->ref20[n] > 0.08) imgr2->aeromask[n] = MOSTLY_ASH;
      
      if (certain_ash_flg[n] == YES) {
        if (imgr1->ref20[n] > 0.1 && imgr1->bt31[n] < 243.0 && imgr1->ref1[n] < 0.70 &&
            rat3865 > 0.2 && fabs(imgr1->lat[n] < 90.0) && imgr1->satzen[n] < 58.0 && imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = ASH_ICE;
      
        /*if (imgr1->ref20[n] > 0.12 && imgr1->bt31[n] < 243.0 && imgr1->ref1[n] < 0.70 &&
            rat3865 > 0.4 && imgr1->satzen[n] < 55.0 && imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = ASH_ICE;*/
	    
	if (imgr1->ref1[n] < 0.40 && imgr1->bt31[n] < 210.0 && imgr1->ref20[n] > 0.06) imgr2->aeromask[n] = MOSTLY_ASH;
	if (imgr1->ref1[n] < 0.50 && imgr1->bt31[n] < 200.0 && imgr1->ref20[n] > 0.06) imgr2->aeromask[n] = MOSTLY_ASH;
      }
      
      /****************************GROUP IV TESTS****************************/
      /*Restoral tests*/
      
      if (certain_ash_flg[n] == NO) {
        if ((imgr1->ref1[n] > 0.12 && imgr1->bt31[n] > foo && rat3865 < 0.70) || 
            (imgr1->ref1[n] > 0.105 && imgr1->bt31[n] > foo+3.5 && rat3865 < 0.85) ||
	    (imgr1->ref1[n] > 0.10 && imgr1->bt31[n] > foo+5.0)) imgr2->aeromask[n] = NO_AEROSOL;    
      
        /*Ash should not be bright, yet very transmissive in IR*/
        if (imgr1->ref1[n] > 0.20 && imgr1->bt31[n] > 280.0 && imgr1->sfc_type[n] != WATER_IGBP &&
	    imgr1->sfc_type[n] != DESERT_IGBP) imgr2->aeromask[n] = NO_AEROSOL;
      }
      
      /*Look for saturated pixels-hot spots*/
      if (imgr1->bt20[n] > 350.0 && imgr1->sfc_type[n] != WATER_IGBP && imgr1->solzen[n] < 75.0) {
	x = n % imgr1->ncol;
	y = n/imgr1->ncol;
	start_col = ((x-25) >= 0) ? x-25 : 0;
	end_col = ((x+25) < imgr1->ncol) ? x+25 : imgr1->ncol-1;
	start_row = ((y-25) >= 0) ? y-25 : 0;
	end_row = ((y+25) < imgr1->nrow) ? y+25 : imgr1->nrow-1;
	/*printf("n=%d, y=%d, x=%d, start_row=%d, end_row=%d, start_col=%d, end_col=%d\n",
	       n,y,x,start_row,end_row,start_col,end_col);*/
	printf("HOT SPOT FOUND: %f  %f  %f  %f\n",imgr1->bt31[n],imgr1->bt20[n],imgr1->lat[n],imgr1->lon[n]);
	for (irow=start_row; irow<end_row; irow++) {
	  for (icol=start_col; icol<end_col; icol++) {
	    hot_flg[(irow*imgr1->ncol)+icol] = YES;
	  }
	}
      }
	    
    } /*End of Daytime ash detection*/ 
      
    /*--------------------------------------------------------------------------*/
    /* If nighttime, use the ash detection algorithm.                           */
    /*--------------------------------------------------------------------------*/
      
    else {
    
      /*BTD[11,12] test ==> TEST #1*/
      
      if (btd1112[n] < -3.0 && imgr1->bt31[n] < 270.0) {
        imgr2->aeromask[n] = MOSTLY_ASH;
      }
      
      if (fabs(imgr1->lat[n]) < 20.0 && btd1112[n] < -1.0) {
        imgr2->aeromask[n] = MOSTLY_ASH;
      }
      
      if (fabs(imgr1->lat[n]) < 20.0 && btd1112[n] < 0.0 &&
          imgr1->bt31[n] > 280.0) {
        imgr2->aeromask[n] = MOSTLY_ASH;
      }
      
      if (btd1112[n] < -1.0 && imgr1->ems20[n] > 1.2 && 
          imgr1->bt31[n] < 275.0 && imgr1->bt31[n] > 230.0 &&
	  fabs(imgr1->lat[n]) < 50.0) {
	imgr2->aeromask[n] = MOSTLY_ASH;
      }
      
      if (btd1112[n] < -0.5 && imgr1->ems20[n] > 1.2 &&
          imgr1->bt31[n] < 278.0 && imgr1->bt31[n] > 230.0 &&
	  fabs(imgr1->lat[n]) < 50.0) {
        imgr2->aeromask[n] = MOSTLY_ASH;
      }
      if (btd1112[n] < -2.0 && imgr1->ems20[n] > 1.0 &&
          imgr1->bt31[n] < 278.0 && imgr1->bt31[n] > 230.0 &&
	  fabs(imgr1->lat[n]) > 50.0) {
        imgr2->aeromask[n] = MOSTLY_ASH;
      }
      if (btd1112[n] < -1.0 && imgr1->ems20[n] > 1.1 &&
          imgr1->bt31[n] < 278.0 && imgr1->bt31[n] > 230.0 &&
	  fabs(imgr1->lat[n]) > 50.0) {
        imgr2->aeromask[n] = MOSTLY_ASH;
      }
      if (btd1112[n] < -0.2 && imgr1->ems20[n] > 1.2 &&
          imgr1->bt31[n] < 278.0 && imgr1->bt31[n] > 230.0 &&
	  fabs(imgr1->lat[n]) > 50.0) {
        imgr2->aeromask[n] = MOSTLY_ASH;
      }
	 
    } /*End of nighttime ash detection*/  
    
  }	/*End of for loop*/
  
  if (DAY == YES) {
  for (n=0; n<imgr1->npts; n++) {
    /*x = n % imgr1->ncol;
    y = n/imgr1->ncol;
    start_col = ((x-25) >= 0) ? x-25 : 0;
    end_col = ((x+25) < imgr1->ncol) ? x+25 : imgr1->ncol-1;
    start_row = ((y-25) >= 0) ? y-25 : 0;
    end_row = ((y+25) < imgr1->nrow) ? y+25 : imgr1->nrow-1;
    ash_count = 0;
    count = 0;
    for (irow=start_row; irow<end_row; irow++) {
      for (icol=start_col; icol<end_col; icol++) {
	if (imgr2->aeromask[n] >=1) ash_count++;
	count++;
      }
    }
    ash_frac[n] = (float)ash_count/(float)count;*/
    if (hot_flg[n] == YES && imgr1->cldmask[n] <= 1 && imgr2->aeromask[n] == NO_AEROSOL) {
      
      rat3865 = imgr1->ref20[n]/imgr1->ref1[n];
      index3 = min(NSCT-1,max(0,(int)(imgr1->scatzen[n]/10.0)));
    
      /*--------------------------------------------------------------------------*/
      /* Get the threshold RAT(3.75,0.65) value for daytime ash detection.        */
      /*--------------------------------------------------------------------------*/ 
      
      if (imgr1->sfc_type[n] == WATER_IGBP) {
        RAT3865_THRES = poly4(A_rat3865_water[index3],B_rat3865_water[index3],
	                      C_rat3865_water[index3],D_rat3865_water[index3],
	         	      E_rat3865_water[index3],imgr1->ref1[n]) + 0.025;
	RAT3865_THRES = RAT3865_THRES-0.2;
	foo = 2.0;
      }
      else if (imgr1->sfc_type[n] == DESERT_IGBP) {
        /*RAT3865_THRES = 999.0;*/
	RAT3865_THRES = poly4(A_rat3865_water[index3],B_rat3865_water[index3],
	                      C_rat3865_water[index3],D_rat3865_water[index3],
	         	      E_rat3865_water[index3],imgr1->ref1[n]) + 0.025;
	RAT3865_THRES = RAT3865_THRES-0.1;
	foo = 1.0;
      }
      else {
	RAT3865_THRES = poly4(A_rat3865_water[index3],B_rat3865_water[index3],
	                      C_rat3865_water[index3],D_rat3865_water[index3],
	         	      E_rat3865_water[index3],imgr1->ref1[n]) + 0.025;
	RAT3865_THRES = RAT3865_THRES-0.1;
	foo = 1.0;
      }
      
      
      if (btd1112[n] < foo && imgr1->bt31[n] < 295.0 && imgr1->ref1[n] > 0.05 && imgr1->ref1[n] < 0.40 &&
          rat3865 > RAT3865_THRES) {
        imgr2->aeromask[n] = MOSTLY_ASH;
        /*if ((imgr1->glintzen[n] < 150.0 && imgr1->sfc_type[n] == WATER_IGBP) || imgr1->sfc_type[n] != WATER_IGBP &&
	     imgr1->sfc_type[n] != DESERT_IGBP) ash_quality[n] = 1;*/
      }
      if ((imgr1->ref1[n] > 0.15 && btd1112[n] > 1.0) || (imgr1->ref1[n] > 0.12 && imgr1->bt31[n] > 290.0)) {
        imgr2->aeromask[n] = NO_AEROSOL;
	/*ash_quality[n] = 0;*/
      }
      if ((imgr1->ref1[n] > 0.10 && imgr1->bt31[n] > 288.0 && imgr1->sfc_type[n] != WATER_IGBP &&
	   imgr1->sfc_type[n] != DESERT_IGBP) || 
          (imgr1->ref1[n] > 0.08 && imgr1->bt31[n] > 292.0 && imgr1->sfc_type[n] != WATER_IGBP &&
	   imgr1->sfc_type[n] != DESERT_IGBP)) {
        imgr2->aeromask[n] = NO_AEROSOL;
	/*ash_quality[n] = 0;*/
      }
      /*should only be used over snow in practice*/
      if (fabs(imgr1->lat[n]) > 50.0 && btd1112[n] < 1.0 && imgr1->bt31[n] < 245.0 && 
          imgr1->ref20[n] > 0.05 && imgr1->ref20[n] < 0.20 && imgr1->ref1[n] > 0.30) {
        imgr2->aeromask[n] = MOSTLY_ASH;
	/*ash_quality[n] = 1;*/
      }
    }
  }
  }
  
  for (n=0; n<imgr1->npts; n++) {
    current_line = n/imgr1->ncol;
    current_pix = n % imgr1->ncol;
    if (imgr2->aeromask[n] == MOSTLY_ASH || imgr2->aeromask[n] == ASH_ICE) {
      ash_check = 0;
      warm_check = 0;
      splitwin_check = 0;
      avg_bt31 = 0.0;
      avg_btd1112 = 0.0;
      avg_count = 0;
      start_row = max(0,current_line-10);
      end_row = min(imgr1->nrow-1,current_line+10);
      start_col = max(0,current_pix-10);
      end_col = min(imgr1->ncol-1,current_pix+10);
      for (y=start_row; y<=end_row; y++) {
        for (x=start_col; x<=end_col; x++) {
	  i = (y*imgr1->ncol) + x;
	  if (imgr2->aeromask[i] == MOSTLY_ASH || imgr2->aeromask[i] == ASH_ICE) ash_check++;
	  if (imgr1->bt31[i] > 293.0 && imgr2->aeromask[n] > 0) warm_check++;
	  if (btd1112[i] > 1.9 && imgr2->aeromask[n] > 0) splitwin_check++;
	  if (imgr1->bt31[i] > 180.0 && imgr1->bt31[i] < 330.0 &&
	      imgr1->bt32[i] > 180.0 && imgr1->bt32[i] < 330.0 && imgr1->cldtype[i] != CIRRUS_TYPE) {
	    avg_bt31 += imgr1->bt31[i];
	    avg_btd1112 += btd1112[i];
	    avg_count++;
	  }
	}
      }
      avg_bt31 /= (float)avg_count;
      avg_btd1112 /= (float)avg_count;
      /*if ((ash_check < 50 && avg_bt31 > 290.0 && avg_btd1112 > 1.5) || ash_check < 20 ||
          (ash_check <= 100 && avg_btd1112 > 1.7 && avg_bt31 > 290.0) ||
	  ((float)warm_check/(float)ash_check > 0.98 && (float)splitwin_check/(float)ash_check > 0.98)) mod.ash_mask[n] = 0;*/
      if (((float)warm_check/(float)ash_check > 0.99 && (float)splitwin_check/(float)ash_check > 0.99 &&
            ash_check < 95) ||
	  ((float)warm_check/(float)ash_check > 0.50 && (float)splitwin_check/(float)ash_check > 0.50 &&
            ash_check < 50) ||
	  (ash_check < 20)) imgr2->aeromask[n] = NO_AEROSOL;
    }
  }
  
  free(hot_flg);
  free(ash_frac);
  free(btd1112);
  free(certain_ash_flg);
  /*free(ash_quality);*/
  
}	/*End of Routine*/
