/*$Id: ph_viirs_cloud_typing.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

/****************************************************************************/
/* VIIRS cloud typing routine (03/02/2004)                                  */
/* Authors:                                                                 */
/* Michael J. Pavolonis, CIMSS/SSEC/UW (mpav@ssec.wisc.edu)                 */
/* Andrew K. Heidinger, NOAA/NESDIS/ORA                                     */
/*                                                                          */
/* This routine contains a series of spectral tests aimed at identifying    */
/* cloud phase and cloud overlap during the day and at night.  Most of the  */
/* tests used here are mature (i.e. daytime cloud overlap) however, a few   */
/* tests may be subject to slight modification in the near future.  This    */
/* routine will loop through pointers of MODIS spectral data.  The MODIS    */
/* cloud mask is used to identify pixels which are cloudy, only cloudy      */
/* pixels are assigned a cloud type.  The following are valid cloud types:  */
/* 2: WATER                                                                 */
/* 3: SUPERCOOLED WATER/MIXED PHASE                                         */
/* 4: OPAQUE ICE                                                            */
/* 5: CIRRUS (non-opaque ice)                                               */
/* 6: MULTI-LAYERED CLOUD                                                   */
/* The values of 0 and 1 are assigned to clear and partly cloudy pixels     */
/* respectively.                                                            */
/*                                                                          */
/* This routine is part of a larger MODIS processing system, if simulated   */
/* VIIRS data is to be used, some modification may be needed (i.e. cloud    */
/* mask codes may be different, variables names will likely need to be      */
/* changed, etc...).                                                        */
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

#include "ph_cld_type.h"

void ph_viirs_cloud_typing (unsigned char verbose,
                            imagerL1* imgr1, imagerL2* imgr2,
			    sounderL1* sndr1, sounderL2* sndr2,
			    nwp_params nwp, rtm_profiles** rtm,
			    rtm_toa rclr, rad_utils rutil)

{

  extern double
    A_win_over[][NVZA],
    B_win_over[][NVZA],
    C_win_over[][NVZA],
    D_win_over[][NVZA],
    E_win_over[][NVZA],
    MIN_win_over[][NVZA],
    A_nir_over_water[],
    B_nir_over_water[],
    C_nir_over_water[],
    D_nir_over_water[],
    E_nir_over_water[],
    A_nir_over_land[],
    B_nir_over_land[],
    C_nir_over_land[],
    D_nir_over_land[],
    E_nir_over_land[],
    A_cirrus[],
    B_cirrus[],
    C_cirrus[],
    D_cirrus[],
    E_cirrus[],
    A_3811[],
    B_3811[],
    C_3811[],
    D_3811[],
    E_3811[],
    A_8511[],
    B_8511[],
    C_8511[],
    D_8511[],
    E_8511[];

  char *rout = {"ph_viirs_cloud_typing.c"};

  int
    DAY,
    n,
    index1,
    index2,
    index3,
    wflg;

  int
    cld_check,
    water_check,
    start_line,
    end_line,
    start_pix,
    end_pix,
    x,
    y,
    current_line,
    current_pix;

  long
    i;

  float
    watfrac;


  float
    BAD = 9999.0;

  float
    WIN_OVER_THRES,
    CIRRUS_THRES,
    NIR_THRES6,
    NIR_THRES20,
    NIR_CIRRUS_THRES,
    MIN_REF26_OVER,
    BTD8511_THRES,
    NIR_OVERLAP_THRES,
    BTD3811_THRES,
    SNOW_REF6_THRES,
    MAX_SNOW_REF6_THRES,


    REF26_WIN_CHECK_THRES,
    NIR_CIRRUS_THRES20,
    BTD1112_N_OVER_L,
    BTD1112_N_OVER_H,
    EMS20_N_OVER_L,
    EMS20_N_OVER_H;

  unsigned char
    *cirrus_quality;

  float
    *nir_ref,
    *ratio1665,
    *btd3811,
    *btd8511,
    *btd1112;

  if ((cirrus_quality = (unsigned char *) malloc(imgr1->npts*sizeof(unsigned char))) == NULL)
    error_allo(rout,"cirrus_quality");

  if ((nir_ref = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"nir_ref");
  if ((ratio1665 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"ratio1665");
  if ((btd3811 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"btd3811");
  if ((btd8511 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"btd8511");
  if ((btd1112 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"btd1112");

  if (verbose == YES)
    fprintf(stdout,"%sApplying PH cloud typing spectral tests\n",EXE_PROMPT);

  /*--------------------------------------------------------------------------*/
  /* Loop over all pixels stored in the pointers.                             */
  /*--------------------------------------------------------------------------*/

  for (n=0; n < imgr1->npts; n++) {
//for (n=1273092; n < 1273112; n++) {

//  if (imgr1->badmask[n]) {
//    imgr2->cldtype[n] = UNKNOWN_TYPE;
//    imgr2->cldphase[n] = UNKNOWN_TYPE;
//    continue;
//  }

    if (imgr1->satid == TERRA)
      nir_ref[n] = imgr1->ref6[n];
    else if (imgr1->satid == AQUA)
      nir_ref[n] = imgr1->ref7[n];
    else
      nir_ref[n] = imgr1->ref6[n];

    ratio1665[n] = nir_ref[n]/imgr1->ref1[n];
    btd3811[n] = imgr1->bt20[n] - imgr1->bt31[n];
    btd8511[n] = imgr1->bt29[n] - imgr1->bt31[n];
    btd1112[n] = imgr1->bt31[n] - imgr1->bt32[n];

    /*--------------------------------------------------------------------------*/
    /* Classify pixels that are 95% or 99% clear according to MODIS cloud       */
    /* mask as clear.                                                           */
    /*--------------------------------------------------------------------------*/

    if (imgr1->cldmask[n] == PROB_CLEAR || imgr1->cldmask[n] == CLEAR) {

      imgr2->cldtype[n] = CLEAR_TYPE;
      imgr2->cldphase[n] = CLEAR_PHASE;

    }

    /*--------------------------------------------------------------------------*/
    /* Determine the cloud type for pixels that are cloudy according to the     */
    /* MODIS cloud mask.*/
    /*--------------------------------------------------------------------------*/

    else {

      /*--------------------------------------------------------------------------*/
      /* Check solar zenith angle to determine if it is day or night.             */
      /*--------------------------------------------------------------------------*/

      DAY = YES;
      if (imgr1->solzen[n] >= 88.0)
      DAY = NO;

      /*--------------------------------------------------------------------------*/
      /* Determine which viewing zenith, solar zenith, and scattering angle bin   */
      /* the current pixel resides in.                                            */
      /*--------------------------------------------------------------------------*/

      index1 = min(NVZA-1,max(0,(int)(imgr1->satzen[n]/10.0)));
      index2 = min(NSZA-1,max(0,(int)(imgr1->solzen[n]/10.0)));
      index3 = min(NSCT-1,max(0,(int)(imgr1->scatzen[n]/10.0)));

      /*--------------------------------------------------------------------------*/
      /* Get the threshold BT(11um) - BT(12um) value for cirrus detection.        */
      /*--------------------------------------------------------------------------*/

      CIRRUS_THRES = max(min(poly4(A_cirrus[index1],B_cirrus[index1],C_cirrus[index1],
                                   D_cirrus[index1],E_cirrus[index1],imgr1->bt31[n]),MAX_CIRRUS),
		         MIN_CIRRUS);

      /*--------------------------------------------------------------------------*/
      /* Get the threshold BT(3.8um) - BT(11um) value for nighttime typing.       */
      /*--------------------------------------------------------------------------*/

      BTD3811_THRES = max(min(poly4(A_3811[index1],B_3811[index1],C_3811[index1],
                                    D_3811[index1],E_3811[index1],btd1112[n]),MAX_3811),
			  MIN_3811);

      /*--------------------------------------------------------------------------*/
      /* Get the threshold BT(8.5um) - BT(11um) value for day or night typing.    */
      /*--------------------------------------------------------------------------*/

      if (imgr1->bt31[n] > 320.0) {
        BTD8511_THRES = BAD;
      }
      else {
        BTD8511_THRES = poly4(A_8511[index1],B_8511[index1],C_8511[index1],
	                      D_8511[index1],E_8511[index1],imgr1->bt31[n]);
      }

      /*--------------------------------------------------------------------------*/
      /* Get the threshold BT(11um) - BT(12um) value for daytime cloud overlap.   */
      /*--------------------------------------------------------------------------*/

      if (imgr1->ref1[n] >= MIN_REF1_OVER && imgr1->ref1[n] <= 0.60) {
        WIN_OVER_THRES = max(poly4(A_win_over[index2][index1],B_win_over[index2][index1],
	                           C_win_over[index2][index1],D_win_over[index2][index1],
				   E_win_over[index2][index1],imgr1->ref1[n]),
			     MIN_win_over[index2][index1]) - 0.25;

      }
      else if (imgr1->ref1[n] > 0.60 && imgr1->ref1[n] <= MAX_REF1_OVER) {
        WIN_OVER_THRES = MIN_win_over[index2][index1] - 0.25;
      }
      else {
        WIN_OVER_THRES = BAD;
      }

      if (imgr1->sfc_type[n] == DESERT_IGBP && imgr1->ref8[n] < MIN_REF04_OVER)
        WIN_OVER_THRES = BAD;

      /*--------------------------------------------------------------------------*/
      /* Get the threshold REF(1.65um) value for daytime cloud overlap.           */
      /*--------------------------------------------------------------------------*/

      if (imgr1->ref26[n] >= MAX_REF26_OVER || imgr1->sfc_type[n] == DESERT_IGBP) {
        NIR_OVERLAP_THRES = BAD;
	REF26_WIN_CHECK_THRES = REF26_WIN_CHECK_THRES_LAND;
      }
      else {
	if (imgr1->sfc_type[n] == WATER_IGBP) {
	  NIR_OVERLAP_THRES = poly4(A_nir_over_water[index3],B_nir_over_water[index3],
	                            C_nir_over_water[index3],D_nir_over_water[index3],
				    E_nir_over_water[index3],imgr1->ref26[n]) + 0.03;
	  if (imgr1->satid == AQUA) NIR_OVERLAP_THRES = NIR_OVERLAP_THRES - 0.09;
	  REF26_WIN_CHECK_THRES = REF26_WIN_CHECK_THRES_WATER;
	}
	else {
	  NIR_OVERLAP_THRES = max(poly4(A_nir_over_land[index3],B_nir_over_land[index3],
	                                C_nir_over_land[index3],D_nir_over_land[index3],
					E_nir_over_land[index3],imgr1->ref26[n]),0.25) + 0.05;
	  if (imgr1->satid == AQUA) NIR_OVERLAP_THRES = NIR_OVERLAP_THRES - 0.09;
	  REF26_WIN_CHECK_THRES = REF26_WIN_CHECK_THRES_LAND;
       	}
      }

      /*--------------------------------------------------------------------------*/
      /* Set some NIR cloud typing thresholds.                                    */
      /*--------------------------------------------------------------------------*/

      if (imgr1->sfc_type[n] == WATER_IGBP) {
        NIR_THRES6 = REF6_THRES_OCEAN;
	NIR_THRES20 = REF20_THRES_OCEAN;
	NIR_CIRRUS_THRES = NIR_CIRRUS_THRES_WATER;
	NIR_CIRRUS_THRES20 = NIR_CIRRUS_THRES_WATER20;
	if (imgr1->lat[n] >= 50.0 || imgr1->lat[n] <= -50.0) {
	  MIN_REF26_OVER = MIN_REF26_OVER_WATER_HIGH;
	} else {
	  MIN_REF26_OVER = MIN_REF26_OVER_WATER_LOW;
	}
	if (imgr1->lat[n] > -30.0 && imgr1->lat[n] < 30.0) {
	  BTD1112_N_OVER_L = BTD1112_N_OVER_L_TROPWATER;
	  BTD1112_N_OVER_H = BTD1112_N_OVER_H_TROPWATER;
	  EMS20_N_OVER_L = EMS20_N_OVER_L_TROPWATER;
	  EMS20_N_OVER_H = EMS20_N_OVER_H_TROPWATER;
        } else {
	  BTD1112_N_OVER_L = BTD1112_N_OVER_L_MIDWATER;
	  BTD1112_N_OVER_H = BTD1112_N_OVER_H_MIDWATER;
	  EMS20_N_OVER_L = EMS20_N_OVER_L_MIDWATER;
	  EMS20_N_OVER_H = EMS20_N_OVER_H_MIDWATER;
        }
      }
      else if (imgr1->sfc_type[n] == DESERT_IGBP) {
	BTD1112_N_OVER_L = BTD1112_N_OVER_L_LAND;
	BTD1112_N_OVER_H = BTD1112_N_OVER_H_LAND;
	EMS20_N_OVER_L = EMS20_N_OVER_L_LAND;
	EMS20_N_OVER_H = EMS20_N_OVER_H_LAND;
	if (imgr1->lat[n] >= 50.0 || imgr1->lat[n] <= -50.0) {
	  MIN_REF26_OVER = MIN_REF26_OVER_LAND_HIGH;
	} else {
	  MIN_REF26_OVER = MIN_REF26_OVER_LAND_LOW;
	}
	if (imgr1->lat[n] > 12.0 && imgr1->lat[n] < 32.0 &&
	    imgr1->lon[n] < 45.0 && imgr1->lon[n] > -20.0) {
	  NIR_CIRRUS_THRES = NIR_CIRRUS_THRES_DESERT;
	  NIR_CIRRUS_THRES20 = NIR_CIRRUS_THRES_DESERT20;
	  BTD1112_N_OVER_L = BAD;
	  BTD1112_N_OVER_H = BAD;
	  EMS20_N_OVER_L = BAD;
	  EMS20_N_OVER_H = BAD;
	}
	else if (imgr1->lat[n] >= 60.0 || imgr1->lat[n] <= -60.0) {
	  NIR_CIRRUS_THRES = NIR_CIRRUS_THRES_WATER;
	  NIR_CIRRUS_THRES20 = NIR_CIRRUS_THRES_WATER20;
	  NIR_THRES6 = REF6_THRES_OCEAN;
	  NIR_THRES20 = REF20_THRES_OCEAN;
	}
	else {
	  NIR_THRES6 = REF6_THRES_LAND;
	  NIR_THRES20 = REF20_THRES_LAND;
	  NIR_CIRRUS_THRES = NIR_CIRRUS_THRES_LAND;
	  NIR_CIRRUS_THRES20 = NIR_CIRRUS_THRES_LAND20;
	}
      }
      else {
        BTD1112_N_OVER_L = BTD1112_N_OVER_L_LAND;
	BTD1112_N_OVER_H = BTD1112_N_OVER_H_LAND;
	EMS20_N_OVER_L = EMS20_N_OVER_L_LAND;
	EMS20_N_OVER_H = EMS20_N_OVER_H_LAND;
        if (imgr1->lat[n] >= 40.0 || imgr1->lat[n] <= -40.0) {
	  MIN_REF26_OVER = MIN_REF26_OVER_LAND_HIGH;
	} else {
	  MIN_REF26_OVER = MIN_REF26_OVER_LAND_LOW;
	}
	if (imgr1->lat[n] >= 60.0 || imgr1->lat[n] <= -60.0) {
	  NIR_THRES6 = REF6_THRES_OCEAN;
	  NIR_THRES20 = REF20_THRES_OCEAN;
	  NIR_CIRRUS_THRES = NIR_CIRRUS_THRES_WATER;
	  NIR_CIRRUS_THRES20 = NIR_CIRRUS_THRES_WATER20;
	}
	else {
	  NIR_THRES6 = REF6_THRES_LAND;
	  NIR_THRES20 = REF20_THRES_LAND;
	  NIR_CIRRUS_THRES = NIR_CIRRUS_THRES_LAND;
	  NIR_CIRRUS_THRES20 = NIR_CIRRUS_THRES_LAND20;
	}
      }

      if (imgr1->lat[n] >= 60.0 || imgr1->lat[n] <= -60.0) {
	SNOW_REF6_THRES = SNOW_REF6_THRES_HIGH;
	MAX_SNOW_REF6_THRES = MAX_SNOW_REF6_THRES_HIGH;
      }
      else {
	SNOW_REF6_THRES = SNOW_REF6_THRES_LOW;
	MAX_SNOW_REF6_THRES = BAD;
      }

      /*--------------------------------------------------------------------------*/
      /* Classify each pixel according to 11 um brightness temperature.           */
      /*--------------------------------------------------------------------------*/

      wflg = 1;

      if (imgr1->bt31[n] <= CERTAIN_ICE_BT11) {
        imgr2->cldtype[n] = TICE_TYPE;
	imgr2->cldphase[n] = ICE_PHASE;
        wflg = 0;
      }
      else if (imgr1->bt31[n] > CERTAIN_ICE_BT11 && imgr1->bt31[n] <= 253.16) {
        imgr2->cldtype[n] = TICE_TYPE;
	imgr2->cldphase[n] = ICE_PHASE;
      }
      else if (imgr1->bt31[n] > 253.16 && imgr1->bt31[n] <= 273.16) {
        imgr2->cldtype[n] = MIXED_TYPE;
	imgr2->cldphase[n] = MIXED_PHASE;
      }
      else {
        imgr2->cldtype[n] = WATER_TYPE;
	imgr2->cldphase[n] = WATER_PHASE;
      }

      /*--------------------------------------------------------------------------*/
      /* If daytime, use the daytime typing algorithm.                            */
      /*--------------------------------------------------------------------------*/

      if (DAY == YES) {

        cirrus_quality[n] = 1;

	/*--------------------------------------------------------------------------*/
	/* Check for cloud overlap with the NIR reflectance test.                   */
	/*--------------------------------------------------------------------------*/

	if (imgr1->ref26[n] < REF26_WIN_CHECK_THRES) {

          if (imgr1->ref26[n] > MIN_REF26_OVER && nir_ref[n] > NIR_OVERLAP_THRES &&
	      ratio1665[n] < 1.0 && imgr1->bt31[n] <  MAX_BT11_NIR_OVER &&
	      btd1112[n] > WIN_OVER_THRES && imgr1->bt31[n] > 220.0) {
            imgr2->cldtype[n] = OVERLAP_TYPE;
	  }

        }
        else {

          if (imgr1->ref26[n] > MIN_REF26_OVER && nir_ref[n] > NIR_OVERLAP_THRES &&
	      ratio1665[n] < 1.0 && imgr1->bt31[n] <  MAX_BT11_NIR_OVER && imgr1->bt31[n] > 220.0) {
            imgr2->cldtype[n] = OVERLAP_TYPE;
	  }

        }

        /*--------------------------------------------------------------------------*/
	/* Check for cloud overlap using the split-window brightness temperature    */
	/* test.                                                                    */
	/*--------------------------------------------------------------------------*/

	if (btd1112[n] > WIN_OVER_THRES && imgr1->bt31[n] < MAX_BT11_WIN_OVER &&
            nir_ref[n] > SNOW_REF6_THRES && imgr1->bt31[n] > 220.0) {
          imgr2->cldtype[n] = OVERLAP_TYPE;
	}

        /*--------------------------------------------------------------------------*/
	/* Check for cirrus clouds using the split window/NIR test.                 */
	/*--------------------------------------------------------------------------*/

	if ((imgr2->cldtype[n] != OVERLAP_TYPE && btd1112[n] > CIRRUS_THRES &&
	     imgr1->ref26[n] > 0.025 && imgr1->ref20[n] < NIR_CIRRUS_THRES20 && imgr1->bt31[n] > 230.0) ||
	    (imgr2->cldtype[n] != OVERLAP_TYPE && imgr1->ref1[n] < 0.40 &&
	     (imgr1->ref26[n]/imgr1->ref1[n]) > 0.12) && imgr1->bt31[n] > 230.0) {
          imgr2->cldtype[n] = CIRRUS_TYPE;
	  imgr2->cldphase[n] = ICE_PHASE;
	}

        /*--------------------------------------------------------------------------*/
	/* Re-classify each pixel if certain spectral conditions are met.           */
	/*--------------------------------------------------------------------------*/

	/*December 9, 2003*/
        /*Added nir ratio to SWBTD overlap test and TICE clouds*/
	/*cannot have a 1.65/0.65 ratio greater than 1.0*/

        if (BTD8511_THRES != BAD) {
          if (btd8511[n] >= BTD8511_THRES && imgr1->bt31[n] < 263.16 && ratio1665[n] < 1.0) {
            if (imgr2->cldtype[n] == MIXED_TYPE) imgr2->cldtype[n] = TICE_TYPE;
	    imgr2->cldphase[n] = ICE_PHASE;
	  }
          if ((wflg == 1 && btd8511[n] < BTD8511_THRES+0.2)  ||
	      (ratio1665[n] > 1.0)) {
            if (imgr2->cldtype[n] == TICE_TYPE) imgr2->cldtype[n] = MIXED_TYPE;
	    imgr2->cldphase[n] = MIXED_PHASE;
	  }
          if (btd8511[n] >= 0.5 && imgr1->ref26[n] > 0.01 &&
	      imgr1->ref20[n] < NIR_CIRRUS_THRES20) {
            if (imgr2->cldtype[n] == WATER_TYPE) imgr2->cldtype[n] = CIRRUS_TYPE;
	    imgr2->cldphase[n] = ICE_PHASE;
	  }
          if (btd8511[n] >= BTD8511_THRES &&
	      imgr1->ref26[n] > 0.025 && imgr1->ref20[n] < NIR_CIRRUS_THRES20) {
            if (imgr2->cldtype[n] == MIXED_TYPE) imgr2->cldtype[n] = CIRRUS_TYPE;
	    imgr2->cldphase[n] = ICE_PHASE;
	  }
        }

	if (imgr2->cldphase[n] == ICE_PHASE && btd1112[n] > CIRRUS_THRES+2.0) imgr2->cldphase[n] = MIXED_PHASE;

      }	/*End of Daytime ctyping*/

      /*--------------------------------------------------------------------------*/
      /* If nighttime, use the nighttime typing algorithm.                        */
      /*--------------------------------------------------------------------------*/

      else {

        cirrus_quality[n] = 0;

	/*--------------------------------------------------------------------------*/
	/* Check for cloud overlap with the nighttime window test.                  */
	/*--------------------------------------------------------------------------*/

	if (btd1112[n] > BTD1112_N_OVER_L && btd1112[n] < BTD1112_N_OVER_H &&
            btd3811[n] > BTD3811_N_OVER_L && btd3811[n] < BTD3811_N_OVER_H &&
            imgr1->ems20[n] > EMS20_N_OVER_L && imgr1->ems20[n] < EMS20_N_OVER_H &&
            imgr1->bt31[n] < 290.0 && imgr1->bt31[n] > 220.0) {
          imgr2->cldtype[n] = OVERLAP_TYPE;
	}

        /*--------------------------------------------------------------------------*/
	/* Check for cirrus clouds using the split window/NIR test.                 */
	/*--------------------------------------------------------------------------*/

	if (imgr2->cldtype[n] != OVERLAP_TYPE && btd1112[n] > CIRRUS_THRES && imgr1->ems20[n] > 1.3 && imgr1->bt31[n] > 230.0) {
          imgr2->cldtype[n] = CIRRUS_TYPE;
	  imgr2->cldphase[n] = ICE_PHASE;
	  if ((btd8511[n] >= BTD8511_THRES && imgr1->bt31[n] < 263.16) ||
	      (imgr1->ems20[n] > 1.9 && imgr1->bt31[n] < 300.0))
	    cirrus_quality[n] = 1;
	}

	if (imgr2->cldtype[n] != OVERLAP_TYPE && imgr2->cldtype[n] != TICE_TYPE && imgr1->ems20[n] > 1.1 && imgr1->bt31[n] < 300.0 && imgr1->bt31[n] > 230.0) {
	  imgr2->cldtype[n] = CIRRUS_TYPE;
	  imgr2->cldphase[n] = ICE_PHASE;
	  if ((btd8511[n] >= BTD8511_THRES && imgr1->bt31[n] < 263.16) ||
	      (imgr1->ems20[n] > 1.9 && imgr1->bt31[n] < 300.0))
	    cirrus_quality[n] = 1;
	}

        /*--------------------------------------------------------------------------*/
	/* Re-classify each pixel if certain spectral conditions are met.           */
	/*--------------------------------------------------------------------------*/

	if (BTD8511_THRES != BAD) {
          if (btd8511[n] >= BTD8511_THRES && imgr1->bt31[n] < 263.16) {
            if (imgr2->cldtype[n] == MIXED_TYPE) imgr2->cldtype[n] = TICE_TYPE;
	    imgr2->cldphase[n] = ICE_PHASE;
	  }
          if (wflg == 1 && btd8511[n] < BTD8511_THRES+0.2) {
            if (imgr2->cldtype[n] == TICE_TYPE) imgr2->cldtype[n] = MIXED_TYPE;
	    imgr2->cldphase[n] = MIXED_PHASE;
	  }
          if (btd8511[n] >= 0.5) {
            if (imgr2->cldtype[n] == WATER_TYPE) imgr2->cldtype[n] = CIRRUS_TYPE;
	    imgr2->cldphase[n] = ICE_PHASE;
	  }
        }

	if (imgr2->cldphase[n] == ICE_PHASE && btd1112[n] > CIRRUS_THRES+2.0) imgr2->cldphase[n] = MIXED_PHASE;

      }	/*End of nighttime typing*/

    }	/*End of cloudy pixel tests*/

  }	/*End of for loop*/

  if (verbose == YES)
    fprintf(stdout,"%sApplying PH cirrus filter\n",EXE_PROMPT);

  /*--------------------------------------------------------------------------*/
  /* Apply spatial filtering to cirrus type.                                  */
  /*--------------------------------------------------------------------------*/

  for (n=0; n<imgr1->npts; n++) {
//for (n=1273092; n < 1273112; n++) {
    current_line = n/imgr1->ncol;
    current_pix = n % imgr1->ncol;
    if (imgr2->cldtype[n] == CIRRUS_TYPE && cirrus_quality[n] == 0) {
      water_check = 0;
      cld_check = 0;
      start_line = max(0,current_line-DY_CIRRUS);
      end_line = min(imgr1->nrow-1,current_line+DY_CIRRUS);
      start_pix = max(0,current_pix-DX_CIRRUS);
      end_pix = min(imgr1->ncol-1,current_pix+DX_CIRRUS);
      for (y=start_line; y<=end_line; y++) {
        for (x=start_pix; x<=end_pix; x++) {
	  i = (y*imgr1->ncol) + x;
          if (imgr2->cldtype[i] == TICE_TYPE || imgr2->cldtype[i] == OVERLAP_TYPE ||
	      (imgr2->cldtype[i] == CIRRUS_TYPE && cirrus_quality[i] == 1)) {
	    cld_check++;
	  }
	  else if (imgr2->cldtype[i] == MIXED_TYPE || imgr2->cldtype[i] == WATER_TYPE) {
	    water_check++;
	    cld_check++;
	  }
	  else if (imgr2->cldtype[i] == CIRRUS_TYPE && cirrus_quality[i] == 0) {
	    cld_check++;
	  }
	}
      }

      watfrac = (float)water_check/(float)cld_check;

      if (water_check > cld_check*0.35) {
        if (imgr1->bt31[n] <= 273.16) {
	  imgr2->cldtype[n] = MIXED_TYPE;
	  imgr2->cldphase[n] = MIXED_PHASE;
	}
	else {
	  imgr2->cldtype[n] = WATER_TYPE;
	  imgr2->cldphase[n] = WATER_PHASE;
	}
	cirrus_quality[n] = 2;
      }
    }

  }

  free(cirrus_quality);
  free(nir_ref);
  free(ratio1665);
  free(btd3811);
  free(btd8511);
  free(btd1112);

}	/*End of Routine*/
