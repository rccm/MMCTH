/*
!C *********************************************************************
!Description:

  Integer function perform_cloud_tests.c
  Generates clear sky confidence and associated output data for current
  pixel.
  Called from create_cloud_mask.c

!Input arguments:
  none

  Inputs stored in structure variable 'pxin' of type "pixel_in" defined
  in pixel.h..

!Output arguments:
  none

  Output cloud mask and associated products are stored in structure
  variable 'pxout' defined in pixel.h.

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include "pixel.h"


int perform_cloud_tests()


{

/*  Declarations */

    extern int water_day();
    extern int land_day();
    extern int water_night();
    extern int land_night();
    extern int polar_day();
    extern int polar_night();
    extern float NightSnow_Inv_check();
    extern float chk_spatial_var();
    extern float chk_sunglint();
    extern float chk_shallow_water();
    extern float chk_land();
    extern float chk_land_night();
    extern float chk_coast();

    int return_code;
    int ret_code;
    int k;

    float conf;

/*  Initializations */

    k = 0;
    while(k < 6) {
      pxout.testbits[k] = 0;
      k++;
    }
    k = 0;
    while(k < 10) {
      pxout.qabits[k] = 0;
      k++;
    }

    conf = -99.0;

/**********************************************************************

    Perform cloud tests.

**********************************************************************/

    if(pxin.polar) {
      if(pxin.day) ret_code = polar_day();
      if(pxin.night) ret_code = polar_night();
    }
    else if (pxin.water) {
      if(pxin.day) ret_code = water_day();
      if(pxin.night) ret_code = water_night();
    }
    else if (pxin.land) {
      if(pxin.day) ret_code = land_day();
      if(pxin.night) ret_code = land_night();
     }
/*
    printf("\nvalue of cm byte 0-2 %d\n", pxout.testbits[0]);
    printf("value of cm byte 1-2 %d\n", pxout.testbits[1]);
    printf("value of cm byte 2-2 %d\n", pxout.testbits[2]);
    printf("value of cm byte 3-2 %d\n", pxout.testbits[3]);
    printf("value of cm byte 4-2 %d\n", pxout.testbits[4]);
    printf("value of cm byte 5-2 %d\n", pxout.testbits[5]);
*/
//  printf("Initial clear sky confidence: %f\n", pxout.init_conf);
//  (void)fflush(stdout);

/**********************************************************************

    Confidence confirmation (clear-sky restoral) tests.

**********************************************************************/

    pxout.intermediate_conf = pxout.init_conf;

//  Night, snow/ice temperature inversion tests.
    if(pxin.night && (pxin.snow || pxin.ice) ) {
      conf = NightSnow_Inv_check();
      pxout.intermediate_conf = conf;
    }

//  Water.
//  printf("calling chk_spatial_var: %d %d %f \n", pxin.water, pxin.uniform, pxout.intermediate_conf);
    if(pxin.water && pxin.uniform && (!pxin.ice) && pxout.intermediate_conf <= 0.95 &&
    		       pxout.intermediate_conf >= 0.05) {
      conf = chk_spatial_var();
      pxout.intermediate_conf = conf;
//    printf("after chk_spatial_var: %f \n", pxout.intermediate_conf);
    }

//  Water day.
    if(pxin.day && pxin.water && pxin.uniform && pxin.sunglint &&
                pxout.intermediate_conf <= 0.95) {
      conf = chk_sunglint();
      pxout.intermediate_conf = conf;
    }

//  Daytime shallow water test.
    if(pxin.day && pxin.water && pxin.sh_ocean && (!pxin.ice) ) {
//    if(pxin.day && pxin.water && (!pxin.ice) ) {
      conf = chk_shallow_water();
      pxout.intermediate_conf = conf;
    }

//  Daytime land.
    if(pxin.day && pxin.land && !pxin.snow && !pxin.ice &&
    		pxout.intermediate_conf <= 0.95) {
      conf = chk_land();
      pxout.intermediate_conf = conf;
    }

//  Daytime coast (land) test.
    if(pxin.day && pxin.land && !pxin.snow && !pxin.ice &&
    	    pxin.coast) {
      conf = chk_coast();
      pxout.intermediate_conf = conf;
    }

//  Nighttime land test.
    if(pxin.night && pxin.land && !pxin.snow && !pxin.ice &&
    		pxout.intermediate_conf <= 0.95) {
      conf = chk_land_night();
      pxout.intermediate_conf = conf;
    }

/*
    printf("\nvalue of cm byte 0-4 %d\n", pxout.testbits[0]);
    printf("value of cm byte 1-4 %d\n", pxout.testbits[1]);
    printf("value of cm byte 2-4 %d\n", pxout.testbits[2]);
    printf("value of cm byte 3-4 %d\n", pxout.testbits[3]);
    printf("value of cm byte 4-4 %d\n", pxout.testbits[4]);
    printf("value of cm byte 5-4 %d\n", pxout.testbits[5]);
*/
/*********************************************************************/

    pxout.confidence = pxout.intermediate_conf;
//  printf("Final confidence: %f \n", pxout.confidence);
//  (void)fflush(stdout);

/*********************************************************************/

    return_code = ret_code;
    return return_code;
}
