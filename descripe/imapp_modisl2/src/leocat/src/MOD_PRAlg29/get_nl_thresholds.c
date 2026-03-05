/*
!C ********************************************************************
!Description:

  Function get_nl_thresholds.c

  Routine for setting 3.9 minus 11 um cloud test thresholds for use
  in night land conditions.
 
!Input parameters:
  btdiff        11-12 um BTD for current pixel

!Output Parameters:
  locut         Zero confidence value of clear sky confidence interval
  hicut         1.0 confidence value of clear sky confidence interval
  midpt         0.5 confidence value of clear sky confidence interval
  power         Power of clear sky confidence curve

!Revision History:
  Converted to C.
  R. Frey           10/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include "thresholds.h"

/**********************************************************************/


void get_nl_thresholds(float btdiff, float *locut, float *hicut,
		               float *midpt, float *power)

/**********************************************************************/

{

/*  Declarations */

    float lo_val;
    float hi_val;
    float lo_val_thr;
    float hi_val_thr;
    float conf_range;
    float a;
    
/**********************************************************************/

/*  Initializations */

    lo_val = 0.0;
    hi_val = 0.0;
    lo_val_thr = 0.0;
    hi_val_thr = 0.0;
    conf_range = 0.0;
    a = 0.0;

/**********************************************************************/

/*  Compute clear sky confidence interval boundaries and mid-point. */

    if(btdiff > bt_diff_bounds[0]) {

      *(hicut) = nl_11_4l[2];
      *(locut) = nl_11_4l[0];
      *(midpt) = nl_11_4l[1];
      *(power) = nl_11_4l[3];

    }

    else if(btdiff < bt_diff_bounds[1]) {

      *(hicut) = nl_11_4h[2];
      *(locut) = nl_11_4h[0];
      *(midpt) = nl_11_4h[1];
      *(power) = nl_11_4h[3];

    }

    else {

      lo_val = bt_diff_bounds[0];
      hi_val = bt_diff_bounds[1];
      lo_val_thr = nl_11_4m[0];
      hi_val_thr = nl_11_4m[1];
      *(power) = nl_11_4m[3];
      conf_range = nl_11_4m[2];

      a = (btdiff - lo_val) / (hi_val - lo_val);
      *(midpt) = lo_val_thr + (a*(hi_val_thr - lo_val_thr));
      *(hicut) = *(midpt) - conf_range;
      *(locut) = *(midpt) + conf_range;

    }
/*
    printf("\nget_nl_threshlds: \n");
    printf("%f %f %f \n", btdiff, bt_diff_bounds[0], bt_diff_bounds[1]);
    printf("%f %f %f %f %f \n", lo_val, hi_val, lo_val_thr, hi_val_thr, conf_range);
    printf("%f %f %f %f %f \n", a, *(midpt), *(locut), *(hicut), *(power));
    printf("\n");
*/
    return;

}
