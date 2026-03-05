/*
!C ********************************************************************
!Description:

  Function get_pn_thresholds.c

  Routine for setting 3.9 minus 11 um cloud test thresholds for use
  in polar night conditions.
 
!Input parameters:
  bt_11         11 um brightness temperature for current pixel
  bt_bnds       11 um brightness temperature bounds for current test
  th_low        Test thresholds for 'bt_11' lt 'bt_bnds(1)'.
  th_mid1       Test thresholds for 'bt_11' between 'bt_bnds(1)' and
                     'bt_bnds(2)'.
  th_mid2       Test thresholds for 'bt_11' between 'bt_bnds(2)' and
                     'bt_bnds(3)'.
  th_mid3       Test thresholds for 'bt_11' between 'bt_bnds(3)' and
                     'bt_bnds(4)'.
  th_hi         Test thresholds for 'bt_11' gt 'bt_bnds(4)'.
 
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


void get_pn_thresholds(float bt_11, float *bt_bnds, float *th_low,
                    float *th_mid1, float *th_mid2, float *th_mid3, 
                    float *th_hi, float *locut, float *hicut, float *midpt,
                    float *power)


/**********************************************************************/

{

/*  Declarations */

    float lo_tmp;
    float hi_tmp;
    float lo_tmp_thr;
    float hi_tmp_thr;
    float conf_range;
    float a;
    
/**********************************************************************/

/*  Initializations */

    lo_tmp = 0.0;
    hi_tmp = 0.0;
    lo_tmp_thr = 0.0;
    hi_tmp_thr = 0.0;
    conf_range = 0.0;
    a = 0.0;

/**********************************************************************/

/*  Compute clear sky confidence interval boundaries and mid-point. */

    if(bt_11 < bt_bnds[0]) {

      *(hicut) = th_low[2];
      *(locut) = th_low[0];
      *(midpt) = th_low[1];
      *(power) = th_low[3];

    }

    else if(bt_11 > bt_bnds[3]) {

      *(hicut) = th_hi[2];
      *(locut) = th_hi[0];
      *(midpt) = th_hi[1];
      *(power) = th_hi[3];

    }

    else {

      if(bt_bnds[1] == 0.0 && bt_bnds[2] == 0.0) {

        lo_tmp = bt_bnds[0];
        hi_tmp = bt_bnds[3];
        lo_tmp_thr = th_mid1[0];
        hi_tmp_thr = th_mid1[1];
        *(power) = th_mid1[3];
        conf_range = th_mid1[2];

      }

      else if(bt_11 >= bt_bnds[0] && bt_11 < bt_bnds[1]) {

        lo_tmp = bt_bnds[0];
        hi_tmp = bt_bnds[1];
        lo_tmp_thr = th_mid1[0];
        hi_tmp_thr = th_mid1[1];
        *(power) = th_mid1[3];
        conf_range = th_mid1[2];

      }

      else if(bt_11 >= bt_bnds[1] && bt_11 < bt_bnds[2]) {

        lo_tmp = bt_bnds[1];
        hi_tmp = bt_bnds[2];
        lo_tmp_thr = th_mid2[0];
        hi_tmp_thr = th_mid2[1];
        *(power) = th_mid2[3];
        conf_range = th_mid2[2];

      }

      else if(bt_11 >= bt_bnds[2] && bt_11 <= bt_bnds[3]) {

        lo_tmp = bt_bnds[2];
        hi_tmp = bt_bnds[3];
        lo_tmp_thr = th_mid3[0];
        hi_tmp_thr = th_mid3[1];
        *(power) = th_mid3[3];
        conf_range = th_mid3[2];

      }

      a = (bt_11 - lo_tmp) / (hi_tmp - lo_tmp);
      *(midpt) = lo_tmp_thr + (a*(hi_tmp_thr - lo_tmp_thr));
      *(hicut) = *(midpt) - conf_range;
      *(locut) = *(midpt) + conf_range;

    }
/*
    printf("\nget_pn_threshlds: \n");
    printf("%f %f %f %f \n", bt_bnds[0], bt_bnds[1], bt_bnds[2], bt_bnds[3]);
    printf("%f %f %f %f %f %f \n", bt_11, lo_tmp, hi_tmp, lo_tmp_thr, hi_tmp_thr, conf_range);
    printf("%f %f %f %f %f \n", a, *(midpt), *(locut), *(hicut), *(power));
    printf("\n");
*/
    return;

}
