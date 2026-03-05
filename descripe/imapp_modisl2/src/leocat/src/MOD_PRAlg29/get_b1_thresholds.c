/*
!C ********************************************************************
!Description:

  Integer function get_b1_thresholds.c

  Computes band 1 (0.66 micron) cloud test thresholds for daytime land
  conditions.
  Computes band 8 (0.412 micron) cloud test thresholds for daytime
  desert conditions.

!Input arguments:
  none

  Input values stored in variables found in pixel.h

!Output arguments:
  locut         Zero confidence value of clear sky confidence interval
  hicut         1.0 confidence value of clear sky confidence interval
  midpt         0.5 confidence value of clear sky confidence interval
  power         Power of clear sky confidence curve

!Revision History:
  Converted to C from FORTRAN
  R. Frey           10/2008


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>
#include "thresholds.h"
#include "pixel.h"

/**********************************************************************/


int get_b1_thresholds(float *locut, float *hicut, float *midpt,
                      float *power)


/**********************************************************************/

{

/*  Declarations */

    int interp=0;
    int a, i, j, k, indvi, icnf, irtn;

    static float delta_ndvi_bin = 0.1;
    float x, x1, x2, y1, y2;
    float coeff[10][3][4], thr[3], thr_adj[3];

//    Transfer all coefficients to one array.

      i = 0;

      for (j=0; j<3; j++) {
        for( k=0; k<4; k++) {

    	  if(pxin.ndvibk < des_ndvi[0]) {

            coeff[0][j][k] = b8ndvi0[i];
            coeff[1][j][k] = b8ndvi1[i];
            coeff[2][j][k] = b8ndvi2[i];
//          printf("coeffs8: %d %d %d %f %f\n",i, j, k, coeff[2][j][k], b8ndvi2[i]);

    	  }
    	  else {

            coeff[0][j][k] = b1ndvi0[i];
            coeff[1][j][k] = b1ndvi1[i];
            coeff[2][j][k] = b1ndvi2[i];
            coeff[3][j][k] = b1ndvi3[i];
            coeff[4][j][k] = b1ndvi4[i];
            coeff[5][j][k] = b1ndvi5[i];
            coeff[6][j][k] = b1ndvi6[i];
            coeff[7][j][k] = b1ndvi7[i];
            coeff[8][j][k] = b1ndvi8[i];
            coeff[9][j][k] = b1ndvi9[i];
//          printf("coeffs1: %d %d %d %f %f\n",i, j, k, coeff[2][j][k], b1ndvi2[i]);

    	  }

          i++;

        }
      }

//    printf("\nband 1, 8 coefficients initialized \n");

/*    Compute thresholds for current pixel as functions of background
      NDVI and scattering angle.
*/

//    Check for special cases and fill values ('ndvibk'=32.767).
      if(pxin.ndvibk < fill_ndvi1[0] && pxin.ndvibk > fill_ndvi2[0]) {

        if(pxin.ndvibk < ndvi_bnd1[0]) {
          indvi = 0;
          x = 0.0;
          y2 = 0.0;
          interp = 0;
        }
        else if(pxin.ndvibk >= ndvi_bnd2[0]) {
          indvi = 9;
          x = 0.0;
          y2 = 0.0;
          interp = 0;
        }
        else {
          interp = 1;
        }

        if(interp) {
          a = ((pxin.ndvibk / delta_ndvi_bin) - 0.5);
          if(a < 0) {
            indvi = 0;
          }
          else {
            indvi = a;
          }
//        indvi = min(nbin_ndvi, max(0, int((pxin.ndvibk / delta_ndvi_bin) - 0.5)) );

          x1 = delta_ndvi_bin*indvi + delta_ndvi_bin/2.0;
          x2 = x1 + delta_ndvi_bin;
          x = (pxin.ndvibk - x1) / ( x2 - x1 );
          if(x < 0.0) x = 0.0;
          if(x > 1.0) x = 1.0;

        }

        for( icnf=0; icnf<3; icnf++) {

          y1 = coeff[indvi][icnf][0] +
               (coeff[indvi][icnf][1] * pxin.sctang) +
               (coeff[indvi][icnf][2] * powf(pxin.sctang,2)) +
               (coeff[indvi][icnf][3] * powf(pxin.sctang,3));

          if(interp) {
            y2 =  coeff[indvi+1][icnf][0] +
                  (coeff[indvi+1][icnf][1] * pxin.sctang) +
                  (coeff[indvi+1][icnf][2] * powf(pxin.sctang,2)) +
                  (coeff[indvi+1][icnf][3] * powf(pxin.sctang,3));
          }

          thr[icnf] = ((1.0 - x) * y1) + (x * y2);
          if(pxin.ndvibk < des_ndvi[0]) {
        	  thr_adj[icnf] = thr[icnf] * thr_adj_fac_desert[0];
          }
          else {
            thr_adj[icnf] = thr[icnf] * thr_adj_fac_land[0];
          }

//        printf("ndvi thr: %f %f %f %f %f %f %f\n", des_ndvi[0], fill_ndvi1[0], fill_ndvi2[0],
//                    ndvi_bnd1[0], ndvi_bnd2[0], thr_adj_fac_desert[0], thr_adj_fac_land[0]);
//        printf("coeffs: %f %f %f %f %f %f %f %f\n", coeff[indvi][icnf][0], coeff[indvi][icnf][1],
//        		  coeff[indvi][icnf][2], coeff[indvi][icnf][3], coeff[indvi+1][icnf][0], coeff[indvi+1][icnf][1],
//        		  coeff[indvi+1][icnf][2], coeff[indvi+1][icnf][3]);
//        printf("interpolation: %d %d %f %f %f %f %f\n", indvi, icnf, pxin.sctang, y1, y2, thr[icnf], thr_adj[icnf]);

        }

        *(hicut) = (thr[0] + thr_adj[0]) / 100.0;
        *(midpt) = (thr[1] + thr_adj[1]) / 100.0;
        *(locut) = (thr[2] + thr_adj[2]) / 100.0;
        *(power) = 2.0;

        irtn = 0;

      }
      else {

//      No NDVI-dependent thresholds can be derived for the current pixel.
    	irtn = -1;

      }

//    printf("\nSubroutine get_b1_thresholds\n");
//    printf("%d %d %d %f %f %f %f %f %f %f %f \n", irtn, interp, indvi, pxin.ndvibk, pxin.sctang,
//            x1, x2, x, *(locut), *(midpt), *(hicut));

      return irtn;

}
