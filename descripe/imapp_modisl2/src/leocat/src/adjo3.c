#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void adjo3(double *p_prof, double *o3_tmp, float *o3_prof, float *adjo3_prof)
{

//    Define stratospheric interpolation points (pressure levels).
      const int gdaspp[] = {18, 20, 26, 29, 35, 39, 44, 46};

//    Define stratosphere boundaries (pressire levels).
      const int strat_top = 20;
      const int strat_bot = 44;

      static int init = 0;
      int levc, intpt, toplev, botlev, kk, jj, m, n;
      float do3, dptot, dp;

/***************************************************************************************/

//    Copy original O3 profile to adjusted array.
      for(kk=0; kk<101; kk++) {
        adjo3_prof[kk] = o3_prof[kk];
      }

/*    if(init == 0) {
        for(kk=0; kk<101; kk++) {
          printf("climo o3: %d %f %f\n", kk, o3_prof[kk], adjo3_prof[kk]);
        }
        for(kk=0; kk<6; kk++) {
          printf("GDAS o3: %d %f\n", kk, o3_tmp[kk]);
        }
      }
*/
/***************************************************************************************/

//    Adjust climatology stratospheric ozone profile.
//    First, replace climatological values with those from GDAS (6 values) and convert
//    from kg/kg dry air to ppmv.
      adjo3_prof[20] = o3_tmp[0] * 1000000 * (28.97 / 48.00);
      adjo3_prof[26] = o3_tmp[1] * 1000000 * (28.97 / 48.00);
      adjo3_prof[29] = o3_tmp[2] * 1000000 * (28.97 / 48.00);
      adjo3_prof[35] = o3_tmp[3] * 1000000 * (28.97 / 48.00);
      adjo3_prof[39] = o3_tmp[4] * 1000000 * (28.97 / 48.00);
      adjo3_prof[44] = o3_tmp[5] * 1000000 * (28.97 / 48.00);

//    for(kk=0; kk<102; kk++) {
//      printf("o3 profiles 2: %d %f %f\n", kk, o3_prof[kk], adjo3_prof[kk]);
//    }

//    Interpolate stratospheric values.
      levc = 1;
      intpt = 0;
      toplev = strat_top - 2;
      botlev = strat_bot + 2;

      for (kk=toplev; kk<botlev+1; kk++) {

//      printf("kk=%d, levc=%d, intpt=%d, O3=%f %f\n", kk,levc,intpt,o3_prof[gdaspp[levc]],
//              o3_prof[gdaspp[intpt]]);

        if(gdaspp[levc] == kk) {

          do3 = adjo3_prof[gdaspp[levc]] - adjo3_prof[gdaspp[intpt]];
          dptot = log(p_prof[gdaspp[levc]]) - log(p_prof[gdaspp[intpt]]);

//        printf("kk=%d, levc=%d, intpt=%d, adjo3_prof=%f, %f, do3=%f, dptot=%f\n",
//                kk,levc,intpt,adjo3_prof[gdaspp[levc]], adjo3_prof[gdaspp[intpt]],
//                do3,dptot);

          m = gdaspp[intpt] + 1;
          n = gdaspp[levc] - 1;

          for(jj=m; jj<n+1; jj++) {

            dp = log(p_prof[jj]) - log(p_prof[gdaspp[intpt]]);
            adjo3_prof[jj] = adjo3_prof[gdaspp[intpt]] + (do3 * (dp / dptot));
//          if(init ==0)
//            printf("jj=%d, o3_prof=%f, adjo3_prof=%f, dp=%f, dp/dptot=%f\n", jj,o3_prof[jj],
//                adjo3_prof[jj],dp,dp/dptot);

          }

          levc += 1;
          intpt += 1;

        }

      }

/*
      if(init == 0) {
        for(kk=0; kk<101; kk++) {
          printf("adjusted o3: %d %f %f\n", kk, o3_prof[kk], adjo3_prof[kk]);
        }
      }
      init = 1;
*/

}
