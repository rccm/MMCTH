#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int adjo3_(double *p_prof, double *o3_tmp, float *o3_prof, float *adjo3_prof)
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
	
	
//    Interpolate stratospheric values.
      levc = 1;
      intpt = 0;
      toplev = strat_top - 2;
      botlev = strat_bot + 2;

      for (kk=toplev; kk<botlev+1; kk++) {


        if(gdaspp[levc] == kk) {

          do3 = adjo3_prof[gdaspp[levc]] - adjo3_prof[gdaspp[intpt]];
          dptot = log(p_prof[gdaspp[levc]]) - log(p_prof[gdaspp[intpt]]);


          m = gdaspp[intpt] + 1;
          n = gdaspp[levc] - 1;

          for(jj=m; jj<n+1; jj++) {

            dp = log(p_prof[jj]) - log(p_prof[gdaspp[intpt]]);
            adjo3_prof[jj] = adjo3_prof[gdaspp[intpt]] + (do3 * (dp / dptot));

          }

          levc += 1;
          intpt += 1;

        }

      }
	return 0;

}
