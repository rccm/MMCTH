/*
************************************************************************
!F77 
 
!Description:
  Adaptation of  MOD_PR35 code Prepare_ECS_Metadata

!Input Parameters:
 stats_out - the cloud mask output statistics

!Output Parameters:
  None

!Revision History:
  K Baggett    09  Converted to C

!Team-Unique Header:

!References and Credits:

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stats.h"

void prepareMetadata(STATS* stats_out)
      
{      
   float PSAValue[numPSA];
   int no250m = 0;
   int i;
/*
   char PSAName[numPSA][30] = {
     "SuccessfulRetrievalPct",
     "VeryHighConfidentClearPct",
     "HighConfidentClearPct",
     "UncertainConfidentClearPct",
     "LowConfidentClearPct",
     "CloudCoverPct250m",
     "ClearPct250m",
     "DayProcessedPct",
     "NightProcessedPct",
     "SunglintProcessedPct",
     "Snow_IceSurfaceProcessedPct",
     "LandProcessedPct",
     "WaterProcessedPct",
     "ThinCirrusSolarFoundPct",
     "ThinCirrusIR_FoundPct",
     "NonCloudObstructionFoundPct",
     "MinSolarZenithAngle",
     "MaxSolarZenithAngle"
   };
*/

   PSAValue[0] = 100.0 - 100.0*(stats_out->numBad/(float)stats_out->numPixels);
   PSAValue[1] = 100.0*(stats_out->numVeryHigh/(float)stats_out->numPixels);
   PSAValue[2] = 100.0*(stats_out->numHigh/(float)stats_out->numPixels);
   PSAValue[3] = 100.0*(stats_out->numUncertain/(float)stats_out->numPixels);
   PSAValue[4] = 100.0*(stats_out->numLow/(float)stats_out->numPixels);

   if (stats_out->num250Pixels > 0){
      PSAValue[5] = 100.0*(stats_out->num250Cloudy/(float)stats_out->num250Pixels);
      PSAValue[6] = 100.0*(stats_out->num250Clear/(float)stats_out->num250Pixels);
   }
   else {
      PSAValue[5] = 0; 
      PSAValue[6] = 0; 
      no250m = 1;
   }
   
   PSAValue[7] = 100.0*(stats_out->numDay/(float)stats_out->numPixels);
   PSAValue[8] = 100.0*(stats_out->numNight/(float)stats_out->numPixels);
   PSAValue[9] = 100.0*(stats_out->numSunglint/(float)stats_out->numPixels);
   PSAValue[10] = 100.0*(stats_out->numIce/(float)stats_out->numPixels);
   PSAValue[11] = 100.0*(stats_out->numLand/(float)stats_out->numPixels);
   PSAValue[12] = 100.0*(stats_out->numWater/(float)stats_out->numPixels);
   PSAValue[13] = 100.0*(stats_out->numThinCirrusSolar/(float)stats_out->numPixels);
   PSAValue[14] = 100.0*(stats_out->numThinCirrusIR/(float)stats_out->numPixels);
   PSAValue[15] = 100.0*(stats_out->numNonCloudObstruction/(float)stats_out->numPixels);
   PSAValue[16] = stats_out->min_sza;
   PSAValue[17] = stats_out->max_sza;

// printf("stats: %d %d %d %d \n", stats_out->numPixels, stats_out->numThinCirrusSolar, stats_out->numThinCirrusIR,
//         stats_out->numNonCloudObstruction);

/*
   fprintf(stdout,"\nnumPixels=%d,num250Pixels=%d\n",stats_out->numPixels,stats_out->num250Pixels);
   fprintf(stdout,"pctSuccessfulRetrieval=%6.2f,pctDay=%6.2f,pctNight=%6.2f,numSunglint=%6.2f,numIce=%6.2f,numLand=%6.2f,numWater=%6.2f\n",
      PSAValue[0],PSAValue[7],PSAValue[8],PSAValue[9],PSAValue[10],
      PSAValue[11],PSAValue[12]);
   fprintf(stdout,"pct250Cloudy=%6.2f,pct250Clear=%6.2f\n",
      PSAValue[5],PSAValue[6]);
   fprintf(stdout,"pctVeryHigh=%6.2f,pctHigh=%6.2f,pctUncertain=%6.2f,pctLow=%6.2f\n",
      PSAValue[1],PSAValue[2],PSAValue[3],PSAValue[4]);
   fprintf(stdout,"pctThinCirrusSolar=%6.2f,pctThinCirrusIR=%6.2f,pctNonCloudObstruction=%6.2f\n",
      PSAValue[13],PSAValue[14],PSAValue[15]);
   fprintf(stdout,"MinSolarZenithAngle=%7.2f,MaxSolarZenithAngle=%7.2f\n",PSAValue[16],PSAValue[17]);
*/

   for (i=0; i<numPSA;++i)
     g_stats[i]=PSAValue[i]; 
}
