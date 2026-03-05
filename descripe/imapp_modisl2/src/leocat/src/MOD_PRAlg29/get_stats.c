/*
************************************************************************
!F77 
 
!Description:
  Routine for summing up the various qa granule statistics to be
    filed as PSA's (Product Specific Attributes)

!Input Parameters:
  None

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
#include "pixel.h"
#include "stats.h"
#include "mask_processing_constants.h"
      
void get_stats(STATS* stats_out)
      
{      
   stats_out->numPixels++;

   if (!pxout.qa1km)
     stats_out->numBad++;

   if (pxin.day)
     stats_out->numDay++;
   
   if (pxin.night)
     stats_out->numNight++;

   if (pxin.sunglint)
     stats_out->numSunglint++;

   if (pxin.snow || pxin.ice)
     stats_out->numIce++;

   if (pxin.land || pxin.coast || pxin.desert)
     stats_out->numLand++;

   if (pxin.water)
     stats_out->numWater++;

   if (pxout.confidence > 0.99)
      stats_out->numVeryHigh++;
   else if (pxout.confidence > 0.95)
      stats_out->numHigh++;
   else if (pxout.confidence > 0.66)
      stats_out->numUncertain++;
   else
      stats_out->numLow++;
   
   if (pxin.thinCirrusSolar)
      stats_out->numThinCirrusSolar++;

   if (pxin.thinCirrusIR)
      stats_out->numThinCirrusIR++;
  
   if (pxin.nonCloudObstruction)
      stats_out->numNonCloudObstruction++;

   /*fprintf(stdout,"numPixels=%d,numDay=%d,numNight=%d,numIce=%d,numLand=%d,numWater=%d\n", 
      stats_out->numPixels,
      stats_out->numDay,stats_out->numNight,stats_out->numIce,
      stats_out->numLand,stats_out->numWater);*/
}
