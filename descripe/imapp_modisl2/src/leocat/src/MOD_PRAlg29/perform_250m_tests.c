/*
************************************************************************
!F77 

!Description:
  Sets up calls to routines that perform 250-m cloud tests based on 
  land/water flag.  Each call to this routine processes a 4x4 pixel 
  group of 250m pixels.

!Input parameters: 

  qkm_ref1                    band 1 250-m reflecatances
  qkm_ref2                    band 2 250-m reflecatances
  qkm_cm                      1-km cloud mask result
  qkm_ice                     surface ice     
  qkm_snow                    surface snow
  qkm_desert                  arid or semi-arid surface conditions
  qkm_coast                   coastal region  
  qkm_land                    land surface  
  qkm_day                     sza < 85 degrees 
  qkm_sunglint                geometric sunglint region
  qkm_visusd                  visible/NIR data used
  qkm_qa1km                   valid 1-km cloud mask

!Output parameters found in pixel.h.

!Revision History:
  R. Frey      10-28-09        Created

!Team-Unique Header:

!References and Credits:

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include "granule.h"
#include "pixel.h"
#include "mask_processing_constants.h"
      
void perform_250m_tests(float qkm_ref1[16], float qkm_ref2[16], int qkm_cm[16],
                        int qkm_ice[16], int qkm_snow[16], int qkm_desert[16],
                        int qkm_coast[16], int qkm_land[16], int qkm_day[16],
                        int qkm_sunglint[16], int qkm_visusd[16], int qkm_qa1km[16])

{

    extern void land_day_250m(int, float[], float[], int[], int[], int[], int[], int[], int[], int[]);
    extern void water_day_250m(int, float[], float[], int[], int[], int[], int[], int[], int[]);

    const int n250p = 16;

    int i;
 
//  Perform tests.

    for(i=0; i<n250p; i++) {

//    printf("Performing 250m tests: %d %d %d %d %d %d %d %d\n", i, qkm_land[i], qkm_qa1km[i], qkm_day[i], qkm_visusd[i],
//             qkm_sunglint[i], qkm_ice[i], qkm_cm[i]);

      if(qkm_land[i]) {

        (void)land_day_250m(i, qkm_ref1, qkm_ref2, qkm_cm, qkm_snow, qkm_day, qkm_visusd,
                            qkm_coast, qkm_desert, qkm_qa1km);

      }
      else{

        (void)water_day_250m(i, qkm_ref1, qkm_ref2, qkm_sunglint, qkm_ice, qkm_cm, qkm_day,
                             qkm_visusd, qkm_qa1km);

      }
       
    }


}
