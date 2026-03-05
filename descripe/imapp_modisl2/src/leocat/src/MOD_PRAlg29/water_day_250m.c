/*
************************************************************************
!F77 

!Description:
  Performs cloud tests for 250-m pixels with water surfaces. 

!Input parameters: 

  pixel index                 0-15
  qkm_ref1                    band 1 250-m reflecatances
  qkm_ref2                    band 2 250-m reflecatances
  qkm_sunglint                geometric sunglint region
  qkm_ice                     surface ice     
  qkm_cm                      1-km cloud mask result
  qkm_day                     sza < 85 degrees 
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
#include <math.h>
#include "pixel.h"
#include "thresholds.h"
#include "mask_processing_constants.h"

      
void water_day_250m(int npx, float qkm_ref1[16], float qkm_ref2[16], int qkm_sunglint[16],
                    int qkm_ice[16], int qkm_cm[16], int qkm_day[16],
                    int qkm_visusd[16], int qkm_qa1km[16])

{

    extern void set_bit(int, unsigned char[]);

    int i;
    int bitno;

    i = npx;
    bitno = 32 + i;

    if(qkm_qa1km[i] && qkm_day[i] && qkm_visusd[i]) {

      if( rint(qkm_ref2[i]) != rint(bad_data) ) {  

        (void) set_bit(bitno, pxout.temp_qabits);

        if( (qkm_ref2[i] >= doref2[1] && qkm_ref2[i] <= doref2[0]) ||
            (qkm_sunglint[i] && qkm_ref2[i] >= doref2[1]) ||
             qkm_ice[i]) {

//        Don't perform 250m tests - copy 1-km result.
          if(qkm_cm[i]) (void) set_bit(bitno, pxout.temp_testbits);

        }
        else if(qkm_ref2[i] < doref2[1]) {

          (void) set_bit(bitno, pxout.temp_testbits);

        }

      }

    }

//  printf("250m water: %d %d %d %d %d %d %d %f %f %f %f\n\n", i, qkm_qa1km[i], qkm_day[i], qkm_visusd[i],
//         qkm_sunglint[i], qkm_ice[i], qkm_cm[i], qkm_ref2[i], doref2[0], doref2[1], doref2[2]);



}
