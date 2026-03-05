/*
************************************************************************
!F77 

!Description:
  Performs cloud tests for 250-m pixels with land surfaces. 

!Input parameters: 

  pixel index                 0-15
  qkm_ref1                    band 1 250-m reflecatances
  qkm_ref2                    band 2 250-m reflecatances
  qkm_cm                      1-km cloud mask result
  qkm_snow                    surface snow    
  qkm_day                     sza < 85 degrees 
  qkm_visusd                  visible/NIR data used
  qkm_coast                   coastal region       
  qkm_desert                  desert region       
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
      
void land_day_250m(int npx, float qkm_ref1[16], float qkm_ref2[16], int qkm_cm[16],
                    int qkm_snow[16], int qkm_day[16], int qkm_visusd[16],
                    int qkm_coast[16], int qkm_desert, int qkm_qa1km[16])

{

    extern void set_bit(int, unsigned char[]);

    int i;
    int bitno;

    float m01;
    float m02;
    float s1, s2, etan, etad, eta, vrat;

    i = npx;
    bitno = 32 + i;

    if(qkm_qa1km[i] && qkm_day[i] && qkm_visusd[i]) {

      if( qkm_cm[i]) {

//      Assume 250 m pixel is clear - fill bit with 1 km result.
        if( rint(qkm_ref1[i]) != rint(bad_data) ) {
          (void) set_bit(bitno, pxout.temp_qabits);
          (void) set_bit(bitno, pxout.temp_testbits);
        }

      }
      else if( qkm_snow[i]) {

//      Assume 250 m pixel is cloudy - but do not perform individual
//      pixel tests, so qa bit is set and test bit stays unset (0).
//      Bit agrees with 1-km result.
        if( rint(qkm_ref1[i]) != rint(bad_data) ) {
          (void) set_bit(bitno, pxout.temp_qabits);
        }

      }
      else {

//      Perform tests on current 250-m pixel. Look for clear sky.

        m01 = qkm_ref1[i];
        m02 = qkm_ref2[i];
        s1 = m01 * 100.0;
        s2 = m02 * 100.0;

        if( rint(m01) != rint(bad_data) && rint(m02) != rint(bad_data) ) {

//        Perform modified GEMI test.

          (void) set_bit(bitno, pxout.temp_qabits);

          etan = 2.0 * (s2-s1) + 1.5*s2 + 0.5*s1;
          etad = s2 + s1 + 0.5;
          eta = etan / etad;
          vrat = eta * (1.0-0.25*eta) - ((s1-0.125) / (1.0-s1));

//        printf("Mod-GEMI: %d %f %f %f %f\n", i, m01, m02, vrat, dlvrat[2]);

          if(vrat > dlvrat[2]) {
//          Clear land.
            (void) set_bit(bitno, pxout.temp_testbits);
          }
          else if(qkm_coast[i] == 0) {
//          Can't be sure of clear, so result is same as 1-km mask (cloudy).
//          Bit stays unset.
          }
          else {
//          Region is coastal. Perform visible threshold test using water threshold.
            if(m02 < doref2[2]) {
              (void) set_bit(bitno, pxout.temp_testbits);
            }
            else {
//            Can't be sure of clear, so result is same as 1-km mask (cloudy).
//            Bit stays unset.
            }  
          }
            
        }

      }

    }

//  printf("250m land: %d %d %d %d %d %d %d %f %f\n\n", i, qkm_qa1km[i], qkm_day[i], qkm_visusd[i],
//         qkm_snow[i], qkm_cm[i], qkm_coast[i], m02, doref2[2]);



}
