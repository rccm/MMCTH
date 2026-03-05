/*
************************************************************************
!F77 

!Description:
  Copies cloud mask result and processing flags for use in subsequent 
  250-m cloud mask processing.

!Input parameters found in pixel.h for current pixel:
  pxout.confidence   cloud mask result
  pxin.land          land background
  pxin.day           sza < 85
  pxin.ice           ice background
  pxin.snow          snow background
  pxin.snglnt        sun glint
  pxin.coast         coast background
  pxin.desert        desert background
  pxin.visusd        visible/NIR data used
  pxout.qabits[10]   QA bit array
  pxout.testbits[6]  cloud mask output bit array

!Output parameters found in pixel.h for use in the 250-m cloud mask
  algorithm after processing the NEXT 1-km pixel.

!Revision History:
  R. Frey      10-26-09        Created

!Team-Unique Header:

!References and Credits:

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include "pixel.h"
      
void save_250m_info()
      
{      

    const int num_cm_bytes = 6;
    const int num_qa_bytes = 10;
    int outbyte;

    pxin.prev_confidence = pxout.confidence;
    pxin.prev_sza = pxin.sza;
    pxin.prev_land = pxin.land;
    pxin.prev_day = pxin.day;
    pxin.prev_ice = pxin.ice;
    pxin.prev_snow = pxin.snow;
    pxin.prev_sunglint = pxin.sunglint;
    pxin.prev_coast = pxin.coast;
    pxin.prev_desert = pxin.desert;
    pxin.prev_visusd = pxin.visusd;
    pxin.prev_qa1km = pxout.qa1km;

    for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
      pxin.prev_testbits[outbyte] = pxout.testbits[outbyte];
    }
    for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
      pxin.prev_qabits[outbyte] = pxout.qabits[outbyte];
    }
//  printf("Pos 2: %d %d %d\n", pxin.prev_qabits[0], pxin.prev_qabits[4], pxin.prev_qabits[5]);
   
}
