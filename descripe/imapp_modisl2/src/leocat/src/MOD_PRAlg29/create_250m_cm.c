/*
************************************************************************
!F77 

!Description:
  Creates 250-m cloud mask. 

!Input parameters found in pixel.h for the PREVIOUS 1-km pixel. In the case 
 of the last 1-km pixel in a scan, also need these values for the CURRENT 
 pixel (without "prev" prefix).
  pxin.prev_confidence        cloud mask result
  pxin.prev_land              land background
  pxin.prev_day               sza < 85
  pxin.prev_ice               ice background
  pxin.prev_snow              snow background
  pxin.prev_snglnt            sun glint
  pxin.prev_coast             coast background
  pxin.prev_desert            desert background
  pxin.prev_visusd            visible/NIR data used
  pxout.prev_qa1km            valid 1-km cloud mask

 Inputs from granule.h:
  g_qkm1                      pointer to 250-m band 1 reflectances
  g_qkm2                      pointer to 250-m band 2 reflectances

 Inputs from calling arguments:
  end_elem_index              number of elements to process in a 1-km scan
  elem_index                  current 1-km pixel number
  scan_index                  current 1-km scan number

!Output parameters found in pixel.h.

!Revision History:
  R. Frey      10-26-09        Created

!Team-Unique Header:

!References and Credits:

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "granule.h"
#include "pixel.h"
#include "mask_processing_constants.h"
#include "stats.h"
      
void create_250m_cm(int idx, int end_elem_index, int elem_index, int scan_index,
                    STATS* stats_out)
      
{      

    extern void perform_250m_tests(float[], float[], int[], int[], int[], int[],
                                   int[], int[], int[], int[], int[], int[]);

    const int n250p = 16;
    const int num_cm_bytes = 6;
    const int num_qa_bytes = 10;

    float cos_sza;
    float qkm_ref1[16];
    float qkm_ref2[16];

    double den;
    double j1;
    double j2;
    double num;
    double qkm_sum[2];
    double qkm_sumsq[2];
    double sqsum;

    int qkm_cm[16];
    int qkm_ice[16];
    int qkm_snow[16];
    int qkm_desert[16];
    int qkm_coast[16];
    int qkm_land[16];
    int qkm_day[16];
    int qkm_sunglint[16];
    int qkm_visusd[16];
    int qkm_qa1km[16];

    int proc_index;
    int i;
    int ll;
    int kk;
    int outbyte;

//*********************************************************************************************************

//  Process first through next-to-last group (4x4) of 250-m pixels associated with current 1-km scan.

    proc_index = elem_index - 1;

//  printf("Processing 250-m cloud mask %d\n", proc_index);

//*********************************************************************************************************

//  Get 250-m reflectances associated with the previously-processed 1-km pixel.

    cos_sza = cos(pxin.prev_sza * dtr);
//  printf("cos_sza: %f %f %f \n", dtr, pxin.prev_sza, cos_sza);

    i = 0;
    j1 = 0.0;
    j2 = 0.0;
    qkm_sum[0] = 0.0;
    qkm_sum[1] = 0.0;
    qkm_sumsq[0] = 0.0;
    qkm_sumsq[1] = 0.0;

    for (kk=0; kk < 4; kk++){ 
      for (ll=0; ll < 4; ll++) {

        qkm_ref1[i] = g_qkm1[ ((scan_index*4 + ll) * g_qkm_ncol) + (proc_index*4 + kk) ] / cos_sza;
        if( qkm_ref1[i] <= 0.0 || (qkm_ref1[i] > 2.0) ) qkm_ref1[i] = bad_data;

        qkm_ref2[i] = g_qkm2[ ((scan_index*4 + ll) * g_qkm_ncol) + (proc_index*4 + kk) ] / cos_sza;
        if( qkm_ref2[i] <= 0.0 || (qkm_ref2[i] > 1.3 &&
            qkm_ref1[i] == bad_data) ) qkm_ref2[i] = bad_data;
        if( qkm_ref2[i] > 1.3 && qkm_ref1[i] != bad_data) qkm_ref2[i] = 1.0;

        if(qkm_ref1[i] != bad_data) {
          j1 = j1 + 1.0;
          qkm_sum[0] = qkm_sum[0] + (double)qkm_ref1[i];
          qkm_sumsq[0] = qkm_sumsq[0] + ((double)qkm_ref1[i] * (double)qkm_ref1[i]);
        }

        if(qkm_ref2[i] != bad_data) {
          j2 = j2 + 1.0;
          qkm_sum[1] = qkm_sum[1] + (double)qkm_ref2[i];
          qkm_sumsq[1] = qkm_sumsq[1] + ((double)qkm_ref2[i] * (double)qkm_ref2[i]);
        }
  
        i++;

      }
    }

//*********************************************************************************************************

//  Get means and standard deviations over the 16-pixel area for bands 1 and 2.

    if(j1 > 1.0) {
      sqsum = qkm_sum[0] * qkm_sum[0];
      num =  (j1 * qkm_sumsq[0]) - sqsum;
      den = j1 * (j1 - 1.0);
      pxout.qkm_std[0] = (float)sqrt(num / den);
      pxout.qkm_mean[0] = (float)qkm_sum[0] / j1;
    }
    else if(j1 == 1.0) {
      pxout.qkm_std[0] = 0.0;
      pxout.qkm_mean[0] = (float)qkm_sum[0];
    }
    else {
      pxout.qkm_std[0] = bad_data;
      pxout.qkm_mean[0] = bad_data;
    }
//  printf("Band 1 250-m stats %f %f %f %f %f %f %f %f\n", j1,qkm_sum[0],sqsum,qkm_sumsq[0],num,den,
//           pxout.qkm_std[0],pxout.qkm_mean[0]);

    if(j2 > 1.0) {
      sqsum = qkm_sum[1] * qkm_sum[1];
      num =  (j2 * qkm_sumsq[1]) - sqsum;
      den = j2 * (j2 - 1.0);
      pxout.qkm_std[1] = (float)sqrt(num / den);
      pxout.qkm_mean[1] = (float)qkm_sum[1] / j2;
    }
    else if(j2 == 1.0) {
      pxout.qkm_std[1] = 0.0;
      pxout.qkm_mean[1] = (float)qkm_sum[1];
    }
    else {
      pxout.qkm_std[1] = bad_data;
      pxout.qkm_mean[1] = bad_data;
    }
//  printf("Band 2 250-m stats %f %f %f %f %f %f %f %f\n", j2,qkm_sum[1],sqsum,qkm_sumsq[1],num,den,
//              pxout.qkm_std[1],pxout.qkm_mean[1]);

//  printf("250-m refs: %f %f %f %f %f %f %f %f\n", qkm_ref2[0], qkm_ref2[1], qkm_ref2[2], qkm_ref2[3],
//                                                  qkm_ref2[4], qkm_ref2[5], qkm_ref2[6], qkm_ref2[7]);
//  printf("250-m refs: %f %f %f %f %f %f %f %f\n\n", qkm_ref2[8], qkm_ref2[9], qkm_ref2[10], qkm_ref2[11],
//                                                  qkm_ref2[12], qkm_ref2[13], qkm_ref2[14], qkm_ref2[15]);

//  Output 250-m stats.
    granule_250m_stats[0][idx-1] = pxout.qkm_mean[0];
    granule_250m_stats[1][idx-1] = pxout.qkm_std[0];
    granule_250m_stats[2][idx-1] = pxout.qkm_mean[1];
    granule_250m_stats[3][idx-1] = pxout.qkm_std[1];

//*********************************************************************************************************

//  Copy the 1-km cloud mask results for the previous 1-km pixel into temporary arrays.

    for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
      pxout.temp_testbits[outbyte] = pxin.prev_testbits[outbyte];
    }
    for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
      pxout.temp_qabits[outbyte] = pxin.prev_qabits[outbyte];
    }
//  printf("Pos 3: %d %d %d\n", pxout.temp_qabits[0], pxout.temp_qabits[4], pxout.temp_qabits[5]);
    
//*********************************************************************************************************

//  Get 1-km information for each 250-m pixel. The first three columns are associated with the
//  previous 1-km pixel and the last column with the current one.

    for(i=0; i< n250p; i++) {

      stats_out->num250Pixels++;

      if(i < 12) {

        if(pxin.prev_confidence > 0.95) {
          qkm_cm[i] = 1;
          stats_out->num250Clear++;
        }
        else {
          qkm_cm[i] = 0;
          stats_out->num250Cloudy++;
        }
        
        qkm_ice[i] = pxin.prev_ice;
        qkm_snow[i] = pxin.prev_snow;
        qkm_desert[i] = pxin.prev_desert;
        qkm_coast[i] = pxin.prev_coast;
        qkm_land[i] = pxin.prev_land;
        qkm_day[i] = pxin.prev_day;
        qkm_sunglint[i] = pxin.prev_sunglint;
        qkm_visusd[i] = pxin.prev_visusd;
        qkm_qa1km[i] = pxin.prev_qa1km;

      }
      else {

        if(pxout.confidence > 0.95) {
          qkm_cm[i] = 1;
          stats_out->num250Clear++;
        }
        else {
          qkm_cm[i] = 0;
          stats_out->num250Cloudy++;
        }
        
        qkm_ice[i] = pxin.ice;
        qkm_snow[i] = pxin.snow;
        qkm_desert[i] = pxin.desert;
        qkm_coast[i] = pxin.coast;
        qkm_land[i] = pxin.land;
        qkm_day[i] = pxin.day;
        qkm_sunglint[i] = pxin.sunglint;
        qkm_visusd[i] = pxin.visusd;
        qkm_qa1km[i] = pxout.qa1km;

      }

    }

//*********************************************************************************************************

//  Calculate 250-m mask for the current 16 pixels.

    (void)perform_250m_tests(qkm_ref1, qkm_ref2, qkm_cm, qkm_ice, qkm_snow, qkm_desert,
                             qkm_coast, qkm_land, qkm_day, qkm_sunglint, qkm_visusd,
                             qkm_qa1km);

//  Copy the data containing the 250-m results back to the output arrays.
  
    for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
      granule_testbits[outbyte][idx-1] = pxout.temp_testbits[outbyte];
    }
    for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
      granule_qabits[outbyte][idx-1] = pxout.temp_qabits[outbyte];
    }
    
//*********************************************************************************************************
//*********************************************************************************************************

//  In the case of the last group of 250-m pixels on a scan, use 1-km info for all pixels.
   
    if( (end_elem_index - elem_index) == 1) {

      proc_index = elem_index;

//    printf("Processing 250-m cloud mask %d\n", proc_index);

//*********************************************************************************************************

//    Get 250-m reflectances associated with the current 1-km pixel.

      cos_sza = cos(pxin.sza * dtr);
//    printf("cos_sza: %f %f %f \n", dtr, pxin.prev_sza, cos_sza);

      i = 0;
      j1 = 0.0;
      j2 = 0.0;
      qkm_sum[0] = 0.0;
      qkm_sum[1] = 0.0;
      qkm_sumsq[0] = 0.0;
      qkm_sumsq[1] = 0.0;

      for (kk=0; kk < 4; kk++){ 
        for (ll=0; ll < 4; ll++) {

          qkm_ref1[i] = g_qkm1[ ((scan_index*4 + ll) * g_qkm_ncol) + (proc_index*4 + kk) ] / cos_sza;
          if( qkm_ref1[i] <= 0.0 || (qkm_ref1[i] > 2.0) ) qkm_ref1[i] = bad_data;

          qkm_ref2[i] = g_qkm2[ ((scan_index*4 + ll) * g_qkm_ncol) + (proc_index*4 + kk) ] / cos_sza;
          if( qkm_ref2[i] <= 0.0 || (qkm_ref2[i] > 1.3 &&
              qkm_ref1[i] == bad_data) ) qkm_ref2[i] = bad_data;
          if( qkm_ref2[i] > 1.3 && qkm_ref1[i] != bad_data) qkm_ref2[i] = 1.0;

          if(qkm_ref1[i] != bad_data) {
            j1 = j1 + 1.0;
            qkm_sum[0] = qkm_sum[0] + (double)qkm_ref1[i];
            qkm_sumsq[0] = qkm_sumsq[0] + ((double)qkm_ref1[i] * (double)qkm_ref1[i]);
          }

          if(qkm_ref2[i] != bad_data) {
            j2 = j2 + 1.0;
            qkm_sum[1] = qkm_sum[1] + (double)qkm_ref2[i];
            qkm_sumsq[1] = qkm_sumsq[1] + ((double)qkm_ref2[i] * (double)qkm_ref2[i]);
          }

          i++;

        }
      }

//*********************************************************************************************************

//    Get means and standard deviations over the 16-pixel area for bands 1 and 2.

      if(j1 > 1.0) {
        sqsum = qkm_sum[0] * qkm_sum[0];
        num =  (j1 * qkm_sumsq[0]) - sqsum;
        den = j1 * (j1 - 1.0);
        pxout.qkm_std[0] = (float)sqrt(num / den);
        pxout.qkm_mean[0] = (float)qkm_sum[0] / j1;
      }   
      else if(j1 == 1.0) {
        pxout.qkm_std[0] = 0.0;
        pxout.qkm_mean[0] = (float)qkm_sum[0];
      }   
      else {
        pxout.qkm_std[0] = bad_data;
        pxout.qkm_mean[0] = bad_data;
      }   
//    printf("Band 1 250-m stats %f %f %f %f %f %f %f %f\n", j1,qkm_sum[0],sqsum,qkm_sumsq[0],num,den,
//                pxout.qkm_std[0],pxout.qkm_mean[0]);
        
      if(j2 > 1.0) { 
        sqsum = qkm_sum[1] * qkm_sum[1];
        num =  (j2 * qkm_sumsq[1]) - sqsum;
        den = j2 * (j2 - 1.0);
        pxout.qkm_std[1] = (float)sqrt(num / den);
        pxout.qkm_mean[1] = (float)qkm_sum[1] / j2;
      }
      else if(j2 == 1.0) {
        pxout.qkm_std[1] = 0.0;
        pxout.qkm_mean[1] = (float)qkm_sum[1]; 
      }                        
      else {                   
        pxout.qkm_std[1] = bad_data;
        pxout.qkm_mean[1] = bad_data;
      }
//    printf("Band 2 250-m stats %f %f %f %f %f %f %f %f\n", j2,qkm_sum[1],sqsum,qkm_sumsq[1],num,den,
//                    pxout.qkm_std[1],pxout.qkm_mean[1]);
 
//    printf("250-m refs: %f %f %f %f %f %f %f %f\n", qkm_ref2[0], qkm_ref2[1], qkm_ref2[2], qkm_ref2[3],
//                                                    qkm_ref2[4], qkm_ref2[5], qkm_ref2[6], qkm_ref2[7]);
//    printf("250-m refs: %f %f %f %f %f %f %f %f\n\n", qkm_ref2[8], qkm_ref2[9], qkm_ref2[10], qkm_ref2[11],
//                                                    qkm_ref2[12], qkm_ref2[13], qkm_ref2[14], qkm_ref2[15]);

//    Output 250-m stats.
      granule_250m_stats[0][idx] = pxout.qkm_mean[0];
      granule_250m_stats[1][idx] = pxout.qkm_std[0];
      granule_250m_stats[2][idx] = pxout.qkm_mean[1];
      granule_250m_stats[3][idx] = pxout.qkm_std[1];

//*********************************************************************************************************

//    Copy the 1-km cloud mask results for the current 1-km pixel into temporary arrays.

      for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
        pxout.temp_testbits[outbyte] = pxout.testbits[outbyte];
      }
      for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
        pxout.temp_qabits[outbyte] = pxout.qabits[outbyte];
      }
    
//*********************************************************************************************************

//    Get 1-km information for each 250-m pixel.

      for(i=0; i< n250p; i++) {

        stats_out->num250Pixels++;
        if(pxout.confidence > 0.95) {
          qkm_cm[i] = 1;
          stats_out->num250Clear++;
        }
        else {
          qkm_cm[i] = 0;
          stats_out->num250Cloudy++;
        }

        qkm_ice[i] = pxin.ice;
        qkm_snow[i] = pxin.snow;
        qkm_desert[i] = pxin.desert;
        qkm_coast[i] = pxin.coast;
        qkm_land[i] = pxin.land;
        qkm_day[i] = pxin.day;
        qkm_sunglint[i] = pxin.sunglint;
        qkm_visusd[i] = pxin.visusd;
        qkm_qa1km[i] = pxout.qa1km;

      }

//*********************************************************************************************************

//    Calculate the 250-m mask for the current 16 pixels.

      (void)perform_250m_tests(qkm_ref1, qkm_ref2, qkm_cm, qkm_ice, qkm_snow, qkm_desert,
                               qkm_coast, qkm_land, qkm_day, qkm_sunglint, qkm_visusd,
                               qkm_qa1km);
   
    }

//*********************************************************************************************************

//  Copy the data containing the 250-m results back to the output arrays.
  
    for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
      granule_testbits[outbyte][idx] = pxout.temp_testbits[outbyte];
    }
    for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
      granule_qabits[outbyte][idx] = pxout.temp_qabits[outbyte];
    }
    
//*********************************************************************************************************
//*********************************************************************************************************

}
