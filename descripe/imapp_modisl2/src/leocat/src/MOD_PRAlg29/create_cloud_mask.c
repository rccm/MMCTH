/*
!C ********************************************************************
!Description:

  Integer function create_cloud_mask.c
  Creates the MODIS cloud mask (MOD35).
  Called from modis_cm_main.c

!Input arguments:
  none

  Inputs from granule.h, pixel.h

!Output arguments:
  none

  int return_code   successful completion is zero, otherwise non-zero
                    (returned through function call)

  Output cloud mask and associated products and arrays in ...

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "granule.h"
#include "pixel.h"
#include "thresholds.h"
#include "stats.h"


int create_cloud_mask()


{

/*  Declarations */

    extern int get_pxldat(STATS *);
    extern int perform_cloud_tests();
    extern void thin_cirrus_tests();
    extern void noncloud_obs_tests();
    extern void set_confdnc();
    extern void set_proc_path();
    extern void set_bit(int, unsigned char[]);
    extern void get_cloud_adj();
    extern void save_250m_info();
    extern void create_250m_cm(int, int, int, int, STATS *);
    extern void get_stats(STATS *);
    extern void prepareMetadata(STATS *);
    extern void set_qa_flags(int, unsigned char[]);

    long int idx;
    long int mem;

    int return_code;
    int scan_index;
    int elem_index;
    int beg_elem_index;
    int end_elem_index;
    int ret_code;
    int num_tests_performed;
    int outbyte;
    int out_indx;
    int num_cm_bytes = 6;
    int num_qa_bytes = 10;
    int num_250m_stats = 4;
    int BadGeoCnt = 0;

    STATS stats_out = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

/*  Initializations */

    printf("Executing create_cloud_mask 29\n");
    return_code = 0;

//  Allocate memory for cloud mask output.
    mem = neles * nscans * sizeof(unsigned char);
    for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
      granule_testbits[outbyte] = (unsigned char *) malloc(mem);
    }
    for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
      granule_qabits[outbyte] = (unsigned char *) malloc(mem);
    }
//  Allocate memnory for 250-m stats output.
    mem = neles * nscans * sizeof(float);
    for (out_indx=0; out_indx<num_250m_stats; ++out_indx) {
      granule_250m_stats[out_indx] = (float *) malloc(mem);
    }

/*  Compute cloud mask */

/*  Loop over scan lines in data granule */

    scan_index = 0;
//  scan_index =  290;
    idx = 0;
    BadGeoCnt = 0;
    while ( scan_index < nscans) {
//  while ( scan_index <  291) {

      pxin.scan = scan_index;

/*    Loop over elements in scan line */

      elem_index = 0;
//    elem_index = 645;
      beg_elem_index = elem_index;

      end_elem_index = neles;
//    end_elem_index = 1200;

      while ( elem_index < end_elem_index) {

        pxin.elem = elem_index;

/*      Fill input data structure for current pixel. */
        ret_code = get_pxldat(&stats_out);

/*      Perform cloud tests. */
        num_tests_performed = perform_cloud_tests();

/*      Perform thin cirrus tests. These tests do not contribute to the final
        confidence of clear sky. */
        (void) thin_cirrus_tests();

/*      Perform non-cloud obstruction tests.  These tests do not contribute to the
        confidence of clear sky. */
        if(!pxin.snow && !pxin.ice) 
          (void) noncloud_obs_tests();

/*      Set output for final confidence of clear sky. */
        (void) set_confdnc();

/*      Set output for processing path information. */
        (void) set_proc_path();

/*      Set cloud mask "valid" flag. */
        if(num_tests_performed >= 1) {
          (void) set_bit(0, pxout.testbits);
          pxout.qa1km = 1;
        }
        else {
          pxout.qa1km = 0;
        }

/*      Set cloud mask QA flags. */
        (void) set_qa_flags(num_tests_performed, pxout.qabits);

/*      Zero out results of the cloud mask if geolocation data are missing. */
        if(pxin.bad_geo) {

	  /*Added by G. Britzolakis on 12/08/11*/
	  if(BadGeoCnt<100){ 
	    printf("Bad geo-location found at %d %d\n", scan_index, elem_index);
	  }
	  BadGeoCnt++;
	  /*Added by G. Britzolakis on 12/08/11*/

          pxout.qa1km = 0;
          for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
            pxout.testbits[outbyte] = 0;
          }
          for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
            pxout.qabits[outbyte] = 0;
          }
        }

/*      Get the output statistics, */
        (void) get_stats(&stats_out);

/*      Store cloud mask and qa info. */
        for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
          granule_testbits[outbyte][idx] = pxout.testbits[outbyte];
        }
        for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
          granule_qabits[outbyte][idx] = pxout.qabits[outbyte];
        }

/*************************************************************************************************/

/*      Process 250-m data. */

        if(proc_qkm && elem_index > beg_elem_index)
          (void)create_250m_cm(idx, end_elem_index, elem_index, scan_index, &stats_out);

        if(proc_qkm) (void)save_250m_info();

//*************************************************************************************************

        idx++;
        elem_index++;
/*
        printf("\nvalue of cm byte 0 %d\n", pxout.testbits[0]);
        printf("value of cm byte 1 %d\n", pxout.testbits[1]);
        printf("value of cm byte 2 %d\n", pxout.testbits[2]);
        printf("value of cm byte 3 %d\n", pxout.testbits[3]);
        printf("value of cm byte 4 %d\n", pxout.testbits[4]);
        printf("value of cm byte 5 %d\n", pxout.testbits[5]);
*/
      }

      scan_index++;

    }
    /* Added by G. Britzolakis on 12/08/11*/
    printf("Number of Bad geo-location found: %d\n", BadGeoCnt);
    /* Added by G. Britzolakis on 12/08/11*/

/*  Prepare cloud mask statistics for output metadata. */
    (void) prepareMetadata(&stats_out);

/*  Copy cloud mask and qa info to output arrays. */
    for (outbyte=0; outbyte<num_cm_bytes; ++outbyte) {
      g_cm[outbyte] = granule_testbits[outbyte];
    }
    for (outbyte=0; outbyte<num_qa_bytes; ++outbyte) {
      g_qa[outbyte] = granule_qabits[outbyte];
    }
    if(proc_qkm) {
      for (out_indx=0; out_indx<num_250m_stats; ++out_indx) {
        g_250m_stats[out_indx] = granule_250m_stats[out_indx];
      }
    }

/*  Get cloud adjacency flag. */
    (void) get_cloud_adj();
/*
    printf("\nvalue of cm byte 0 %d\n", *(g_cm[0]));
    printf("value of cm byte 1 %d\n", *(g_cm[1]));
    printf("value of cm byte 2 %d\n", *(g_cm[2]));
    printf("value of cm byte 3 %d\n", *(g_cm[3]));
    printf("value of cm byte 4 %d\n", *(g_cm[4]));
    printf("value of cm byte 5 %d\n", *(g_cm[5]));
*/
/*
    printf("\nvalue of qa byte 0 %d\n", *(g_qa[0]));
    printf("value of qa byte 1 %d\n", *(g_qa[1]));
    printf("value of qa byte 2 %d\n", *(g_qa[2]));
    printf("value of qa byte 3 %d\n", *(g_qa[3]));
    printf("value of qa byte 4 %d\n", *(g_qa[4]));
    printf("value of qa byte 5 %d\n", *(g_qa[5]));
    printf("value of qa byte 6 %d\n", *(g_qa[6]));
    printf("value of qa byte 7 %d\n", *(g_qa[7]));
    printf("value of qa byte 8 %d\n", *(g_qa[8]));
    printf("value of qa byte 9 %d\n", *(g_qa[9]));
*/
    return return_code;

}
