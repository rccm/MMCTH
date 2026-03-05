/*
!F77

!DESCRIPTION:
  Finds if current pixel is adjacent to a cloudy pixel.   


!Input parameters:
   none

!Output Parameters:
   none    
    
   Output stored in output "testbit" array; bit #12

!Revision History:
   original version 12/07

!Team-unique Header:

!References and Credits:

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>
#include "pixel.h"
#include "mask_processing_constants.h"
#include "granule.h"


void get_cloud_adj()


{

/**********************************************************************/
	
/*  Declarations */
	
    extern int check_line_edge(int, int, int);
    extern int check_elem_edge(int, int, int);
    extern int check_bits(int, unsigned char[]);
    extern void set_bit(int, unsigned char[]);
    
	int imvl, imve, nl, ne;
	int i,k;
	int line, elem, idx;
	int line_edge, elem_edge;
	int bit_test2;
	int bit_test0;
	int cloud;
	
	long int j;
	long int ide;

	unsigned char adj_testbits[6];
	unsigned char current_testbits[6];
	unsigned char currentqa_testbits[10];
	
/**********************************************************************/	
	
/*  Initializations */
	
	line = 0;
	elem = 0;
	idx = 0;
	
/**********************************************************************/    

//  Cycle through the granule.
	
    for (line=0; line<nscans; ++line) {
      for (elem=0; elem<neles; ++elem) {
	
        current_testbits[0] = *(g_cm[0] + idx);
        current_testbits[1] = *(g_cm[1] + idx);
        currentqa_testbits[1] = *(g_qa[1] + idx);

//      Check for valid cloud mask.
        bit_test0 = check_bits(0, current_testbits);
        if(bit_test0) {

          (void) set_bit(12, currentqa_testbits); 
		  
//        Check for data border pixels.
	  line_edge = check_line_edge(line, nscans, lines_in_edge);
	  elem_edge = check_elem_edge(elem, neles, elems_in_edge);
		  
          if(line_edge && line == 0) {
            nl = 2;
            imvl = 0;
          }
          else if(line_edge && line > 0) {
            nl = 2;
            imvl = 1;
          }
          else {
            nl = nlcntx;
            imvl = 1;
          }
    
          if(elem_edge && elem == 0) {
            ne = 2;
            imve = 0;
          }
          else if(elem_edge && elem > 0) {
            ne = 2;
            imve = 1;
          }
          else {
            ne = necntx;
            imve = 1;
          }
        
//        Get cloud adjacency for current pixel.
    
          cloud = 0;
          for(i=0; i<nl; ++i) {
            for(k=0; k<ne; ++k) {
    	
              j = ( (line + (i - imvl)) * neles) + elem;
              ide = j + (k - imve);
            
              adj_testbits[0] = *(g_cm[0] + ide);
            
              bit_test0 = check_bits(0, adj_testbits);
              bit_test2 = check_bits(2, adj_testbits);
              if(bit_test0 && !bit_test2) cloud = 1;
          
            }
          }
        
          if(!cloud) (void) set_bit(12, current_testbits);
/*
          printf("\nvalue of cm byte 0-7 %d\n", current_testbits[0]);
          printf("value of cm byte 1-7 %d\n", current_testbits[1]);
          printf("value of cm byte 2-7 %d\n", current_testbits[2]);
          printf("value of cm byte 3-7 %d\n", current_testbits[3]);
          printf("value of cm byte 4-7 %d\n", current_testbits[4]);
          printf("value of cm byte 5-7 %d\n", current_testbits[5]);
*/        
          g_cm[1][idx] = current_testbits[1];
          g_qa[1][idx] = currentqa_testbits[1];

        }
        
        idx++;
        
      }
    }	  
    
    return;
        
}
