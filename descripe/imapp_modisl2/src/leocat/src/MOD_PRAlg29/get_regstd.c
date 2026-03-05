/*
!F77

!DESCRIPTION:
  Computes standard deviation of (center) pixel of interest and the     
  surrounding eight pixel values for indicated channel.   


!Input parameters:
   band          MODIS band number (1-36)

!Output Parameters:
   sigma         Standard deviation of indicated MODIS channel
 
   Other outputs stored in structure variable 'rg_var' of type "regional_var"
   in pixel.h (rg_var.mean).
    

!Revision History:
   10/04 Collection 5b   R. Frey
   Added calculation and output of mean value.
   
   11/07 Converted to C  R. Frey

!Team-unique Header:

!References and Credits:
   See Cloud Mask ATBD.

!END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>
#include "pixel.h"
#include "mask_processing_constants.h"
#include "granule.h"


float get_regstd(int band)


{

/**********************************************************************/
	
/*  Declarations */
	
	int imvl, imve, nl, ne;
	int i,k;
	
	long int j;
	long int ide;
	
	float sigma;
	
	double n, sum, sumsq, tbb;
	double a, sqsum, num, den, sig, mn;

/**********************************************************************/	
	
/*  Initializations */
 
    n = 0.0;
    sum = 0.0;
    sumsq = 0.0;
    sqsum = 0.0;
    num = 0.0;
    den = 0.0;
    
/**********************************************************************/    

    if(pxin.line_edge && pxin.scan == 0) {
      nl = 2;
      imvl = 0;
    }
    else if(pxin.line_edge && pxin.scan > 0) {
      nl = 2;
      imvl = 1;
    }
    else {
      nl = nlcntx;
      imvl = 1;
    }
    
    if(pxin.elem_edge && pxin.elem == 0) {
      ne = 2;
      imve = 0;
    }
    else if(pxin.elem_edge && pxin.elem > 0) {
      ne = 2;
      imve = 1;
    }
    else {
      ne = necntx;
      imve = 1;
    }

//  Get standard deviation about central pixel.
    
    for(i=0; i<nl; ++i) {
      for(k=0; k<ne; ++k) {
    	
        j = ( (pxin.scan + (i - imvl)) * neles) + pxin.elem;
        ide = j + (k - imve);  
	    tbb = (double)(*( ( g_rad[band-1] ) + ide ));
//      printf("tbb: %ld %ld %f \n", j, ide, tbb);
	    
	    if( rintf(tbb) != rintf(bad_data) ) {
	      n = n + (double)(1.0);
	      sum = sum + tbb;          
	      sumsq = sumsq + (tbb * tbb);
	    }
          
      }
    } 

    if(n > 1.0) {
      sqsum = sum * sum;
      num = (n * sumsq) - sqsum;
      den = n * (n - (double)(1.0));
      a = num / den;
      sig = sqrt(a);
      mn = sum / n;
    }
    else if(n == 1.0) {
      sig = 0.0;
      mn = sum;
    }
    else {
      sig = (double)bad_data;
      mn = (double)bad_data;
    } 
    sigma = (float)sig;
    rg_var.mean = (float)mn;
    
//  printf("std dev: %d %d %d %d %d %d \n", pxin.line_edge, pxin.elem_edge, imvl, imve, ne, nl);
//  printf("std dev: %f %f %f %f %f %f %f %f \n", n, sum, sumsq, sqsum, num, den, sig, mn);

    return sigma;
        
}
