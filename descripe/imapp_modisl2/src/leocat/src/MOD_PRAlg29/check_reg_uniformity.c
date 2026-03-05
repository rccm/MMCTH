/*
!C ********************************************************************
!Description:

  Integer function check_reg_uniformity.c
  Evaluates surface uniformity.
  Called from create_cloud_mask.c

!Input arguments:
  none

  Inputs from mask_processing_constants.h

!Output arguments:
  none

  int uniform       value of 1 if uniform, otherwise 0
                    (returned through function call)

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include "mask_processing_constants.h"
#include "pixel.h"
#include "granule.h"

int check_reg_uniformity()


{

/*  Declarations */

    int loc_uniform;
    int nland;
    int nwater;
    int ncoast;
    int ntotal;
    int nl, ne;
    int imve, imvl;
    int pxsnow;
    
    long int j, ide;
    
    unsigned char nise_val;
    unsigned char lsf;

/*  Initializations */

    nwater = 0;
    nland = 0;
    ncoast = 0;
    ntotal = nlcntx * necntx;
    loc_uniform = 1;
    imve = (necntx-1) / 2;
    imvl = (nlcntx-1) / 2;

//  Check for edges of granule.    
    if(pxin.line_edge == 1 || pxin.elem_edge == 1) {
    	
      loc_uniform = 0;
      
    }
//  Check for snow/ice in current pixel.
/*  else if(pxin.snow || pxin.ice) {
    	  
      loc_uniform = 0;  

    }  */
    else {
    	
//    Check region (3x3 pixel area centered on current pixel) for scene uniformity.
    
      for ( nl=0; nl<nlcntx; ++nl) {
        for ( ne=0; ne<necntx; ++ne) {

    	  j = ( (pxin.scan + (nl - imvl)) * neles) + pxin.elem;
    	  ide = j + (ne - imve);
    	
//        Ecosystem consistency.
    	  if ( (pxin.eco_type - *(g_eco+ide)) != 0) loc_uniform = 0;
    	
//        Snow/ice consistency.
    	  nise_val = *(g_nise+ide);
    	  pxsnow = 0;
    	  if (nise_val == 103 || nise_val == 104) pxsnow = 1;
    	  if(pxsnow == 1) loc_uniform = 0;
//  	  printf("loc_uniform: %d %d %d \n", nise_val, pxsnow, loc_uniform);
    	  
//        Count land, water, and coast pixels in region for consistency calculation
//        below.
    	  lsf = *(g_lsf+ide);
    	  if(lsf == 1 || lsf == 4) {
            nland++;
    	  }
    	  else if(lsf == 2 || lsf == 3) {
    	    ncoast++;
    	  }
    	  else {
    	    nwater++;
    	  }
    	
        }
      }
                
//    Check for coastlines in current region.
      if( (nwater + ncoast) == ntotal) {
    	
        if(nwater != ntotal) {
          loc_uniform = 0;
//        Provide "double coastlines".
          pxin.coast = 1;
          pxin.land = 1;
          pxin.water = 0;
        }
    
      }
      else if( (nland + ncoast) == ntotal) {
     
        if(nland != ntotal) loc_uniform = 0;
    	
      }
      else {
    	
        loc_uniform = 0;
        pxin.coast = 1;
        pxin.land = 1;
        pxin.water = 0;
    	
      }
      
    }      
/*    
    printf("\nFinal uniformity: %d \n", loc_uniform);
    printf("edges, snow: %d %d %d %d \n", pxin.line_edge, pxin.elem_edge, pxin.snow, pxin.ice);
    printf("lwc: %d %d %d \n", nland, nwater, ncoast);
*/    
    return loc_uniform;

}
