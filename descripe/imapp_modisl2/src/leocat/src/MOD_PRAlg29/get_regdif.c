/*
************************************************************************
!F77

 !DESCRIPTION: 
   Computes reflectance or brightness temperature differences
   between the (center) pixel of interest and the surrounding
   eight pixel values.

 !Input parameters:
   None

   Inputs stored in structure variable 'pxin' of type "pixel_in" defined
   in pixel.h

 !Output Parameters:
   None

   Outputs stored in structure variable 'rg_var' of type "regional_var"
   defined in pixel.h

 !Revision History:
   Converted to C     11/07        R. Frey

 !Team-unique Header:

 !References and Credits:
 See Cloud Mask ATBD

 !END
***********************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>
#include "mask_processing_constants.h"
#include "pixel.h"
#include "granule.h"


void get_regdif(int band)


{

/*  Declarations */

    int center_line;
    int center_pix;
    int imve;
    int imvl;
    int nl;
    int ne;
    int idx;

    long int j, ide;

    float center_val;
    float surr_val;
    
/*  Initializations */

//  printf("Executing get_regdif \n");
    
    imve = 1;
    imvl = 1;
    center_line = (nlcntx / 2) + 1;
    center_pix = (necntx / 2) + 1;

/**********************************************************************/

//  Get differences about pixel of interest for input band.
 
    center_val = pxin.rad[band-1];
//  printf("center_val %f \n", center_val);

    idx = 0;
    for ( nl=0; nl<nlcntx; ++nl) {
      for ( ne=0; ne<necntx; ++ne) {

        j = ( (pxin.scan + (nl - imvl)) * neles) + pxin.elem;
        ide = j + (ne - imve);
        surr_val = *( ( g_rad[band-1] ) + ide );

        if(nl != (center_line - 1) || ne != (center_pix - 1) ) {
          if(rintf(center_val) != rintf(bad_data) &&
             rintf(surr_val) != rintf(bad_data)) {
            
            rg_var.diff[idx] = surr_val - center_val;

          }
          else {
            rg_var.diff[idx] = bad_data;
          }
//        printf("%d %d %ld %ld %d %f %f \n", nl,ne,j,ide,idx,surr_val,rg_var.diff[idx]);
          
          idx++;

        }

      }
    }

    return;

}
