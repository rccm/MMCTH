/*
!C ********************************************************************
!Description:

  Integer function polar_day.c
  Computes the clear sky confidence for daytime polar (poleward of 60
  degrees latitude) pixels.
  Called from perform_cloud_tests.c

!Input arguments:
  none

  Inputs from structure variable 'pxin' of type "pixel_in" defined in
  pixel.h.

!Output arguments:
  none

  int return_code   returns number of cloud tests performed
                    (returned through function call)

  Output initial clear sky confidence and associated data stored in
  structure variable 'pxout' of type "pixel_out" defined in pixel.h.

!Revision History:
  Original version taken from MODIS cloud mask FORTRAN code developed
  at CIMSS, UW-Madison.
  Converted to C               R. Frey           10/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

// Includes

#include <stdio.h>
#include "pixel.h"


int polar_day()


{

//  Declarations


    extern int Antarctic_day();
    extern int polarday_snow();
    extern int polarday_ocean();
    extern int polarday_land();
    extern int polarday_coast();
    extern int polarday_desert();
    extern int polarday_desert_c();

    int return_code;
    int ret_code;

//  Initializations

    return_code = 0;

//  Get initial clear sky confidence for current pixel.

    if(pxin.Antarctic && pxin.land) {
    	
      ret_code = Antarctic_day();
    }
    else if(pxin.snow || pxin.ice) {

      ret_code = polarday_snow();
    }
    else if(pxin.land) {    	
    
      if(pxin.desert && pxin.coast) {

        ret_code = polarday_desert_c();
      }
      else if(pxin.coast) {

        ret_code = polarday_coast();
      }
      else if(pxin.desert) {

        ret_code = polarday_desert();
      }
      else {

        ret_code = polarday_land();
  
      }

    }
    else {
    	
      ret_code = polarday_ocean();
    	
    }
      
    return_code = ret_code;
    return return_code;

}
