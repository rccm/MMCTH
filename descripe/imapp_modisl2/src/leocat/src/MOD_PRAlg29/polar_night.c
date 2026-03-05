/*
!C ********************************************************************
!Description:

  Integer function polar_night.c
  Computes the clear sky confidence for nighttime polar (poleward of 60
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


int polar_night()


{

	
//  Declarations

    extern int polarnight_snow();
    extern int polarnight_ocean();
    extern int polarnight_land();

    int return_code;
    int ret_code;

//  Initializations

    return_code = 0;

//  Get initial clear sky confidence for current pixel.

    if(pxin.snow || pxin.ice) {

      ret_code = polarnight_snow();
      
    }
    
    else if(pxin.land) {    	
    
      ret_code = polarnight_land();

    }
    
    else {
    	
      ret_code = polarnight_ocean();
    	
    }
      
    return_code = ret_code;
    return return_code;

}
