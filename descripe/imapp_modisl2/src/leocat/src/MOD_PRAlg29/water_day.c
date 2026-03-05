/*
!C ********************************************************************
!Description:

  Integer function water_day.c
  Computes the initial clear sky confidence for daytime water pixels.
  Called from perform_cloud_tests.c

!Input arguments:
  none

  Inputs from structure variable 'pxin' of type "pixel_in" defined in
  pixel.h.

!Output arguments:
  none

  int return_code   successful completion is zero, otherwise non-zero
                    (returned through function call)

  Output initial clear sky confidence and associated data stored in
  structure variable 'pxout' of type "pixel_out" defined in pixel.h.

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include "pixel.h"


int water_day()


{

/*  Declarations */

    extern int ocean_day();
    extern int day_snow();

    int return_code;
    int ret_code;

/*  Initializations */

    return_code = 0;

/*  Get initial clear sky confidence for current pixel.
    Day_snow is called in the case of ice-covered ocean.
*/

    if(pxin.ice) {

      ret_code = day_snow(); }

    else {

      ret_code = ocean_day();
  
    }

    return_code = ret_code;
    return return_code;

}
