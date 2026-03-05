/*
!C ********************************************************************
!Description:

  Integer function water_night.c
  Computes the initial clear sky confidence for nighttime water pixels.
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
  R. Frey           07/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include "pixel.h"


int water_night()


{

/*  Declarations */

    extern int ocean_night();
    extern int night_snow(int lnd);
  

    int return_code;
    int ret_code;
    int lnd;

/*  Initializations */

    return_code = 0;

/*  Get initial clear sky confidence for current pixel.
    night_snow is called in the case of ice-covered water.
*/

    if(pxin.ice) {

      lnd = 0;
      ret_code = night_snow(lnd);

    }

    else {

      ret_code = ocean_night();
  
    }

    return_code = ret_code;
    return return_code;

}
