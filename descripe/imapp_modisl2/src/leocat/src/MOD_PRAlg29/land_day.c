/*
!C ********************************************************************
!Description:

  Integer function land_day.c
  Computes the clear sky confidence for daytime land pixels.
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
  Original version taken from MODIS cloud mask FORTRAN code developed
  at CIMSS, UW-Madison.
  Converted to C               R. Frey           06/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include "pixel.h"


int land_day()


{

/*  Declarations */


    extern int day_snow();
    extern int landday();
    extern int landday_coast();
    extern int landday_desert();
    extern int landday_desert_c();

    int return_code;
    int ret_code;

/*  Initializations */

    return_code = 0;

/*  Get initial clear sky confidence for current pixel.
*/

    if(pxin.snow || pxin.ice) {

      ret_code = day_snow(); }

    else if(pxin.desert && pxin.coast) {

      ret_code = landday_desert_c(); }
      
    else if(pxin.coast) {

      ret_code = landday_coast(); }

    else if(pxin.desert) {

      ret_code = landday_desert(); }
  
    else {

      ret_code = landday();
  
    }

    return_code = ret_code;
    return return_code;

}
