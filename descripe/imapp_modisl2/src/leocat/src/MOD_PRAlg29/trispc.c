/*
!C ********************************************************************
!Description:

  Float function trispc.c    

  Regression function relating the pw values to the 8-11um btdif.
  The difference is related to the amount of pw in the atmosphere
  due to the weak water vapor lines present in this spectral
  region.  The regressions were determined from actual HIRS data
  and provide thresholds that can be compared to measured 8-11
  differences.

!Input arguments:
  tdf1      measured 8-11 micron BTD

!Output arguments:
  none

  Output value returned through function call.

!Revision History:
  Original from MODIS Cloud Mask FORTRAN code developed at CIMSS,
  UW-Madison (trispc.f).

  Converted to C
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>


float trispc(float tdf1)


{

    float a1, a2, a3, a4;
    float tri_thr;
    float x;

/*  Regression coefficients */
    a1 =  2.7681;
    a2 = -3.7290;
    a3 =  1.0540;
    a4 = -0.1020;

    x = tdf1;
    tri_thr = a1 + x * (a2 + x * (a3 + x * a4));

    return tri_thr;

}
