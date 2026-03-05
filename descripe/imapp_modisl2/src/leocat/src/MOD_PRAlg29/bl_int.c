#include <stdlib.h>
#include <stdio.h>
#include <math.h> 

#include "sst.h"

int bilinearInterpSST(int, int, float, float, const SST *, float *);

/* ===== bilinearInterpSST ================================================== */
/*
*| Name:
*|    bilinearInterpSST - Bi-linear interpolation scheme used with Reynolds
*|                     blended SST data
*|
*| Interface:
*|    int
*|    bilinearInterpSST (int xindx, int yindx, float lat, float lon, const SST *sstGrid float returnValue)
*|
*| Input:
*|   i         Grid index in the x direction
*|   j         Grid index in the y direction
*|   lat       Pixel latitude value
*|   lon       Pixel longitude value
*|   sstGrid   The SST grid from the Reynold's SST file
*|
*|
*| Input and Output:
*|    none
*|
*| Output:
*|    t    Interpolated SST Value
*|
*| Return values:
*|    0        - Success.
*|   -1        - Unsuccessful
*|
*|
*/

int
  bilinearInterpSST (int i, int j, float lat, float lon, const SST *sstGrid, float * t)
{
  int npoints_x = 360;
  int npoints_y = 180;
  float dx=1.0;
  float dy=1.0;
  int doYInterp = 1;
  int i1, j1;
  float p, ip0, ip1, jp0, jp1, t00, t01, t0, t10, t11, t1, ldi, ldj;

  if ( j < 0 || j > (npoints_y -1))
    return(-1);
  if ( i < 0 || i > (npoints_x -1))
    return(-1);
  
  ldi = lon - i; 
  if (ldi >= 0.5) {
     i1 = i + 1;
     if ( i == (npoints_x - 1))
        i1 = 0;
  }
  else {
     i1 = i - 1;
     if ( i == 0)
        i1 =  npoints_x - 1;
  }
  
  ldj = lat - (j - 90.0);

  if (ldj >= 0.5) {
     j1 = j + 1;
     if ( j == (npoints_y - 1))
        doYInterp = 0;
  }
  else {
     j1 = j - 1;
     if ( j == 0)
        doYInterp = 0;
  }

  p = lon - (i1 + 0.5);
  if ( p < -1.0 ) 
     p = p + 360.0;
  if ( p > 1.0) 
     p = p - 360.0;
  ip0 = fabs( p / dx);
  ip1 = dx - ip0;

  t00 = sstGrid->values[j][i];
  t01 = sstGrid->values[j][i1];
  t0 = (ip0 * t00) + (ip1 * t01);

  if( doYInterp ) {
    p = lat - (j1 - 90.0 + 0.5);
    jp0 = fabs( p / dy );
    jp1 = dy - jp0;

    t10 = sstGrid->values[j1][i];
    t11 = sstGrid->values[j1][i1];
    t1 = (ip0 * t10) + (ip1 * t11);

    *t = (jp0 * t0) + (jp1 * t1);
  }
  else
    *t = t0;

  return (0);
}

