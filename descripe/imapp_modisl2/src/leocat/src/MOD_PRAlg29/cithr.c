/*
!C ********************************************************************
!Description:

  Float function cithr.c      
  Computes the bi-dimensional linear or quadratic interpolation
  (lagrange form) of a given value to a table of scan angle and 11 um
  brightness temperature dependent 11-12 micron brightness temperature
  differences (thin cirrus cloud test thresholds).  Current look-up
  table from J. Key (NOAA, CIMSS, UW-Madison).

!Input arguments:
  int key            linear (1) or quadratic (2) interpolation
  float sec_vza      secant of the viewing zenith angle
  float bt11         11 micron measured brightness temperature (K)

!Output arguments:
  none

  float tci_thr  thin cirrus cloud test threshold
                 (returned through function call)

!Revision History:

  Original from MODIS Cloud Mask FORTRAN code developed at CIMSS,
  UW-Madison (tview.f).

  Converted to C
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>

/**********************************************************************/


float cithr(int key, float sec_vza, float bt11)


/**********************************************************************/

{

/*  Declarations */

    int i,i0,i1,i2,ii;
    int j,j0,j1,j2,jj;
    int Iflag;

    float lt0,lt1,lt2,lu0,lu1,lu2;
    float p,p0,p1,p2;
    float t,t0,t1,t2;
    float tci_thr;
    float u,u0,u1,u2;

    float ttab[13] = {190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0,
                      260.0, 270.0, 280.0, 290.0, 300.0 ,310.0 };

    float utab[5] = { 2.00, 1.75, 1.50, 1.25, 1.00 };

    float tab[5][13] =
     {.559,.424,.286,.137,.123,.198,.333,.696,1.217,3.184,5.178,8.269,12.452,
      .542,.416,.294,.162,.156,.240,.366,.704,1.184,2.926,4.854,7.885,11.985,
      .520,.405,.305,.194,.199,.294,.409,.715,1.141,2.591,4.433,7.389,11.381,
      .491,.391,.319,.238,.257,.367,.467,.729,1.083,2.140,3.866,6.720,10.567,
      .450,.370,.340,.300,.340,.470,.550,.750,1.000,1.500,3.060,5.770, 9.410 };
  
      
/*  Initializations */

    lt0 = 0.0;
    lt1 = 0.0;
    lt2 = 0.0;
    lu0 = 0.0;
    lu1 = 0.0;
    lu2 = 0.0;
    p = 0.0;
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    t = 0.0;
    t0 = 0.0;
    t1 = 0.0;
    t2 = 0.0;
    u = 0.0;
    u1 = 0.0;
    u2 = 0.0;
    i0 = 0;
    i1 = 0;
    i2 = 0;
    j0 = 0;
    j1 = 0;
    j2 = 0;
    jj = 0;
    
/*  Bounds check. */

    u = sec_vza;
    t = bt11;
    if (u > utab[0]) u = utab[0];
    if (u < utab[4]) u = utab[4];
    if (t < ttab[0]) t = ttab[0];
    if (t > ttab[12]) t = ttab[12];

    Iflag = 0;
    i = 1;
    while(i <= 4) {

      ii = i;
      if (u >= utab[i] && Iflag == 0) {
        Iflag = 1;
        if (key == 1) {
          i0 = ii - 1;
          i1 = ii; }
        else {
          if (ii == 4) {
            i0 = ii - 2;
            i1 = ii - 1;
            i2 = ii; }
          else {
            i0 = ii - 1;
            i1 = ii;
            i2 = ii + 1;
          }
        }
      }

      i++;

    }

    Iflag = 0;
    j = 1;
    while(j <= 12) {

      jj = j;
      if (t <= ttab[j] && Iflag == 0) {
        Iflag = 1;
        if (key == 1) {
          j0 = jj - 1;
          j1 = jj; }
        else {
          if (jj == 12) {
            j0 = jj - 2;
            j1 = jj - 1;
            j2 = jj; }
          else {
            j0 = jj - 1;
            j1 = jj;
            j2 = jj + 1;
          }
        }
      }

      j++;

    }

/*  Branch on interpolation method. */

    if (key == 1) {

/*    linear scheme */
/*    designate index values */

      u0 = utab[i0];
      u1 = utab[i1];
      t0 = ttab[j0];
      t1 = ttab[j1];

/*    lagrange polynomials */
      lu0 = (u-u1) / (u0-u1);
      lu1 = (u-u0) / (u1-u0);
      lt0 = (t-t1) / (t0-t1);
      lt1 = (t-t0) / (t1-t0);

/*    interpolating polynomials for the first dimension */
      p0 = tab[i0][j0]*lu0 + tab[i1][j0]*lu1;
      p1 = tab[i0][j1]*lu0 + tab[i1][j1]*lu1;

/*    interpolating polynomial for second dimension */
      p = p0*lt0 + p1*lt1;

//    printf("lin cithr: %f %f %f %f %f %f %f\n", u0,u1,t0,t1,p0,p1,p);

    }

    else {

/*    quadratic scheme */
/*    designate index values */

      u0 = utab[i0];
      u1 = utab[i1];
      u2 = utab[i2];
      t0 = ttab[j0];
      t1 = ttab[j1];
      t2 = ttab[j2];
      
/*    lagrange polynomials */
      lu0 = (u-u1) * (u-u2) / (u0-u1) / (u0-u2);
      lu1 = (u-u0) * (u-u2) / (u1-u0) / (u1-u2);
      lu2 = (u-u0) * (u-u1) / (u2-u0) / (u2-u1);
      lt0 = (t-t1) * (t-t2) / (t0-t1) / (t0-t2);
      lt1 = (t-t0) * (t-t2) / (t1-t0) / (t1-t2);
      lt2 = (t-t0) * (t-t1) / (t2-t0) / (t2-t1);
      
/*    interpolating polynomials for the first dimension */
      p0 = tab[i0][j0] * lu0 + tab[i1][j0] * lu1 + tab[i2][j0] * lu2;
      p1 = tab[i0][j1] * lu0 + tab[i1][j1] * lu1 + tab[i2][j1] * lu2;
      p2 = tab[i0][j2] * lu0 + tab[i1][j2] * lu1 + tab[i2][j2] * lu2;
      
/*    interpolating polynomial for second dimension */
      p = p0*lt0 + p1*lt1 + p2*lt2;

/*    printf("quad cithr: %f %f %f %f\n", p0,p1,p2,p); */

    }

    tci_thr = p;

    return tci_thr;

}
