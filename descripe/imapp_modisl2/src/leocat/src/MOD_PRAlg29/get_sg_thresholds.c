/*
!C ********************************************************************
!Description:

  Function get_sg_thresholds.c
  Computes 0.86 micron cloud test thresholds for water surfaces 
    characterized by sun-glint (reflectance angle <= 36 degrees).

!Input arguments:
  float refang            reflectance or "sun-glint" angle.

!Output arguments:
  float * dosgref2        pointer to output threshold array

!Revision History:
  Original from MODIS cloud mask FORTRAN code developed by CIMSS,
    UW-Madison.
  Converted to C.
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include "thresholds.h"

/**********************************************************************/


void get_sg_thresholds(float refang, float* dosgref2)


/**********************************************************************/

{

/*  Declarations */

    float a;
    float conf_range;
    float hi_ang;
    float hi_ang_val;
    float lo_ang;
    float lo_ang_val;
    float midpt;
    float power;
    
/**********************************************************************/

/*  Initializations */

    a = 0.0;

/**********************************************************************/

/*  Compute clear sky confidence mid-point and boundaries. */

    if(refang > snglnt_bounds[3]) {

      /*printf("Problem in get_sg_thresholds.c\n");*/
      /*printf("Reflectance angle > %f\n", snglnt_bounds[3]);*/
      dosgref2[2] = doref2[2];
      dosgref2[0] = doref2[0];
      dosgref2[1] = doref2[1];
      dosgref2[3] = doref2[3]; }

    else if(refang <= snglnt_bounds[1]) {

      dosgref2[2] = snglnt0[2];
      dosgref2[0] = snglnt0[0];
      dosgref2[1] = snglnt0[1];
      dosgref2[3] = snglnt0[3]; }

    else {

      if(refang > snglnt_bounds[1] && refang <= snglnt_bounds[2]) {

        lo_ang = snglnt_bounds[1];
        hi_ang = snglnt_bounds[2];
        lo_ang_val = snglnt10[0];
        hi_ang_val = snglnt10[1];
        power = snglnt10[3];
        conf_range = snglnt10[2]; }

      else if(refang > snglnt_bounds[2] && refang <= snglnt_bounds[3]) {

        lo_ang = snglnt_bounds[2];
        hi_ang = snglnt_bounds[3];
        lo_ang_val = snglnt20[0];
        hi_ang_val = snglnt20[1];
        power = snglnt20[3];
        conf_range = snglnt20[2];

      }

      a = (refang - lo_ang) / (hi_ang - lo_ang);
      midpt = lo_ang_val + a*(hi_ang_val - lo_ang_val);
      dosgref2[1] = midpt;
      dosgref2[2] = midpt - conf_range;
      dosgref2[0] = midpt + conf_range;
      dosgref2[3] = power;

    }

/*  printf("\n%f %f %f %f %f\n", snglnt_bounds[0],snglnt_bounds[1],
          snglnt_bounds[2],snglnt_bounds[3],refang);
    printf("%f %f %f %f %f\n", a, dosgref2[0],dosgref2[1],dosgref2[2],
          dosgref2[3]);
*/

    return;

}
