/*$Id: fast_planck_routines.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <hdf.h>
#include "common_leocat.h"
#include "radutils_leocat.h"

void load_fast_planck_modis(rad_utils *rutil)
{
  char *rout = {"load_fast_planck_modis.c"};
  int i;

  rutil->B_table = allocate_2d_float_ptr("rutil->B_table", 17, nplanck);
  if ((rutil->T_planck = (float *) malloc(nplanck*sizeof(float))) == NULL)
    error_allo(rout,"rutil->T_planck");
    
  for (i=0; i<nplanck; i++) {
    rutil->T_planck[i] = 180.0 + i;
    rutil->B_table[0][i] = rutil->planck_rad_ptr(20, rutil->T_planck[i]);
    rutil->B_table[1][i] = rutil->planck_rad_ptr(21, rutil->T_planck[i]);
    rutil->B_table[2][i] = rutil->planck_rad_ptr(22, rutil->T_planck[i]);
    rutil->B_table[3][i] = rutil->planck_rad_ptr(23, rutil->T_planck[i]);
    rutil->B_table[4][i] = rutil->planck_rad_ptr(24, rutil->T_planck[i]);
    rutil->B_table[5][i] = rutil->planck_rad_ptr(25, rutil->T_planck[i]);
    rutil->B_table[6][i] = 0.0;
    rutil->B_table[7][i] = rutil->planck_rad_ptr(27, rutil->T_planck[i]);
    rutil->B_table[8][i] = rutil->planck_rad_ptr(28, rutil->T_planck[i]);
    rutil->B_table[9][i] = rutil->planck_rad_ptr(29, rutil->T_planck[i]);
    rutil->B_table[10][i] = rutil->planck_rad_ptr(30, rutil->T_planck[i]);
    rutil->B_table[11][i] = rutil->planck_rad_ptr(31, rutil->T_planck[i]);
    rutil->B_table[12][i] = rutil->planck_rad_ptr(32, rutil->T_planck[i]);
    rutil->B_table[13][i] = rutil->planck_rad_ptr(33, rutil->T_planck[i]);
    rutil->B_table[14][i] = rutil->planck_rad_ptr(34, rutil->T_planck[i]);
    rutil->B_table[15][i] = rutil->planck_rad_ptr(35, rutil->T_planck[i]);
    rutil->B_table[16][i] = rutil->planck_rad_ptr(36, rutil->T_planck[i]);
  }
  
}

void destroy_fast_planck_modis(rad_utils *rutil)
{
//  char *rout = {"destroy_fast_planck_modis.c"};
  
  destroy_2d_float_ptr(17, rutil->B_table);
  free(rutil->T_planck);
}

float planck_rad_fast_modis(int ch, float temp, float *T_planck, float **B_table)
{
  int index, l;
  float rad;
  
  index = ch - 20;
  
  l = max(0,min(nplanck-2,locate(T_planck, nplanck, temp)));
  
  rad = B_table[index][l] + (temp-T_planck[l])*(B_table[index][l+1]-B_table[index][l])/
    (T_planck[l+1]-T_planck[l]);
  return (rad);
}

float planck_bt_fast_modis(int ch, float rad, float *T_planck, float **B_table)
{
  int index, l;
  float dB_dT, bt;
  
  index = ch - 20;
  
  l = max(0,min(nplanck-2,locate(B_table[index], nplanck, rad)));
  dB_dT = (B_table[index][l+1]-B_table[index][l])/(T_planck[l+1]-T_planck[l]);
  bt = T_planck[l] + (rad - B_table[index][l])/dB_dT;
  
  return (bt);
}

float planck_rad_fast_index_modis(int ch, float temp, float *T_planck, float **B_table, float *dB_dT)
{  
  int l, index;
  
  index = ch - 20;
  
  l = max(0,min(nplanck-2,locate(T_planck, nplanck, temp)));
  *dB_dT = (B_table[index][l+1]-B_table[index][l])/(T_planck[l+1]-T_planck[l]);
  return B_table[index][l] + (temp-T_planck[l])*dB_dT[0];
  
}

float planck_bt_fast_index_modis(int ch, float rad, float *T_planck, float **B_table, float *dB_dT)
{
  int l, index;
  
  index = ch - 20;
  
  l= max(0,min(nplanck-2,locate(B_table[index], nplanck, rad)));
  *dB_dT = (B_table[index][l+1]-B_table[index][l])/(T_planck[l+1]-T_planck[l]);
  return T_planck[l] + (rad - B_table[index][l])/(dB_dT[0]);

}
