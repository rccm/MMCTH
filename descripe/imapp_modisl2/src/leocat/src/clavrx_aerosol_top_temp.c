/*$Id: clavrx_aerosol_top_temp.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <hdf.h>
#include "common_leocat.h"
#include "nwp_leocat.h"
#include "rtm_leocat.h"
#include "radutils_leocat.h"
#include "imagerL1_leocat.h"
#include "imagerL2_leocat.h"
#include "imagerL3_leocat.h"
#include "sounderL1_leocat.h"
#include "sounderL2_leocat.h"
#include "sounderL3_leocat.h"


void clavrx_aerosol_top_temp (unsigned char verbose,
                              imagerL1* imgr1, imagerL2* imgr2, 
			      sounderL1* sndr1, sounderL2* sndr2,
			      nwp_params nwp, rtm_profiles** rtm, 
			      rtm_toa rclr, rad_utils rutil)
{
//  char *rout = {"clavrx_cloud_top_temp.c"};
//  unsigned char DEBUG = 0;
  
//  const size_t nrow = 2, ncol = 2;
//  const int iter_max = 10;
//  const float conv_crit = 1.0;
//  int iter, ifail, pp;
//  float conv_test;
//  float y[ncol], f[ncol];
//  float x[nrow], x_ap[nrow];
//  
//  long i;
//  int current_pix, current_line, i1, i2, j1, j2, i_nwp;
//  int lev_bnd, l, ilev;
//  float *btd1112, bt11_uni, btd1112_uni, mintmp, maxtmp, determinant, tau_ap, Tc_temp;
//  float tsfc_est, tpw_est, z_bnd, trans_ac_11, rad_ac_11, trans_ac_12, rad_ac_12,
//        rad_bnd_11, rad_clear_11, rad_bnd_12, rad_clear_12, sfc_emiss_11, sfc_emiss_12;
//  float T_planck[nplanck], B11[nplanck], B12[nplanck];
//  float layer_lapse_rate, layer_lapse_rate_p, *level_temperature, beta;
//  float emiss_11, emiss_12, demiss12_demiss11, xx, dB_dTc_11, dB_dTc_12, Bc_11, Bc_12, rad_11, rad_12,
//        dB_dT_11, dB_dT_12, sign1, sign2, iter_avg=0.0;
  
  if (verbose == YES) fprintf(stdout,"%sIn clavrx_cloud_top_temp.c\n",EXE_PROMPT);
  
  system("date");
  
 
  
  system("date");

}
