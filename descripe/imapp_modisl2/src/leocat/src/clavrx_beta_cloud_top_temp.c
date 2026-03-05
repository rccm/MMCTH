/*$Id: clavrx_beta_cloud_top_temp.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

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


void clavrx_beta_cloud_top_temp (unsigned char verbose,
                                 imagerL1* imgr1, imagerL2* imgr2, 
			         sounderL1* sndr1, sounderL2* sndr2,
			         nwp_params nwp, rtm_profiles** rtm, 
			         rtm_toa rclr, rad_utils rutil)
{
  char *rout = {"clavrx_beta_cloud_top_temp.c"};
  unsigned char DEBUG = 0, two_layer_flg = 0;
  
  const size_t nrow = 2, ncol = 2;
  const int iter_max = 10;
  const double conv_crit = 1.0;
  const float z_low = 2000.0;
  int iter, ifail;
  double conv_test;
  double y[ncol], f[ncol];
  double x[nrow], x_ap[nrow];
  double **mdiff, **xdiff, **tmp_mat2, **tmp_mat3, **delta_x;
  double **K, **Dy, **tmp_mat1;
  double **Sa, **Sa_inv, **Sx, **Sx_inv, **A, **E;
  double **Sy, **Sy_inv, **mat_trans22, **mat_add22, **mat_add21;
  
  long i;
  int current_pix, current_line, i1, i2, j1, j2, i_nwp, ivza;
  int l, ilev, imin, nprof, index;
  float *btd1112, bt11_uni, btd1112_uni, mintmp, maxtmp, tau_ap, Tc_temp, xmin, PTOP = 100.0, Tmean;
  float tsfc_est, trans_ac_11, rad_ac_11, trans_ac_12, rad_ac_12,
         rad_clear_11,  rad_clear_12;
  float layer_lapse_rate, layer_lapse_rate_p, *tt, *pp, *zz, *rad_atm31, *rad_atm32, *trans_atm31, *trans_atm32;
  float emiss_11, emiss_12, demiss12_demiss11, dB_dTc_11, dB_dTc_12, Bc_11, Bc_12, rad_11, rad_12,
        dB_dT_11, dB_dT_12, sign1, sign2, iter_avg=0.0, delL_11, delL_12;
  
  imin=0;
  mdiff = calloc_2d_double_ptr(rout, 2, 1);
  xdiff = calloc_2d_double_ptr(rout, 2, 1);  
  tmp_mat2 = calloc_2d_double_ptr(rout, 2, 1);
  tmp_mat3 = calloc_2d_double_ptr(rout, 2, 1);
  delta_x = calloc_2d_double_ptr(rout, 2, 1);
  K = calloc_2d_double_ptr(rout, 2, 2);
  Dy = calloc_2d_double_ptr(rout, 2, 2);
  tmp_mat1 = calloc_2d_double_ptr(rout, 2, 2);
  Sa = calloc_2d_double_ptr(rout, 2, 2);
  Sa_inv = calloc_2d_double_ptr(rout, 2, 2);
  Sx = calloc_2d_double_ptr(rout, 2, 2);
  Sx_inv = calloc_2d_double_ptr(rout, 2, 2);
  A = calloc_2d_double_ptr(rout, 2, 2);
  E = calloc_2d_double_ptr(rout, 2, 2);
  Sy = calloc_2d_double_ptr(rout, 2, 2);
  Sy_inv = calloc_2d_double_ptr(rout, 2, 2);
  mat_trans22 = calloc_2d_double_ptr(rout, 2, 2);
  mat_add22 = calloc_2d_double_ptr(rout, 2, 2);
  mat_add21 = calloc_2d_double_ptr(rout, 2, 1);
  
  mat_identity(E, 2, 2);
  
  /*----------------------------------------------------------------------------
    Calculate the split window brightness temperature difference.
  ----------------------------------------------------------------------------*/
  
  if ((btd1112 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"btd1112");
    
  for (i=0; i<imgr1->npts; i++) btd1112[i] = imgr1->bt31[i] - imgr1->bt32[i];
  
  if ((tt = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"tt");
  if ((pp = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"pp");
  if ((zz = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"zz");
  if ((rad_atm31 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"rad_atm31");
  if ((rad_atm32 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"rad_atm32");
  if ((trans_atm31 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"trans_atm31");
  if ((trans_atm32 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"trans_atm32");
  
  /*----------------------------------------------------------------------------
    Loop over all satellite pixels.
  ----------------------------------------------------------------------------*/
  
  for (i=0; i<imgr1->npts; i++) {
  
    if (DEBUG) exit(1);
    DEBUG = 0;
    if (imgr1->cldtype[i] == 5 && imgr1->bt31[i] > 275.0 && imgr1->ref26[i] > 0.02) DEBUG = 0;
    
    imgr2->cldbeta[i] = 0.0;
    imgr2->cldbeta_high[i] = 0.0;
    imgr2->cldbeta_mid[i] = 0.0;
    imgr2->cldbeta_low[i] = 0.0;
    imgr2->emiss11_high[i] = 0.0;
    imgr2->emiss11_mid[i] = 0.0;
    imgr2->emiss11_low[i] = 0.0;
    
    /*if (imgr1->cldmask[i] == PROB_CLEAR || imgr1->cldmask[i] == CLEAR || imgr1->cldmask[i] == PROB_CLOUDY) continue;*/
    
    /*----------------------------------------------------------------------------
      Compute spatial uniformity of observations for assigning obs. uncertainty.
    ----------------------------------------------------------------------------*/
    
    current_pix = i % imgr1->ncol;
    current_line = i/imgr1->ncol;
    i1 = max(0,current_pix-1);
    i2 = min(imgr1->ncol-1,current_pix+1);
    j1 = max(0,current_line-1);
    j2 = min(imgr1->nrow-1,current_line+1);    
    array_minmax_sub_float(i1, i2, j1, j2, imgr1->ncol, btd1112, &mintmp, &maxtmp);
    btd1112_uni = maxtmp - mintmp;
    array_minmax_sub_float(i1, i2, j1, j2, imgr1->ncol, imgr1->bt31, &mintmp, &maxtmp);
    bt11_uni = maxtmp - mintmp;  
    if (maxtmp > 380.0) {
      bt11_uni = 0.2;
      btd1112_uni = 0.1;
    }    
    
    /*----------------------------------------------------------------------------
      For convenience assign nwp indice to static value for given pixel and
      set the temperature profile to a 1D array.
    ----------------------------------------------------------------------------*/
    
    i_nwp = imgr1->index_nwp[i];
    ivza = imgr1->index_vza[i];
    
    xmin = DBL_MAX;
    for (ilev=0; ilev<nwp.nlevels-5; ilev++) {
      if (nwp.map[i_nwp].t_lev[ilev] < xmin && nwp.map[i_nwp].p_lev[ilev] >= PTOP) {
        xmin = nwp.map[i_nwp].t_lev[ilev];
        imin = ilev;
      }
    }
    
    l = 0;
    nprof = 0;
    for (ilev=imin; ilev<=nwp.map[i_nwp].sfc_level; ilev++) {
      tt[l] = nwp.map[i_nwp].t_lev[ilev];
      pp[l] = nwp.map[i_nwp].p_lev[ilev];
      zz[l] = nwp.map[i_nwp].z_lev[ilev];
      rad_atm31[l] = rtm[i_nwp][ivza].rad_atm31_clr[ilev];
      rad_atm32[l] = rtm[i_nwp][ivza].rad_atm32_clr[ilev];
      trans_atm31[l] = rtm[i_nwp][ivza].trans_atm31_clr[ilev];
      trans_atm32[l] = rtm[i_nwp][ivza].trans_atm32_clr[ilev];
      l++;
    }
    nprof = l++;      
    
    /*----------------------------------------------------------------------------
      Assign values to measurement vector, y.
    ----------------------------------------------------------------------------*/
    
    y[0] = imgr1->bt31[i];
    y[1] = btd1112[i];
    
    /*----------------------------------------------------------------------------
      Parameterize forward model error which is based on scene uniformity.
    ----------------------------------------------------------------------------*/
    
    mat_zero(Sy, 2, 2);
    mat_zero(Sy_inv, 2, 2);
    
    Sy[0][0] = 2.0 + min(2.0,bt11_uni/2.0);
    Sy[1][1] = 0.5 + min(1.0,btd1112_uni/2.0);    
    mat_inv(Sy, Sy_inv, 2, 2); 
    
    /*----------------------------------------------------------------------------
      Measurement error is assumed to be 0.0?
    ----------------------------------------------------------------------------*/
    
    /*----------------------------------------------------------------------------
      Determine a priori conditions.
    ----------------------------------------------------------------------------*/
    
    /*Default values*/
    x_ap[0] = imgr1->bt31[i];
    tau_ap = 1.0;
    
    switch (imgr1->cldtype[i]) {
      case FOG_TYPE:
        x_ap[0] = imgr1->bt31[i];
	tau_ap = 1.2; /*ec=0.70*/
	Sa[0][0] = (10.0*10.0);
	Sa[1][1] = (0.5*0.5);
	imgr2->cldbeta[i] = 1.3;
        imgr2->cldbeta_high[i] = 0.9;
        imgr2->cldbeta_mid[i] = 1.1;
        imgr2->cldbeta_low[i] = 1.3;
        imgr2->emiss11_high[i] = -99.0;
        imgr2->emiss11_mid[i] = -99.0;
        imgr2->emiss11_low[i] = -99.0;
        break;
      case WATER_TYPE:
        x_ap[0] = imgr1->bt31[i];
	tau_ap = 2.3; /*ec=0.90*/
	Sa[0][0] = (10.0*10.0);
	Sa[1][1] = (0.2*0.2);
	imgr2->cldbeta[i] = 1.3;
        imgr2->cldbeta_high[i] = 0.9;
        imgr2->cldbeta_mid[i] = 1.1;
        imgr2->cldbeta_low[i] = 1.3;
        imgr2->emiss11_high[i] = -99.0;
        imgr2->emiss11_mid[i] = -99.0;
        imgr2->emiss11_low[i] = -99.0;
	break;
      case MIXED_TYPE:
        x_ap[0] = imgr1->bt31[i];
	tau_ap = 2.3; /*ec=0.90*/
	Sa[0][0] = (10.0*10.0);
	Sa[1][1] = (0.2*0.2);
	imgr2->cldbeta[i] = 1.3;
        imgr2->cldbeta_high[i] = 0.9;
        imgr2->cldbeta_mid[i] = 1.1;
        imgr2->cldbeta_low[i] = 1.3;
        imgr2->emiss11_high[i] = -99.0;
        imgr2->emiss11_mid[i] = -99.0;
        imgr2->emiss11_low[i] = -99.0;
	break;
      case TICE_TYPE:
	x_ap[0] = imgr1->bt31[i];
	tau_ap = 2.3; /*ec=0.90*/
	Sa[0][0] = (10.0*10.0);
	Sa[1][1] = (0.2*0.2);
	imgr2->cldbeta[i] = 1.08;
        imgr2->cldbeta_high[i] = 1.16;
        imgr2->cldbeta_mid[i] = 1.10;
        imgr2->cldbeta_low[i] = 1.08;
        imgr2->emiss11_high[i] = -99.0;
        imgr2->emiss11_mid[i] = -99.0;
        imgr2->emiss11_low[i] = -99.0;
	break;
      case CIRRUS_TYPE:
        x_ap[0] = nwp.map[i_nwp].t_tropo + 40.0;
	tau_ap = 0.7; /*ec=0.50*/
	Sa[0][0] = (25.0*25.0);
	Sa[1][1] = (0.5*0.5);
	imgr2->cldbeta[i] = 1.08;
        imgr2->cldbeta_high[i] = 1.00;
        imgr2->cldbeta_mid[i] = 1.16;
        imgr2->cldbeta_low[i] = 1.3;
        imgr2->emiss11_high[i] = -99.0;
        imgr2->emiss11_mid[i] = -99.0;
        imgr2->emiss11_low[i] = -99.0;
	break;
      case OVERLAP_TYPE:
        x_ap[0] = nwp.map[i_nwp].t_tropo + 40.0;
	tau_ap = 1.2; /*ec=0.70*/
	Sa[0][0] = (25.0*25.0);
	Sa[1][1] = (0.5*0.5);
	imgr2->cldbeta[i] = 1.08;
        imgr2->cldbeta_high[i] = 1.00;
        imgr2->cldbeta_mid[i] = 1.16;
        imgr2->cldbeta_low[i] = 1.3;
        imgr2->emiss11_high[i] = -99.0;
        imgr2->emiss11_mid[i] = -99.0;
        imgr2->emiss11_low[i] = -99.0;
        break;
      default:
        x_ap[0] = max(imgr1->bt31[i] - 10.0,200.0);
	tau_ap = 0.7; /*ec=0.50*/
	Sa[0][0] = (50.0*50.0);
	Sa[1][1] = (0.5*0.5);
	imgr2->cldbeta[i] = 1.3;
        imgr2->cldbeta_high[i] = 0.9;
        imgr2->cldbeta_mid[i] = 1.1;
        imgr2->cldbeta_low[i] = 1.3;
        imgr2->emiss11_high[i] = -99.0;
        imgr2->emiss11_mid[i] = -99.0;
        imgr2->emiss11_low[i] = -99.0;
	break;
    }
    
    /*----------------------------------------------------------------------------
      Compute first guess emissivity from optical depth guess.
    ----------------------------------------------------------------------------*/
    
    x_ap[1] = 1.0 - exp((-1.0)*tau_ap/imgr1->cos_satzen[i]);
    
    /*----------------------------------------------------------------------------
      Determine the inverse of the first guess weight matrix (Sa).
    ----------------------------------------------------------------------------*/
    
    mat_inv(Sa, Sa_inv, 2, 2);
    
    /*----------------------------------------------------------------------------
      If in DEBUG mode, print some information.
    ----------------------------------------------------------------------------*/
    
    if (DEBUG) {
      fprintf(stdout,"%s%s - new retrieval parameters\n",EXE_PROMPT,rout);
      fprintf(stdout,"%sy[0] = %f, y[1] = %f\n",EXE_PROMPT,y[0],y[1]);
      fprintf(stdout,"%scld_type = %d, bt11_uni = %f, btd1112_uni = %f\n",
        EXE_PROMPT,imgr1->cldtype[i],bt11_uni,btd1112_uni);
      fprintf(stdout,"%sx_ap[0] = %f, x_ap[1] = %f, tau_ap = %f, T_tropo = %f\n",
        EXE_PROMPT,x_ap[0],x_ap[1],tau_ap,nwp.map[i_nwp].t_tropo);
      
      fprintf(stdout,"%sirow=1, icol=1, Sy=%f\n",EXE_PROMPT,Sy[0][0]);
      fprintf(stdout,"%sirow=1, icol=2, Sy=%f\n",EXE_PROMPT,Sy[0][1]);
      fprintf(stdout,"%sirow=2, icol=1, Sy=%f\n",EXE_PROMPT,Sy[1][0]);
      fprintf(stdout,"%sirow=2, icol=2, Sy=%f\n",EXE_PROMPT,Sy[1][1]);
      
      fprintf(stdout,"%sirow=1, icol=1, Sy_inv=%f\n",EXE_PROMPT,Sy_inv[0][0]);
      fprintf(stdout,"%sirow=1, icol=2, Sy_inv=%f\n",EXE_PROMPT,Sy_inv[0][1]);
      fprintf(stdout,"%sirow=2, icol=1, Sy_inv=%f\n",EXE_PROMPT,Sy_inv[1][0]);
      fprintf(stdout,"%sirow=2, icol=2, Sy_inv=%f\n",EXE_PROMPT,Sy_inv[1][1]);
      
      fprintf(stdout,"%sirow=1, icol=1, Sa=%f\n",EXE_PROMPT,Sa[0][0]);
      fprintf(stdout,"%sirow=1, icol=2, Sa=%f\n",EXE_PROMPT,Sa[0][1]);
      fprintf(stdout,"%sirow=2, icol=1, Sa=%f\n",EXE_PROMPT,Sa[1][0]);
      fprintf(stdout,"%sirow=2, icol=2, Sa=%f\n",EXE_PROMPT,Sa[1][1]);
      
      fprintf(stdout,"%sirow=1, icol=1, Sa_inv=%f\n",EXE_PROMPT,Sa_inv[0][0]);
      fprintf(stdout,"%sirow=1, icol=2, Sa_inv=%f\n",EXE_PROMPT,Sa_inv[0][1]);
      fprintf(stdout,"%sirow=2, icol=1, Sa_inv=%f\n",EXE_PROMPT,Sa_inv[1][0]);
      fprintf(stdout,"%sirow=2, icol=2, Sa_inv=%f\n\n",EXE_PROMPT,Sa_inv[1][1]);
    }
    
    if (imgr1->cldtype[i] == OVERLAP_TYPE && two_layer_flg) {
      ilev = max(0,min(nprof-1,locate(zz,nprof,z_low)));
      tsfc_est = tt[ilev];
      rad_clear_11 = (rutil.planck_rad_fast_ptr(31, tsfc_est, rutil.T_planck, rutil.B_table))*trans_atm31[ilev] +
        rad_atm31[ilev];
      rad_clear_12 = (rutil.planck_rad_fast_ptr(32, tsfc_est, rutil.T_planck, rutil.B_table))*trans_atm32[ilev] +
        rad_atm32[ilev];
    }
    else {
    
      /*----------------------------------------------------------------------------
        Define some parameters used in constructing the forward model, f.
      ----------------------------------------------------------------------------*/
    
      tsfc_est = nwp.map[i_nwp].t_sfc;
    
      /*----------------------------------------------------------------------------
        Compute 11 and 12 micron clear sky TOA radiance including sfc emissivity,
        but not surface reflection.  The Planck radiance is linearly interpolated
        between two values in the planck radiance table already in memory.
      ----------------------------------------------------------------------------*/
    
      rad_clear_11 = rtm[i_nwp][ivza].rad31_clr;
      rad_clear_12 = rtm[i_nwp][ivza].rad32_clr;
    }
    
    /*----------------------------------------------------------------------------
      Parameterize 11 um/12 um emissivity ratio (beta).
    ----------------------------------------------------------------------------*/
    
    /*----------------------------------------------------------------------------
      Get beta for highest valid level.
    ----------------------------------------------------------------------------*/    
        
    for (ilev=0; ilev<nprof; ilev++) {
      Tmean = tt[ilev];
    
      Bc_11 = rutil.planck_rad_fast_ptr(31, Tmean, rutil.T_planck, rutil.B_table);
      Bc_12 = rutil.planck_rad_fast_ptr(32, Tmean, rutil.T_planck, rutil.B_table);
		
      trans_ac_11 = trans_atm31[ilev];
      rad_ac_11 = rad_atm31[ilev];
      trans_ac_12 = trans_atm32[ilev];
      rad_ac_12 = rad_atm32[ilev];
    
      rad_11 = rutil.planck_rad_fast_ptr(31, imgr1->bt31[i], rutil.T_planck, rutil.B_table);
      delL_11 = rad_11 - rad_clear_11;
      rad_12 = rutil.planck_rad_fast_ptr(32, imgr1->bt32[i], rutil.T_planck, rutil.B_table);
      delL_12 = rad_12 - rad_clear_12;
    
      emiss_11 = (delL_11)/(rad_ac_11 + trans_ac_11*Bc_11 - rad_clear_11);
      emiss_12 = (delL_12)/(rad_ac_12 + trans_ac_12*Bc_12 - rad_clear_12);
            
      if (emiss_11 > 0.0 && emiss_11 < 1.0 && emiss_12 > 0.0 && emiss_12 < 1.0) {
        imgr2->cldbeta_high[i] = log(1.0 - emiss_12)/log(1.0-emiss_11);
	imgr2->emiss11_high[i] = emiss_11;
        break;
      }
      else {
        continue;
      }
    }
    
    /*----------------------------------------------------------------------------
      Get beta for middle level.
    ----------------------------------------------------------------------------*/
    
    ilev = max(0,min(nprof-1,locate(tt,nprof,(imgr1->bt31[i]+tt[0])/2.0)));
    /*ilev = ilev/2;*/
    Tmean = tt[ilev];
    
    Bc_11 = rutil.planck_rad_fast_ptr(31, Tmean, rutil.T_planck, rutil.B_table);
    Bc_12 = rutil.planck_rad_fast_ptr(32, Tmean, rutil.T_planck, rutil.B_table);
    
    trans_ac_11 = trans_atm31[ilev];
    rad_ac_11 = rad_atm31[ilev];
    trans_ac_12 = trans_atm32[ilev];
    rad_ac_12 = rad_atm32[ilev];
    
    rad_11 = rutil.planck_rad_fast_ptr(31, imgr1->bt31[i], rutil.T_planck, rutil.B_table);
    delL_11 = rad_11 - rad_clear_11;
    rad_12 = rutil.planck_rad_fast_ptr(32, imgr1->bt32[i], rutil.T_planck, rutil.B_table);
    delL_12 = rad_12 - rad_clear_12;
    
    emiss_11 = (delL_11)/(rad_ac_11 + trans_ac_11*Bc_11 - rad_clear_11);
    emiss_12 = (delL_12)/(rad_ac_12 + trans_ac_12*Bc_12 - rad_clear_12);
            
    if (emiss_11 > 0.0 && emiss_11 < 1.0 && emiss_12 > 0.0 && emiss_12 < 1.0) {
      imgr2->cldbeta_mid[i] = log(1.0 - emiss_12)/log(1.0-emiss_11);
      imgr2->emiss11_mid[i] = emiss_11;
    }
    
    /*----------------------------------------------------------------------------
      Get beta for the lowest valid level.
    ----------------------------------------------------------------------------*/
    
    index = max(0,min(nprof-1,locate(tt,nprof,imgr1->bt31[i])));
    index += 0;
    
    for (ilev=index; ilev>=0; ilev--) {
      
      Tmean = tt[ilev];
    
      Bc_11 = rutil.planck_rad_fast_ptr(31, Tmean, rutil.T_planck, rutil.B_table);
      Bc_12 = rutil.planck_rad_fast_ptr(32, Tmean, rutil.T_planck, rutil.B_table);
    
      trans_ac_11 = trans_atm31[ilev];
      rad_ac_11 = rad_atm31[ilev];
      trans_ac_12 = trans_atm32[ilev];
      rad_ac_12 = rad_atm32[ilev];
    
      rad_11 = rutil.planck_rad_fast_ptr(31, imgr1->bt31[i], rutil.T_planck, rutil.B_table);
      delL_11 = rad_11 - rad_clear_11;
      rad_12 = rutil.planck_rad_fast_ptr(32, imgr1->bt32[i], rutil.T_planck, rutil.B_table);
      delL_12 = rad_12 - rad_clear_12;
    
      emiss_11 = (delL_11)/(rad_ac_11 + trans_ac_11*Bc_11 - rad_clear_11);
      emiss_12 = (delL_12)/(rad_ac_12 + trans_ac_12*Bc_12 - rad_clear_12);
            
      if (emiss_11 > 0.0 && emiss_11 < 1.0 && emiss_12 > 0.0 && emiss_12 < 1.0) {
        imgr2->cldbeta_low[i] = log(1.0 - emiss_12)/log(1.0-emiss_11);
        imgr2->emiss11_low[i] = emiss_11;
	break;
      }
      else {
        continue;
      }
    }
    
    /*----------------------------------------------------------------------------
      Determine the best beta to use in the cloud retrieval.
    ----------------------------------------------------------------------------*/
    
    switch (imgr1->cldtype[i]) {
      case FOG_TYPE:
        imgr2->cldbeta[i] = imgr2->cldbeta_low[i];
        break;
      case WATER_TYPE:
        imgr2->cldbeta[i] = imgr2->cldbeta_low[i];
	break;
      case MIXED_TYPE:
        imgr2->cldbeta[i] = imgr2->cldbeta_low[i];
	break;
      case TICE_TYPE:
        imgr2->cldbeta[i] = imgr2->cldbeta_mid[i];
	imgr2->cldbeta[i] = imgr2->cldbeta_mid[i];
	break;
      case CIRRUS_TYPE:
        imgr2->cldbeta[i] = (imgr2->cldbeta_mid[i] + imgr2->cldbeta_high[i])/2.0;
	imgr2->cldbeta[i] = imgr2->cldbeta_mid[i];
	break;
      case OVERLAP_TYPE:
        imgr2->cldbeta[i] = (imgr2->cldbeta_mid[i] + imgr2->cldbeta_high[i])/2.0;
	imgr2->cldbeta[i] = imgr2->cldbeta_mid[i];
        break;
      default:
        imgr2->cldbeta[i] = imgr2->cldbeta_mid[i];
	break;
    }
    
    if (imgr1->aeromask[i] == MOSTLY_ASH) x_ap[0] = imgr1->bt31[i]-10.0;
    
    if (DEBUG) {
      printf("beta = %f\n",imgr2->cldbeta[i]);
      printf("Tmean = %f\n",(imgr1->bt31[i]+tt[0])/2.0);
    }
    
    /*trans_ac_11 = trans_atm31[ilev] +
                   (Tmean - tt[ilev]) *
                   (trans_atm31[ilev+1] - trans_atm31[ilev]) /
                   (tt[ilev+1] - tt[ilev]);
		 
      trans_ac_12 = trans_atm32[ilev] +
                   (Tmean - tt[ilev]) *
                   (trans_atm32[ilev+1] - trans_atm32[ilev]) /
                   (tt[ilev+1] - tt[ilev]);
		 
      rad_ac_11 = rad_atm31[ilev] +
                  (Tmean - tt[ilev]) *
                  (rad_atm31[ilev+1] - rad_atm31[ilev]) /
                  (tt[ilev+1] - tt[ilev]);
		 
      rad_ac_12 = rad_atm32[ilev] +
                  (Tmean - tt[ilev]) *
                  (rad_atm32[ilev+1] - rad_atm32[ilev]) /
                  (tt[ilev+1] - tt[ilev]);*/
    
    /*----------------------------------------------------------------------------
      Start the retrieval loop.
    ----------------------------------------------------------------------------*/
    
    iter = 0;
    conv_test = FLT_MAX;
    ifail = 0;
    
    x[0] = x_ap[0];
    x[1] = x_ap[1];
    
    while(1) {
      iter++;
      if (DEBUG) {
        fprintf(stdout,"%sTesting for convergence for iter, %d, y[0]-f[0] = %f, y[1]-f[1] = %f, conv_test = %f, conv_crit = %f\n",
	 EXE_PROMPT,iter,y[0]-f[0],y[1]-f[1],conv_test,conv_crit); 
      }
      
      /*----------------------------------------------------------------------------
        Test for convergence.
      ----------------------------------------------------------------------------*/
      
      if (conv_test < conv_crit) {
        if (DEBUG) fprintf(stdout,"%sConverged - iter=%d, conv_test=%f, conv_crit=%f\n",EXE_PROMPT,iter,conv_test,conv_crit);
	ifail = 0;
	break;
      }
      
      /*----------------------------------------------------------------------------
        Exit if no convergence after iter_max iterations.
      ----------------------------------------------------------------------------*/
      
      if (iter > iter_max) {
        if (DEBUG) fprintf(stdout,"%sFailed converegence\n",EXE_PROMPT);
	ifail = 1;
	break;
      }
      
      Tc_temp = x[0];
      
      /*----------------------------------------------------------------------------
        Estimate cloud altitude based on current estimate of cloud temperature.
      ----------------------------------------------------------------------------*/
      
      ilev = max(0,min(nprof-1,locate(tt,nprof,Tc_temp)));

      layer_lapse_rate = (tt[ilev+1] - tt[ilev]) / 
        (zz[ilev+1] - zz[ilev]);
                 
      if (layer_lapse_rate != 0.0)
        imgr2->cldz[i] = zz[ilev] + 
          (Tc_temp - tt[ilev]) / layer_lapse_rate;
      else
        imgr2->cldz[i] = zz[ilev];
     
      /*----------------------------------------------------------------------------
        Estimate above cloud radiance and transmittance by interpolating
	precomputed radiance and transmittance in Z.
      ----------------------------------------------------------------------------*/

      trans_ac_11 = trans_atm31[ilev] +
        (imgr2->cldz[i] - zz[ilev]) *
        (trans_atm31[ilev+1] - trans_atm31[ilev]) /
        (zz[ilev+1] - zz[ilev]);

      rad_ac_11 = rad_atm31[ilev] +
        (imgr2->cldz[i] - zz[ilev]) *
        (rad_atm31[ilev+1] - rad_atm31[ilev]) /
        (zz[ilev+1] - zz[ilev]);

      trans_ac_12 = trans_atm32[ilev] +
        (imgr2->cldz[i] - zz[ilev]) *
        (trans_atm32[ilev+1] - trans_atm32[ilev]) /
        (zz[ilev+1] - zz[ilev]);

      rad_ac_12 = rad_atm32[ilev] +
        (imgr2->cldz[i] - zz[ilev]) *
        (rad_atm32[ilev+1] - rad_atm32[ilev]) /
        (zz[ilev+1] - zz[ilev]);       
      
      /*----------------------------------------------------------------------------
        Call forward model.
      ----------------------------------------------------------------------------*/
      
      emiss_11 = x[1];
      demiss12_demiss11 = imgr2->cldbeta[i] * pow((1.0-emiss_11),(imgr2->cldbeta[i]-1.0));
      emiss_12 = 1.0 - pow((1.0-emiss_11),imgr2->cldbeta[i]);

      /*----------------------------------------------------------------------------
        Compute Planck emission for current cloud temperature.
      ----------------------------------------------------------------------------*/
      
      Bc_11 = rutil.planck_rad_fast_index_ptr(31, Tc_temp, rutil.T_planck, rutil.B_table, &dB_dTc_11);
      Bc_12 = rutil.planck_rad_fast_index_ptr(32, Tc_temp, rutil.T_planck, rutil.B_table, &dB_dTc_12);
      
      /*----------------------------------------------------------------------------
        Compute forward model for 11 um channel (ignoring scattering).
	Equation below is equivalent to:
	rad = sfc(toa) + below_cld_atmos(toa) + cld(toa) + above_cld_atmos(toa)
      ----------------------------------------------------------------------------*/
      
      rad_11 = emiss_11*rad_ac_11 + trans_ac_11 * emiss_11 * Bc_11 +
        (1.0 - emiss_11) * rad_clear_11;
      
      f[0] = rutil.planck_bt_fast_index_ptr(31, rad_11, rutil.T_planck, rutil.B_table, &dB_dT_11);
      
      /*----------------------------------------------------------------------------
        Compute forward model for 12 um channel.
	Equation below is equivalent to:
	rad = sfc(toa) + below_cld_atmos(toa) + cld(toa) + above_cld_atmos(toa)
      ----------------------------------------------------------------------------*/
      
      rad_12 = emiss_12*rad_ac_12 + trans_ac_12 * emiss_12 * Bc_12 +
        (1.0 - emiss_12) * rad_clear_12;
     
      f[1] = f[0] - rutil.planck_bt_fast_index_ptr(32, rad_12, rutil.T_planck, rutil.B_table, &dB_dT_12);
      
      /*----------------------------------------------------------------------------
        Compute kernel matrix - NOTE: index convention is NOT opposite of Fortran 
	since the Fortran90 procedure matmul assumes a matrix A(nrow,ncol) which is
	counter-intuitive considering the normal column major fortran convention.
      ----------------------------------------------------------------------------*/
      
      K[0][0] = (trans_ac_11 * emiss_11 * dB_dTc_11) / dB_dT_11;
      K[0][1] = (rad_ac_11 + trans_ac_11*Bc_11 - rad_clear_11)/dB_dT_11;
      K[1][0] = K[0][0] - trans_ac_12 * emiss_12 * dB_dTc_12 / dB_dT_12;
      K[1][1] = K[0][1] - (rad_ac_12+trans_ac_12*Bc_12-rad_clear_12)*(demiss12_demiss11)/dB_dT_12;
      /*printf("dT_11/dT_cld=%f, dbdt_1112/dT_cld=%f\n",K[0][0],K[1][0]);
      printf("dT_11/demiss_cld=%f, dbdt_1112/demiss_cld=%f\n",K[0][1],K[1][1]);
      if (imgr1->cldtype[i] = OVERLAP_TYPE) exit(1);*/
      
      /*----------------------------------------------------------------------------
        Compute next step as in Rodgers, 1976 - Eq. 102.
      ----------------------------------------------------------------------------*/
      
      mat_zero(tmp_mat1,2, 2);
      mat_mul(Sy_inv, 2, 2, K, 2, 2, tmp_mat1);
      mat_trans(K, mat_trans22, 2, 2);
      mat_mul(mat_trans22, 2, 2, tmp_mat1, 2, 2, mat_add22);
      mat_add(mat_add22, Sa_inv, Sx_inv, 2, 2);
      
      /*----------------------------------------------------------------------------
        Continue calculation of Eq. 102.
      ----------------------------------------------------------------------------*/
      
      mat_inv(Sx_inv, Sx, 2, 2);
      
     /*----------------------------------------------------------------------------
        Determine pertubation in state vector.
      ----------------------------------------------------------------------------*/
      
      mat_zero(tmp_mat2, 2, 1);
      mat_zero(tmp_mat3, 2, 1);
      
      mdiff[0][0] = y[0] - f[0];
      mdiff[1][0] = y[1] - f[1];
      xdiff[0][0] = x_ap[0]-x[0];
      xdiff[1][0] = x_ap[1]-x[1];
      
      mat_mul(Sy_inv, 2, 2, mdiff, 2, 1, tmp_mat2);
      mat_trans(K, mat_trans22, 2, 2);
      mat_mul(mat_trans22, 2, 2, tmp_mat2, 2, 1, mat_add21);     
      mat_zero(tmp_mat2, 2, 1);
      mat_mul(Sa_inv, 2, 2, xdiff, 2, 1, tmp_mat2);
      mat_add(mat_add21, tmp_mat2, tmp_mat3, 2, 1);
      mat_mul(Sx, 2, 2, tmp_mat3, 2, 1, delta_x);
      
      /*----------------------------------------------------------------------------
        Control the step size.
      ----------------------------------------------------------------------------*/
      
      mat_zero(tmp_mat2, 2, 1);
      mat_mul(Sx_inv, 2, 2, delta_x, 2, 1, tmp_mat2);
      tmp_mat2[0][0] *= delta_x[0][0];
      tmp_mat2[1][0] *= delta_x[1][0];
      
      conv_test = fabs(tmp_mat2[0][0] + tmp_mat2[1][0]);
      
      sign1 = 1.0;
      sign2 = 1.0;
      if (delta_x[0][0] < 0.0) sign1 = -1.0;
      if (delta_x[1][0] < 0.0) sign2 = -1.0;
      
      if (fabs(delta_x[0][0]) > fabs(x[0])/2.0)
        delta_x[0][0] = fabs(x[0]/2.0)*sign1;
	
      if (fabs(delta_x[1][0]) > fabs(x[1])/2.0)
        delta_x[1][0] = fabs(x[1]/2.0)*sign2;
	
      if (fabs(delta_x[0][0]) > 20.0)
        delta_x[0][0] = 20.0*sign1;
	
      if (fabs(delta_x[1][0]) > 0.2)
        delta_x[1][0] = 0.2*sign2;
	
      /*----------------------------------------------------------------------------
        Update state or retrieval vector.
      ----------------------------------------------------------------------------*/
      
       x[0] += delta_x[0][0];
       x[1] += delta_x[1][0];
       
       /*----------------------------------------------------------------------------
         Constrain state vector to reasonable vales - may be questionable.
       ----------------------------------------------------------------------------*/
       
       x[0] = max(180.0,min(tsfc_est,x[0]));
       x[1] = max(0.0,min(x[1],1.0));
       
       if (DEBUG) {
         printf("x/2.0 = %f  %f\n",x[0]/2.0,x[1]/2.0);
	 printf("delta_x = %f  %f\n",delta_x[0][0],delta_x[1][0]);
	 printf("x = %f  %f\n",x[0],x[1]);
	 printf("conv_test = %f\n",conv_test);
       }
      
    }
    
    /*----------------------------------------------------------------------------
       End of retrieval loop.
     ----------------------------------------------------------------------------*/
     
    /*----------------------------------------------------------------------------
       Determine final retrieval and corresponding error estimate and quality.
     ----------------------------------------------------------------------------*/
     
     if (ifail == 0) {
       mat_zero(tmp_mat1, 2, 2);
       mat_trans(K, mat_trans22, 2, 2);
       mat_mul(mat_trans22, 2, 2, Sy_inv, 2, 2, tmp_mat1);
       mat_mul(Sx, 2, 2, tmp_mat1, 2, 2, Dy);
       mat_mul(Dy, 2, 2, K, 2, 2, A);
       
       /*if (imgr1->cldtype[i] == CIRRUS_TYPE && imgr1->bt31[i] < 275.0) {
         printf("\n%5.1f  %5.1f\n",A[0][0],A[0][1]);
         printf("%5.1f  %5.1f\n",A[1][0],A[1][1]);
         exit(1);
       }*/
             
       imgr2->cldt[i] = x[0];
       if (x[1] < 1.00) {
         imgr2->cod_ir[i] = (-1.0)*imgr1->cos_satzen[i]*log(1.0 - x[1]);
	 imgr2->cldemiss[i] = 1.0 - exp((-1.0)*imgr2->cod_ir[i]);
	 imgr2->cod_vis[i] = 2.0*imgr2->cod_ir[i];
	 imgr2->cld_solm[i] = 1;
       }
       else {
         imgr2->cod_ir[i] = 10.0;
	 imgr2->cldemiss[i] = 1.0;
	 imgr2->cod_vis[i] = 20.0;
	 imgr2->cld_solm[i] = 1;
       }
       
     }
     else {
       imgr2->cldt[i] = imgr1->bt31[i];
       imgr2->cod_ir[i] = 10.0;
       imgr2->cldemiss[i] = 1.0;
       imgr2->cod_vis[i] = 20.0;
       imgr2->cld_solm[i] = 0;
       A[0][0] = -1.0;
       A[1][0] = -1.0;
       A[0][1] = -1.0;
       A[1][1] = -1.0;
     }
     
     /*----------------------------------------------------------------------------
       Find corresponding pressure and height for retrieved cloud temperature.
     ----------------------------------------------------------------------------*/
     
    if (imgr2->cldt[i] > 0.0) {
      ilev = max(0,min(nprof-1,locate(tt,nprof,imgr2->cldt[i])));

      layer_lapse_rate = (tt[ilev+1] - tt[ilev]) / 
        (zz[ilev+1] - zz[ilev]);
      layer_lapse_rate_p = (tt[ilev+1] - tt[ilev]) / 
        (pp[ilev+1] - pp[ilev]);
                 
      if (layer_lapse_rate != 0.0) {
        imgr2->cldz[i] = zz[ilev] + 
          (imgr2->cldt[i] - tt[ilev]) / layer_lapse_rate;
	imgr2->cldp[i] = pp[ilev] + 
          (imgr2->cldt[i] - tt[ilev]) / layer_lapse_rate_p;
      }
      else {
        imgr2->cldz[i] = zz[ilev];
	imgr2->cldp[i] = pp[ilev];
      }
    }
     
    iter_avg += (float)(iter-1);
    
  }
  
  /*----------------------------------------------------------------------------
     End of loop over satellite pixels.
   ----------------------------------------------------------------------------*/
   
   iter_avg /= (float)(imgr1->npts);
   printf("Average number of iterations = %f\n",iter_avg);
  
  /*----------------------------------------------------------------------------
    Free up local memory variables.
  ----------------------------------------------------------------------------*/
  
  destroy_2d_double_ptr(2, mdiff);
  destroy_2d_double_ptr(2, xdiff);
  destroy_2d_double_ptr(2, tmp_mat2);
  destroy_2d_double_ptr(2, tmp_mat3);
  destroy_2d_double_ptr(2, delta_x);
  destroy_2d_double_ptr(2, K);
  destroy_2d_double_ptr(2, Dy);
  destroy_2d_double_ptr(2, tmp_mat1);
  destroy_2d_double_ptr(2, Sa);
  destroy_2d_double_ptr(2, Sa_inv);
  destroy_2d_double_ptr(2, Sx);
  destroy_2d_double_ptr(2, Sx_inv);
  destroy_2d_double_ptr(2, A);
  destroy_2d_double_ptr(2, E);
  destroy_2d_double_ptr(2, Sy);
  destroy_2d_double_ptr(2, Sy_inv);
  destroy_2d_double_ptr(2, mat_trans22);
  destroy_2d_double_ptr(2, mat_add22);
  destroy_2d_double_ptr(2, mat_add21);
  
  free(btd1112);
  free(tt);
  free(pp);
  free(zz);
  free(rad_atm31);
  free(rad_atm32);
  free(trans_atm31);
  free(trans_atm32);

}
