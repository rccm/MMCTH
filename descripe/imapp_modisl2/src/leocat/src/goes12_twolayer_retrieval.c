/*$Id: goes12_twolayer_retrieval.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

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


void create_algorithm_pointers(int16, char *[], imagerL2 *);
void destroy_algorithm_pointers(int16, char *[], imagerL2 *);
void modis06_cloud_prop(unsigned char, imagerL1 *, imagerL2 *, sounderL1 *, sounderL2 *, nwp_params, rtm_profiles **, rtm_toa, rad_utils);

void goes12_twolayer_retrieval (unsigned char verbose,
                                imagerL1* imgr1, imagerL2* imgr2, 
			        sounderL1* sndr1, sounderL2* sndr2,
			        nwp_params nwp, rtm_profiles** rtm, 
			        rtm_toa rclr, rad_utils rutil)
{
  char *rout = {"goes12_twolayer_retrieval.c"};
  char *sds_name[MAX_ALGO_OUTPUT] = {"cloud_top_temperature", "cloud_top_pressure", "cloud_top_height",
  "cloud_emissivity", "cloud_optical_depth_ir", "cloud_optical_depth_vis", "cloud_top_method"};
  
  unsigned char DEBUG = 0;
  const size_t nrow = 3, ncol = 3;
  const int iter_max = 10, DX=20, DY=20;
  const double conv_crit = 1.0;
  const float z_low = 3000.0;
  int iter, ifail;
  double conv_test;
  double y[ncol], f[ncol];
  double x[nrow], x_ap[nrow];
  double **mdiff, **xdiff, **tmp_mat2, **tmp_mat3, **delta_x;
  double **K, **Dy, **tmp_mat1;
  double **Sa, **Sa_inv, **Sx, **Sx_inv, **A, **E;
  double **Sy, **Sy_inv, **mat_trans33, **mat_add33, **mat_add31;
  
  int inwp, ivza;
  long i;
  int current_pix, current_line, i1, i2, j1, j2, imin=0, nprof, l, ilev;
  int nlow, start_row, end_row, start_col, end_col, irow, icol, index;
  float *cod_total_vis,   *bt33;
  float *btd1113, bt11_uni, btd1113_uni, bt67_uni, mintmp, maxtmp, xmin, PTOP = 100.0;
  float tau_high_ap, tau_low, layer_lapse_rate, layer_lapse_rate_p, *tt, *pp, *zz, junk, cod_ir;
  float *rad_atm27, *rad_atm31, *rad_atm33, *trans_atm27, *trans_atm31, *trans_atm33;
  float tsfc_est, rad_clear_67, rad_clear_11, rad_clear_13, Tc_temp_high, Tc_temp_low;
  float trans_ahc_67, rad_ahc_67, trans_ahc_11, rad_ahc_11, trans_ahc_13, rad_ahc_13;
  float trans_alc_67, rad_alc_67, trans_alc_11, rad_alc_11, trans_alc_13, rad_alc_13;
  float emiss_11_high, emiss_13_high, emiss_67_high, beta1_high, beta2_high, demiss13_demiss11_high, demiss67_demiss11_high;
  float emiss_11_low, emiss_13_low, emiss_67_low, beta1_low, beta2_low, demiss13_demiss11_low, demiss67_demiss11_low;
  float Bc_67_high, Bc_11_high, Bc_13_high, dB_dTc_67_high, dB_dTc_11_high, dB_dTc_13_high;
  float Bc_67_low, Bc_11_low, Bc_13_low, dB_dTc_67_low, dB_dTc_11_low, dB_dTc_13_low;
  float dB_dT_11,dB_dT_13, dB_dT_67, rad_11, rad_13, rad_67;
  float sign1, sign2, sign3, *btd1167, avg_lowt;
  
  if (verbose == YES) fprintf(stdout,"%sIn goes12_twolayer_retrieval.c\n",EXE_PROMPT);
  
  system("date");
 
  imgr2->cod_vis = NULL;
  imgr2->cldreff = NULL;
  modis06_cloud_prop(1, imgr1, imgr2, sndr1, sndr2, nwp, rtm, rclr, rutil);
  
  if ((cod_total_vis = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
    error_allo(rout,"cod_total_vis");
  for (i=0; i<imgr1->npts; i++) cod_total_vis[i] = imgr2->cod_vis[i];
  free(imgr2->cod_vis);
  
  create_algorithm_pointers(7, sds_name, imgr2);
  co2singlech_cloud_top_temp(1, imgr1, imgr2, sndr1, sndr2, nwp, rtm, rclr, rutil);
  
  /*nbias = 0;
  bias = 0;
  for (i=0; i<imgr1->npts; i++) {
    inwp = imgr1->index_nwp[i];
    ivza = imgr1->index_vza[i];
    if (imgr1->cldmask[i] == CLEAR) {
      bt33 = rutil.planck_bt_fast_ptr(33, imgr1->rad33[i], rutil.T_planck, rutil.B_table);
      bias += (bt33 - rtm[inwp][ivza].bt33_clr);
      bias_abs += fabs(bt33 - rtm[inwp][ivza].bt33_clr);
      nbias++;
    }
  }
  
  bias /= (float)nbias;
  bias_abs /= (float)nbias;
  printf("bias=%f, bias_abs=%f\n",bias,bias_abs);*/
  
  mdiff = calloc_2d_double_ptr(rout, 3, 1);
  xdiff = calloc_2d_double_ptr(rout, 3, 1);  
  tmp_mat2 = calloc_2d_double_ptr(rout, 3, 1);
  tmp_mat3 = calloc_2d_double_ptr(rout, 3, 1);
  delta_x = calloc_2d_double_ptr(rout, 3, 1);
  K = calloc_2d_double_ptr(rout, 3, 3);
  Dy = calloc_2d_double_ptr(rout, 3, 3);
  tmp_mat1 = calloc_2d_double_ptr(rout, 3, 3);
  Sa = calloc_2d_double_ptr(rout, 3, 3);
  Sa_inv = calloc_2d_double_ptr(rout, 3, 3);
  Sx = calloc_2d_double_ptr(rout, 3, 3);
  Sx_inv = calloc_2d_double_ptr(rout, 3, 3);
  A = calloc_2d_double_ptr(rout, 3, 3);
  E = calloc_2d_double_ptr(rout, 3, 3);
  Sy = calloc_2d_double_ptr(rout, 3, 3);
  Sy_inv = calloc_2d_double_ptr(rout, 3, 3);
  mat_trans33 = calloc_2d_double_ptr(rout, 3, 3);
  mat_add33 = calloc_2d_double_ptr(rout, 3, 3);
  mat_add31 = calloc_2d_double_ptr(rout, 3, 1);
  
  if ((btd1113 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"btd1113");
  if ((btd1167 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"btd1167");
  if ((bt33 = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
    error_allo(rout,"bt33");
  if ((tt = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"tt");
  if ((pp = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"pp");
  if ((zz = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"zz");
  if ((rad_atm27 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"rad_atm27");
  if ((rad_atm31 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"rad_atm31");
  if ((rad_atm33 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"rad_atm33");
  if ((trans_atm27 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"trans_atm27");
  if ((trans_atm31 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"trans_atm31");
  if ((trans_atm33 = (float *) malloc(nwp.nlevels*sizeof(float))) == NULL)
    error_allo(rout,"trans_atm33");
    
  for (i=0; i<imgr1->npts; i++) {
    bt33[i] = rutil.planck_bt_fast_ptr(33, imgr1->rad33[i], rutil.T_planck, rutil.B_table);
    btd1113[i] = imgr1->bt31[i] - bt33[i];
    btd1167[i] = imgr1->bt31[i] - imgr1->bt27[i];
  }
  
  for (i=0; i<imgr1->npts; i++) {
    imgr2->cldp_high[i] = imgr2->cldp[i];
    imgr2->cldt_high[i] = imgr2->cldt[i];
    imgr2->cldz_high[i] = imgr2->cldz[i];
    imgr2->cldemiss_high[i] = imgr2->cldemiss[i];
    imgr2->cod_vis_high[i] = imgr2->cod_vis[i];
    imgr2->cldp_low[i] = MISSING_FLOAT;
    imgr2->cldt_low[i] = MISSING_FLOAT;
    imgr2->cldz_low[i] = MISSING_FLOAT;
    imgr2->cldemiss_low[i] = MISSING_FLOAT;
    imgr2->cod_vis_low[i] = MISSING_FLOAT;
    
    if (imgr1->cldtype[i] != OVERLAP_TYPE || imgr2->cldp[i] > 600.0 || imgr2->cod_vis[i] > 4.0) continue;
    
    /*----------------------------------------------------------------------------
      Compute spatial uniformity of observations for assigning obs. uncertainty.
    ----------------------------------------------------------------------------*/
    
    current_pix = i % imgr1->ncol;
    current_line = i/imgr1->ncol;
    i1 = max(0,current_pix-1);
    i2 = min(imgr1->ncol-1,current_pix+1);
    j1 = max(0,current_line-1);
    j2 = min(imgr1->nrow-1,current_line+1);    
    array_minmax_sub_float(i1, i2, j1, j2, imgr1->ncol, btd1167, &mintmp, &maxtmp);
    bt67_uni = maxtmp - mintmp;
    array_minmax_sub_float(i1, i2, j1, j2, imgr1->ncol, btd1113, &mintmp, &maxtmp);
    btd1113_uni = maxtmp - mintmp;
    array_minmax_sub_float(i1, i2, j1, j2, imgr1->ncol, imgr1->bt31, &mintmp, &maxtmp);
    bt11_uni = maxtmp - mintmp;    
    if (maxtmp > 380.0) {
      bt67_uni = 0.1;
      bt11_uni = 0.2;
      btd1113_uni = 0.1;
    }
    
    avg_lowt = 0.0;
    nlow = 0;
    start_row = max(0, current_line-DX);
    end_row = min(current_line+DX, imgr1->nrow);
    start_col = max(0, current_pix-DY);
    end_col = min(current_pix+DY, imgr1->ncol);
    for (irow=start_row; irow<=end_row; irow++) {
      for (icol=end_col; icol<=end_col; icol++) {
        index = (irow*imgr1->ncol) + icol;
	if ((imgr1->cldtype[index] == WATER_TYPE || imgr1->cldtype[index] == MIXED_TYPE) && imgr2->cldp[index] > 500.0) {
	  avg_lowt += imgr2->cldt[index];
	  nlow++;
	}
      }
    }
    
    if (nlow > 0)
      avg_lowt /= (float)nlow;
    
    /*----------------------------------------------------------------------------
      For convenience assign nwp indice to static value for given pixel and
      set the temperature profile to a 1D array.
    ----------------------------------------------------------------------------*/
    
    inwp = imgr1->index_nwp[i];
    ivza = imgr1->index_vza[i];
    
    xmin = DBL_MAX;
    for (ilev=0; ilev<nwp.nlevels-20; ilev++) {
      if (nwp.map[inwp].t_lev[ilev] < xmin && nwp.map[inwp].p_lev[ilev] >= PTOP) {
        xmin = nwp.map[inwp].t_lev[ilev];
        imin = ilev;
      }
    }
    
    l = 0;
    nprof = 0;
    for (ilev=imin; ilev<=nwp.map[inwp].sfc_level; ilev++) {
      tt[l] = nwp.map[inwp].t_lev[ilev];
      pp[l] = nwp.map[inwp].p_lev[ilev];
      zz[l] = nwp.map[inwp].z_lev[ilev];
      rad_atm27[l] = rtm[inwp][ivza].rad_atm27_clr[ilev];
      rad_atm31[l] = rtm[inwp][ivza].rad_atm31_clr[ilev];
      rad_atm33[l] = rtm[inwp][ivza].rad_atm33_clr[ilev];
      trans_atm27[l] = rtm[inwp][ivza].trans_atm27_clr[ilev];
      trans_atm31[l] = rtm[inwp][ivza].trans_atm31_clr[ilev];
      trans_atm33[l] = rtm[inwp][ivza].trans_atm33_clr[ilev];
      l++;
    }
    nprof = l++;
    
    /*----------------------------------------------------------------------------
      Assign values to measurement vector, y.
    ----------------------------------------------------------------------------*/
    
    y[0] = imgr1->bt31[i];
    y[1] = btd1113[i];
    y[2] = btd1167[i];
    
    /*----------------------------------------------------------------------------
      Parameterize forward model error which is based on scene uniformity.
    ----------------------------------------------------------------------------*/
    
    mat_zero(Sy, 3, 3);
    mat_zero(Sy_inv, 3, 3);
    
    Sy[0][0] = 2.0 + min(2.0,bt11_uni/2.0);
    Sy[1][1] = 0.5 + min(1.0,btd1113_uni/2.0);
    Sy[2][2] = 0.5 + min(1.0,bt67_uni/2.0);   
    mat3d_inv(Sy, Sy_inv);
    
    /*----------------------------------------------------------------------------
      Measurement error is assumed to be 0.0?
    ----------------------------------------------------------------------------*/
    
    /*----------------------------------------------------------------------------
      Determine a priori conditions.
    ----------------------------------------------------------------------------*/
    
    ilev = max(0,min(nprof-1,locate(zz,nprof,z_low)));
    tau_high_ap = min(imgr2->cod_vis[i],4.0);
    tau_low = cod_total_vis[i] - imgr2->cod_vis[i];
    
    x_ap[0] = imgr2->cldt[i];
    x_ap[1] = 1.0 - exp((-1.0)*tau_high_ap);
    
    if (nlow > 0)
      x_ap[2] = avg_lowt;
    else
      x_ap[2] = tt[ilev];
    
    Sa[0][0] = (1000.0*1000.0);
    
    if (tau_high_ap > 4.0)
      Sa[1][1] = (1000.0*1000.0);
    else 
      Sa[1][1] = (1000.0*1000.0);
    
    if (nlow > 0)
      Sa[2][2] = (0.0*0.0);
    else
      Sa[2][2] = (0.0*0.0);
    
    /*----------------------------------------------------------------------------
      Determine the inverse of the first guess weight matrix (Sa).
    ----------------------------------------------------------------------------*/
    
    mat3d_inv(Sa, Sa_inv);
    
    /*----------------------------------------------------------------------------
      Define some parameters used in constructing the forward model, f.
    ----------------------------------------------------------------------------*/
    
    tsfc_est = nwp.map[inwp].t_sfc;
    
    /*----------------------------------------------------------------------------
      Compute 11 and 12 micron clear sky TOA radiance including sfc emissivity,
      but not surface reflection.  The Planck radiance is linearly interpolated
      between two values in the planck radiance table already in memory.
    ----------------------------------------------------------------------------*/
    
    rad_clear_67 = rtm[inwp][ivza].rad27_clr;
    rad_clear_11 = rtm[inwp][ivza].rad31_clr;
    rad_clear_13 = rtm[inwp][ivza].rad33_clr;
    
    /*----------------------------------------------------------------------------
      Start the retrieval loop.
    ----------------------------------------------------------------------------*/
    
    iter = 0;
    conv_test = FLT_MAX;
    ifail = 0;
    
    x[0] = x_ap[0];
    x[1] = x_ap[1];
    x[2] = x_ap[2];
    
    while(1) {
      iter++;
      
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
      
      Tc_temp_high = x[0];
      Tc_temp_low = x[2];
      
      /*----------------------------------------------------------------------------
        Estimate cloud altitude based on current estimate of cloud temperature.
      ----------------------------------------------------------------------------*/
      
      /*----------------------------------------------------------------------------
        Top cloud layer.
      ----------------------------------------------------------------------------*/
      
      ilev = max(0,min(nprof-1,locate(tt,nprof,Tc_temp_high)));

      layer_lapse_rate = (tt[ilev+1] - tt[ilev]) / 
        (zz[ilev+1] - zz[ilev]);
                 
      if (layer_lapse_rate != 0.0)
        imgr2->cldz_high[i] = zz[ilev] + 
          (Tc_temp_high - tt[ilev]) / layer_lapse_rate;
      else
        imgr2->cldz_high[i] = zz[ilev];
	
      /*----------------------------------------------------------------------------
        Estimate above high cloud radiance and transmittance by interpolating
	precomputed radiance and transmittance in Z.
      ----------------------------------------------------------------------------*/
      
      trans_ahc_67 = trans_atm27[ilev] +
        (imgr2->cldz_high[i] - zz[ilev]) *
        (trans_atm27[ilev+1] - trans_atm27[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      rad_ahc_67 = rad_atm27[ilev] +
        (imgr2->cldz_high[i] - zz[ilev]) *
        (rad_atm27[ilev+1] - rad_atm27[ilev]) /
        (zz[ilev+1] - zz[ilev]);
      
      trans_ahc_11 = trans_atm31[ilev] +
        (imgr2->cldz_high[i] - zz[ilev]) *
        (trans_atm31[ilev+1] - trans_atm31[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      rad_ahc_11 = rad_atm31[ilev] +
        (imgr2->cldz_high[i] - zz[ilev]) *
        (rad_atm31[ilev+1] - rad_atm31[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      trans_ahc_13 = trans_atm33[ilev] +
        (imgr2->cldz_high[i] - zz[ilev]) *
        (trans_atm33[ilev+1] - trans_atm33[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      rad_ahc_13 = rad_atm33[ilev] +
        (imgr2->cldz_high[i] - zz[ilev]) *
        (rad_atm33[ilev+1] - rad_atm33[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      /*----------------------------------------------------------------------------
        Bottom cloud layer.
      ----------------------------------------------------------------------------*/
      
      ilev = max(0,min(nprof-1,locate(tt,nprof,Tc_temp_low)));

      layer_lapse_rate = (tt[ilev+1] - tt[ilev]) / 
        (zz[ilev+1] - zz[ilev]);
      layer_lapse_rate_p = (tt[ilev+1] - tt[ilev]) / 
        (pp[ilev+1] - pp[ilev]);
                 
      if (layer_lapse_rate != 0.0) {
        imgr2->cldz_low[i] = zz[ilev] + 
          (Tc_temp_low - tt[ilev]) / layer_lapse_rate;
	imgr2->cldp_low[i] = pp[ilev] + 
          (Tc_temp_low - tt[ilev]) / layer_lapse_rate_p;
      }
      else {
        imgr2->cldz_low[i] = zz[ilev];
	imgr2->cldp_low[i] = pp[ilev];
      }
	
      /*----------------------------------------------------------------------------
        Estimate above low cloud radiance and transmittance by interpolating
	precomputed radiance and transmittance in Z.
      ----------------------------------------------------------------------------*/
      
      trans_alc_67 = trans_atm27[ilev] +
        (imgr2->cldz_low[i] - zz[ilev]) *
        (trans_atm27[ilev+1] - trans_atm27[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      rad_alc_67 = rad_atm27[ilev] +
        (imgr2->cldz_low[i] - zz[ilev]) *
        (rad_atm27[ilev+1] - rad_atm27[ilev]) /
        (zz[ilev+1] - zz[ilev]);
      
      trans_alc_11 = trans_atm31[ilev] +
        (imgr2->cldz_low[i] - zz[ilev]) *
        (trans_atm31[ilev+1] - trans_atm31[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      rad_alc_11 = rad_atm31[ilev] +
        (imgr2->cldz_low[i] - zz[ilev]) *
        (rad_atm31[ilev+1] - rad_atm31[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      trans_alc_13 = trans_atm33[ilev] +
        (imgr2->cldz_low[i] - zz[ilev]) *
        (trans_atm33[ilev+1] - trans_atm33[ilev]) /
        (zz[ilev+1] - zz[ilev]);
	
      rad_alc_13 = rad_atm33[ilev] +
        (imgr2->cldz_low[i] - zz[ilev]) *
        (rad_atm33[ilev+1] - rad_atm33[ilev]) /
        (zz[ilev+1] - zz[ilev]);
      
      /*----------------------------------------------------------------------------
        Call forward model.
      ----------------------------------------------------------------------------*/
      
      /*based on hexagonal solid columns with reff=30 um*/
      if (Tc_temp_high < 253.0) {
        beta1_high = 1.05;
        beta2_high = 0.89;
      }
      else {
        beta1_high = 1.34;
        beta2_high = 1.14;
      }
      
      /*based on spheres with reff=10 um*/
      beta1_low = 1.34;
      beta2_low = 1.14;
      
      /*----------------------------------------------------------------------------
        Get emissivities for top cloud layer.
      ----------------------------------------------------------------------------*/
      
      emiss_11_high = x[1];
      
      demiss13_demiss11_high = beta1_high * pow((1.0-emiss_11_high),(beta1_high-1.0));
      emiss_13_high = 1.0 - pow((1.0-emiss_11_high),beta1_high);
      
      demiss67_demiss11_high = beta2_high * pow((1.0-emiss_11_high),(beta2_high-1.0));
      emiss_67_high = 1.0 - pow((1.0-emiss_11_high),beta2_high);
      
      /*----------------------------------------------------------------------------
        Get emissivities for bottom cloud layer.
      ----------------------------------------------------------------------------*/
      
      emiss_11_low = 1.0 - exp((-1.0)*tau_low);
      
      demiss13_demiss11_low = beta1_low * pow((1.0-emiss_11_low),(beta1_low-1.0));
      emiss_13_low = 1.0 - pow((1.0-emiss_11_low),beta1_low);
      
      demiss67_demiss11_low = beta2_low * pow((1.0-emiss_11_low),(beta2_low-1.0));
      emiss_67_low = 1.0 - pow((1.0-emiss_11_low),beta2_low);
      
      /*----------------------------------------------------------------------------
        Compute Planck emission for current high cloud temperature.
      ----------------------------------------------------------------------------*/
      
      Bc_67_high = rutil.planck_rad_fast_index_ptr(27, Tc_temp_high, rutil.T_planck, rutil.B_table, &dB_dTc_67_high);
      Bc_11_high = rutil.planck_rad_fast_index_ptr(31, Tc_temp_high, rutil.T_planck, rutil.B_table, &dB_dTc_11_high);
      Bc_13_high = rutil.planck_rad_fast_index_ptr(33, Tc_temp_high, rutil.T_planck, rutil.B_table, &dB_dTc_13_high);
      
      /*----------------------------------------------------------------------------
        Compute Planck emission for current low cloud temperature.
      ----------------------------------------------------------------------------*/
      
      Bc_67_low = rutil.planck_rad_fast_index_ptr(27, Tc_temp_low, rutil.T_planck, rutil.B_table, &dB_dTc_67_low);
      Bc_11_low = rutil.planck_rad_fast_index_ptr(31, Tc_temp_low, rutil.T_planck, rutil.B_table, &dB_dTc_11_low);
      Bc_13_low = rutil.planck_rad_fast_index_ptr(33, Tc_temp_low, rutil.T_planck, rutil.B_table, &dB_dTc_13_low);
      
      /*----------------------------------------------------------------------------
        Compute forward model for 11 um channel (ignoring scattering).
      ----------------------------------------------------------------------------*/
      
      rad_11 = (1.0-emiss_11_low)*(1.0-emiss_11_high)*(rad_clear_11 - rad_alc_11) +
               (1.0-emiss_11_high)*(emiss_11_low*Bc_11_low*trans_alc_11 + rad_alc_11) +
	       emiss_11_high*rad_ahc_11 + 
	       emiss_11_high*Bc_11_high*trans_ahc_11;
	       
      f[0] = rutil.planck_bt_fast_index_ptr(31, rad_11, rutil.T_planck, rutil.B_table, &dB_dT_11);
      
      /*----------------------------------------------------------------------------
        Compute forward model for 13 um channel (ignoring scattering).
      ----------------------------------------------------------------------------*/
      
      rad_13 = (1.0-emiss_13_low)*(1.0-emiss_13_high)*(rad_clear_13 - rad_alc_13) +
               (1.0-emiss_13_high)*(emiss_13_low*Bc_13_low*trans_alc_13 + rad_alc_13) +
	       emiss_13_high*rad_ahc_13 + 
	       emiss_13_high*Bc_13_high*trans_ahc_13;
	       
      f[1] = f[0] - rutil.planck_bt_fast_index_ptr(33, rad_13, rutil.T_planck, rutil.B_table, &dB_dT_13);
      
      /*----------------------------------------------------------------------------
        Compute forward model for 6.7 um channel (ignoring scattering).
      ----------------------------------------------------------------------------*/
      
      rad_67 = (1.0-emiss_67_low)*(1.0-emiss_67_high)*(rad_clear_67 - rad_alc_67) +
               (1.0-emiss_67_high)*(emiss_67_low*Bc_67_low*trans_alc_67 + rad_alc_67) +
	       emiss_67_high*rad_ahc_67 + 
	       emiss_67_high*Bc_67_high*trans_ahc_67;
	       
      f[2] = f[0] - rutil.planck_bt_fast_index_ptr(27, rad_67, rutil.T_planck, rutil.B_table, &dB_dT_67);
      
      if (imgr2->cldemiss[i] < 0.90 && DEBUG) {
      printf("\n11 um\n");
      printf("tau_total_vis=%f\n",cod_total_vis[i]);
      printf("tcld_low=%f, pcld_low=%f, cldtau_low=%f\n",Tc_temp_low,imgr2->cldp_low[i],tau_low);
      printf("tcld_high=%f, pcld_high=%f, cldtau_high=%f\n",Tc_temp_high,imgr2->cldp[i],imgr2->cod_vis[i]);
      printf("rclr=%f\n",rad_clear_11);
      printf("trans_alc=%f\n",trans_alc_11);
      printf("rad_alc=%f\n",rad_alc_11);
      printf("rad_low=%f, emiss_low=%f\n",Bc_11_low,emiss_11_low);
      printf("trans_ahc=%f\n",trans_ahc_11);
      printf("rad_ahc=%f\n",rad_ahc_11);
      printf("rad_high=%f, emiss_high=%f\n",Bc_11_high,emiss_11_high);
      printf("rad_calc=%f, bt_calc=%f\n",rad_11, f[0]);
      printf("rad_obs=%f, bt_obs=%f\n",rutil.planck_rad_fast_ptr(31, imgr1->bt31[i], rutil.T_planck, rutil.B_table),imgr1->bt31[i]);
      
      printf("\n13 um\n");
      printf("tau_total_vis=%f\n",cod_total_vis[i]);
      printf("tcld_low=%f, pcld_low=%f, cldtau_low=%f\n",Tc_temp_low,imgr2->cldp_low[i],tau_low);
      printf("tcld_high=%f, pcld_high=%f, cldtau_high=%f\n",Tc_temp_high,imgr2->cldp[i],imgr2->cod_vis[i]);
      printf("rclr=%f\n",rad_clear_13);
      printf("trans_alc=%f\n",trans_alc_13);
      printf("rad_alc=%f\n",rad_alc_13);
      printf("rad_low=%f, emiss_low=%f\n",Bc_13_low,emiss_13_low);
      printf("trans_ahc=%f\n",trans_ahc_13);
      printf("rad_ahc=%f\n",rad_ahc_13);
      printf("rad_high=%f, emiss_high=%f\n",Bc_13_high,emiss_13_high);
      printf("rad_calc=%f, bt_calc=%f, bt11-bt13_calc=%f\n",rad_13,
        rutil.planck_bt_fast_ptr(33, rad_13, rutil.T_planck, rutil.B_table),
	f[1]);
      printf("rad_obs=%f, bt_obs=%f, bt11-bt13_obs=%f\n",imgr1->rad33[i],bt33[i],imgr1->bt31[i]-bt33[i]);
      
      printf("\n6.7 um\n");
      printf("tau_total_vis=%f\n",cod_total_vis[i]);
      printf("tcld_low=%f, pcld_low=%f, cldtau_low=%f\n",Tc_temp_low,imgr2->cldp_low[i],tau_low);
      printf("tcld_high=%f, pcld_high=%f, cldtau_high=%f\n",Tc_temp_high,imgr2->cldp[i],imgr2->cod_vis[i]);
      printf("rclr=%f\n",rad_clear_67);
      printf("trans_alc=%f\n",trans_alc_67);
      printf("rad_alc=%f\n",rad_alc_67);
      printf("rad_low=%f, emiss_low=%f\n",Bc_67_low,emiss_67_low);
      printf("trans_ahc=%f\n",trans_ahc_67);
      printf("rad_ahc=%f\n",rad_ahc_67);
      printf("rad_high=%f, emiss_high=%f\n",Bc_67_high,emiss_67_high);
      printf("rad_calc=%f, bt_calc=%f, bt11-bt67_calc=%f\n",rad_67,
        rutil.planck_bt_fast_ptr(27, rad_67, rutil.T_planck, rutil.B_table),
        f[2]);
      printf("rad_obs=%f, bt_obs=%f, bt11-bt67_obs=%f\n",rutil.planck_rad_fast_ptr(27, imgr1->bt27[i], rutil.T_planck, rutil.B_table),
        imgr1->bt27[i],btd1167[i]);
      }
           
     /*----------------------------------------------------------------------------
        Compute kernel matrix - NOTE: index convention is NOT opposite of Fortran 
	since the Fortran90 procedure matmul assumes a matrix A(nrow,ncol) which is
	counter-intuitive considering the normal column major fortran convention.
      ----------------------------------------------------------------------------*/
      
      /*df/dT_hc*/
      K[0][0] = (trans_ahc_11 * emiss_11_high * dB_dTc_11_high) / dB_dT_11;
      K[1][0] = K[0][0] - ((trans_ahc_13 * emiss_13_high * dB_dTc_13_high) / dB_dT_13);
      K[2][0] = K[0][0] - ((trans_ahc_67 * emiss_67_high * dB_dTc_67_high) / dB_dT_67);
      
      /*df/demiss_hc*/
      K[0][1] = ((-1.0)*(1.0-emiss_11_low)*(rad_clear_11-rad_alc_11) - 
        emiss_11_low*Bc_11_low*trans_alc_11 - rad_alc_11 + 
	rad_ahc_11 + Bc_11_high*trans_ahc_11)/dB_dT_11;
      
      junk = ((-1.0)*(1.0-emiss_13_low)*(rad_clear_13-rad_alc_13) - 
        emiss_13_low*Bc_13_low*trans_alc_13 - rad_alc_13 + 
	rad_ahc_13 + Bc_13_high*trans_ahc_13)/dB_dT_13;
      K[1][1] = K[0][1] - (junk)*(demiss13_demiss11_high);
      
      junk = ((-1.0)*(1.0-emiss_67_low)*(rad_clear_67-rad_alc_67) - 
        emiss_67_low*Bc_67_low*trans_alc_67 - rad_alc_67 + 
	rad_ahc_67 + Bc_67_high*trans_ahc_67)/dB_dT_67;
      K[2][1] = K[0][1] - (junk)*(demiss67_demiss11_high);
      
      /*df/dT_lc*/
      K[0][2] = ((trans_alc_11 * emiss_11_low * dB_dTc_11_low)*(1.0 - emiss_11_high)) / dB_dT_11;
      K[1][2] = K[0][2] - (((trans_alc_13 * emiss_13_low * dB_dTc_13_low)*(1.0 - emiss_13_high)) / dB_dT_13);
      K[2][2] = K[0][2] - (((trans_alc_67 * emiss_67_low * dB_dTc_67_low)*(1.0 - emiss_67_high)) / dB_dT_67);
      
      if (imgr2->cldemiss[i] < 0.90 && DEBUG) {
      printf("\ndT_11/dT_hc=%f, dbtd_1113/dT_hc=%f, dbtd_1167/dT_hc=%f\n",K[0][0],K[1][0],K[2][0]);
      printf("dT_11/demiss_hc=%f, dbtd_1113/demiss_hc=%f, dbtd_1167/demiss_hc=%f\n",K[0][1],K[1][1],K[2][1]);
      printf("dT_11/dT_lc=%f, dbtd_1113/dT_lc=%f, dbtd_1167/dT_lc=%f\n",K[0][2],K[1][2],K[2][2]);
      }
      
      /*----------------------------------------------------------------------------
        Compute next step as in Rodgers, 1976 - Eq. 102.
      ----------------------------------------------------------------------------*/
      
      mat_zero(tmp_mat1,3, 3);
      mat_mul(Sy_inv, 3, 3, K, 3, 3, tmp_mat1);
      mat_trans(K, mat_trans33, 3, 3);
      mat_mul(mat_trans33, 3, 3, tmp_mat1, 3, 3, mat_add33);
      mat_add(mat_add33, Sa_inv, Sx_inv, 3, 3);
      
      /*----------------------------------------------------------------------------
        Continue calculation of Eq. 102.
      ----------------------------------------------------------------------------*/
      
      mat3d_inv(Sx_inv, Sx);
      
      /*----------------------------------------------------------------------------
        Determine pertubation in state vector.
      ----------------------------------------------------------------------------*/
      
      mat_zero(tmp_mat2, 3, 1);
      mat_zero(tmp_mat3, 3, 1);
      
      mdiff[0][0] = y[0] - f[0];
      mdiff[1][0] = y[1] - f[1];
      mdiff[2][0] = y[2] - f[2];
      xdiff[0][0] = x_ap[0]-x[0];
      xdiff[1][0] = x_ap[1]-x[1];
      xdiff[2][0] = x_ap[2]-x[2];
      
      mat_mul(Sy_inv, 3, 3, mdiff, 3, 1, tmp_mat2);
      mat_trans(K, mat_trans33, 3, 3);
      mat_mul(mat_trans33, 3, 3, tmp_mat2, 3, 1, mat_add31);     
      mat_zero(tmp_mat2, 3, 1);
      mat_mul(Sa_inv, 3, 3, xdiff, 3, 1, tmp_mat2);
      mat_add(mat_add31, tmp_mat2, tmp_mat3, 3, 1);
      mat_mul(Sx, 3, 3, tmp_mat3, 3, 1, delta_x);
      
      /*----------------------------------------------------------------------------
        Control the step size.
      ----------------------------------------------------------------------------*/
      
      mat_zero(tmp_mat2, 3, 1);
      mat_mul(Sx_inv, 3, 3, delta_x, 3, 1, tmp_mat2);
      tmp_mat2[0][0] *= delta_x[0][0];
      tmp_mat2[1][0] *= delta_x[1][0];
      tmp_mat2[2][0] *= delta_x[2][0];
      
      conv_test = fabs(tmp_mat2[0][0] + tmp_mat2[1][0] + tmp_mat2[2][0]);
      
      if (imgr2->cldemiss[i] < 0.90 && DEBUG) {
      printf("\ndelta_x[0]=%f, delta_x[1]=%f, delta_x[2]=%f\n",delta_x[0][0],delta_x[1][0],delta_x[2][0]);
      printf("x[0]=%f, x[1]=%f, x[2]=%f\n",x[0],x[1],x[2]);
      }
      
      sign1 = 1.0;
      sign2 = 1.0;
      sign3 = 1.0;
      if (delta_x[0][0] < 0.0) sign1 = -1.0;
      if (delta_x[1][0] < 0.0) sign2 = -1.0;
      if (delta_x[2][0] < 0.0) sign3 = -1.0;
      
      if (fabs(delta_x[0][0]) > fabs(x[0])/2.0)
        delta_x[0][0] = fabs(x[0]/2.0)*sign1;
	
      if (fabs(delta_x[1][0]) > fabs(x[1])/2.0)
        delta_x[1][0] = fabs(x[1]/2.0)*sign2;
	
      if (fabs(delta_x[2][0]) > fabs(x[2])/2.0)
        delta_x[2][0] = fabs(x[2]/2.0)*sign3;
	
      if (fabs(delta_x[0][0]) > 20.0)
        delta_x[0][0] = 20.0*sign1;
	
      if (fabs(delta_x[1][0]) > 0.2)
        delta_x[1][0] = 0.2*sign2;
	
      if (fabs(delta_x[2][0]) > 20.0)
        delta_x[2][0] = 20.0*sign3;
	
      /*----------------------------------------------------------------------------
        Update state or retrieval vector.
      ----------------------------------------------------------------------------*/
      
       x[0] += delta_x[0][0];
       x[1] += delta_x[1][0];
       x[2] += delta_x[2][0];
       
       /*----------------------------------------------------------------------------
         Constrain state vector to reasonable vales - may be questionable.
       ----------------------------------------------------------------------------*/
       
       x[0] = max(180.0,min(tt[nprof-1],x[0]));
       x[1] = max(0.0,min(x[1],0.99));
       x[2] = max(x[0]-5.0,min(tt[nprof-1],x[2]));
       
       tau_low = cod_total_vis[i] - 2.0*((-1.0)*imgr1->cos_satzen[i]*log(1.0 - x[1]));
      
    }
    
    /*----------------------------------------------------------------------------
       End of retrieval loop.
     ----------------------------------------------------------------------------*/
     
     
    /*----------------------------------------------------------------------------
       Determine final retrieval and corresponding error estimate and quality.
     ----------------------------------------------------------------------------*/
     
    if (ifail == 0) {
      mat_zero(tmp_mat1, 3, 3);
      mat_trans(K, mat_trans33, 3, 3);
      mat_mul(mat_trans33, 3, 3, Sy_inv, 3, 3, tmp_mat1);
      mat_mul(Sx, 3, 3, tmp_mat1, 3, 3, Dy);
      mat_mul(Dy, 3, 3, K, 3, 3, A);
             
      imgr2->cldt_high[i] = x[0];
      imgr2->cldt_low[i] = x[2];
      if (x[1] < 1.00) {
        cod_ir = (-1.0)*imgr1->cos_satzen[i]*log(1.0 - x[1]);
	imgr2->cldemiss_high[i] = 1.0 - exp((-1.0)*cod_ir);
	imgr2->cod_vis_high[i] = 2.0*cod_ir;
	imgr2->cod_vis_low[i] = tau_low;
	cod_ir = tau_low/2.0;
	imgr2->cldemiss_low[i] = 1.0 - exp((-1.0)*cod_ir);
      }
      else {
	imgr2->cldemiss_high[i] = imgr2->cldemiss[i];
	imgr2->cod_vis_high[i] = imgr2->cod_vis[i];
	imgr2->cldemiss_low[i] = MISSING_FLOAT;
	imgr2->cod_vis_low[i] = MISSING_FLOAT;
	imgr2->cldt_low[i] = MISSING_FLOAT;
	imgr2->cldp_low[i] = MISSING_FLOAT;
	imgr2->cldz_low[i] = MISSING_FLOAT;
      }
       
    }
    else {
      imgr2->cldt_high[i] = imgr2->cldt[i];
      imgr2->cldemiss_high[i] = imgr2->cldemiss[i];
      imgr2->cod_vis_high[i] = imgr2->cod_vis[i];
      if (nlow == 0) {
        imgr2->cldt_low[i] = MISSING_FLOAT;
        imgr2->cldemiss_low[i] = MISSING_FLOAT;
        imgr2->cod_vis_low[i] = MISSING_FLOAT;
        imgr2->cldp_low[i] = MISSING_FLOAT;
        imgr2->cldz_low[i] = MISSING_FLOAT;
      }
      else {
        imgr2->cldt_low[i] = avg_lowt;
        imgr2->cldemiss_low[i] = 1.0;
        imgr2->cod_vis_low[i] = cod_total_vis[i] - imgr2->cod_vis[i];
      }
      A[0][0] = -1.0;
      A[1][0] = -1.0;
      A[2][0] = -1.0;
      A[0][1] = -1.0;
      A[1][1] = -1.0;
      A[2][1] = -1.0;
      A[0][2] = -1.0;
      A[1][2] = -1.0;
      A[2][2] = -1.0;
    }
    
    if (imgr2->cldemiss[i] < 0.90 && DEBUG) {
      printf("\n%5.1f  %5.1f  %5.1f\n",A[0][0],A[0][1],A[0][2]);
      printf("%5.1f  %5.1f  %5.1f\n",A[1][0],A[1][1],A[1][2]);
      printf("%5.1f  %5.1f  %5.1f\n",A[2][0],A[2][1],A[2][2]);
    }
    if (cod_total_vis[i] > 50.0 && DEBUG) exit(1);
     
    /*----------------------------------------------------------------------------
      Find corresponding pressure and height for retrieved cloud temperature.
    ----------------------------------------------------------------------------*/
     
    if (imgr2->cldt_high[i] > 0.0) {
      ilev = max(0,min(nprof-1,locate(tt,nprof,imgr2->cldt_high[i])));

      layer_lapse_rate = (tt[ilev+1] - tt[ilev]) / 
        (zz[ilev+1] - zz[ilev]);
      layer_lapse_rate_p = (tt[ilev+1] - tt[ilev]) / 
        (pp[ilev+1] - pp[ilev]);
                 
      if (layer_lapse_rate != 0.0) {
        imgr2->cldz_high[i] = zz[ilev] + 
          (imgr2->cldt_high[i] - tt[ilev]) / layer_lapse_rate;
        imgr2->cldp_high[i] = pp[ilev] + 
          (imgr2->cldt_high[i] - tt[ilev]) / layer_lapse_rate_p;
      }
      else {
        imgr2->cldz_high[i] = zz[ilev];
        imgr2->cldp_high[i] = pp[ilev];
      }
    }
    
    if (imgr2->cldt_low[i] > 0.0) {
      ilev = max(0,min(nprof-1,locate(tt,nprof,imgr2->cldt_low[i])));

      layer_lapse_rate = (tt[ilev+1] - tt[ilev]) / 
        (zz[ilev+1] - zz[ilev]);
      layer_lapse_rate_p = (tt[ilev+1] - tt[ilev]) / 
        (pp[ilev+1] - pp[ilev]);
                 
      if (layer_lapse_rate != 0.0) {
        imgr2->cldz_low[i] = zz[ilev] + 
          (imgr2->cldt_low[i] - tt[ilev]) / layer_lapse_rate;
        imgr2->cldp_low[i] = pp[ilev] + 
          (imgr2->cldt_low[i] - tt[ilev]) / layer_lapse_rate_p;
      }
      else {
        imgr2->cldz_low[i] = zz[ilev];
        imgr2->cldp_low[i] = pp[ilev];
      }
    }
  }
  
  destroy_algorithm_pointers(7, sds_name, imgr2);
  free(cod_total_vis);
  free(btd1113);
  free(bt33);
  free(tt);
  free(pp);
  free(zz);
  free(rad_atm27);
  free(rad_atm31);
  free(rad_atm33);
  free(trans_atm27);
  free(trans_atm31);
  free(trans_atm33);
  
  destroy_2d_double_ptr(3, mdiff);
  destroy_2d_double_ptr(3, xdiff);
  destroy_2d_double_ptr(3, tmp_mat2);
  destroy_2d_double_ptr(3, tmp_mat3);
  destroy_2d_double_ptr(3, delta_x);
  destroy_2d_double_ptr(3, K);
  destroy_2d_double_ptr(3, Dy);
  destroy_2d_double_ptr(3, tmp_mat1);
  destroy_2d_double_ptr(3, Sa);
  destroy_2d_double_ptr(3, Sa_inv);
  destroy_2d_double_ptr(3, Sx);
  destroy_2d_double_ptr(3, Sx_inv);
  destroy_2d_double_ptr(3, A);
  destroy_2d_double_ptr(3, E);
  destroy_2d_double_ptr(3, Sy);
  destroy_2d_double_ptr(3, Sy_inv);
  destroy_2d_double_ptr(3, mat_trans33);
  destroy_2d_double_ptr(3, mat_add33);
  destroy_2d_double_ptr(3, mat_add31);
  
  system("date");
}
