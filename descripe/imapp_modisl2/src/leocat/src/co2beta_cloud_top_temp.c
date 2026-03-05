/*$Id: co2beta_cloud_top_temp.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

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


void co2beta_cloud_top_temp (unsigned char verbose,
                             imagerL1* imgr1, imagerL2* imgr2, 
			     sounderL1* sndr1, sounderL2* sndr2,
			     nwp_params nwp, rtm_profiles** rtm, 
			     rtm_toa rclr, rad_utils rutil)
{
  FILE
   *fp;
  char *rout = {"co2beta_cloud_top_temp.c"};
  char beta_file[MAX_STR_LEN];
  unsigned  wv_corr = 0, write_beta_prof = 0;
  const int nchan = 6, it[] = {5, 5, 4, 4, 3}, ib[] = {4, 3, 3, 2, 2}, chan_num[] = {31, 32, 33, 34, 35, 36};
  const float PBOT = 700, PTOP = 100;
  float noise[nchan], rad27, rad_clr27, rad_cld27;
  double emiss11, emiss12, trans_ac11, trans_ac12, rad_ac11, rad_ac12, rad_cld11, rad_cld12, min_emiss, max_emiss;
  double *trans_atm31, *trans_atm32, *rad_atm31, *rad_atm32, *emiss_prof31, *emiss_prof32, *beta_prof, *emiss_prof27;
  double delL[nchan], lhs, rad[nchan], rad_clr[nchan];
  double xmin, *pp, *zz, *tt, a, b, double_cldp, sos, eca_sol;
  long i, ii, current_pix, current_line;
  int i_nwp, ivza,  ilev, isol, iprof, iprof_win, il, l, ret_val, imin, j;
  size_t nsol;
  typedef struct Rcld_Calculations {unsigned char flag; int itop, ibot, nprof, nprof_win; double **delR, **rhs, *Iwin;} rcld_params;
  rcld_params **rcld;
  imin=0;
  nsol = sizeof(it)/sizeof(int);
  printf("%d\n",nsol);
  
  if (verbose == YES) fprintf(stdout,"%sIn co2standard_cloud_top_temp.c\n",EXE_PROMPT);
  
  system("date");
  
  for (i=0; i<nchan; i++) noise[i] = rutil.rad_v2w_ptr(chan_num[i], rutil.NEDR[chan_num[i]-1]);
  
  if ((rcld = (rcld_params **) malloc(sizeof(rcld_params *) * nwp.ncells)) == NULL)
    error_allo(rout,"rcld**");
  for (i=0; i<nwp.ncells; i++) {
    if ((rcld[i] = (rcld_params *) malloc(sizeof(rcld_params) * nwp.rtm_nvzen)) == NULL)
      error_allo(rout,"rcld*");
    for (ii=0; ii<nwp.rtm_nvzen; ii++) rcld[i][ii].flag = NO;
  }
  
  /*----------------------------------------------------------------------------
    Calculate the split window brightness temperature difference.
  ----------------------------------------------------------------------------*/
    
  if ((pp = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"pp");
    
  if ((zz = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"zz");
    
  if ((tt = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"tt");
    
  if ((trans_atm31 = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"trans_atm31");
    
  if ((trans_atm32 = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"trans_atm32");
    
  if ((rad_atm31 = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"rad_atm31");
    
  if ((rad_atm32 = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"rad_atm32");
    
  if ((emiss_prof31 = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"emiss_prof31");
    
  if ((emiss_prof32 = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"emiss_prof32");
    
  if ((emiss_prof27 = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"emiss_prof27");
    
  if ((beta_prof = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
    error_allo(rout,"beta_prof");
  
  /*----------------------------------------------------------------------------
    Loop over all satellite pixels.
  ----------------------------------------------------------------------------*/
  
  sprintf(beta_file,"../l2_output/%s.beta",imgr1->l1b_filename2);
  fp = fopen(beta_file, "w");
  
  for (i=0; i<imgr1->npts; i++) {
  
    /*if (imgr1->cldmask[i] == PROB_CLEAR || imgr1->cldmask[i] == CLEAR || imgr1->cldmask[i] == PROB_CLOUDY) continue;*/
    
    i_nwp = imgr1->index_nwp[i];
    ivza = imgr1->index_vza[i];
      
    rad_clr[0] = rtm[i_nwp][ivza].rad31_clr;
    rad_clr[1] = rtm[i_nwp][ivza].rad32_clr;
    rad_clr[2] = rtm[i_nwp][ivza].rad33_clr;
    rad_clr[3] = rtm[i_nwp][ivza].rad34_clr;
    rad_clr[4] = rtm[i_nwp][ivza].rad35_clr;
    rad_clr[5] = rtm[i_nwp][ivza].rad36_clr;
    
    rad[0] = rutil.planck_rad_fast_ptr(31, imgr1->bt31[i], rutil.T_planck, rutil.B_table);
    delL[0] = rad[0] - rad_clr[0];
    
    rad[1] = rutil.planck_rad_fast_ptr(32, imgr1->bt32[i], rutil.T_planck, rutil.B_table);
    delL[1] = rad[1] - rad_clr[1];
    
    rad[2] = imgr1->rad33[i];
    delL[2] = rad[2] - rad_clr[2];
    
    rad[3] = imgr1->rad34[i];
    delL[3] = rad[3] - rad_clr[3];
    
    rad[4] = imgr1->rad35[i];
    delL[4] = rad[4] - rad_clr[4];
    
    rad[5] = imgr1->rad36[i];
    delL[5] = rad[5] - rad_clr[5];
    
    rad27 = rutil.planck_rad_fast_ptr(29, imgr1->bt29[i], rutil.T_planck, rutil.B_table);
    rad_clr27 = rtm[i_nwp][ivza].rad29_clr;
    
    if (rcld[i_nwp][ivza].flag == NO) {
      
      if ((rcld[i_nwp][ivza].Iwin = (double *) malloc(nwp.nlevels*sizeof(double))) == NULL)
        error_allo(rout,"Iwin");
      rcld[i_nwp][ivza].delR = allocate_2d_double_ptr("delR", nchan, nwp.nlevels);
      rcld[i_nwp][ivza].rhs = allocate_2d_double_ptr("rhs", nchan, nwp.nlevels);
            
      iprof = 0;
      iprof_win = 0;
      xmin = DBL_MAX;
      for (ilev=0; ilev<=nwp.map[i_nwp].sfc_level-5; ilev++) {
        if (nwp.map[i_nwp].t_lev[ilev] < xmin && nwp.map[i_nwp].p_lev[ilev] >= PTOP) {
          xmin = nwp.map[i_nwp].t_lev[ilev];
          imin = ilev;
        }
      }
      rcld[i_nwp][ivza].itop = imin;
      
      for (ilev=imin; ilev<=nwp.map[i_nwp].sfc_level; ilev++) {
        
	if (nwp.map[i_nwp].p_lev[ilev] <= PBOT) {
	  	    
	  rcld[i_nwp][ivza].delR[0][iprof] = rtm[i_nwp][ivza].rad_atm31_clr[ilev] + 
            rutil.planck_rad_fast_ptr(31, nwp.map[i_nwp].t_lev[ilev], rutil.T_planck, rutil.B_table)*rtm[i_nwp][ivza].trans_atm31_clr[ilev] -
            rad_clr[0];
      
          rcld[i_nwp][ivza].delR[1][iprof] = rtm[i_nwp][ivza].rad_atm32_clr[ilev] + 
            rutil.planck_rad_fast_ptr(32, nwp.map[i_nwp].t_lev[ilev], rutil.T_planck, rutil.B_table)*rtm[i_nwp][ivza].trans_atm32_clr[ilev] -
            rad_clr[1];
      
          rcld[i_nwp][ivza].delR[2][iprof] = rtm[i_nwp][ivza].rad_atm33_clr[ilev] + 
            rutil.planck_rad_fast_ptr(33, nwp.map[i_nwp].t_lev[ilev], rutil.T_planck, rutil.B_table)*rtm[i_nwp][ivza].trans_atm33_clr[ilev] -
            rad_clr[2];
	
          rcld[i_nwp][ivza].delR[3][iprof] = rtm[i_nwp][ivza].rad_atm34_clr[ilev] + 
            rutil.planck_rad_fast_ptr(34, nwp.map[i_nwp].t_lev[ilev], rutil.T_planck, rutil.B_table)*rtm[i_nwp][ivza].trans_atm34_clr[ilev] -
            rad_clr[3];
	
          rcld[i_nwp][ivza].delR[4][iprof] = rtm[i_nwp][ivza].rad_atm35_clr[ilev] + 
            rutil.planck_rad_fast_ptr(35, nwp.map[i_nwp].t_lev[ilev], rutil.T_planck, rutil.B_table)*rtm[i_nwp][ivza].trans_atm35_clr[ilev] -
            rad_clr[4];
	
          rcld[i_nwp][ivza].delR[5][iprof] = rtm[i_nwp][ivza].rad_atm36_clr[ilev] + 
            rutil.planck_rad_fast_ptr(36, nwp.map[i_nwp].t_lev[ilev], rutil.T_planck, rutil.B_table)*rtm[i_nwp][ivza].trans_atm36_clr[ilev] -
            rad_clr[5];
	  
          for (isol=0; isol<nsol; isol++)
	    rcld[i_nwp][ivza].rhs[isol][iprof] = rcld[i_nwp][ivza].delR[it[isol]][iprof]/rcld[i_nwp][ivza].delR[ib[isol]][iprof];
	  
	  iprof++;
	}
	rcld[i_nwp][ivza].Iwin[iprof_win] = rtm[i_nwp][ivza].rad_atm31_clr[ilev] + 
            rutil.planck_rad_fast_ptr(31, nwp.map[i_nwp].t_lev[ilev], rutil.T_planck, rutil.B_table)*rtm[i_nwp][ivza].trans_atm31_clr[ilev];
	iprof_win++;
      }
      rcld[i_nwp][ivza].nprof = iprof++;
      rcld[i_nwp][ivza].nprof_win = iprof_win++;
      rcld[i_nwp][ivza].ibot = (imin + rcld[i_nwp][ivza].nprof_win) - 1;
      
      rcld[i_nwp][ivza].flag = YES;
    }
    
    for (l=0; l<nwp.nlevels; l++) {
      emiss_prof31[l] = MISSING_FLOAT;
      emiss_prof32[l] = MISSING_FLOAT;
      emiss_prof27[l] = MISSING_FLOAT;
      beta_prof[l] = MISSING_FLOAT;
    }
    l = 0;
    for (ilev=rcld[i_nwp][ivza].itop; ilev<=rcld[i_nwp][ivza].ibot; ilev++) {
      pp[l] = nwp.map[i_nwp].p_lev[ilev];
      zz[l] = nwp.map[i_nwp].z_lev[ilev];
      tt[l] = nwp.map[i_nwp].t_lev[ilev];
      trans_atm31[l] = rtm[i_nwp][ivza].trans_atm31_clr[ilev];
      trans_atm32[l] = rtm[i_nwp][ivza].trans_atm32_clr[ilev];
      rad_atm31[l] = rtm[i_nwp][ivza].rad_atm31_clr[ilev];
      rad_atm32[l] = rtm[i_nwp][ivza].rad_atm32_clr[ilev];
      
      if (write_beta_prof) {
	rad_cld27 = rutil.planck_rad_fast_ptr(29, tt[l], rutil.T_planck, rutil.B_table);
        emiss_prof27[ilev] = (rad27 - rad_clr27)/(rtm[i_nwp][ivza].rad_atm29_clr[ilev] + rtm[i_nwp][ivza].trans_atm29_clr[ilev]*rad_cld27 - rad_clr27);
      
        rad_cld11 = rutil.planck_rad_fast_ptr(31, tt[l], rutil.T_planck, rutil.B_table);
        emiss_prof31[ilev] = (delL[0])/(rad_atm31[l] + trans_atm31[l]*rad_cld11 - rad_clr[0]);
        rad_cld12 = rutil.planck_rad_fast_ptr(32, tt[l], rutil.T_planck, rutil.B_table);
        emiss_prof32[ilev] = (delL[1])/(rad_atm32[l] + trans_atm32[l]*rad_cld12 - rad_clr[1]);
            
        if (emiss_prof31[ilev] > 0.0 && emiss_prof31[ilev] < 1.0 && emiss_prof32[ilev] > 0.0 && emiss_prof32[ilev] < 1.0)
          beta_prof[ilev] = log(1.0 - emiss_prof32[ilev])/log(1.0-emiss_prof31[ilev]);
        else {
          beta_prof[ilev] = MISSING_FLOAT;
	  emiss_prof31[ilev] = MISSING_FLOAT;
	  emiss_prof32[ilev] = MISSING_FLOAT;
        }
      }
      l++;
    }
    
    array_minmax_double (nwp.nlevels, emiss_prof31, &min_emiss, &max_emiss);
    
    if (wv_corr)
      ret_val = lin_interp_sorted(rcld[i_nwp][ivza].Iwin, rad[0], pp, rcld[i_nwp][ivza].nprof_win, &il, &a, &double_cldp);
    else
      ret_val = lin_interp_sorted(tt, imgr1->bt31[i], pp, rcld[i_nwp][ivza].nprof_win, &il, &a, &double_cldp);
    if (ret_val < 0) {
      imgr2->cldz[i] = (float)zz[0];
      imgr2->cldp[i] = (float)pp[0];
      imgr2->cldt[i] = (float)tt[0];
      trans_ac11 = trans_atm31[0];
      trans_ac12 = trans_atm32[0];
      rad_ac11 = rad_atm31[0];
      rad_ac12 = rad_atm32[0];
    }
    else
    if (ret_val > 0) {
      imgr2->cldz[i] = (float)zz[rcld[i_nwp][ivza].nprof_win-1];
      imgr2->cldp[i] = (float)pp[rcld[i_nwp][ivza].nprof_win-1];
      imgr2->cldt[i] = (float)tt[rcld[i_nwp][ivza].nprof_win-1];
      trans_ac11 = trans_atm31[rcld[i_nwp][ivza].nprof_win-1];
      trans_ac12 = trans_atm32[rcld[i_nwp][ivza].nprof_win-1];
      rad_ac11 = rad_atm31[rcld[i_nwp][ivza].nprof_win-1];
      rad_ac12 = rad_atm32[rcld[i_nwp][ivza].nprof_win-1];
    }
    else {
      imgr2->cldp[i] = (float)double_cldp;
      imgr2->cldz[i] = (float)(zz[il] + a * (zz[il+1] - zz[il]));
      imgr2->cldt[i] = (float)(tt[il] + a * (tt[il+1] - tt[il]));
      trans_ac11 = (trans_atm31[il] + a * (trans_atm31[il+1] - trans_atm31[il]));
      trans_ac12 = (trans_atm32[il] + a * (trans_atm32[il+1] - trans_atm32[il]));
      rad_ac11 = (rad_atm31[il] + a * (rad_atm31[il+1] - rad_atm31[il]));
      rad_ac12 = (rad_atm32[il] + a * (rad_atm32[il+1] - rad_atm32[il]));
    }
    imgr2->cld_solm[i] = 0;
    imgr2->cldemiss[i] = 1.0;
    imgr2->cod_ir[i] = 10.0;
    imgr2->cod_vis[i] = 20.0;
    
    xmin = DBL_MAX;
    for (isol=0; isol<nsol; isol++){ 
      if ((delL[it[isol]] <= -noise[it[isol]]) &&
          (delL[ib[isol]] <= -noise[ib[isol]])) {
	
	lhs = delL[it[isol]]/delL[ib[isol]];
	
	
	if (lin_interp_sorted(rcld[i_nwp][ivza].rhs[isol], lhs, pp, rcld[i_nwp][ivza].nprof, &il, &a, &double_cldp) <= 0) {
	  b = rcld[i_nwp][ivza].Iwin[il] + a * (rcld[i_nwp][ivza].Iwin[il+1] - rcld[i_nwp][ivza].Iwin[il]);
	  
	  eca_sol = delL[0]/(b - rad_clr[0]);
	  if (eca_sol < 0.01 || eca_sol > 1.0) continue;
	  
	  sos = 0.0;
	  for (j=2; j<nchan; ++j) {
	    b = rcld[i_nwp][ivza].delR[j][il] + a * (rcld[i_nwp][ivza].delR[j][il+1] - rcld[i_nwp][ivza].delR[j][il]);
	    b = delL[j] - eca_sol * b;
	    sos += b*b;
	  }
	  if (sos < xmin) {
	    xmin = sos;
	    imgr2->cldp[i] = (float)double_cldp;
            imgr2->cldz[i] = (float)(zz[il] + a * (zz[il+1] - zz[il]));
            imgr2->cldt[i] = (float)(tt[il] + a * (tt[il+1] - tt[il]));
	    imgr2->cldemiss[i] = eca_sol;
	    imgr2->cld_solm[i] = isol;
	    /*imgr2->cldt[i] = (imgr1->bt31[i]+tt[0])/2.0; 
	    lin_interp_sorted(tt, imgr2->cldt[i], pp, rcld[i_nwp][ivza].nprof_win, &il, &a, &double_cldp);*/
	    trans_ac11 = (trans_atm31[il] + a * (trans_atm31[il+1] - trans_atm31[il]));
            trans_ac12 = (trans_atm32[il] + a * (trans_atm32[il+1] - trans_atm32[il]));
            rad_ac11 = (rad_atm31[il] + a * (rad_atm31[il+1] - rad_atm31[il]));
            rad_ac12 = (rad_atm32[il] + a * (rad_atm32[il+1] - rad_atm32[il]));
	    if (eca_sol < 1.0) {
	      imgr2->cod_ir[i] = (-1.0)*imgr1->cos_satzen[i]*log(1.0 - eca_sol);
	      imgr2->cod_vis[i] = 2.0*imgr2->cod_ir[i];
	    }
	    else {
	      imgr2->cod_ir[i] = 10.0;
	      imgr2->cod_vis[i] = 20.0;
	    }
	  }
	}
      }
    }
    rad_cld11 = rutil.planck_rad_fast_ptr(31, imgr2->cldt[i], rutil.T_planck, rutil.B_table);
    rad_cld12 = rutil.planck_rad_fast_ptr(32, imgr2->cldt[i], rutil.T_planck, rutil.B_table);
    emiss11 = (delL[0])/(rad_ac11 + trans_ac11*rad_cld11 - rad_clr[0]);
    emiss12 = (delL[1])/(rad_ac12 + trans_ac12*rad_cld12 - rad_clr[1]);
    
    /*if (1) {
      l = 0;
      printf("emiss11=%f, emiss12=%f, beta=%f, cldt=%f, cldemiss=%f, btd1112=%f, ref1=%f\n",
        emiss11,emiss12,log(1.0 - emiss12)/log(1.0-emiss11),imgr2->cldt[i],imgr2->cldemiss[i],imgr1->bt31[i]-imgr1->bt32[i],imgr1->ref1[i]);
      printf("lat=%f, lon=%f, bt31=%f\n",imgr1->lat[i],imgr1->lon[i],imgr1->bt31[i]);
      for (ilev=rcld[i_nwp][ivza].itop; ilev<=rcld[i_nwp][ivza].ibot; ilev++) {
	printf("l=%d, emiss11=%f, emiss12=%f, beta=%f, cldt=%f\n",l,emiss_prof31[ilev],emiss_prof32[ilev],beta_prof[ilev],tt[l]);
	l++;
      }
      exit(1);
    }*/
      
    imgr2->cldbeta[i] = 0.0;
    if (imgr1->cldmask[i] == CLOUDY || imgr1->cldmask[i] == PROB_CLOUDY || imgr1->aeromask[i] == MOSTLY_ASH) {
      if (emiss11 < 1.0 && emiss12 < 1.0) {
        imgr2->cldbeta[i] = log(1.0 - emiss12)/log(1.0-emiss11);
      }
      else {
        imgr2->cldbeta[i] = 1.0;
      }
    }
      
    /*imgr2->cld_solm[i] = CLOUDY;
    if (min_emiss == 9999.0 || max_emiss == -9999.0) {
      imgr2->cld_solm[i] = CLEAR;
    }
    else {
      if ((min_emiss < 0.025 && max_emiss < 0.2) ||
          (emiss_prof31[(rcld[i_nwp][ivza].itop+rcld[i_nwp][ivza].ibot)/2] < 0.10 && 
	   emiss_prof31[(rcld[i_nwp][ivza].itop+rcld[i_nwp][ivza].ibot)/2] > 0.0)) imgr2->cld_solm[i] = CLEAR;
    }*/
    
    current_pix = i % imgr1->ncol;
    current_line = i/imgr1->ncol;
    if ((current_line % 5) == 0 && (current_pix % 5) == 0 && write_beta_prof) {
      /*printf("Writing pixel # %d\n",i);*/
      fwrite(&i, sizeof(long), 1, fp);
      fwrite(&imgr1->lat[i], sizeof(float), 1, fp);
      fwrite(&imgr1->lon[i], sizeof(float), 1, fp);
      fwrite(&imgr1->satzen[i], sizeof(float), 1, fp);
      fwrite(&imgr1->cldtype[i], sizeof(unsigned char), 1, fp);
      fwrite(&imgr1->ref1[i], sizeof(float), 1, fp);
      fwrite(&imgr1->ref20[i], sizeof(float), 1, fp);
      fwrite(&imgr1->ref26[i], sizeof(float), 1, fp);
      fwrite(&imgr1->bt31[i], sizeof(float), 1, fp);
      fwrite(&imgr1->bt32[i], sizeof(float), 1, fp);
      fwrite(&imgr2->cldt[i], sizeof(float), 1, fp);
      fwrite(&imgr2->cldp[i], sizeof(float), 1, fp);
      fwrite(&imgr2->cldz[i], sizeof(float), 1, fp);
      fwrite(&imgr2->cldemiss[i], sizeof(float), 1, fp);
      fwrite(&imgr2->cldbeta[i], sizeof(float), 1, fp);
      fwrite(&imgr2->cld_solm[i], sizeof(unsigned char), 1, fp);
      fwrite(&rcld[i_nwp][ivza].itop, sizeof(int), 1, fp);
      fwrite(&rcld[i_nwp][ivza].ibot, sizeof(int), 1, fp);
      fwrite(nwp.map[i_nwp].t_lev, sizeof(float), nwp.nlevels, fp);
      fwrite(nwp.map[i_nwp].p_lev, sizeof(float), nwp.nlevels, fp);
      fwrite(nwp.map[i_nwp].z_lev, sizeof(float), nwp.nlevels, fp);
      fwrite(emiss_prof27, sizeof(double), nwp.nlevels, fp);
      fwrite(emiss_prof31, sizeof(double), nwp.nlevels, fp);
      fwrite(emiss_prof32, sizeof(double), nwp.nlevels, fp);
      fwrite(beta_prof, sizeof(double), nwp.nlevels, fp);
    }
  
  }
  
  fclose(fp);
  
  free(pp);
  free(zz);
  free(tt);
  free(trans_atm31);
  free(trans_atm32);
  free(rad_atm31);
  free(rad_atm32);
  free(emiss_prof27);
  free(emiss_prof31);
  free(emiss_prof32);
  free(beta_prof);
  
  for (i=0; i<imgr1->npts; i++) {
    i_nwp = imgr1->index_nwp[i];
    ivza = imgr1->index_vza[i];
    if (rcld[i_nwp][ivza].flag == YES) {
      destroy_2d_double_ptr(nchan, rcld[i_nwp][ivza].delR);
      destroy_2d_double_ptr(nchan, rcld[i_nwp][ivza].rhs);
      free(rcld[i_nwp][ivza].Iwin);
      rcld[i_nwp][ivza].flag = NO;  
    }
  }
  for (i=0; i<nwp.ncells; i++) free(rcld[i]);
  free(rcld);
  
  /*fp = fopen("../l2_output/ref20", "w");
  fwrite(&imgr1->npts, sizeof(long), 1, fp);
  fwrite(imgr1->ref20, sizeof(float), imgr1->npts, fp);
  fclose(fp);*/
  
  system("date");

}
