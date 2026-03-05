/*$Id: co2standard_cloud_top_temp.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

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


void co2standard_cloud_top_temp (unsigned char verbose,
                                 imagerL1* imgr1, imagerL2* imgr2, 
			         sounderL1* sndr1, sounderL2* sndr2,
			         nwp_params nwp, rtm_profiles** rtm, 
			         rtm_toa rclr, rad_utils rutil)
{
  char *rout = {"co2standard_cloud_top_temp.c"};
  unsigned char  wv_corr = 1;
  const int nchan = 6, it[] = {5, 5, 4, 4, 3}, ib[] = {4, 3, 3, 2, 2}, chan_num[] = {31, 32, 33, 34, 35, 36};
  
/* Added by raf */
  const float FILLF = -999.0;
  const unsigned char FILLB = 99;
 
  const float PBOT = 700, PTOP = 100;
  float noise[nchan];
  double delL[nchan], lhs, rad[nchan], rad_clr[nchan];
  double xmin, *pp, *zz, *tt, a, b, double_cldp, sos, eca_sol;
  long i, ii;
  int i_nwp, ivza, iclr, ilev, isol, iprof, iprof_win, il, l, ret_val, imin=0, j;
  size_t nsol;
  typedef struct Rcld_Calculations {unsigned char flag; int itop, ibot, nprof, nprof_win; double **delR, **rhs, *Iwin;} rcld_params;
  rcld_params **rcld;
  
  nsol = sizeof(it)/sizeof(int);
  
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
  
  /*----------------------------------------------------------------------------
    Loop over all satellite pixels.
  ----------------------------------------------------------------------------*/
  
  for (i=0; i<imgr1->npts; i++) {
    
/* Added by raf */
	imgr2->cldz[i] = FILLF;
	imgr2->cldp[i] = FILLF;
	imgr2->cldt[i] = FILLF;
	imgr2->cldemiss[i] = FILLF;
	imgr2->cod_ir[i] = FILLF;
	imgr2->cod_vis[i] = FILLF;
	imgr2->cld_solm[i] = FILLB;
	  
//  if (imgr1->badmask[i]) continue;
//  if (imgr1->cldmask[i] == PROB_CLEAR || imgr1->cldmask[i] == CLEAR || imgr1->cldmask[i] == PROB_CLOUDY) continue;
	if (imgr1->cldmask[i] == PROB_CLEAR || imgr1->cldmask[i] == CLEAR) continue;
    
    i_nwp = imgr1->index_nwp[i];
    ivza = imgr1->index_vza[i];
    iclr = imgr1->index_rclr[i];
      
    /*rad_clr[0] = rtm[i_nwp][ivza].rad31_clr;
    rad_clr[1] = rtm[i_nwp][ivza].rad32_clr;
    rad_clr[2] = rtm[i_nwp][ivza].rad33_clr;
    rad_clr[3] = rtm[i_nwp][ivza].rad34_clr;
    rad_clr[4] = rtm[i_nwp][ivza].rad35_clr;
    rad_clr[5] = rtm[i_nwp][ivza].rad36_clr;*/
    
    rad_clr[0] = rclr.rad31_clr[iclr];
    rad_clr[1] = rclr.rad32_clr[iclr];
    rad_clr[2] = rclr.rad33_clr[iclr];
    rad_clr[3] = rclr.rad34_clr[iclr];
    rad_clr[4] = rclr.rad35_clr[iclr];
    rad_clr[5] = rclr.rad36_clr[iclr];
    
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
	  	    
	  rcld[i_nwp][ivza].delR[0][iprof] = rtm[i_nwp][ivza].cloud_prof31[ilev] - rad_clr[0];
	  rcld[i_nwp][ivza].delR[1][iprof] = rtm[i_nwp][ivza].cloud_prof32[ilev] - rad_clr[1];
	  rcld[i_nwp][ivza].delR[2][iprof] = rtm[i_nwp][ivza].cloud_prof33[ilev] - rad_clr[2];
	  rcld[i_nwp][ivza].delR[3][iprof] = rtm[i_nwp][ivza].cloud_prof34[ilev] - rad_clr[3];
	  rcld[i_nwp][ivza].delR[4][iprof] = rtm[i_nwp][ivza].cloud_prof35[ilev] - rad_clr[4];
	  rcld[i_nwp][ivza].delR[5][iprof] = rtm[i_nwp][ivza].cloud_prof36[ilev] - rad_clr[5];
	  
          for (isol=0; isol<nsol; isol++)
	    rcld[i_nwp][ivza].rhs[isol][iprof] = rcld[i_nwp][ivza].delR[it[isol]][iprof]/rcld[i_nwp][ivza].delR[ib[isol]][iprof];
	  
	  iprof++;
	}
	rcld[i_nwp][ivza].Iwin[iprof_win] = rtm[i_nwp][ivza].cloud_prof31[ilev];
	iprof_win++;
      }
      rcld[i_nwp][ivza].nprof = iprof++;
      rcld[i_nwp][ivza].nprof_win = iprof_win++;
      rcld[i_nwp][ivza].ibot = (imin + rcld[i_nwp][ivza].nprof_win) - 1;
      
      rcld[i_nwp][ivza].flag = YES;
    }
    
    l = 0;
    for (ilev=rcld[i_nwp][ivza].itop; ilev<=rcld[i_nwp][ivza].ibot; ilev++) {
      pp[l] = nwp.map[i_nwp].p_lev[ilev];
      zz[l] = nwp.map[i_nwp].z_lev[ilev];
      tt[l] = nwp.map[i_nwp].t_lev[ilev];
      l++;
    }
    
    if (wv_corr)
      ret_val = lin_interp_sorted(rcld[i_nwp][ivza].Iwin, rad[0], pp, rcld[i_nwp][ivza].nprof_win, &il, &a, &double_cldp);
    else
      ret_val = lin_interp_sorted(tt, imgr1->bt31[i], pp, rcld[i_nwp][ivza].nprof_win, &il, &a, &double_cldp);
    if (ret_val < 0) {
      imgr2->cldz[i] = (float)zz[0];
      imgr2->cldp[i] = (float)pp[0];
      imgr2->cldt[i] = (float)tt[0];
    }
    else
    if (ret_val > 0) {
      imgr2->cldz[i] = (float)zz[rcld[i_nwp][ivza].nprof_win-1];
      imgr2->cldp[i] = (float)pp[rcld[i_nwp][ivza].nprof_win-1];
      imgr2->cldt[i] = (float)tt[rcld[i_nwp][ivza].nprof_win-1];
    }
    else {
      imgr2->cldp[i] = (float)double_cldp;
      imgr2->cldz[i] = (float)(zz[il] + a * (zz[il+1] - zz[il]));
      imgr2->cldt[i] = (float)(tt[il] + a * (tt[il+1] - tt[il]));
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
	    imgr2->cld_solm[i] = isol + 1;
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
  }
  
  free(pp);
  free(zz);
  free(tt);
  
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

}
