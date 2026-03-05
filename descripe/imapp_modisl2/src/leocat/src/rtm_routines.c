/*$Id: rtm_routines.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "nwp_leocat.h"
#include "rtm_leocat.h"
#include "radutils_leocat.h"

void tran_modisd101_(char *,int *,int*,float *,float *,float *,float *,int *,int *,int *,float *,int *);

/*----------------------------------------------------------------------------
  Call the CRTM for MODIS using this routine.
----------------------------------------------------------------------------*/


void run_crtm_modis(rtm_profiles *rtm, profile_params *nwp, rad_utils rutil,
                    int nlevels, float *surface_emissivity, float satzen,
		    int8 *chflg_modis, int platform_id, char *algData_dir, int *pr_year,
                    int *pr_jday, rtm_toa *rclr, int iclr)
{

  extern void
    crtm_wrapper_mp_run_forward_model_(int *, int *, double *, double *, double *,
                                       double *, double *, double *, double *, double *,
				       double *, double *, double *, double *,
				       double *, double *, double *, double *, double *,
				       double *, double (*)[]);


}

/*----------------------------------------------------------------------------
  Call the PLOD for MODIS using this routine.
----------------------------------------------------------------------------*/

void run_plod_modis(rtm_profiles *rtm, profile_params *nwp, rad_utils rutil,
                    int nlevels, float *surface_emissivity, float satzen,
		    int8 *chflg_modis, int platform_id, char *algData_dir, int *pr_year,
                    int *pr_jday, rtm_toa *rclr, int iclr)
{

  char *rout = {"run_plod_modis.c"};
  int idet, iok, il, isfc, ch;
  double a;

  idet = 0;

//printf("RTM routines: algData_dir: %s \n", algData_dir);
  char plod_coeff_dir[256];
  sprintf(plod_coeff_dir,"%s/%s",algData_dir,"MODIS/rtm_plod");
//printf("RTM PLOD coeff directory: %s\n", plod_coeff_dir);
//(void)fflush(stdout);

//printf("year jday: %d %d\n", *pr_year, *pr_jday);
  isfc = nwp->sfc_level;
  il = isfc - 1;
  a = (nwp->p_sfc - nwp->p_lev[il]) / (nwp->p_lev[il+1] - nwp->p_lev[il]);
//nwp->p_lev[isfc] = nwp->p_sfc;
  nwp->t_lev[isfc] = nwp->t_lev[il] + a * (nwp->t_lev[il+1] - nwp->t_lev[il]);

/*
  int j;
  printf("RTM: %d %d %f %f %f\n", isfc, il, a, nwp->p_lev[isfc], nwp->t_lev[isfc]);
  for (j=0; j<101; j++) {
      printf("%d %f %f %f %f \n", j, nwp->p_lev[j], nwp->t_lev[j],nwp->w_lev[j],nwp->adjo3_lev[j]);
  }
*/ 

  if (chflg_modis[19]) {
    ch = 20;
    if ((rtm->trans_atm20_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm20_clr");
    if ((rtm->rad_atm20_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm20_clr");
    if ((rtm->cloud_prof20 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof20");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm20_clr,&iok);
    /*tran_modisd101_(algData_dir,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm20_clr,&iok);*/
	tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm20_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 20 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm20_clr[isfc] = rtm->trans_atm20_clr[il] + a * (rtm->trans_atm20_clr[il+1] - rtm->trans_atm20_clr[il]);

//  printf("Calling clear_atm_rad for band 20\n");
//  (void)fflush(stdout);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss20[iclr], rtm->trans_atm20_clr, rtm->rad_atm20_clr, rtm->cloud_prof20,
		   &rtm->rad20_clr, &rtm->bt20_clr);
  }
//printf("Done with RTM for band 20\n");
//(void)fflush(stdout);

  if (chflg_modis[20]) {
    ch = 21;
    if ((rtm->trans_atm21_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm21_clr");
    if ((rtm->rad_atm21_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm21_clr");
    if ((rtm->cloud_prof21 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof21");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm21_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm21_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 21 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm21_clr[isfc] = rtm->trans_atm21_clr[il] + a * (rtm->trans_atm21_clr[il+1] - rtm->trans_atm21_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss21[iclr], rtm->trans_atm21_clr, rtm->rad_atm21_clr, rtm->cloud_prof21,
		   &rtm->rad21_clr, &rtm->bt21_clr);
  }

  if (chflg_modis[21]) {
    ch = 22;
    if ((rtm->trans_atm22_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm22_clr");
    if ((rtm->rad_atm22_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm22_clr");
    if ((rtm->cloud_prof22 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof22");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm22_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm22_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 22 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm22_clr[isfc] = rtm->trans_atm22_clr[il] + a * (rtm->trans_atm22_clr[il+1] - rtm->trans_atm22_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss22[iclr], rtm->trans_atm22_clr, rtm->rad_atm22_clr, rtm->cloud_prof22,
		   &rtm->rad22_clr, &rtm->bt22_clr);
  }

  if (chflg_modis[22]) {
    ch = 23;
    if ((rtm->trans_atm23_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm23_clr");
    if ((rtm->rad_atm23_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm23_clr");
    if ((rtm->cloud_prof23 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof23");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm23_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm23_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 23 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm23_clr[isfc] = rtm->trans_atm23_clr[il] + a * (rtm->trans_atm23_clr[il+1] - rtm->trans_atm23_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss23[iclr], rtm->trans_atm23_clr, rtm->rad_atm23_clr, rtm->cloud_prof23,
		   &rtm->rad23_clr, &rtm->bt23_clr);
  }

  if (chflg_modis[23]) {
    ch = 24;
    if ((rtm->trans_atm24_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm24_clr");
    if ((rtm->rad_atm24_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm24_clr");
    if ((rtm->cloud_prof24 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof24");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm24_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm24_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 24 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm24_clr[isfc] = rtm->trans_atm24_clr[il] + a * (rtm->trans_atm24_clr[il+1] - rtm->trans_atm24_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss24[iclr], rtm->trans_atm24_clr, rtm->rad_atm24_clr, rtm->cloud_prof24,
		   &rtm->rad24_clr, &rtm->bt24_clr);
  }

  if (chflg_modis[24]) {
    ch = 25;
    if ((rtm->trans_atm25_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm25_clr");
    if ((rtm->rad_atm25_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm25_clr");
    if ((rtm->cloud_prof25 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof25");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm25_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm25_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 25 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm25_clr[isfc] = rtm->trans_atm25_clr[il] + a * (rtm->trans_atm25_clr[il+1] - rtm->trans_atm25_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss25[iclr], rtm->trans_atm25_clr, rtm->rad_atm25_clr, rtm->cloud_prof25,
		   &rtm->rad25_clr, &rtm->bt25_clr);
  }

  if (chflg_modis[26]) {
    ch = 27;
    if ((rtm->trans_atm27_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm27_clr");
    if ((rtm->rad_atm27_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm27_clr");
    if ((rtm->cloud_prof27 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof27");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm27_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm27_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 27 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm27_clr[isfc] = rtm->trans_atm27_clr[il] + a * (rtm->trans_atm27_clr[il+1] - rtm->trans_atm27_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss27[iclr], rtm->trans_atm27_clr, rtm->rad_atm27_clr, rtm->cloud_prof27,
		   &rtm->rad27_clr, &rtm->bt27_clr);
  }

  if (chflg_modis[27]) {
    ch = 28;
    if ((rtm->trans_atm28_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm28_clr");
    if ((rtm->rad_atm28_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm28_clr");
    if ((rtm->cloud_prof28 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof28");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm28_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm28_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 28 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm28_clr[isfc] = rtm->trans_atm28_clr[il] + a * (rtm->trans_atm28_clr[il+1] - rtm->trans_atm28_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss28[iclr], rtm->trans_atm28_clr, rtm->rad_atm28_clr, rtm->cloud_prof28,
		   &rtm->rad28_clr, &rtm->bt28_clr);
  }

  if (chflg_modis[28]) {
    ch = 29;
    if ((rtm->trans_atm29_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm29_clr");
    if ((rtm->rad_atm29_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm29_clr");
    if ((rtm->cloud_prof29 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof29");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm29_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm29_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 29 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm29_clr[isfc] = rtm->trans_atm29_clr[il] + a * (rtm->trans_atm29_clr[il+1] - rtm->trans_atm29_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss29[iclr], rtm->trans_atm29_clr, rtm->rad_atm29_clr, rtm->cloud_prof29,
		   &rtm->rad29_clr, &rtm->bt29_clr);
    /*printf("run_plod_modis: ch: %d t_sfc: %f sfc_emiss%d: %f trans_atm%d_clr[isfc]:%f rad%d_clr:%f bt%d_clr: %f\n", ch, nwp->t_sfc, ch, rclr->sfc_emiss29[iclr], 
          ch, rtm->trans_atm29_clr[isfc], ch, rtm->rad29_clr,ch, rtm->bt29_clr);*/
  }

  if (chflg_modis[29]) {
    ch = 30;
    if ((rtm->trans_atm30_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm30_clr");
    if ((rtm->rad_atm30_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm30_clr");
    if ((rtm->cloud_prof30 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof30");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm30_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm30_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 30 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm30_clr[isfc] = rtm->trans_atm30_clr[il] + a * (rtm->trans_atm30_clr[il+1] - rtm->trans_atm30_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss30[iclr], rtm->trans_atm30_clr, rtm->rad_atm30_clr, rtm->cloud_prof30,
		   &rtm->rad30_clr, &rtm->bt30_clr);
  }

  if (chflg_modis[30]) {
    ch = 31;
    if ((rtm->trans_atm31_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm31_clr");
    if ((rtm->rad_atm31_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm31_clr");
    if ((rtm->cloud_prof31 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof31");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm31_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm31_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 31 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm31_clr[isfc] = rtm->trans_atm31_clr[il] + a * (rtm->trans_atm31_clr[il+1] - rtm->trans_atm31_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss31[iclr], rtm->trans_atm31_clr, rtm->rad_atm31_clr, rtm->cloud_prof31,
		   &rtm->rad31_clr, &rtm->bt31_clr);
    /*printf("run_plod_modis: ch: %d t_sfc: %f sfc_emiss%d: %f trans_atm%d_clr[isfc]:%f rad%d_clr:%f bt%d_clr: %f\n", ch, nwp->t_sfc, ch, rclr->sfc_emiss31[iclr], 
          ch, rtm->trans_atm31_clr[isfc], ch, rtm->rad31_clr,ch, rtm->bt31_clr);*/

  }

  if (chflg_modis[31]) {
    ch = 32;
    if ((rtm->trans_atm32_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm32_clr");
    if ((rtm->rad_atm32_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm32_clr");
    if ((rtm->cloud_prof32 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof32");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm32_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm32_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 32 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm32_clr[isfc] = rtm->trans_atm32_clr[il] + a * (rtm->trans_atm32_clr[il+1] - rtm->trans_atm32_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss32[iclr], rtm->trans_atm32_clr, rtm->rad_atm32_clr, rtm->cloud_prof32,
		   &rtm->rad32_clr, &rtm->bt32_clr);
   /* printf("run_plod_modis: ch: %d t_sfc: %f sfc_emiss%d: %f trans_atm%d_clr[isfc]:%f rad%d_clr:%f bt%d_clr: %f\n", ch, nwp->t_sfc, ch, rclr->sfc_emiss32[iclr], 
          ch, rtm->trans_atm32_clr[isfc], ch, rtm->rad32_clr,ch, rtm->bt32_clr);*/
  }

  if (chflg_modis[32]) {
    ch = 33;
    if ((rtm->trans_atm33_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm33_clr");
    if ((rtm->rad_atm33_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm33_clr");
    if ((rtm->cloud_prof33 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof33");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm33_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm33_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 33 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm33_clr[isfc] = rtm->trans_atm33_clr[il] + a * (rtm->trans_atm33_clr[il+1] - rtm->trans_atm33_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss33[iclr], rtm->trans_atm33_clr, rtm->rad_atm33_clr, rtm->cloud_prof33,
		   &rtm->rad33_clr, &rtm->bt33_clr);
    /*printf("run_plod_modis: ch: %d t_sfc: %f sfc_emiss%d: %f trans_atm%d_clr[isfc]:%f rad%d_clr:%f bt%d_clr: %f\n", ch, nwp->t_sfc, ch, rclr->sfc_emiss33[iclr], 
          ch, rtm->trans_atm33_clr[isfc], ch, rtm->rad33_clr,ch, rtm->bt33_clr);*/
  }

  if (chflg_modis[33]) {
    ch = 34;
    if ((rtm->trans_atm34_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm34_clr");
    if ((rtm->rad_atm34_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm34_clr");
    if ((rtm->cloud_prof34 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof34");
//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm34_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm34_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 34 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm34_clr[isfc] = rtm->trans_atm34_clr[il] + a * (rtm->trans_atm34_clr[il+1] - rtm->trans_atm34_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss34[iclr], rtm->trans_atm34_clr, rtm->rad_atm34_clr, rtm->cloud_prof34,
		   &rtm->rad34_clr, &rtm->bt34_clr);
    /* printf("run_plod_modis: ch: %d t_sfc: %f sfc_emiss%d: %f trans_atm%d_clr[isfc]:%f rad%d_clr:%f bt%d_clr: %f\n", ch, nwp->t_sfc, ch, rclr->sfc_emiss34[iclr], 
          ch, rtm->trans_atm34_clr[isfc], ch, rtm->rad34_clr,ch, rtm->bt34_clr);*/
  }

  if (chflg_modis[34]) {
    ch = 35;
    if ((rtm->trans_atm35_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm35_clr");
    if ((rtm->rad_atm35_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm35_clr");
    if ((rtm->cloud_prof35 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof35");

//  printf("Calling tran_modisd101 for band 35\n");
//  (void)fflush(stdout);
//  (void) sleep(10);

//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm35_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm35_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 35 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm35_clr[isfc] = rtm->trans_atm35_clr[il] + a * (rtm->trans_atm35_clr[il+1] - rtm->trans_atm35_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss35[iclr], rtm->trans_atm35_clr, rtm->rad_atm35_clr, rtm->cloud_prof35,
		   &rtm->rad35_clr, &rtm->bt35_clr);
    /*printf("run_plod_modis: ch: %d t_sfc: %f sfc_emiss%d: %f trans_atm%d_clr[isfc]:%f rad%d_clr:%f bt%d_clr: %f\n", ch, nwp->t_sfc, ch, rclr->sfc_emiss35[iclr], 
          ch, rtm->trans_atm35_clr[isfc], ch, rtm->rad35_clr,ch, rtm->bt35_clr);*/
  }

  if (chflg_modis[35]) {
    ch = 36;
    if ((rtm->trans_atm36_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->trans_atm36_clr");
    if ((rtm->rad_atm36_clr = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->rad_atm36_clr");
    if ((rtm->cloud_prof36 = (float *) calloc(nlevels,sizeof(float))) == NULL)
      error_allo(rout,"rtm->cloud_prof36");

//  printf("Calling tran_modisd101 for band 36\n");
//  (void)fflush(stdout);
//  sleep(10);

//  tran_modisd101_(nwp->t_lev,nwp->w_lev,nwp->o3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm36_clr,&iok);
    tran_modisd101_(algData_dir,pr_year,pr_jday,nwp->t_lev,nwp->w_lev,nwp->adjo3_lev,&satzen,&platform_id,&ch,&idet,rtm->trans_atm36_clr,&iok);
    if (iok != 0) {
      fprintf(stderr,"%sError calculating channel 36 PLOD transmittances - aborting\n",EXE_PROMPT);
      exit(EXIT_FAILURE);
    }
    rtm->trans_atm36_clr[isfc] = rtm->trans_atm36_clr[il] + a * (rtm->trans_atm36_clr[il+1] - rtm->trans_atm36_clr[il]);
    clear_atm_rad (rutil, nwp->t_lev, nwp->t_sfc, nwp->sfc_level+1,
                   ch, rclr->sfc_emiss36[iclr], rtm->trans_atm36_clr, rtm->rad_atm36_clr, rtm->cloud_prof36,
		   &rtm->rad36_clr, &rtm->bt36_clr);
    /*printf("run_plod_modis: ch: %d t_sfc: %f sfc_emiss%d: %f trans_atm%d_clr[isfc]:%f rad%d_clr:%f bt%d_clr: %f\n", ch, nwp->t_sfc, ch, rclr->sfc_emiss36[iclr], 
          ch, rtm->trans_atm36_clr[isfc], ch, rtm->rad36_clr,ch, rtm->bt36_clr);*/
  }


}

void destroy_rtm_ptrs_modis(rtm_profiles * rtm, nwp_params nwp, int8 *chflg_modis)
{
  if (chflg_modis[19]) {
    free(rtm->trans_atm20_clr);
    free(rtm->rad_atm20_clr);
    free(rtm->cloud_prof20);
  }

  if (chflg_modis[20]) {
    free(rtm->trans_atm21_clr);
    free(rtm->rad_atm21_clr);
    free(rtm->cloud_prof21);
  }

  if (chflg_modis[21]) {
    free(rtm->trans_atm22_clr);
    free(rtm->rad_atm22_clr);
    free(rtm->cloud_prof22);
  }

  if (chflg_modis[22]) {
    free(rtm->trans_atm23_clr);
    free(rtm->rad_atm23_clr);
    free(rtm->cloud_prof23);
  }

  if (chflg_modis[23]) {
    free(rtm->trans_atm24_clr);
    free(rtm->rad_atm24_clr);
    free(rtm->cloud_prof24);
  }

  if (chflg_modis[24]) {
    free(rtm->trans_atm25_clr);
    free(rtm->rad_atm25_clr);
    free(rtm->cloud_prof25);
  }

  if (chflg_modis[26]) {
    free(rtm->trans_atm27_clr);
    free(rtm->rad_atm27_clr);
    free(rtm->cloud_prof27);
  }

  if (chflg_modis[27]) {
    free(rtm->trans_atm28_clr);
    free(rtm->rad_atm28_clr);
    free(rtm->cloud_prof28);
  }

  if (chflg_modis[28]) {
    free(rtm->trans_atm29_clr);
    free(rtm->rad_atm29_clr);
    free(rtm->cloud_prof29);
  }

  if (chflg_modis[29]) {
    free(rtm->trans_atm30_clr);
    free(rtm->rad_atm30_clr);
    free(rtm->cloud_prof30);
  }

  if (chflg_modis[30]) {
    free(rtm->trans_atm31_clr);
    free(rtm->rad_atm31_clr);
    free(rtm->cloud_prof31);
  }

  if (chflg_modis[31]) {
    free(rtm->trans_atm32_clr);
    free(rtm->rad_atm32_clr);
    free(rtm->cloud_prof32);
  }

  if (chflg_modis[32]) {
    free(rtm->trans_atm33_clr);
    free(rtm->rad_atm33_clr);
    free(rtm->cloud_prof33);
  }

  if (chflg_modis[33]) {
    free(rtm->trans_atm34_clr);
    free(rtm->rad_atm34_clr);
    free(rtm->cloud_prof34);
  }

  if (chflg_modis[34]) {
    free(rtm->trans_atm35_clr);
    free(rtm->rad_atm35_clr);
    free(rtm->cloud_prof35);
  }

  if (chflg_modis[35]) {
    free(rtm->trans_atm36_clr);
    free(rtm->rad_atm36_clr);
    free(rtm->cloud_prof36);
  }
}

void get_clear_toa_rad (rtm_profiles rtm, profile_params nwp, rad_utils rutil,
                        int nlevels, rtm_toa *rclr, int8 *chflg, int iclr)
{

//  fprintf(stdout,"%sRTM 0\n",EXE_PROMPT);
//	printf("chflg: %d %d \n", chflg[19],nwp.sfc_level);
//	printf("chflg: %f %f \n", nwp.t_sfc,rclr->sfc_emiss20[iclr]);
//	printf("chflg: %f \n", rtm.trans_atm20_clr[nwp.sfc_level]);
//	printf("chflg: %f \n", rtm.rad_atm20_clr[nwp.sfc_level]);

    if (chflg[19]) {
    clear_toa_rad (20, rtm.rad_atm20_clr[nwp.sfc_level], rtm.trans_atm20_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss20[iclr], rutil, &rclr->rad20_clr[iclr], &rclr->bt20_clr[iclr]);
  }


  if (chflg[20]) {
    clear_toa_rad (21, rtm.rad_atm21_clr[nwp.sfc_level], rtm.trans_atm21_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss21[iclr], rutil, &rclr->rad21_clr[iclr], &rclr->bt21_clr[iclr]);
  }
  if (chflg[21]) {
    clear_toa_rad (22, rtm.rad_atm22_clr[nwp.sfc_level], rtm.trans_atm22_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss22[iclr], rutil, &rclr->rad22_clr[iclr], &rclr->bt22_clr[iclr]);
  }
  if (chflg[22]) {
    clear_toa_rad (23, rtm.rad_atm23_clr[nwp.sfc_level], rtm.trans_atm23_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss23[iclr], rutil, &rclr->rad23_clr[iclr], &rclr->bt23_clr[iclr]);
  }
  if (chflg[23]) {
    clear_toa_rad (24, rtm.rad_atm24_clr[nwp.sfc_level], rtm.trans_atm24_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss24[iclr], rutil, &rclr->rad24_clr[iclr], &rclr->bt24_clr[iclr]);
  }
  if (chflg[24]) {
    clear_toa_rad (25, rtm.rad_atm25_clr[nwp.sfc_level], rtm.trans_atm25_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss25[iclr], rutil, &rclr->rad25_clr[iclr], &rclr->bt25_clr[iclr]);
  }
  if (chflg[26]) {
    clear_toa_rad (27, rtm.rad_atm27_clr[nwp.sfc_level], rtm.trans_atm27_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss27[iclr], rutil, &rclr->rad27_clr[iclr], &rclr->bt27_clr[iclr]);
  }

  if (chflg[27]) {
    clear_toa_rad (28, rtm.rad_atm28_clr[nwp.sfc_level], rtm.trans_atm28_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss28[iclr], rutil, &rclr->rad28_clr[iclr], &rclr->bt28_clr[iclr]);
  }
  if (chflg[28]) {
    clear_toa_rad (29, rtm.rad_atm29_clr[nwp.sfc_level], rtm.trans_atm29_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss29[iclr], rutil, &rclr->rad29_clr[iclr], &rclr->bt29_clr[iclr]);
  }
  if (chflg[29]) {
    clear_toa_rad (30, rtm.rad_atm30_clr[nwp.sfc_level], rtm.trans_atm30_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss30[iclr], rutil, &rclr->rad30_clr[iclr], &rclr->bt30_clr[iclr]);
  }

  if (chflg[30]) {
    clear_toa_rad (31, rtm.rad_atm31_clr[nwp.sfc_level], rtm.trans_atm31_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss31[iclr], rutil, &rclr->rad31_clr[iclr], &rclr->bt31_clr[iclr]);
  }

//  fprintf(stdout,"%sRTM 0\n",EXE_PROMPT);
//      printf("chflg: %d %d \n", chflg[30],nwp.sfc_level);
//      printf("chflg: %f %f \n", nwp.t_sfc,rclr->sfc_emiss31[iclr]);
//      printf("chflg: %f \n", rtm.trans_atm31_clr[nwp.sfc_level]);
//      printf("chflg: %f \n", rtm.rad_atm31_clr[nwp.sfc_level]);
//      printf("chflg: %f %f \n", rclr->rad31_clr[iclr], rclr->bt31_clr[iclr]);
  
  if (chflg[31]) {
    clear_toa_rad (32, rtm.rad_atm32_clr[nwp.sfc_level], rtm.trans_atm32_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss32[iclr], rutil, &rclr->rad32_clr[iclr], &rclr->bt32_clr[iclr]);
  }
  if (chflg[32]) {
    clear_toa_rad (33, rtm.rad_atm33_clr[nwp.sfc_level], rtm.trans_atm33_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss33[iclr], rutil, &rclr->rad33_clr[iclr], &rclr->bt33_clr[iclr]);
  }

  if (chflg[33]) {
    clear_toa_rad (34, rtm.rad_atm34_clr[nwp.sfc_level], rtm.trans_atm34_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss34[iclr], rutil, &rclr->rad34_clr[iclr], &rclr->bt34_clr[iclr]);
  }
  if (chflg[34]) {
    clear_toa_rad (35, rtm.rad_atm35_clr[nwp.sfc_level], rtm.trans_atm35_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss35[iclr], rutil, &rclr->rad35_clr[iclr], &rclr->bt35_clr[iclr]);
  }
  if (chflg[35]) {
    clear_toa_rad (36, rtm.rad_atm36_clr[nwp.sfc_level], rtm.trans_atm36_clr[nwp.sfc_level], nwp.t_sfc,
                   rclr->sfc_emiss36[iclr], rutil, &rclr->rad36_clr[iclr], &rclr->bt36_clr[iclr]);
  }
}

void clear_atm_rad (rad_utils rutil, float *tt, float tsfc, int nlevels,
                    int ich, float surface_emissivity,
		    float *trn, float *rad_atm_clr, float *cloud_prof,
		    float *rad_clr, float *bt_clr)
{

  int ilevel;
  float B1, B2, dtrn;

  B1 = rutil.planck_rad_fast_ptr(ich, tt[0], rutil.T_planck, rutil.B_table);
  rad_atm_clr[0] = 0.0;
  cloud_prof[0] = B1*trn[0];

  for (ilevel=1; ilevel<nlevels; ilevel++) {
    B2 = rutil.planck_rad_fast_ptr(ich, tt[ilevel], rutil.T_planck, rutil.B_table);
    dtrn = -(trn[ilevel] - trn[ilevel-1]);
    rad_atm_clr[ilevel] = rad_atm_clr[ilevel-1] +
      (B1 + B2)/2. * dtrn;
    cloud_prof[ilevel] = rad_atm_clr[ilevel] + B2*trn[ilevel];
    B1 = B2;
  }

  clear_toa_rad (ich, rad_atm_clr[nlevels-1], trn[nlevels-1], tsfc,
                 surface_emissivity, rutil, rad_clr, bt_clr);

}

void clear_toa_rad (int ich, float rad_atm, float tau_atm, float tsfc,
                    float esfc, rad_utils rutil, float *rad_clr, float *bt_clr)
{

  *rad_clr = rad_atm +
    esfc*rutil.planck_rad_fast_ptr(ich,tsfc, rutil.T_planck, rutil.B_table)*tau_atm;

  *bt_clr = rutil.planck_bt_fast_ptr(ich,*rad_clr, rutil.T_planck, rutil.B_table);
}

void create_rclr_ptrs(int8 * chflg, rtm_toa *rclr)
{

  char *rout = {"create_rclr_ptrs"};

  if (chflg[19]) {
    if ((rclr->rad20_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad20_clr");
    if ((rclr->bt20_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt20_clr");
  }
  if (chflg[20]) {
    if ((rclr->rad21_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad21_clr");
    if ((rclr->bt21_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt21_clr");
  }
  if (chflg[21]) {
    if ((rclr->rad22_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad22_clr");
    if ((rclr->bt22_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt22_clr");
  }
  if (chflg[22]) {
    if ((rclr->rad23_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad23_clr");
    if ((rclr->bt23_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt23_clr");
  }
  if (chflg[23]) {
    if ((rclr->rad24_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad24_clr");
    if ((rclr->bt24_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt24_clr");
  }
  if (chflg[24]) {
    if ((rclr->rad25_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad25_clr");
    if ((rclr->bt25_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt25_clr");
  }
  if (chflg[26]) {
    if ((rclr->rad27_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad27_clr");
    if ((rclr->bt27_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt27_clr");
  }
  if (chflg[27]) {
    if ((rclr->rad28_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad28_clr");
    if ((rclr->bt28_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt28_clr");
  }
  if (chflg[28]) {
    if ((rclr->rad29_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad29_clr");
    if ((rclr->bt29_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt29_clr");
  }
  if (chflg[29]) {
    if ((rclr->rad30_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad30_clr");
    if ((rclr->bt30_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt30_clr");
  }
  if (chflg[30]) {
    if ((rclr->rad31_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad31_clr");
    if ((rclr->bt31_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt31_clr");
  }
  if (chflg[31]) {
    if ((rclr->rad32_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad32_clr");
    if ((rclr->bt32_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt32_clr");
  }
  if (chflg[32]) {
    if ((rclr->rad33_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad33_clr");
    if ((rclr->bt33_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt33_clr");
  }
  if (chflg[33]) {
    if ((rclr->rad34_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad34_clr");
    if ((rclr->bt34_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt34_clr");
  }
  if (chflg[34]) {
    if ((rclr->rad35_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad35_clr");
    if ((rclr->bt35_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt35_clr");
  }
  if (chflg[35]) {
    if ((rclr->rad36_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->rad36_clr");
    if ((rclr->bt36_clr = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->bt36_clr");
  }
}

void destroy_rclr_ptrs(int8 * chflg, rtm_toa *rclr)
{



  free(rclr->flag);

  if (chflg[19]) {
    free(rclr->rad20_clr), free(rclr->bt20_clr), free(rclr->sfc_emiss20);
  }
  if (chflg[20]) {
    free(rclr->rad21_clr), free(rclr->bt21_clr), free(rclr->sfc_emiss21);
  }
  if (chflg[21]) {
    free(rclr->rad22_clr), free(rclr->bt22_clr), free(rclr->sfc_emiss22);
  }
  if (chflg[22]) {
    free(rclr->rad23_clr), free(rclr->bt23_clr), free(rclr->sfc_emiss23);
  }
  if (chflg[23]) {
    free(rclr->rad24_clr), free(rclr->bt24_clr), free(rclr->sfc_emiss24);
  }
  if (chflg[24]) {
    free(rclr->rad25_clr), free(rclr->bt25_clr), free(rclr->sfc_emiss25);
  }
  if (chflg[26]) {
    free(rclr->rad27_clr), free(rclr->bt27_clr), free(rclr->sfc_emiss27);
  }
  if (chflg[27]) {
    free(rclr->rad28_clr), free(rclr->bt28_clr), free(rclr->sfc_emiss28);
  }
  if (chflg[28]) {
    free(rclr->rad29_clr), free(rclr->bt29_clr), free(rclr->sfc_emiss29);
  }
  if (chflg[29]) {
    free(rclr->rad30_clr), free(rclr->bt30_clr), free(rclr->sfc_emiss30);
  }
  if (chflg[30]) {
    free(rclr->rad31_clr), free(rclr->bt31_clr), free(rclr->sfc_emiss31);
  }
  if (chflg[31]) {
    free(rclr->rad32_clr), free(rclr->bt32_clr), free(rclr->sfc_emiss32);
  }
  if (chflg[32]) {
    free(rclr->rad33_clr), free(rclr->bt33_clr), free(rclr->sfc_emiss33);
  }
  if (chflg[33]) {
    free(rclr->rad34_clr), free(rclr->bt34_clr), free(rclr->sfc_emiss34);
  }
  if (chflg[34]) {
    free(rclr->rad35_clr), free(rclr->bt35_clr), free(rclr->sfc_emiss35);
  }
  if (chflg[35]) {
    free(rclr->rad36_clr), free(rclr->bt36_clr), free(rclr->sfc_emiss36);
  }
}
