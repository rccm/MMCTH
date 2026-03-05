/*$Id: create_algorithm_pointers.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

#include <stdio.h>
#include <stdlib.h>
#include <hdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "imagerL2_leocat.h"

void create_algorithm_pointers(int16 nout, char *sds_name[MAX_ALGO_OUTPUT],
                               imagerL2 *imgr2)
{
char *rout = {"create_algorithm_pointers.c"};
int i;

  for (i=0; i<nout; i++) {
    if (strcmp(sds_name[i],"cloud_mask") == 0) {
      printf("Creating Algorithm pointer for cloud mask of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->cldmask = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->cldmask");
    }
    else if (strcmp(sds_name[i],"land_mask") == 0) {
      printf("Creating Algorithm pointer for land mask of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->landmask = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->landmask");
	  }
    else if (strcmp(sds_name[i],"cloud_type") == 0) {
      printf("Creating Algorithm pointer for cloud type of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->cldtype = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->cldtype");
    }
    else if (strcmp(sds_name[i],"cloud_phase") == 0) {
      printf("Creating Algorithm pointer for cloud phase of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->cldphase = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->cldphase");
    }
    else if (strcmp(sds_name[i],"aerosol_mask") == 0) {
      if ((imgr2->aeromask = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->aeromask");
    }
    else if (strcmp(sds_name[i],"so2mask") == 0) {
      if ((imgr2->so2mask = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->so2mask");
    }
    else if (strcmp(sds_name[i],"cop_qf0") == 0) {	// GPC
      if ((imgr2->cop_qf0 = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->cop_qf0");
    }
    else if (strcmp(sds_name[i],"cop_qf1") == 0) {	// GPC
      if ((imgr2->cop_qf1 = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->cop_qf1");
    }
    else if (strcmp(sds_name[i],"inwctt_qf") == 0) {	// GPC
      if ((imgr2->inwctt_qf = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->inwctt_qf");
    }
    else if (strcmp(sds_name[i],"flag_arr1") == 0) {
      printf("Creating Algorithm pointer for flag_arr1 of %d bytes times %ld points\n",10,imgr2->npts);	// GPC
      imgr2->flag_arr1 = calloc_2d_uchar8_ptr("imgr2->flag_arr1", 10, imgr2->npts);
      if ((imgr2->stats_arr = (float *) malloc(18*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->stats_arr");
    }
    else if (strcmp(sds_name[i],"flag_arr2") == 0) {
      printf("Creating Algorithm pointer for flag_arr2 of %d bytes times %ld points\n",10,imgr2->npts);	// GPC
      imgr2->flag_arr2 = calloc_2d_uchar8_ptr("imgr2->flag_arr2", imgr2->npts, 10);
    }
//  Also check whether or not 250-m data will be output.  RAF
    else if (strcmp(sds_name[i],"stats_250m") == 0 && imgr2->output_250 == YES) {
      printf("Creating Algorithm pointer for stats_250m of %d times %ld points\n",4,imgr2->npts);	// RAF
      imgr2->stats_250m = calloc_2d_float_ptr("imgr2->stats_250m", imgr2->npts, 4);
    }
    else if (strcmp(sds_name[i],"ctParmQf") == 0) {
//    printf("Creating Algorithm pointer for ctParmQf of %d bytes times %ld points\n",3,imgr2->npts);	// GPC
      imgr2->ctParmQf = calloc_2d_uchar8_ptr("imgr2->ctParmQf", 3, imgr2->npts);
    }
    else if (strcmp(sds_name[i],"sst") == 0) {
      if ((imgr2->sst = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->sst");
    }
    else if (strcmp(sds_name[i],"lst") == 0) {
      if ((imgr2->lst = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->lst");
    }
    else if (strcmp(sds_name[i],"ist") == 0) {
      if ((imgr2->ist = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->ist");
    }
    else if (strcmp(sds_name[i],"ndvi") == 0) {
      if ((imgr2->ndvi = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->ndvi");
    }
    else if (strcmp(sds_name[i],"tpw") == 0) {
      if ((imgr2->tpw = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->tpw");
    }
    else if (strcmp(sds_name[i],"cloud_top_temperature") == 0) {
      if ((imgr2->cldt = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldt");
    }
    else if (strcmp(sds_name[i],"cloudTopTemp_Ice_NightWater") == 0) {	//GPC
      printf("Creating Algorithm pointer for cloudTopTemp_Ice_NightWater of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->cldt_IceNW = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)		// GPC
        error_allo(rout,"imgr2->cldt_IceNW");
    }
    else if (strcmp(sds_name[i],"aerosol_top_temperature") == 0) {
      if ((imgr2->aerot = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aerot");
    }
    else if (strcmp(sds_name[i],"cloud_top_temperature_high") == 0) {
      if ((imgr2->cldt_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldt_high");
    }
    else if (strcmp(sds_name[i],"cloud_top_temperature_low") == 0) {
      if ((imgr2->cldt_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldt_low");
    }
    else if (strcmp(sds_name[i],"cloud_top_pressure") == 0) {
      if ((imgr2->cldp = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldp");
    }
    else if (strcmp(sds_name[i],"aerosol_top_pressure") == 0) {
      if ((imgr2->aerop = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aerop");
    }
    else if (strcmp(sds_name[i],"cloud_top_method") == 0) {
      if ((imgr2->cld_solm = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->cld_solm");
    }
    else if (strcmp(sds_name[i],"aerosol_top_method") == 0) {
      if ((imgr2->aero_solm = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->aero_solm");
    }
    else if (strcmp(sds_name[i],"cloud_top_pressure_high") == 0) {
      if ((imgr2->cldp_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldp_high");
    }
    else if (strcmp(sds_name[i],"cloud_top_pressure_low") == 0) {
      if ((imgr2->cldp_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldp_low");
    }
    else if (strcmp(sds_name[i],"cloud_top_height") == 0) {
      if ((imgr2->cldz = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldz");
    }
// Added RAF
    else if (strcmp(sds_name[i],"clear_bt31") == 0) {
      if ((imgr2->clrbt31 = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->clrbt");
    }
// Added RAF
    else if (strcmp(sds_name[i],"surface_temperature") == 0) {
      printf("Creating Algorithm float pointer for imgr2->sfctmp of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->sfctmp = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->sfctmp");
    }
// Added RAF
    else if (strcmp(sds_name[i],"cld_emiss11") == 0) {
      if ((imgr2->cld_emiss11 = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cld_emiss11");
    }
// Added RAF
    else if (strcmp(sds_name[i],"cld_emiss12") == 0) {
      if ((imgr2->cld_emiss12 = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cld_emiss12");
    }
// Added RAF
    else if (strcmp(sds_name[i],"cld_emiss13") == 0) {
      if ((imgr2->cld_emiss13 = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cld_emiss13");
    }
// Added RAF
    else if (strcmp(sds_name[i],"cld_emiss85") == 0) {
      if ((imgr2->cld_emiss85 = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cld_emiss85");
    }
// Added RAF
    else if (strcmp(sds_name[i],"modis_C6_IRP") == 0) {
      if ((imgr2->modis_C6_IRP = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->modis_C6_IRP");
    }
// Added RAF
    else if (strcmp(sds_name[i],"IRP_CTH_Consistency_Flag") == 0) {
      if ((imgr2->IRP_CTH_Consistency_Flag = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->IRP_CTH_Consistency_Flag");
    }
// Added RAF 
    else if (strcmp(sds_name[i],"os_top_flag") == 0) {
      if ((imgr2->os_top_flag = (unsigned char *) malloc(imgr2->npts*sizeof(unsigned char))) == NULL)
        error_allo(rout,"imgr2->os_top_flag");
    }
    else if (strcmp(sds_name[i],"aerosol_top_height") == 0) {
      if ((imgr2->aeroz = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aeroz");
    }
    else if (strcmp(sds_name[i],"cloud_top_height_high") == 0) {
      if ((imgr2->cldz_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldz_high");
    }
    else if (strcmp(sds_name[i],"cloud_top_height_low") == 0) {
      if ((imgr2->cldz_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldz_low");
    }
    else if (strcmp(sds_name[i],"cloud_particle_effective_radius") == 0) {
      printf("Creating Algorithm pointer for cloud_particle_effective_radius of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->cldreff = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldreff");
    }
    else if (strcmp(sds_name[i],"cloud_particle_effective_radius_high") == 0) {
      if ((imgr2->cldreff_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldreff_high");
    }
    else if (strcmp(sds_name[i],"cloud_particle_effective_radius_low") == 0) {
      if ((imgr2->cldreff_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldreff_low");
    }
    else if (strcmp(sds_name[i],"cloud_particle_effective_diameter") == 0) {
      if ((imgr2->clddeff = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->clddeff");
    }
    else if (strcmp(sds_name[i],"cloud_particle_effective_diameter_high") == 0) {
      if ((imgr2->clddeff_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->clddeff_high");
    }
    else if (strcmp(sds_name[i],"cloud_particle_effective_diameter_low") == 0) {
      if ((imgr2->clddeff_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->clddeff_low");
    }
    else if (strcmp(sds_name[i],"cloud_emissivity") == 0) {
      if ((imgr2->cldemiss = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldemiss");
    }
    else if (strcmp(sds_name[i],"aerosol_emissivity") == 0) {
      if ((imgr2->aeroemiss = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aeroemiss");
    }
    else if (strcmp(sds_name[i],"cloud_emissivity_high") == 0) {
      if ((imgr2->cldemiss_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldemiss_high");
    }
    else if (strcmp(sds_name[i],"cloud_emissivity_low") == 0) {
      if ((imgr2->cldemiss_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldemiss_low");
    }
    else if (strcmp(sds_name[i],"cloud_optical_depth_ir") == 0) {
      printf("Creating Algorithm pointer for cloud_optical_depth_ir of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->cod_ir = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cod_ir");
    }
    else if (strcmp(sds_name[i],"cloud_optical_depth_ir_high") == 0) {
      if ((imgr2->cod_ir_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cod_ir_high");
    }
    else if (strcmp(sds_name[i],"cloud_optical_depth_ir_low") == 0) {
      if ((imgr2->cod_ir_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cod_ir_low");
    }
    else if (strcmp(sds_name[i],"cloud_optical_depth_vis") == 0) {
      printf("Creating Algorithm pointer for cloud_optical_depth_vis of %ld points\n",imgr2->npts);	// GPC
      if ((imgr2->cod_vis = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cod_vis");
    }
    else if (strcmp(sds_name[i],"cloud_optical_depth_vis_high") == 0) {
      if ((imgr2->cod_vis_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cod_vis_high");
    }
    else if (strcmp(sds_name[i],"cloud_optical_depth_vis_low") == 0) {
      if ((imgr2->cod_vis_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cod_vis_low");
    }
    else if (strcmp(sds_name[i],"cloud_bottom") == 0) {
      if ((imgr2->cldbot = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldbot");
    }
    else if (strcmp(sds_name[i],"aerosol_bottom") == 0) {
      if ((imgr2->aerobot = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aerobot");
    }
    else if (strcmp(sds_name[i],"cloud_bottom_high") == 0) {
      if ((imgr2->cldbot_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldbot_high");
    }
    else if (strcmp(sds_name[i],"cloud_bottom_low") == 0) {
      if ((imgr2->cldbot_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldbot_low");
    }
    else if (strcmp(sds_name[i],"cloud_thickness") == 0) {
      if ((imgr2->clddz = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->clddz");
    }
    else if (strcmp(sds_name[i],"aerosol_thickness") == 0) {
      if ((imgr2->aerodz = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aerodz");
    }
    else if (strcmp(sds_name[i],"cloud_thickness_high") == 0) {
      if ((imgr2->clddz_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->clddz_high");
    }
    else if (strcmp(sds_name[i],"cloud_thickness_low") == 0) {
      if ((imgr2->clddz_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->clddz_low");
    }
    else if (strcmp(sds_name[i],"cloud_fraction") == 0) {
      if ((imgr2->cldfrac = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldfrac");
    }
    else if (strcmp(sds_name[i],"aerosol_fraction") == 0) {
      if ((imgr2->aerofrac = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aerofrac");
    }
    else if (strcmp(sds_name[i],"cloud_liquid_water_path") == 0) {
      if ((imgr2->cldlwp = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldlwp");
    }
    else if (strcmp(sds_name[i],"cloud_ice_water_path") == 0) {
      if ((imgr2->cldiwp = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldiwp");
    }
    else if (strcmp(sds_name[i],"cloud_liquid_water_content") == 0) {
      if ((imgr2->cldlwc = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldlwc");
    }
    else if (strcmp(sds_name[i],"cloud_ice_water_content") == 0) {
      if ((imgr2->cldiwc = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldiwc");
    }
    else if (strcmp(sds_name[i],"cloud_beta") == 0) {
      if ((imgr2->cldbeta = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldbeta");
    }
    else if (strcmp(sds_name[i],"cloud_beta_high") == 0) {
      if ((imgr2->cldbeta_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldbeta_high");
    }
    else if (strcmp(sds_name[i],"cloud_beta_mid") == 0) {
      if ((imgr2->cldbeta_mid = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldbeta_mid");
    }
    else if (strcmp(sds_name[i],"cloud_beta_low") == 0) {
      if ((imgr2->cldbeta_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->cldbeta_low");
    }
    else if (strcmp(sds_name[i],"emiss11_high") == 0) {
      if ((imgr2->emiss11_high = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->emiss11_high");
    }
    else if (strcmp(sds_name[i],"emiss11_mid") == 0) {
      if ((imgr2->emiss11_mid = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->emiss11_mid");
    }
    else if (strcmp(sds_name[i],"emiss11_low") == 0) {
      if ((imgr2->emiss11_low = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->emiss11_low");
    }
    else if (strcmp(sds_name[i],"column_aerosol_amount") == 0) {
      if ((imgr2->col_aer = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->col_aer");
    }
    else if (strcmp(sds_name[i],"aerosol_optical_depth_vis") == 0) {
      if ((imgr2->aod_vis = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aod_vis");
    }
    else if (strcmp(sds_name[i],"aerosol_optical_depth_ir") == 0) {
      if ((imgr2->aod_ir = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aod_ir");
    }
    else if (strcmp(sds_name[i],"aerosol_particle_effective_radius") == 0) {
      if ((imgr2->aerreff = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aerreff");
    }
    else if (strcmp(sds_name[i],"aerosol_particle_effective_diameter") == 0) {
      if ((imgr2->aerdeff = (float *) malloc(imgr2->npts*sizeof(float))) == NULL)
        error_allo(rout,"imgr2->aerdeff");
    }
  }
}
