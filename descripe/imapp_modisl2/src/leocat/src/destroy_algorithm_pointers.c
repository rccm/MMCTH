/*$Id: destroy_algorithm_pointers.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

#include <hdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "imagerL2_leocat.h"

void destroy_algorithm_pointers(int16 nout, char *sds_name[MAX_ALGO_OUTPUT], imagerL2 *imgr2)
{
	int i;

	for (i=0; i<nout; i++) {
//              printf("freeing SDS %d: %s\n",i,sds_name[i]);

		if (strcmp(sds_name[i],"cloud_mask") == 0) {
			free(imgr2->cldmask);
		}
		else if (strcmp(sds_name[i],"land_mask") == 0) {
			free(imgr2->landmask);
		}
		else if (strcmp(sds_name[i],"cloud_type") == 0) {
			free(imgr2->cldtype);
		}
		else if (strcmp(sds_name[i],"cloud_phase") == 0) {
			free(imgr2->cldphase);
		}
		else if (strcmp(sds_name[i],"aerosol_mask") == 0) {
			free(imgr2->aeromask);
		}
		else if (strcmp(sds_name[i],"so2mask") == 0) {
			free(imgr2->so2mask);
		}
		else if (strcmp(sds_name[i],"cop_qf0") == 0) {	// GPC
			free(imgr2->cop_qf0);
		}
		else if (strcmp(sds_name[i],"cop_qf1") == 0) {	// GPC
			free(imgr2->cop_qf1);
		}
		else if (strcmp(sds_name[i],"inwctt_qf") == 0) {	// GPC
			free(imgr2->inwctt_qf);
		}
		else if (strcmp(sds_name[i],"flag_arr1") == 0) {
			destroy_2d_uchar8_ptr(10, imgr2->flag_arr1);
		}
		else if (strcmp(sds_name[i],"flag_arr2") == 0) {
			destroy_2d_uchar8_ptr(imgr2->npts, imgr2->flag_arr2);
		}
//              Also check whether or not 250-m data was output.    RAF
		else if (strcmp(sds_name[i],"stats_250m") == 0 && imgr2->output_250 == YES) {
			destroy_2d_float_ptr(imgr2->npts, imgr2->stats_250m);
		}
		else if (strcmp(sds_name[i],"sst") == 0) {
			free(imgr2->sst);
		}
		else if (strcmp(sds_name[i],"lst") == 0) {
			free(imgr2->lst);
		}
		else if (strcmp(sds_name[i],"ist") == 0) {
			free(imgr2->ist);
		}
		else if (strcmp(sds_name[i],"ndvi") == 0) {
			free(imgr2->ndvi);
		}
		else if (strcmp(sds_name[i],"tpw") == 0) {
			free(imgr2->tpw);
		}
		else if (strcmp(sds_name[i],"cloud_top_temperature") == 0) {
			free(imgr2->cldt);
		}
		else if (strcmp(sds_name[i],"cloudTopTemp_Ice_NightWater") == 0) {	// GPC
			free(imgr2->cldt_IceNW);											// GPC
		}																	// GPC
		else if (strcmp(sds_name[i],"ctParmQf") == 0) {						// GPC
			destroy_2d_uchar8_ptr(3, imgr2->ctParmQf);						// GPC
		}																	// GPC
		else if (strcmp(sds_name[i],"aerosol_top_temperature") == 0) {
			free(imgr2->aerot);
		}
		else if (strcmp(sds_name[i],"cloud_top_temperature_high") == 0) {
			free(imgr2->cldt_high);
		}
		else if (strcmp(sds_name[i],"cloud_top_temperature_low") == 0) {
			free(imgr2->cldt_low);
		}
		else if (strcmp(sds_name[i],"cloud_top_pressure") == 0) {
			free(imgr2->cldp);
		}
		else if (strcmp(sds_name[i],"aerosol_top_pressure") == 0) {
			free(imgr2->aerop);
		}
		else if (strcmp(sds_name[i],"cloud_top_method") == 0) {
			free(imgr2->cld_solm);
		}
		else if (strcmp(sds_name[i],"aerosol_top_method") == 0) {
			free(imgr2->aero_solm);
		}
		else if (strcmp(sds_name[i],"cloud_top_pressure_high") == 0) {
			free(imgr2->cldp_high);
		}
		else if (strcmp(sds_name[i],"cloud_top_pressure_low") == 0) {
			free(imgr2->cldp_low);
		}
		else if (strcmp(sds_name[i],"cloud_top_height") == 0) {
			free(imgr2->cldz);
		}
		// Added RAF
		else if (strcmp(sds_name[i],"clear_bt31") == 0) {
			free(imgr2->clrbt31);
		}
		// Added RAF
		else if (strcmp(sds_name[i],"surface_temperature") == 0) {
			free(imgr2->sfctmp);
		}
		// Added RAF
		else if (strcmp(sds_name[i],"cld_emiss11") == 0) {
			free(imgr2->cld_emiss11);
		}
		// Added RAF
		else if (strcmp(sds_name[i],"cld_emiss12") == 0) {
			free(imgr2->cld_emiss12);
		}
		// Added RAF
		else if (strcmp(sds_name[i],"cld_emiss13") == 0) {
			free(imgr2->cld_emiss13);
		}
		// Added RAF
		else if (strcmp(sds_name[i],"cld_emiss85") == 0) {
			free(imgr2->cld_emiss85);
		}
		// Added RAF
		else if (strcmp(sds_name[i],"modis_C6_IRP") == 0) {
			free(imgr2->modis_C6_IRP);
		}
		// Added RAF
		else if (strcmp(sds_name[i],"IRP_CTH_Consistency_Flag") == 0) {
			free(imgr2->IRP_CTH_Consistency_Flag);
		}
		// Added RAF (test)
		else if (strcmp(sds_name[i],"os_top_flag") == 0) {
			free(imgr2->os_top_flag);
		}
		else if (strcmp(sds_name[i],"aerosol_top_height") == 0) {
			free(imgr2->aeroz);
		}
		else if (strcmp(sds_name[i],"cloud_top_height_high") == 0) {
			free(imgr2->cldz_high);
		}
		else if (strcmp(sds_name[i],"cloud_top_height_low") == 0) {
			free(imgr2->cldz_low);
		}
		else if (strcmp(sds_name[i],"cloud_particle_effective_radius") == 0) {
			free(imgr2->cldreff);
		}
		else if (strcmp(sds_name[i],"cloud_particle_effective_radius_high") == 0) {
			free(imgr2->cldreff_high);
		}
		else if (strcmp(sds_name[i],"cloud_particle_effective_radius_low") == 0) {
			free(imgr2->cldreff_low);
		}
		else if (strcmp(sds_name[i],"cloud_particle_effective_diameter") == 0) {
			free(imgr2->clddeff);
		}
		else if (strcmp(sds_name[i],"cloud_particle_effective_diameter_high") == 0) {
			free(imgr2->clddeff_high);
		}
		else if (strcmp(sds_name[i],"cloud_particle_effective_diameter_low") == 0) {
			free(imgr2->clddeff_low);
		}
		else if (strcmp(sds_name[i],"cloud_emissivity") == 0) {
			free(imgr2->cldemiss);
		}
		else if (strcmp(sds_name[i],"aerosol_emissivity") == 0) {
			free(imgr2->aeroemiss);
		}
		else if (strcmp(sds_name[i],"cloud_emissivity_high") == 0) {
			free(imgr2->cldemiss_high);
		}
		else if (strcmp(sds_name[i],"cloud_emissivity_low") == 0) {
			free(imgr2->cldemiss_low);
		}
		else if (strcmp(sds_name[i],"cloud_optical_depth_ir") == 0) {
			free(imgr2->cod_ir);
		}
		else if (strcmp(sds_name[i],"cloud_optical_depth_ir_high") == 0) {
			free(imgr2->cod_ir_high);
		}
		else if (strcmp(sds_name[i],"cloud_optical_depth_ir_low") == 0) {
			free(imgr2->cod_ir_low);
		}
		else if (strcmp(sds_name[i],"cloud_optical_depth_vis") == 0) {
			free(imgr2->cod_vis);
		}
		else if (strcmp(sds_name[i],"cloud_optical_depth_vis_high") == 0) {
			free(imgr2->cod_vis_high);
		}
		else if (strcmp(sds_name[i],"cloud_optical_depth_vis_low") == 0) {
			free(imgr2->cod_vis_low);
		}
		else if (strcmp(sds_name[i],"cloud_bottom") == 0) {
			free(imgr2->cldbot);
		}
		else if (strcmp(sds_name[i],"aerosol_bottom") == 0) {
			free(imgr2->aerobot);
		}
		else if (strcmp(sds_name[i],"cloud_bottom_high") == 0) {
			free(imgr2->cldbot_high);
		}
		else if (strcmp(sds_name[i],"cloud_bottom_low") == 0) {
			free(imgr2->cldbot_low);
		}
		else if (strcmp(sds_name[i],"cloud_thickness") == 0) {
			free(imgr2->clddz);
		}
		else if (strcmp(sds_name[i],"aerosol_thickness") == 0) {
			free(imgr2->aerodz);
		}
		else if (strcmp(sds_name[i],"cloud_thickness_high") == 0) {
			free(imgr2->clddz_high);
		}
		else if (strcmp(sds_name[i],"cloud_thickness_low") == 0) {
			free(imgr2->clddz_low);
		}
		else if (strcmp(sds_name[i],"cloud_fraction") == 0) {
			free(imgr2->cldfrac);
		}
		else if (strcmp(sds_name[i],"aerosol_fraction") == 0) {
			free(imgr2->aerofrac);
		}
		else if (strcmp(sds_name[i],"cloud_liquid_water_path") == 0) {
			free(imgr2->cldlwp);
		}
		else if (strcmp(sds_name[i],"cloud_ice_water_path") == 0) {
			free(imgr2->cldiwp);
		}
		else if (strcmp(sds_name[i],"cloud_liquid_water_content") == 0) {
			free(imgr2->cldlwc);
		}
		else if (strcmp(sds_name[i],"cloud_ice_water_content") == 0) {
			free(imgr2->cldiwc);
		}
		else if (strcmp(sds_name[i],"cloud_beta") == 0) {
			free(imgr2->cldbeta);
		}
		else if (strcmp(sds_name[i],"cloud_beta_high") == 0) {
			free(imgr2->cldbeta_high);
		}
		else if (strcmp(sds_name[i],"cloud_beta_mid") == 0) {
			free(imgr2->cldbeta_mid);
		}
		else if (strcmp(sds_name[i],"cloud_beta_low") == 0) {
			free(imgr2->cldbeta_low);
		}
		else if (strcmp(sds_name[i],"emiss11_high") == 0) {
			free(imgr2->emiss11_high);
		}
		else if (strcmp(sds_name[i],"emiss11_mid") == 0) {
			free(imgr2->emiss11_mid);
		}
		else if (strcmp(sds_name[i],"emiss11_low") == 0) {
			free(imgr2->emiss11_low);
		}
		else if (strcmp(sds_name[i],"column_aerosol_amount") == 0) {
			free(imgr2->col_aer);
		}
		else if (strcmp(sds_name[i],"aerosol_optical_depth_vis") == 0) {
			free(imgr2->aod_vis);
		}
		else if (strcmp(sds_name[i],"aerosol_optical_depth_ir") == 0) {
			free(imgr2->aod_ir);
		}
		else if (strcmp(sds_name[i],"aerosol_particle_effective_radius") == 0) {
			free(imgr2->aerreff);
		}
		else if (strcmp(sds_name[i],"aerosol_particle_effective_diameter") == 0) {
			free(imgr2->aerdeff);
		}
	}
}
