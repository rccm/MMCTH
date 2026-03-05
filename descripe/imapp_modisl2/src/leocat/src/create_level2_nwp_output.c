/*$Id: create_level2_nwp_output.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <hdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "nwp_leocat.h"
#include "rtm_leocat.h"
#include "output_leocat.h"

void write_hdf (hdf_output, void *);

void create_level2_nwp_output(int32 sd_id, int8 nwp_source, nwp_params nwp, imagerL1 imgr1)
{
  char *rout = {"create_level2_nwp_output.c"};
  char tsurf_units[MAX_STR_LEN], tpw_units[MAX_STR_LEN],
       snow_units[MAX_STR_LEN], tropo_t_units[MAX_STR_LEN];

  hdf_output hdf;
  int i;
  float minval, maxval;
  float *pix_temp;
  
  hdf.sd_id = sd_id;
  strcpy(hdf.reference,"NA");
  hdf.rank = 2;
  hdf.dimen[0] = imgr1.nrow;
  hdf.dimen[1] = imgr1.ncol;
  strcpy(hdf.dim_name[0],"NPIXELS_ALONG_TRACK");
  strcpy(hdf.dim_name[1],"NPIXELS_ACROSS_TRACK");
  
  if ((pix_temp = (float *) malloc(imgr1.npts*sizeof(float))) == NULL)
    error_allo(rout,"pix_temp");  
    
  switch (nwp_source) {
    case (GFS_NWP):
      strcpy(tsurf_units, "K");
      strcpy(tpw_units, "Kg/m^2");
      strcpy(snow_units,"Kg/m^2");
      strcpy(tropo_t_units,"K");
      break;
    case (NCEP_NWP):
      strcpy(tsurf_units, "K");
      strcpy(tpw_units, "Kg/m^2");
      strcpy(snow_units,"Kg/m^2");
      strcpy(tropo_t_units,"K");
      break;
    case (ECMWF_NWP):
      strcpy(tsurf_units, "K");
      strcpy(tpw_units, "Kg/m^2");
      strcpy(snow_units,"Kg/m^2");
      strcpy(tropo_t_units,"K");
      break;
    default:
      fprintf(stderr,"Specified NWP type in not supported - see help\n");
      exit(EXIT_FAILURE);
      break;
  }
  
  strcpy(hdf.sds_name,"Surface_Temperature_NWP");
  hdf.range_min = 190.0;
  hdf.range_max = 360.0;
  strcpy(hdf.units,tsurf_units);
  hdf.scaled_flg = linear_scale;
  hdf.scaled_type = DFNT_INT16;
  for (i=0; i<imgr1.npts; i++) pix_temp[i] = nwp.map[imgr1.index_nwp[i]].t_sfc;
  array_minmax_float(imgr1.npts, pix_temp, &minval, &maxval);
  write_hdf(hdf, (void *) pix_temp);
  
  strcpy(hdf.sds_name,"TPW_NWP");
  hdf.range_min = 0.0;
  hdf.range_max = 10.0;
  strcpy(hdf.units,tpw_units);
  hdf.scaled_flg = linear_scale;
  hdf.scaled_type = DFNT_INT8;
  for (i=0; i<imgr1.npts; i++) pix_temp[i] = nwp.map[imgr1.index_nwp[i]].tpw;
  array_minmax_float(imgr1.npts, pix_temp, &minval, &maxval);
  write_hdf(hdf, (void *) pix_temp);
  
  strcpy(hdf.sds_name,"Snow_Depth_NWP");
  hdf.range_min = 0.0;
  hdf.range_max = 10.0;
  strcpy(hdf.units,snow_units);
  hdf.scaled_flg = linear_scale;
  hdf.scaled_type = DFNT_INT8;
  for (i=0; i<imgr1.npts; i++) pix_temp[i] = nwp.map[imgr1.index_nwp[i]].snow_sfc;
  array_minmax_float(imgr1.npts, pix_temp, &minval, &maxval);
  write_hdf(hdf, (void *) pix_temp);
  
  strcpy(hdf.sds_name,"Tropopause_Temperature_NWP");
  hdf.range_min = 180.0;
  hdf.range_max = 250.0;
  strcpy(hdf.units,tropo_t_units);
  hdf.scaled_flg = linear_scale;
  hdf.scaled_type = DFNT_INT8;
  for (i=0; i<imgr1.npts; i++) pix_temp[i] = nwp.map[imgr1.index_nwp[i]].t_tropo;
  array_minmax_float(imgr1.npts, pix_temp, &minval, &maxval);
  write_hdf(hdf, (void *) pix_temp);
  
  free(pix_temp);
}
