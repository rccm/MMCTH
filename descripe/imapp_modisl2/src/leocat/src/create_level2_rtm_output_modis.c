/*$Id: create_level2_rtm_output_modis.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <hdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "nwp_leocat.h"
#include "rtm_leocat.h"
#include "output_leocat.h"

void write_hdf (hdf_output, void *);

void create_level2_rtm_output_modis(int32 sd_id, imagerL1 imgr1, nwp_params nwp,
                                    rtm_profiles rtm, int8 *chflg_modis)
{
   
  /*hdf.sd_id = sd_id;
  strcpy(hdf.reference,"NA");
  hdf.rank = 2;
  hdf.dimen[0] = imgr1.nrow;
  hdf.dimen[1] = imgr1.ncol;
  strcpy(hdf.dim_name[0],"NPIXELS_ALONG_TRACK");
  strcpy(hdf.dim_name[1],"NPIXELS_ACROSS_TRACK");
  
  if ((pix_temp = (float *) malloc(imgr1.npts*sizeof(float))) == NULL)
    error_allo(rout,"pix_temp");
    
  hdf.range_min = range_min_bt;
  hdf.range_max = range_max_bt;
  strcpy(hdf.units,units_bt);
  hdf.scaled_flg = scaled_flg_bt;
  hdf.scaled_type = scaled_type_bt_pixel;  
    
  if (chflg_modis[19]) {
    strcpy(hdf.sds_name,"BT20_Clear");
    write_hdf(hdf, (void *) rtm.bt20_clr);
  }
  
  if (chflg_modis[20]) {
    strcpy(hdf.sds_name,"BT21_Clear");
    write_hdf(hdf, (void *) rtm.bt21_clr);
  }
  
  if (chflg_modis[21]) {
    strcpy(hdf.sds_name,"BT22_Clear");
    write_hdf(hdf, (void *) rtm.bt22_clr);
  }
  
  if (chflg_modis[22]) {
    strcpy(hdf.sds_name,"BT23_Clear");
    write_hdf(hdf, (void *) rtm.bt23_clr);
  }
  
  if (chflg_modis[23]) {
    strcpy(hdf.sds_name,"BT24_Clear");
    write_hdf(hdf, (void *) rtm.bt24_clr);
  }
  
  if (chflg_modis[24]) {
    strcpy(hdf.sds_name,"BT25_Clear");
    write_hdf(hdf, (void *) rtm.bt25_clr);
  }
  
  if (chflg_modis[26]) {
    strcpy(hdf.sds_name,"BT27_Clear");
    write_hdf(hdf, (void *) rtm.bt27_clr);
  }
  
  if (chflg_modis[27]) {
    strcpy(hdf.sds_name,"BT28_Clear");
    write_hdf(hdf, (void *) rtm.bt28_clr);
  }
  
  if (chflg_modis[28]) {
    strcpy(hdf.sds_name,"BT29_Clear");
    write_hdf(hdf, (void *) rtm.bt29_clr);
  }
  
  if (chflg_modis[29]) {
    strcpy(hdf.sds_name,"BT30_Clear");
    write_hdf(hdf, (void *) rtm.bt30_clr);
  }
  
  if (chflg_modis[30]) {
    strcpy(hdf.sds_name,"BT31_Clear");
    write_hdf(hdf, (void *) rtm.bt31_clr);
  }
  
  if (chflg_modis[31]) {
    strcpy(hdf.sds_name,"BT32_Clear");
    write_hdf(hdf, (void *) rtm.bt32_clr);
  }
  
  if (chflg_modis[32]) {
    strcpy(hdf.sds_name,"BT33_Clear");
    write_hdf(hdf, (void *) rtm.bt33_clr);
  }
  
  if (chflg_modis[33]) {
    strcpy(hdf.sds_name,"BT34_Clear");
    write_hdf(hdf, (void *) rtm.bt34_clr);
  }
  
  if (chflg_modis[34]) {
    strcpy(hdf.sds_name,"BT35_Clear");
    write_hdf(hdf, (void *) rtm.bt35_clr);
  }
  
  if (chflg_modis[35]) {
    strcpy(hdf.sds_name,"BT36_Clear");
    write_hdf(hdf, (void *) rtm.bt36_clr);
  }
  
  free(pix_temp);*/
}
