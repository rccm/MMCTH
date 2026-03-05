/*$Id: create_level2_default_output.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "output_leocat.h"

void write_hdf (hdf_output, void *);

void create_level2_default_output(int32 sd_id, imagerL1 imgr1, char *nwp_source_name)
{

	
	printf("Writing level 2 default attributes to HDF...\n"); //GPC
  
	char l1bdata_type[MAX_STR_LEN];
	int32 status;
	float lat_region[2], lon_region[2];
	hdf_output hdf;

	if (imgr1.directbc_flg == YES)
	  strcpy(l1bdata_type,"Direct Broadcast");
	else
	  strcpy(l1bdata_type,"Granule");
	  
	lat_region[0] = imgr1.bounds[0], lat_region[1] = imgr1.bounds[2];
	lon_region[0] = imgr1.bounds[1], lon_region[1] = imgr1.bounds[3];

	status = SDsetattr(sd_id,"LEOCAT_VERSION",DFNT_CHAR8,strlen(PVERSION),(void *) &PVERSION[0]);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - LEOCAT_VERSION - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	// GPCstatus = SDsetattr(sd_id,"L1B_FILENAME",DFNT_CHAR8,MAX_STR_LEN,(void *) &imgr1.l1b_filename1[0]);
	status = SDsetattr(sd_id,"L1B_FILENAME",DFNT_CHAR8,strlen(&imgr1.l1b_filename1[0]),(void *) &imgr1.l1b_filename1[0]); //GPC
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - L1B_FILENAME - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	// GPC status = SDsetattr(sd_id,"SATNAME",DFNT_CHAR8,MAX_STR_LEN,(void *) &imgr1.satname[0]);
	status = SDsetattr(sd_id,"SATNAME",DFNT_CHAR8,strlen(&imgr1.satname[0]),(void *) &imgr1.satname[0]);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - SATNAME - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	// GPC status = SDsetattr(sd_id,"INTRUMENT_NAME",DFNT_CHAR8,MAX_STR_LEN,(void *) &imgr1.instrumentname[0]);
	status = SDsetattr(sd_id,"INTRUMENT_NAME",DFNT_CHAR8,strlen(&imgr1.instrumentname[0]),(void *) &imgr1.instrumentname[0]);//GPC
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - INSTRUMENT_NAME - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	status = SDsetattr(sd_id,"L1B_RESOLUTION",DFNT_FLOAT32,1,(void *) &imgr1.l1b_data_res);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - L1B_RESOLUTION - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	// GPC status = SDsetattr(sd_id,"L1B_RESOLUTION_UNITS",DFNT_CHAR8,MAX_STR_LEN,(void *) &imgr1.l1b_data_res_units[0]);
	status = SDsetattr(sd_id,"L1B_RESOLUTION_UNITS",DFNT_CHAR8,strlen(&imgr1.l1b_data_res_units[0]),(void *) &imgr1.l1b_data_res_units[0]); //GPC
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - L1B_RESOLUTION_UNITS - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	// GPC status = SDsetattr(sd_id,"SATNODE",DFNT_CHAR8,MAX_STR_LEN,(void *) &imgr1.asc_des_attr[0]);
	status = SDsetattr(sd_id,"SATNODE",DFNT_CHAR8,strlen(&imgr1.asc_des_attr[0]),(void *) &imgr1.asc_des_attr[0]); //GPC
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - SATNODE - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	// GPC status = SDsetattr(sd_id,"L1BDATA_TYPE",DFNT_CHAR8,MAX_STR_LEN,(void *) &l1bdata_type[0]);
	status = SDsetattr(sd_id,"L1BDATA_TYPE",DFNT_CHAR8,strlen(&l1bdata_type[0]),(void *) &l1bdata_type[0]); //GPC
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - L1BDATA_TYPE - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	status = SDsetattr(sd_id,"YEAR",DFNT_INT16,1,(void *) &imgr1.year);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - YEAR - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	status = SDsetattr(sd_id,"DAY",DFNT_INT16,1,(void *) &imgr1.jday);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - DAY - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	status = SDsetattr(sd_id,"TIME",DFNT_FLOAT,1,(void *) &imgr1.time);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - TIME - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	status = SDsetattr(sd_id,"EARTH-SUN_DISTANCE",DFNT_FLOAT32,1,(void *) &imgr1.earth2sun_dist);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - EARTH-SUN_DISTANCE - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	status = SDsetattr(sd_id,"LAT_REGION",DFNT_FLOAT32,2,(void *) lat_region);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - LAT_REGION - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	status = SDsetattr(sd_id,"LON_REGION",DFNT_FLOAT32,2,(void *) lon_region);
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - LON_REGION - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}
	
	// GPC status = SDsetattr(sd_id,"NWP_SOURCE",DFNT_CHAR8,MAX_STR_LEN,(void *) &nwp_source_name[0]);
	status = SDsetattr(sd_id,"NWP_SOURCE",DFNT_CHAR8,strlen(&nwp_source_name[0]),(void *) &nwp_source_name[0]); //GPC
	if (status == FAIL) {
	  fprintf(stderr,"%sCannot create L1b global attribute - NWP_SOURCE - aborting\n",EXE_PROMPT);
	  exit(EXIT_FAILURE);
	}


	hdf.sd_id = sd_id;
	strcpy(hdf.reference,"NA");
	hdf.rank = 2;
	hdf.dimen[0] = imgr1.nrow;
	hdf.dimen[1] = imgr1.ncol;
	strcpy(hdf.dim_name[0],"NPIXELS_ALONG_TRACK");
	strcpy(hdf.dim_name[1],"NPIXELS_ACROSS_TRACK");

	strcpy(hdf.sds_name,"Latitude");
	hdf.range_min = range_min_lat;
	hdf.range_max = range_max_lat;
	strcpy(hdf.units,units_lat);
	hdf.scaled_flg = scaled_flg_lat;
	hdf.scaled_type = scaled_type_lat_pixel;
	write_hdf(hdf, (void *) imgr1.lat);

	strcpy(hdf.sds_name,"Longitude");
	hdf.range_min = range_min_lon;
	hdf.range_max = range_max_lon;
	strcpy(hdf.units,units_lon);
	hdf.scaled_flg = scaled_flg_lon;
	hdf.scaled_type = scaled_type_lon_pixel;
	write_hdf(hdf, (void *) imgr1.lon);

	strcpy(hdf.sds_name,"Solar_Zenith_Angle");
	hdf.range_min = range_min_solzen;
	hdf.range_max = range_max_solzen;
	strcpy(hdf.units,units_solzen);
	hdf.scaled_flg = scaled_flg_solzen;
	hdf.scaled_type = scaled_type_solzen_pixel;
	write_hdf(hdf, (void *) imgr1.solzen);

	strcpy(hdf.sds_name,"Satellite_Zenith_Angle");
	hdf.range_min = range_min_satzen;
	hdf.range_max = range_max_satzen;
	strcpy(hdf.units,units_satzen);
	hdf.scaled_flg = scaled_flg_satzen;
	hdf.scaled_type = scaled_type_satzen_pixel;
	write_hdf(hdf, (void *) imgr1.satzen);

	strcpy(hdf.sds_name,"Relative_Azimuth");
	hdf.range_min = range_min_relaz;
	hdf.range_max = range_max_relaz;
	strcpy(hdf.units,units_relaz);
	hdf.scaled_flg = scaled_flg_relaz;
	hdf.scaled_type = scaled_type_relaz_pixel;
	write_hdf(hdf, (void *) imgr1.relaz);

        strcpy(hdf.sds_name,"IGBP_Surface_Type");
        hdf.range_min = range_min_sfc_type;
        hdf.range_max = range_max_sfc_type;
        strcpy(hdf.units,units_sfc_type);
        hdf.scaled_flg = scaled_flg_sfc_type;
        hdf.scaled_type = scaled_type_sfc_type_pixel;
        write_hdf(hdf, (void *) imgr1.sfc_type);

}
