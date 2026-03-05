/*$Id: create_level2_output.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "imagerL2_leocat.h"
#include "output_leocat.h"

#define numPSA 18

void write_hdf(hdf_output, void *);
void write_hdf_multidim(hdf_output, int32 *, int32 *, int32 *, void *);

void create_level2_output(int32 sd_id, int16 algo_num, int16 nout,
		char *keyword, char *reference, char *sds_name[MAX_ALGO_OUTPUT],
		imagerL2 imgr2) {
	int i, n, status, x, y, iout;
	int32 start[3], edge[3], sds_id;
	hdf_output hdf;
        char PSAName[numPSA][30] = {
          "SuccessfulRetrievalPct",
          "VeryHighConfidentClearPct",
          "HighConfidentClearPct",
          "UncertainConfidentClearPct",
          "LowConfidentClearPct",
          "CloudCoverPct250m",
          "ClearPct250m",
          "DayProcessedPct",
          "NightProcessedPct",
          "SunglintProcessedPct",
          "Snow_IceSurfaceProcessedPct",
          "LandProcessedPct",
          "WaterProcessedPct",
          "ThinCirrusSolarFoundPct",
          "ThinCirrusIR_FoundPct",
          "NonCloudObstructionFoundPct",
          "MinSolarZenithAngle",
          "MaxSolarZenithAngle"
        };

	hdf.sd_id = sd_id;
	hdf.algo_num = algo_num;
	strcpy(hdf.reference, reference);
	hdf.rank = 2;
	hdf.dimen[0] = imgr2.nrow;
	hdf.dimen[1] = imgr2.ncol;
	strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
	strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");

        if(imgr2.output_cm_stats_flag == 1) {

//        Write out the statistics   K. Baggett 
          printf("\nWriting out the statistics to global attributes section\n");

          for (i=0; i<numPSA; i++) {
            printf("parameter=%s,val=%7.2f\n",PSAName[i],imgr2.stats_arr[i]);
            status = SDsetattr(sd_id,PSAName[i],DFNT_FLOAT,1,(void *) &imgr2.stats_arr[i]);
            if (status == FAIL) {
              fprintf(stderr,"%sCannot create L1b global attribute - %s - aborting\n",EXE_PROMPT,PSAName[i]);
              exit(EXIT_FAILURE);
            }
          }

//        Write out the NDVI file name as a global attribute
          status = SDsetattr(sd_id,"NDVIFile",DFNT_CHAR8,strlen(&imgr2.ndviFile[0]),(void *) &imgr2.ndviFile[0]);
          printf("NDVIfile = %s\n\n",imgr2.ndviFile);
//        Write out the thresholds file name as a global attribute
          status=SDsetattr(sd_id,"ThresholdFile",DFNT_CHAR8,strlen(&imgr2.threshFile[0]),(void *)&imgr2.threshFile[0]);
          printf("Thesholdfile = %s\n\n",imgr2.threshFile);
  
        }

//RAF   Added this because LEOCAT inflexible as to number of output SDSs - cannot change on the fly.  There will
//      be no output 250-m stats if no 250-m data processed (imgr2.output_250).  Night MODIS data has no 250-m
//      data. So when no 250-m output, hdf.sds_name contains one more item than available data outputs.  See
//      "250m_stats" below.
        if(algo_num == 29 && imgr2.output_250 == NO) {
          iout = nout - 1;
        }
        else {
          iout = nout;
        }
//      printf("output: %d %d %d %d \n", imgr2.output_250, NO, iout, nout);
//      (void)fflush(stdout);

	for (i=0; i<iout; i++) {

		sprintf(hdf.sds_name, "%s_%s", keyword, sds_name[i]);
		printf("hdf.sds_name[%d]: %s \n", i, hdf.sds_name);

		if (strcmp(sds_name[i], "cloud_mask") == 0) {
			printf("Writing Cloud Mask to hdf...\n"); //GPC
			hdf.range_min = range_min_cldmask;
			hdf.range_max = range_max_cldmask;
			strcpy(hdf.units, units_cldmask);
			hdf.scaled_flg = scaled_flg_cldmask;
			hdf.scaled_type = scaled_type_cldmask_pixel;
			write_hdf(hdf, (void *) imgr2.cldmask);
			printf("\nFinished writing imgr2.cldmask to HDF\n");
		} else if (strcmp(sds_name[i], "land_mask") == 0) {
			hdf.range_min = range_min_landmask;
			hdf.range_max = range_max_landmask;
			strcpy(hdf.units, units_landmask);
			hdf.scaled_flg = scaled_flg_landmask;
			hdf.scaled_type = scaled_type_landmask_pixel;
			write_hdf(hdf, (void *) imgr2.landmask);
		} else if (strcmp(sds_name[i], "cloud_type") == 0) {
			hdf.range_min = range_min_cldtype;
			hdf.range_max = range_max_cldtype;
			strcpy(hdf.units, units_cldtype);
			hdf.scaled_flg = scaled_flg_cldtype;
			hdf.scaled_type = scaled_type_cldtype_pixel;
			write_hdf(hdf, (void *) imgr2.cldtype);
		} else if (strcmp(sds_name[i], "cloud_phase") == 0) {
			printf("Writing Cloud Phase to hdf...\n"); //GPC
			///////// GPC
			hdf.rank = 2;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			///////// GPC
			hdf.range_min = range_min_cldphase;
			hdf.range_max = range_max_cldphase;
			strcpy(hdf.units, units_cldphase);
			hdf.scaled_flg = scaled_flg_cldphase;
			hdf.scaled_type = scaled_type_cldphase_pixel;
			write_hdf(hdf, (void *) imgr2.cldphase);
			printf("\nFinished writing imgr2.cldphase to HDF\n");
		} else if (strcmp(sds_name[i], "aerosol_mask") == 0) {
			hdf.range_min = range_min_aeromask;
			hdf.range_max = range_max_aeromask;
			strcpy(hdf.units, units_aeromask);
			hdf.scaled_flg = scaled_flg_aeromask;
			hdf.scaled_type = scaled_type_aeromask_pixel;
			write_hdf(hdf, (void *) imgr2.aeromask);
		} else if (strcmp(sds_name[i], "so2mask") == 0) {
			hdf.range_min = range_min_so2mask;
			hdf.range_max = range_max_so2mask;
			strcpy(hdf.units, units_so2mask);
			hdf.scaled_flg = scaled_flg_so2mask;
			hdf.scaled_type = scaled_type_so2mask_pixel;
			write_hdf(hdf, (void *) imgr2.so2mask);
		} else if (strcmp(sds_name[i], "cop_qf0") == 0) {        // GPC
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");
			hdf.range_min = range_min_byteArray;
			hdf.range_max = range_max_byteArray;
			strcpy(hdf.units, units_byteArray);
			hdf.scaled_flg = scaled_flg_byteArray;
			hdf.scaled_type = scaled_type_byteArray_pixel;
			write_hdf(hdf, (void *) imgr2.cop_qf0);
		} else if (strcmp(sds_name[i], "cop_qf1") == 0) {        // GPC
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");
			hdf.range_min = range_min_byteArray;
			hdf.range_max = range_max_byteArray;
			strcpy(hdf.units, units_byteArray);
			hdf.scaled_flg = scaled_flg_byteArray;
			hdf.scaled_type = scaled_type_byteArray_pixel;
			write_hdf(hdf, (void *) imgr2.cop_qf1);
		} else if (strcmp(sds_name[i], "inwctt_qf") == 0) {         // GPC
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");
			hdf.range_min = range_min_byteArray;
			hdf.range_max = range_max_byteArray;
			strcpy(hdf.units, units_byteArray);
			hdf.scaled_flg = scaled_flg_byteArray;
			hdf.scaled_type = scaled_type_byteArray_pixel;
			write_hdf(hdf, (void *) imgr2.inwctt_qf);
		} else if (strcmp(sds_name[i], "flag_arr1") == 0) {

			hdf.rank = 3;
			//if (imgr2.numberFlagsArray1Buffers == 0)
				imgr2.numberFlagsArray1Buffers = 6;
			hdf.dimen[0] = imgr2.numberFlagsArray1Buffers;
			hdf.dimen[1] = imgr2.nrow;
			hdf.dimen[2] = imgr2.ncol;

			sprintf(hdf.dim_name[0],"%s_%s",keyword, "FLAG1_BYTE_SEGMENT");
			strcpy(hdf.dim_name[1], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[2], "NPIXELS_ACROSS_TRACK");
			hdf.range_min = range_min_cldmask;
			hdf.range_max = range_max_cldmask;
			strcpy(hdf.units, units_cldmask);
			hdf.scaled_flg = scaled_flg_cldmask;
			hdf.scaled_type = scaled_type_cldmask_pixel;
			start[1] = 0;
			start[2] = 0;
			edge[0] = 1;
			edge[1] = imgr2.nrow;
			edge[2] = imgr2.ncol;
			sds_id = -999.0;

//	                printf("Uninitialised %s: %d %d %s \n",sds_name[i], sds_id, hdf.dimen[0],
//			             hdf.dim_name[0]); //GPC

			for (n=0; n<imgr2.numberFlagsArray1Buffers; n++) {

				if (imgr2.flag_arr1[n] != NULL) {
					start[0] = n;

					write_hdf_multidim(hdf, start, edge, &sds_id,
							(void *) imgr2.flag_arr1[n]);
				}
			}
			status = SDendaccess(sds_id);
			if (status == FAIL) {
				fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
				exit(EXIT_FAILURE);
			}
			hdf.rank = 2;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");

//                      printf("Initialised %s: %d %d %s \n",sds_name[i], sds_id, hdf.dimen[0],
//                              hdf.dim_name[0]); //GPC

		} else if (strcmp(sds_name[i], "flag_arr2") == 0) {

			hdf.rank = 3;
			imgr2.numberFlagsArray2Buffers = 10;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			hdf.dimen[2] = imgr2.numberFlagsArray2Buffers;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");
			strcpy(hdf.dim_name[2], "FLAG2_BYTE_SEGMENT");
			hdf.range_min = range_min_cldmask;
			hdf.range_max = range_max_cldmask;
			strcpy(hdf.units, units_cldmask);
			hdf.scaled_flg = scaled_flg_cldmask;
			hdf.scaled_type = scaled_type_cldmask_pixel;
			start[2] = 0;
			edge[0] = 1;
			edge[1] = 1;
			edge[2] = imgr2.numberFlagsArray2Buffers;
			sds_id = -999.0;

//                      printf("Uninitialised %s: %d %d %s \n",sds_name[i], sds_id, hdf.dimen[0],
//                                  hdf.dim_name[0]); //GPC

			for (n=0; n<imgr2.npts; n++) {
				i2xy(n, imgr2.ncol, &x, &y);
				start[0] = y;
				start[1] = x;
				write_hdf_multidim(hdf, start, edge, &sds_id,
						(void *) imgr2.flag_arr2[n]);
			}

			status = SDendaccess(sds_id);
			if (status == FAIL) {
				fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
				exit(EXIT_FAILURE);
			}
			hdf.rank = 2;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");

//                      printf("Initialised %s: %d %d %s \n",sds_name[i], sds_id, hdf.dimen[0],
//                                   hdf.dim_name[0]); //GPC

		} else if (strcmp(sds_name[i], "stats_250m") == 0 && imgr2.output_250 == YES) {

			hdf.rank = 3;
			imgr2.numberStats250Buffers = 4;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			hdf.dimen[2] = imgr2.numberStats250Buffers;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");
			strcpy(hdf.dim_name[2], "NUM_250M_STATS");
			hdf.range_min = range_min_stats250;
			hdf.range_max = range_max_stats250;
			strcpy(hdf.units, units_stats250);
			hdf.scaled_flg = scaled_flg_stats250;
			hdf.scaled_type = scaled_type_stats250_pixel;
			start[2] = 0;
			edge[0] = 1;
			edge[1] = 1;
			edge[2] = imgr2.numberStats250Buffers;
			sds_id = -999.0;

//                      printf("Uninitialised %s: %d %d %s \n",sds_name[i], sds_id, hdf.dimen[0],
//                                  hdf.dim_name[0]);

			for (n=0; n<imgr2.npts; n++) {
				i2xy(n, imgr2.ncol, &x, &y);
				start[0] = y;
				start[1] = x;
				write_hdf_multidim(hdf, start, edge, &sds_id,
						(void *) imgr2.stats_250m[n]);
			}

			status = SDendaccess(sds_id);
			if (status == FAIL) {
				fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
				exit(EXIT_FAILURE);
			}

			hdf.rank = 2;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");

//                      printf("Initialised %s: %d %d %s \n",sds_name[i], sds_id, hdf.dimen[0],
//                                    hdf.dim_name[0]);

		} else if (strcmp(sds_name[i], "sst") == 0) {
			hdf.range_min = range_min_sst;
			hdf.range_max = range_max_sst;
			strcpy(hdf.units, units_sst);
			hdf.scaled_flg = scaled_flg_sst;
			hdf.scaled_type = scaled_type_sst_pixel;
			write_hdf(hdf, (void *) imgr2.sst);
		} else if (strcmp(sds_name[i], "lst") == 0) {
			hdf.range_min = range_min_lst;
			hdf.range_max = range_max_lst;
			strcpy(hdf.units, units_lst);
			hdf.scaled_flg = scaled_flg_lst;
			hdf.scaled_type = scaled_type_lst_pixel;
			write_hdf(hdf, (void *) imgr2.lst);
		} else if (strcmp(sds_name[i], "ist") == 0) {
			hdf.range_min = range_min_ist;
			hdf.range_max = range_max_ist;
			strcpy(hdf.units, units_ist);
			hdf.scaled_flg = scaled_flg_ist;
			hdf.scaled_type = scaled_type_ist_pixel;
			write_hdf(hdf, (void *) imgr2.ist);
		} else if (strcmp(sds_name[i], "ndvi") == 0) {
			hdf.range_min = range_min_ndvi;
			hdf.range_max = range_max_ndvi;
			strcpy(hdf.units, units_ndvi);
			hdf.scaled_flg = scaled_flg_ndvi;
			hdf.scaled_type = scaled_type_ndvi_pixel;
			write_hdf(hdf, (void *) imgr2.ndvi);
		} else if (strcmp(sds_name[i], "tpw") == 0) {
			hdf.range_min = range_min_tpw;
			hdf.range_max = range_max_tpw;
			strcpy(hdf.units, units_tpw);
			hdf.scaled_flg = scaled_flg_tpw;
			hdf.scaled_type = scaled_type_tpw_pixel;
			write_hdf(hdf, (void *) imgr2.tpw);
		} else if (strcmp(sds_name[i], "cloud_top_temperature") == 0) {
			hdf.range_min = range_min_cldt;
			hdf.range_max = range_max_cldt;
			strcpy(hdf.units, units_cldt);
			hdf.scaled_flg = scaled_flg_cldt;
			hdf.scaled_type = scaled_type_cldt_pixel;
			write_hdf(hdf, (void *) imgr2.cldt);

		} else if (strcmp(sds_name[i], "cloudTopTemp_Ice_NightWater") == 0) {
			printf("Writing Cloud Top Temperature (Ice/NightWater) to hdf...\n"); //GPC
			///////// GPC
			hdf.rank = 2;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			///////// GPC
			hdf.range_min = range_min_cldt;
			hdf.range_max = range_max_cldt;
			strcpy(hdf.units, units_cldt);
			hdf.scaled_flg = scaled_flg_cldt;
			hdf.scaled_type = scaled_type_cldt_pixel;
			printf("hdf.scaled_flg = %d\n",hdf.scaled_flg);	// GPC
			printf("hdf.scaled_type = %d\n",hdf.scaled_type);	// GPC
			write_hdf(hdf, (void *) imgr2.cldt_IceNW);
			printf("\nFinished writing imgr2.cldt_IceNW to HDF\n");

		} else if (strcmp(sds_name[i], "ctParmQf") == 0) {		// GPC

			hdf.rank = 3;
//			imgr2.numberFlagsArray1Buffers = 3;
//			hdf.dimen[0] = imgr2.numberFlagsArray1Buffers;
			int number_ctParmQfArrayBuffers = 3;
			hdf.dimen[0] = number_ctParmQfArrayBuffers;
			hdf.dimen[1] = imgr2.nrow;
			hdf.dimen[2] = imgr2.ncol;

			sprintf(hdf.dim_name[0],"%s_%s",keyword, "FLAG1_BYTE_SEGMENT");
			strcpy(hdf.dim_name[1], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[2], "NPIXELS_ACROSS_TRACK");
			hdf.range_min = range_min_cldmask;
			hdf.range_max = range_max_cldmask;
			strcpy(hdf.units, units_cldmask);
			hdf.scaled_flg = scaled_flg_cldmask;
			hdf.scaled_type = scaled_type_cldmask_pixel;
			start[1] = 0;
			start[2] = 0;
			edge[0] = 1;
			edge[1] = imgr2.nrow;
			edge[2] = imgr2.ncol;
			sds_id = -999.0;

			printf("Uninitialised %s: %d %d %s \n",sds_name[i], sds_id, hdf.dimen[0],
					hdf.dim_name[0]); //GPC

			for (n=0; n < number_ctParmQfArrayBuffers; n++) {

				if (imgr2.ctParmQf[n] != NULL) {
					start[0] = n;

					write_hdf_multidim(hdf, start, edge, &sds_id,
							(void *) imgr2.ctParmQf[n]);
				}
			}
			status = SDendaccess(sds_id);
			if (status == FAIL) {
				fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
				exit(EXIT_FAILURE);
			}
			hdf.rank = 2;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			strcpy(hdf.dim_name[0], "NPIXELS_ALONG_TRACK");
			strcpy(hdf.dim_name[1], "NPIXELS_ACROSS_TRACK");

			printf("Initialised %s: %d %d %s \n",sds_name[i], sds_id, hdf.dimen[0],
					hdf.dim_name[0]); //GPC

		} else if (strcmp(sds_name[i], "aerosol_top_temperature") == 0) {
			hdf.range_min = range_min_cldt;
			hdf.range_max = range_max_cldt;
			strcpy(hdf.units, units_cldt);
			hdf.scaled_flg = scaled_flg_cldt;
			hdf.scaled_type = scaled_type_cldt_pixel;
			write_hdf(hdf, (void *) imgr2.aerot);
		} else if (strcmp(sds_name[i], "cloud_top_temperature_high") == 0) {
			hdf.range_min = range_min_cldt;
			hdf.range_max = range_max_cldt;
			strcpy(hdf.units, units_cldt);
			hdf.scaled_flg = scaled_flg_cldt;
			hdf.scaled_type = scaled_type_cldt_pixel;
			write_hdf(hdf, (void *) imgr2.cldt_high);
		} else if (strcmp(sds_name[i], "cloud_top_temperature_low") == 0) {
			hdf.range_min = range_min_cldt;
			hdf.range_max = range_max_cldt;
			strcpy(hdf.units, units_cldt);
			hdf.scaled_flg = scaled_flg_cldt;
			hdf.scaled_type = scaled_type_cldt_pixel;
			write_hdf(hdf, (void *) imgr2.cldt_low);
		} else if (strcmp(sds_name[i], "cloud_top_pressure") == 0) {
			hdf.range_min = range_min_cldp;
			hdf.range_max = range_max_cldp;
			strcpy(hdf.units, units_cldp);
			hdf.scaled_flg = scaled_flg_cldp;
			hdf.scaled_type = scaled_type_cldp_pixel;
			write_hdf(hdf, (void *) imgr2.cldp);
		} else if (strcmp(sds_name[i], "aerosol_top_pressure") == 0) {
			hdf.range_min = range_min_cldp;
			hdf.range_max = range_max_cldp;
			strcpy(hdf.units, units_cldp);
			hdf.scaled_flg = scaled_flg_cldp;
			hdf.scaled_type = scaled_type_cldp_pixel;
			write_hdf(hdf, (void *) imgr2.aerop);
		} else if (strcmp(sds_name[i], "cloud_top_method") == 0) {
			hdf.range_min = range_min_solm;
			hdf.range_max = range_max_solm;
			strcpy(hdf.units, units_solm);
			hdf.scaled_flg = scaled_flg_solm;
			hdf.scaled_type = scaled_type_solm_pixel;
			write_hdf(hdf, (void *) imgr2.cld_solm);
		} else if (strcmp(sds_name[i], "aerosol_top_method") == 0) {
			hdf.range_min = range_min_solm;
			hdf.range_max = range_max_solm;
			strcpy(hdf.units, units_solm);
			hdf.scaled_flg = scaled_flg_solm;
			hdf.scaled_type = scaled_type_solm_pixel;
			write_hdf(hdf, (void *) imgr2.aero_solm);
		} else if (strcmp(sds_name[i], "cloud_top_pressure_high") == 0) {
			hdf.range_min = range_min_cldp;
			hdf.range_max = range_max_cldp;
			strcpy(hdf.units, units_cldp);
			hdf.scaled_flg = scaled_flg_cldp;
			hdf.scaled_type = scaled_type_cldp_pixel;
			write_hdf(hdf, (void *) imgr2.cldp_high);
		} else if (strcmp(sds_name[i], "cloud_top_pressure_low") == 0) {
			hdf.range_min = range_min_cldp;
			hdf.range_max = range_max_cldp;
			strcpy(hdf.units, units_cldp);
			hdf.scaled_flg = scaled_flg_cldp;
			hdf.scaled_type = scaled_type_cldp_pixel;
			write_hdf(hdf, (void *) imgr2.cldp_low);
		} else if (strcmp(sds_name[i], "cloud_top_height") == 0) {
			hdf.range_min = range_min_cldz;
			hdf.range_max = range_max_cldz;
			strcpy(hdf.units, units_cldz);
			hdf.scaled_flg = scaled_flg_cldz;
			hdf.scaled_type = scaled_type_cldz_pixel;
			write_hdf(hdf, (void *) imgr2.cldz);
// Added RAF
		} else if (strcmp(sds_name[i], "clear_bt31") == 0) {
			hdf.range_min = range_min_clrbt31;
			hdf.range_max = range_max_clrbt31;
			strcpy(hdf.units, units_clrbt31);
			hdf.scaled_flg = scaled_flg_clrbt31;
			hdf.scaled_type = scaled_type_clrbt31_pixel;
			write_hdf(hdf, (void *) imgr2.clrbt31);
// Added RAF
		} else if (strcmp(sds_name[i], "surface_temperature") == 0) {
			hdf.range_min = range_min_sfctmp;
			hdf.range_max = range_max_sfctmp;
			strcpy(hdf.units, units_sfctmp);
			hdf.scaled_flg = scaled_flg_sfctmp;
			hdf.scaled_type = scaled_type_sfctmp_pixel;
			write_hdf(hdf, (void *) imgr2.sfctmp);
// Added RAF
		} else if (strcmp(sds_name[i], "cld_emiss11") == 0) {
			hdf.range_min = range_min_cld_emiss11;
			hdf.range_max = range_max_cld_emiss11;
			strcpy(hdf.units, units_cld_emiss11);
			hdf.scaled_flg = scaled_flg_cld_emiss11;
			hdf.scaled_type = scaled_type_cld_emiss11_pixel;
			write_hdf(hdf, (void *) imgr2.cld_emiss11);
// Added RAF
		} else if (strcmp(sds_name[i], "cld_emiss12") == 0) {
			hdf.range_min = range_min_cld_emiss12;
			hdf.range_max = range_max_cld_emiss12;
			strcpy(hdf.units, units_cld_emiss12);
			hdf.scaled_flg = scaled_flg_cld_emiss12;
			hdf.scaled_type = scaled_type_cld_emiss12_pixel;
			write_hdf(hdf, (void *) imgr2.cld_emiss12);
// Added RAF
		} else if (strcmp(sds_name[i], "cld_emiss13") == 0) {
			hdf.range_min = range_min_cld_emiss13;
			hdf.range_max = range_max_cld_emiss13;
			strcpy(hdf.units, units_cld_emiss13);
			hdf.scaled_flg = scaled_flg_cld_emiss13;
			hdf.scaled_type = scaled_type_cld_emiss13_pixel;
			write_hdf(hdf, (void *) imgr2.cld_emiss13);
// Added RAF
		} else if (strcmp(sds_name[i], "cld_emiss85") == 0) {
			hdf.range_min = range_min_cld_emiss85;
			hdf.range_max = range_max_cld_emiss85;
			strcpy(hdf.units, units_cld_emiss85);
			hdf.scaled_flg = scaled_flg_cld_emiss85;
			hdf.scaled_type = scaled_type_cld_emiss85_pixel;
			write_hdf(hdf, (void *) imgr2.cld_emiss85);
// Added RAF
		} else if (strcmp(sds_name[i], "modis_C6_IRP") == 0) {
			hdf.range_min = range_min_modis_C6_IRP;
			hdf.range_max = range_max_modis_C6_IRP;
			strcpy(hdf.units, units_modis_C6_IRP);
			hdf.scaled_flg = scaled_flg_modis_C6_IRP;
			hdf.scaled_type = scaled_type_modis_C6_IRP_pixel;
			write_hdf(hdf, (void *) imgr2.modis_C6_IRP);
// Added RAF
		} else if (strcmp(sds_name[i], "IRP_CTH_Consistency_Flag") == 0) {
			hdf.range_min = range_min_IRP_CTH_Consistency_Flag;
			hdf.range_max = range_max_IRP_CTH_Consistency_Flag;
			strcpy(hdf.units, units_IRP_CTH_Consistency_Flag);
			hdf.scaled_flg = scaled_flg_IRP_CTH_Consistency_Flag;
			hdf.scaled_type = scaled_type_IRP_CTH_Consistency_Flag_pixel;
			write_hdf(hdf, (void *) imgr2.IRP_CTH_Consistency_Flag);
// Added RAF (test)
		} else if (strcmp(sds_name[i], "os_top_flag") == 0) {
			hdf.range_min = range_min_os_top_flag;
			hdf.range_max = range_max_os_top_flag;
			strcpy(hdf.units, units_os_top_flag);
			hdf.scaled_flg = scaled_flg_os_top_flag;
			hdf.scaled_type = scaled_type_os_top_flag_pixel;
			write_hdf(hdf, (void *) imgr2.os_top_flag);

		} else if (strcmp(sds_name[i], "aerosol_top_height") == 0) {
			hdf.range_min = range_min_cldz;
			hdf.range_max = range_max_cldz;
			strcpy(hdf.units, units_cldz);
			hdf.scaled_flg = scaled_flg_cldz;
			hdf.scaled_type = scaled_type_cldz_pixel;
			write_hdf(hdf, (void *) imgr2.aeroz);
		} else if (strcmp(sds_name[i], "cloud_top_height_high") == 0) {
			hdf.range_min = range_min_cldz;
			hdf.range_max = range_max_cldz;
			strcpy(hdf.units, units_cldz);
			hdf.scaled_flg = scaled_flg_cldz;
			hdf.scaled_type = scaled_type_cldz_pixel;
			write_hdf(hdf, (void *) imgr2.cldz_high);
		} else if (strcmp(sds_name[i], "cloud_top_height_low") == 0) {
			hdf.range_min = range_min_cldz;
			hdf.range_max = range_max_cldz;
			strcpy(hdf.units, units_cldz);
			hdf.scaled_flg = scaled_flg_cldz;
			hdf.scaled_type = scaled_type_cldz_pixel;
			write_hdf(hdf, (void *) imgr2.cldz_low);
		} else if (strcmp(sds_name[i], "cloud_particle_effective_radius") == 0) {
			printf("Writing Cloud Effective Particle Radius to hdf...\n"); //GPC
			///////// GPC
			hdf.rank = 2;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			///////// GPC
			hdf.range_min = range_min_cldreff;
			hdf.range_max = range_max_cldreff;
			strcpy(hdf.units, units_cldreff);
			hdf.scaled_flg = scaled_flg_cldreff;
			hdf.scaled_type = scaled_type_cldreff_pixel;
			printf("hdf.scaled_flg = %d\n",hdf.scaled_flg);	// GPC
			printf("hdf.scaled_type = %d\n",hdf.scaled_type);	// GPC
			write_hdf(hdf, (void *) imgr2.cldreff);
			printf("\nFinished writing imgr2.cldreff to HDF\n");
		} else if (strcmp(sds_name[i], "cloud_particle_effective_radius_high")
				== 0) {
			hdf.range_min = range_min_cldreff;
			hdf.range_max = range_max_cldreff;
			strcpy(hdf.units, units_cldreff);
			hdf.scaled_flg = scaled_flg_cldreff;
			hdf.scaled_type = scaled_type_cldreff_pixel;
			write_hdf(hdf, (void *) imgr2.cldreff_high);
		} else if (strcmp(sds_name[i], "cloud_particle_effective_radius_low")
				== 0) {
			hdf.range_min = range_min_cldreff;
			hdf.range_max = range_max_cldreff;
			strcpy(hdf.units, units_cldreff);
			hdf.scaled_flg = scaled_flg_cldreff;
			hdf.scaled_type = scaled_type_cldreff_pixel;
			write_hdf(hdf, (void *) imgr2.cldreff_low);
		} else if (strcmp(sds_name[i], "cloud_particle_effective_diameter")
				== 0) {
			hdf.range_min = range_min_clddeff;
			hdf.range_max = range_max_clddeff;
			strcpy(hdf.units, units_clddeff);
			hdf.scaled_flg = scaled_flg_clddeff;
			hdf.scaled_type = scaled_type_clddeff_pixel;
			write_hdf(hdf, (void *) imgr2.clddeff);
		} else if (strcmp(sds_name[i], "cloud_particle_effective_diameter_high")
				== 0) {
			hdf.range_min = range_min_clddeff;
			hdf.range_max = range_max_clddeff;
			strcpy(hdf.units, units_clddeff);
			hdf.scaled_flg = scaled_flg_clddeff;
			hdf.scaled_type = scaled_type_clddeff_pixel;
			write_hdf(hdf, (void *) imgr2.clddeff_high);
		} else if (strcmp(sds_name[i], "cloud_particle_effective_diameter_low")
				== 0) {
			hdf.range_min = range_min_clddeff;
			hdf.range_max = range_max_clddeff;
			strcpy(hdf.units, units_clddeff);
			hdf.scaled_flg = scaled_flg_clddeff;
			hdf.scaled_type = scaled_type_clddeff_pixel;
			write_hdf(hdf, (void *) imgr2.clddeff_low);
		} else if (strcmp(sds_name[i], "cloud_emissivity") == 0) {
			hdf.range_min = range_min_cldemiss;
			hdf.range_max = range_max_cldemiss;
			strcpy(hdf.units, units_cldemiss);
			hdf.scaled_flg = scaled_flg_cldemiss;
			hdf.scaled_type = scaled_type_cldemiss_pixel;
			write_hdf(hdf, (void *) imgr2.cldemiss);
		} else if (strcmp(sds_name[i], "aerosol_emissivity") == 0) {
			hdf.range_min = range_min_cldemiss;
			hdf.range_max = range_max_cldemiss;
			strcpy(hdf.units, units_cldemiss);
			hdf.scaled_flg = scaled_flg_cldemiss;
			hdf.scaled_type = scaled_type_cldemiss_pixel;
			write_hdf(hdf, (void *) imgr2.aeroemiss);
		} else if (strcmp(sds_name[i], "cloud_emissivity_high") == 0) {
			hdf.range_min = range_min_cldemiss;
			hdf.range_max = range_max_cldemiss;
			strcpy(hdf.units, units_cldemiss);
			hdf.scaled_flg = scaled_flg_cldemiss;
			hdf.scaled_type = scaled_type_cldemiss_pixel;
			write_hdf(hdf, (void *) imgr2.cldemiss_high);
		} else if (strcmp(sds_name[i], "cloud_emissivity_low") == 0) {
			hdf.range_min = range_min_cldemiss;
			hdf.range_max = range_max_cldemiss;
			strcpy(hdf.units, units_cldemiss);
			hdf.scaled_flg = scaled_flg_cldemiss;
			hdf.scaled_type = scaled_type_cldemiss_pixel;
			write_hdf(hdf, (void *) imgr2.cldemiss_low);
		} else if (strcmp(sds_name[i], "cloud_optical_depth_ir") == 0) {
			hdf.range_min = range_min_cod_ir;
			hdf.range_max = range_max_cod_ir;
			strcpy(hdf.units, units_cod_ir);
			hdf.scaled_flg = scaled_flg_cod_ir;
			hdf.scaled_type = scaled_type_cod_ir_pixel;
			write_hdf(hdf, (void *) imgr2.cod_ir);
			printf("\nFinished writing imgr2.cod_ir to HDF\n");
		} else if (strcmp(sds_name[i], "cloud_optical_depth_ir_high") == 0) {
			hdf.range_min = range_min_cod_ir;
			hdf.range_max = range_max_cod_ir;
			strcpy(hdf.units, units_cod_ir);
			hdf.scaled_flg = scaled_flg_cod_ir;
			hdf.scaled_type = scaled_type_cod_ir_pixel;
			write_hdf(hdf, (void *) imgr2.cod_ir_high);
		} else if (strcmp(sds_name[i], "cloud_optical_depth_ir_low") == 0) {
			hdf.range_min = range_min_cod_ir;
			hdf.range_max = range_max_cod_ir;
			strcpy(hdf.units, units_cod_ir);
			hdf.scaled_flg = scaled_flg_cod_ir;
			hdf.scaled_type = scaled_type_cod_ir_pixel;
			write_hdf(hdf, (void *) imgr2.cod_ir_low);
		} else if (strcmp(sds_name[i], "cloud_optical_depth_vis") == 0) {
			printf("Writing Cloud Optical Depth (visible) to hdf...\n"); //GPC
			///////// GPC
			hdf.rank = 2;
			hdf.dimen[0] = imgr2.nrow;
			hdf.dimen[1] = imgr2.ncol;
			///////// GPC
			hdf.range_min = range_min_cod_vis;
			hdf.range_max = range_max_cod_vis;
			strcpy(hdf.units, units_cod_vis);
			hdf.scaled_flg = scaled_flg_cod_vis;
			hdf.scaled_type = scaled_type_cod_vis_pixel;
			printf("hdf.scaled_flg = %d\n",hdf.scaled_flg);	// GPC
			printf("hdf.scaled_type = %d\n",hdf.scaled_type);	// GPC
			write_hdf(hdf, (void *) imgr2.cod_vis);
			printf("\nFinished writing imgr2.cod_vis to HDF\n");
		} else if (strcmp(sds_name[i], "cloud_optical_depth_vis_high") == 0) {
			hdf.range_min = range_min_cod_vis;
			hdf.range_max = range_max_cod_vis;
			strcpy(hdf.units, units_cod_vis);
			hdf.scaled_flg = scaled_flg_cod_vis;
			hdf.scaled_type = scaled_type_cod_vis_pixel;
			write_hdf(hdf, (void *) imgr2.cod_vis_high);
		} else if (strcmp(sds_name[i], "cloud_optical_depth_vis_low") == 0) {
			hdf.range_min = range_min_cod_vis;
			hdf.range_max = range_max_cod_vis;
			strcpy(hdf.units, units_cod_vis);
			hdf.scaled_flg = scaled_flg_cod_vis;
			hdf.scaled_type = scaled_type_cod_vis_pixel;
			write_hdf(hdf, (void *) imgr2.cod_vis_low);
		} else if (strcmp(sds_name[i], "cloud_bottom") == 0) {
			hdf.range_min = range_min_cldbot;
			hdf.range_max = range_max_cldbot;
			strcpy(hdf.units, units_cldbot);
			hdf.scaled_flg = scaled_flg_cldbot;
			hdf.scaled_type = scaled_type_cldbot_pixel;
			write_hdf(hdf, (void *) imgr2.cldbot);
		} else if (strcmp(sds_name[i], "aerosol_bottom") == 0) {
			hdf.range_min = range_min_cldbot;
			hdf.range_max = range_max_cldbot;
			strcpy(hdf.units, units_cldbot);
			hdf.scaled_flg = scaled_flg_cldbot;
			hdf.scaled_type = scaled_type_cldbot_pixel;
			write_hdf(hdf, (void *) imgr2.aerobot);
		} else if (strcmp(sds_name[i], "cloud_bottom_high") == 0) {
			hdf.range_min = range_min_cldbot;
			hdf.range_max = range_max_cldbot;
			strcpy(hdf.units, units_cldbot);
			hdf.scaled_flg = scaled_flg_cldbot;
			hdf.scaled_type = scaled_type_cldbot_pixel;
			write_hdf(hdf, (void *) imgr2.cldbot_high);
		} else if (strcmp(sds_name[i], "cloud_bottom_low") == 0) {
			hdf.range_min = range_min_cldbot;
			hdf.range_max = range_max_cldbot;
			strcpy(hdf.units, units_cldbot);
			hdf.scaled_flg = scaled_flg_cldbot;
			hdf.scaled_type = scaled_type_cldbot_pixel;
			write_hdf(hdf, (void *) imgr2.cldbot_low);
		} else if (strcmp(sds_name[i], "cloud_thickness") == 0) {
			hdf.range_min = range_min_clddz;
			hdf.range_max = range_max_clddz;
			strcpy(hdf.units, units_clddz);
			hdf.scaled_flg = scaled_flg_clddz;
			hdf.scaled_type = scaled_type_clddz_pixel;
			write_hdf(hdf, (void *) imgr2.clddz);
		} else if (strcmp(sds_name[i], "aerosol_thickness") == 0) {
			hdf.range_min = range_min_clddz;
			hdf.range_max = range_max_clddz;
			strcpy(hdf.units, units_clddz);
			hdf.scaled_flg = scaled_flg_clddz;
			hdf.scaled_type = scaled_type_clddz_pixel;
			write_hdf(hdf, (void *) imgr2.aerodz);
		} else if (strcmp(sds_name[i], "cloud_thickness_high") == 0) {
			hdf.range_min = range_min_clddz;
			hdf.range_max = range_max_clddz;
			strcpy(hdf.units, units_clddz);
			hdf.scaled_flg = scaled_flg_clddz;
			hdf.scaled_type = scaled_type_clddz_pixel;
			write_hdf(hdf, (void *) imgr2.clddz_high);
		} else if (strcmp(sds_name[i], "cloud_thickness_low") == 0) {
			hdf.range_min = range_min_clddz;
			hdf.range_max = range_max_clddz;
			strcpy(hdf.units, units_clddz);
			hdf.scaled_flg = scaled_flg_clddz;
			hdf.scaled_type = scaled_type_clddz_pixel;
			write_hdf(hdf, (void *) imgr2.clddz_low);
		} else if (strcmp(sds_name[i], "cloud_fraction") == 0) {
			hdf.range_min = range_min_cldfrac;
			hdf.range_max = range_max_cldfrac;
			strcpy(hdf.units, units_cldfrac);
			hdf.scaled_flg = scaled_flg_cldfrac;
			hdf.scaled_type = scaled_type_cldfrac_pixel;
			write_hdf(hdf, (void *) imgr2.cldfrac);
		} else if (strcmp(sds_name[i], "aerosol_fraction") == 0) {
			hdf.range_min = range_min_cldfrac;
			hdf.range_max = range_max_cldfrac;
			strcpy(hdf.units, units_cldfrac);
			hdf.scaled_flg = scaled_flg_cldfrac;
			hdf.scaled_type = scaled_type_cldfrac_pixel;
			write_hdf(hdf, (void *) imgr2.aerofrac);
		} else if (strcmp(sds_name[i], "cloud_liquid_water_path") == 0) {
			hdf.range_min = range_min_cldlwp;
			hdf.range_max = range_max_cldlwp;
			strcpy(hdf.units, units_cldlwp);
			hdf.scaled_flg = scaled_flg_cldlwp;
			hdf.scaled_type = scaled_type_cldlwp_pixel;
			write_hdf(hdf, (void *) imgr2.cldlwp);
		} else if (strcmp(sds_name[i], "cloud_ice_water_path") == 0) {
			hdf.range_min = range_min_cldiwp;
			hdf.range_max = range_max_cldiwp;
			strcpy(hdf.units, units_cldiwp);
			hdf.scaled_flg = scaled_flg_cldiwp;
			hdf.scaled_type = scaled_type_cldiwp_pixel;
			write_hdf(hdf, (void *) imgr2.cldiwp);
		} else if (strcmp(sds_name[i], "cloud_liquid_water_content") == 0) {
			hdf.range_min = range_min_cldlwc;
			hdf.range_max = range_max_cldlwc;
			strcpy(hdf.units, units_cldlwc);
			hdf.scaled_flg = scaled_flg_cldlwc;
			hdf.scaled_type = scaled_type_cldlwc_pixel;
			write_hdf(hdf, (void *) imgr2.cldlwc);
		} else if (strcmp(sds_name[i], "cloud_ice_water_content") == 0) {
			hdf.range_min = range_min_cldiwc;
			hdf.range_max = range_max_cldiwc;
			strcpy(hdf.units, units_cldiwc);
			hdf.scaled_flg = scaled_flg_cldiwc;
			hdf.scaled_type = scaled_type_cldiwc_pixel;
			write_hdf(hdf, (void *) imgr2.cldiwc);
		} else if (strcmp(sds_name[i], "cloud_beta") == 0) {
			hdf.range_min = range_min_cldbeta;
			hdf.range_max = range_max_cldbeta;
			strcpy(hdf.units, units_cldbeta);
			hdf.scaled_flg = scaled_flg_cldbeta;
			hdf.scaled_type = scaled_type_cldbeta_pixel;
			write_hdf(hdf, (void *) imgr2.cldbeta);
		} else if (strcmp(sds_name[i], "cloud_beta_high") == 0) {
			hdf.range_min = range_min_cldbeta;
			hdf.range_max = range_max_cldbeta;
			strcpy(hdf.units, units_cldbeta);
			hdf.scaled_flg = scaled_flg_cldbeta;
			hdf.scaled_type = scaled_type_cldbeta_pixel;
			write_hdf(hdf, (void *) imgr2.cldbeta_high);
		} else if (strcmp(sds_name[i], "cloud_beta_mid") == 0) {
			hdf.range_min = range_min_cldbeta;
			hdf.range_max = range_max_cldbeta;
			strcpy(hdf.units, units_cldbeta);
			hdf.scaled_flg = scaled_flg_cldbeta;
			hdf.scaled_type = scaled_type_cldbeta_pixel;
			write_hdf(hdf, (void *) imgr2.cldbeta_mid);
		} else if (strcmp(sds_name[i], "cloud_beta_low") == 0) {
			hdf.range_min = range_min_cldbeta;
			hdf.range_max = range_max_cldbeta;
			strcpy(hdf.units, units_cldbeta);
			hdf.scaled_flg = scaled_flg_cldbeta;
			hdf.scaled_type = scaled_type_cldbeta_pixel;
			write_hdf(hdf, (void *) imgr2.cldbeta_low);
		} else if (strcmp(sds_name[i], "emiss11_high") == 0) {
			hdf.range_min = range_min_cldemiss;
			hdf.range_max = range_max_cldemiss;
			strcpy(hdf.units, units_cldemiss);
			hdf.scaled_flg = scaled_flg_cldemiss;
			hdf.scaled_type = scaled_type_cldemiss_pixel;
			write_hdf(hdf, (void *) imgr2.emiss11_high);
		} else if (strcmp(sds_name[i], "emiss11_mid") == 0) {
			hdf.range_min = range_min_cldemiss;
			hdf.range_max = range_max_cldemiss;
			strcpy(hdf.units, units_cldemiss);
			hdf.scaled_flg = scaled_flg_cldemiss;
			hdf.scaled_type = scaled_type_cldemiss_pixel;
			write_hdf(hdf, (void *) imgr2.emiss11_mid);
		} else if (strcmp(sds_name[i], "emiss11_low") == 0) {
			hdf.range_min = range_min_cldemiss;
			hdf.range_max = range_max_cldemiss;
			strcpy(hdf.units, units_cldemiss);
			hdf.scaled_flg = scaled_flg_cldemiss;
			hdf.scaled_type = scaled_type_cldemiss_pixel;
			write_hdf(hdf, (void *) imgr2.emiss11_low);
		} else if (strcmp(sds_name[i], "column_aerosol_amount") == 0) {
			hdf.range_min = range_min_colaer;
			hdf.range_max = range_max_colaer;
			strcpy(hdf.units, units_colaer);
			hdf.scaled_flg = scaled_flg_colaer;
			hdf.scaled_type = scaled_type_colaer_pixel;
			write_hdf(hdf, (void *) imgr2.col_aer);
		} else if (strcmp(sds_name[i], "aerosol_optical_depth_vis") == 0) {
			hdf.range_min = range_min_aod_vis;
			hdf.range_max = range_max_aod_vis;
			strcpy(hdf.units, units_aod_vis);
			hdf.scaled_flg = scaled_flg_aod_vis;
			hdf.scaled_type = scaled_type_aod_vis_pixel;
			write_hdf(hdf, (void *) imgr2.aod_vis);
		} else if (strcmp(sds_name[i], "aerosol_optical_depth_ir") == 0) {
			hdf.range_min = range_min_aod_ir;
			hdf.range_max = range_max_aod_ir;
			strcpy(hdf.units, units_aod_ir);
			hdf.scaled_flg = scaled_flg_aod_ir;
			hdf.scaled_type = scaled_type_aod_ir_pixel;
			write_hdf(hdf, (void *) imgr2.aod_ir);
		} else if (strcmp(sds_name[i], "aerosol_particle_effective_radius")
				== 0) {
			hdf.range_min = range_min_aerreff;
			hdf.range_max = range_max_aerreff;
			strcpy(hdf.units, units_aerreff);
			hdf.scaled_flg = scaled_flg_aerreff;
			hdf.scaled_type = scaled_type_aerreff_pixel;
			write_hdf(hdf, (void *) imgr2.aerreff);
		} else if (strcmp(sds_name[i], "aerosol_particle_effective_diameter")
				== 0) {
			hdf.range_min = range_min_aerdeff;
			hdf.range_max = range_max_aerdeff;
			strcpy(hdf.units, units_aerdeff);
			hdf.scaled_flg = scaled_flg_aerdeff;
			hdf.scaled_type = scaled_type_aerdeff_pixel;
			write_hdf(hdf, (void *) imgr2.aerdeff);
		} else {
			fprintf(stderr,"%sUndefined algorithm output variable name - consult documentation-aborting\n",EXE_PROMPT);
			exit(EXIT_FAILURE);
		}
	}

}
