/*$Id: leocat.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/
/*-----------------------------------------------------------------
 LEOCAT

 Copyright (C) 2006  Michael J. Pavolonis
 National Oceanic and Atmospheric Administration

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.
 -----------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 Necessary header files
 ----------------------------------------------------------------------------*/
#define SFC_DIR ""
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <hdf.h>
#include <mfhdf.h>
#include "type_kind.h"
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "imagerL2_leocat.h"
#include "sounderL1_leocat.h"
#include "sounderL2_leocat.h"
#include "rtm_leocat.h"
#include "nwp_leocat.h"
#include "radutils_leocat.h"
#include "algorithm_interface.h"
#include "algorithms_leocat.h"
#include "input_leocat.h"
#include "ph_cld_type_threshold_coeffs.h"
#include "clavrx_volash_mask_threshold_coeffs.h"
#include "utils.h"

/*----------------------------------------------------------------------------
 Define the NEDR for each platform
 ----------------------------------------------------------------------------*/

const double NEDR_terra[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 8.00, 100.0, 
				1.0, 1.25 };
const double NEDR_aqua[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 8.00, 8.0, 1.0, 1.25 };
const double NEDR_npp[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0};
 /*const double NEDR_npoess[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,*/
 /*0.0, 0.0, 0.0, 0.0, 0.0, 0.0,*/
 /*0.0, 0.0, 0.0, 0.0, 0.0, 0.0,*/
 /*0.0, 0.0, 0.0, 0.0, 0.0, 0.0,*/
 /*0.0, 0.0, 0.0, 0.0, 0.0, 0.0,*/
 /*0.5, 0.5, 1.0, 1.0, 1.0, 1.0};*/

/*----------------------------------------------------------------------------
 Planck Function routine prototypes
 ----------------------------------------------------------------------------*/

// MODIS Planck functions
float modis_planck_terra(int, float);
float modis_planck_terra_shf(int, float);
float modis_planck_aqua(int, float);
float modis_planck_aqua_shf(int, float);
float planck_btemp_terra(int, float);
float planck_btemp_terra_shf(int, float);
float planck_btemp_aqua(int, float);
float planck_btemp_aqua_shf(int, float);
double rad_v2w_terra(int, double);
double rad_v2w_terra_shf(int, double);
double rad_v2w_aqua(int, double);
double rad_v2w_aqua_shf(int, double);
void load_fast_planck_modis(rad_utils *);
void destroy_fast_planck_modis(rad_utils *);

// VIIRS Planck functions
/*float viirs_planck_npp(int, float);*/
/*float viirs_planck_npoess(int, float);*/
/*float planck_btemp_npp(int, float);*/
/*float planck_btemp_npoess(int, float);*/
/*double rad_v2w_npp(int, double);*/
/*double rad_v2w_npoess(int, double);*/
/*void load_fast_planck_viirs(rad_utils *);*/
/*void destroy_fast_planck_viirs(rad_utils *);*/

/*----------------------------------------------------------------------------
 Instrument read routine prototypes
 ----------------------------------------------------------------------------*/

void read_modis(int8*, imagerL1*, int8, int8);
void read_modisqkm(int8*, imagerL1*, int8);
void read_viirs(int8*, imagerL1*, int8);
void airs_level1b_read(char *, sounderL1 *);

/*----------------------------------------------------------------------------
 Level 1b pointer destruction prototypes
 ----------------------------------------------------------------------------*/

void deallocate_modis(int8*, imagerL1 *);
void deallocate_viirs(int8*, imagerL1 *);

/*----------------------------------------------------------------------------
 NWP read routine prototypes
 ----------------------------------------------------------------------------*/

void main_gfs(char *, int, int, float, float *, nwp_params *);
void main_gdas(char *, int, int, float, float *, nwp_params *);
/*void main_ncep(char *, int, int, float, float *, nwp_params *);
 void main_ecmwf(char *, int, int, float, float *, nwp_params *);*/

/*----------------------------------------------------------------------------
 Level 2 output routine prototypes
 ----------------------------------------------------------------------------*/

void create_algorithm_pointers(int16, char *[], imagerL2 *);
void create_level2_default_output(int32, imagerL1, char *);
void create_level2_nwp_output(int32, int8, nwp_params, imagerL1);
void create_level2_output(int32, int16, int16, char *, char *, char *[], imagerL2);
void create_level2_rtm_output_modis(int32, imagerL1, nwp_params, rtm_profiles, int8 *);
/*void create_level2_rtm_output_viirs(int32, imagerL1, nwp_params, rtm_profiles, int8 *);*/
void create_rclr_ptrs(int8 *, rtm_toa *);

/*----------------------------------------------------------------------------
 Level 2 pointer destruction prototypes
 ----------------------------------------------------------------------------*/

void destroy_algorithm_pointers(int16, char *[], imagerL2 *);
void destroy_nwp_pointers(profile_params *);
void destroy_rtm_ptrs_modis(rtm_profiles *, nwp_params, int8 *);
/*void destroy_rtm_ptrs_viirs(rtm_profiles *, nwp_params, int8 *);*/
void destroy_rclr_ptrs(int8 *, rtm_toa *);

/*----------------------------------------------------------------------------
 RTM routine prototypes
 ----------------------------------------------------------------------------*/

void run_plod_modis(rtm_profiles *, profile_params *, rad_utils, int, float *,
		float, int8 *, int, char *, int *, int *, rtm_toa *, int);
void run_crtm_modis(rtm_profiles *, profile_params *, rad_utils, int, float *,
		float, int8 *, int, char *, int *, int *, rtm_toa *, int);
/*void run_plod_viirs(rtm_profiles *, profile_params *, rad_utils,
 int, float *, float, int8 *, int);
 void run_crtm_viirs(rtm_profiles *, profile_params *, rad_utils,
 int, float *, float, int8 *, int);*/

void get_clear_toa_rad(rtm_profiles, profile_params, rad_utils, int, rtm_toa *,
		int8 *, int);

/*----------------------------------------------------------------------------
 Utility routines just used in this routine
 ----------------------------------------------------------------------------*/

void version();
void help();
void print_alist();
void check_command_line_input(int, int, int, char**, char*);
void print_input_options(input_options, int8, int8, int8);

unsigned char * IGBP_grid2pixel(float *, float *, long, unsigned char *);
float * elev_grid2pixel(float *, float *, long, int16 *);
float * sst_grid2pixel(float *, float *, long, int, float *);
void nwp_grid2pixel(nwp_params *, imagerL1 *);
void read_seebor_emiss_imager(imagerL1 *, rtm_toa *);
void read_constant_emiss_imager(imagerL1 *, rtm_toa *);
void read_IGBP_map(unsigned char *, char *, char *);

/*----------------------------------------------------------------------------
 Start of main processing routine
 ----------------------------------------------------------------------------*/
char gdas_file[200];
char gdas_file1[200];
char gdas_file2[200];
char nise_file[200];
char sst_file[200];

extern void algorithmInformation_15(void *initializeMe);
extern void algorithmInformation_16(void *initializeMe);
extern void algorithmInformation_17(void *initializeMe);
extern void algorithmInformation_18(void *initializeMe);
extern void algorithmInformation_19(void *initializeMe);
extern void algorithmInformation_20(void *initializeMe);
extern void algorithmInformation_21(void *initializeMe);
extern void algorithmInformation_22(void *initializeMe);
extern void algorithmInformation_23(void *initializeMe);
extern void algorithmInformation_24(void *initializeMe);
extern void algorithmInformation_25(void *initializeMe);
extern void algorithmInformation_26(void *initializeMe);
extern void algorithmInformation_27(void *initializeMe);
extern void algorithmInformation_28(void *initializeMe);
extern void algorithmInformation_29(void *initializeMe);

extern void setImgr1Time(char *YYYYDDD, char *HHMM, imagerL1 *imgr1);



int main(int argc, char* argv[]) {

	printf("\n\nExecuting : %s\n\n",argv[0]); //GPC

	/*----------------------------------------------------------------------------
	 Declare needed variables
	 ----------------------------------------------------------------------------*/
	int i;
	int ii;
	FILE *fptr;

	int8 default_cmask_needed = NO, default_ctype_needed = NO,
			default_amask_needed = NO, nwp_needed = NO, rtm_needed = NO;
	int index, ichan, nchan, ifile, nchan_rtm, nwp_rclr_flag;
	int icell;
	int ivza;
	int iclr;
	char comment[MAX_STR_LEN];
	unsigned char *lmask_map;
	int16 *elev_map;
	float *sst_map;
	char *rout = { "leocat.c" }, level2_filename[MAX_STR_LEN],
			  level2_filenamebuff[MAX_STR_LEN],level3_filename[MAX_STR_LEN];
	long nlmask_map, nelev_map, nsst_map;
	float satzen_mid, *surface_emissivity;
	int32 l2_sd_id=0, l3_sd_id=0, status;

	/*----------------------------------------------------------------------------
	 Declare needed type define structures
	 ----------------------------------------------------------------------------*/

	input_options cl, def, opt;
	imagerL1 imgr1;
	imagerL2 imgr2;
	sounderL1 sndr1;
	sounderL2 sndr2;
	nwp_params nwp;
	rtm_profiles **rtm=NULL;
	rtm_toa rclr;
	rad_utils rutil;

	/*----------------------------------------------------------------------------
	 Declare needed pointers to functions
	 ----------------------------------------------------------------------------*/

	void (*rtm_function_ptr)(rtm_profiles *, profile_params *, rad_utils, int,
			float *, float, int8 *, int, char *, int *, int *, rtm_toa *, 
                        int)=NULL;

	void (*destroy_l1b_ptrs)(int8*, imagerL1*)=NULL;

	void (*destroy_fast_planck)(rad_utils *)=NULL;
	void (*load_fast_planck)(rad_utils *)=NULL;

	void (*destroy_rtm_ptrs)(rtm_profiles *, nwp_params, int8 *)=NULL;

	/*----------------------------------------------------------------------------
	 Declare structures used to help keep track of processing time.
	 ----------------------------------------------------------------------------*/

	time_t start_time_total, end_time_total;
	time_t start_time_algo, end_time_algo;
	double tdiff_total, tdiff_algo;
	int dsec_total, dsec_algo;
	int dmin_total, dmin_algo;
	int dhr_total, dhr_algo;
	int dday_total, dday_algo;
	int l2_name_set=0;

	char defaultFilenameWithPath[80];

	/*----------------------------------------------------------------------------
	 Begin timing the total processing time.
	 ----------------------------------------------------------------------------*/

	start_time_total = time(NULL);

	/*----------------------------------------------------------------------------
	 Print software version information.
	 ----------------------------------------------------------------------------*/

	version();

        printf("\nGetting algorithm information:\n");

//  algorithmInformation_15(&algo_info[14]);
//	algorithmInformation_16(&algo_info[15]);
	algorithmInformation_17(&algo_info[16]);
//	algorithmInformation_18(&algo_info[17]);
//	algorithmInformation_19(&algo_info[18]);
	//algorithmInformation_20(&algo_info[19]);	// v1.5.0.18 VCM
//	algorithmInformation_21(&algo_info[20]);	// v1.5.0.48 VCM
//      algorithmInformation_22(&algo_info[21]);
//      algorithmInformation_23(&algo_info[22]);
//      algorithmInformation_24(&algo_info[23]);
//	algorithmInformation_25(&algo_info[24]);
//	algorithmInformation_26(&algo_info[25]);
//	algorithmInformation_27(&algo_info[26]);
//	algorithmInformation_28(&algo_info[27]);
	algorithmInformation_29(&algo_info[28]);

	/*----------------------------------------------------------------------------
	 Initialize input options.
	 ----------------------------------------------------------------------------*/

	imgr2.output_cm_stats_flag = 0;
	imgr1.proc_250 = NO;
	imgr2.output_250 = NO;
	cl.l1b_data_res = -1.0;
	cl.level3_res = -1.0;
	strcpy(cl.nwp_source_name, UNKNOWN_STR);
	strcpy(cl.rtm_source_name, UNKNOWN_STR);
	strcpy(cl.data_source_name, UNKNOWN_STR);
	strcpy(cl.cmask_source_name, UNKNOWN_STR);
	strcpy(cl.ctype_source_name, UNKNOWN_STR);
	strcpy(cl.amask_source_name, UNKNOWN_STR);
	strcpy(cl.l1b_dir_name, UNKNOWN_STR);
	strcpy(cl.cmask_dir_name, UNKNOWN_STR);
	strcpy(cl.cmask_file_name, UNKNOWN_STR);
	strcpy(cl.cprod_dir_name, UNKNOWN_STR);
	strcpy(cl.nwp_dir_name, UNKNOWN_STR);
	strcpy(cl.l2_dir_name, UNKNOWN_STR);
	strcpy(cl.l3_dir_name, UNKNOWN_STR);
	strcpy(cl.geo_file, UNKNOWN_STR);
	strcpy(cl.sfc_dir_name, UNKNOWN_STR);
	strcpy(cl.algData_dir_name, UNKNOWN_STR);
	strcpy(cl.staticData_dir_name, UNKNOWN_STR);
	strcpy(cl.l1b_file_name, UNKNOWN_STR);
	strcpy(cl.qkm_file_name, UNKNOWN_STR);
	strcpy(cl.csrb_file_name, UNKNOWN_STR);
	strcpy(cl.csrb_dir_name, UNKNOWN_STR);
	strcpy(cl.flist_name, default_flist_name);
	cl.make_level2 = NOTKNOWN;
	cl.make_level3 = NOTKNOWN;
	cl.rtm_out = NOTKNOWN;
	cl.nwp_out = NOTKNOWN;
	cl.geo_interp = NO;
	cl.single_file = NO;
	cl.verbose = NO;
        cl.aqua_plod_spectral_shift = NO;
        cl.terra_plod_spectral_shift = NO;
	cl.algo2process = 0;
	for (i=0; i<NALGO; i++)
		cl.algo_flg[i] = NO;

	/*----------------------------------------------------------------------------
	 Parse any command line arguments.
	 ----------------------------------------------------------------------------*/
	strcat(defaultFilenameWithPath, default_file_name);
	memset(gdas_file, 0, 200);
	memset(gdas_file1, 0, 200);
	memset(gdas_file2, 0, 200);
	for (i=1; i<argc; i++) {
		if (strcmp(argv[i], "-help") == 0) {
			help();
			exit(EXIT_SUCCESS);
		} else if (strcmp(argv[i], "-alist") == 0) {
			print_alist();
			exit(EXIT_SUCCESS);
		} else if (strcmp(argv[i], "-wd") == 0) {
			check_command_line_input(++i, argc, YES, argv, "wd");
			defaultFilenameWithPath[0]=0;
			strcat(defaultFilenameWithPath, argv[i]);
			strcat(defaultFilenameWithPath, "/");

			strcat(defaultFilenameWithPath, default_file_name);

		} else if (strcmp(argv[i], "-satType") == 0) {
			check_command_line_input(++i, argc, YES, argv, "satType");

			imgr1.directbc_flg = NO;

			if (strcmp(argv[i], "aqua") == 0) {
				strcpy(imgr1.satname, "Aqua");
				imgr1.satid = AQUA;
			} else if (strcmp(argv[i], "terra") == 0) {
				strcpy(imgr1.satname, "Terra");
				imgr1.satid = TERRA;
			} else if (strcmp(argv[i], "npp") == 0) {
				strcpy(imgr1.satname, "NPP");
				imgr1.satid = NPP;
			} else if (strcmp(argv[i], "npoess") == 0) {
				strcpy(imgr1.satname, "NPOESS");
				imgr1.satid = NPOESS;
			} else {
				printf("Error need satellite id");
			}

		} else if (strcmp(argv[i], "-satYYYDDD") == 0) {
			check_command_line_input(++i, argc, YES, argv, "satYYYDDD");

			setImgr1Time(argv[i], NULL, &imgr1);

		} else if (strcmp(argv[i], "-satHHMM") == 0) {
			check_command_line_input(++i, argc, YES, argv, "satHHMM");

			setImgr1Time(NULL, argv[i], &imgr1);

		}

		else if (strcmp(argv[i], "-geo_file") == 0) {
			check_command_line_input(++i, argc, YES, argv, "geo_file");
			strcpy(imgr1.geo_filename1, argv[i]);
			strcpy(cl.geo_file, argv[i]);
			printf("GeoFile: %s %s\n", imgr1.geo_filename1, cl.geo_file);

		} else if (strcmp(argv[i], "-gdas_file") == 0) {
			check_command_line_input(++i, argc, YES, argv, "gdas_file");
			strcpy(gdas_file, argv[i]);
			printf("Gdasfile: %s\n", gdas_file);

                } else if (strcmp(argv[i], "-gdas_file1") == 0) {
                        check_command_line_input(++i, argc, YES, argv, "gdas_file1");
                        strcpy(gdas_file1, argv[i]);
                        printf("Gdasfile1: %s\n", gdas_file1);

                } else if (strcmp(argv[i], "-gdas_file2") == 0) {
                        check_command_line_input(++i, argc, YES, argv, "gdas_file2");
                        strcpy(gdas_file2, argv[i]);
                        printf("Gdasfile2: %s\n", gdas_file2);

		} else if (strcmp(argv[i], "-nise_file") == 0) {
			check_command_line_input(++i, argc, YES, argv, "nise_file");
			strcpy(nise_file, argv[i]);
			printf("nise_file: %s\n", nise_file);

                } else if (strcmp(argv[i], "-sst_file") == 0) {
                        check_command_line_input(++i, argc, YES, argv, "sst_file");
                        strcpy(sst_file, argv[i]);
                        printf("sst_file: %s\n", sst_file);

		} else if (strcmp(argv[i], "-l2_name") == 0) {
			check_command_line_input(++i, argc, YES, argv, "l2_name");
			strcpy(level2_filename, argv[i]);
			l2_name_set=1;
		}

		else if (strcmp(argv[i], "-data") == 0) {
			check_command_line_input(++i, argc, YES, argv, "data");
			strcpy(cl.data_source_name, argv[i]);
		} else if (strcmp(argv[i], "-data_res") == 0) {
			check_command_line_input(++i, argc, NO, argv, "data_res");
			cl.l1b_data_res = atof(argv[i]);
		} else if (strcmp(argv[i], "-nwp") == 0) {
			check_command_line_input(++i, argc, YES, argv, "nwp");
			strcpy(cl.nwp_source_name, argv[i]);

			printf("nwp_source_name:%s\n", argv[i]);

		} else if (strcmp(argv[i], "-rtm") == 0) {
			check_command_line_input(++i, argc, YES, argv, "rtm");
			strcpy(cl.rtm_source_name, argv[i]);
		} else if (strcmp(argv[i], "-nol2") == 0) {
			cl.make_level2 = NO;
		} else if (strcmp(argv[i], "-nol3") == 0) {
			cl.make_level3 = NO;
		} else if (strcmp(argv[i], "-rtm_out") == 0) {
			cl.rtm_out = YES;
		} else if (strcmp(argv[i], "-nwp_out") == 0) {
			cl.nwp_out = YES;
		} else if (strcmp(argv[i], "-geo_interp") == 0) {
			cl.geo_interp = YES;
		} else if (strcmp(argv[i], "-l3res") == 0) {
			check_command_line_input(++i, argc, NO, argv, "l3res");
			cl.level3_res = atof(argv[i]);
		} else if (strcmp(argv[i], "-cmask") == 0) {
			check_command_line_input(++i, argc, YES, argv, "cmask");
			strcpy(cl.cmask_source_name, argv[i]);
		} else if (strcmp(argv[i], "-ctype") == 0) {
			check_command_line_input(++i, argc, YES, argv, "ctype");
			strcpy(cl.ctype_source_name, argv[i]);
		} else if (strcmp(argv[i], "-amask") == 0) {
			check_command_line_input(++i, argc, YES, argv, "amask");
			strcpy(cl.amask_source_name, argv[i]);
		} else if (strcmp(argv[i], "-a") == 0) {
			if (++i == argc) {
				fprintf(
						stderr,
						"%sInvalid command line input option following [-a], displaying help...\n",
						EXE_PROMPT);
				help();
				exit(EXIT_FAILURE);
			}
			while (argv[i][0] != '-') {
				index = atoi(argv[i])-1;
				printf("\nAlgorithm number: %d\n",index+1);
				if (index < 0 || index >= NALGO) {
					fprintf(
							stderr,
							"%sInvalid algorithm index found, displaying algorithm information...\n",
							EXE_PROMPT);
					print_alist();
					exit(EXIT_FAILURE);
				}
				algo_order[cl.algo2process] = index+1;
				cl.algo_flg[index] = YES;
				cl.algo2process++;
				if (++i == argc)
					break;
			}
			i--;
		} else if (strcmp(argv[i], "-l1b_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "l1b_dir");
			strcpy(cl.l1b_dir_name, argv[i]);
		} else if (strcmp(argv[i], "-l1b_qkm_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "l1b_qkm_dir");
			strcpy(cl.qkm_dir_name, argv[i]);
		} else if (strcmp(argv[i], "-csrb_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "csrb_dir");
			strcpy(cl.csrb_dir_name, argv[i]);
                        strcpy(imgr1.csrb_dir_name, cl.csrb_dir_name);
			printf("CSRB directory: %s \n", cl.csrb_dir_name);
		} else if (strcmp(argv[i], "-cmask_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "cmask_dir");
			strcpy(cl.cmask_dir_name, argv[i]);
		} else if (strcmp(argv[i], "-cmask_file") == 0) {
			check_command_line_input(++i, argc, YES, argv, "cmask_file");
			strcpy(cl.cmask_file_name, argv[i]);
		} else if (strcmp(argv[i], "-cprod_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "cprod_dir");
			strcpy(cl.cprod_dir_name, argv[i]);
		} else if (strcmp(argv[i], "-sfc_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "sfc_dir");
			strcpy(cl.sfc_dir_name, argv[i]);
			strcpy(imgr1.sfc_dir_name, argv[i]);
			printf("sfc_data directory: %s \n", imgr1.sfc_dir_name);
		} else if (strcmp(argv[i], "-algData_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "algData_dir");
			strcpy(cl.algData_dir_name, argv[i]);
			strcpy(imgr1.algData_dir_name, argv[i]);
			printf("algData_dir_name: %s \n", imgr1.algData_dir_name);
		} else if (strcmp(argv[i], "-staticData_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "staticData_dir");
			strcpy(cl.staticData_dir_name, argv[i]);
			strcpy(imgr1.staticData_dir_name, argv[i]);
			printf("staticData_dir_name: %s \n", imgr1.staticData_dir_name);
		} else if (strcmp(argv[i], "-nwp_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "nwp_dir");
			strcpy(cl.nwp_dir_name, argv[i]);
		} else if (strcmp(argv[i], "-l2_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "l2_dir");
			strcpy(cl.l2_dir_name, argv[i]);
		} else if (strcmp(argv[i], "-l3_dir") == 0) {
			check_command_line_input(++i, argc, YES, argv, "l3_dir");
			strcpy(cl.l3_dir_name, argv[i]);
		} else if (strcmp(argv[i], "-qkm_file") == 0) {
			check_command_line_input(++i, argc, YES, argv, "qkm_file");
			strcpy(cl.qkm_file_name, argv[i]);
			imgr1.proc_250 = YES;
			imgr2.output_250 = YES;
		} else if (strcmp(argv[i], "-csrb_file") == 0) {
			check_command_line_input(++i, argc, YES, argv, "csrb_file");
			strcpy(cl.csrb_file_name, argv[i]);
                        strcpy(imgr1.csrb_filename, cl.csrb_file_name);
			printf("CSRB file: %s \n", cl.csrb_file_name);
		} else if (strcmp(argv[i], "-f") == 0) {
			check_command_line_input(++i, argc, YES, argv, "f");
			strcpy(cl.l1b_file_name, argv[i]);
			cl.single_file = YES;
		} else if (strcmp(argv[i], "-flist_name") == 0) {
			check_command_line_input(++i, argc, YES, argv, "flist_name");
			strcpy(cl.flist_name, argv[i]);
		} else if (strcmp(argv[i], "-verbose") == 0) {
			cl.verbose = YES;
		} else if (strcmp(argv[i], "-aqua_plod_spectral_shift") == 0) {
			cl.aqua_plod_spectral_shift= YES;
		} else if (strcmp(argv[i], "-terra_plod_spectral_shift") == 0) {
			cl.terra_plod_spectral_shift= YES;
		} else if (strcmp(argv[i], "-version") == 0 || strcmp(argv[i], "-v")
				== 0) {
			exit(EXIT_SUCCESS);
		} else {
			fprintf(stderr,"%sUnknown input option [%s], displaying help...\n",
					EXE_PROMPT, argv[i]);
			help();
			exit(EXIT_FAILURE);
		}
	}

	/*----------------------------------------------------------------------------
	 Open the default file.
	 ----------------------------------------------------------------------------*/

	if ((fptr = fopen(defaultFilenameWithPath, "r")) == NULL) {
		fprintf(
				stderr,
				"%sCannot open the default file %s\nUSER:check file location and format - aborting\n",
				EXE_PROMPT, defaultFilenameWithPath);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Read in the default values.
	 ----------------------------------------------------------------------------*/

	fscanf(fptr, "%s", &def.data_source_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%f", &def.l1b_data_res);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.nwp_source_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.rtm_source_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.rtm_out);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.nwp_out);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%d", &def.make_level2);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%d", &def.make_level3);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%f", &def.level3_res);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.cmask_source_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.ctype_source_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.amask_source_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	def.algo2process = 0;
	for (i=0; i<NALGO; i++) {
		fscanf(fptr, "%d", &def.algo_flg[i]);
		if (def.algo_flg[i] == 1)
			def.algo2process++;
		fgets(comment, MAX_STR_LEN, fptr);
	}
	fscanf(fptr, "%s", &def.l1b_dir_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.geo_file[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.cmask_dir_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.cmask_file_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.cprod_dir_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.nwp_dir_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.l2_dir_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.l3_dir_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);
	fscanf(fptr, "%s", &def.sfc_dir_name[0]);
	fgets(comment, MAX_STR_LEN, fptr);

	/*----------------------------------------------------------------------------
	 Close the default file.
	 ----------------------------------------------------------------------------*/

	fclose(fptr);

	/*----------------------------------------------------------------------------
	 Determine actual options to be used based on combined command line and
	 default values.
	 ----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------
	 Take care of the strings first.
	 ----------------------------------------------------------------------------*/

	if (strcmp(cl.data_source_name, UNKNOWN_STR) == 0){
		strcpy(opt.data_source_name, def.data_source_name);
	}else{
		strcpy(opt.data_source_name, cl.data_source_name);
	}
	if (strcmp(cl.nwp_source_name, UNKNOWN_STR) == 0)
		strcpy(opt.nwp_source_name, def.nwp_source_name);
	else
		strcpy(opt.nwp_source_name, cl.nwp_source_name);

	if (strcmp(cl.rtm_source_name, UNKNOWN_STR) == 0)
		strcpy(opt.rtm_source_name, def.rtm_source_name);
	else
		strcpy(opt.rtm_source_name, cl.rtm_source_name);

	if (strcmp(cl.cmask_source_name, UNKNOWN_STR) == 0)
		strcpy(opt.cmask_source_name, def.cmask_source_name);
	else
		strcpy(opt.cmask_source_name, cl.cmask_source_name);

	if (strcmp(cl.ctype_source_name, UNKNOWN_STR) == 0)
		strcpy(opt.ctype_source_name, def.ctype_source_name);
	else
		strcpy(opt.ctype_source_name, cl.ctype_source_name);

	if (strcmp(cl.amask_source_name, UNKNOWN_STR) == 0)
		strcpy(opt.amask_source_name, def.amask_source_name);
	else
		strcpy(opt.amask_source_name, cl.amask_source_name);

	if (strcmp(cl.l1b_dir_name, UNKNOWN_STR) == 0)
		strcpy(opt.l1b_dir_name, def.l1b_dir_name);
	else
		strcpy(opt.l1b_dir_name, cl.l1b_dir_name);

	if (strcmp(cl.geo_file, UNKNOWN_STR) == 0)
		strcpy(opt.geo_file, def.geo_file);
	else
		strcpy(opt.geo_file, cl.geo_file);

	if (strcmp(cl.sfc_dir_name, UNKNOWN_STR) == 0)
		strcpy(opt.sfc_dir_name, def.sfc_dir_name);
	else
		strcpy(opt.sfc_dir_name, cl.sfc_dir_name);

	if (strcmp(cl.cmask_dir_name, UNKNOWN_STR) == 0)
//              strcpy(opt.cmask_dir_name, def.cmask_dir_name);
		strcpy(opt.cmask_dir_name, cl.cmask_dir_name);
	else
		strcpy(opt.cmask_dir_name, cl.cmask_dir_name);

	if (strcmp(cl.cmask_file_name, UNKNOWN_STR) == 0)
                strcpy(opt.cmask_file_name, def.cmask_file_name);
	else
		strcpy(opt.cmask_file_name, cl.cmask_file_name);

	if (strcmp(cl.cprod_dir_name, UNKNOWN_STR) == 0)
		strcpy(opt.cprod_dir_name, def.cprod_dir_name);
	else
		strcpy(opt.cprod_dir_name, cl.cprod_dir_name);

	if (strcmp(cl.nwp_dir_name, UNKNOWN_STR) == 0)
		strcpy(opt.nwp_dir_name, def.nwp_dir_name);
	else
		strcpy(opt.nwp_dir_name, cl.nwp_dir_name);

	if (strcmp(cl.l2_dir_name, UNKNOWN_STR) == 0)
		strcpy(opt.l2_dir_name, def.l2_dir_name);
	else
		strcpy(opt.l2_dir_name, cl.l2_dir_name);

	if (strcmp(cl.l3_dir_name, UNKNOWN_STR) == 0)
		strcpy(opt.l3_dir_name, def.l3_dir_name);
	else
		strcpy(opt.l3_dir_name, cl.l3_dir_name);

	strcpy(imgr1.sfc_dir_name, opt.sfc_dir_name);
        strcpy(imgr1.cmask_dir_name, opt.cmask_dir_name);
        strcpy(imgr1.cmask_file_name, opt.cmask_file_name);
	strcpy(imgr1.cprod_dir_name, opt.cprod_dir_name);

	/*----------------------------------------------------------------------------
	 Now the numerical values.
	 ----------------------------------------------------------------------------*/

	opt.l1b_data_res = (cl.l1b_data_res == (float)NOTKNOWN) ? def.l1b_data_res
			: cl.l1b_data_res;
	opt.level3_res = (cl.level3_res == (float)NOTKNOWN) ? def.level3_res
			: cl.level3_res;
	opt.rtm_out = (cl.rtm_out == NOTKNOWN) ? def.rtm_out : cl.rtm_out;
	opt.nwp_out = (cl.nwp_out == NOTKNOWN) ? def.nwp_out : cl.nwp_out;
	opt.make_level2 = (cl.make_level2 == NOTKNOWN) ? def.make_level2
			: cl.make_level2;
	opt.make_level3 = (cl.make_level3 == NOTKNOWN) ? def.make_level3
			: cl.make_level3;
	opt.verbose = cl.verbose;
        opt.aqua_plod_spectral_shift = cl.aqua_plod_spectral_shift;
        opt.terra_plod_spectral_shift = cl.terra_plod_spectral_shift;
	opt.geo_interp = cl.geo_interp;

	if (cl.algo2process == 0) {
		for (i=0; i<NALGO; i++)
			opt.algo_flg[i] = def.algo_flg[i];
		opt.algo2process = def.algo2process;
	} else {
		for (i=0; i<NALGO; i++)
			opt.algo_flg[i] = cl.algo_flg[i];
		opt.algo2process = cl.algo2process;
	}

	if (opt.algo2process == 0) {
		fprintf(
				stderr,
				"%sNo algorithms were specified, the default algorithms are not excuted if\n\tadditional algorithms are not specified.\n",
				EXE_PROMPT);
		fprintf(
				stderr,
				"%sChoose at least 1 algorithm from list (default cloud mask/cloud type can be chosen from list).\n",
				EXE_PROMPT);
		print_alist();
		exit(EXIT_FAILURE);
	}

	strcpy(opt.flist_name, cl.flist_name);
	opt.single_file = cl.single_file;

	/*----------------------------------------------------------------------------
	 Convert certain string options to useable integers.
	 ----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------
	 Data source
	 ----------------------------------------------------------------------------*/

	if (strcmp(opt.data_source_name, "modis") == 0) {
		opt.data_source = MODIS_DAT;
		opt.l1b_data_res = 1.0;
		nchan = NMODIS_CHAN;
		nchan_rtm = NMODIS_CHAN_CRTM;
		rutil.planck_rad_fast_ptr = planck_rad_fast_modis;
		rutil.planck_bt_fast_ptr = planck_bt_fast_modis;
		rutil.planck_rad_fast_index_ptr = planck_rad_fast_index_modis;
		rutil.planck_bt_fast_index_ptr = planck_bt_fast_index_modis;
		destroy_l1b_ptrs = deallocate_modis;
		load_fast_planck = load_fast_planck_modis;
		destroy_fast_planck = destroy_fast_planck_modis;
		destroy_rtm_ptrs = destroy_rtm_ptrs_modis;
	} else if (strcmp(opt.data_source_name, "modis5") == 0) {
		opt.data_source = MODIS5_DAT;
		opt.l1b_data_res = 5.0;
		nchan = NMODIS_CHAN;
		nchan_rtm = NMODIS_CHAN_CRTM;
		rutil.planck_rad_fast_ptr = planck_rad_fast_modis;
		rutil.planck_bt_fast_ptr = planck_bt_fast_modis;
		rutil.planck_rad_fast_index_ptr = planck_rad_fast_index_modis;
		rutil.planck_bt_fast_index_ptr = planck_bt_fast_index_modis;
		destroy_l1b_ptrs = deallocate_modis;
		load_fast_planck = load_fast_planck_modis;
		destroy_fast_planck = destroy_fast_planck_modis;
		destroy_rtm_ptrs = destroy_rtm_ptrs_modis;
	} else if (strcmp(opt.data_source_name, "modis_avgx") == 0) {
		opt.data_source = MODIS_AVGX_DAT;
		sprintf(opt.data_source_name, "modis_avg_%3.1f", opt.l1b_data_res);
		nchan = NMODIS_CHAN;
		nchan_rtm = NMODIS_CHAN_CRTM;
		rutil.planck_rad_fast_ptr = planck_rad_fast_modis;
		rutil.planck_bt_fast_ptr = planck_bt_fast_modis;
		rutil.planck_rad_fast_index_ptr = planck_rad_fast_index_modis;
		rutil.planck_bt_fast_index_ptr = planck_bt_fast_index_modis;
		destroy_l1b_ptrs = deallocate_modis;
		load_fast_planck = load_fast_planck_modis;
		destroy_fast_planck = destroy_fast_planck_modis;
		destroy_rtm_ptrs = destroy_rtm_ptrs_modis;
	} else if (strcmp(opt.data_source_name, "viirs") == 0) {
		opt.data_source = VIIRS_DAT;
		opt.make_level2 = NO;
		opt.l1b_data_res = 0.75;
		nchan = NVIIRS_CHAN;
		nchan_rtm = NVIIRS_CHAN_CRTM;
		/*rutil.planck_rad_fast_ptr = planck_rad_fast_viirs;*/
		/*rutil.planck_bt_fast_ptr = planck_bt_fast_viirs;*/
		/*rutil.planck_rad_fast_index_ptr = planck_rad_fast_index_viirs;*/
		/*rutil.planck_bt_fast_index_ptr = planck_bt_fast_index_viirs;*/
		destroy_l1b_ptrs = deallocate_viirs;
		/*load_fast_planck = load_fast_planck_viirs;*/
		/*destroy_fast_planck = destroy_fast_planck_viirs;*/
		/*destroy_rtm_ptrs = destroy_rtm_ptrs_viirs;*/
	} else {
		fprintf(
				stderr,
				"%sUnknown data source, [%s], check command line option [-data] or %s, displaying help...\n",
		EXE_PROMPT,opt.data_source_name,default_file_name);
		exit(EXIT_FAILURE);
	}

	strcpy(imgr1.instrumentname, opt.data_source_name);
	imgr1.data_source = opt.data_source;
	imgr1.l1b_data_res = opt.l1b_data_res;
	strcpy(imgr1.l1b_data_res_units, "km");

	/*----------------------------------------------------------------------------
	 Default nwp
	 ----------------------------------------------------------------------------*/

	if (strcmp(opt.nwp_source_name, "gfs") == 0)
		opt.nwp_source = GFS_NWP;
	else if (strcmp(opt.nwp_source_name, "gdas") == 0)
		opt.nwp_source = GDAS_NWP;
	else if (strcmp(opt.nwp_source_name, "ncep") == 0)
		opt.nwp_source = NCEP_NWP;
	else if (strcmp(opt.nwp_source_name, "ecmwf") == 0)
		opt.nwp_source = ECMWF_NWP;
	else {
		fprintf(
				stderr,
				"%sUnknown nwp source, [%s], check command line option [-nwp] or %s, displaying help...\n",
				EXE_PROMPT, opt.nwp_source_name, default_file_name);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Default rtm
	 ----------------------------------------------------------------------------*/

	if (strcmp(opt.rtm_source_name, "crtm") == 0) {
		opt.rtm_source = CRTM_RTM;
		if (opt.data_source == MODIS_DAT || opt.data_source == MODIS5_DAT
				|| opt.data_source == MODIS_AVGX_DAT) {
			rtm_function_ptr = run_crtm_modis;
		} else if (opt.data_source == VIIRS_DAT) {
			/*rtm_function_ptr = run_crtm_viirs;*/
		}
	} else if (strcmp(opt.rtm_source_name, "plod") == 0) {
		opt.rtm_source = PLOD_RTM;
		if (opt.data_source == MODIS_DAT || opt.data_source == MODIS5_DAT
				|| opt.data_source == MODIS_AVGX_DAT) {
			rtm_function_ptr = run_plod_modis;
		} else if (opt.data_source == VIIRS_DAT) {
			/*rtm_function_ptr = run_plod_viirs;*/
		}
	} else {
		fprintf(
				stderr,
				"%sUnknown nwp source, [%s], check command line option [-rtm] or %s, displaying help...\n",
		EXE_PROMPT,opt.rtm_source_name,default_file_name);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Default cloud mask
	 ----------------------------------------------------------------------------*/

	if (strcmp(opt.cmask_source_name, "modis35") == 0)
		opt.cmask_source = MODIS35_CMASK;
	else {
		fprintf(
				stderr,
				"%sUnknown default cloud mask, [%s], check command line option [-cmask] or %s, displaying help...\n",
		EXE_PROMPT,opt.cmask_source_name,default_file_name);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Default cloud type
	 ----------------------------------------------------------------------------*/

	if (strcmp(opt.ctype_source_name, "ph05") == 0)
		opt.ctype_source = PH05_CTYPE;
	else if (strcmp(opt.ctype_source_name, "modis_ir") == 0)
		opt.ctype_source = MODIS_IR_CTYPE;
	else {
		fprintf(
				stderr,
				"%sUnknown default cloud type, [%s], check command line option [-ctype] or %s, displaying help...\n",
		EXE_PROMPT,opt.ctype_source_name,default_file_name);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Default aerosol mask
	 ----------------------------------------------------------------------------*/

	if (strcmp(opt.amask_source_name, "volash_4ch") == 0)
		opt.amask_source = VOLASH_4CH;
	else {
		fprintf(
				stderr,
				"%sUnknown default aerosol mask, [%s], check command line option [-amask] or %s, displaying help...\n",
		EXE_PROMPT,opt.amask_source_name,default_file_name);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Determine if the default cloud mask, cloud type, aerosol mask, nwp, and rtm
	 are needed by other algorithms.
	 ----------------------------------------------------------------------------*/

	for (i=0; i<NALGO; i++) {
		if ((opt.algo_flg[i] == YES && algo_info[i].cmask_needed == YES)
				|| (opt.algo_flg[i] == YES && i == opt.cmask_source-1))
			default_cmask_needed = YES;
		if ((opt.algo_flg[i] == YES && algo_info[i].ctype_needed == YES)
				|| (opt.algo_flg[i] == YES && i == opt.ctype_source-1))
			default_ctype_needed = YES;
		if ((opt.algo_flg[i] == YES && algo_info[i].amask_needed == YES)
				|| (opt.algo_flg[i] == YES && i == opt.amask_source-1))
			default_amask_needed = YES;
		if (opt.algo_flg[i] == YES && algo_info[i].rtm_needed == YES)
			rtm_needed = YES;
		if (opt.algo_flg[i] == YES && algo_info[i].nwp_needed == YES)
			nwp_needed = YES;
	}

	/*----------------------------------------------------------------------------
	 Make sure that the default cloud mask, cloud type, and aerosol mask
	 algorithms are not counted twice.
	 ----------------------------------------------------------------------------*/

	if (default_cmask_needed == YES) {
		if (opt.algo_flg[opt.cmask_source-1] == NO)
			opt.algo2process++;
		else
			opt.algo_flg[opt.cmask_source-1] = NO;
	}

	if (default_ctype_needed == YES) {
		if (opt.algo_flg[opt.ctype_source-1] == NO)
			opt.algo2process++;
		else
			opt.algo_flg[opt.ctype_source-1] = NO;
	}

	if (default_amask_needed == YES) {
		if (opt.algo_flg[opt.amask_source-1] == NO)
			opt.algo2process++;
		else
			opt.algo_flg[opt.amask_source-1] = NO;
	}

	/*----------------------------------------------------------------------------
	 Print the input options.
	 ----------------------------------------------------------------------------*/

	print_input_options(opt, default_cmask_needed, default_ctype_needed,
			default_amask_needed);

	/*----------------------------------------------------------------------------
	 Determine which instrument channels are needed.
	 ----------------------------------------------------------------------------*/

        printf("Determining which instrument channels are required...\n");
	if (strcmp(opt.data_source_name, "modis") == 0) {
		printf("Determining MODIS channels required...\n");
		if ((imgr1.chflg = (int8 *) malloc(sizeof(int8) * NMODIS_CHAN)) == NULL)error_allo(rout,"imgr1.chflg");
        if(imgr1.proc_250) {
			if ((imgr1.chflg_qkm = (int8 *) malloc(sizeof(int8) * N250M_CHAN)) == NULL)error_allo(rout,"imgr1.chflg_qkm");
        }
		if ((imgr1.chflg_rtm = (int8 *) malloc(sizeof(int8) * NMODIS_CHAN)) == NULL)error_allo(rout,"imgr1.chflg_rtm");
		if ((surface_emissivity = (float *) malloc(sizeof(float) * NMODIS_CHAN))
				== NULL)error_allo(rout,"surface_emissivity");
		if ((rutil.NEDR = (double *) malloc(sizeof(double) * NMODIS_CHAN)) == NULL)error_allo(rout,"rutil.NEDR");

		for (ichan=0; ichan<NMODIS_CHAN; ichan++) {
			imgr1.chflg[ichan] = 0;
			for (i=0; i<NALGO; i++) {
				imgr1.chflg[ichan] += algo_info[i].chflg[ichan]
						*opt.algo_flg[i];
			}
			imgr1.chflg[ichan]
				+= algo_info[opt.cmask_source - 1].chflg[ichan];
			if (default_ctype_needed == YES)
				imgr1.chflg[ichan]
						+= algo_info[opt.ctype_source-1].chflg[ichan];
						
		}

		if(imgr1.proc_250) {
			for (ichan=0; ichan<N250M_CHAN; ichan++) {
				imgr1.chflg_qkm[ichan] = 1;
			}
		}
		printf("\nFinal MODIS channel status...\n");
		for (ichan=0; ichan<NMODIS_CHAN; ichan++)
			printf("%3d",ichan+1);
		printf("\n");
		for (ichan=0; ichan<NMODIS_CHAN; ichan++)
			printf("%3d",imgr1.chflg[ichan]);

	} else if (strcmp(opt.data_source_name, "viirs") == 0) {
		// If we are calling multiple algorithms, run through each and determine channels required 
		// by entire suite.
		printf("Determining VIIRS channels required...\n");
		if ((imgr1.chflg = (int8 *) malloc(sizeof(int8) * NVIIRS_CHAN)) == NULL)error_allo(rout,"imgr1.chflg");
		if ((imgr1.chflg_rtm = (int8 *) malloc(sizeof(int8) * NVIIRS_CHAN)) == NULL)error_allo(rout,"imgr1.chflg_rtm");
		if ((surface_emissivity = (float *) malloc(sizeof(float) * NVIIRS_CHAN))
				== NULL)error_allo(rout,"surface_emissivity");
		if ((rutil.NEDR = (double *) malloc(sizeof(double) * NVIIRS_CHAN)) == NULL)error_allo(rout,"rutil.NEDR");

		// For every channel, determine if there is at least one algorithm that requires it.
		for (ichan=0; ichan<NVIIRS_CHAN; ichan++) {
			imgr1.chflg[ichan] = 0;
			for (i=0; i<NALGO; i++) {
				imgr1.chflg[ichan] += algo_info[i].chflg[ichan]*opt.algo_flg[i];
			}
		}
		printf("Final VIIRS channel status...\n");
		for (ichan=0; ichan<NVIIRS_CHAN; ichan++)
			printf("%d ",imgr1.chflg[ichan]);

	} else {
		fprintf(
				stderr,
				"%sUnknown data source, [%s], check command line option [-data] or %s, displaying help...\n",
		EXE_PROMPT,opt.data_source_name,default_file_name);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Determine the number of viewing angle bins used in the RTM.
	 ----------------------------------------------------------------------------*/

	if (rtm_needed == YES)
		nwp.rtm_nvzen = ((int) ceil(1.0 / RTM_VZA_BINSIZE)) + 1;

	/*----------------------------------------------------------------------------
	 Read in land mask.
	 ----------------------------------------------------------------------------*/

	if (opt.verbose == YES){
		fprintf(stdout,"\n%sReading land surface type map from opt.sfc_dir_name = %s\n",EXE_PROMPT,opt.sfc_dir_name);
	}

	nlmask_map = (long) nlon_IGBP * (long) nlat_IGBP;
	if ((lmask_map = (unsigned char *) malloc(sizeof(unsigned char)
	  		* nlmask_map)) == NULL)error_allo(rout,"lmask_map");
	read_IGBP_map(lmask_map, opt.sfc_dir_name,SFC_DIR);

	/*----------------------------------------------------------------------------
	 Read in surface elevation.
	 ----------------------------------------------------------------------------*/

	if (opt.verbose == YES){
		fprintf(stdout,"\n%sReading surface elevation map from opt.sfc_dir_name = %s\n",EXE_PROMPT,opt.sfc_dir_name);
	}
	nelev_map = (long) nlon_elev * (long) nlat_elev;
	if ((elev_map = (int16 *) malloc(sizeof(int16) * nelev_map)) == NULL)error_allo(rout,"elev_map");
	read_elev_map(elev_map, opt.sfc_dir_name,SFC_DIR);

	/*----------------------------------------------------------------------------
	 Read in climatological sst.
	 ----------------------------------------------------------------------------*/

	if (opt.verbose == YES){
		fprintf(stdout,"\n%sReading climatological sst map from opt.sfc_dir_name = %s\n",EXE_PROMPT,opt.sfc_dir_name);
	}
	nsst_map = (long) nlon_sst * (long) nlat_sst;
	if ((sst_map = (float *) malloc(sizeof(float) * nsst_map * 12)) == NULL)error_allo(rout,"sst_map");
	read_sst_map(sst_map, opt.sfc_dir_name,SFC_DIR);

	/*----------------------------------------------------------------------------
	 Open file list file, if needed.
	 ----------------------------------------------------------------------------*/

	if (opt.single_file == NO) {
		if ((fptr = fopen(opt.flist_name, "r")) == NULL) {
			fprintf(
					stderr,
					"%sCannot open the L1b list file %s\nUSER:check file location and format - aborting\n",
					EXE_PROMPT, opt.flist_name);
			exit(EXIT_FAILURE);
		}
	}

	/*----------------------------------------------------------------------------
	 Loop over each L1b file to process.
	 ----------------------------------------------------------------------------*/

	ifile = 0;

	for (;;) {

		/*----------------------------------------------------------------------------
		 Grab the current Level 1b radiance filename.
		 ----------------------------------------------------------------------------*/

		strcpy(imgr1.l1b_dir_name, opt.l1b_dir_name);
		if (opt.single_file == YES) {
			strcpy(imgr1.l1b_filename1, cl.l1b_file_name);
		} else {
			fscanf(fptr, "%s", opt.l1b_file_name);
			fgets(comment, MAX_STR_LEN, fptr);
			strcpy(imgr1.l1b_filename1, opt.l1b_file_name);
			if (feof(fptr))
				break;
		}

		ifile++;
		fprintf(stdout,"\n%sProcessing granule #%d\n",EXE_PROMPT,ifile);
		if (opt.verbose == YES)
			fprintf(stdout,"%sfilename = %s\n",EXE_PROMPT,imgr1.l1b_filename1);

		/*----------------------------------------------------------------------------
		 Select the appropriate read routine.
		 ----------------------------------------------------------------------------*/

		printf("Reading in imager data based on opt.data_source...\n");

		switch (opt.data_source) {

			/*----------------------------------------------------------------------------
			 Read in MODIS 1 km data, if used.
			 ----------------------------------------------------------------------------*/

			/*if ((imgr1.landmask = (unsigned char *) malloc(imgr1.npts*/
					/**sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1.landmask");*/

			case MODIS_DAT:
				if ((imgr1.landmask = (unsigned char *) malloc(imgr1.npts
						*sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1.landmask");
				read_modis(imgr1.chflg, &imgr1, opt.geo_interp, opt.verbose);
				break;

			/*----------------------------------------------------------------------------
			 Read in MODIS 5 km data, if used.
			 ----------------------------------------------------------------------------*/

			case MODIS5_DAT:
				/*read_modis5(imgr1.chflg, &imgr1, opt.verbose);*/
				break;

			/*----------------------------------------------------------------------------
			 Read in MODIS 1 km data and degrade spatial resolution, if used.
			 ----------------------------------------------------------------------------*/

			case MODIS_AVGX_DAT:
				/*read_modis_avgx(imgr1.chflg, opt.modis_avgx_res, &imgr1, opt.verbose);*/
				break;

			/*----------------------------------------------------------------------------
			 Read in VIIRS data, if used.
			 ----------------------------------------------------------------------------*/

			case VIIRS_DAT:
				printf("\n\nCalling read_viirs()\n\n"); 
				read_viirs(imgr1.chflg, &imgr1, opt.verbose);
				break;

		}

		if(imgr1.proc_250) {
		  strcpy(imgr1.qkm_dir_name, cl.qkm_dir_name);
		  strcpy(imgr1.qkm_filename, cl.qkm_file_name);
		  printf("\n\nCalling read_modisqkm() %s %s \n\n", imgr1.qkm_dir_name, imgr1.qkm_filename);
		  read_modisqkm(imgr1.chflg_qkm, &imgr1, opt.verbose);
		}
  

		/*----------------------------------------------------------------------------
		 Based on the spacecraft ID, grab the proper Planck routines.
		 ----------------------------------------------------------------------------*/

		switch (imgr1.satid) {
		case TERRA:
			printf("Getting TERRA Planck routines\n");
			printf("shifted Planck routines: %d \n", opt.terra_plod_spectral_shift);

                        if(opt.terra_plod_spectral_shift == YES) {
			  rutil.planck_rad_ptr = modis_planck_terra_shf;
			  rutil.planck_bt_ptr = planck_btemp_terra_shf;
                          rutil.rad_v2w_ptr = rad_v2w_terra_shf;
                        }
                        else {
			  rutil.planck_rad_ptr = modis_planck_terra;
			  rutil.planck_bt_ptr = planck_btemp_terra;
			  rutil.rad_v2w_ptr = rad_v2w_terra;
                        }
			for (ichan=0; ichan<NMODIS_CHAN; ichan++)
				rutil.NEDR[ichan] = NEDR_terra[ichan];
			break;
		case AQUA:
			printf("Getting AQUA Planck routines\n");
			printf("shifted Planck routines: %d \n", opt.aqua_plod_spectral_shift);

                        if(opt.aqua_plod_spectral_shift == YES) {
			  rutil.planck_rad_ptr = modis_planck_aqua_shf;
			  rutil.planck_bt_ptr = planck_btemp_aqua_shf;
                          rutil.rad_v2w_ptr = rad_v2w_aqua_shf;
                        }
                        else {
			  rutil.planck_rad_ptr = modis_planck_aqua;
			  rutil.planck_bt_ptr = planck_btemp_aqua;
			  rutil.rad_v2w_ptr = rad_v2w_aqua;
                        }
			for (ichan=0; ichan<NMODIS_CHAN; ichan++)
				rutil.NEDR[ichan] = NEDR_aqua[ichan];
			break;
		case NPP:
			printf("Getting NPP Planck routines\n"); 

			/*rutil.planck_rad_ptr = viirs_planck_npp;*/
			/*rutil.planck_bt_ptr = planck_btemp_npp;*/
			/*rutil.rad_v2w_ptr = rad_v2w_npp;*/

			// TODO : Do we need NEDR for NPP ?
			/*for (ichan=0; ichan<NVIIRS_CHAN; ichan++)*/
				/*rutil.NEDR[ichan] = NEDR_npp[ichan];*/
			break;
		case NPOESS:
			printf("Getting NPOESS Planck routines\n"); // TODO

			break;
		default:
			fprintf(stderr,"%sUnknown satellite imager - aborting\n",EXE_PROMPT);
			exit(EXIT_FAILURE);
			break;
		}

		/*----------------------------------------------------------------------------
		 Load the Planck tables into memory to allow for fast calculations.
		 ----------------------------------------------------------------------------*/

		if (opt.data_source == (MODIS_DAT || MODIS5_DAT || MODIS_AVGX_DAT)) {
			printf("Loading fast Planck routines\n");(void)fflush(stdout);
			load_fast_planck(&rutil);
		}

		/*----------------------------------------------------------------------------
		 Copy some dimension information.
		 ----------------------------------------------------------------------------*/

		imgr2.ncol = imgr1.ncol;
		imgr2.nrow = imgr1.nrow;
		imgr2.npts = imgr1.npts;

		/*----------------------------------------------------------------------------
		 Determine the pixel level surface type from database.
		 ----------------------------------------------------------------------------*/

		switch (opt.data_source) {
			case MODIS_DAT:
				if (opt.verbose == YES)
					fprintf(stdout,"%sConverting surface type from grid to pixel level.\n",
							EXE_PROMPT);

                                imgr1.sfc_type = IGBP_grid2pixel(imgr1.lat, imgr1.lon, imgr1.npts, lmask_map);
                                free(lmask_map);
				break;

			case VIIRS_DAT:
				// For VIIRS we get this handled in read_viirs().
				break;
		}

		/*----------------------------------------------------------------------------
		 Determine the pixel level surface elevation from database.
		 ----------------------------------------------------------------------------*/

		if (opt.verbose == YES)
			fprintf(
					stdout,
					"%sConverting surface elevation from grid to pixel level.\n",
					EXE_PROMPT);
		switch (opt.data_source) {
			case MODIS_DAT:
				imgr1.zsfc = elev_grid2pixel(imgr1.lat, imgr1.lon, imgr1.npts, elev_map);
				break;
			case VIIRS_DAT:
				// For VIIRS we get this from geolocation in read_viirs().
				break;
		}

		/*----------------------------------------------------------------------------
		 Determine the pixel level sst from database.
		 ----------------------------------------------------------------------------*/

		if (opt.verbose == YES)
			fprintf(stdout,"%sConverting sst from grid to pixel level.\n",EXE_PROMPT);
		imgr1.sst = sst_grid2pixel(imgr1.lat, imgr1.lon, imgr1.npts,
				imgr1.month, sst_map);

		/*----------------------------------------------------------------------------
		 Gather RTM parameters.
		 ----------------------------------------------------------------------------*/

		if (nwp_needed == YES) {

			/*----------------------------------------------------------------------------
			 Read in NWP data, if needed.
			 ----------------------------------------------------------------------------*/

			switch (opt.nwp_source) {
			case GFS_NWP:
				if (opt.verbose == YES)
					fprintf(stdout,"\n%sReading in GFS data.\n",EXE_PROMPT);
				main_gfs(opt.nwp_dir_name, imgr1.year, imgr1.jday, imgr1.time,
						imgr1.bounds, &nwp);
				break;
			case GDAS_NWP:
				if (opt.verbose == YES)
					fprintf(stdout,"\n%sReading in GDAS data.\n",EXE_PROMPT);
				main_gdas(opt.nwp_dir_name, imgr1.year, imgr1.jday, imgr1.time,
						imgr1.bounds, &nwp);
				break;
			case NCEP_NWP:
				if (opt.verbose == YES)
					fprintf(stdout,"\n%sReading in NCEP-Reanalysis data.\n",EXE_PROMPT);
				/*main_ncep(opt.nwp_dir_name, imgr1.year, imgr1.jday, imgr1.time, imgr1.bounds, &nwp);*/
				break;
			case ECMWF_NWP:
				if (opt.verbose == YES)
					fprintf(stdout,"\n%sReading in ECMWF data.\n",EXE_PROMPT);
				/*main_ecmwf(opt.nwp_dir_name, imgr1.year, imgr1.jday, imgr1.time, imgr1.bounds, &nwp);*/
				break;
			}

			/*----------------------------------------------------------------------------
			 Determine the matching nwp grid point for each pixel and determine the
			 viewing zenith and solar zenith angle associated with each nwp grid cell.
			 ----------------------------------------------------------------------------*/

			if (opt.verbose == YES)
				fprintf(
						stdout,
						"%sMapping NWP model output to satellite observations\n",
						EXE_PROMPT);
			nwp_grid2pixel(&nwp, &imgr1);

		}

		if (rtm_needed == YES) {
			printf("Using Radiative Transfer Model...\n");

			/*----------------------------------------------------------------------------
			 Read in surface emissivity from database.
			 ----------------------------------------------------------------------------*/

			if (opt.verbose == YES)
				fprintf(stdout,"%sReading in Seebor surface emissivity.\n",EXE_PROMPT);
			read_seebor_emiss_imager(&imgr1, &rclr);

			/*read_constant_emiss_imager(&imgr1, &rclr);*/

			/*----------------------------------------------------------------------------
			 Allocate memory for RTM output structures.
			 ----------------------------------------------------------------------------*/

			if ((rtm = (rtm_profiles **) malloc(sizeof(rtm_profiles *)
					* nwp.ncells)) == NULL)error_allo(rout,"rtm**");

			for (i=0; i<nwp.ncells; i++) {
				if ((rtm[i] = (rtm_profiles *) malloc(sizeof(rtm_profiles)
						* nwp.rtm_nvzen)) == NULL)error_allo(rout,"rtm *");
				for (ii=0; ii<nwp.rtm_nvzen; ii++)
					rtm[i][ii].flag = NO;
			}

			create_rclr_ptrs(imgr1.chflg, &rclr);

			/*----------------------------------------------------------------------------
			 Perform all required RTM calculations.
			 ----------------------------------------------------------------------------*/

			if (opt.verbose == YES)
				fprintf(stdout,"%sRunning RTM N_CHAN=%d imgPts=%ld\n",	
						EXE_PROMPT,nchan,imgr1.npts);
						/*EXE_PROMPT,NMODIS_CHAN,imgr1.npts);*/

                        for (i=0; i<nchan; i++)
  		            surface_emissivity[i] = 0.990;

                        for (i=0; i<imgr1.npts; i++) {

                           if(imgr1.index_nwp[i] != -999) {

				icell = imgr1.index_nwp[i];
				ivza = imgr1.index_vza[i];
				iclr = imgr1.index_rclr[i];

				if (rtm[icell][ivza].flag == NO) {
					satzen_mid = vza_mid(imgr1.cos_satzen[i]);
					rtm_function_ptr(&rtm[icell][ivza], &nwp.map[icell], rutil,
							nwp.nlevels, surface_emissivity, satzen_mid,
							imgr1.chflg, (int)imgr1.satid, imgr1.algData_dir_name,
                                                        &imgr1.year, &imgr1.jday, &rclr, iclr);

                                        rtm[icell][ivza].flag = YES;
					nwp_rclr_flag = NO;


				}

                                if (rclr.flag[iclr] == NO || rclr.flag[iclr] == YES) {
					get_clear_toa_rad(rtm[icell][ivza], nwp.map[icell], rutil,
							nwp.nlevels, &rclr, imgr1.chflg, iclr);
					rclr.flag[iclr] = YES;
					nwp_rclr_flag = YES;
				}   
                          }

			}

			if (opt.verbose == YES)
				fprintf(stdout,"%sRunning RTM-done\n",EXE_PROMPT);(void)fflush(stdout);
		}

		/*----------------------------------------------------------------------------
		 If writing a Level 2 file, open it and write global attributes and other
		 relevant SDS's.
		 ----------------------------------------------------------------------------*/

		if (opt.make_level2 == YES) {

			if (l2_name_set==0) {
				sprintf(level2_filename,
						"%s/leocatL2.%s_%s.%d%3.3d.%2.2d%2.2d.hdf",
						opt.l2_dir_name, imgr1.satname, opt.data_source_name,
						imgr1.year, imgr1.jday, imgr1.hour, imgr1.minute);
                        }else if(l2_name_set==1){
                          strcpy(level2_filenamebuff, level2_filename);
                          strcpy(level2_filename, "");
                          sprintf(level2_filename, "%s%s%s", opt.l2_dir_name, "/" ,level2_filenamebuff);
                        }

			if (opt.verbose == YES)
				fprintf(stdout,"%sLevel 2 output filename: %s\n", EXE_PROMPT,
						level2_filename);

			l2_sd_id = SDstart(level2_filename, DFACC_CREATE);
			if (l2_sd_id == FAIL) {
				fprintf(
						stderr,
						"%sCannot create Level 2 output file:\n%s%s - aborting\n",
						EXE_PROMPT, EXE_PROMPT,level2_filename);
				exit(EXIT_FAILURE);
			}
			if (opt.verbose == YES)
				fprintf(
						stdout,"%sWriting default output (l2_sd_id = %d) to Level 2 file...\n",
						EXE_PROMPT,l2_sd_id); //GPC

			// Output non-product specific default data
			create_level2_default_output(l2_sd_id, imgr1, opt.nwp_source_name);
			printf("\nCreated defult output to Level 2 file... "); //GPC
			(void)fflush(stdout);

			if (opt.rtm_out == YES && rtm_needed == YES) {
				if (opt.verbose == YES)
					fprintf(stdout,"%sWriting RTM output to Level 2 file...\n",
							EXE_PROMPT);
				/*create_level2_rtm_output_modis(l2_sd_id, imgr1, nwp, rtm, imgr1.chflg);*/
			}
			if (opt.nwp_out == YES && nwp_needed == YES) {
				if (opt.verbose == YES)
					fprintf(stdout,"%sWriting NWP output to Level 2 file...\n",
							EXE_PROMPT);
				create_level2_nwp_output(l2_sd_id, opt.nwp_source, nwp, imgr1);
			}
		}

		/*----------------------------------------------------------------------------
		 If writing a Level 3 file, open it and write global attributes and other
		 relevant SDS's.
		 ----------------------------------------------------------------------------*/

		if (opt.make_level3 == YES) {
			sprintf(level3_filename,
					"%s/leocatL3.%s_%s.%d%3.3d.%2.2d%2.2d.hdf",
					opt.l3_dir_name, imgr1.satname, opt.data_source_name,
					imgr1.year, imgr1.jday, imgr1.hour, imgr1.minute);
			if (opt.verbose == YES)
				fprintf(stdout,"%sLevel 3 output filename: %s\n",
						EXE_PROMPT,level3_filename);
			l3_sd_id = SDstart(level3_filename, DFACC_CREATE);
			if (l3_sd_id == FAIL) {
				fprintf(
						stderr,"%sCannot create Level 3 output file:\n%s%s - aborting\n",
						EXE_PROMPT,EXE_PROMPT,level3_filename);
				exit(EXIT_FAILURE);
			}
		}

		/*----------------------------------------------------------------------------
		 Get default cloud mask.
		 ----------------------------------------------------------------------------*/

		printf("Checking if default cloud mask needed ... "); //GPC
		if (default_cmask_needed == YES) {
			printf("YES\n"); //GPC
			/*----------------------------------------------------------------------------
			 Call the default cloud mask subroutine.
			 ----------------------------------------------------------------------------*/

			if (opt.verbose == YES) {
				fprintf(stdout,"\n%s/***************************************************************/\n",EXE_PROMPT);
				fprintf(stdout,"%sGenerating default cloud mask, %s...\n",EXE_PROMPT,algo_info[opt.cmask_source-1].algo_name);
				fprintf(stdout,"%s/***************************************************************/\n",EXE_PROMPT);
			}

			create_algorithm_pointers(algo_info[opt.cmask_source-1].nout,
					algo_info[opt.cmask_source-1].sds_name, &imgr2);
			if ((imgr1.cldmask = (unsigned char *) malloc(imgr1.npts
					*sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1.cldmask");

			/*----------------------------------------------------------------------------
			 The pixel level land mask must first be initialized using the grid-based
			 pixel land mask in case the default cloud mask does produce a land mask.
			 ----------------------------------------------------------------------------*/


			for (i=0; i<imgr1.npts; i++) {
				if (imgr1.sfc_type[i] == WATER_IGBP)
					imgr2.landmask[i] = WATER_PLM;
				else if (imgr1.sfc_type[i] == DESERT_IGBP)
					imgr2.landmask[i] = DESERT_PLM;
				else
					imgr2.landmask[i] = LAND_PLM;
			}

			start_time_algo = time(NULL);
			algo_info[opt.cmask_source-1].algo_function_ptr(opt.verbose,
					&imgr1, &imgr2, &sndr1, &sndr2, nwp, rtm, rclr, rutil);

			/*----------------------------------------------------------------------------
			 Determine total processing time for each algorithm.
			 ----------------------------------------------------------------------------*/

			end_time_algo = time(NULL);
			tdiff_algo = difftime(end_time_algo, start_time_algo);
			time_elapsed(tdiff_algo, &dday_algo, &dhr_algo, &dmin_algo,
					&dsec_algo);
			print_time_elapsed(dday_algo, dhr_algo, dmin_algo, dsec_algo, 0);

			/*----------------------------------------------------------------------------
			 Copy the cloud mask to another structure.
			 ----------------------------------------------------------------------------*/

			for (i=0; i<imgr1.npts; i++) {
				imgr1.cldmask[i] = imgr2.cldmask[i];
			}

			/*----------------------------------------------------------------------------
			 Write cloud mask to Level 2 file, if needed.
			 ----------------------------------------------------------------------------*/

			if (opt.make_level2 == YES) {
				printf("\n\nWriting default cloudmask to HDF\n\n");
				create_level2_output(l2_sd_id,
						algo_info[opt.cmask_source-1].algo_num,
						algo_info[opt.cmask_source-1].nout,
						algo_info[opt.cmask_source-1].keyword,
						algo_info[opt.cmask_source-1].reference,
						algo_info[opt.cmask_source-1].sds_name, imgr2);
			}

			/*----------------------------------------------------------------------------
			 Write cloud mask to Level 3 file, if needed.
			 ----------------------------------------------------------------------------*/

			/*----------------------------------------------------------------------------
			 If MODIS data are averaged, compute cloud fraction.
			 ----------------------------------------------------------------------------*/

			if (opt.data_source == MODIS_AVGX_DAT && opt.cmask_source
					== MODIS35_CMASK)
				;

			/*----------------------------------------------------------------------------
			 Write cloud fraction to Level 2 file, if needed.
			 ----------------------------------------------------------------------------*/

			/*----------------------------------------------------------------------------
			 Copy the cloud fraction to another structure.
			 ----------------------------------------------------------------------------*/

			destroy_algorithm_pointers(algo_info[opt.cmask_source-1].nout,
					algo_info[opt.cmask_source-1].sds_name, &imgr2);

		}
		if (default_cmask_needed != YES) printf("NO\n"); //GPC

		/*----------------------------------------------------------------------------
		 Get default cloud type, if needed.
		 ----------------------------------------------------------------------------*/

		printf("Checking if default cloud type needed ... ");
		if (default_ctype_needed == YES) {

			/*----------------------------------------------------------------------------
			 Call the default cloud type subroutine.
			 ----------------------------------------------------------------------------*/

			if (opt.verbose == YES) {
				fprintf(stdout,"\n%s/***************************************************************/\n",EXE_PROMPT);
				fprintf(stdout,"%sGenerating default cloud type, %s...\n",EXE_PROMPT,algo_info[opt.ctype_source-1].algo_name);
				fprintf(stdout,"%s/***************************************************************/\n",EXE_PROMPT);
			}

			create_algorithm_pointers(algo_info[opt.ctype_source-1].nout,
					algo_info[opt.ctype_source-1].sds_name, &imgr2);
			if ((imgr1.cldtype = (unsigned char *) malloc(imgr1.npts
					*sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1.cldtype");
			if ((imgr1.cldphase = (unsigned char *) malloc(imgr1.npts
					*sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1.cldphase");

			start_time_algo = time(NULL);
			algo_info[opt.ctype_source-1].algo_function_ptr(opt.verbose,
					&imgr1, &imgr2, &sndr1, &sndr2, nwp, rtm, rclr, rutil);

			/*----------------------------------------------------------------------------
			 Determine total processing time for each algorithm.
			 ----------------------------------------------------------------------------*/

			end_time_algo = time(NULL);
			tdiff_algo = difftime(end_time_algo, start_time_algo);
			time_elapsed(tdiff_algo, &dday_algo, &dhr_algo, &dmin_algo,
					&dsec_algo);
			print_time_elapsed(dday_algo, dhr_algo, dmin_algo, dsec_algo, 0);

			/*----------------------------------------------------------------------------
			 Copy the cloud type to another structure.
			 ----------------------------------------------------------------------------*/

			for (i=0; i<imgr1.npts; i++) {
				imgr1.cldtype[i] = imgr2.cldtype[i];
				imgr1.cldphase[i] = imgr2.cldphase[i];
			}

			/*----------------------------------------------------------------------------
			 Write cloud type to Level 2 file, if needed.
			 ----------------------------------------------------------------------------*/

			if (opt.make_level2 == YES) {
				create_level2_output(l2_sd_id,
						algo_info[opt.ctype_source-1].algo_num,
						algo_info[opt.ctype_source-1].nout,
						algo_info[opt.ctype_source-1].keyword,
						algo_info[opt.ctype_source-1].reference,
						algo_info[opt.ctype_source-1].sds_name, imgr2);
			}

			/*----------------------------------------------------------------------------
			 Write cloud type to Level 3 file, if needed.
			 ----------------------------------------------------------------------------*/

			destroy_algorithm_pointers(algo_info[opt.ctype_source-1].nout,
					algo_info[opt.ctype_source-1].sds_name, &imgr2);
		}

		/*----------------------------------------------------------------------------
		 Get default aerosol mask.
		 ----------------------------------------------------------------------------*/

		printf("Checking if default aerosol mask needed ... "); 
		if (default_amask_needed == YES) {

			/*----------------------------------------------------------------------------
			 Call the default cloud mask subroutine.
			 ----------------------------------------------------------------------------*/

			if (opt.verbose == YES) {
				fprintf(stdout,"\n%s/***************************************************************/\n",EXE_PROMPT);
				fprintf(stdout,"%sGenerating default aerosol mask, %s...\n",EXE_PROMPT,algo_info[opt.amask_source-1].algo_name);
				fprintf(stdout,"%s/***************************************************************/\n",EXE_PROMPT);
			}

			create_algorithm_pointers(algo_info[opt.amask_source-1].nout,
					algo_info[opt.amask_source-1].sds_name, &imgr2);
			if ((imgr1.aeromask = (unsigned char *) malloc(imgr1.npts
					*sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1.aeromask");
			start_time_algo = time(NULL);
			algo_info[opt.amask_source-1].algo_function_ptr(opt.verbose,
					&imgr1, &imgr2, &sndr1, &sndr2, nwp, rtm, rclr, rutil);

			/*----------------------------------------------------------------------------
			 Determine total processing time for each algorithm.
			 ----------------------------------------------------------------------------*/

			end_time_algo = time(NULL);
			tdiff_algo = difftime(end_time_algo, start_time_algo);
			time_elapsed(tdiff_algo, &dday_algo, &dhr_algo, &dmin_algo,
					&dsec_algo);
			print_time_elapsed(dday_algo, dhr_algo, dmin_algo, dsec_algo, 0);

			/*----------------------------------------------------------------------------
			 Copy the aerosol mask to another structure.
			 ----------------------------------------------------------------------------*/

			for (i=0; i<imgr1.npts; i++)
				imgr1.aeromask[i] = imgr2.aeromask[i];

			/*----------------------------------------------------------------------------
			 Write aerosol mask to Level 2 file, if needed.
			 ----------------------------------------------------------------------------*/

			if (opt.make_level2 == YES) {
				create_level2_output(l2_sd_id,
						algo_info[opt.amask_source-1].algo_num,
						algo_info[opt.amask_source-1].nout,
						algo_info[opt.amask_source-1].keyword,
						algo_info[opt.amask_source-1].reference,
						algo_info[opt.amask_source-1].sds_name, imgr2);
			}

			/*----------------------------------------------------------------------------
			 Write aerosol mask to Level 3 file, if needed.
			 ----------------------------------------------------------------------------*/

			/*----------------------------------------------------------------------------
			 If MODIS data are averaged, compute aerosol fraction.
			 ----------------------------------------------------------------------------*/

			if (opt.data_source == MODIS_AVGX_DAT)
				;

			/*----------------------------------------------------------------------------
			 Write aerosol fraction to Level 2 file, if needed.
			 ----------------------------------------------------------------------------*/

			/*----------------------------------------------------------------------------
			 Copy the aerosol fraction to another structure.
			 ----------------------------------------------------------------------------*/

			destroy_algorithm_pointers(algo_info[opt.amask_source-1].nout,
					algo_info[opt.amask_source-1].sds_name, &imgr2);

		}
		if (default_amask_needed != YES) printf("NO\n"); //GPC

		/*----------------------------------------------------------------------------
		 Run through all other cloud algorithms.
		 ----------------------------------------------------------------------------*/

		int algIndex = 0;
		while (algo_order[algIndex] != 0) {
			i = algo_order[algIndex]-1;
			algIndex++;

			/*----------------------------------------------------------------------------
			 Process only desired algorithms.
			 ----------------------------------------------------------------------------*/

			if (opt.algo_flg[i] == YES) {
				/*
				 if (algo_info[i].algo_function_info_ptr != NULL )
				 {
				 algo_info[i].algo_function_info_ptr(algo_info[i].algo_name,algo_info[i].keyword);
				 }
				 */

				if (opt.verbose == YES) {
					fprintf(stdout,"\n%s/***************************************************************/\n",EXE_PROMPT);
					fprintf(stdout,"%sProcessing %s algorithm...\n",EXE_PROMPT,algo_info[i].algo_name);
					fprintf(stdout,"%s/***************************************************************/\n",EXE_PROMPT);
				}

				/*----------------------------------------------------------------------------
				 Create algorithm pointers.
				 ----------------------------------------------------------------------------*/

				if (opt.verbose == YES) {
					fprintf(stdout,"%sCreating algorithm pointers\n",EXE_PROMPT);	//GPC
					fprintf(stdout,"\talgo_info[%d].nout=%d\n",i,algo_info[i].nout);//GPC
					int alg_SDS;
					for (alg_SDS=0; alg_SDS<algo_info[i].nout; alg_SDS++){
						fprintf(stdout,"\talgo_info[%d].sds_name[%d]=%s\n",i,alg_SDS,algo_info[i].sds_name[alg_SDS]);//GPC
					}
				}
				create_algorithm_pointers(algo_info[i].nout,
						algo_info[i].sds_name, &imgr2);
				//GPC
//				if (opt.verbose == YES) {
//					fprintf(stdout,"Part C %d %s\n",i, algo_info[i].algo_name);
//				}

				/*----------------------------------------------------------------------------
				 Call the algorithm subroutine.
				 ----------------------------------------------------------------------------*/

                                printf("Calling %s subroutine\n",algo_info[i].algo_name);

				// Initialise structure containing channels
//                              printf("\tInitialising %d  channels\n",nchan);
				for (ichan=0; ichan<nchan; ichan++) {
					imgr1.chflg_rtm[ichan] = algo_info[i].chflg[ichan];
				}

				start_time_algo = time(NULL);
                                printf("\tCalling pointer to %s\n",algo_info[i].algo_name);
				algo_info[i].algo_function_ptr(opt.verbose, &imgr1, &imgr2,
						&sndr1, &sndr2, nwp, rtm, rclr, rutil);

				printf("\n\nReturned from Algorithm %s\n",algo_info[i].algo_name);
                                (void)fflush(stdout);

				/* if provided by algorithm ask for name and version */
//				if (opt.verbose == YES) {
//					fprintf(stdout,"Part A\n");
//				}

				/*----------------------------------------------------------------------------
				 Determine total processing time for each algorithm.
				 ----------------------------------------------------------------------------*/

				end_time_algo = time(NULL);
				tdiff_algo = difftime(end_time_algo, start_time_algo);
				time_elapsed(tdiff_algo, &dday_algo, &dhr_algo, &dmin_algo,
						&dsec_algo);
				print_time_elapsed(dday_algo, dhr_algo, dmin_algo, dsec_algo, 0);

				/*----------------------------------------------------------------------------
				 Write algorithm results to a Level 2 file, if needed.
				 ----------------------------------------------------------------------------*/

				printf("\nAlgorithm Level2 output...");//GPC
				if (opt.make_level2 == YES) {
					printf("YES\n");//GPC
					if (opt.verbose == YES)
						fprintf(stdout,"%sWriting %s algorithm results to Level 2 file...\nFile=%s\n",
						EXE_PROMPT,algo_info[i].algo_name,level2_filename);

					printf("Calling create_level2_output() for l2_sd_id = %d\n",l2_sd_id);//GPC
                                        (void)fflush(stdout);
					create_level2_output(l2_sd_id, algo_info[i].algo_num,
							algo_info[i].nout, algo_info[i].keyword,
							algo_info[i].reference, algo_info[i].sds_name,
							imgr2);
					printf("\nFinished writing to HDF file with l2_sd_id = %d\n",l2_sd_id);//GPC
				}
				if (opt.make_level2 != YES) printf("NO\n");//GPC

				/*----------------------------------------------------------------------------
				 Create Level 3 output from algorithm results and write to a Level 3 file,
				 if needed.
				 ----------------------------------------------------------------------------*/

				if (opt.make_level3 == YES) {
					if (opt.verbose == YES)
						fprintf(stdout,"%sWriting %s algorithm results to Level 3 file...\n",EXE_PROMPT,algo_info[i].algo_name);
				}

				/*----------------------------------------------------------------------------
				 Destroy algorithm pointers.
				 ----------------------------------------------------------------------------*/

				if (opt.verbose == YES)
					fprintf(stdout,"%sDestroying algorithm pointers\n",EXE_PROMPT);
				destroy_algorithm_pointers(algo_info[i].nout, algo_info[i].sds_name, &imgr2);
			}

		}

		/*----------------------------------------------------------------------------
		 If a Level 2 file was created, close it.
		 ----------------------------------------------------------------------------*/

		if (opt.make_level2 == YES) {

			status = SDend(l2_sd_id);
			printf("Closing L2 file with id %d returns status %d\n",l2_sd_id,status);
			if (status == FAIL) {
				HEprint(stdout, 0);
				fprintf(stderr,"%sCannot close Level 2 output file:\n%s%s - aborting\n",EXE_PROMPT,EXE_PROMPT,level2_filename);
				exit(EXIT_FAILURE);
			}
		}

		/*----------------------------------------------------------------------------
		 If a Level 3 file was created, close it.
		 ----------------------------------------------------------------------------*/

		if (opt.make_level3 == YES) {
			status = SDend(l3_sd_id);
			if (status == FAIL) {
				fprintf(stderr,"%sCannot close Level 3 output file:\n%s%s - aborting\n",EXE_PROMPT,EXE_PROMPT,level3_filename);
				exit(EXIT_FAILURE);
			}
		}

		/*----------------------------------------------------------------------------
		 Destroy satellite data pointers.
		 ----------------------------------------------------------------------------*/

		if (opt.verbose == YES)
			fprintf(stdout,"\n%sDestroying satellite data pointers\n",EXE_PROMPT);
		destroy_l1b_ptrs(imgr1.chflg, &imgr1);
		printf("Finished destroying satellite data pointers\n");

		printf("Destroying fast planck pointers\n");
		if (opt.data_source == (MODIS_DAT || MODIS5_DAT || MODIS_AVGX_DAT)) {
			destroy_fast_planck(&rutil);
		}
		printf("Finished destroying fast planck pointers\n");

		/*----------------------------------------------------------------------------
		 Destroy default cloud/aerosol mask and type pointers.
		 ----------------------------------------------------------------------------*/

		if (default_cmask_needed == YES) {
			free(imgr1.cldmask);
			free(imgr1.landmask);
		}
		if (default_amask_needed == YES)
			free(imgr1.aeromask);
		if (default_ctype_needed == YES)
			free(imgr1.cldtype);

		if (rtm_needed == YES || nwp_needed == YES) {

			/*----------------------------------------------------------------------------
			 Destroy nwp pointers, if needed.
			 ----------------------------------------------------------------------------*/

			if (opt.verbose == YES)
				fprintf(stdout,"%sDestroying NWP data pointers\n",EXE_PROMPT);
			for (i=0; i<nwp.ncells; i++) {
				destroy_nwp_pointers(&nwp.map[i]);
			}
			free(nwp.map);

			/*----------------------------------------------------------------------------
			 Destroy rtm pointers, if needed.
			 ----------------------------------------------------------------------------*/

			if (opt.verbose == YES)
				fprintf(stdout,"%sDestroying RTM data pointers\n",EXE_PROMPT);
			if (rtm_needed == YES) {
				for (i=0; i<imgr1.npts; i++) {
                                   if(imgr1.index_nwp[i] != -999) {
					icell = imgr1.index_nwp[i];
					ivza = imgr1.index_vza[i];
					if (rtm[icell][ivza].flag == YES) {
						destroy_rtm_ptrs(&rtm[icell][ivza], nwp, imgr1.chflg);
						rtm[icell][ivza].flag = NO;
					}
                                   }
				}
				for (i=0; i<nwp.ncells; i++)
					free(rtm[i]);
				free(rtm);
				destroy_rclr_ptrs(imgr1.chflg, &rclr);
			}
			free(imgr1.index_vza);
			free(imgr1.index_nwp);
			free(imgr1.index_rclr);
		}

		/*----------------------------------------------------------------------------
		 Exit loop if no more files to process.
		 ----------------------------------------------------------------------------*/

		if (opt.single_file == YES) {
			break;
		}
	}

	/*----------------------------------------------------------------------------
	 Close file list file, if open.
	 ----------------------------------------------------------------------------*/

	if (opt.single_file == NO)
		fclose(fptr);

	/*----------------------------------------------------------------------------
	 Deallocate some remaining memory.
	 ----------------------------------------------------------------------------*/

	free(elev_map);
	free(sst_map);
	free(imgr1.chflg);
	free(imgr1.chflg_rtm);
	free(surface_emissivity);
	free(rutil.NEDR);

	/*----------------------------------------------------------------------------
	 Print finishing statement.
	 ----------------------------------------------------------------------------*/

	fprintf(stdout,"%sFinished with processing.\n\n",EXE_PROMPT);

	/*----------------------------------------------------------------------------
	 Determine total processing time.
	 ----------------------------------------------------------------------------*/

	end_time_total = time(NULL);
	tdiff_total = difftime(end_time_total, start_time_total);
	time_elapsed(tdiff_total, &dday_total, &dhr_total, &dmin_total, &dsec_total);
	print_time_elapsed(dday_total, dhr_total, dmin_total, dsec_total, 1);

	return (0);
}

/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/******************************************************************************/

void version() {
	fprintf(stdout,"\n  %s, v%s\n",EXE_NAME,PVERSION);
}

void print_alist() {
	int ialgo, iout;
	fprintf(stdout,"\n%sPrinting algorithm information...\n",EXE_PROMPT);
	for (ialgo=0; ialgo<NALGO; ialgo++) {
		fprintf(stdout,"Algorithm index = %d:\n\tname = %s\n\treference = %s\n",
		algo_info[ialgo].algo_num,algo_info[ialgo].algo_name,
		algo_info[ialgo].reference);
		for (iout=0; iout<algo_info[ialgo].nout; iout++) {
			fprintf(stdout,"\tOutput index = %d   Output name = %s_%s\n", iout
					+ 1, algo_info[ialgo].keyword,
					algo_info[ialgo].sds_name[iout]);
		}
		fprintf(stdout,"\n");
	}
}

/******************************************************************************/
/******************************************************************************/

void help() {
	fprintf(stdout,"\n%shelp: option list\n\n",EXE_PROMPT);

	fprintf(stdout,"  -a[1, 2, ...NALGO]\n");
	fprintf(stdout,
			"    Use this option to specify which algorithms are executed, not including\n");
	fprintf(
			stdout,
			"    dependencies. Specify the algorithms by number (see -alist to determine algorithm\n");
	fprintf(
			stdout,
			"    number). Multiple algorithm numbers can be specified by leaving a space between each\n");
	fprintf(stdout,"    number.\n\n");

	fprintf(stdout,"  -alist\n");
	fprintf(stdout,
			"    List information about each available algorithm then exit.\n\n");

	fprintf(stdout,"  -amask[volash_4ch]\n");
	fprintf(stdout,
			"    Use this option to specify the default aerosol mask algorithm.\n\n");

	fprintf(stdout,"  -cmask[modis35]\n");
	fprintf(stdout,
			"    Use this option to specify the default cloud mask algorithm.\n\n");

	fprintf(stdout,"  -cmask_dir[directory name]\n");
	fprintf(stdout,
			"    Use this option to specify the directory where the pre-processed\n");
	fprintf(stdout,"    cloud mask file(s) are located.\n\n");

	fprintf(stdout,"  -cprod_dir[directory name]\n");
	fprintf(stdout,
			"    Use this option to specify the directory where the pre-processed\n");
	fprintf(stdout,"    cloud product file(s) are located.\n\n");

	fprintf(stdout,"  -ctype[ph05, modis_ir]\n");
	fprintf(stdout,
			"    Use this option to specify the default cloud type algorithm.\n\n");

	fprintf(stdout,"  -data [modis, modis5, modis_avgx, viirs]\n");
	fprintf(stdout,"    Use this option to specify the imager data source.\n\n");

	fprintf(stdout,"  -data_res [1.0 - ?]\n");
	fprintf(stdout,
			"    If the modis_avgx option is set, set this option to the\n");
	fprintf(stdout,
			"    desired degraded spatial resolution (must be >= 1.0).\n\n");

	fprintf(stdout,"  -f[Level1b filename]\n");
	fprintf(stdout,
			"    Use this option to specify the name of the Level 1B file to process.\n\n");

	fprintf(stdout,"  -flist_name[Level1b filelist name]\n");
	fprintf(
			stdout,
			"    Use this option to specify the name of a text file that has a list of Level\n");
	fprintf(stdout,"    1B file(s) to process.\n\n");

	fprintf(stdout,"  -l1b_dir[directory name]\n");
	fprintf(
			stdout,
			"    Use this option to specify the directory where the Level 1b file(s) are located.\n\n");

	fprintf(stdout,"  -l2_dir[directory name]\n");
	fprintf(
			stdout,
			"    Use this option to specify the directory where the Level 2 output file(s) are to\n");
	fprintf(stdout,"    be written.\n\n");

	fprintf(stdout,"  -l3res[0.125 - 2.5]\n");
	fprintf(
			stdout,
			"    Use this option to specify the Level 3 grid cell resolution in degrees.\n\n");

	fprintf(stdout,"  -l3_dir[directory name]\n");
	fprintf(
			stdout,
			"    Use this option to specify the directory where the Level 3 output file(s) are to\n");
	fprintf(stdout,"    be written.\n\n");

	fprintf(stdout,"  -nol2\n");
	fprintf(stdout,
			"    Set this option if no Level 2 data are to be created.\n\n");

	fprintf(stdout,"  -nol3\n");
	fprintf(stdout,
			"    Set this option if no Level 3 data are to be created.\n\n");

	fprintf(stdout,"  -nwp [gfs, gdas, ncep, ecmwf]\n");
	fprintf(stdout,"    Use this option to specify the nwp data source.\n");
	fprintf(stdout,
			"    This option is only valid if nwp is needed by at least 1 algorithm.\n\n");

	fprintf(stdout,"  -nwp_dir[directory name]\n");
	fprintf(
			stdout,
			"    Use this option to specify the directory where the NWP file(s) are located.\n\n");

	fprintf(stdout,"  -nwp_out\n");
	fprintf(
			stdout,
			"    Set this option if certain NWP parameters are to be added to Level 2 file.\n\n");

	fprintf(stdout,"  -rtm [plod, crtm]\n");
	fprintf(stdout,
			"    Use this option to specify the fast radiative transfer model to use.\n");
	fprintf(
			stdout,
			"    This option is only valid if an RTM is needed by at least 1 algorithm.\n\n");

	fprintf(stdout,"  -rtm_out\n");
	fprintf(
			stdout,
			"    Set this option if certain clear sky radiative transfer calculations are\n");
	fprintf(stdout,"    to be added to Level 2 file.\n\n");

	fprintf(stdout,"  -verbose\n");
	fprintf(stdout,
			"    Output detailed processing information to the standard output device.\n\n");

	fprintf(stdout,"  -version or -v\n");
	fprintf(stdout,"    Display the code version number then exit.\n\n");
}

/******************************************************************************/
/******************************************************************************/

// If this routine is called we are expecting command line input of the form "-switch value".
// If the input does not conform to the above, we try to diagnose the problem, and EXIT.

void check_command_line_input(int iarg, int narg, int string_check, char** arg,
		char* opt) {

	// We've specified a switch without a corresponding value...
	if (iarg == narg) {
		fprintf(
				stderr,
				"%sInvalid command line input option following [-%s],(iarg == narg) displaying help...\n",
				EXE_PROMPT, opt);
		help();
		exit(EXIT_FAILURE);
	}

	// Instead of the corresponding value, our switch is followed by another switch...
	if (arg[iarg][0] == '-') {
		fprintf(
				stderr,
				"%sInvalid command line input option following [-%s],(arg[iarg][0] == '-') displaying help...\n",
				EXE_PROMPT, opt);
		help();
		exit(EXIT_FAILURE);
	}

	if (string_check == YES) {
		if (strlen(arg[iarg]) > MAX_STR_LEN) {
			fprintf(
					stderr,
					"%s String length following input option [-%s] exceeds the maximum length of %d, aborting\n",
					EXE_PROMPT, opt, MAX_STR_LEN);
			exit(EXIT_FAILURE);
		}
	}
}

/******************************************************************************/
/******************************************************************************/

void print_input_options(input_options o, int8 default_cmask_needed,
		int8 default_ctype_needed, int8 default_amask_needed) {
	int i;

	fprintf(stdout,"\n%sPrinting input options...\n", EXE_PROMPT);

	if (o.data_source == MODIS_AVGX_DAT)
		fprintf(stdout,"DATA TYPE:  %s - averaged to: %5.2f km\n",
				o.data_source_name, o.l1b_data_res);
	else
		fprintf(stdout,"DATA TYPE:  %s\n",o.data_source_name);
	fprintf(stdout,"DEFAULT NWP:  %s\n",o.nwp_source_name);
	fprintf(stdout,"DEFAULT RTM:  %s\n",o.rtm_source_name);
	fprintf(stdout,"DEFAULT CLOUD MASK:  %s\n",o.cmask_source_name);
	fprintf(stdout,"DEFAULT CLOUD TYPE:  %s\n",o.ctype_source_name);
	fprintf(stdout,"DEFAULT AEROSOL MASK:  %s\n",o.amask_source_name);
	fprintf(stdout,"LEVEL 1b DATA DIRECTORY:  %s\n",o.l1b_dir_name);
	fprintf(stdout,"GEOLOCATION FILE:  %s\n",o.geo_file);
	fprintf(stdout,"CLOUD MASK DIRECTORY (MAY BE NA):  %s\n",o.cmask_dir_name);
	fprintf(stdout,"CLOUD MASK FILE (MAY BE NA):  %s\n",o.cmask_file_name);
	fprintf(stdout,"CLOUD PRODUCT DIRECTORY (MAY BE NA):  %s\n",o.cprod_dir_name);
	fprintf(stdout,"NWP DIRECTORY (MAY BE NA):  %s\n",o.nwp_dir_name);
	if (o.rtm_out == YES)
		fprintf(stdout,"OUTPUT RTM OUTPUT TO HDF: YES\n");
	else
		fprintf(stdout,"OUTPUT RTM OUTPUT TO HDF: NO\n");
	if (o.nwp_out == YES)
		fprintf(stdout,"OUTPUT NWP OUTPUT TO HDF: YES\n");
	else
		fprintf(stdout,"OUTPUT NWP OUTPUT TO HDF: NO\n");
	if (o.make_level2 == YES) {
		fprintf(stdout,"MAKE LEVEL 2 PRODUCTS:  YES\n");
		fprintf(stdout,"LEVEL 2 OUTPUT DIRECTORY:  %s\n",o.l2_dir_name);
	} else {
		fprintf(stdout,"MAKE LEVEL 2 PRODUCTS:  NO\n");
		fprintf(stdout,"LEVEL 2 OUTPUT DIRECTORY:  NA\n");
	}
	if (o.make_level3 == YES) {
		fprintf(stdout,"MAKE LEVEL 3 PRODUCTS:  YES\n");
		fprintf(stdout,"LEVEL 3 PRODUCT RESOLUTION:  %f\n",o.level3_res);
		fprintf(stdout,"LEVEL 3 OUTPUT DIRECTORY:  %s\n",o.l3_dir_name);
	} else {
		fprintf(stdout,"MAKE LEVEL 3 PRODUCTS:  NO\n");
		fprintf(stdout,"LEVEL 3 PRODUCT RESOLUTION:  NA\n");
		fprintf(stdout,"LEVEL 3 OUTPUT DIRECTORY:  NA\n");
	}
	fprintf(stdout,"ALGORITHM(S) TO BE PROCESSED, INCLUDING DEFAULT ALGORITHMS:\n");
	if (default_cmask_needed == YES)
		fprintf(stdout,"\tAlgorithm index = %d, Algorithm name = %s\n",
				algo_info[o.cmask_source - 1].algo_num,
				algo_info[o.cmask_source - 1].algo_name);
	if (default_ctype_needed == YES)
		fprintf(stdout,"\tAlgorithm index = %d, Algorithm name = %s\n",
				algo_info[o.ctype_source - 1].algo_num,
				algo_info[o.ctype_source - 1].algo_name);
	if (default_amask_needed == YES)
		fprintf(stdout,"\tAlgorithm index = %d, Algorithm name = %s\n",
				algo_info[o.amask_source - 1].algo_num,
				algo_info[o.amask_source - 1].algo_name);
	for (i=0; i<NALGO; i++) {
		if (o.algo_flg[i] == YES && i != o.cmask_source-1 && i
				!= o.ctype_source-1)
			fprintf(stdout,"\tAlgorithm index = %d, Algorithm name = %s\n",
					algo_info[i].algo_num, algo_info[i].algo_name);
	}
	fprintf(
			stdout,
			"\n%sTotal number of algorithms to process (including default) = %d\n",
			EXE_PROMPT, o.algo2process);

}
