/*$Id: read_satdata.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <glob.h>

#include <hdf.h>
#include <mfhdf.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <common_leocat.h>
#include <imagerL1_leocat.h>

void modis_level1b_read(char *, int, int, int, float, int, void **,
		unsigned char *, int *, int *);
void modis_level03_read(char *, int, int, int, void **, unsigned char *, int *,
		int *);
float modis_planck_terra(int, float);
float modis_planck_aqua(int, float);
float * modis_geo_interp_1000(int, int, float *);
void get_geo_filename(int8, int8 *, imagerL1 *);
unsigned char * unpack_detector_flags(int, int, unsigned char *);

void setImgr1Time(char *YYYYDDD, char *HHMM, imagerL1 *imgr1) {
	int ichar=0;
	char year_str[] = { '\0', '\0', '\0', '\0', '\0' };
	char jday_str[] = { '\0','\0', '\0', '\0' };
	char hh_str[] = { '\0', '\0', '\0' };
	char mm_str[] = {'\0', '\0', '\0' };

	if (YYYYDDD != NULL) {
		for (ichar = 0; ichar < 4; ichar++)
			year_str[ichar] = YYYYDDD[ichar];
		imgr1->year = atoi(year_str);

		for (ichar = 0; ichar < 3; ichar++)
			jday_str[ichar] = YYYYDDD[ichar + 4];
		imgr1->jday = atoi(jday_str);
	}

	if (HHMM != NULL) {

		for (ichar = 0; ichar < 2; ichar++)
			hh_str[ichar] = HHMM[ichar];
		imgr1->hour = atoi(hh_str);
		for (ichar = 0; ichar < 2; ichar++)
			mm_str[ichar] = HHMM[2 + ichar];
		imgr1->minute = atoi(mm_str);
		imgr1->time = (float) imgr1->hour + imgr1->minute / 60.0;
	}



	imgr1->ileap = leap_year_check(imgr1->year);
	imgr1->month = jday2month(imgr1->jday, imgr1->ileap);
	imgr1->day = jday2day(imgr1->jday, imgr1->ileap);


}

void read_modis(int8 *chflg, imagerL1 *imgr1, int8 geo_interp, int8 verbose) {

	char *stmp1, *rout = { "read_modis" },
			full_filename[MAX_STR_LEN],
			att_name[MAX_STR_LEN];
	char year_str[] = { '\0', '\0', '\0', '\0', '\0' }, jday_str[] = { '\0',
			'\0', '\0', '\0' }, hh_str[] = { '\0', '\0', '\0' }, mm_str[] = {
			'\0', '\0', '\0' };
	int8 des_node, asc_node, scale_lon_flg;
	int ichar, i, nchan, i_nadir;
	uint16 *ui16junk;
	int32 id, geo_id, status, att_index, data_type, count;
	int ncol5, nrow5, npts5;
	float sin_solzen, sin_satzen, cos_relaz, ch20_ems, minlat, maxlat, minlon,
			maxlon;
	float32 *earth2sun_temp, *f32junk, SOLAR_20;
	void *attr;
	uint8 *detector_qc;

	int offsetForTime = 0;

	/*----------------------------------------------------------------------------
	 Use the filename to determine the satname and if direct broadcast.
	 ----------------------------------------------------------------------------*/

	// SSEC PEATE FILENAME aqua0.MODIS.GEO.2008012.0010.005.2008014.122632.hdf

//      printf("read_modis: %s\n", imgr1->l1b_filename1);
	char *filename = imgr1->l1b_filename1;


	/*

	if ((stmp1 = strstr(filename, "aqua0.MODIS.KM_RAD")) != NULL) {
		strcpy(imgr1->satname, "Aqua");
		imgr1->satid = AQUA;
		imgr1->directbc_flg = NO;
//              printf("found %s\n", stmp1);
		offsetForTime = strlen("aqua0.MODIS.LAC.");

	}

	else if ((stmp1 = strstr(imgr1->l1b_filename1, "MOD021KM.A")) != NULL) {
		strcpy(imgr1->satname, "Terra");
		imgr1->satid = TERRA;
		imgr1->directbc_flg = NO;
		offsetForTime = strlen("MOD021KM.A");

	} else if ((stmp1 = strstr(imgr1->l1b_filename1, "t1.")) != NULL) {
		strcpy(imgr1->satname, "Terra");
		imgr1->satid = TERRA;
		imgr1->directbc_flg = YES;
	} else if ((stmp1 = strstr(imgr1->l1b_filename1, "MYD021KM.A")) != NULL) {
		strcpy(imgr1->satname, "Aqua");
		imgr1->satid = AQUA;
		imgr1->directbc_flg = NO;
		offsetForTime = strlen("MOD021KM.A");

	} else if ((stmp1 = strstr(imgr1->l1b_filename1, "a1.")) != NULL) {
		strcpy(imgr1->satname, "Aqua");
		imgr1->satid = AQUA;
		imgr1->directbc_flg = YES;
	} else {
		fprintf(
				stderr,"%sData file name, %s, does not conform to DAAC or direct broadcast naming standards.\n",
				EXE_PROMPT,imgr1->l1b_filename1);
		fprintf(
				stderr,"%sThe expected naming convention is used to parse information about granule-aborting\n",
				EXE_PROMPT);
		exit(EXIT_FAILURE);
	}

*/
	/*----------------------------------------------------------------------------
	 Use the DAAC filename to determine the year, Julian day, and time.
	 ----------------------------------------------------------------------------*/
/*
//      printf("Using %s\n", stmp1);

	if (imgr1->directbc_flg == NO) {
		for (ichar = 0; ichar < 4; ichar++)
			year_str[ichar] = stmp1[ichar + offsetForTime];
		imgr1->year = atoi(year_str);
		for (ichar = 0; ichar < 3; ichar++)
			jday_str[ichar] = stmp1[ichar + offsetForTime + 4];
		imgr1->jday = atoi(jday_str);
		for (ichar = 0; ichar < 2; ichar++)
			hh_str[ichar] = stmp1[ichar + offsetForTime + 8];
		imgr1->hour = atoi(hh_str);
		for (ichar = 0; ichar < 2; ichar++)
			mm_str[ichar] = stmp1[ichar + offsetForTime + 10];
		imgr1->minute = atoi(mm_str);
		imgr1->time = (float) imgr1->hour + imgr1->minute / 60.0;
	}
*/

	/*----------------------------------------------------------------------------
	 Use the direct broadcast filename to determine the year, Julian day,
	 and time.
	 ----------------------------------------------------------------------------*/
/*
	else {
		if (strstr(imgr1->l1b_filename1, ".1000m.") == NULL) {
			fprintf(
					stderr,"%sDirect broadcast MODIS file, %s, is not 1 km resolution, cannot continue.\n",
					EXE_PROMPT,imgr1->l1b_filename1);
			exit(EXIT_FAILURE);
		}
		for (ichar = 0; ichar < 2; ichar++)
			year_str[ichar] = stmp1[ichar + 3];
		imgr1->year = atoi(year_str);
		imgr1->year = (imgr1->year < 60) ? imgr1->year + 2000 : imgr1->year
				+ 1900;
		for (ichar = 0; ichar < 3; ichar++)
			jday_str[ichar] = stmp1[ichar + 5];
		imgr1->jday = atoi(jday_str);
		for (ichar = 0; ichar < 2; ichar++)
			hh_str[ichar] = stmp1[ichar + 9];
		imgr1->hour = atoi(hh_str);
		for (ichar = 0; ichar < 2; ichar++)
			mm_str[ichar] = stmp1[ichar + 11];
		imgr1->minute = atoi(mm_str);
		imgr1->time = (float) imgr1->hour + imgr1->minute / 60.0;
	}

*/
/*
	imgr1->ileap = leap_year_check(imgr1->year);
	imgr1->month = jday2month(imgr1->jday, imgr1->ileap);
	imgr1->day = jday2day(imgr1->jday, imgr1->ileap);
*/

	if (verbose == YES)
		fprintf(stdout,"%ssatname=%s, year=%d, jday=%d, time=%f\n",
				EXE_PROMPT,imgr1->satname, imgr1->year, imgr1->jday,
				imgr1->time);

	if (geo_interp == NO)
//              printf("Calling get_geo_filename()\n");
		get_geo_filename(verbose, &geo_interp, imgr1);

	/*----------------------------------------------------------------------------
	 Open the hdf file.
	 ----------------------------------------------------------------------------*/

	sprintf(full_filename, "%s/%s", imgr1->l1b_dir_name, imgr1->l1b_filename1);
	id = SDstart(full_filename, DFACC_READ);
	if (id == FAIL) {
		fprintf(stderr,"%sInvalid HDF file, %s\n", EXE_PROMPT,full_filename);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Read in the Earth-Sun Distance global attribute.
	 ----------------------------------------------------------------------------*/

	strcpy(att_name, "Earth-Sun Distance");
	att_index = SDfindattr(id, att_name);
	status = SDattrinfo(id, att_index, att_name, &data_type, &count);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDattrinfo - Earth-Sun Distance.");

	if ((attr = (void *) malloc(count * DFKNTsize(data_type))) == NULL)error_allo(rout,"attr-Earth-Sun Distance");
	status = SDreadattr(id, att_index, attr);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDreadattr - Earth-Sun Distance.");

	earth2sun_temp = (float32 *) attr;
	imgr1->earth2sun_dist = earth2sun_temp[0];
	free(earth2sun_temp);

	/*----------------------------------------------------------------------------
	 Read in detector quality flags.
	 ----------------------------------------------------------------------------*/

	strcpy(att_name, "Detector Quality Flag");
	att_index = SDfindattr(id, att_name);
	status = SDattrinfo(id, att_index, att_name, &data_type, &count);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDattrinfo - Detector Quality Flag.");

	if ((attr = (void *) malloc(count * DFKNTsize(data_type))) == NULL)error_allo(rout,"attr-Detector Quality Flag");
	status = SDreadattr(id, att_index, attr);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDreadattr - Detector Quality Flag.");

	detector_qc = (uint8 *) attr;

	/*----------------------------------------------------------------------------
	 Unpack detector quality flags.
	 ----------------------------------------------------------------------------*/

	imgr1->noise_det = unpack_detector_flags(1, count, detector_qc);
	imgr1->dead_det = unpack_detector_flags(2, count, detector_qc);
	imgr1->gain_det = unpack_detector_flags(3, count, detector_qc);
	imgr1->range_det = unpack_detector_flags(4, count, detector_qc);
	imgr1->source_det = unpack_detector_flags(5, count, detector_qc);
	imgr1->residuals_det = unpack_detector_flags(6, count, detector_qc);
	imgr1->crosstalk_det = unpack_detector_flags(7, count, detector_qc);
	free(detector_qc);

	/*----------------------------------------------------------------------------
	 Read in number of scans.
	 ----------------------------------------------------------------------------*/

	int32 *numScans_temp;
	strcpy(att_name, "Number of Scans");
	att_index = SDfindattr(id, att_name);
	status = SDattrinfo(id, att_index, att_name, &data_type, &count);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDattrinfo - Number_of_Scans.");

	if ((attr = (void *) malloc(count * DFKNTsize(data_type))) == NULL)error_allo(rout,"attr-Number_of_Scans");
	status = SDreadattr(id, att_index, attr);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDreadattr - Number_of_Scans.");

	numScans_temp = (int32 *) attr;
	imgr1->nScans = numScans_temp[0];
	free(numScans_temp);

	printf("\nNumber of scans = %d\n",imgr1->nScans);
	imgr1->nDetectors = 10;
//      printf("\nNumber of detectors = %d\n\n",imgr1->nDetectors);

	/*----------------------------------------------------------------------------
	 Read in the 1000 m dimensions.
	 ----------------------------------------------------------------------------*/

	modis_level1b_read("OPEN_1000", imgr1->satid, id, 1, 1.0, 0,
			(void **) &ui16junk, imgr1->badmask, &imgr1->nrow, &imgr1->ncol);
	free(ui16junk);

	imgr1->npts = (long) imgr1->nrow * (long) imgr1->ncol;

	if ((imgr1->badmask = (unsigned char *) malloc(imgr1->npts
			* sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1->badmask");
	for (i = 0; i < imgr1->npts; i++)
		imgr1->badmask[i] = NO;

	/*----------------------------------------------------------------------------
	 If no geolocation file, then read in 5 km geo from L1b and interpolate.
	 ----------------------------------------------------------------------------*/

	if (geo_interp) {

		/*----------------------------------------------------------------------------
		 Read in the latitude and interpolate to 1000m.
		 ----------------------------------------------------------------------------*/

		if (verbose == YES)
			fprintf(stdout,"\n%sReading in Latitude.\n", EXE_PROMPT);
		modis_level1b_read("OPEN_1000", imgr1->satid, id, 6, 0.0, 0,
				(void **) &f32junk, imgr1->badmask, &nrow5, &ncol5);
		npts5 = ncol5 * nrow5;
		array_minmax_float(npts5, f32junk, &minlat, &maxlat);
		if (verbose == YES)
			fprintf(stdout,"%sInterpolating Latitude.\n", EXE_PROMPT);
		imgr1->lat = modis_geo_interp_1000(ncol5, nrow5, f32junk);
		free(f32junk);

		/*----------------------------------------------------------------------------
		 Read in the longitude and interpolate to 1000m.
		 ----------------------------------------------------------------------------*/

		if (verbose == YES)
			fprintf(stdout,"%sReading in Longitude.\n", EXE_PROMPT);
		modis_level1b_read("OPEN_1000", imgr1->satid, id, 7, 0.0, 0,
				(void **) &f32junk, imgr1->badmask, &nrow5, &ncol5);
		array_minmax_float(npts5, f32junk, &minlon, &maxlon);
		if (verbose == YES)
			fprintf(stdout,"%sInterpolating in Longitude.\n", EXE_PROMPT);

		/*----------------------------------------------------------------------------
		 If not in the polar regions, but possibly near the date line scale the
		 longitude to 0 - 360.
		 ----------------------------------------------------------------------------*/

		if (minlat < 80.0 && maxlat > -80.0 && minlon < -60.0 && maxlon > 60.0) {
			for (i = 0; i < npts5; i++) {
				if (f32junk[i] < 0.0 && f32junk[i] > -180.0)
					f32junk[i] += 360.0;
			}
		}
		imgr1->lon = modis_geo_interp_1000(ncol5, nrow5, f32junk);
		free(f32junk);
		array_minmax_float(imgr1->npts, imgr1->lat, &minlat, &maxlat);
		array_minmax_float(imgr1->npts, imgr1->lon, &minlon, &maxlon);
		for (i = 0; i < imgr1->npts; i++) {
			if (imgr1->lon[i] > 180.0)
				imgr1->lon[i] -= 360.0;
			if (imgr1->lat[i] == MISSING_FLOAT)
				imgr1->badmask[i] = YES;
		}

		/*----------------------------------------------------------------------------
		 Read in the solar zenith angle and interpolate to 1000m.
		 ----------------------------------------------------------------------------*/

		if (verbose == YES)
			fprintf(stdout,"%sReading in Solar Zenith Angle.\n", EXE_PROMPT);
		modis_level1b_read("OPEN_1000", imgr1->satid, id, 8, 0.0, 0,
				(void **) &f32junk, imgr1->badmask, &nrow5, &ncol5);
		if (verbose == YES)
			fprintf(stdout,"%sInterpolating Solar Zenith Angle.\n", EXE_PROMPT);
		imgr1->solzen = modis_geo_interp_1000(ncol5, nrow5, f32junk);
		free(f32junk);

		/*----------------------------------------------------------------------------
		 Read in the satellite zenith angle and interpolate to 1000m.
		 ----------------------------------------------------------------------------*/

		if (verbose == YES)
			fprintf(stdout,"%sReading in Satellite Zenith Angle.\n", EXE_PROMPT);
		modis_level1b_read("OPEN_1000", imgr1->satid, id, 9, 0.0, 0,
				(void **) &f32junk, imgr1->badmask, &nrow5, &ncol5);
		if (verbose == YES)
			fprintf(
					stdout,"%sInterpolating Satellite Zenith Angle Latitude.\n",
					EXE_PROMPT);
		imgr1->satzen = modis_geo_interp_1000(ncol5, nrow5, f32junk);
		free(f32junk);

		/*----------------------------------------------------------------------------
		 Read in the relative azimuth and interpolate to 1000m.
		 ----------------------------------------------------------------------------*/

		if (verbose == YES)
			fprintf(stdout,"%sReading in Relative Azimuth.\n", EXE_PROMPT);
		modis_level1b_read("OPEN_1000", imgr1->satid, id, 10, 0.0, 0,
				(void **) &f32junk, imgr1->badmask, &nrow5, &ncol5);
		if (verbose == YES)
			fprintf(stdout,"%sInterpolating Relative Azimuth.\n", EXE_PROMPT);
		imgr1->relaz = modis_geo_interp_1000(ncol5, nrow5, f32junk);
		free(f32junk);

	}

	/*----------------------------------------------------------------------------
	 If geolocation file is available, then read in 1 km geo.
	 ----------------------------------------------------------------------------*/

	else {
		/*----------------------------------------------------------------------------
		 Open the hdf file.
		 ----------------------------------------------------------------------------*/

		fprintf(stdout,"Geo name %s\n", imgr1->geo_filename1);

		geo_id = SDstart(imgr1->geo_filename1, DFACC_READ);
		if (geo_id == FAIL) {
			fprintf(stderr,"%sInvalid HDF file, %s\n",
					EXE_PROMPT,imgr1->geo_filename1);
			exit(EXIT_FAILURE);
		}
		if (verbose == YES)
			fprintf(stdout,"%sReading in Latitude.\n", EXE_PROMPT);
		modis_level03_read("OPEN_1000", geo_id, 1, 0, (void **) &imgr1->lat,
				imgr1->badmask, &imgr1->nrow, &imgr1->ncol);
		if (verbose == YES)
			fprintf(stdout,"%sReading in Longitude.\n", EXE_PROMPT);
		modis_level03_read("OPEN_1000", geo_id, 2, 0, (void **) &imgr1->lon,
				imgr1->badmask, &imgr1->nrow, &imgr1->ncol);
		if (verbose == YES)
			fprintf(stdout,"%sReading in Solar Zenith Angle.\n", EXE_PROMPT);
		modis_level03_read("OPEN_1000", geo_id, 3, 0, (void **) &imgr1->solzen,
				imgr1->badmask, &imgr1->nrow, &imgr1->ncol);
		if (verbose == YES)
			fprintf(stdout,"%sReading in Sensor Zenith Angle.\n", EXE_PROMPT);
		modis_level03_read("OPEN_1000", geo_id, 4, 0, (void **) &imgr1->satzen,
				imgr1->badmask, &imgr1->nrow, &imgr1->ncol);
		if (verbose == YES)
			fprintf(stdout,"%sReading in Relative Azimuth.\n", EXE_PROMPT);
		modis_level03_read("OPEN_1000", geo_id, 5, 0, (void **) &imgr1->relaz,
				imgr1->badmask, &imgr1->nrow, &imgr1->ncol);
		if (verbose == YES)
			fprintf(stdout,"%sReading in land/sea tag.\n", EXE_PROMPT);
		modis_level03_read("OPEN_1000", geo_id, 6, 0,
				(void **) &imgr1->landmask, imgr1->badmask, &imgr1->nrow,
				&imgr1->ncol);
		if (verbose == YES)
			fprintf(stdout,"%sReading in Solar Azimuth Angle.\n", EXE_PROMPT);
		modis_level03_read("OPEN_1000", geo_id, 7, 0, (void **) &imgr1->solaz,
				imgr1->badmask, &imgr1->nrow, &imgr1->ncol);
		if (verbose == YES)
			fprintf(stdout,"%sReading in Sensor Azimuth Angle.\n", EXE_PROMPT);
		modis_level03_read("OPEN_1000", geo_id, 8, 0, (void **) &imgr1->sataz,
				imgr1->badmask, &imgr1->nrow, &imgr1->ncol);

		/*----------------------------------------------------------------------------
		 Close the HDF file.
		 ----------------------------------------------------------------------------*/

		status = SDend(geo_id);
		if (status == FAIL) {
			fprintf(stderr,"%s%s-cannot close %s\n", EXE_PROMPT,rout,
					imgr1->geo_filename1);
			exit(EXIT_FAILURE);
		}

		array_minmax_float(imgr1->npts, imgr1->lat, &minlat, &maxlat);
		array_minmax_float(imgr1->npts, imgr1->lon, &minlon, &maxlon);

//      RAF added code.

//              printf("RAF lons %f %f \n", minlon, maxlon);

		if(minlon < -179.0 && maxlon > 179.0 && maxlat < 75.0 && minlat > -75.0) {
			float splon=1000.0, lnlon=-1000.0;
			for (i = 0; i < imgr1->npts; i++) {
			   if(imgr1->lon[i] > 0.0) {
			      if(imgr1->lon[i] < splon) splon = imgr1->lon[i];
			   }
			   if(imgr1->lon[i] < 0.0) {
			      if(imgr1->lon[i] > lnlon) lnlon = imgr1->lon[i];
			   }
			}
			maxlon = splon;
			minlon = lnlon;
//                      printf("new lons: %f %f\n", minlon, maxlon);
		}


		scale_lon_flg = NO;
		if (minlat < 80.0 && maxlat > -80.0 && minlon < -60.0 && maxlon > 60.0) {
			scale_lon_flg = YES;
			for (i = 0; i < imgr1->npts; i++) {
				if (imgr1->lon[i] < 0.0)
					imgr1->lon[i] += 360.0;
			}
		}
		if (scale_lon_flg == YES) {
			for (i = 0; i < imgr1->npts; i++) {
				if (imgr1->lon[i] > 180.0)
					imgr1->lon[i] -= 360.0;
			}
		}
	}

	/*----------------------------------------------------------------------------
	 Set min/max lat and lon limit array.
	 ----------------------------------------------------------------------------*/

	imgr1->bounds[0] = minlat, imgr1->bounds[1] = minlon;
	imgr1->bounds[2] = maxlat, imgr1->bounds[3] = maxlon;
	if (minlon > 180.0)
		minlon -= 360.0;
	if (maxlon > 180.0)
		maxlon -= 360.0;

	if (verbose == YES)
		fprintf(
				stdout,"\n%sMODIS 1 km dimensions: npixels=%d, nlines=%d, ntotal=%ld\n",
				EXE_PROMPT, imgr1->ncol, imgr1->nrow, imgr1->npts);

	/*----------------------------------------------------------------------------
	 Allocate memory for various parameters that are to be calculated.
	 ----------------------------------------------------------------------------*/

	if ((imgr1->cos_solzen = (float *) malloc(imgr1->npts * sizeof(float)))
			== NULL)error_allo(rout,"imgr1->cos_solzen");
	if ((imgr1->cos_satzen = (float *) malloc(imgr1->npts * sizeof(float)))
			== NULL)error_allo(rout,"imgr1->cos_satzen");

	if ((imgr1->scatzen = (float *) malloc(imgr1->npts * sizeof(float)))
			== NULL)error_allo(rout,"imgr1->scatzen");
	if ((imgr1->glintzen = (float *) malloc(imgr1->npts * sizeof(float)))
			== NULL)error_allo(rout,"imgr1->glintzen");

	/*----------------------------------------------------------------------------
	 Calculate trigometric parameters.
	 ----------------------------------------------------------------------------*/

//      if (verbose == YES)
//              fprintf(stdout,"%sCalculating scattering and glint angles\n",
//		EXE_PROMPT);

	for (i = 0; i < imgr1->npts; i++) {
		imgr1->cos_solzen[i] = cos(DTOR*imgr1->solzen[i]);
		imgr1->cos_satzen[i] = cos(DTOR*imgr1->satzen[i]);

		sin_solzen = sin(DTOR*imgr1->solzen[i]);
		sin_satzen = sin(DTOR*imgr1->satzen[i]);
		cos_relaz = cos(DTOR*imgr1->relaz[i]);

		imgr1->scatzen[i] = (acos(imgr1->cos_solzen[i] * imgr1->cos_satzen[i]
				+ sin_solzen * sin_satzen * cos_relaz)) / DTOR;
		imgr1->glintzen[i] = (acos(-1.0 * imgr1->cos_solzen[i]
				* imgr1->cos_satzen[i] - sin_solzen * sin_satzen * cos_relaz))
				/ DTOR;

	}

	/*----------------------------------------------------------------------------
	 Determine ascending/descending flag.
	 ----------------------------------------------------------------------------*/

	des_node = NO;
	asc_node = NO;
	if ((imgr1->asc_des_flg = (unsigned char *) malloc(imgr1->nrow
			* sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1->asc_des_flg");
	i_nadir = (int) (imgr1->ncol / 2.0);
	for (i = 0; i < imgr1->nrow - 1; i++) {
		imgr1->asc_des_flg[i] = MISSING_NODE;
		if (imgr1->badmask[(i * imgr1->ncol) + i_nadir] == NO
				&& imgr1->badmask[((i + 1) * imgr1->ncol) + i_nadir] == NO) {
			if (imgr1->lat[(i * imgr1->ncol) + i_nadir] >= imgr1->lat[((i + 1)
					* imgr1->ncol) + i_nadir]) {
				imgr1->asc_des_flg[i] = DES_NODE;
				des_node = YES;
			} else {
				imgr1->asc_des_flg[i] = ASC_NODE;
				asc_node = YES;
			}
		}
	}
	imgr1->asc_des_flg[imgr1->ncol - 1] = imgr1->asc_des_flg[imgr1->ncol - 2];
	if (des_node == YES && asc_node == NO)
		strcpy(imgr1->asc_des_attr, "Descending");
	else if (des_node == NO && asc_node == YES)
		strcpy(imgr1->asc_des_attr, "Ascending");
	else if (des_node == YES && asc_node == YES)
		strcpy(imgr1->asc_des_attr, "Both");
	else {
		fprintf(stderr,"%sCannot determine satellite node - aborting\n",
				EXE_PROMPT);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Determine how many spectral channels are to be read in.
	 ----------------------------------------------------------------------------*/

	if (chflg[19] && !chflg[30])
		chflg[30] = 1;
	nchan = 0;
	for (i = 0; i < NMODIS_CHAN; i++) {
		if (chflg[i])
			nchan++;
	}
	if (verbose == YES)
		fprintf(stdout,"%s%d channels are to be read.\n\n", EXE_PROMPT,nchan);

	/*----------------------------------------------------------------------------
	 Read in data from needed channels.
	 ----------------------------------------------------------------------------*/

	if (nchan > 0) {
		if (chflg[0]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 1 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 1.0, 0,
					(void **) &imgr1->ref1, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref1[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[1]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 2 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 2.0, 0,
					(void **) &imgr1->ref2, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref2[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[2]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 3 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 3.0, 0,
					(void **) &imgr1->ref3, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref3[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[3]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 4 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 4.0, 0,
					(void **) &imgr1->ref4, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref4[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[4]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 5 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 5.0, 0,
					(void **) &imgr1->ref5, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref5[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[5]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 6 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 6.0, 0,
					(void **) &imgr1->ref6, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref6[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[6]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 7 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 7.0, 0,
					(void **) &imgr1->ref7, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref7[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[7]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 8 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 8.0, 0,
					(void **) &imgr1->ref8, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref8[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[8]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 9 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 9.0, 0,
					(void **) &imgr1->ref9, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref9[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[9]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 10 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 10.0, 0,
					(void **) &imgr1->ref10, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref10[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[10]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 11 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 11.0, 0,
					(void **) &imgr1->ref11, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref11[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[11]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 12 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 12.0, 0,
					(void **) &imgr1->ref12, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref12[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[12]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 13lo Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 13.0, 0,
					(void **) &imgr1->ref13l, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref13l[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[13]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 14lo Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 14.0, 0,
					(void **) &imgr1->ref14l, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref14l[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[14]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 15 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 15.0, 0,
					(void **) &imgr1->ref15, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref15[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[15]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 16 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 16.0, 0,
					(void **) &imgr1->ref16, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref16[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[16]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 17 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 17.0, 0,
					(void **) &imgr1->ref17, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref17[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[17]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 18 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 18.0, 0,
					(void **) &imgr1->ref18, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref18[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[18]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 19 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 19.0, 0,
					(void **) &imgr1->ref19, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref19[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[19]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 20 Radiance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 20.0, 0,
					(void **) &imgr1->rad20, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			if (verbose == YES)
				fprintf(
						stdout,"%sReading in Channel 20 Brightness Temperature.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 20.0, 0,
					(void **) &imgr1->bt20, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[20]) {
			if (verbose == YES)
				fprintf(
						stdout,"%sReading in Channel 21 Brightness Temperature.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 21.0, 0,
					(void **) &imgr1->bt21, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[21]) {
			if (verbose == YES)
				fprintf(
						stdout,"%sReading in Channel 22 Brightness Temperature.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 22.0, 0,
					(void **) &imgr1->bt22, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[22]) {
			if (verbose == YES)
				fprintf(
						stdout,"%sReading in Channel 23 Brightness Temperature.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 23.0, 0,
					(void **) &imgr1->bt23, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[23]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 24 Radiance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 24.0, 0,
					(void **) &imgr1->rad24, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[24]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 25 Radiance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 25.0, 0,
					(void **) &imgr1->rad25, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[25]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 26 Reflectance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 4, 26.0, 0,
					(void **) &imgr1->ref26, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			for (i = 0; i < imgr1->npts; i++)
				imgr1->ref26[i] /= imgr1->cos_solzen[i];
		}
		if (chflg[26]) {
			if (verbose == YES)
				fprintf(
						stdout,"%sReading in Channel 27 Brightness Temperature.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 27.0, 0,
					(void **) &imgr1->bt27, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[27]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 28 Radiance.\n",EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 28.0, 0,
					(void **) &imgr1->rad28, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
				fprintf(stdout,"%sReading in Channel 28 Brightness Temperature.\n",EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 28.0, 0,
					(void **) &imgr1->bt28, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[28]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 29 Radiance.\n",EXE_PROMPT);	// GPC
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 29.0, 0,
					(void **) &imgr1->rad29, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);	// GPC
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 29 Brightness Temperature.\n",EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 29.0, 0,
					(void **) &imgr1->bt29, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[29]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 30 Brightness Temperature.\n",EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 30.0, 0,
					(void **) &imgr1->bt30, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[30]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 31 Radiance.\n",EXE_PROMPT);	// GPC
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 31.0, 0,
					(void **) &imgr1->rad31, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);	// GPC
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 31 Brightness Temperature.\n",EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 31.0, 0,
					(void **) &imgr1->bt31, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[31]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 32 Radiance.\n",EXE_PROMPT);	// GPC
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 32.0, 0,
					(void **) &imgr1->rad32, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);	// GPC
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 32 Brightness Temperature.\n",EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 32.0, 0,
					(void **) &imgr1->bt32, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[32]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 33 Radiance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 33.0, 0,
					(void **) &imgr1->rad33, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			if (verbose == YES) fprintf(stdout,"%sReading in Channel 33 Brightness Temperature.\n",EXE_PROMPT);
			 modis_level1b_read ("OPEN_1000", imgr1->satid, id, 5, 33.0, 0, 
					(void **) &imgr1->bt33, imgr1->badmask, &imgr1->nrow, &imgr1->ncol);
		}
		if (chflg[33]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 34 Radiance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 34.0, 0,
					(void **) &imgr1->rad34, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			if (verbose == YES)
				fprintf(
						stdout,"%sReading in Channel 34 Brightness Temperature.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 34.0, 0,
					(void **) &imgr1->bt34, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[34]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 35 Radiance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 35.0, 0,
					(void **) &imgr1->rad35, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			if (verbose == YES)
				fprintf(
						stdout,"%sReading in Channel 35 Brightness Temperature.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 35.0, 0,
					(void **) &imgr1->bt35, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
		if (chflg[35]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 36 Radiance.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 3, 36.0, 0,
					(void **) &imgr1->rad36, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
			if (verbose == YES)
				fprintf(
						stdout,"%sReading in Channel 36 Brightness Temperature.\n",
						EXE_PROMPT);
			modis_level1b_read("OPEN_1000", imgr1->satid, id, 5, 36.0, 0,
					(void **) &imgr1->bt36, imgr1->badmask, &imgr1->nrow,
					&imgr1->ncol);
		}
	}

	/*----------------------------------------------------------------------------
	 Close the HDF file.
	 ----------------------------------------------------------------------------*/

	status = SDend(id);
	if (status == FAIL) {
		fprintf(stderr,"%s%s-cannot close %s\n", EXE_PROMPT,rout, full_filename);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 If channel 20 is needed, calculate the pseudo-reflectance and emissivity.
	 ----------------------------------------------------------------------------*/

	if (chflg[19]) {
//              if (verbose == YES)
//	                fprintf(
//			stdout,"\n%sCalculating channel 20 pseudo reflectance and emissivity\n",
//			EXE_PROMPT);
		if ((imgr1->ref20 = (float *) malloc(imgr1->npts * sizeof(float)))
				== NULL)error_allo(rout,"imgr1->ref20");
		if ((imgr1->ems20 = (float *) malloc(imgr1->npts * sizeof(float)))
				== NULL)error_allo(rout,"imgr1->ems20");
		for (i = 0; i < imgr1->npts; i++) {
			imgr1->ref20[i] = MISSING_FLOAT;
			imgr1->ems20[i] = MISSING_FLOAT;
			if (imgr1->badmask[i] == NO) {
				if (imgr1->satid == TERRA) {
					ch20_ems = modis_planck_terra(20, imgr1->bt31[i]);
					SOLAR_20 = SOLAR_20_TERRA;
				} else if (imgr1->satid) {
					ch20_ems = modis_planck_aqua(20, imgr1->bt31[i]);
					SOLAR_20 = SOLAR_20_AQUA;
				} else {
					fprintf(stderr,"%s%s - Unknown platform name - aborting\n",
							EXE_PROMPT,rout);
					exit(EXIT_FAILURE);
				}
				imgr1->ref20[i] = PI * (imgr1->rad20[i] - ch20_ems)
						/ (SOLAR_20 * (imgr1->cos_solzen[i] / pow(
								imgr1->earth2sun_dist, 2)) - PI * ch20_ems);
				imgr1->ems20[i] = (imgr1->rad20[i]) / (ch20_ems);
			}
		}
	}

	if (verbose == YES)
		fprintf(
				stdout,"\n%sMODIS: {minlat=%f, maxlat=%f, minlon=%f, maxlon=%f}\n",
				EXE_PROMPT, imgr1->bounds[0], imgr1->bounds[2],
				imgr1->bounds[1], imgr1->bounds[3]);
}

void read_modisqkm(int8 *chflg_qkm, imagerL1 *imgr1, int8 verbose) {


	char *rout = { "read_modisqkm" },
			full_filename[MAX_STR_LEN],
			att_name[MAX_STR_LEN];

	int i, nchan;
	uint16 *ui16junk;
	int32 id, status, att_index, data_type, count;
	float32 *earth2sun_temp;
	void *attr;
	uint8 *detector_qc;


//      printf("read_modisqkm: qkm_filename is: %s %d\n", imgr1->qkm_filename, verbose);
//      (void)fflush(stdout);

	if (verbose == YES)
		fprintf(stdout,"%ssatname=%s, year=%d, jday=%d, time=%f\n",
				EXE_PROMPT,imgr1->satname, imgr1->year, imgr1->jday,
				imgr1->time);

	/*----------------------------------------------------------------------------
	 Open the hdf file.
	 ----------------------------------------------------------------------------*/

	sprintf(full_filename, "%s/%s", imgr1->qkm_dir_name, imgr1->qkm_filename);
	id = SDstart(full_filename, DFACC_READ);
	if (id == FAIL) {
		fprintf(stderr,"%sInvalid HDF file, %s\n", EXE_PROMPT,full_filename);
		exit(EXIT_FAILURE);
	}

	/*----------------------------------------------------------------------------
	 Read in the Earth-Sun Distance global attribute.
	 ----------------------------------------------------------------------------*/

	strcpy(att_name, "Earth-Sun Distance");
	att_index = SDfindattr(id, att_name);
	status = SDattrinfo(id, att_index, att_name, &data_type, &count);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDattrinfo - Earth-Sun Distance.");

	if ((attr = (void *) malloc(count * DFKNTsize(data_type))) == NULL)error_allo(rout,"attr-Earth-Sun Distance");
	status = SDreadattr(id, att_index, attr);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDreadattr - Earth-Sun Distance.");

	earth2sun_temp = (float32 *) attr;
	imgr1->earth2sun_dist = earth2sun_temp[0];
	free(earth2sun_temp);

	/*----------------------------------------------------------------------------
	 Read in detector quality flags.
	 ----------------------------------------------------------------------------*/

	strcpy(att_name, "Detector Quality Flag");
	att_index = SDfindattr(id, att_name);
	status = SDattrinfo(id, att_index, att_name, &data_type, &count);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDattrinfo - Detector Quality Flag.");

	if ((attr = (void *) malloc(count * DFKNTsize(data_type))) == NULL)error_allo(rout,"attr-Detector Quality Flag");
	status = SDreadattr(id, att_index, attr);
	if (status == FAIL)error_exit(rout,"There was an error when calling SDreadattr - Detector Quality Flag.");

	detector_qc = (uint8 *) attr;

	/*----------------------------------------------------------------------------
	 Unpack detector quality flags.
	 ----------------------------------------------------------------------------*/

	imgr1->noise_det = unpack_detector_flags(1, count, detector_qc);
	imgr1->dead_det = unpack_detector_flags(2, count, detector_qc);
	imgr1->gain_det = unpack_detector_flags(3, count, detector_qc);
	imgr1->range_det = unpack_detector_flags(4, count, detector_qc);
	imgr1->source_det = unpack_detector_flags(5, count, detector_qc);
	imgr1->residuals_det = unpack_detector_flags(6, count, detector_qc);
	imgr1->crosstalk_det = unpack_detector_flags(7, count, detector_qc);
	free(detector_qc);

	/*----------------------------------------------------------------------------
	 Read in the 250 m dimensions.
	 ----------------------------------------------------------------------------*/

//      printf("read_modisqkm: get the 250 m dimensions\n");
	modis_level1b_read("OPEN_250", imgr1->satid, id, 1, 1.0, 0,
			(void **) &ui16junk, imgr1->qkm_badmask, &imgr1->qkm_nrow, &imgr1->qkm_ncol);
	free(ui16junk);

	imgr1->qkm_npts = (long) imgr1->qkm_nrow * (long) imgr1->qkm_ncol;

	if ((imgr1->qkm_badmask = (unsigned char *) malloc(imgr1->qkm_npts
			* sizeof(unsigned char))) == NULL)error_allo(rout,"imgr1->qkm_badmask");
	for (i = 0; i < imgr1->qkm_npts; i++)
		imgr1->qkm_badmask[i] = NO;

  
	if (verbose == YES)
		fprintf(
				stdout,"\n%sMODIS 250-m dimensions: npixels=%d, nlines=%d, ntotal=%ld\n",
				EXE_PROMPT, imgr1->qkm_ncol, imgr1->qkm_nrow, imgr1->qkm_npts);

	/*----------------------------------------------------------------------------
	 Determine how many spectral channels are to be read in.
	 ----------------------------------------------------------------------------*/

	nchan = 0;
	for (i = 0; i < N250M_CHAN; i++) {
		if (chflg_qkm[i])
			nchan++;
	}
	if (verbose == YES)
		fprintf(stdout,"%s%d channels are to be read.\n\n", EXE_PROMPT,nchan);

	/*----------------------------------------------------------------------------
	 Read in data from needed channels.
	 ----------------------------------------------------------------------------*/

	if (nchan > 0) {
		if (chflg_qkm[0]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 1 Reflectance.\n",
						EXE_PROMPT);

                        if ((imgr1->qkm_ref1 = (float32 *) malloc(imgr1->qkm_npts*sizeof(float32))) == NULL)
                            error_allo(rout,"imgr1->qkm_ref1");

			modis_level1b_read("OPEN_250", imgr1->satid, id, 4, 1.0, 0,
  			                (void **) &imgr1->qkm_ref1, imgr1->qkm_badmask, &imgr1->qkm_nrow,
					&imgr1->qkm_ncol);

		}

		if (chflg_qkm[1]) {
			if (verbose == YES)
				fprintf(stdout,"%sReading in Channel 2 Reflectance.\n",
						EXE_PROMPT);

                        if ((imgr1->qkm_ref2 = (float32 *) malloc(imgr1->qkm_npts*sizeof(float32))) == NULL)
                            error_allo(rout,"imgr1->qkm_ref2");

			modis_level1b_read("OPEN_250", imgr1->satid, id, 4, 2.0, 0,
  			                (void **) &imgr1->qkm_ref2, imgr1->qkm_badmask, &imgr1->qkm_nrow,
					&imgr1->qkm_ncol);

		}
	}

	/*----------------------------------------------------------------------------
	 Close the HDF file.
	 ----------------------------------------------------------------------------*/

	status = SDend(id);
	if (status == FAIL) {
		fprintf(stderr,"%s%s-cannot close %s\n", EXE_PROMPT,rout, full_filename);
		exit(EXIT_FAILURE);
	}

//      if (verbose == YES)
//		fprintf(
//				stdout,"\n%sMODIS: {minlat=%f, maxlat=%f, minlon=%f, maxlon=%f}\n",
//				EXE_PROMPT, imgr1->bounds[0], imgr1->bounds[2],
//				imgr1->bounds[1], imgr1->bounds[3]);
   
}

void deallocate_modis(int8 *chflg, imagerL1 *imgr1) {
	int i, nchan;

	// Free quantities that are always allocated memory

	free(imgr1->noise_det);
	free(imgr1->dead_det);
	free(imgr1->gain_det);
	free(imgr1->range_det);
	free(imgr1->source_det);
	free(imgr1->residuals_det);
	free(imgr1->crosstalk_det);


	free(imgr1->lat);
	free(imgr1->lon);
	free(imgr1->solzen);
	free(imgr1->satzen);
	free(imgr1->relaz);
	free(imgr1->cos_solzen);
	free(imgr1->cos_satzen);
	free(imgr1->scatzen);
	free(imgr1->glintzen);
	free(imgr1->sfc_type);
	free(imgr1->zsfc);
	free(imgr1->sst);
	free(imgr1->asc_des_flg);
	free(imgr1->badmask);

	// Determine which spectral channels were allocated memory,
	// and free the associated arrays

	nchan = 0;
	for (i = 0; i < NMODIS_CHAN; i++) {
		if (chflg[i])
			nchan++;
	}

	if (nchan > 0) {
		if (chflg[0]) {
			free(imgr1->ref1);
		}
		if (chflg[1]) {
			free(imgr1->ref2);
		}
		if (chflg[2]) {
			free(imgr1->ref3);
		}
		if (chflg[3]) {
			free(imgr1->ref4);
		}
		if (chflg[4]) {
			free(imgr1->ref5);
		}
		if (chflg[5]) {
			free(imgr1->ref6);
		}
		if (chflg[6]) {
			free(imgr1->ref7);
		}
		if (chflg[7]) {
			free(imgr1->ref8);
		}
		if (chflg[8]) {
			free(imgr1->ref9);
		}
		if (chflg[9]) {
			free(imgr1->ref10);
		}
		if (chflg[10]) {
			free(imgr1->ref11);
		}
		if (chflg[11]) {
			free(imgr1->ref12);
		}
		if (chflg[12]) {
			free(imgr1->ref13l);
		}
		if (chflg[13]) {
			free(imgr1->ref14l);
		}
		if (chflg[14]) {
			free(imgr1->ref15);
		}
		if (chflg[15]) {
			free(imgr1->ref16);
		}
		if (chflg[16]) {
			free(imgr1->ref17);
		}
		if (chflg[17]) {
			free(imgr1->ref18);
		}
		if (chflg[18]) {
			free(imgr1->ref19);
		}
		if (chflg[19]) {
			free(imgr1->ref20);
			free(imgr1->rad20);
			free(imgr1->ems20);
			free(imgr1->bt20);
		}
		if (chflg[20]) {
			free(imgr1->bt21);
		}
		if (chflg[21]) {
			free(imgr1->bt22);
		}
		if (chflg[22]) {
			free(imgr1->bt23);
		}
		if (chflg[23]) {
			free(imgr1->rad24);
		}
		if (chflg[24]) {
			free(imgr1->rad25);
		}
		if (chflg[25]) {
			free(imgr1->ref26);
		}
		if (chflg[26]) {
			free(imgr1->bt27);
		}
		if (chflg[27]) {
			free(imgr1->bt28);
		}
		if (chflg[28]) {
			free(imgr1->bt29);
		}
		if (chflg[29]) {
			free(imgr1->bt30);
		}
		if (chflg[20]) {
			free(imgr1->bt31);
		}
		if (chflg[31]) {
			free(imgr1->bt32);
		}
		if (chflg[32]) {
			free(imgr1->rad33);
		}
		if (chflg[33]) {
			free(imgr1->rad34);
		}
		if (chflg[34]) {
			free(imgr1->rad35);
		}
		if (chflg[35]) {
			free(imgr1->rad36);
		}
		
		/*if (chflg[32]) {free(imgr1->bt33);}
		 if (chflg[33]) {free(imgr1->bt34);}
		 if (chflg[34]) {free(imgr1->bt35);}
		 if (chflg[35]) {free(imgr1->bt36);}*/
	}
}

void get_geo_filename(int8 verbose, int8 *geo_interp, imagerL1 *imgr1) {

        char sys_command[MAX_STR_LEN];

	sprintf(sys_command, "echo nofile > temp_flist");
	system(sys_command);

//      if (verbose == YES)
//             fprintf(
//		      stdout,"%sGenerating MODIS geolocation file list using UNIX commands\n",
//		      EXE_PROMPT);
//      if ((imgr1->geo_filename1 != NULL) && (strlen(imgr1->geo_filename1) != 0)) {
//
//              printf("Command line supplied name: %s\n", imgr1->geo_filename1);
//      }

	/*----------------------------------------------------------------------------
	 Check to make sure a file was found.
	 ----------------------------------------------------------------------------*/

	if (strlen(imgr1->geo_filename1) == 0) {
		fprintf(
				stderr,"\n%sCannot find matching MODIS geolocation file - aborting\n",
				EXE_PROMPT);
		exit(EXIT_FAILURE);
	}

//      if (verbose)
//              fprintf(stdout,"%sGeo file = %s\n", EXE_PROMPT,imgr1->geo_filename1);
}
