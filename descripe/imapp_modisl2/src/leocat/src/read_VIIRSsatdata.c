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

void get_VIIRS_geo_filename(int8, int8 *, imagerL1 *);
void viirs_geoloc_read(char *, hid_t, void **, int *, int *, int *);
int viirs_level1_read(char *, hid_t, void **, int *, int *, int *);

void get_VIIRS_geo_filename(int8 verbose, int8 *geo_interp, imagerL1 *imgr1_ptr) {
	//FILE *fptr;
	char sys_command[MAX_STR_LEN];
	//char sys_command1[MAX_STR_LEN], sys_command[MAX_STR_LEN];
	//int yy = 0;

	sprintf(sys_command, "echo nofile > temp_flist");
	system(sys_command);

	if (verbose == YES)
		fprintf(
				stdout,"%sGenerating MODIS geolocation file list using UNIX commands\n",
				EXE_PROMPT);
	if ((imgr1_ptr->geo_filename1 != NULL) && (strlen(imgr1_ptr->geo_filename1) != 0)) {

		printf("Command line supplied name: %s\n", imgr1_ptr->geo_filename1);
	}

	/*----------------------------------------------------------------------------
	 Check to make sure a file was found.
	 ----------------------------------------------------------------------------*/

	if (strlen(imgr1_ptr->geo_filename1) == 0) {
		fprintf(
				stderr,"\n%sCannot find matching VIIRS geolocation file - aborting\n",
				EXE_PROMPT);
		exit(EXIT_FAILURE);
	}

	if (verbose)
		fprintf(stdout,"%sGeo file = %s\n", EXE_PROMPT,imgr1_ptr->geo_filename1);
}

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *	                         VIIRS Level 1 Data                              *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************/

// TODO : These FILL values are from the most recent CodeBook (CDFCB-X_Volume_3) 
//        for the VIIRS SDR. Remove these definitions when we have the 1.5.0.48
//        codebase installed

#define NA_UINT16_FILL          65535
#define MISS_UINT16_FILL        65534
#define ONBOARD_PT_UINT16_FILL  65533
#define ONGROUND_PT_UINT16_FILL 65532
#define ERR_UINT16_FILL         65531
#define VDNE_UINT16_FILL        65529
#define SOUB_UINT16_FILL        65528

#define NA_FLOAT32_FILL          -999.9
#define MISS_FLOAT32_FILL        -999.8
#define ONBOARD_PT_FLOAT32_FILL  -999.7
#define ONGROUND_PT_FLOAT32_FILL -999.6
#define ERR_FLOAT32_FILL         -999.5
#define VDNE_FLOAT32_FILL        -999.3

#define FLOAT64_FILL_TEST       -999.0e0
#define INT64_FILL_TEST         -990
#define UINT64_FILL_TEST         18446744073709551607
#define FLOAT32_FILL_TEST       -999.0
#define INT32_FILL_TEST         -990
#define UINT32_FILL_TEST         424967287
#define INT16_FILL_TEST         -990
#define UINT16_FILL_TEST         65527
#define INT8_FILL_TEST           119
#define UINT8_FILL_TEST          247

void read_viirs(int8 *chflg, imagerL1 *imgr1_ptr, int8 verbose) {

	char *rout = { "read_viirs" }, full_filename[MAX_STR_LEN];

	char dsetName[100];
	int i, j,row,col,scan,detector;
	long idx;
	int ichar, nchan;
	int nRank;
	int nrows,ncols,nscans,ndets;
	long npts;

	void *bufferPtr_void;

	hid_t   fileID, dspaceID, dsetID;          /* Handles */
	herr_t  h5status;
	hsize_t dimSize[2];

	/*****************************************************************************
     *	                         VIIRS Geolocation                               *
     *****************************************************************************/

	/*----------------------------------------------------------------------------
	 Use the filename to determine the satname and if direct broadcast.
	 ----------------------------------------------------------------------------*/

	printf("read_viirs: %s\n", imgr1_ptr->l1b_filename1);

	printf("Geo file name : %s\n", imgr1_ptr->geo_filename1);
	sprintf(full_filename, "%s", imgr1_ptr->geo_filename1);
	printf("Full Geo name %s\n", full_filename);
	
	printf("Level 1 dir name  : %s\n", imgr1_ptr->l1b_dir_name);
	printf("Level 1 file name : %s\n", imgr1_ptr->l1b_filename1);
	sprintf(full_filename, "%s/%s", imgr1_ptr->l1b_dir_name, imgr1_ptr->l1b_filename1);
	printf("Full Level 1 name %s\n", full_filename);

	if (verbose == YES)
		fprintf(stdout,"\n%ssatname=%s, year=%d, jday=%d, time=%f\n",
				EXE_PROMPT,imgr1_ptr->satname, imgr1_ptr->year, imgr1_ptr->jday,
				imgr1_ptr->time);

	/*----------------------------------------------------------------------------
	 Open the geolocation file
	 ----------------------------------------------------------------------------*/

	char ModGeoGroup[] =  "All_Data/VIIRS-MOD-GEO_All/";

	printf("Geo file name : %s\n", imgr1_ptr->geo_filename1);
	sprintf(full_filename, "%s", imgr1_ptr->geo_filename1);

	// Open the file as Read Only, and default File Access property list
	fileID = H5Fopen(full_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (fileID == FAIL) {
		fprintf(stderr,"%sInvalid HDF file, %s\n", EXE_PROMPT,full_filename);
		exit(EXIT_FAILURE);
	}

	// Set the number of VIIRS MOD resolution detectors (TODO : Set from viirs_sdr.h)
	imgr1_ptr->nDetectors = 16;

    /******************************
	 *  Read in the NumberOfScans *
     ******************************/

	sprintf(dsetName, "%s%s", ModGeoGroup,"NumberOfScans");

	// Get the dataset ID
	dsetID = H5Dopen2 (fileID, dsetName, H5P_DEFAULT);

	// Read the ActualScans 
	h5status = H5Dread (dsetID,              // dataset ID
					  H5T_NATIVE_INT,    // ident of memory datatype
					  H5S_ALL,             // ident of memory dataspace
					  H5S_ALL,             // ident of dataset's dataspace in file
					  H5P_DEFAULT,         // File access properties
					  &imgr1_ptr->nScans);        // Pointer to buffer receiving data

	printf("\n%s = %d\n",dsetName,imgr1_ptr->nScans);

    /*****************************************
	 *  Read in the ModeScan dataset  *
     *****************************************/

	sprintf(dsetName, "%s%s", ModGeoGroup,"ModeScan");
	dsetID = H5Dopen2 (fileID, dsetName, H5P_DEFAULT);
	dspaceID = H5Dget_space(dsetID);
	nRank = H5Sget_simple_extent_dims(dspaceID, dimSize, NULL);
	int nViirsScans=dimSize[0];		// FIXME : Get this value from viirsRDR.h
	imgr1_ptr->ModeScan = (int *) malloc(dimSize[0] * sizeof(int));

	// Read the ModeScan dataset into the ModeScan array
	h5status = H5Dread (dsetID,               // dataset ID
					  H5T_NATIVE_INT,       // ident of memory datatype
					  H5S_ALL,              // ident of memory dataspace
					  H5S_ALL,              // ident of dataset's dataspace in file
					  H5P_DEFAULT,          // File access properties
					  imgr1_ptr->ModeScan); // Pointer to buffer receiving data

	printf("\nModeScan =\n");
	printf (" [ ");
	for (scan=0; scan<nViirsScans; scan++){
		printf (" %4d", imgr1_ptr->ModeScan[scan]);
	}
	printf (" ]\n");

	// Get first and last good rows of geolocation

	int firstGoodScan,lastGoodScan,firstGoodRow,lastGoodRow;

	printf("Searching for first good row...\n");
	scan = 0;
	while (scan<nViirsScans){
		printf("\tModeScan[%d] = %d\n",scan,imgr1_ptr->ModeScan[scan]);
		if (imgr1_ptr->ModeScan[scan]<=2){
			firstGoodScan = scan;
			break;
		}
		scan++;
	}
	firstGoodRow = firstGoodScan*imgr1_ptr->nDetectors;

	printf("Searching for last good row...\n");
	scan = nViirsScans-1;
	while (scan>-1){
		printf("\tModeScan[%d] = %d\n",scan,imgr1_ptr->ModeScan[scan]);
		if (imgr1_ptr->ModeScan[scan]<=2){
			lastGoodScan = scan;
			break;
		}
		scan--;
	}
	lastGoodRow = imgr1_ptr->nDetectors*(lastGoodScan+1)-1;

	printf("firstGoodScan/Row = %d,%d\n",firstGoodScan,firstGoodRow);
	printf("lastGoodScan/Row = %d,%d\n",lastGoodScan,lastGoodRow);

    /**********************************
	 *  Read in the Latitude dataset  *
     **********************************/

	sprintf(dsetName, "%s%s", ModGeoGroup,"Latitude");
	printf("\nDataset name is : %s\n",dsetName);

	// Read into the imgr1_ptr structure
	viirs_geoloc_read (dsetName, fileID, (void **) &bufferPtr_void, 
		&nRank, &imgr1_ptr->nrow, &imgr1_ptr->ncol);

	dimSize[0] = imgr1_ptr->nrow;
	dimSize[1] = imgr1_ptr->ncol;

	imgr1_ptr->npts = (long) imgr1_ptr->nrow * (long) imgr1_ptr->ncol;

	printf("imgr1_ptr->nScans = %d\n",imgr1_ptr->nScans);
	printf("imgr1_ptr->nrow = %d\n",imgr1_ptr->nrow);
	printf("imgr1_ptr->ncol = %d\n",imgr1_ptr->ncol);
	printf("imgr1_ptr->npts = %ld\n",imgr1_ptr->npts);
	printf("Returned rank = %d\n",nRank);

	npts = imgr1_ptr->npts;
	nrows = imgr1_ptr->nrow;
	ncols = imgr1_ptr->ncol;
	nscans = imgr1_ptr->nScans;
	ndets = imgr1_ptr->nDetectors;

	imgr1_ptr->lat = (float *) bufferPtr_void;
	bufferPtr_void = 0;

    /***********************************
	 *  Read in the Longitude dataset  *
     ***********************************/

	sprintf(dsetName, "%s%s", ModGeoGroup,"Longitude");
	printf("\nDataset name is : %s\n",dsetName);

	// Read into the imgr1_ptr structure
	viirs_geoloc_read (dsetName, fileID, (void **) &bufferPtr_void, &nRank, &nrows, &ncols);

	imgr1_ptr->lon = (float *) bufferPtr_void;
	bufferPtr_void = 0;

    /**************************************
	 *  Read in the Solar Zenith dataset  *
     **************************************/

	sprintf(dsetName, "%s%s", ModGeoGroup,"SolarZenithAngle");
	printf("\nDataset name is : %s\n",dsetName);

	// Read into the imgr1_ptr structure
	viirs_geoloc_read (dsetName, fileID, (void **) &bufferPtr_void, &nRank, &nrows, &ncols);

	imgr1_ptr->solzen = (float *) bufferPtr_void;
	bufferPtr_void = 0;

    /***************************************
	 *  Read in the Solar Azimuth dataset  *
     ***************************************/

	sprintf(dsetName, "%s%s", ModGeoGroup,"SolarAzimuthAngle");
	printf("\nDataset name is : %s\n",dsetName);

	// Read into the imgr1_ptr structure
	viirs_geoloc_read (dsetName, fileID, (void **) &bufferPtr_void, &nRank, &nrows, &ncols);

	imgr1_ptr->solaz = (float *) bufferPtr_void;
	bufferPtr_void = 0;

    /******************************************
	 *  Read in the Satellite Zenith dataset  *
     ******************************************/

	sprintf(dsetName, "%s%s", ModGeoGroup,"SatelliteZenithAngle");
	printf("\nDataset name is : %s\n",dsetName);

	// Read into the imgr1_ptr structure
	viirs_geoloc_read (dsetName, fileID, (void **) &bufferPtr_void, &nRank, &nrows, &ncols);

	imgr1_ptr->satzen = (float *) bufferPtr_void;
	bufferPtr_void = 0;

    /******************************************
	 *  Read in the Satellite Azimuth dataset  *
     ******************************************/

	sprintf(dsetName, "%s%s", ModGeoGroup,"SatelliteAzimuthAngle");
	printf("\nDataset name is : %s\n",dsetName);

	// Read into the imgr1_ptr structure
	viirs_geoloc_read (dsetName, fileID, (void **) &bufferPtr_void, &nRank, &nrows, &ncols);

	imgr1_ptr->sataz = (float *) bufferPtr_void;
	bufferPtr_void = 0;

    /********************************
	 *  Read in the Terrain Height  *
     ********************************/

	// FIXME : Terrain Height from VIIRS geolocation is buggered
	/*sprintf(dsetName, "%s%s", ModGeoGroup,"Height");*/
	/*printf("\nDataset name is : %s\n",dsetName);*/

	// Read into the imgr1_ptr structure
	/*viirs_geoloc_read (dsetName, fileID, (void **) &bufferPtr_void, &nRank, &nrows, &ncols);*/

	/*imgr1_ptr->zsfc = (float *) bufferPtr_void;*/
	/*bufferPtr_void = 0;*/

	/*----------------------------------------------------------------------------
	 Close the geolocation file
	 ----------------------------------------------------------------------------*/

	// Close and release file.
	h5status = H5Fclose (fileID);

	/*----------------------------------------------------------------------------
	 TODO : Calculate a bunch of stuff for the geolocation.  This is required for 
            the GDAS code to be able to find the correct subset of the gdas files, 
            and grid the associated datasets.

            This is a duplicate of the MODIS code. Fold this into a function to 
            be used for multiple sensors.
	 ----------------------------------------------------------------------------*/


	// Determine the minimum and maximum lat and long

	float minlat, maxlat, minlon, maxlon;
	int8 scale_lon_flg;

	// TODO : Check for other error flags in geolocation

	minlat = 9999.0;
	maxlat = -9999.0;
	for (idx = 0; idx < npts; idx++) {
		if ((imgr1_ptr->lat[idx] - VDNE_FLOAT32_FILL)>1.e-3) {
			minlat = min(imgr1_ptr->lat[idx],minlat);
			maxlat = max(imgr1_ptr->lat[idx],maxlat);
		}
	}

	minlon = 9999.0;
	maxlon = -9999.0;
	for (idx = 0; idx < npts; idx++) {
		if ((imgr1_ptr->lon[idx] - VDNE_FLOAT32_FILL)>1.e-3) {
			minlon = min(imgr1_ptr->lon[idx],minlon);
			maxlon = max(imgr1_ptr->lon[idx],maxlon);
		}
	}

	printf("\nGranule minlat,maxlat %f %f\n", minlat, maxlat);
	printf("Granule minlon,maxlon %f %f\n", minlon, maxlon);


	// This is the well behaved condition
	if(minlon < -179.0 && maxlon > 179.0 && maxlat < 75.0 && minlat > -75.0) {
		float splon=1000.0, lnlon=-1000.0;
		for (i = 0; i < imgr1_ptr->npts; i++) {
		   if(imgr1_ptr->lon[i] > 0.0) {
			  if(imgr1_ptr->lon[i] < splon) splon = imgr1_ptr->lon[i];
		   }
		   if(imgr1_ptr->lon[i] < 0.0) {
			  if(imgr1_ptr->lon[i] > lnlon) lnlon = imgr1_ptr->lon[i];
		   }
		}
		maxlon = splon;
		minlon = lnlon;
	}

	// FIXME : Checks whether this granule crosses the international dateline
	scale_lon_flg = NO;
	if (minlat < 80.0 && maxlat > -80.0 && minlon < -60.0 && maxlon > 60.0) {
		printf("scale_lon_flag triggered.\n");
		scale_lon_flg = YES;
		for (i = 0; i < imgr1_ptr->npts; i++) {
			if (imgr1_ptr->lon[i] < 0.0)
				imgr1_ptr->lon[i] += 360.0;
		}
	}
	if (scale_lon_flg == YES) {
		for (i = 0; i < imgr1_ptr->npts; i++) {
			if (imgr1_ptr->lon[i] > 180.0)
				imgr1_ptr->lon[i] -= 360.0;
		}
	}

	// Set min/max lat and lon limit array.

	imgr1_ptr->bounds[0] = minlat, imgr1_ptr->bounds[1] = minlon;
	imgr1_ptr->bounds[2] = maxlat, imgr1_ptr->bounds[3] = maxlon;
	if (minlon > 180.0)
		minlon -= 360.0;
	if (maxlon > 180.0)
		maxlon -= 360.0;

	printf("New minlat,maxlat %f %f\n", minlat, maxlat);
	printf("New minlon,maxlon %f %f\n", minlon, maxlon);

	/*----------------------------------------------------------------------------
	 Calculate trigometric parameters.
	 ----------------------------------------------------------------------------*/

	if ((imgr1_ptr->cos_satzen = (float *) malloc(imgr1_ptr->npts * sizeof(float)))
			== NULL)error_allo(rout,"imgr1_ptr->cos_satzen");

	if (verbose == YES)
		fprintf(stdout,"\n%sCalculating scattering and glint angles\n\n",
				EXE_PROMPT);

	for (i = 0; i < imgr1_ptr->npts; i++) {
		imgr1_ptr->cos_satzen[i] = cos(DTOR*imgr1_ptr->satzen[i]);
	}

	/*****************************************************************************
     *****************************************************************************
     *	       Digital Elevation Model: Land Sea Mask & Terrain Height           *
     *****************************************************************************
     *****************************************************************************/

	char *DEM_fileNames[6] = {"dem30ARC_W180N90.hdf","dem30ARC_W60N90.hdf","dem30ARC_E60N90.hdf",
                                  "dem30ARC_W180N0.hdf","dem30ARC_W60N0.hdf","dem30ARC_E60N0.hdf"};
	float DEM_startLat[6] = {90.,90.,90.,0.,0.,0.};
	float DEM_endLat[6] = {0.,0.,0.,-90.,-90.,-90.};
	float DEM_startLon[6] = {-180.,-60.,60.,-180.,-60.,60.};
	float DEM_endLon[6] = {-60.,60.,180.,-60.,60.,180.};
	float granuleCnrLats[4];
	float granuleCnrLons[4];

	int DEM_neededTiles[6] = {0,0,0,0,0,0};
	int DEM_pixel2Tile[npts];
	int cnrIdx,tileIdx;
	int DEM_nrows,DEM_ncols;
	long DEM_npts;

	DEM_nrows = 10800;
	DEM_ncols = 14400;
	DEM_npts = DEM_nrows * DEM_ncols;

	granuleCnrLats[0] = imgr1_ptr->lat[firstGoodRow*ncols];
	granuleCnrLats[1] = imgr1_ptr->lat[ncols * (firstGoodRow + 1) - 1];
	granuleCnrLats[2] = imgr1_ptr->lat[lastGoodRow*ncols];
	granuleCnrLats[3] = imgr1_ptr->lat[ncols * (lastGoodRow + 1) - 1];
	granuleCnrLons[0] = imgr1_ptr->lon[firstGoodRow*ncols];
	granuleCnrLons[1] = imgr1_ptr->lon[ncols * (firstGoodRow + 1) - 1];
	granuleCnrLons[2] = imgr1_ptr->lon[lastGoodRow*ncols];
	granuleCnrLons[3] = imgr1_ptr->lon[ncols * (lastGoodRow + 1) - 1];

	// Determine which of the DEM files are needed
	for (cnrIdx=0;cnrIdx<4;cnrIdx++){
		printf("Lat/Long of point %d is: (%f,%f)\n",cnrIdx,granuleCnrLats[cnrIdx],granuleCnrLons[cnrIdx]);
		for (tileIdx=0;tileIdx<6;tileIdx++){
			if ((granuleCnrLats[cnrIdx]<=DEM_startLat[tileIdx]) 
					&& (granuleCnrLats[cnrIdx]>DEM_endLat[tileIdx])
					&& (granuleCnrLons[cnrIdx]>=DEM_startLon[tileIdx]) 
					&& (granuleCnrLons[cnrIdx]<DEM_endLon[tileIdx])){
				printf("\tGranule point is in Tile %d\n",tileIdx);
				DEM_neededTiles[tileIdx] = 1;
			}
		}
	}

	// Assign the appropriate tile to each pixel
	printf("Assign the appropriate tile to each pixel\n");
	for (scan=0;scan<nViirsScans;scan++){
		if (imgr1_ptr->ModeScan[scan] <= 2){
			for (detector=0;detector<ndets;detector++){
				for (col=0;col<ncols;col++){
					idx = scan*ndets*ncols + detector*ncols + col;
					for (tileIdx=0;tileIdx<6;tileIdx++){
						if ((imgr1_ptr->lat[idx]<=DEM_startLat[tileIdx])
						    && (imgr1_ptr->lat[idx]>DEM_endLat[tileIdx])
						    && (imgr1_ptr->lon[idx]>=DEM_startLon[tileIdx])
						    && (imgr1_ptr->lon[idx]<DEM_endLon[tileIdx]))
						{
							DEM_pixel2Tile[idx] = tileIdx;	
						}
					} // End tile loop
				} // End column loop
			} // End detector loop
		}else{
			for (detector=0;detector<ndets;detector++){
				for (col=0;col<ncols;col++){
					idx = scan*ndets*ncols + detector*ncols + col;
					DEM_pixel2Tile[idx] = -1;	
				} // End column loop
			} // End detector loop
		} // End ModeScan if-block
	} // End scan loop

	printf("\nThe needed tiles are : ");
	for (tileIdx=0;tileIdx<6;tileIdx++)printf(" %d ",DEM_neededTiles[tileIdx]);
	printf("\n");

	// Open DEM files and read in contents

	int DEM_fileID[6], LSM_dsetIDx[6], LSM_dsetID[6],TerrHgt_dsetIDx[6], TerrHgt_dsetID[6];
	int start[2], edge[2];
	int h4status;

	// Open the file and get the file id
	for (tileIdx=0;tileIdx<6;tileIdx++){
		DEM_fileID[tileIdx] = -1;
		if (DEM_neededTiles[tileIdx]){
			sprintf(full_filename, "%s/%s/%s", imgr1_ptr->sfc_dir_name,"LSM",DEM_fileNames[tileIdx]);
			printf("\nOpening LSM file %s\n",full_filename);
			DEM_fileID[tileIdx] = SDstart(full_filename, DFACC_READ);
			if (DEM_fileID[tileIdx] == FAIL) {
				fprintf(stderr,"%sInvalid HDF file, %s\n",EXE_PROMPT,DEM_fileNames[tileIdx]);
				exit (EXIT_FAILURE);
			}
		}
	}

	char * LSM_dsetName = "LandWater";
	char * TerrHgt_dsetName = "Elevation";

	// Get the index of the dataset in the file from the dataset name.
	for (tileIdx=0;tileIdx<6;tileIdx++){
		LSM_dsetIDx[tileIdx] = -1;
		if (DEM_fileID[tileIdx] != -1){
			LSM_dsetIDx[tileIdx] = SDnametoindex(DEM_fileID[tileIdx],LSM_dsetName);
			TerrHgt_dsetIDx[tileIdx] = SDnametoindex(DEM_fileID[tileIdx],TerrHgt_dsetName);
		}
	}

	// Connect to dataset and get the dataset id
	for (tileIdx=0;tileIdx<6;tileIdx++){
		LSM_dsetID[tileIdx] = -1;
		TerrHgt_dsetID[tileIdx] = -1;
		if ((LSM_dsetIDx[tileIdx] != -1) && (TerrHgt_dsetIDx[tileIdx] != -1)){
			printf("\nConnecting to datasets in DEM file %s\n",DEM_fileNames[tileIdx]);
			LSM_dsetID[tileIdx] = SDselect(DEM_fileID[tileIdx],LSM_dsetIDx[tileIdx]);
			TerrHgt_dsetID[tileIdx] = SDselect(DEM_fileID[tileIdx],TerrHgt_dsetIDx[tileIdx]);
			if (LSM_dsetID[tileIdx] == FAIL) {
				fprintf(stderr,"%s%s-cannot select specified SDS, =%s\n",EXE_PROMPT,rout,LSM_dsetName);
				exit(EXIT_FAILURE);
			}
			if (TerrHgt_dsetID[tileIdx] == FAIL) {
				fprintf(stderr,"%s%s-cannot select specified SDS, =%s\n",EXE_PROMPT,rout,TerrHgt_dsetName);
				exit(EXIT_FAILURE);
			}
		}
	}

	// Reserve memory for DEM array of tile pointers, and initialise these
	// pointers

	int idx_free;

	int8 **LSM;
	if ((LSM = (int8 **) malloc(6 * sizeof(int8 *)))==NULL)
		error_allo(rout,"LSM");

	int16 **TerrHgt;
	if ((TerrHgt = (int16 **) malloc(6 * sizeof(int16 *)))==NULL)
		error_allo(rout,"TerrHgt");

	for (tileIdx=0;tileIdx<6;tileIdx++){
		LSM[tileIdx] = 0;
		TerrHgt[tileIdx] = 0;
		if (DEM_neededTiles[tileIdx]){
			if((LSM[tileIdx] = (int8 *) malloc(DEM_npts * sizeof(int8)))==NULL)
				printf("\nFailed to allocate memory for LSM tile %d\n",tileIdx);
			if((TerrHgt[tileIdx] = (int16 *) malloc(DEM_npts * sizeof(int16)))==NULL)
				printf("\nFailed to allocate memory for TerrHgt tile %d\n",tileIdx);
		}
	}

	// Offset to start reading dataset
	start[0] = 0;
	start[1] = 0;
	// Number of elements in each dimension to read
	edge[0] = DEM_nrows;
	edge[1] = DEM_ncols;

	for (tileIdx=0;tileIdx<6;tileIdx++){
		if (LSM_dsetID[tileIdx] != -1){
			printf("\tReading in dataset from DEM file %s\n",DEM_fileNames[tileIdx]);
			printf("\tDataset ID for file is %s is %d\n",DEM_fileNames[tileIdx],LSM_dsetID[tileIdx]);
			h4status = SDreaddata(LSM_dsetID[tileIdx],start,NULL,edge,LSM[tileIdx]);
			h4status = SDreaddata(TerrHgt_dsetID[tileIdx],start,NULL,edge,TerrHgt[tileIdx]);
			if (h4status == FAIL) {
				fprintf(stderr,"%s%s-cannot read specified SDS %s, freeing memory.\n",EXE_PROMPT,rout,LSM_dsetName);
				for (idx_free=0;idx_free<6;idx_free++){
					if (LSM[idx_free] != NULL)
						printf("\nTile %d is to be freed\n",idx_free);
						free(LSM[idx_free]);
						free(TerrHgt[idx_free]);
				}
				free(LSM);
				free(TerrHgt);
				exit(EXIT_FAILURE);
			}else{
				printf("\tSuccessfully read in tile %d DEM file %s\n\n",tileIdx,DEM_fileNames[tileIdx]);
			}
		}
	}   

	//
	// Granulate the LSM and TerrHgt datasets
	//

	printf("Granulating the LSM and TerrHgt datasets...\n");

	long gridRow, gridCol, gridIdx;
	float DEM_dlat, DEM_dlon;

	// Allocate memory for the LSM
	if ((imgr1_ptr->landmask = (unsigned char *) malloc(npts*sizeof(unsigned char))) == NULL)
		error_allo(rout,"imgr1_ptr->landmask");
	// Allocate memory for the TerrHgt
	if ((imgr1_ptr->zsfc = (float *) malloc(npts*sizeof(float))) == NULL)
		error_allo(rout,"imgr1_ptr->zsfc");

	// Grid increments in degrees
	DEM_dlat = 0.5/60.;
	DEM_dlon = 0.5/60.;

	// Determine the grid cell for each pixel, and assign the LSM for that cell to the 
	// corrsponding pixel.
	/* TODO : Check if we need to make any tweaks to account for 
	   TODO : dataline crossings and hight latitude granules */

	for (idx=0; idx<npts; idx++) {
		if (DEM_pixel2Tile[idx] != -1){ 	
			tileIdx = DEM_pixel2Tile[idx];

			////// FIXME : Need to account for reduced rows perhaps? ///////////
			gridRow = max(0,min(DEM_nrows-1, (long) (fabs(imgr1_ptr->lat[idx]-DEM_startLat[tileIdx])/DEM_dlat)));
			gridCol = max(0,min(DEM_ncols-1, (long) (fabs(imgr1_ptr->lon[idx]-DEM_startLon[tileIdx])/DEM_dlon)));
			////////////////////////////////////////////////////////////////////

			gridIdx = (DEM_ncols*gridRow) + gridCol;

			imgr1_ptr->landmask[idx] = LSM[tileIdx][gridIdx];
			imgr1_ptr->zsfc[idx] = (float) TerrHgt[tileIdx][gridIdx];
		}else{
			imgr1_ptr->landmask[idx] = UINT8_FILL_TEST;
			imgr1_ptr->zsfc[idx] = FLOAT32_FILL_TEST;
		}
	}
	printf("   success\n");

	////////////////////////////////////////////////////////////////////////
	// Output land/sea mask to a h5 file as a test
	fileID = H5Fcreate ("LEOCAT-DEM_LSM.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	dspaceID = H5Screate_simple (2, dimSize, NULL);
	dsetID = H5Dcreate2 (fileID, "LandSeaMask",H5T_NATIVE_UCHAR, dspaceID, H5P_DEFAULT,
				H5P_DEFAULT, H5P_DEFAULT);
	h5status = H5Dwrite (dsetID, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
				imgr1_ptr->landmask);
	h5status = H5Dclose (dsetID);
	h5status = H5Sclose (dspaceID);
	h5status = H5Fclose (fileID);
	///////////////////////////////////////////////////////////////////////
	// Output terrain height to a h5 file as a test
	fileID = H5Fcreate ("LEOCAT-DEM_TerrHgt.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	dspaceID = H5Screate_simple (2, dimSize, NULL);
	dsetID = H5Dcreate2 (fileID, "TerrainHeight",H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
				H5P_DEFAULT, H5P_DEFAULT);
	h5status = H5Dwrite (dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
				imgr1_ptr->zsfc);
	h5status = H5Dclose (dsetID);
	h5status = H5Sclose (dspaceID);
	h5status = H5Fclose (fileID);
	///////////////////////////////////////////////////////////////////////

	// End access to dataset
	for (tileIdx=0;tileIdx<6;tileIdx++){
		if (LSM_dsetID[tileIdx] != -1){
			printf("\nEnding access to dataset in LSM file %s\n",DEM_fileNames[tileIdx]);
			h4status = SDendaccess(LSM_dsetID[tileIdx]);
			h4status = SDendaccess(TerrHgt_dsetID[tileIdx]);
			if (h4status == FAIL) {
				fprintf(stderr,"%s%s-cannot end access to specified SDS, =%s\n",EXE_PROMPT,rout,LSM_dsetName);
				for (idx_free=0;idx_free<6;idx_free++){
					if (LSM[idx_free] != NULL){
						printf("\nTile %d is to be freed\n",idx_free);
						free(LSM[idx_free]);
						free(TerrHgt[idx_free]);
					}
				}
				free(LSM);
				free(TerrHgt);
				exit(EXIT_FAILURE);
			}
		}
	}

	// Close file(s)
	for (tileIdx=0;tileIdx<6;tileIdx++){
		if (DEM_fileID[tileIdx] != -1){
			printf("\nClosing DEM file %s\n",DEM_fileNames[tileIdx]);
			h4status = SDend(DEM_fileID[tileIdx]);
			if (h4status == FAIL) {
				fprintf(stderr,"%s%s-cannot close %s\n",EXE_PROMPT,rout,DEM_fileNames[tileIdx]);
				for (idx_free=0;idx_free<6;idx_free++){
					if (LSM[idx_free] != NULL)
						printf("\nTile %d is to be freed\n",idx_free);
						free(LSM[idx_free]);
						free(TerrHgt[idx_free]);
				}
				free(LSM);
				free(TerrHgt);
				exit(EXIT_FAILURE);
			}
		}
	}

	// Free memory from LSM and TerrHgt

	for (idx_free=0;idx_free<6;idx_free++){
		if (LSM[idx_free] != NULL)
			printf("\nTile %d is to be freed\n",idx_free);
			free(LSM[idx_free]);
			free(TerrHgt[idx_free]);
	}
	free(LSM);
	free(TerrHgt);

	printf("\nEnding granulation of DEM\n\n");


	/*****************************************************************************
     *****************************************************************************
     *	              Surface type from the IGBP Ecosystem Map                   *
     *****************************************************************************
     *****************************************************************************/

	char * IGBP_filename = "IGBP.EcoMap.v1.0.2004.129.v004.hdf";
	char * IGBP_dsetName = "IGBP_Land_Cover_Type";

	long IGBP_nrows = 10800;
	long IGBP_ncols = 21600;
	double LatScale[IGBP_nrows],LonScale[IGBP_ncols];
	long IGBP_npts = IGBP_nrows * IGBP_ncols;
	long IGBP_granNumRows,IGBP_granNumCols;
	long IGBP_granMinRow, IGBP_granMaxRow;
	long IGBP_granMinCol, IGBP_granMaxCol;
	long IGBP_granNumPts;
	int IGBP_fileID, IGBP_dsetIDx, IGBP_dsetID;

	// Grid increments in degrees
	double IGBP_dlat, IGBP_dlon;
	IGBP_dlat = 1./60.;
	IGBP_dlon = 1./60.;

	// Determine the IGBP indices which bound the granule
	for(gridRow=0;gridRow<IGBP_nrows;gridRow++){
		LatScale[gridRow] = 90. - (IGBP_dlat/2.) - ((double)gridRow) * IGBP_dlat;
	}

	printf("\nMaximum LatScale = %13.8f\n",LatScale[0]);
	printf("Minimum LatScale = %13.8f\n",LatScale[IGBP_nrows-1]);

	printf("\nMaximum Granule Latitude = %13.8f\n",maxlat);
	for(gridRow=0;gridRow<IGBP_nrows;gridRow++){
		if(LatScale[gridRow] < maxlat){
			printf("LatScale dropped below maxlat at %13.8f\n",LatScale[gridRow]);
			IGBP_granMinRow = gridRow-1;
			break;
		}
	}
	printf("Final maxlat gridRow = %ld\n",IGBP_granMinRow);
	printf("Max LatScale = %13.8f\n",LatScale[IGBP_granMinRow]);

	printf("\nMinimum Granule Latitude = %13.8f\n",minlat);
	for(gridRow=start[0]+1;gridRow<IGBP_nrows;gridRow++){
		if(LatScale[gridRow] < minlat){
			printf("LatScale below minlat at %13.8f\n",LatScale[gridRow]);
			IGBP_granMaxRow = gridRow;
			break;
		}
	}
	printf("Final minlat gridRow = %ld\n",IGBP_granMaxRow);
	printf("Min LatScale = %13.8f\n\n",LatScale[IGBP_granMaxRow]);

	//////////

	for(gridCol=0;gridCol<IGBP_ncols;gridCol++){
		LonScale[gridCol] = ((double)gridCol) * IGBP_dlon - 180. + (IGBP_dlon/2.);
	}

	printf("\nMinimum Granule Longitude = %13.8f\n",minlon);
	for(gridCol=0;gridCol<IGBP_ncols;gridCol++){
		if(LonScale[gridCol] > minlon){
			printf("LonScale above minlon at %13.8f\n",LonScale[gridCol]);
			IGBP_granMinCol = gridCol-1;
			break;
		}
	}
	printf("Final minlon gridCol = %ld\n",IGBP_granMinCol);
	printf("Min LonScale = %13.8f\n",LonScale[IGBP_granMinCol]);

	printf("\nMaximum Granule Longitude = %13.8f\n",maxlon);
	for(gridCol=start[1]+1;gridCol<IGBP_ncols;gridCol++){
		if(LonScale[gridCol] > maxlon){
			printf("LonScale above maxlon at %13.8f\n",LonScale[gridCol]);
			IGBP_granMaxCol = gridCol;
			break;
		}
	}
	printf("Final maxlon gridCol = %ld\n",IGBP_granMaxCol);
	printf("Max LonScale = %13.8f\n",LonScale[IGBP_granMaxCol]);

	// Reserve memory for IGBP array
	int8 *IGBP;
	IGBP_granNumRows = IGBP_granMaxRow - IGBP_granMinRow + 1;
	IGBP_granNumCols = IGBP_granMaxCol - IGBP_granMinCol + 1;
	IGBP_granNumPts = IGBP_granNumRows * IGBP_granNumCols;

	if((IGBP = (int8 *) malloc(IGBP_granNumPts * sizeof(int8)))==NULL)
		printf("\nFailed to allocate memory for IGBP file\n");

	// Set dataset boundaries and extents
	start[0] = IGBP_granMinRow;
	edge[0] = IGBP_granNumRows;
	start[1] = IGBP_granMinCol;
	edge[1] = IGBP_granNumCols;

	// Open the file and get the file id
	IGBP_fileID = -1;
	sprintf(full_filename, "%s/%s", imgr1_ptr->sfc_dir_name,IGBP_filename);
	printf("\nOpening IGBP file %s\n",full_filename);
	IGBP_fileID = SDstart(full_filename, DFACC_READ);
	if (IGBP_fileID == FAIL) {
		fprintf(stderr,"%sInvalid HDF file, %s\n",EXE_PROMPT,IGBP_filename);
		exit (EXIT_FAILURE);
	}
	printf("\nIGBP_fileID = %d\n",IGBP_fileID);

	// Get the index of the dataset in the file from the dataset name.
	IGBP_dsetIDx = -1;
	if (IGBP_fileID != -1){
		IGBP_dsetIDx = SDnametoindex(IGBP_fileID,IGBP_dsetName);
		printf("\nIGBP_dsetIDx = %d\n",IGBP_dsetIDx);
	}

	// Connect to dataset and get the dataset id
	IGBP_dsetID = -1;
	if (IGBP_dsetIDx != -1){
		printf("\nConnecting to dataset in IGBP file %s\n",IGBP_filename);
		IGBP_dsetID = SDselect(IGBP_fileID,IGBP_dsetIDx);
		if (IGBP_dsetID == FAIL) {
			fprintf(stderr,"%s%s-cannot select specified SDS, =%s\n",EXE_PROMPT,rout,IGBP_dsetName);
			exit(EXIT_FAILURE);
		}
		printf("\nIGBP_dsetID = %d\n",IGBP_dsetID);
	}

	// Read the IGBP dataset
	printf("\nDataset ID for file is %s is %d\n",IGBP_filename,IGBP_dsetID);
	h4status = SDreaddata(IGBP_dsetID,start,NULL,edge, (VOIDP) IGBP);
	if (h4status == FAIL) {
		fprintf(stderr,"%s%s-cannot read specified SDS %s, freeing memory.\n",EXE_PROMPT,rout,IGBP_dsetName);
		free(IGBP);
		exit(EXIT_FAILURE);
	}

	// End access to dataset
	printf("\nEnding access to dataset in IGBP file %s\n",IGBP_filename);
	h4status = SDendaccess(IGBP_dsetID);
	if (h4status == FAIL) {
		fprintf(stderr,"%s%s-cannot end access to specified SDS, =%s\n",EXE_PROMPT,rout,IGBP_dsetName);
		free(IGBP);
		exit(EXIT_FAILURE);
	}

	// Close file(s)
	printf("\nClosing IGBP file %s\n",IGBP_filename);
	h4status = SDend(IGBP_fileID);
	if (h4status == FAIL) {
		fprintf(stderr,"%s%s-cannot close %s\n",EXE_PROMPT,rout,IGBP_filename);
		free(IGBP);
		exit(EXIT_FAILURE);
	}

	//
	// Granulate the IGBP file
	//

	// Allocate memory for the surface type
	if ((imgr1_ptr->sfc_type = (unsigned char *) malloc(npts*sizeof(unsigned char))) == NULL)
		error_allo(rout,"imgr1_ptr->sfc_type");

	// Determine the grid cell for each pixel, and assign the LSM for that cell to the 
	// corrsponding pixel.
	/* TODO : Check if we need to make any tweaks to account for
	   TODO : dataline crossings and hight latitude granules
	   TODO : for import of IGBP surface type.               */

	for (idx=0; idx<npts; idx++) {
		if (fabs(imgr1_ptr->lat[idx]-VDNE_FLOAT32_FILL) < 0.1
		 || fabs(imgr1_ptr->lon[idx]-VDNE_FLOAT32_FILL) < 0.1){

			imgr1_ptr->sfc_type[idx] = UINT8_FILL_TEST;

		}else{
			gridRow = max(0,min(IGBP_granNumRows-1, (long) (fabs(imgr1_ptr->lat[idx]-LatScale[IGBP_granMinRow])/IGBP_dlat)));
			gridCol = max(0,min(IGBP_granNumCols-1, (long) (fabs(imgr1_ptr->lon[idx]-LonScale[IGBP_granMinCol])/IGBP_dlon)));

			gridIdx = IGBP_granNumCols*gridRow + gridCol;
			imgr1_ptr->sfc_type[idx] = IGBP[gridIdx];
		}
	}

	////////////////////////////////////////////////////////////////////////
	// Output IGBP surface type to a h5 file as a test
	fileID = H5Fcreate ("LEOCAT-IGBP_QST.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	dspaceID = H5Screate_simple (2, dimSize, NULL);
	dsetID = H5Dcreate2 (fileID, IGBP_dsetName,H5T_NATIVE_UCHAR, dspaceID, H5P_DEFAULT,
				H5P_DEFAULT, H5P_DEFAULT);
	h5status = H5Dwrite (dsetID, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
				imgr1_ptr->sfc_type);
	h5status = H5Dclose (dsetID);
	h5status = H5Sclose (dspaceID);
	h5status = H5Fclose (fileID);
	///////////////////////////////////////////////////////////////////////

	// Free memory from IGBP
	free(IGBP);

	printf("\nEnding granulation of IGBP\n\n");

	/*****************************************************************************
     *****************************************************************************
     *	      VIIRS Radiances, Reflectances, and Brightness Temperatures         *
     *****************************************************************************
     *****************************************************************************/

	/*----------------------------------------------------------------------------
	 Determine how many spectral channels are to be read in.
	 ----------------------------------------------------------------------------*/

	int channel;

	nchan = 0;
	for (channel = 0; channel < NVIIRS_CHAN; channel++) {
		if (chflg[channel])
			nchan++;
	}
	if (verbose == YES)
		fprintf(stdout,"\n\n%s%d channels are to be read.\n\n", EXE_PROMPT,nchan);

	/*----------------------------------------------------------------------------
	 Construct the full set of level 1 filenames for the VIIRS granules

	 TODO : Add ability to specify granule filenames from command line
	 ----------------------------------------------------------------------------*/

	char *level1ChanPrefix[NVIIRS_CHAN] = { "SVDNB_npp_","SVI01_npp_","SVI02_npp_","SVI03_npp_","SVI04_npp_",
		"SVI05_npp_","SVM01_npp_","SVM02_npp_","SVM03_npp_","SVM04_npp_","SVM05_npp_",
		"SVM06_npp_","SVM07_npp_","SVM08_npp_","SVM09_npp_","SVM10_npp_","SVM11_npp_",
		"SVM12_npp_","SVM13_npp_","SVM14_npp_","SVM15_npp_","SVM16_npp_"};

	char level1FileNames[NVIIRS_CHAN][256] = {"\0","\0","\0","\0","\0","\0","\0","\0",
		"\0","\0","\0","\0","\0","\0","\0","\0","\0","\0","\0","\0","\0","\0"};

	int imageryType[NVIIRS_CHAN] = {DAYNIGHT,
		IMAGERY,IMAGERY,IMAGERY,IMAGERY,IMAGERY,
		MODERATE,MODERATE,MODERATE,MODERATE,MODERATE,
		MODERATE,MODERATE,MODERATE,MODERATE,MODERATE,
		MODERATE,MODERATE,MODERATE,MODERATE,MODERATE,
		MODERATE};

	char *fileString, *dateString, *createString, *level1Suffix;
	unsigned int underScores[8];
	char * char_ptr = 0;
	int charCount;
	char grpName[100];
	float scaleFactors[2];
	hid_t group_ID;

	// Reserve memory for strings
	fileString = (char *) malloc(100 * sizeof(char));
	level1Suffix = (char *) malloc(20 * sizeof(char));

	printf("Level 1 dir name  : %s\n", imgr1_ptr->l1b_dir_name);
	printf("Level 1 file name : %s\n", imgr1_ptr->l1b_filename1);
	sprintf(full_filename, "%s/%s", imgr1_ptr->l1b_dir_name, imgr1_ptr->l1b_filename1);
	printf("Full Level 1 name %s\n", full_filename);

	// Copy the level1 filename from the imgr1 structure.	
	charCount = sprintf(fileString, "%s",imgr1_ptr->l1b_filename1);

	// Determine where the underscores are in the filename
	ichar = 0;
	char_ptr = strchr(fileString,'_');
	while (char_ptr != NULL)
	{
		underScores[ichar] = (char_ptr - fileString);
		ichar++;
		char_ptr=strchr(char_ptr+1,'_');
	}

	// Determine the file suffix
	ichar = 0;
	char_ptr = fileString + underScores[6];
	charCount = sprintf(level1Suffix, "%s",char_ptr);
	printf("level1Suffix = %s\n",level1Suffix);

	// Get the creation date string from the filename
	createString = fileString + underScores[5] + 1;
	*(fileString + underScores[6]) = (int)NULL;

	printf("Creation Dates string = %s\n",createString);

	// Get the granule date string from the filename
	dateString = fileString + underScores[1] + 1;
	*(dateString + underScores[5]-underScores[1]-1) = (int)NULL;

	printf("Datestring = %s\n",dateString);


	/*----------------------------------------------------------------------------
	 Determine the filenames for each channel file.         
	 ----------------------------------------------------------------------------*/

	char globString[256];
	glob_t globbuf;

	for (channel=0;channel<NVIIRS_CHAN;channel++){
		if (chflg[channel]){		// This channel required?

			// Construct the file pattern glob, for this date and channel, and find any matches.
			sprintf(globString,"%s/%s%s%s%s",imgr1_ptr->l1b_dir_name,level1ChanPrefix[channel],
				dateString,"*",level1Suffix);

			if ((glob(globString, GLOB_MARK, NULL, &globbuf)) == GLOB_NOMATCH){
				printf("\tThere were no matches found for channel %d in %s\n",channel,imgr1_ptr->l1b_dir_name);
				return;
			}

			/*printf("\tNumber of Matches: %zd\n",globbuf.gl_pathc);*/

			// Take the first filename match
			sprintf(level1FileNames[channel],"%s",globbuf.gl_pathv[0]);
			globfree (&globbuf);
		}
	}
	// Free the allocated memory for the filename arrays.
	free(level1Suffix);
	free(fileString);

	/*----------------------------------------------------------------------------
	    Read in the required channels

		TODO : Change algorithm structure to allow specification of which data
               (radiance, reflectance, BT) are required for each channel.
	 ----------------------------------------------------------------------------*/

	int scaleDataSet = 0;
	long npts_channel;
	unsigned short *bufferPtr_ushort;

	// Define pointer to an array of channel structure pointers
	if ((imgr1_ptr->channel_ptr = (imagerL1_chanType **) malloc( NVIIRS_CHAN * sizeof(imagerL1_chanType *))) == NULL)
		error_allo(rout,"imgr1_ptr->channel_ptr");

	// Set all channel pointers to NULL
	for (channel=0;channel<NVIIRS_CHAN;channel++)
		imgr1_ptr->channel_ptr[channel] = 0;

	// Initialise the number of points in a channel dataset to the number of MOD
	// resolution points
	npts_channel = npts;

	// Initialise a channel structure and return a pointer to it for each valid channel, 
	// and return a null pointer for the rest.
	for (channel=0;channel<NVIIRS_CHAN;channel++){

		printf("\nChannel %3d: \n",channel);
		if (chflg[channel]){		// This channel required?

			// Initialise a channel structure for each channel, and point to it with imgr1_ptr->channel_ptr[channel]
			if ((imgr1_ptr->channel_ptr[channel] = (imagerL1_chanType *) malloc(sizeof(imagerL1_chanType))) == NULL)
				error_allo(rout,"imgr1_ptr->channel_ptr[channel]");

			// Figure out the imagery type...
			switch (imageryType[channel]){
				case DAYNIGHT:
					/*printf("\tThis is the Day-Night Band\n");*/
					sprintf(grpName,"%s","/All_Data/VIIRS-DNB-SDR_All");
					break;
				case IMAGERY:
					/*printf("\tThis is the Imagery Resolution Band\n");*/
					sprintf(grpName,"%s%d%s","/All_Data/VIIRS-I",channel,"-SDR_All");
					break;
				case MODERATE:
					/*printf("\tThis is the Moderate Resolution Band\n");*/
					sprintf(grpName,"%s%d%s","/All_Data/VIIRS-M",channel-5,"-SDR_All");
					break;
			}
			/*printf("\tDataset pathname: %s\n",grpName);*/

			// Open the file as Read Only, and default File Access property list
			/*printf("\tOpening the file %s\n",level1FileNames[channel]);*/
			fileID = H5Fopen(level1FileNames[channel], H5F_ACC_RDONLY, H5P_DEFAULT);
			if (fileID == FAIL) {
				fprintf(stderr,"%sInvalid HDF file, %s\n", EXE_PROMPT,level1FileNames[channel]);
				exit(EXIT_FAILURE);
			}

			// Get the group id of the VIIRS SDR group
			if ((group_ID = H5Gopen2 (fileID, grpName, H5P_DEFAULT)) < 0){
				printf(">>> Unable to find the group %s in file %s\n",grpName,level1FileNames[channel]);
				h5status = H5Fclose (fileID);
				if (h5status == 0) printf("\tSuccessfully closed file.\n");
				return;
			}

			/////////////////////////////////////////////////////////////////////
			//  Try getting the radiance, and the number of rows, columns and  //
			//  number of points for this dataset.                             //
			/////////////////////////////////////////////////////////////////////

			sprintf(dsetName,"%s%s",grpName,"/Radiance");
			printf("\tReading in channel %d radiance...\n",channel);

			h5status = H5LTfind_dataset ( group_ID, "Radiance" );
			/*printf("\tThe return value for finding Radiance is %d\n",h5status);*/

			if (h5status){
				h5status = viirs_level1_read (
					dsetName,
					fileID, 
					(void **) &bufferPtr_void,
					&nRank,
					&imgr1_ptr->channel_ptr[channel]->nrows,
					&imgr1_ptr->channel_ptr[channel]->ncols
				);
				/*printf("\tThe return value for reading the Radiance is %d\n",h5status);*/
				imgr1_ptr->channel_ptr[channel]->npts = 
					imgr1_ptr->channel_ptr[channel]->nrows
					* imgr1_ptr->channel_ptr[channel]->ncols;

				npts_channel = imgr1_ptr->channel_ptr[channel]->npts;

				// If the radiance scale factors exist, get the factors
				scaleDataSet = 0;
				h5status = H5LTfind_dataset ( group_ID, "RadianceFactors" );
				if (h5status){
					scaleDataSet = 1;
					sprintf(dsetName,"%s%s",grpName,"/RadianceFactors");
					dsetID = H5Dopen2 (fileID, dsetName, H5P_DEFAULT);
					h5status = H5Dread (dsetID,              // dataset ID
						H5T_NATIVE_FLOAT,    // ident of memory datatype
						H5S_ALL,             // ident of memory dataspace
						H5S_ALL,             // ident of dataset's dataspace in file
						H5P_DEFAULT,         // File access properties
						scaleFactors);        // Pointer to buffer receiving data

						/*printf("\tRadiance Scale Factor = %f\n",scaleFactors[0]);*/
						/*printf("\tRadiance Offset = %f\n",scaleFactors[1]);*/
				}

				// Allocate memory for the radiance if we need to unscale, or point the radiance
				// pointer to the void pointer.
				bufferPtr_ushort = 0;
				if (scaleDataSet){
					if ((imgr1_ptr->channel_ptr[channel]->radiance = (float *) malloc(npts_channel*sizeof(float))) == NULL)
						error_allo(rout,"radiance allocation failed");
					bufferPtr_ushort = (unsigned short *) bufferPtr_void;
					bufferPtr_void = 0;
				}else{
					imgr1_ptr->channel_ptr[channel]->radiance = (float *) bufferPtr_void;
					bufferPtr_void = 0;
				}
				
				// Print out the initial read of the radiance dataset.
				/*printf("\n\tOriginal Radiance =\n");*/
				/*for (i=0; i<10; i++) {*/
					/*printf ("\t [ ");*/
					/*for (j=0; j<10; j++){*/
						/*idx = i * imgr1_ptr->channel_ptr[channel]->ncols + j ;*/
						/*if (scaleDataSet){*/
							/*printf (" %u", bufferPtr_ushort[idx]);*/
						/*}else{*/
							/*printf (" %11.6f", imgr1_ptr->channel_ptr[channel]->radiance[idx]);*/
						/*}*/
					/*}*/
					/*printf ("\t ]\n");*/
				/*}*/
				/*printf("\n\n");*/

				// Scale the radiance, and free the pointer to the scaled Radiance, if we need to
				if (scaleDataSet){
					for (idx = 0; idx < npts_channel; idx++){
						if (bufferPtr_ushort[idx] > UINT16_FILL_TEST){
							imgr1_ptr->channel_ptr[channel]->radiance[idx] = FLOAT32_FILL_TEST;
						}else{
							imgr1_ptr->channel_ptr[channel]->radiance[idx] = 
								(float) bufferPtr_ushort[idx] * scaleFactors[0] + scaleFactors[1];
						}
					}
					free(bufferPtr_ushort);
				}

				// Print out the final radiance dataset.
				/*printf("\n\tFinal Radiance =\n");*/
				/*for (i=0; i<10; i++) {*/
					/*printf ("\t [ ");*/
					/*for (j=0; j<10; j++){*/
						/*idx = i * imgr1_ptr->channel_ptr[channel]->ncols + j ;*/
						/*printf (" %11.6f", imgr1_ptr->channel_ptr[channel]->radiance[idx]);*/
					/*}*/
					/*printf ("\t ]\n");*/
				/*}*/
				/*printf("\n\n");*/


			}else{   // Radiance not here, zero radiance pointer
				imgr1_ptr->channel_ptr[channel]->radiance = 0;
			}

			////////////////////////////////////////////////////////////////////////
			//  Try getting the reflectance, and the number of rows, columns and  //
			//  number of points for this dataset.                                //
			////////////////////////////////////////////////////////////////////////

			sprintf(dsetName,"%s%s",grpName,"/Reflectance");
			printf("\tReading in channel %d reflectance...\n",channel);

			h5status = H5LTfind_dataset ( group_ID, "Reflectance" );
			/*printf("\tThe return value for finding Reflectance is %d\n",h5status);*/

			if (h5status){
				h5status = viirs_level1_read (dsetName, fileID, 
					(void **) &bufferPtr_void,
					&nRank,
					&imgr1_ptr->channel_ptr[channel]->nrows,
					&imgr1_ptr->channel_ptr[channel]->ncols
				);
				/*printf("\tThe return value for reading the Reflectance is %d\n",h5status);*/

				// If the reflectance scale factors exist, get the factors
				scaleDataSet = 0;
				h5status = H5LTfind_dataset ( group_ID, "ReflectanceFactors" );
				if (h5status){
					scaleDataSet = 1;
					sprintf(dsetName,"%s%s",grpName,"/ReflectanceFactors");
					dsetID = H5Dopen2 (fileID, dsetName, H5P_DEFAULT);
					h5status = H5Dread (dsetID,              // dataset ID
						H5T_NATIVE_FLOAT,    // ident of memory datatype
						H5S_ALL,             // ident of memory dataspace
						H5S_ALL,             // ident of dataset's dataspace in file
						H5P_DEFAULT,         // File access properties
						scaleFactors);        // Pointer to buffer receiving data

						/*printf("\tReflectance Scale Factor = %f\n",scaleFactors[0]);*/
						/*printf("\tReflectance Offset = %f\n",scaleFactors[1]);*/
				}

				// Allocate memory for the reflectance if we need to unscale, or point the reflectance
				// pointer to the void pointer.
				bufferPtr_ushort = 0;
				if (scaleDataSet){
					if ((imgr1_ptr->channel_ptr[channel]->reflectance = (float *) malloc(npts_channel*sizeof(float))) == NULL)
						error_allo(rout,"reflectance allocation failed");
					bufferPtr_ushort = (unsigned short *) bufferPtr_void;
					bufferPtr_void = 0;
				}else{
					imgr1_ptr->channel_ptr[channel]->reflectance = (float *) bufferPtr_void;
					bufferPtr_void = 0;
				}
				
				// Print out the initial read of the reflectance dataset.
				/*printf("\n\tOriginal Reflectance =\n");*/
				/*for (i=0; i<10; i++) {*/
					/*printf ("\t [ ");*/
					/*for (j=0; j<10; j++){*/
						/*idx = i * imgr1_ptr->channel_ptr[channel]->ncols + j ;*/
						/*if (scaleDataSet){*/
							/*printf (" %u", bufferPtr_ushort[idx]);*/
						/*}else{*/
							/*printf (" %11.6f", imgr1_ptr->channel_ptr[channel]->reflectance[idx]);*/
						/*}*/
					/*}*/
					/*printf ("\t ]\n");*/
				/*}*/
				/*printf("\n\n");*/

				// Scale the reflectance, and free the pointer to the scaled Reflectance, if we need to
				if (scaleDataSet){
					for (idx = 0; idx < npts_channel; idx++){
						if (bufferPtr_ushort[idx] > UINT16_FILL_TEST){
							imgr1_ptr->channel_ptr[channel]->reflectance[idx] = FLOAT32_FILL_TEST;
						}else{
							imgr1_ptr->channel_ptr[channel]->reflectance[idx] = 
								(float) bufferPtr_ushort[idx] * scaleFactors[0] + scaleFactors[1];
						}
					}
					free(bufferPtr_ushort);
				}

				// Print out the final reflectance dataset.
				/*printf("\n\tFinal Reflectance =\n");*/
				/*for (i=0; i<10; i++) {*/
					/*printf ("\t [ ");*/
					/*for (j=0; j<10; j++){*/
						/*idx = i * imgr1_ptr->channel_ptr[channel]->ncols + j ;*/
						/*printf (" %11.6f", imgr1_ptr->channel_ptr[channel]->reflectance[idx]);*/
					/*}*/
					/*printf ("\t ]\n");*/
				/*}*/
				/*printf("\n\n");*/

			}else{
				imgr1_ptr->channel_ptr[channel]->reflectance = 0;
			}

			///////////////////////////////////////////////////////////////////////////////////
			//  Try getting the brightness temperature, and the number of rows, columns and  //
			//  number of points for this dataset.                                           //
			///////////////////////////////////////////////////////////////////////////////////

			sprintf(dsetName,"%s%s",grpName,"/BrightnessTemperature");
			printf("\tReading in channel %d brightnessTemp...\n",channel);

			h5status = H5LTfind_dataset ( group_ID, "BrightnessTemperature" );
			/*printf("\tThe return value for finding BrightnessTemperature is %d\n",h5status);*/
			if (h5status){
				h5status = viirs_level1_read (dsetName, fileID, 
					(void **) &bufferPtr_void,
					&nRank,
					&imgr1_ptr->channel_ptr[channel]->nrows,
					&imgr1_ptr->channel_ptr[channel]->ncols
				);
				/*printf("\tThe return value for reading the BrightnessTemperature is %d\n",h5status);*/

				// If the brightness temperature scale factors exist, get the factors
				scaleDataSet = 0;
				h5status = H5LTfind_dataset ( group_ID, "BrightnessTemperatureFactors" );
				if (h5status){
					scaleDataSet = 1;
					sprintf(dsetName,"%s%s",grpName,"/BrightnessTemperatureFactors");
					dsetID = H5Dopen2 (fileID, dsetName, H5P_DEFAULT);
					h5status = H5Dread (dsetID,              // dataset ID
						H5T_NATIVE_FLOAT,    // ident of memory datatype
						H5S_ALL,             // ident of memory dataspace
						H5S_ALL,             // ident of dataset's dataspace in file
						H5P_DEFAULT,         // File access properties
						scaleFactors);        // Pointer to buffer receiving data

						/*printf("\tBrightnessTemperature Scale Factor = %f\n",scaleFactors[0]);*/
						/*printf("\tBrightnessTemperature Offset = %f\n",scaleFactors[1]);*/
				}

				// Allocate memory for the brightnessTemp if we need to unscale, or point the brightnessTemp
				// pointer to the void pointer.
				bufferPtr_ushort = 0;
				if (scaleDataSet){
					if ((imgr1_ptr->channel_ptr[channel]->brightnessTemp = (float *) malloc(npts_channel*sizeof(float))) == NULL)
						error_allo(rout,"brightnessTemp allocation failed");
					bufferPtr_ushort = (unsigned short *) bufferPtr_void;
					bufferPtr_void = 0;
				}else{
					imgr1_ptr->channel_ptr[channel]->brightnessTemp = (float *) bufferPtr_void;
					bufferPtr_void = 0;
				}
				
				// Print out the initial read of the brightnessTemp dataset.
				/*printf("\n\tOriginal BrightnessTemperature =\n");*/
				/*for (i=0; i<10; i++) {*/
					/*printf ("\t [ ");*/
					/*for (j=0; j<10; j++){*/
						/*idx = i * imgr1_ptr->channel_ptr[channel]->ncols + j ;*/
						/*if (scaleDataSet){*/
							/*printf (" %u", bufferPtr_ushort[idx]);*/
						/*}else{*/
							/*printf (" %11.6f", imgr1_ptr->channel_ptr[channel]->brightnessTemp[idx]);*/
						/*}*/
					/*}*/
					/*printf ("\t ]\n");*/
				/*}*/
				/*printf("\n\n");*/

				// Scale the brightnessTemp, and free the pointer to the scaled BrightnessTemperature, if we need to
				if (scaleDataSet){
					for (idx = 0; idx < npts_channel; idx++){
						if (bufferPtr_ushort[idx] > UINT16_FILL_TEST){
							imgr1_ptr->channel_ptr[channel]->brightnessTemp[idx] = FLOAT32_FILL_TEST;
						}else{
							imgr1_ptr->channel_ptr[channel]->brightnessTemp[idx] = 
								(float) bufferPtr_ushort[idx] * scaleFactors[0] + scaleFactors[1];
						}
					}
					free(bufferPtr_ushort);
				}

				// Print out the final brightnessTemp dataset.
				/*printf("\n\tFinal BrightnessTemperature =\n");*/
				/*for (i=0; i<10; i++) {*/
					/*printf ("\t [ ");*/
					/*for (j=0; j<10; j++){*/
						/*idx = i * imgr1_ptr->channel_ptr[channel]->ncols + j ;*/
						/*printf (" %11.6f", imgr1_ptr->channel_ptr[channel]->brightnessTemp[idx]);*/
					/*}*/
					/*printf ("\t ]\n");*/
				/*}*/
				/*printf("\n\n");*/

			}else{
				imgr1_ptr->channel_ptr[channel]->brightnessTemp = 0;
			}

			// Close and release file.
			h5status = H5Fclose (fileID);
			if (h5status == 0) printf("\tSuccessfully closed file.\n");

		}else{
			// We don't need this channel, skip and set channel pointer to NULL
			/*printf("\tis to be skipped.\n");*/
			imgr1_ptr->channel_ptr[channel] = 0;
			continue;
		}
	}
	printf("\n\n");		

	return ;

}

/****************************************************************************************
 ****************************************************************************************
 *	      Cleanup of VIIRS Radiances, Reflectances, and Brightness Temperatures         *
 *        and geolocation                                                               *
 ****************************************************************************************
 ****************************************************************************************/

void deallocate_viirs(int8 *dummy, imagerL1 *imgr1_ptr) {

	int channel;

	// Free quantities that are always allocated memory

	free(imgr1_ptr->ModeScan);
	free(imgr1_ptr->lat);
	free(imgr1_ptr->lon);
	free(imgr1_ptr->solzen);
	free(imgr1_ptr->solaz);
	free(imgr1_ptr->satzen);
	free(imgr1_ptr->cos_satzen);
	free(imgr1_ptr->sataz);
	free(imgr1_ptr->zsfc);

	// Determine which spectral channels were allocated memory,
	// free the associated arrays, then free the channel structure
	// pointer.

	printf("Freeing VIIRS channel arrays...");

	for (channel=0;channel<NVIIRS_CHAN;channel++){

		/*printf("Channel %3d : ",channel);*/
		if (imgr1_ptr->channel_ptr[channel] != NULL){

			/*printf(" is valid, freeing array pointers...\n");*/
			if (imgr1_ptr->channel_ptr[channel]->radiance != NULL){
				/*printf("\tFreeing the radiance pointer \n");*/
				free(imgr1_ptr->channel_ptr[channel]->radiance);
			}
			if (imgr1_ptr->channel_ptr[channel]->reflectance != NULL){
				/*printf("\tFreeing the reflectance pointer \n");*/
				free(imgr1_ptr->channel_ptr[channel]->reflectance);
			}
			if (imgr1_ptr->channel_ptr[channel]->brightnessTemp != NULL){
				/*printf("\tFreeing the brightness temperature pointer \n");*/
				free(imgr1_ptr->channel_ptr[channel]->brightnessTemp);
			}

			/*printf("\tfreeing the channel structure pointer\n");	*/
			free(imgr1_ptr->channel_ptr[channel]);
		}else{
			/*printf(" is NULL, ignoring\n");	*/
		}
	}
	free(imgr1_ptr->channel_ptr);
	printf(" done\n");

}
