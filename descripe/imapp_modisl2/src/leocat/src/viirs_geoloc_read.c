/******************************************************************************/
/* This routine reads in VIIRS geolocation data from an HDF5 file             */
/*                                                                            */
/* Input:                                                                     */
/*         dsetName: HDF5 path of the desired dataset.                        */
/*         fileId  : HDF5 file handle.                                        */
/*         buf_ptr : pointer to the dataset array.                            */
/*         nrank   : the rank of the dataset (0==scalar,1==vector,2==matrix)  */
/*         nrows   : the number of rows of the dataset.                       */
/*         ncols   : the number of columns of the dataset.                    */
/*                                                                            */
/* Author:                                                                    */
/*         Geoff Cureton                                                      */
/*         CIMSS/SSEC                                                         */
/******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <common_leocat.h>

int viirs_geoloc_read (char *dsetName, hid_t fileId, void **buf_ptr, int *nrank, int *nrows, int *ncols)
{
	hid_t       dsetID, dspaceID, dTypeID;          /* Handles */
	herr_t      h5_dsetStatus,h5_dspaceStatus,h5_dsetReadStatus;
	hsize_t     dimSize[2];
	int rank, npts;
	int rows,cols;
	float *buf_float;
	unsigned short *buf_ushort;
	char *rout = {" viirs_geoloc_read"};

    /*************************
	 *  Read in the dataset  *
     *************************/

	// Get the dataset ID for the dataset
	if ((dsetID = H5Dopen2 (fileId, dsetName, H5P_DEFAULT)) < 0){
		printf("\n\t>>> Error opening VIIRS L1 dataset %s\n",dsetName);
		*buf_ptr = 0;
		return dsetID;
	}

	// Get the datatype ID for the desired dataset
	dTypeID = H5Dget_type(dsetID);
	dTypeID = H5Tget_native_type(dTypeID,H5T_DIR_ASCEND);

	// Get the datatype description.
	char *str;
	size_t len;
	herr_t status;
	status = H5LTdtype_to_text(dTypeID, NULL, H5LT_DDL, &len);
	str = (char *) malloc(len * sizeof(char));
	status = H5LTdtype_to_text(dTypeID, str, H5LT_DDL, &len);
	printf("\tString description for the type : %s\n",str);
	free(str);

    // Get dataspace ID  
	if ((dspaceID = H5Dget_space(dsetID)) < 0){
		printf("\n\t>>> Error getting dataspace for VIIRS L1 dataset %s\n",dsetName);
		h5_dsetStatus = H5Dclose (dsetID);
		*buf_ptr = 0;
		return dspaceID;
	}

	// Get the number of dimensions, dimension sizes and maximum dimension
	// sizes for the dataset (through the associated dataspace)
	if ((rank = H5Sget_simple_extent_dims(dspaceID, dimSize, NULL)) < 0){
		printf("\n\t>>> Error getting dataspace dimensions for VIIRS L1 dataset %s\n",dsetName);
		h5_dspaceStatus = H5Sclose (dspaceID);
		h5_dsetStatus = H5Dclose (dsetID);
		*buf_ptr = 0;
		return rank;
	}
	h5_dspaceStatus = H5Sclose (dspaceID);

	*nrank = rank;
    
	rows = dimSize[0];
	cols = dimSize[1];
	npts = rows * cols; 

	buf_ushort = 0;
	buf_float = 0;

	/*H5T_NATIVE_CHAR  	,     // char*/
	/*H5T_NATIVE_SCHAR 	,     // signed char*/
	/*H5T_NATIVE_UCHAR 	,     // unsigned char*/
	/*H5T_NATIVE_SHORT 	,     // short*/
	/*H5T_NATIVE_USHORT 	,     // unsigned short*/
	/*H5T_NATIVE_INT 		,     // int*/
	/*H5T_NATIVE_UINT 	,     // unsigned*/
	/*H5T_NATIVE_LONG 	,     // long*/
	/*H5T_NATIVE_ULONG 	,     // unsigned long*/
	/*H5T_NATIVE_LLONG 	,     // long long*/
	/*H5T_NATIVE_ULLONG 	,     // unsigned long long*/
	/*H5T_NATIVE_FLOAT 	,     // float*/
	/*H5T_NATIVE_DOUBLE 	,     // double*/
	/*H5T_NATIVE_LDOUBLE 	,     // long double*/
	/*H5T_NATIVE_HSIZE 	,     // hsize_t*/
	/*H5T_NATIVE_HSSIZE 	,     // hssize_t*/
	/*H5T_NATIVE_HERR 	,     // herr_t*/
	/*H5T_NATIVE_HBOOL 	      // hbool_t*/

	// Reserve memory for the array, returning a pointer to that array.
	/*if ((buf_float = (float *) malloc(npts*sizeof(float))) == NULL)*/
		/*error_allo(rout,"buf_float");*/
	/**buf_ptr = buf_float;*/
	if (H5Tequal(dTypeID,H5T_NATIVE_USHORT)){
		if ((buf_ushort = (unsigned short *) malloc(npts*sizeof(unsigned short))) == NULL)
			error_allo(rout,"buf_ushort");
		*buf_ptr = buf_ushort;
	}else if(H5Tequal(dTypeID,H5T_NATIVE_FLOAT)){
		if ((buf_float = (float *) malloc(npts*sizeof(float))) == NULL)
			error_allo(rout,"buf_float");
		*buf_ptr = buf_float;
	}else{
		printf("\n\t>>> Error allocating C-type from HDF5 type for VIIRS L1 dataset %s\n",dsetName);
		h5_dspaceStatus = H5Sclose (dspaceID);
		h5_dsetStatus = H5Dclose (dsetID);
		*buf_ptr = 0;
	}

	// Read the dataset into the array
	/*h5status = H5Dread (dsetID,            // dataset ID*/
					  /*dTypeID         ,    // ident of memory datatype*/
					  /*H5S_ALL,             // ident of memory dataspace*/
					  /*H5S_ALL,             // ident of dataset's dataspace in file*/
					  /*H5P_DEFAULT,         // File access properties*/
					  /*buf_float);          // Pointer to buffer receiving data*/
	h5_dsetReadStatus = 0;
	if (H5Tequal(dTypeID,H5T_NATIVE_USHORT)){
		h5_dsetReadStatus = H5Dread (dsetID,     // dataset ID
							dTypeID,             // ident of memory datatype
							H5S_ALL,             // ident of memory dataspace
							H5S_ALL,             // ident of dataset's dataspace in file
							H5P_DEFAULT,         // File access properties
							buf_ushort);           // Pointer to buffer receiving data
	}else if(H5Tequal(dTypeID,H5T_NATIVE_FLOAT)){
		h5_dsetReadStatus = H5Dread (dsetID,     // dataset ID
							dTypeID,             // ident of memory datatype
							H5S_ALL,             // ident of memory dataspace
							H5S_ALL,             // ident of dataset's dataspace in file
							H5P_DEFAULT,         // File access properties
							buf_float);           // Pointer to buffer receiving data
	}else{
		printf("\n\t>>> Datatype Error reading into C-type from HDF5 type for VIIRS L1 dataset %s\n",dsetName);
		h5_dspaceStatus = H5Sclose (dspaceID);
		h5_dsetStatus = H5Dclose (dsetID);
		*buf_ptr = 0;
	}

	if (h5_dsetReadStatus < 0){		// There was an error, cleaning up...
		printf("\n\t>>> Error reading VIIRS L1 dataset %s into buf_(ushort/float)\n",dsetName);
		h5_dsetStatus = H5Dclose (dsetID);
		if (buf_ushort != NULL) free(buf_ushort);
		if (buf_float != NULL) free(buf_float);
		*buf_ptr = 0;
		return h5_dsetReadStatus;
	}

	*nrows = rows;
	*ncols = cols;

	// Close and release dataset and dataspace
	/*h5status = H5Dclose (dsetID);*/
	h5_dsetStatus = H5Dclose (dsetID);
	/*h5status = H5Sclose (dspaceID);*/

	/*return;*/
	return h5_dsetStatus;
	

}
