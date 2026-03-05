/*----------------------------------------------------------------------
#    Copyright (C) 2011,  Space Science and Engineering Center,
#    University C  of Wisconsin-Madison, Madison WI.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------*/
/******************************************************************************/
/* Read in MODIS Level 1B 1 km Data from HDF Scientific Data Sets.  Every     */
/* every line and pixel is read into a pointer for a given band.              */
/*    						                              */
/* Input:								      */
/*        filename: the name of a DAAC or direct broadcast MODIS level 1b     */
/*                  HDF file                                                  */				                     
/*          option: 							      */
/*                  1: Output raw data straight from HDF file                 */
/*                  2: Output corrected counts                                */
/*                  3: Output radiance                                        */
/*                  4: Output reflectance (SZA uncorrected) (1-19, 16 only)   */
/*                  5: Output brightness temperature (K) (20-25, 27-36 only)  */
/*                  6: Output latitude (deg.)	                              */
/*                  7: Output longitude (deg.)                                */
/*                  8: Output solar zenith angle (deg.)                       */
/*                  9: Output satellite zenith angle (deg.)                   */
/*                 10: Output relative azimuth (deg.)                         */
/*									      */
/*            band: the MODIS band number of interest (1-36) use 13.0 for     */
/*                  13lo and 13.5 for 13hi and 14.0 for 14lo and 14.5 for     */
/*                  14hi						      */
/* Author:								      */
/*         Michael Pavolonis					              */
/*         CIMSS/SSEC			  			              */
/******************************************************************************/


#include <math.h>
#include <hdf.h>
#include "subpixel_cloudfrac.h"

float32 planck_btemp(float32, float32);

/*----------------------------------------------------------------------------*/

void modis_level1b_read (char *filename, int option, float band, void **buf,
                         int *nlines, int *npixels, float *earth2sun,
			 int SUB_FLG)
  
{

  void
    *attr;
    
  char8
    sds_name[MAX_NC_NAME],
    sds_name_check[MAX_NC_NAME],
    att_name[MAX_NC_NAME],
    scale_string[MAX_NC_NAME],
    offset_string[MAX_NC_NAME],
    *units;
  
  int16
    *buf_i16,
    *buf_i16_2,
    pix_offset;
  
  uint16
    *buf_u16,
    *valid_range,
    *fill_val;
    
  float32
    *scale_fac,
    *offset,
    *earth2sun_temp,
    *buf_f32,
    wl,
    wl_meters,
    rad;
    
  float64
    *buf_f64;
  
  int32
    status,
    id,
    sd_id,
    sds_index,
    att_index,
    band_index, 
    data_type,
    count,
    npts,
    n,
    rank,
    nattr,
    dimen[MAX_VAR_DIMS],
    start[3],
    edge[3];	
    
    const float32 cwn[] = {2.641775E+03, 2.505277E+03, 2.518028E+03, 2.465428E+03,
                           2.235815E+03, 2.200346E+03, 0.0000000000, 1.477967E+03,
		           1.362737E+03, 1.173190E+03, 1.027715E+03, 9.080884E+02,
		           8.315399E+02, 7.483394E+02, 7.308963E+02, 7.188681E+02,
		           7.045367E+02};
  
    const float32 tcs[] = {9.993411E-01, 9.998646E-01, 9.998584E-01, 9.998682E-01,
                           9.998819E-01, 9.998845E-01, 0.0000000000, 9.994877E-01, 
		           9.994918E-01, 9.995495E-01, 9.997398E-01, 9.995608E-01,
		           9.997256E-01, 9.999160E-01, 9.999167E-01, 9.999191E-01,
		           9.999281E-01};

    const float32 tci[] = {4.770532E-01, 9.262664E-02, 9.757996E-02, 8.929242E-02,
                           7.310901E-02, 7.060415E-02, 0.0000000000, 2.204921E-01,
		           2.046087E-01, 1.599191E-01, 8.253401E-02, 1.302699E-01,
		           7.181833E-02, 1.972608E-02, 1.913568E-02, 1.817817E-02,
		           1.583042E-02};
  
  pix_offset = 0;
  if (SUB_FLG == 2) pix_offset = 1;
  
  /*--------------------------------------------------------------------------*/
  /* Open the HDF file.                                                       */
  /*--------------------------------------------------------------------------*/
  
  id = SDstart(filename, DFACC_READ);
  if (id == -1) {
    fprintf(stderr,"ERROR: Invalid HDF file: %s\n", filename);
    exit(EXIT_FAILURE);
  }
 
  /*--------------------------------------------------------------------------*/
  /* Check for improper input.						      */
  /*--------------------------------------------------------------------------*/
  
  if (option == 4 && (band >= 20. && band != 26.))
    error_exit("Reflectance units are only valid for bands 1-19 and 26 only-exiting.");
  
  if (option == 5 && (band <= 19. || band == 26.))
    error_exit("Temperature units are only valid for bands 20-25 and 27-36 only-exiting.");
  
  /*--------------------------------------------------------------------------*/
  /* Read data for a given MODIS spectral band.                               */
  /*--------------------------------------------------------------------------*/
  
  if (option >= 1 && option <= 5) {
  
    /*--------------------------------------------------------------------------*/
    /* Determine which SDS the MODIS band of interest resides in.               */
    /*--------------------------------------------------------------------------*/
  
    if (band >= 1. && band <= 2.) {
      if ((strstr(filename,"02QKM.A") != NULL) || (strstr(filename,".250m.hdf") != NULL))
        strcpy(sds_name,"EV_250_RefSB");
      else
        strcpy(sds_name,"EV_250_Aggr1km_RefSB");
      band_index = (int32) band - 1;  
    }
    else if (band >= 3. && band <= 7.) {
      strcpy(sds_name,"EV_500_Aggr1km_RefSB");
      band_index = (int32) band - 3;
    }
    else if (band >= 8. && band <= 12.) {
      strcpy(sds_name,"EV_1KM_RefSB");
      band_index = (int32) band - 8;
    }
    else if (band >= 13. && band < 15.) {
      strcpy(sds_name,"EV_1KM_RefSB");
      band_index = (int32) 2.*band - 21.;
    }
    else if (band >= 15. && band <= 19.) {
      strcpy(sds_name,"EV_1KM_RefSB");
      band_index = (int32) band - 6.;
    }
    else if (band == 26.) {
      if (!SUB_FLG) {
        strcpy(sds_name,"EV_Band26");
        band_index = 0;
      }
      else {
        strcpy(sds_name,"EV_1KM_RefSB");
        band_index = 14;
      }
    }
    else if (band >= 20. && band <= 25.) {
      strcpy(sds_name,"EV_1KM_Emissive");
      band_index = (int32) band - 20;
    }
    else if (band >= 27. && band <= 36.) {
      strcpy(sds_name,"EV_1KM_Emissive");
      band_index = (int32) band - 21;
    }
    else {
      fprintf(stderr,"The band number specified (%d) is not available-exiting.\n",band);
      exit(EXIT_FAILURE);
    }
  
    /*--------------------------------------------------------------------------*/
    /* Convert the name of the SDS of interest to an index number.              */
    /*--------------------------------------------------------------------------*/
  
    sds_index = SDnametoindex(id,sds_name);
  
    /*--------------------------------------------------------------------------*/
    /* Decide which units for the MODIS data are desired.			*/
    /*--------------------------------------------------------------------------*/
  
    switch (option) {
      case 1:
        /*printf("\n\tRetrieving raw unscaled data for band %4.1f\n.",band);*/
        break;
      case 2:
        /*printf("\n\tRetrieving corrected counts for band %4.1f\n.",band);*/
        strcpy(scale_string,"corrected_counts_scales");
        strcpy(offset_string,"corrected_counts_offsets");
        break;
      case 3:
        /*printf("\n\tRetrieving radiance data for band %4.1f\n.",band);*/
        strcpy(scale_string,"radiance_scales");
        strcpy(offset_string,"radiance_offsets");
        break;
      case 4:
        /*printf("\n\tRetrieving reflectance data for band %4.1f\n.",band);*/
        strcpy(scale_string,"reflectance_scales");
        strcpy(offset_string,"reflectance_offsets");
        break;
      case 5:
        /*printf("\n\tRetrieving brightness temperature data for band %4.1f\n.",band);*/
        strcpy(scale_string,"radiance_scales");
        strcpy(offset_string,"radiance_offsets");
        break;
      default:
        error_exit("Option input is not valid-exiting");
    }
    
    /*--------------------------------------------------------------------------*/
    /* Read in the Earth-sun-distance attribute.                                */
    /*--------------------------------------------------------------------------*/
    
    strcpy(att_name,"Earth-Sun Distance");
    att_index = SDfindattr(id,att_name);
    status = SDattrinfo(id,att_index,att_name,&data_type,&count);
    if (status == -1)
      error_exit("There was an error when reading attribute Earth-Sun Distance, exiting.");
   
    attr = (void *) malloc(count*DFKNTsize(data_type));
    status = SDreadattr(id,att_index,attr);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDreadattr-exiting.\n");
      exit (1);
    }
    earth2sun_temp = (float32 *) attr;
    *earth2sun = earth2sun_temp[0];
    
    /*------------------------------------------------------------------------*/
    /* Attach to the desired scientific data set in the HDF file.             */
    /*------------------------------------------------------------------------*/
    
    sd_id = SDselect(id,sds_index);
  
    /*--------------------------------------------------------------------------*/
    /* Read in the valid_range attribute.                                       */
    /*--------------------------------------------------------------------------*/
  
    strcpy(att_name,"valid_range");
    att_index = SDfindattr(sd_id,att_name);
    status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDfindattr-exiting.\n");
      exit (1);
    }
    attr = (void *) malloc(count*DFKNTsize(data_type));
    status = SDreadattr(sd_id,att_index,attr);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDreadattr-exiting.\n");
      exit (1);
    }
    valid_range = (uint16 *) attr;
  
    /*--------------------------------------------------------------------------*/
    /* Read in the fill_value attribute.					*/
    /*--------------------------------------------------------------------------*/
  
    strcpy(att_name,"_FillValue");
    att_index = SDfindattr(sd_id,att_name);
    status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDfindattr-exiting.\n");
      exit (1);
    }
    attr = (void *) malloc(count*DFKNTsize(data_type));
    status = SDreadattr(sd_id,att_index,attr);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDreadattr-exiting.\n");
      exit (1);
    }
    fill_val = (uint16 *) attr;
  
    /*--------------------------------------------------------------------------*/
    /* Unless raw data is desired, read in the appropriate data scale factor    */ 
    /* and offset.                                                              */
    /*--------------------------------------------------------------------------*/
  
    if (option != 1) {
  
      strcpy(att_name,scale_string);
      att_index = SDfindattr(sd_id,att_name);
      status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
      if (status == -1) {
        printf("\n\tThere was an error when calling SDfindattr-exiting.\n");
        exit (1);
      }
      attr = (void *) malloc(count*DFKNTsize(data_type));
      status = SDreadattr(sd_id,att_index,attr);
      if (status == -1) {
        printf("\n\tThere was an error when calling SDreadattr-exiting.\n");
        exit (1);
      }
      scale_fac = (float32 *) attr;
  
      strcpy(att_name,offset_string);
      att_index = SDfindattr(sd_id,att_name);
      status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
      if (status == -1) {
        printf("\n\tThere was an error when calling SDfindattr-exiting.\n");
        exit (1);
      }
      attr = (void *) malloc(count*DFKNTsize(data_type));
      status = SDreadattr(sd_id,att_index,attr);
      if (status == -1) {
        printf("\n\tThere was an error when calling SDreadattr-exiting.\n");
        exit (1);
      }
      offset = (float32 *) attr;
    }

    /*--------------------------------------------------------------------------*/
    /* Read in information about the dimensions of the data of interest.        */
    /*--------------------------------------------------------------------------*/
  
    status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDgetinfo-exiting.\n");
      exit (1);
    }
  
    /*--------------------------------------------------------------------------*/
    /* Set arrays that are needed when reading in the actual data with          */
    /* dimension information.						        */
    /*--------------------------------------------------------------------------*/
    
    if (band == 26. && SUB_FLG == 0) {
    
      start[0] = 0;
      start[1] = 0;
  
      edge[0] = dimen[0];
      edge[1] = dimen[1]-pix_offset;
    
      npts = edge[0]*edge[1];
      *nlines = dimen[0];
      *npixels = dimen[1]-pix_offset;
    } 
    else {
    
      start[0] = band_index;
      start[1] = 0;
      start[2] = 0;
  
      edge[0] = 1;
      edge[1] = dimen[1];
      edge[2] = dimen[2]-pix_offset;
    
      npts = edge[1]*edge[2];
      *nlines = dimen[1];
      *npixels = dimen[2]-pix_offset;
    }
  
    /*------------------------------------------------------------------------*/
    /* Allocate memory for the data to be read in from the file and the       */ 
    /* unscaled data.						              */
    /*------------------------------------------------------------------------*/
  
    buf_u16 = (uint16 *) malloc(npts*sizeof(uint16));
    if (option == 1) {
      *buf = buf_u16;
    }
    else {
      buf_f32 = (float32 *) malloc(npts*sizeof(float32));
      *buf = buf_f32;
    }
  
    /*------------------------------------------------------------------------*/
    /* Read in the data to a buffer.		                              */
    /*------------------------------------------------------------------------*/
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_u16);
  
    /*------------------------------------------------------------------------*/
    /* Print the status of the read: status 0 if successful, -1 if failed.    */
    /*------------------------------------------------------------------------*/

    if (status != 0) {
      printf("\n\tThere was an error in reading the HDF data-exiting.\n");
      exit (1);
    }

    /*------------------------------------------------------------------------*/
    /* Detach from the scientific data set.                                   */
    /*------------------------------------------------------------------------*/
    
    status = SDendaccess(sd_id);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDendaccess-exiting.\n");
      exit (1);
    }
  
    /*--------------------------------------------------------------------------*/
    /* Scale the MODIS data, if desired and convert to brightnes temperature if */
    /* desired.								        */
    /*--------------------------------------------------------------------------*/
  
    if (option != 1) {
    
      for (n=0; n<npts; n++) {
   
        if (buf_u16[n] == fill_val[0]) {
          buf_f32[n] = -999.0;
        }
        else {
          buf_f32[n] = scale_fac[band_index] * (buf_u16[n] - offset[band_index]);
    
          if (option == 5) {
            /*buf_f32[n] = (planck_btemp(1.0e+4/cwn[(int32)band-20],buf_f32[n]) -
                                       tci[(int32)band-20])/tcs[(int32)band-20];*/
				       
	   wl = 1.0e+4/cwn[(int32)band-20];
	   wl_meters = 1.0e-6 * wl;
           buf_f32[n] = ((c2_planck / (wl_meters * log(c1_planck / 
	               (1.0e+6 * buf_f32[n] * pow(wl_meters,5)) + 1.0))) - 
		       tci[(int32)band-20])/tcs[(int32)band-20];
          }  
        }
    
      }
      
      free(buf_u16);
    }
    
    free(valid_range);
    free(fill_val);
    free(scale_fac);
    free(offset);
    
  }
  
  /*--------------------------------------------------------------------------*/
  /* Read latitude or longitude.`	                                      */
  /*--------------------------------------------------------------------------*/
  
  if (option == 6 || option == 7) {
  
    switch (option) {
      case 6:
        strcpy(sds_name,"Latitude");
        break;
      case 7:
        strcpy(sds_name,"Longitude");  
	break;
    }
  
    sds_index = SDnametoindex(id,sds_name);
    sd_id = SDselect(id,sds_index);
    status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDgetinfo-exiting.\n");
      exit (1);
    }
    
    start[0] = 0;
    start[1] = 0;
  
    edge[0] = dimen[0];
    edge[1] = dimen[1]-pix_offset;
    
    npts = edge[0]*edge[1];
    *nlines = dimen[0];
    *npixels = dimen[1]-pix_offset;
    
    buf_f32 = (float32 *) malloc(npts*sizeof(float32));
    *buf = buf_f32;
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_f32);
    
    if (status != 0) {
      printf("\n\tThere was an error in reading the HDF data-exiting.\n");
      exit (1);
    }   
    
    status = SDendaccess(sd_id);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDendaccess-exiting.\n");
      exit (1);
    }
    
  }
  
  /*--------------------------------------------------------------------------*/
  /* Read SZA or VZA.	                                                      */
  /*--------------------------------------------------------------------------*/
  
  if (option == 8 || option == 9) {
  
    switch (option) {
      case 8:
        strcpy(sds_name,"SolarZenith"); 
	break;
      case 9:
        strcpy(sds_name,"SensorZenith");
	break;  
    }
          
    sds_index = SDnametoindex(id,sds_name);
    sd_id = SDselect(id,sds_index);
    status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDgetinfo-exiting.\n");
      exit (1);
    }
    
    start[0] = 0;
    start[1] = 0;
  
    edge[0] = dimen[0];
    edge[1] = dimen[1]-pix_offset;
    
    npts = edge[0]*edge[1];
    *nlines = dimen[0];
    *npixels = dimen[1]-pix_offset;
    
    buf_i16 = (int16 *) malloc(npts*sizeof(int16));
    buf_f32 = (float32 *) malloc(npts*sizeof(float32));
    *buf = buf_f32;
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
    
    if (status != 0) {
      printf("\n\tThere was an error in reading the HDF data-exiting.\n");
      exit (1);
    }   
    
    status = SDendaccess(sd_id);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDendaccess-exiting.\n");
      exit (1);
    }
    
    for (n=0; n<npts; n++) {
      buf_f32[n] = (float32) buf_i16[n] * 0.01;
    }
    
    free(buf_i16);
    
  }
  
  /*--------------------------------------------------------------------------*/
  /* Read Solar and Sensor azimuth and compute relative azimuth.              */
  /*--------------------------------------------------------------------------*/
  
  if (option == 10) {
          
    sds_index = SDnametoindex(id,"SolarAzimuth");
    sd_id = SDselect(id,sds_index);
    status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDgetinfo-exiting.\n");
      exit (1);
    }
    
    start[0] = 0;
    start[1] = 0;
  
    edge[0] = dimen[0];
    edge[1] = dimen[1]-pix_offset;
    
    npts = edge[0]*edge[1];
    *nlines = dimen[0];
    *npixels = dimen[1]-pix_offset;
    
    buf_i16 = (int16 *) malloc(npts*sizeof(int16));
    buf_i16_2 = (int16 *) malloc(npts*sizeof(int16));
    buf_f32 = (float32 *) malloc(npts*sizeof(float32));
    *buf = buf_f32; 
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
    
    if (status != 0) {
      printf("\n\tThere was an error in reading the HDF data-exiting.\n");
      exit (1);
    }   
    
    status = SDendaccess(sd_id);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDendaccess-exiting.\n");
      exit (1);
    }
    
    sds_index = SDnametoindex(id,"SensorAzimuth");
    sd_id = SDselect(id,sds_index);
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16_2);
    
    if (status != 0) {
      printf("\n\tThere was an error in reading the HDF data-exiting.\n");
      exit (1);
    }   
    
    status = SDendaccess(sd_id);
    if (status == -1) {
      printf("\n\tThere was an error when calling SDendaccess-exiting.\n");
      exit (1);
    }
    
    /* Relative azimuth angle, abs(solar_az - sensor_az), so that 0
       is looking into the sun and 180 is looking away.  Range: 0..180. */
    
    for (n=0; n<npts; n++) {
//      buf_f32[n] = 180.0 - (((float32) buf_i16[n] * 0.01) - ((float32) buf_i16_2[n] * 0.01));
      buf_f32[n] = fabsf(((float32) buf_i16[n] * 0.01) - ((float32) buf_i16_2[n] * 0.01));
      if (buf_f32[n] > 180.) buf_f32[n] = 360. - buf_f32[n];
    }
    
    free(buf_i16);
    free(buf_i16_2);
    
  }    
  
  /*--------------------------------------------------------------------------*/
  /* Close the HDF file.                                                      */
  /*--------------------------------------------------------------------------*/
  
  status = SDend(id);
  if (status == -1) {
    printf("\n\tThere was an error when calling SDend-exiting.\n");
    exit (1);
  }
  
  /*--------------------------------------------------------------------------*/
  /* DONE with the routine so exit.                                           */
  /*--------------------------------------------------------------------------*/

}  

/******************************************************************************/


/******************************************************************************/
/* This function is used to compute brightness temperature in Kelvin.  Input  */
/* is as follows: wavelength (wl) in um and radiance in W*m-2*str-1*um-1      */
/******************************************************************************/
 
float32 planck_btemp (float32 wl,float32 rad)

{

  float32
    planckT,
    wl_meters;
    
  wl_meters = 1.0e-6 * wl;
  
  planckT = c2_planck / (wl_meters * log(c1_planck / (1.0e+6 * rad * pow(wl_meters,5)) + 1.0));
  
  return planckT;
  
}
  
/******************************************************************************/
/******************************************************************************/
