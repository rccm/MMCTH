/*$Id: modis_level1b_read.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

/******************************************************************************/
/* Read in MODIS Level 1B 1 km Data from HDF Scientific Data Sets.  Every     */
/* every line and pixel is read into a pointer for a given band.              */
/*                                                                            */
/* Input:                                                                     */
/*        filename: the name of a DAAC or direct broadcast MODIS level 1b     */
/*                  HDF file                                                  */				                     
/*          option:                                                           */
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
/*                                                                            */
/*            band: the MODIS band number of interest (1-36) use 13.0 for     */
/*                  13lo and 13.5 for 13hi and 14.0 for 14lo and 14.5 for     */
/*                  14hi                                                      */
/* Author:                                                                    */
/*         Michael Pavolonis                                                  */
/*         CIMSS/SSEC                                                         */
/******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"

float planck_btemp_terra (int, float);
float planck_btemp_aqua (int, float);

/*----------------------------------------------------------------------------*/

void modis_level1b_read (char *filename, int satname, int id, int option, float band, int SUB_FLG,
                         void **buf, unsigned char *badmask, int *nlines, int *npixels)
  
{

  void
    *attr;
    
  char8
    sds_name[MAX_NC_NAME],
    sds_name_check[MAX_NC_NAME],
    att_name[MAX_NC_NAME],
    scale_string[MAX_NC_NAME],
    offset_string[MAX_NC_NAME],
    *rout = {"modis_level1b_read"};

  
//  uint8
//    *buf_ui8;
  
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
    *buf_f32,
//    wl,
//    wl_meters,
//    rad,
    range[2];
    
  float64
//    *buf_f64,
    *scale_fac_f64;
  
  int32
    status,
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
  
   scale_fac=NULL;
   offset=NULL;
   buf_f32=NULL;
  
  pix_offset = 0;
  if (SUB_FLG == 2) pix_offset = 1;
  
  /*--------------------------------------------------------------------------*/
  /* Open the HDF file.                                                       */
  /*--------------------------------------------------------------------------*/
  
  if (strcmp(filename,"OPEN_250") != 0 && strcmp(filename,"OPEN_500") != 0 && strcmp(filename,"OPEN_1000") != 0) {
    id = SDstart(filename, DFACC_READ);
    if (id == FAIL) {
      fprintf(stderr,"%sInvalid HDF file, %s\n",EXE_PROMPT,filename);
      exit(EXIT_FAILURE);
    }
  }
 
  /*--------------------------------------------------------------------------*/
  /* Check for improper input.						      */
  /*--------------------------------------------------------------------------*/
  
  if (option == 4 && (band >= 20. && band != 26.))
    error_exit(rout,"Reflectance units are only valid for bands 1-19 and 26 only.");
  
  if (option == 5 && (band <= 19. || band == 26.))
    error_exit(rout,"Temperature units are only valid for bands 20-25 and 27-36 only.");
  
  /*--------------------------------------------------------------------------*/
  /* Read data for a given MODIS spectral band.                               */
  /*--------------------------------------------------------------------------*/
  
  if (option >= 1 && option <= 5) {
  
    /*--------------------------------------------------------------------------*/
    /* Determine which SDS the MODIS band of interest resides in.               */
    /*--------------------------------------------------------------------------*/
  
    if (band >= 1. && band <= 2.) {
      if ((strstr(filename,"02QKM.A") != NULL) || (strstr(filename,".250m.hdf") != NULL) || (strstr(filename,"OPEN_250") != NULL))
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
      fprintf(stderr,"%s%s-The band number specified (%f) is not available-exiting.\n",EXE_PROMPT,rout,band);
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
        fprintf(stderr,"%s%s-invalid option input, option=%d\n",EXE_PROMPT,rout,option);
        exit(EXIT_FAILURE);
    }
    
    /*------------------------------------------------------------------------*/
    /* Attach to the desired scientific data set in the HDF file.             */
    /*------------------------------------------------------------------------*/
    
    sd_id = SDselect(id,sds_index);
    if (sd_id == FAIL) {
      fprintf(stderr,"%s%s-cannot select specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
  
    /*--------------------------------------------------------------------------*/
    /* Read in the valid_range attribute.                                       */
    /*--------------------------------------------------------------------------*/
  
    strcpy(att_name,"valid_range");
    att_index = SDfindattr(sd_id,att_name);
    status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDattrinfo - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }
    
    if ((attr = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
      error_allo(rout,"attr-valid_range");
    status = SDreadattr(sd_id,att_index,attr);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDreadattr - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }
    
    valid_range = (uint16 *) attr;
  
    /*--------------------------------------------------------------------------*/
    /* Read in the fill_value attribute.					*/
    /*--------------------------------------------------------------------------*/
  
    strcpy(att_name,"_FillValue");
    att_index = SDfindattr(sd_id,att_name);
    status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDattrinfo - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }
      
    if ((attr = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
      error_allo(rout,"attr-_FillValue");
    status = SDreadattr(sd_id,att_index,attr);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDreadattr - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
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
      if (status == FAIL) {
        fprintf(stderr,"%s%s-There was an error when calling SDattrinfo - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
        exit(EXIT_FAILURE);
      }

      if ((attr = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
        error_allo(rout,"attr-scale_factor");
      status = SDreadattr(sd_id,att_index,attr);
      if (status == FAIL) {
        fprintf(stderr,"%s%s-There was an error when calling SDreadattr - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
        exit(EXIT_FAILURE);
      }
	
      scale_fac = (float32 *) attr;
  
      strcpy(att_name,offset_string);
      att_index = SDfindattr(sd_id,att_name);
      status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
      if (status == FAIL) {
        fprintf(stderr,"%s%s-There was an error when calling SDattrinfo - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
        exit(EXIT_FAILURE);
      }
        
      if ((attr = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
        error_allo(rout,"attr-offset");
      status = SDreadattr(sd_id,att_index,attr);
      if (status == FAIL) {
        fprintf(stderr,"%s%s-There was an error when calling SDreadattr - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
        exit(EXIT_FAILURE);
      }
        
      offset = (float32 *) attr;
    }

    /*--------------------------------------------------------------------------*/
    /* Read in information about the dimensions of the data of interest.        */
    /*--------------------------------------------------------------------------*/
  
    status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot get info about specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
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
  
    if ((buf_u16 = (uint16 *) malloc(npts*sizeof(uint16))) == NULL)
      error_allo(rout,"buf_u16");
    if (option == 1) {
      *buf = buf_u16;
    }
    else {
      if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
        error_allo(rout,"buf_f32");
      *buf = buf_f32;
    }
  
    /*------------------------------------------------------------------------*/
    /* Read in the data to a buffer.		                              */
    /*------------------------------------------------------------------------*/
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_u16);
  
    /*------------------------------------------------------------------------*/
    /* Print the status of the read: status 0 if successful, -1 if failed.    */
    /*------------------------------------------------------------------------*/

    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot read specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }

    /*------------------------------------------------------------------------*/
    /* Detach from the scientific data set.                                   */
    /*------------------------------------------------------------------------*/
    
    status = SDendaccess(sd_id);

    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot end access to specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
  
    /*--------------------------------------------------------------------------*/
    /* Scale the MODIS data, if desired and convert to brightnes temperature if */
    /* desired.								        */
    /*--------------------------------------------------------------------------*/
  
    if (option != 1) {
    
//    printf("scaling begin %d %f %f\n", fill_val[0], scale_fac[band_index], offset[band_index]);
//    (void)fflush(stdout);
      
      for (n=0; n<npts; n++) {
   
	if (buf_u16[n] == fill_val[0]) {
          buf_f32[n] = MISSING_FLOAT;

//        if (strcmp(filename,"OPEN_250") == 0) {
//          printf("setting badmask %d\n", n);
//          (void)fflush(stdout);
//        }
      
          badmask[n] = YES;
        }
        else {

//        if (strcmp(filename,"OPEN_250") == 0) {
//          printf("scaling data %f %f %d\n", scale_fac[band_index], offset[band_index], buf_u16[n]);
//          (void)fflush(stdout);
//        }
      
	  buf_f32[n] = scale_fac[band_index] * (buf_u16[n] - offset[band_index]);
    
          if (option == 5) {				       
	   if (satname == TERRA)
	     buf_f32[n] = planck_btemp_terra(band, buf_f32[n]);
	   else if (satname == AQUA)
	     buf_f32[n] = planck_btemp_aqua(band, buf_f32[n]);
	   else {
	     fprintf(stderr,"%s%s - Unknown platform name - aborting\n",EXE_PROMPT,rout);
	     exit(EXIT_FAILURE);
	   }
	   
          }  
        }
    
      }

      free(buf_u16);
      free(scale_fac);
      free(offset);
    }
    
    free(valid_range);
    free(fill_val);
    
  }
  
  /*--------------------------------------------------------------------------*/
  /* Read latitude or longitude.`	                                      */
  /*--------------------------------------------------------------------------*/
  
  if (option == 6 || option == 7) {
  
    switch (option) {
      case 6:
        strcpy(sds_name,"Latitude");
	range[0] = -90.0;
	range[1] = 90.0;
        break;
      case 7:
        strcpy(sds_name,"Longitude");
	range[0] = -180.0;
	range[1] = 180.0;  
	break;
    }
  
    sds_index = SDnametoindex(id,sds_name);
    sd_id = SDselect(id,sds_index);
    if (sd_id == FAIL) {
      fprintf(stderr,"%s%s-cannot select specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot get info on specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    start[0] = 0;
    start[1] = 0;
  
    edge[0] = dimen[0];
    edge[1] = dimen[1]-pix_offset;
    
    npts = edge[0]*edge[1];
    *nlines = dimen[0];
    *npixels = dimen[1]-pix_offset;
    
    if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
      error_allo(rout,"buf_f32");
    *buf = buf_f32;
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_f32);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot read specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }   
    
    status = SDendaccess(sd_id);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot end access to specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    for (n=0; n<npts; n++) {
      if (buf_f32[n] <  range[0] || buf_f32[n] >  range[1]) {
	buf_f32[n] = MISSING_FLOAT;
      }
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
    if (sd_id == FAIL) {
      fprintf(stderr,"%s%s-cannot select specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    strcpy(att_name,"scale_factor");
    att_index = SDfindattr(sd_id,att_name);
    status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDattrinfo - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }

    if ((attr = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
      error_allo(rout,"attr-scale_factor");
    status = SDreadattr(sd_id,att_index,attr);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDreadattr - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }
	
    scale_fac_f64 = (float64 *) attr;
    
    status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot get info on specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    start[0] = 0;
    start[1] = 0;
  
    edge[0] = dimen[0];
    edge[1] = dimen[1]-pix_offset;
    
    npts = edge[0]*edge[1];
    *nlines = dimen[0];
    *npixels = dimen[1]-pix_offset;
    
    if ((buf_i16 = (int16 *) malloc(npts*sizeof(int16))) == NULL)
      error_allo(rout,"buf_i16");
    if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
      error_allo(rout,"buf_f32");
    *buf = buf_f32;
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot read specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }  
    
    status = SDendaccess(sd_id);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot end access to specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    for (n=0; n<npts; n++) {
      buf_f32[n] = (float32) buf_i16[n] * scale_fac_f64[0];
    }
    
    free(buf_i16);
    free(scale_fac_f64);
    
  }
  
  /*--------------------------------------------------------------------------*/
  /* Read Solar and Sensor azimuth and compute relative azimuth.              */
  /*--------------------------------------------------------------------------*/
  
  if (option == 10) {
          
    sds_index = SDnametoindex(id,"SolarAzimuth");
    sd_id = SDselect(id,sds_index);
    if (sd_id == FAIL) {
      fprintf(stderr,"%s%s-cannot select specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    strcpy(att_name,"scale_factor");
    att_index = SDfindattr(sd_id,att_name);
    status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDattrinfo - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }

    if ((attr = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
      error_allo(rout,"attr-scale_factor");
    status = SDreadattr(sd_id,att_index,attr);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDreadattr - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }
	
    scale_fac_f64 = (float64 *) attr;
    
    status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot get info on specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    start[0] = 0;
    start[1] = 0;
  
    edge[0] = dimen[0];
    edge[1] = dimen[1]-pix_offset;
    
    npts = edge[0]*edge[1];
    *nlines = dimen[0];
    *npixels = dimen[1]-pix_offset;
    
    if ((buf_i16 = (int16 *) malloc(npts*sizeof(int16))) == NULL)
      error_allo(rout,"buf_i16");
    if ((buf_i16_2 = (int16 *) malloc(npts*sizeof(int16))) == NULL)
      error_allo(rout,"buf_i16_2");
    if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
      error_allo(rout,"buf_f32");
    *buf = buf_f32; 
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot read specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }  
    
    status = SDendaccess(sd_id);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot end access to specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    sds_index = SDnametoindex(id,"SensorAzimuth");
    sd_id = SDselect(id,sds_index);
    sd_id = SDselect(id,sds_index);
    if (sd_id == FAIL) {
      fprintf(stderr,"%s%s-cannot select specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16_2);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot read specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }   
    
    status = SDendaccess(sd_id);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot end access to specified SDS, option=%d\n",EXE_PROMPT,rout,option);
      exit(EXIT_FAILURE);
    }
    
    
    for (n=0; n<npts; n++) {
      /*buf_f32[n] = 180.0 - (((float32) buf_i16[n] * scale_fac_f64[0]) - ((float32) buf_i16_2[n] * scale_fac_f64[0]));*/
      buf_f32[n] = fabs(180.0 - fabs(((float32) buf_i16_2[n] * scale_fac_f64[0]) - ((float32) buf_i16[n] * scale_fac_f64[0])));

     }
    
    free(buf_i16);
    free(buf_i16_2);
    free(scale_fac_f64);
    
  }
  
  /*--------------------------------------------------------------------------*/
  /* Close the HDF file.                                                      */
  /*--------------------------------------------------------------------------*/
  
  if (strcmp(filename,"OPEN_250") != 0 && strcmp(filename,"OPEN_500") != 0 && strcmp(filename,"OPEN_1000") != 0) {
    status = SDend(id);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-cannot close %s\n",EXE_PROMPT,rout,filename);
      exit(EXIT_FAILURE);
    }
  }
  
  /*--------------------------------------------------------------------------*/
  /* DONE with the routine so exit.                                           */
  /*--------------------------------------------------------------------------*/

}

unsigned char * unpack_detector_flags(int flag_id, int ndet, unsigned char *packed_flags)
{

  char *rout = {"unpack_detector_flags"};
  unsigned char m_bits, *flag;
  int n_shift, i;
  
  m_bits = pow(2,flag_id-1);
  n_shift = flag_id - 1;
  
  /*printf("flag_id=%d, n_bits=%d, n_shift=%d\n",flag_id,m_bits,n_shift);*/
  
  if ((flag = (unsigned char *) malloc(ndet*sizeof(unsigned char))) == NULL)
    error_allo(rout,"flag");
    
  for (i=0; i<ndet; i++) {
    flag[i] = (packed_flags[i] & m_bits) >> n_shift;
  }
  
  return(flag);

}
  
/******************************************************************************/
/******************************************************************************/
