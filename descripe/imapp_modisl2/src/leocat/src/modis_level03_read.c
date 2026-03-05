/*$Id: modis_level03_read.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

/******************************************************************************/
/* This routine reads in 1 km MODIS geolocation data from an HDF.             */
/*                                                                            */
/* Input:                                                                     */
/*         filename: the name of the L03 geolocation file.                    */
/*                                                                            */
/* Author:                                                                    */
/*         Michael Pavolonis                                                  */
/*         CIMSS/SSEC                                                         */
/******************************************************************************/


#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"

void modis_level03_read (char *filename, int id, int option, int SUB_FLG, void **buf,
                         unsigned char *badmask, int *nlines, int *npixels)
			 
{

  void
    *attr,
    *attr1,
    *attr2;
    
  char8
    sds_name[MAX_NC_NAME],
    sds_name_check[MAX_NC_NAME],
    att_name[MAX_NC_NAME],
  

   
    *rout = {"modis_level03_read"};
  
  unsigned char
    *bufc8;
  
  int16
    *buf_i16,
    pix_offset,
    *valid_range_i16,
    *fill_val_i16;
    
  float32
    *buf_f32,
    *buf_f32_2,
    *valid_range_f32,
    *fill_val_f32;
    
  float64
    
    *scale_fac_f64;
  
  int32
    status,
    sd_id,
    sds_index,
    att_index,

    data_type,
    count,
    npts,
    n,
    rank,
    nattr,
    dimen[MAX_VAR_DIMS],
    start[2],
    edge[2];
    
  pix_offset = 0;
  if (SUB_FLG == 2) pix_offset = 1;
  
  /*--------------------------------------------------------------------------*/
  /* Open the HDF file.                                                       */
  /*--------------------------------------------------------------------------*/
  
  if (strcmp(filename,"OPEN_1000") != 0) {
    id = SDstart(filename, DFACC_READ);
    if (id == FAIL) {
      fprintf(stderr,"%sInvalid HDF file, %s\n",EXE_PROMPT,filename);
      exit (EXIT_FAILURE);
    }
  }
  
  /*--------------------------------------------------------------------------*/
  /* Check for improper input.                                                */
  /*--------------------------------------------------------------------------*/
  
  if (option < 1 || option > 8) {
    fprintf(stderr,"%smodis_level03_read - Invalid option, %d\n",EXE_PROMPT,option);
    exit (EXIT_FAILURE);
  }
  
  switch (option) {
    case 1:
      strcpy(sds_name,"Latitude");
      break;
    case 2:
      strcpy(sds_name,"Longitude");  
      break;
    case 3:
      strcpy(sds_name,"SolarZenith");
      break;
    case 4:
      strcpy(sds_name,"SensorZenith");
      break;
    case 5:
      strcpy(sds_name,"SolarAzimuth");
      break;
    case 6:
      strcpy(sds_name,"Land/SeaMask");
      break;
    case 7:
      strcpy(sds_name,"SolarAzimuth");
      break;
    case 8:
      strcpy(sds_name,"SensorAzimuth");
      break;
  }
  
  /*--------------------------------------------------------------------------*/
  /* Convert the name of the SDS of interest to an index number.              */
  /*--------------------------------------------------------------------------*/
  
  sds_index = SDnametoindex(id,sds_name);
    
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
  
  if( option < 6 || option > 6) {
    strcpy(att_name,"valid_range");
    att_index = SDfindattr(sd_id,att_name);
    status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDattrinfo - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }
    
    if ((attr1 = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
    error_allo(rout,"attr-valid_range");
    status = SDreadattr(sd_id,att_index,attr1);
    if (status == FAIL) {
      fprintf(stderr,"%s%s-There was an error when calling SDreadattr - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
      exit(EXIT_FAILURE);
    }
    
  }  
  
  /*--------------------------------------------------------------------------*/
  /* Read in the fill_value attribute.					      */
  /*--------------------------------------------------------------------------*/
  
  strcpy(att_name,"_FillValue");
  att_index = SDfindattr(sd_id,att_name);
  status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-There was an error when calling SDattrinfo - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
    exit(EXIT_FAILURE);
  }
      
  if ((attr2 = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
    error_allo(rout,"attr-_FillValue");
  status = SDreadattr(sd_id,att_index,attr2);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-There was an error when calling SDreadattr - %s, option=%d.",EXE_PROMPT,rout,att_name,option);
    exit(EXIT_FAILURE);
  }
    
  /*--------------------------------------------------------------------------*/
  /* Set arrays that are needed when reading in the actual data with          */
  /* dimension information.						      */
  /*--------------------------------------------------------------------------*/
  
  status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot get info about specified SDS, option=%d\n",EXE_PROMPT,rout,option);
    exit(EXIT_FAILURE);
  }
  
  start[0] = 0;
  start[1] = 0;
  
  edge[0] = dimen[0];
  edge[1] = dimen[1]-pix_offset;
    
  npts = edge[0]*edge[1];
  *nlines = dimen[0];
  *npixels = dimen[1]-pix_offset;
  
  /*--------------------------------------------------------------------------*/
  /* Read latitude or longitude.`	                                      */
  /*--------------------------------------------------------------------------*/
  
  if (option == 1 || option == 2) {
    
    valid_range_f32 = (float32 *) attr1;
    fill_val_f32 = (float32 *) attr2;
    
    if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
      error_allo(rout,"buf_f32");
    *buf = buf_f32;
    
    /*------------------------------------------------------------------------*/
    /* Read in the data to a buffer.		                              */
    /*------------------------------------------------------------------------*/
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_f32);
    
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
    
    for (n=0; n<npts; n++) {
      if (buf_f32[n] <  valid_range_f32[0] || buf_f32[n] >  valid_range_f32[1]) {
	buf_f32[n] = MISSING_FLOAT;
	badmask[n] = YES;
      }
    }
    
  }
  
  /*--------------------------------------------------------------------------*/
  /* Read SZA or VZA.	                                                      */
  /*--------------------------------------------------------------------------*/
  
  else if (option >= 3 && option <= 5) {
    
    valid_range_i16 = (int16 *) attr1;
    fill_val_i16 = (int16 *) attr2;
    
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
    
    if ((buf_i16 = (int16 *) malloc(npts*sizeof(int16))) == NULL)
      error_allo(rout,"buf_i16");
    if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
      error_allo(rout,"buf_f32");
    *buf = buf_f32;
    
    /*------------------------------------------------------------------------*/
    /* Read in the data to a buffer.		                              */
    /*------------------------------------------------------------------------*/
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
    
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
    
//  printf("Scale Factor = %f\n",scale_fac_f64);

    for (n=0; n<npts; n++) {
      buf_f32[n] = (float32) buf_i16[n] * scale_fac_f64[0];
    }
    
    free(buf_i16);
    free(scale_fac_f64);
    
  }
  
  /*--------------------------------------------------------------------------*/
  /* Read Solar and Sensor azimuth and compute relative azimuth.              */
  /*--------------------------------------------------------------------------*/
  
  if (option == 5) {
          
    sds_index = SDnametoindex(id,"SensorAzimuth");
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
    
    if ((buf_i16 = (int16 *) malloc(npts*sizeof(int16))) == NULL)
      error_allo(rout,"buf_i16");
    if ((buf_f32_2 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
      error_allo(rout,"buf_f32_2");
    
    /*------------------------------------------------------------------------*/
    /* Read in the data to a buffer.		                              */
    /*------------------------------------------------------------------------*/
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
    
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
    
    for (n=0; n<npts; n++) {
      buf_f32_2[n] = (float32) buf_i16[n] * scale_fac_f64[0];
/*    float32 a = buf_f32[n]; */
      
      /* Relative azimuth calculation as used in MOD35 */ 
      buf_f32[n] = fabs(180.0 - fabs((buf_f32_2[n] - buf_f32[n])));
/*      printf("buf: %d %f %f %f\n", n, a, buf_f32_2[n], buf_f32[n]); */
    }
    
    free(buf_i16);
    free(buf_f32_2);
    free(scale_fac_f64);
    
  }    
  
  /* Read land/sea tag. */
  
  if( option == 6) {
	  
	  	  
	  if ((bufc8 = (unsigned char *) malloc(npts*sizeof(unsigned char))) == NULL)
	        error_allo(rout,"bufc8");
	    *buf = bufc8;
	  
	    /*------------------------------------------------------------------------*/
	    /* Read in the data to a buffer.		                              */
	    /*------------------------------------------------------------------------*/
	    
	    status = SDreaddata(sd_id,start,NULL,edge,bufc8);
	    	    
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
	    
/*      for (n=0; n<100; n++) {
          printf("bufc8: %d \n", bufc8[n]);
	    }
*/   
	  
  }
  
  /*--------------------------------------------------------------------------*/
  /* Read Solar azimuth                                                       */
  /*--------------------------------------------------------------------------*/
  
  if( option == 7) {

    valid_range_i16 = (int16 *) attr1;
    fill_val_i16 = (int16 *) attr2;
    
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
    
    if ((buf_i16 = (int16 *) malloc(npts*sizeof(int16))) == NULL)
      error_allo(rout,"buf_i16");
    if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
      error_allo(rout,"buf_f32");
    *buf = buf_f32;
    
    /*------------------------------------------------------------------------*/
    /* Read in the data to a buffer.		                              */
    /*------------------------------------------------------------------------*/
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
    
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
    
//  printf("Scale Factor = %f\n",scale_fac_f64);

    for (n=0; n<npts; n++) {
      buf_f32[n] = (float32) buf_i16[n] * scale_fac_f64[0];
    }
    
    free(buf_i16);
    free(scale_fac_f64);    
  }    

  /*--------------------------------------------------------------------------*/
  /* Read Sensor azimuth                                                      */
  /*--------------------------------------------------------------------------*/

  if( option == 8) {

    valid_range_i16 = (int16 *) attr1;
    fill_val_i16 = (int16 *) attr2;
    
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
    
    if ((buf_i16 = (int16 *) malloc(npts*sizeof(int16))) == NULL)
      error_allo(rout,"buf_i16");
    if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
      error_allo(rout,"buf_f32");
    *buf = buf_f32;
    
    /*------------------------------------------------------------------------*/
    /* Read in the data to a buffer.		                              */
    /*------------------------------------------------------------------------*/
    
    status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
    
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
    
//  printf("Scale Factor = %f\n",scale_fac_f64);

    for (n=0; n<npts; n++) {
      buf_f32[n] = (float32) buf_i16[n] * scale_fac_f64[0];
    }
    
    free(buf_i16);
    free(scale_fac_f64);
  }    

  /*--------------------------------------------------------------------------*/
  /* Close the HDF file.                                                      */
  /*--------------------------------------------------------------------------*/
  
  if (strcmp(filename,"OPEN_1000") != 0) {
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

/******************************************************************************/
/******************************************************************************/
