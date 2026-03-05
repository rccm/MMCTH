/*$Id: modis_level06_read.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

/******************************************************************************/
/* This routine reads in the MODIS cloud mask from an HDF.                    */
/*                                                                            */
/* Input:                                                                     */
/*         filename: the name of the 06 cloud product file.                     */
/*                                                                            */
/* Author:                                                                    */
/*         Michael Pavolonis					              */
/*         CIMSS/SSEC			  			              */
/******************************************************************************/


#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"

/*----------------------------------------------------------------------------*/

void modis_level06_read (char *filename, int option, void **buf, int *nlines, 
                         int *npixels)
  
{
    
  void
    *attr;
    
  int16
    *buf_i16;
  
  char8
    sds_name[MAX_NC_NAME],
    sds_name_check[MAX_NC_NAME],
    att_name[MAX_NC_NAME],

    *rout = {"modis_level06_read"};
  
  uint16
 
    *valid_range,
    *fill_val;
    
  int32
    status,
    id,
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
    
  float32
    *buf_f32;
    
  float64
    *scale_fac,
    *offset;
  
  switch (option) {
    case 1:
      strcpy(sds_name,"Effective_Particle_Radius");
      break;
    case 2:
      strcpy(sds_name,"Cloud_Optical_Thickness");
      break; 
    default:
        fprintf(stderr,"%s%s - Option, %d is not valid-aborting\n",EXE_PROMPT,rout,option);
        exit(EXIT_FAILURE);
  }
  
  /*--------------------------------------------------------------------------*/
  /* Open the HDF file.                                                       */
  /*--------------------------------------------------------------------------*/
  
  id = SDstart(filename, DFACC_READ);
  if (id == FAIL) {
    fprintf(stderr,"%s%s - Invalid HDF file:\n %s\n",EXE_PROMPT,rout,filename);
    exit(EXIT_FAILURE);
  }
 
  /*--------------------------------------------------------------------------*/
  /* Convert the name of the SDS of interest to an index number.              */
  /*--------------------------------------------------------------------------*/
  
  sds_index = SDnametoindex(id,sds_name);
  sd_id = SDselect(id,sds_index);
  if (sd_id == FAIL) {
    fprintf(stderr,"%s%s - Invalid SDS, %s - aborting\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }
  status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot get info about specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }
  if (data_type != DFNT_INT16) {
    fprintf(stderr,"%s%s-SDS, %s is not of type INT16 as expected\n",EXE_PROMPT,rout,sds_name);
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
	
  scale_fac = (float64 *) attr;
  
  strcpy(att_name,"add_offset");
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
        
  offset = (float64 *) attr;
  
  
  /*--------------------------------------------------------------------------*/
  /* Set arrays that are needed when reading in the actual data with          */
  /* dimension information.						      */
  /*--------------------------------------------------------------------------*/
  
  start[0] = 0;
  start[1] = 0;
  
  edge[0] = dimen[0];
  edge[1] = dimen[1];
    
  npts = edge[0]*edge[1];
  *nlines = dimen[0];
  *npixels = dimen[1];
  
  /*------------------------------------------------------------------------*/
  /* Allocate memory for the data to be read in from the file.	            */
  /*------------------------------------------------------------------------*/
  
  if ((buf_i16 = (int16 *) malloc(npts*sizeof(int16))) == NULL)
    error_allo(rout,"buf_i16");
  
  if ((buf_f32 = (float32 *) malloc(npts*sizeof(float32))) == NULL)
    error_allo(rout,"buf_f32");
  *buf = buf_f32;
  
  /*------------------------------------------------------------------------*/
  /* Read in the data to a buffer.		                            */
  /*------------------------------------------------------------------------*/
  
  status = SDreaddata(sd_id,start,NULL,edge,buf_i16);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot read specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }   
    
  /*------------------------------------------------------------------------*/
  /* Detach from the scientific data set.                                   */
  /*------------------------------------------------------------------------*/
  
  status = SDendaccess(sd_id);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot end access to specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }
  
  /*--------------------------------------------------------------------------*/
  /* Close the HDF file.                                                      */
  /*--------------------------------------------------------------------------*/
  
  status = SDend(id);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot close %s\n",EXE_PROMPT,rout,filename);
    exit(EXIT_FAILURE);
  }
  
  /*------------------------------------------------------------------------*/
  /* Perform the necessary unscaling.                                    */
  /*------------------------------------------------------------------------*/
  
  for (n=0; n<npts; n++) {
    buf_f32[n] = scale_fac[0] * (buf_i16[n] - offset[0]);
  }
  
  free(buf_i16);
  free(scale_fac);
  free(offset);
  free(valid_range);
  free(fill_val);
  
  /*--------------------------------------------------------------------------*/
  /* DONE with the routine so exit.                                           */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
