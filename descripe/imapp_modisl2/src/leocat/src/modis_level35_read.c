/*$Id: modis_level35_read.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

/******************************************************************************/
/* This routine reads in the MODIS cloud mask from an HDF.                    */
/*                                                                            */
/* Input:                                                                     */
/*         filename: the name of the L35 cloud mask file.                     */
/*                                                                            */
/* Author:                                                                    */
/*         Michael Pavolonis					              */
/*         CIMSS/SSEC			  			              */
/******************************************************************************/


#include <math.h>
#include <hdf.h>
#include <mfhdf.h>

#include "common_leocat.h"

typedef unsigned int short int48;

#define M_BITS_CMASK 0006
#define M_BITS_LMASK 0300
#define M_BITS_SMASK 0040
#define M_BITS_GMASK 0020

#define N_BIT_SHIFT_CMASK 1
#define N_BIT_SHIFT_LMASK 6
#define N_BIT_SHIFT_SMASK 5
#define N_BIT_SHIFT_GMASK 4

/*----------------------------------------------------------------------------*/

void modis_level35_read (char *filename, int option, void **buf, int *nlines, 
                         int *npixels)
  
{
    
  uchar8
    *buf_uc8;
  
  char8
  
    sds_name_check[MAX_NC_NAME],

  

   
    *rout = {"modis_level35_read"};
  
  int
    M_BITS,
    N_BIT_SHIFT;
  
  int32
    status,
    id,
    sd_id,
    sds_index,
    data_type,
    npts,
    n,
    rank,
    nattr,
    dimen[MAX_VAR_DIMS],
    start[3],
    edge[3];		 
  
  switch (option) {
    case 1:
      M_BITS = M_BITS_CMASK;
      N_BIT_SHIFT = N_BIT_SHIFT_CMASK;
      break;
    case 2:
      M_BITS = M_BITS_LMASK;
      N_BIT_SHIFT = N_BIT_SHIFT_LMASK;
      break; 
    case 3:
      M_BITS = M_BITS_SMASK;
      N_BIT_SHIFT = N_BIT_SHIFT_SMASK;
      break;
    case 4:
      M_BITS = M_BITS_GMASK;
      N_BIT_SHIFT = N_BIT_SHIFT_GMASK;
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
  
  sds_index = SDnametoindex(id,"Cloud_Mask");
  sd_id = SDselect(id,sds_index);
  if (sd_id == FAIL) {
    fprintf(stderr,"%s%s - Invalid SDS, Cloud_Mask - aborting\n",EXE_PROMPT,rout);
    exit(EXIT_FAILURE);
  }
  status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot get info about specified SDS, Cloud_Mask\n",EXE_PROMPT,rout);
    exit(EXIT_FAILURE);
  }
  
  /*--------------------------------------------------------------------------*/
  /* Set arrays that are needed when reading in the actual data with          */
  /* dimension information.						      */
  /*--------------------------------------------------------------------------*/
  
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  
  edge[0] = 1;
  edge[1] = dimen[1];
  edge[2] = dimen[2];
    
  npts = edge[1]*edge[2];
  *nlines = dimen[1];
  *npixels = dimen[2];
  
  /*------------------------------------------------------------------------*/
  /* Allocate memory for the data to be read in from the file.	            */
  /*------------------------------------------------------------------------*/
  
  if ((buf_uc8 = (uchar8 *) malloc(npts*sizeof(uchar8))) == NULL)
    error_allo(rout,"buf_uc8");
  *buf = buf_uc8;
  
  /*------------------------------------------------------------------------*/
  /* Read in the data to a buffer.		                            */
  /*------------------------------------------------------------------------*/
  
  status = SDreaddata(sd_id,start,NULL,edge,buf_uc8);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot read specified SDS, Cloud_Mask\n",EXE_PROMPT,rout);
    exit(EXIT_FAILURE);
  }   
    
  /*------------------------------------------------------------------------*/
  /* Detach from the scientific data set.                                   */
  /*------------------------------------------------------------------------*/
  
  status = SDendaccess(sd_id);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot end access to specified SDS, Cloud_Mask\n",EXE_PROMPT,rout);
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
  /* Perform the necessary bit shifting.                                    */
  /*------------------------------------------------------------------------*/
  
  for (n=0; n<npts; n++) {
    buf_uc8[n] = (buf_uc8[n] & M_BITS) >> N_BIT_SHIFT;
  } 
  
  /*--------------------------------------------------------------------------*/
  /* DONE with the routine so exit.                                           */
  /*--------------------------------------------------------------------------*/

}

/******************************************************************************/
/******************************************************************************/
