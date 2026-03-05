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
    sds_name[MAX_NC_NAME],
    sds_name_check[MAX_NC_NAME],
    att_name[MAX_NC_NAME],
    scale_string[MAX_NC_NAME],
    offset_string[MAX_NC_NAME],
    *units;
  
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
        printf("\n\tOption input is not valid-exiting\n");
        exit(1);
  }
  
  /*--------------------------------------------------------------------------*/
  /* Open the HDF file.                                                       */
  /*--------------------------------------------------------------------------*/
  
  id = SDstart(filename, DFACC_READ);
  if (id == -1) {
    printf("\n\tERROR: Invalid HDF file: %s\n", filename);
    exit (1);
  }
  
  /*--------------------------------------------------------------------------*/
  /* Convert the name of the SDS of interest to an index number.              */
  /*--------------------------------------------------------------------------*/
  
  sds_index = SDnametoindex(id,"Cloud_Mask");
  sd_id = SDselect(id,sds_index);
  status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
  if (status == -1) {
    printf("\n\tThere was an error when calling SDgetinfo-exiting.\n");
    exit (1);
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
  
  buf_uc8 = (uchar8 *) malloc(npts*sizeof(uchar8));
  *buf = buf_uc8;
  
  /*------------------------------------------------------------------------*/
  /* Read in the data to a buffer.		                            */
  /*------------------------------------------------------------------------*/
  
  status = SDreaddata(sd_id,start,NULL,edge,buf_uc8);
  
  /*------------------------------------------------------------------------*/
  /* Print the status of the read: status 0 if successful, -1 if failed.    */
  /*------------------------------------------------------------------------*/
  
  if (status == 0) {
    /*printf("\n\tThe HDF scientific data was read successfully.\n");*/
  }
  else {
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
  /* Close the HDF file.                                                      */
  /*--------------------------------------------------------------------------*/
  
  status = SDend(id);
  if (status == -1) {
    printf("\n\tThere was an error when calling SDend-exiting.\n");
    exit (1);
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
