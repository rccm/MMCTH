/* leocat_to_mod35.c */
/* Copies Collection 6 LEOCAT Cloud Mask to the Collection 5 Cloud Mask Output file */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"
#include "output_leocat.h"

void write_hdf_multidim(hdf_output, int32 *, int32 *, int32 *, void *);

int main(int argc, char* argv[]) {

    hdf_output hdf;
    char leocatFilename[256];
    char modisFilename[256];
    char sds_name[30];
    char sds_name_check[30];
    int32 sd_id, sds_id, sds_index, status;
    int32 rank, data_type, nattr;
    int32 start[3], stride[3], edge[3], dimen[3];
    int32 start1[3], stride1[3], edge1[3], dimen1[3];
    int32 flag1_byte_segment = 6;
    int32 flag2_byte_segment = 10;
    int32 num_250m_stats = 4;
    int32 numBandSPI = 2;
    int32 i,j,k,indx,count;
    int32 attrId;
    double scale_factor;
    double add_offset;

    unsigned char *cldmask_arr;
    unsigned char *qa_arr;
    float32 *stats_arr;
    float32 *cloud_mask_SPI;
    short *cloud_mask_SPI_Out;


    if (argc != 3) {
       fprintf(stdout,"To execute, type: leocat_to_mod35 <LEOCAT Collection 6 Cloud Mask HDF file> <Collection 5 Cloud Mask File>\n");
       return(0);
    }

    leocatFilename[0]=0;
    modisFilename[0]=0;

    strcat(leocatFilename,argv[1]);
    strcat(modisFilename,argv[2]);

    sd_id = SDstart(leocatFilename, DFACC_READ);
    if (sd_id <= 0) {
      fprintf(stderr, "Could not open LEOCAT HDF file %s\n", leocatFilename);
      exit (-1);
    }

    /* Select proper SDS and read the data */
   (void) strcpy(sds_name, "AlgKey29_flag_arr1");

   sds_index = SDnametoindex(sd_id, sds_name);
   sds_id = SDselect(sd_id, sds_index);
   status = SDgetinfo(sds_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
   if (status == FAIL) {
      fprintf(stderr,"Cannot get info about %s SDS from %s\n",sds_name,leocatFilename);
      exit(EXIT_FAILURE);
   }
   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   stride[0] = 1;
   stride[1] = 1;
   stride[2] = 1;
   edge[2] = dimen[2];
   edge[1] = dimen[1];
   edge[0] = flag1_byte_segment;
   cldmask_arr = (unsigned char*) malloc(flag1_byte_segment * edge[1] * edge[2]);

   status = SDreaddata(sds_id, start, stride, edge, cldmask_arr);
   if (status == FAIL){
      fprintf(stdout,"Collection 6 AlgKey29_flag_arr1 (cloud mask) not read from %s\n",leocatFilename);
      exit(EXIT_FAILURE);
   }
      
   fprintf(stdout,"Collection 6 AlgKey29_flag_arr1 (cloud mask) read successfully from %s\n",leocatFilename);

   /* Write Collection 6 cloud mask out to existing Collection 5 HDF file */
   hdf.sd_id = SDstart(modisFilename,DFACC_WRITE);
   if (hdf.sd_id == FAIL) {
          fprintf(stderr,"%sCannot open Level 2 output file:\n%s%s - aborting\n",
             EXE_PROMPT, EXE_PROMPT,modisFilename);
          exit(EXIT_FAILURE);
   }

   hdf.scaled_flg = scaled_flg_cldmask;
   hdf.scaled_type = scaled_type_cldmask_pixel;
   hdf.rank = 3;
   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   (void) strcpy(hdf.sds_name,"Cloud_Mask");
   sds_index = SDnametoindex(hdf.sd_id,hdf.sds_name);
   sds_id = SDselect(hdf.sd_id,sds_index);


   write_hdf_multidim(hdf, start, edge, &sds_id,
            (void *) cldmask_arr);
   status = SDendaccess(sds_id);
   if (status == FAIL) {
       fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
       exit(EXIT_FAILURE);
   }
   fprintf(stdout,"%s overwritten successfully in %s\n",hdf.sds_name,modisFilename);

   /* Read in Collection 6 QA array */
   (void) strcpy(sds_name, "AlgKey29_flag_arr2");
   
   sds_index = SDnametoindex(sd_id, sds_name);
   sds_id = SDselect(sd_id, sds_index);
   status = SDgetinfo(sds_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
   if (status == FAIL) {
      fprintf(stderr,"Cannot get info about %s SDS from %s\n",sds_name,leocatFilename);
      exit(EXIT_FAILURE);
   }
   edge[2] = flag2_byte_segment;
   edge[1] = dimen[1];
   edge[0] = dimen[0];
    
   qa_arr = (unsigned char*) malloc(flag2_byte_segment * edge[0] * edge[1]);

   status = SDreaddata(sds_id, start, stride, edge, qa_arr);
   if (status == FAIL){
      fprintf(stdout,"Collection 6 AlgKey29_flag_arr2 (QA) not read from %s\n",leocatFilename);
      exit(EXIT_FAILURE);
   }
      
   fprintf(stdout,"Collection 6 AlgKey29_flag_arr2 (QA) read successfully from %s\n",leocatFilename);

   /* Write Collection 6 QA out to existing Collection 5 HDF file */

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   (void) strcpy(hdf.sds_name,"Quality_Assurance");
   sds_index = SDnametoindex(hdf.sd_id,hdf.sds_name);
   sds_id = SDselect(hdf.sd_id,sds_index);

   write_hdf_multidim(hdf, start, edge, &sds_id, (void *) qa_arr);

   status = SDendaccess(sds_id);
   if (status == FAIL) {
       fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
       exit(EXIT_FAILURE);
   }
   fprintf(stdout,"%s overwritten successfully in %s\n",hdf.sds_name,modisFilename);
   /*fprintf(stdout, "Initialized %s: %d %d %s \n",hdf.sds_name, sds_id, hdf.dimen[0],
        hdf.dim_name[0]);*/ //GPC 
   
   /* Read in Collection 6 statistics array */
   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   (void) strcpy(sds_name, "AlgKey29_stats_250m");
   
   sds_index = SDnametoindex(sd_id, sds_name);
   sds_id = SDselect(sd_id, sds_index);
   status = SDgetinfo(sds_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
   if (status == FAIL) {
      fprintf(stderr,"Cannot get info about %s SDS from %s\n",sds_name,leocatFilename);
      /* exit(EXIT_FAILURE);*/
   }

   edge[2] = num_250m_stats;
   edge[1] = dimen[1];
   edge[0] = dimen[0];
   stats_arr = (float32 *) malloc(num_250m_stats * edge[0] * edge[1] * sizeof(float32));

   edge1[0] = dimen[0];
   edge1[1] = dimen[1];
   edge1[2] = numBandSPI;
   cloud_mask_SPI = (float32 *) malloc(numBandSPI * edge1[0] * edge1[1] * sizeof(float32));
   cloud_mask_SPI_Out = (short *) malloc(numBandSPI * edge1[0] * edge1[1] * sizeof(short));

   status = SDreaddata(sds_id, start, stride, edge, stats_arr);
   /* Should not exit(EXIT_FAILURE) here because Ref_250m_stats not available for night time case */
   if (status == FAIL){
      fprintf(stdout,"Collection 6 AlgKey29_stats_250m not read from %s\n",leocatFilename);
      /*exit(EXIT_FAILURE);*/
   }
   else {
      /* Write Collection 6 250 m statistics out to existing Collection 5 HDF file */
      start[0] = 0;
      start[1] = 0;
      start[2] = 0;
   
      (void) strcpy(hdf.sds_name,"Cloud_Mask_SPI");
      sds_index = SDnametoindex(hdf.sd_id,hdf.sds_name);
      sds_id = SDselect(hdf.sd_id,sds_index);
      
      attrId = SDfindattr(sds_id, "scale_factor");
      SDreadattr(sds_id, attrId, &scale_factor);

      attrId = SDfindattr(sds_id, "add_offset");
      SDreadattr(sds_id, attrId, &add_offset);
      
      start1[0] = 0;
      start1[1] = 0;
      start1[2] = 0;

      count = 0;
      indx = 0;
      for(i=0;i<edge[0];i++){
	for(j=0;j<edge[1];j++){

	  cloud_mask_SPI[count] = (stats_arr[indx+1]/stats_arr[indx])*100;
	  cloud_mask_SPI[count+1] = (stats_arr[indx+3]/stats_arr[indx+2])*100;
	  
	  
	  if(cloud_mask_SPI[count]<0){
	    cloud_mask_SPI_Out[count] = -9999;
	  }else{
	    cloud_mask_SPI_Out[count] = nint(cloud_mask_SPI[count]/scale_factor+add_offset);
	  }
	  
	  if(cloud_mask_SPI[count+1]<0){
	    cloud_mask_SPI_Out[count+1] = -9999;
	  }else{
	    cloud_mask_SPI_Out[count+1] = nint(cloud_mask_SPI[count+1]/scale_factor+add_offset);
	  }
	  
	  count = count+2;
	  indx = indx+4;
	}
      }

      write_hdf_multidim(hdf, start1, edge1, &sds_id, (void *) cloud_mask_SPI_Out);
      
      status = SDendaccess(sds_id);
      if (status == FAIL) {
          fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
          exit(EXIT_FAILURE);
      }
      fprintf(stdout,"%s overwritten successfully in %s\n",hdf.sds_name,modisFilename);
   }

   /* Close LEOCAT file. */
   sd_id = SDend(sd_id);
   /* Close MODIS file */
   hdf.sd_id = SDend(hdf.sd_id);
   fprintf(stdout, "Closing HDF access to %s and %s successful\n",leocatFilename,modisFilename);

   free(cldmask_arr);
   free(qa_arr);
   free(stats_arr);
   free(cloud_mask_SPI);
   free(cloud_mask_SPI_Out);
   return(0);
}
