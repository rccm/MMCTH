/* leocat_to_mod06.c */
/* Copies Collection 6 LEOCAT Cloud Top (MOD06) data to the Collection 5 Cloud Top Output file */

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
   int32 start[2], edge[2], dimen[2];

   /* Matlab doesn't like names that start with numbers, so the reason for 2 varNames */
   char * varNames1[9] = {"1km_cloud_top_temperature",
                            "1km_cloud_top_pressure",
                            "1km_cloud_top_height",
                            "1km_cloud_emissivity",
                            "1km_surface_temperature",
                            "1km_cld_emiss11",
                            "1km_cld_emiss12",
                            "1km_cld_emiss13",
                            "1km_cld_emiss85"};
   char * varNames2[9] = {"cloud_top_temperature_1km",
                            "cloud_top_pressure_1km",
                            "cloud_top_height_1km",
                            "cloud_emissivity_1km",
                            "surface_temperature_1km",
                            "cloud_emiss11_1km",
                            "cloud_emiss12_1km",
                            "cloud_emiss13_1km",
                            "cloud_emiss85_1km"};
   /* scales and offsets match up with the varNames above */
   float64 scale[9]={0.01,0.1,1.0,0.01,0.01,0.0001,0.0001,0.0001,0.0001};
   float64 offset[9]={-15000.0,0.0,0.0,0.0,-15000.0,0.0,0.0,0.0,0.0};
   int fillValues[9] = {-999,-999,-999,127,-999,-999,-999,-999,-999};
   int isShort[9]={1,1,1,0,1,1,1,1,1};
   char * byteNames1[4] = {"1km_cloud_top_method",
                          "1km_modis_C6_IRP",
                          "1km_os_top_flag",
                          "1km_IRP_CTH_Consistency_Flag"};
   char * byteNames2[4] = {"cloud_top_method_1km",
                          "Cloud_Phase_Infrared_1km",
                          "os_top_flag_1km",
                          "IRP_CTH_Consistency_Flag_1km"};
   unsigned char *byte_arr;
   float32 *float_arr;
   int16 *short_arr;
   int i, n;


   if (argc != 3) {
       fprintf(stdout,"To execute, type: MOD_PRL2M06.exe <LEOCAT Collection 6 Cloud Top HDF file> <Collection 5 Cloud Top File>\n");
       return(0);
   }

   leocatFilename[0]=0;
   modisFilename[0]=0;

   strcat(leocatFilename,argv[1]);
   strcat(modisFilename,argv[2]);

    /* Select proper SDS and read the data */
   sd_id = SDstart(leocatFilename, DFACC_READ);
   if (sd_id <= 0) {
      fprintf(stderr, "Could not open LEOCAT HDF file %s\n", leocatFilename);
      exit (-1);
   }
   /* Open Collection 5 cloud top data file for writing */
   hdf.sd_id = SDstart(modisFilename,DFACC_WRITE);
   if (hdf.sd_id == FAIL) {
        fprintf(stderr,"%sCannot open Level 2 output file:\n%s%s - aborting\n",
           EXE_PROMPT, EXE_PROMPT,modisFilename);
        exit(EXIT_FAILURE);
   }

   for (i=0; i<9; i++)
   {
     (void) strcpy(sds_name, varNames1[i]);
     fprintf(stdout,"Processing %s array\n",sds_name);

     sds_index = SDnametoindex(sd_id, sds_name);
     sds_id = SDselect(sd_id, sds_index);
     status = SDgetinfo(sds_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
     if (status == FAIL) {
      fprintf(stderr,"Cannot get info about %s SDS from %s\n",sds_name,leocatFilename);
      exit(EXIT_FAILURE);
     }
     start[0] = 0;
     start[1] = 0;
     edge[0] = dimen[0];
     edge[1] = dimen[1];
     fprintf(stderr,"edge=%d %d, rank=%d\n",edge[0],edge[1],rank);

     float_arr = (float32*) malloc(edge[0] * edge[1] *sizeof(float32));
     status = SDreaddata(sds_id, start, NULL, edge, float_arr);
     if (status == FAIL){
      fprintf(stdout,"Collection 6 %s not read from %s\n",sds_name,leocatFilename);
      exit(EXIT_FAILURE);
     }
     status = SDendaccess(sds_id);
     if (status == FAIL) {
       fprintf(stderr,"%sError ending access to Collection 6 SDS, %s - aborting\n",EXE_PROMPT,sds_name);
       exit(EXIT_FAILURE);
     }
      
     fprintf(stdout,"Collection 6 %s read successfully from %s\n",sds_name,leocatFilename);

     /* edge values should be correct */
     (void) strcpy(hdf.sds_name,varNames2[i]);
     sds_index = SDnametoindex(hdf.sd_id,hdf.sds_name);
     sds_id = SDselect(hdf.sd_id,sds_index);

     /* write output to the Collection 5 Output file - short or a byte */
     /* LEOCAT fill values in this case are -999, but need to be different for Collection 5 Cloud Top Output file */ 
     hdf.scaled_flg = no_scale;
     hdf.rank = rank;
     if (isShort[i])
     {
       short_arr = (int16*) malloc(edge[0]*edge[1]*sizeof(int16));
       hdf.scaled_type = DFNT_INT16;
       for (n=0;  n<edge[0] * edge[1]; n++)
          if (nint(float_arr[n])!=-999)
            short_arr[n]=nint(float_arr[n]/scale[i] + offset[i]);
          else
            short_arr[n]=fillValues[i];
       write_hdf_multidim(hdf, start, edge, &sds_id, (void *) short_arr);
       free(short_arr);
     }
     else
     {
       byte_arr = (unsigned char*) malloc(edge[0] * edge[1]);
       hdf.scaled_type = DFNT_INT8;
       for (n=0;  n<edge[0] * edge[1]; n++)
          if (nint(float_arr[n])!=-999)
            byte_arr[n]=nint(float_arr[n]/scale[i] + offset[i]);
          else
            byte_arr[n]=fillValues[i];

       write_hdf_multidim(hdf, start, edge, &sds_id, (void *) byte_arr);
       free(byte_arr);
     }


     status = SDendaccess(sds_id);
     if (status == FAIL) {
       fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
       exit(EXIT_FAILURE);
     }
     fprintf(stdout,"%s overwritten successfully in %s\n",hdf.sds_name,modisFilename);
     free(float_arr);
   }

   /* Read in Collection 6 byte array variables*/
   for (i=0; i<4; i++)
   {
     (void) strcpy(sds_name, byteNames1[i]);
     fprintf(stdout,"Processing %s array\n",sds_name);
   
     sds_index = SDnametoindex(sd_id, sds_name);
     sds_id = SDselect(sd_id, sds_index);
     status = SDgetinfo(sds_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
     if (status == FAIL) {
      fprintf(stderr,"Cannot get info about %s SDS from %s\n",sds_name,leocatFilename);
      exit(EXIT_FAILURE);
     }
     start[0] = 0;
     start[1] = 0;
     edge[0] = dimen[0];
     edge[1] = dimen[1];
     fprintf(stderr,"edge=%d %d\n",edge[0],edge[1]);

     byte_arr = (unsigned char*) malloc(edge[0] * edge[1]);
     status = SDreaddata(sds_id, start, NULL, edge, byte_arr);
     if (status == FAIL){
      fprintf(stdout,"Collection 6 %s field not read from %s\n",sds_name,leocatFilename);
      exit(EXIT_FAILURE);
     }
      
     fprintf(stdout,"Collection 6 %s field read successfully from %s\n",sds_name,leocatFilename);

     /* Write Collection 6 byte array variable out to existing Collection 5 HDF file */
     (void) strcpy(hdf.sds_name,byteNames2[i]);
     sds_index = SDnametoindex(hdf.sd_id,hdf.sds_name);
     sds_id = SDselect(hdf.sd_id,sds_index);
     hdf.scaled_flg = scaled_flg_solm;
     hdf.scaled_type = scaled_type_solm_pixel;
     hdf.rank = rank;

     write_hdf_multidim(hdf, start, edge, &sds_id, (void *) byte_arr);

     status = SDendaccess(sds_id);
     if (status == FAIL) {
       fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
       exit(EXIT_FAILURE);
     }
     fprintf(stdout,"%s overwritten successfully in %s\n",hdf.sds_name,modisFilename);
     free(byte_arr);
   }

   /* Close LEOCAT file. */
   sd_id = SDend(sd_id);
   /* Close MODIS file */
   hdf.sd_id = SDend(hdf.sd_id);

   return(0);
}
