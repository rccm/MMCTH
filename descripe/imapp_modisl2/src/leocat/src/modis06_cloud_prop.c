/*$Id: modis06_cloud_prop.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <hdf.h>
#include "common_leocat.h"
#include "nwp_leocat.h"
#include "rtm_leocat.h"
#include "radutils_leocat.h"
#include "imagerL1_leocat.h"
#include "imagerL2_leocat.h"
#include "imagerL3_leocat.h"
#include "sounderL1_leocat.h"
#include "sounderL2_leocat.h"
#include "sounderL3_leocat.h"


void modis_level06_read (char *, int, void **, int *, int *);

void modis06_cloud_prop (unsigned char verbose,
                         imagerL1* imgr1, imagerL2* imgr2, 
			 sounderL1* sndr1, sounderL2* sndr2,
			 nwp_params nwp, rtm_profiles** rtm, 
			 rtm_toa rclr, rad_utils rutil)
{
//  char *rout = {"modis06_cloud_prop.c"};
  
  FILE *fptr;
//  char8
//    sds_name[MAX_NC_NAME],
//    sds_name_check[MAX_NC_NAME],
//    att_name[MAX_NC_NAME],
//    scale_string[MAX_NC_NAME],
//    offset_string[MAX_NC_NAME],
//    *units;
  char filename[MAX_STR_LEN],
       sys_command1[MAX_STR_LEN],sys_command2[MAX_STR_LEN],sys_command[MAX_STR_LEN];
  int yy=0;
//  int32 id, status;
  
  sprintf(sys_command,"echo nofile > temp_flist");
  system(sys_command);
  
  if (verbose == YES)
    fprintf(stdout,"%sGenerating MODIS cloud product file list using UNIX commands\n",EXE_PROMPT);

  /*----------------------------------------------------------------------------
    Deal with direct broadcast data.
  ----------------------------------------------------------------------------*/
  
  if (imgr1->directbc_flg == YES) {
    if (imgr1->year > 100) yy = imgr1->year - ((imgr1->year/100)*100);
    switch(imgr1->satid) {
    
    /*----------------------------------------------------------------------------
      Terra
    ----------------------------------------------------------------------------*/
    
    case TERRA:
      sprintf(sys_command1,"%s/*%2.2d%3.3d.%2.2d%2.2d.mod06.hdf*",imgr1->cprod_dir_name,yy,imgr1->jday,imgr1->hour,imgr1->minute);
      sprintf(sys_command,"ls %s > temp_flist",sys_command1);
      system(sys_command);
      if ((fptr = fopen("temp_flist","r")) == NULL) {
        fprintf(stderr,"%sCannot open cloud product filename text file - temp_flist - aborting\n",EXE_PROMPT);
	exit(EXIT_FAILURE);
      }
      fscanf(fptr,"%s",&filename[0]);
      fclose(fptr);
      break;
    
    /*----------------------------------------------------------------------------
      Aqua
    ----------------------------------------------------------------------------*/
    
    case AQUA:
      sprintf(sys_command1,"%s/*%2.2d%3.3d.%2.2d%2.2d.mod06.hdf*",imgr1->cprod_dir_name,yy,imgr1->jday,imgr1->hour,imgr1->minute);
      sprintf(sys_command,"ls %s > temp_flist",sys_command1);
      system(sys_command);
      if ((fptr = fopen("temp_flist","r")) == NULL) {
        fprintf(stderr,"%sCannot open cloud mask filename text file - temp_flist - aborting\n",EXE_PROMPT);
	exit(EXIT_FAILURE);
      }
      fscanf(fptr,"%s",&filename[0]);
      fclose(fptr);
      break;
    default:   
      fprintf(stderr,"%sMODIS 06 cloud product can only be read if satellite = Terra or Aqua, not %s - aborting\n",EXE_PROMPT,imgr1->satname);
      exit(EXIT_FAILURE);
      break;
    }
  }
  
  /*----------------------------------------------------------------------------
    Deal with DAAC data.
  ----------------------------------------------------------------------------*/
  
  else {
    switch(imgr1->satid) {
    
    /*----------------------------------------------------------------------------
      Terra
    ----------------------------------------------------------------------------*/
    
    case TERRA:
      
      /*----------------------------------------------------------------------------
        Subsetted (5 km)
      ----------------------------------------------------------------------------*/
      
      if (imgr1->data_source == MODIS5_DAT) {
        sprintf(sys_command1,"%s/MODATML2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cprod_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
	sprintf(sys_command2,"%s/MODATML2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cmask_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
        sprintf(sys_command,"ls %s %s > temp_flist",sys_command1,sys_command2);
        system(sys_command);
      }
      
      /*----------------------------------------------------------------------------
        Standard 1 km
      ----------------------------------------------------------------------------*/
      
      else {
        sprintf(sys_command1,"%s/MOD06_L2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cprod_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
	sprintf(sys_command,"ls %s > temp_flist",sys_command1);
        system(sys_command);
      }
      if ((fptr = fopen("temp_flist","r")) == NULL) {
        fprintf(stderr,"%sCannot open cloud product filename text file - temp_flist - aborting\n",EXE_PROMPT);
	exit(EXIT_FAILURE);
      }
      fscanf(fptr,"%s",&filename[0]);
      fclose(fptr);
      break;
    
    /*----------------------------------------------------------------------------
      Aqua
    ----------------------------------------------------------------------------*/
    
    case AQUA:
      
      /*----------------------------------------------------------------------------
        Subsetted (5 km)
      ----------------------------------------------------------------------------*/
      
      if (imgr1->data_source == MODIS5_DAT) {
        sprintf(sys_command1,"%s/MYDATML2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cprod_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
	sprintf(sys_command2,"%s/MYDATML2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cmask_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
        sprintf(sys_command,"ls %s %s > temp_flist",sys_command1,sys_command2);
        system(sys_command);
      }
      
      /*----------------------------------------------------------------------------
        Standard 1 km
      ----------------------------------------------------------------------------*/
      
      else {
        sprintf(sys_command1,"%s/MYD06_L2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cprod_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
	sprintf(sys_command,"ls %s > temp_flist",sys_command1);
        system(sys_command);
      }
      if ((fptr = fopen("temp_flist","r")) == NULL) {
        fprintf(stderr,"%sCannot open cloud product filename text file - temp_flist - aborting\n",EXE_PROMPT);
	exit(EXIT_FAILURE);
      }
      fscanf(fptr,"%s",&filename[0]);
      fclose(fptr);
      break;
    default:   
      fprintf(stderr,"%sMODIS 06 cloud product can only be read if satellite = Terra or Aqua, not %s - aborting\n",EXE_PROMPT,imgr1->satname);
      exit(EXIT_FAILURE);
      break;
    }
  }

   /*----------------------------------------------------------------------------
    Check to make sure a file was found.
  ----------------------------------------------------------------------------*/
  
  if (strlen(filename) == 0) {
    fprintf(stderr,"\n%sCannot find matching MODIS cloud product file - aborting\n",EXE_PROMPT);
    exit(EXIT_FAILURE);
  }  
  
  /*----------------------------------------------------------------------------
    Read in cloud products.
  ----------------------------------------------------------------------------*/  
  
  free(imgr2->cod_vis);
  free(imgr2->cldreff);
  
  if (verbose == YES)
    fprintf(stdout,"%sReading cloud product file:\n%s\n",EXE_PROMPT,filename);
  modis_level06_read(filename, 1, (void **) &imgr2->cldreff, &imgr2->nrow, &imgr2->ncol);
  modis_level06_read(filename, 2, (void **) &imgr2->cod_vis, &imgr2->nrow, &imgr2->ncol);
  
  /*--------------------------------------------------------------------------*/
  /* DONE with the routine so exit.                                           */
  /*--------------------------------------------------------------------------*/
  
}
