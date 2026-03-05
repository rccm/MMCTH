/*$Id: modis35_cloud_mask.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

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


void modis_level35_read (char *, int, void **, int *, int *);

void modis35_cloud_mask (unsigned char verbose,
                         imagerL1* imgr1, imagerL2* imgr2,
			 sounderL1* sndr1, sounderL2* sndr2,
			 nwp_params nwp, rtm_profiles** rtm,
			 rtm_toa rclr, rad_utils rutil)
{
  FILE *fptr;

  char filename[MAX_STR_LEN],
       sys_command1[MAX_STR_LEN], sys_command2[MAX_STR_LEN],
       sys_command3[MAX_STR_LEN], sys_command[MAX_STR_LEN],
       mask_dir[MAX_STR_LEN], mask_file[MAX_STR_LEN];
  int yy=0;

  sprintf(sys_command,"echo nofile > temp_flist");
  system(sys_command);

  if (verbose == YES)
    fprintf(stdout,"%sGenerating MODIS cloud mask file name\n",EXE_PROMPT);

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
      sprintf(sys_command1,"%s/*%2.2d%3.3d.%2.2d%2.2d.mod35.hdf*",imgr1->cmask_dir_name,yy,imgr1->jday,imgr1->hour,imgr1->minute);
      sprintf(sys_command2,"%s/*%2.2d%3.3d.%2.2d%2.2d.mod06.hdf*",imgr1->cmask_dir_name,yy,imgr1->jday,imgr1->hour,imgr1->minute);
      sprintf(sys_command3,"%s/*%2.2d%3.3d.%2.2d%2.2d.mod06.hdf*",imgr1->cprod_dir_name,yy,imgr1->jday,imgr1->hour,imgr1->minute);
      sprintf(sys_command,"ls %s %s %s > temp_flist",sys_command1,sys_command2,sys_command3);
      system(sys_command);
      if ((fptr = fopen("temp_flist","r")) == NULL) {
        fprintf(stderr,"%sCannot open cloud mask filename text file - temp_flist - aborting\n",EXE_PROMPT);
	exit(EXIT_FAILURE);
      }
      fscanf(fptr,"%s",&filename[0]);
      fclose(fptr);
      break;

    /*----------------------------------------------------------------------------
      Aqua
    ----------------------------------------------------------------------------*/

    case AQUA:
      sprintf(sys_command1,"%s/*%2.2d%3.3d.%2.2d%2.2d.mod35.hdf*",imgr1->cmask_dir_name,yy,imgr1->jday,imgr1->hour,imgr1->minute);
      sprintf(sys_command2,"%s/*%2.2d%3.3d.%2.2d%2.2d.mod06.hdf*",imgr1->cmask_dir_name,yy,imgr1->jday,imgr1->hour,imgr1->minute);
      sprintf(sys_command3,"%s/*%2.2d%3.3d.%2.2d%2.2d.mod06.hdf*",imgr1->cprod_dir_name,yy,imgr1->jday,imgr1->hour,imgr1->minute);
      sprintf(sys_command,"ls %s %s %s > temp_flist",sys_command1,sys_command2,sys_command3);
      system(sys_command);
      if ((fptr = fopen("temp_flist","r")) == NULL) {
        fprintf(stderr,"%sCannot open cloud mask filename text file - temp_flist - aborting\n",EXE_PROMPT);
	exit(EXIT_FAILURE);
      }
      fscanf(fptr,"%s",&filename[0]);
      fclose(fptr);
      break;
    default:
      fprintf(stderr,"%sMODIS L35 cloud mask can only be read if satellite = Terra or Aqua, not %s - aborting\n",EXE_PROMPT,imgr1->satname);
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
        sprintf(sys_command1,"%s/MODATML2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cmask_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
	sprintf(sys_command2,"%s/MODATML2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cprod_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
        sprintf(sys_command,"ls %s %s > temp_flist",sys_command1,sys_command2);
        system(sys_command);
      }

      /*----------------------------------------------------------------------------
        Standard 1 km
      ----------------------------------------------------------------------------*/

      else {

//      Old way was to search 'cmask_dir_name' for correct date and time.  For MODAPS, input file name comes directly
//      from command line, stored in 'imgr1->cmask_file_name'.

/*
        sprintf(sys_command1,"%s/MOD35_L2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cmask_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
        sprintf(sys_command2,"%s/MOD06.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cmask_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
        sprintf(sys_command3,"%s/MOD06.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cprod_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
//      sprintf(sys_command,"ls %s %s %s > temp_flist",sys_command1,sys_command2,sys_command3);
        sprintf(sys_command,"ls %s > temp_flist",sys_command1);
        system(sys_command);
*/

      }

/*
      if ((fptr = fopen("temp_flist","r")) == NULL) {
        fprintf(stderr,"%sCannot open cloud mask filename text file - temp_flist - aborting\n",EXE_PROMPT);
	exit(EXIT_FAILURE);
      }
      fscanf(fptr,"%s",&filename[0]);
      fclose(fptr);
*/

      printf("cloud mask directory: %s\n", imgr1->cmask_dir_name);
      printf("cloud mask file: %s\n", imgr1->cmask_file_name);
      if(strcmp(imgr1->cmask_dir_name, UNKNOWN_STR) == 0) {
        strcpy(filename, imgr1->cmask_file_name);
      }
      else {
        strcpy(mask_dir, imgr1->cmask_dir_name);
        strcpy(mask_file, imgr1->cmask_file_name);
        strcat(mask_dir, "/");
        strcat(mask_dir, mask_file);
        strcpy(filename, mask_dir);
      }

      break;

    /*----------------------------------------------------------------------------
      Aqua
    ----------------------------------------------------------------------------*/

    case AQUA:

      /*----------------------------------------------------------------------------
        Subsetted (5 km)
      ----------------------------------------------------------------------------*/

      if (imgr1->data_source == MODIS5_DAT) {
        sprintf(sys_command1,"%s/MYDATML2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cmask_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
	sprintf(sys_command2,"%s/MYDATML2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cprod_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
        sprintf(sys_command,"ls %s %s > temp_flist",sys_command1,sys_command2);
        system(sys_command);
      }

      /*----------------------------------------------------------------------------
        Standard 1 km
      ----------------------------------------------------------------------------*/

      else {

//      Old way was to search 'cmask_dir_name' for correct date and time.  For MODAPS, input file name comes directly
//      from command line, stored in 'imgr1->cmask_file_name'.

/*
        sprintf(sys_command1,"%s/MYD35_L2.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cmask_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
        sprintf(sys_command2,"%s/MYD06.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cmask_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
        sprintf(sys_command3,"%s/MYD06.A%d%3.3d.%2.2d%2.2d*.hdf*",imgr1->cprod_dir_name,imgr1->year,imgr1->jday,imgr1->hour,imgr1->minute);
//      sprintf(sys_command,"ls %s %s %s > temp_flist",sys_command1,sys_command2,sys_command3);
        sprintf(sys_command,"ls %s > temp_flist",sys_command1);
        system(sys_command);
*/
      }

/*
      if ((fptr = fopen("temp_flist","r")) == NULL) {
        fprintf(stderr,"%sCannot open cloud mask filename text file - temp_flist - aborting\n",EXE_PROMPT);
	exit(EXIT_FAILURE);
      }
      fscanf(fptr,"%s",&filename[0]);
      fclose(fptr);
*/

      printf("cloud mask directory: %s\n", imgr1->cmask_dir_name);
      printf("cloud mask file: %s\n", imgr1->cmask_file_name);
      if(strcmp(imgr1->cmask_dir_name, UNKNOWN_STR) == 0) {
        strcpy(filename, imgr1->cmask_file_name);
      }
      else {
        strcpy(mask_dir, imgr1->cmask_dir_name);
        strcpy(mask_file, imgr1->cmask_file_name);
        strcat(mask_dir, "/");
        strcat(mask_dir, mask_file);
        strcpy(filename, mask_dir);
      }

      break;

    default:
      fprintf(stderr,"%sMODIS L35 cloud mask can only be read if satellite = Terra or Aqua, not %s - aborting\n",EXE_PROMPT,imgr1->satname);
      exit(EXIT_FAILURE);
      break;

    }
  }

   /*----------------------------------------------------------------------------
    Check to make sure a file was found.
  ----------------------------------------------------------------------------*/

  if (strlen(filename) == 0) {
    fprintf(stderr,"\n%sCannot find matching MODIS cloud mask file - aborting\n",EXE_PROMPT);
    exit(EXIT_FAILURE);
  }

  /*----------------------------------------------------------------------------
    Read in cloud mask.
  ----------------------------------------------------------------------------*/

  free(imgr2->cldmask);
  free(imgr2->landmask);
  if (verbose == YES)
    fprintf(stdout,"%sReading cloud mask file:\n%s\n",EXE_PROMPT,filename);
  modis_level35_read(filename, 1, (void **) &imgr2->cldmask, &imgr2->nrow, &imgr2->ncol);
  modis_level35_read(filename, 2, (void **) &imgr2->landmask, &imgr2->nrow, &imgr2->ncol);

}
