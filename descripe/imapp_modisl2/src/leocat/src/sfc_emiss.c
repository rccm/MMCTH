/*$Id: sfc_emiss.c,v 1.1.1.1 2006/12/05 14:27:49 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "rtm_leocat.h"
#include "nwp_leocat.h"

/*const float fill_emiss[] = {0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
			    0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
			    0.000, 0.978, 0.978, 0.978, 0.978, 0.978,
			    0.978, 0.000, 0.979, 0.979, 0.983, 0.988,
			    0.993, 0.986, 0.968, 0.968, 0.968. 0.968};*/

const float fill_emiss[] = {0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
			    0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
			    0.000, 0.985, 0.985, 0.985, 0.985, 0.985,
			    0.985, 0.000, 0.985, 0.985, 0.985, 0.985,
			    0.985, 0.985, 0.985, 0.985, 0.985, 0.985};

void seebor_channels(int32, int8 *, int32 *, int32 *, int32 *, rtm_toa *);
void seebor_channels_dateline(int32, int8 *, int32 *, int32 *, int32 *, int32 *, int32 *, int32 *, rtm_toa *);
void read_seebor_hdf(int32, char8 *, int32 *, int32 *, int32 *, float, void **);
void cat_emiss(int, int, int, float *, float *, float *);
void constant_channels(int8 *, rtm_toa *);
void set_emiss_arr(float *, float, long);

void read_seebor_emiss_imager(imagerL1 *imgr1, rtm_toa *rclr)
{
  char *rout = {"read_seebor_emiss_imager"};
  char filename[MAX_STR_LEN];
  int year, jday[] = {1,32,60,91,121,152,182,213,244,274,305,335};
  int32 id, status;
  int8 dateline_flg, polar_flg;
  int n, ilat1, ilat2, ilon1, ilon2, ilon1_2, ilon2_2, index;
  int32 istart[2], istride[2], iedge[2],
        istart_2[2], istride_2[2], iedge_2[2];

  if ((imgr1->index_rclr = (int *) malloc(imgr1->npts*sizeof(int))) == NULL)
    error_allo(rout,"imgr1->index_rclr");

  year = 2005;
  index = imgr1->month - 1 ;

  strcpy(filename, imgr1->sfc_dir_name);
  strcat(filename, "/global_emiss_intABI_2005001.hdf");

  printf("read_seebor_emiss_imager filename=%s\n",filename);

  dateline_flg = 0;
//    if (imgr1->bounds[1] > 180.0 || imgr1->bounds[3] > 180.0) dateline_flg = 1;
      if(imgr1->bounds[1] < 0.0 && imgr1->bounds[3] > 0.0 && imgr1->bounds[1] <= -90.0
		  && imgr1->bounds[3] > 90.0) dateline_flg = 1;

//   RAF added.
      polar_flg = 0;
      if(imgr1->bounds[0] <= -75.0 || imgr1->bounds[2] >= 75.0) polar_flg = 1;

  id = SDstart(filename, DFACC_READ);
  if (id == FAIL) {
    fprintf(stderr,"%sInvalid HDF file, %s\n",EXE_PROMPT,filename);
    exit(EXIT_FAILURE);
  }

  if (dateline_flg == 0 || polar_flg ==1) {
    ilat1 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->bounds[2]-first_lat_sfcemiss)/dlat_sfcemiss)));
    ilat2 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->bounds[0]-first_lat_sfcemiss)/dlat_sfcemiss)));

    ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->bounds[1]-first_lon_sfcemiss)/dlon_sfcemiss)));
    ilon2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->bounds[3]-first_lon_sfcemiss)/dlon_sfcemiss)));

    if (ilat1 > ilat2) swap_int_values(&ilat1, &ilat2);
    if (ilon1 > ilon2) swap_int_values(&ilon1, &ilon2);

    istart[0] = ilat1;
    istart[1] = ilon1;

    istride[0] = 1;
    istride[1] = 1;

    iedge[0] = (ilat2-ilat1)+1;
    iedge[1] = (ilon2-ilon1)+1;

    printf("edge: %f %f %d %d %d %d %d %d\n", imgr1->bounds[1],imgr1->bounds[3],ilat1,ilat2,ilon1,ilon2,iedge[0],iedge[1]);
    seebor_channels(id, imgr1->chflg, istart, istride, iedge, rclr);

    for (n=0; n<imgr1->npts; n++) {
      ilat1 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->lat[n]-first_lat_sfcemiss)/dlat_sfcemiss)));
      ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->lon[n]-first_lon_sfcemiss)/dlon_sfcemiss)));

      ilat2 = max(0,min((ilat1 - istart[0]),iedge[0]-1));
      ilon2 = max(0,min((ilon1 - istart[1]),iedge[1]-1));
      imgr1->index_rclr[n] = (iedge[1]*ilat2) + ilon2;

//    if(n > 2318049 && n < 2319051) {
//      printf("sfc_emiss: %d %d %d %d %d %d %d %d\n", n, ilat1, ilon1, ilat2, ilon2, imgr1->index_rclr[n], 
//             imgr1->index_nwp[n], imgr1->index_vza[n]);
//    }

    }

  }
  else {

    ilat1 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->bounds[2]-first_lat_sfcemiss)/dlat_sfcemiss)));
    ilat2 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->bounds[0]-first_lat_sfcemiss)/dlat_sfcemiss)));

//    ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->bounds[1]-first_lon_sfcemiss)/dlon_sfcemiss)));   //orig
    ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->bounds[3]-first_lon_sfcemiss)/dlon_sfcemiss)));
    ilon2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(180.0-first_lon_sfcemiss)/dlon_sfcemiss)));

    ilon1_2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(-180.0-first_lon_sfcemiss)/dlon_sfcemiss)));
//    ilon2_2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs((imgr1->bounds[3]-360.0)-first_lon_sfcemiss)/dlon_sfcemiss)));    //orig
    ilon2_2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->bounds[1]-first_lon_sfcemiss)/dlon_sfcemiss)));

    printf("first ilons dateline: %d %d %d %d\n", ilon1, ilon2, ilon1_2, ilon2_2);

    if (ilat1 > ilat2) swap_int_values(&ilat1, &ilat2);
    if (ilon1 > ilon2) swap_int_values(&ilon1, &ilon2);
    if (ilon1_2 > ilon2_2) swap_int_values(&ilon1_2, &ilon2_2);

    istart[0] = ilat1;
    istart[1] = ilon1;

    istride[0] = 1;
    istride[1] = 1;

    iedge[0] = (ilat2-ilat1)+1;
    iedge[1] = (ilon2-ilon1)+1;

    istart_2[0] = ilat1;
    istart_2[1] = ilon1_2;

    istride_2[0] = 1;
    istride_2[1] = 1;

    iedge_2[0] = (ilat2-ilat1)+1;
    iedge_2[1] = (ilon2_2-ilon1_2)+1;

    printf("dateline edge: %f %f %d %d %d %d %d %d %d %d %d %d\n", imgr1->bounds[1],imgr1->bounds[3],ilat1,ilat2,ilon1,ilon2,
    		iedge[0],iedge[1],ilon1_2,ilon2_2,iedge_2[0],iedge_2[1]);

    seebor_channels_dateline(id, imgr1->chflg, istart, istride, iedge, istart_2, istride_2, iedge_2, rclr);

    for (n=0; n<imgr1->npts; n++) {
//  for (n=0; n<10; n++) {
      ilat1 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->lat[n]-first_lat_sfcemiss)/dlat_sfcemiss)));
      ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->lon[n]-first_lon_sfcemiss)/dlon_sfcemiss)));

      if (imgr1->lon[n] >= 0.0) {
        ilat2 = max(0,min((ilat1 - istart[0]),iedge[0]-1));
        ilon2 = max(0,min((ilon1 - istart[1]),iedge[1]-1));
        imgr1->index_rclr[n] = iedge_2[1]*ilat2 + ((iedge[1]*ilat2) + ilon2);
      }
      else {
        ilat2 = max(0,min((ilat1 - istart_2[0]),iedge_2[0]-1));
        ilon2 = max(0,min((ilon1 - istart_2[1]),iedge_2[1]-1));
//        imgr1->index_rclr[n] = iedge[1]*ilat2 + ((iedge_2[1]*ilat2) + ilon2);
        imgr1->index_rclr[n] = iedge[1]*(ilat2+1) + ((iedge_2[1]*ilat2) + ilon2);
      }

//    printf("rclr index: %d %d %d %d %d %d\n", n, ilat1, ilat2, ilon1, ilon2, imgr1->index_rclr[n]);

    }
  }

  status = SDend(id);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot close %s\n",EXE_PROMPT,rout,filename);
    exit(EXIT_FAILURE);
  }

  if ((rclr->flag = (unsigned char *) malloc(rclr->npts*sizeof(unsigned char))) == NULL)
    error_allo(rout,"rclr->flag");
  for (n=0; n<rclr->npts; n++) rclr->flag[n] = NO;

}

void seebor_channels(int32 id, int8 *chflg, int32 istart[2], int32 istride[2], int32 iedge[2], rtm_toa *rclr)
{

  long int npts;
  char *rout = {"seebor_channels"};

  rclr->npts = iedge[0]*iedge[1];
  npts = rclr->npts;

  printf("seebor npts: %d\n", rclr->npts);

  if (chflg[19]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[19], (void **) &rclr->sfc_emiss20);
  }
  if (chflg[20]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[20], (void **) &rclr->sfc_emiss21);
  }
  if (chflg[21]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[21], (void **) &rclr->sfc_emiss22);
  }
  if (chflg[22]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[22], (void **) &rclr->sfc_emiss23);
  }
  if (chflg[23]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[23], (void **) &rclr->sfc_emiss24);
  }
  if (chflg[24]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[24], (void **) &rclr->sfc_emiss25);
  }
  if (chflg[26]) {
    read_seebor_hdf(id, "emiss9", istart, istride, iedge, fill_emiss[26], (void **) &rclr->sfc_emiss27);
  }
  if (chflg[27]) {
    read_seebor_hdf(id, "emiss10", istart, istride, iedge, fill_emiss[27], (void **) &rclr->sfc_emiss28);
  }
  if (chflg[28]) {
    read_seebor_hdf(id, "emiss11", istart, istride, iedge, fill_emiss[28], (void **) &rclr->sfc_emiss29);
  }
  if (chflg[29]) {
    read_seebor_hdf(id, "emiss12", istart, istride, iedge, fill_emiss[29], (void **) &rclr->sfc_emiss30);
  }
  if (chflg[30]) {
    read_seebor_hdf(id, "emiss14", istart, istride, iedge, fill_emiss[30], (void **) &rclr->sfc_emiss31);
  }
  if (chflg[31]) {
    read_seebor_hdf(id, "emiss15", istart, istride, iedge, fill_emiss[31], (void **) &rclr->sfc_emiss32);
  }
  if (chflg[32]) {
    read_seebor_hdf(id, "emiss16", istart, istride, iedge, fill_emiss[32], (void **) &rclr->sfc_emiss33);
  }
  if (chflg[33]) {
    read_seebor_hdf(id, "emiss16", istart, istride, iedge, fill_emiss[33], (void **) &rclr->sfc_emiss34);
  }
  if (chflg[34]) {
    read_seebor_hdf(id, "emiss16", istart, istride, iedge, fill_emiss[34], (void **) &rclr->sfc_emiss35);
  }
  if (chflg[35]) {
    read_seebor_hdf(id, "emiss16", istart, istride, iedge, fill_emiss[35], (void **) &rclr->sfc_emiss36);
  }

}

void seebor_channels_dateline(int32 id, int8 *chflg, int32 istart[2], int32 istride[2], int32 iedge[2],
                              int32 istart_2[2], int32 istride_2[2], int32 iedge_2[2], rtm_toa *rclr)
{

  char *rout = {"seebor_channels_dateline"};
  int count1, count2, npts;
  float *temp1, *temp2;

  read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[19], (void **) &temp1);
  read_seebor_hdf(id, "emiss7", istart_2, istride_2, iedge_2, fill_emiss[19], (void **) &temp2);
  free(temp1), free(temp2);
  count1 = iedge[0]*iedge[1];
  count2 = iedge_2[0]*iedge_2[1];
  npts = count1 + count2;
  rclr->npts = npts;
  printf("seebor_channels_dateline: rclr npts = %d\n", npts);

  if (chflg[19]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[19], (void **) &temp1);
    read_seebor_hdf(id, "emiss7", istart_2, istride_2, iedge_2, fill_emiss[19], (void **) &temp2);
    if ((rclr->sfc_emiss20 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss20");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss20);
    free(temp1), free(temp2);
  }
  if (chflg[20]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[20], (void **) &temp1);
    read_seebor_hdf(id, "emiss7", istart_2, istride_2, iedge_2, fill_emiss[20], (void **) &temp2);
    if ((rclr->sfc_emiss21 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss21");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss21);
    free(temp1), free(temp2);
  }
  if (chflg[21]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[21], (void **) &temp1);
    read_seebor_hdf(id, "emiss7", istart_2, istride_2, iedge_2, fill_emiss[21], (void **) &temp2);
    if ((rclr->sfc_emiss22 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss22");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss22);
    free(temp1), free(temp2);
  }
  if (chflg[22]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[22], (void **) &temp1);
    read_seebor_hdf(id, "emiss7", istart_2, istride_2, iedge_2, fill_emiss[22], (void **) &temp2);
    if ((rclr->sfc_emiss23 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss23");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss23);
    free(temp1), free(temp2);
  }
  if (chflg[23]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[23], (void **) &temp1);
    read_seebor_hdf(id, "emiss7", istart_2, istride_2, iedge_2, fill_emiss[23], (void **) &temp2);
    if ((rclr->sfc_emiss24 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss24");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss24);
    free(temp1), free(temp2);
  }
  if (chflg[24]) {
    read_seebor_hdf(id, "emiss7", istart, istride, iedge, fill_emiss[24], (void **) &temp1);
    read_seebor_hdf(id, "emiss7", istart_2, istride_2, iedge_2, fill_emiss[24], (void **) &temp2);
    if ((rclr->sfc_emiss25 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss25");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss25);
    free(temp1), free(temp2);
  }
  if (chflg[26]) {
    read_seebor_hdf(id, "emiss9", istart, istride, iedge, fill_emiss[26], (void **) &temp1);
    read_seebor_hdf(id, "emiss9", istart_2, istride_2, iedge_2, fill_emiss[26], (void **) &temp2);
    if ((rclr->sfc_emiss27 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss27");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss27);
    free(temp1), free(temp2);
  }
  if (chflg[27]) {
    read_seebor_hdf(id, "emiss10", istart, istride, iedge, fill_emiss[27], (void **) &temp1);
    read_seebor_hdf(id, "emiss10", istart_2, istride_2, iedge_2, fill_emiss[27], (void **) &temp2);
    if ((rclr->sfc_emiss28 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss28");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss28);
    free(temp1), free(temp2);
  }
  if (chflg[28]) {
    read_seebor_hdf(id, "emiss11", istart, istride, iedge, fill_emiss[28], (void **) &temp1);
    read_seebor_hdf(id, "emiss11", istart_2, istride_2, iedge_2, fill_emiss[28], (void **) &temp2);
    if ((rclr->sfc_emiss29 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss29");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss29);
    free(temp1), free(temp2);
  }
  if (chflg[29]) {
    read_seebor_hdf(id, "emiss12", istart, istride, iedge, fill_emiss[29], (void **) &temp1);
    read_seebor_hdf(id, "emiss12", istart_2, istride_2, iedge_2, fill_emiss[29], (void **) &temp2);
    if ((rclr->sfc_emiss30 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss30");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss30);
    free(temp1), free(temp2);
  }
  if (chflg[30]) {
    read_seebor_hdf(id, "emiss14", istart, istride, iedge, fill_emiss[30], (void **) &temp1);
    read_seebor_hdf(id, "emiss14", istart_2, istride_2, iedge_2, fill_emiss[30], (void **) &temp2);
    if ((rclr->sfc_emiss31 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss31");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss31);
    free(temp1), free(temp2);
  }
  if (chflg[31]) {
    read_seebor_hdf(id, "emiss15", istart, istride, iedge, fill_emiss[31], (void **) &temp1);
    read_seebor_hdf(id, "emiss15", istart_2, istride_2, iedge_2, fill_emiss[31], (void **) &temp2);
    if ((rclr->sfc_emiss32 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss32");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss32);
    free(temp1), free(temp2);
  }
  if (chflg[32]) {
    read_seebor_hdf(id, "emiss16", istart, istride, iedge, fill_emiss[32], (void **) &temp1);
    read_seebor_hdf(id, "emiss16", istart_2, istride_2, iedge_2, fill_emiss[32], (void **) &temp2);
    if ((rclr->sfc_emiss33 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss33");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss33);
    free(temp1), free(temp2);
  }
  if (chflg[33]) {
    read_seebor_hdf(id, "emiss16", istart, istride, iedge, fill_emiss[33], (void **) &temp1);
    read_seebor_hdf(id, "emiss16", istart_2, istride_2, iedge_2, fill_emiss[33], (void **) &temp2);
    if ((rclr->sfc_emiss34 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss34");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss34);
    free(temp1), free(temp2);
  }
  if (chflg[34]) {
    read_seebor_hdf(id, "emiss16", istart, istride, iedge, fill_emiss[34], (void **) &temp1);
    read_seebor_hdf(id, "emiss16", istart_2, istride_2, iedge_2, fill_emiss[34], (void **) &temp2);
    if ((rclr->sfc_emiss35 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss35");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss35);
    free(temp1), free(temp2);
  }
  if (chflg[35]) {
    read_seebor_hdf(id, "emiss16", istart, istride, iedge, fill_emiss[35], (void **) &temp1);
    read_seebor_hdf(id, "emiss16", istart_2, istride_2, iedge_2, fill_emiss[35], (void **) &temp2);
    if ((rclr->sfc_emiss36 = (float *) malloc(npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss36");
    cat_emiss(iedge[0], iedge[1], iedge_2[1], temp1, temp2, rclr->sfc_emiss36);
    free(temp1), free(temp2);
  }

}

void read_seebor_hdf(int32 id, char8 *sds_name, int32 istart[2], int32 istride[2], int32 iedge[2],
                     float fill_emiss, void **emiss)
{

  char *rout = {"read_seebor_hdf"};
  char8 sds_name_check[MAX_NC_NAME], att_name[MAX_NC_NAME];
  int32 sds_index, sd_id, status, data_type, rank, dimen[MAX_VAR_DIMS],
        nattr, att_index, count;
  void *attr, *bufp;
  int16 *buf_i16;
  float32 *scale_fac, *offset, *buffer;
  int npts, n;

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

  if (iedge[0] < 0) iedge[0] = dimen[0];
  if (iedge[1] < 0) iedge[1] = dimen[1];

  iedge[0] = min(dimen[0] - istart[0],iedge[0]);
  iedge[1] = min(dimen[1] - istart[1],iedge[1]);

  npts = iedge[0]*iedge[1];

  if ((bufp = (void *) malloc(npts*DFKNTsize(data_type))) == NULL)
    error_allo(rout,"bufp");

  if ((buffer = (float32 *) malloc(npts*sizeof(float32))) == NULL)
    error_allo(rout,"buffer");

  *emiss = buffer;

  status = SDreaddata(sd_id,istart,NULL,iedge,bufp);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot read specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }

  buf_i16 = (int16 *)bufp;

  strcpy(att_name,"scale_factor");
  att_index = SDfindattr(sd_id,"scale_factor");
  status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
  if (status == FAIL)
    error_exit(rout,"There was an error when calling SDattrinfo - scale_factor.");

  if ((attr = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
    error_allo(rout,att_name);
  status = SDreadattr(sd_id,att_index,attr);
  if (status == FAIL)
    error_exit(rout,"There was an error when calling SDreadattr - scale_factor.");

  scale_fac = (float32 *) attr;

  strcpy(att_name,"add_offset");
  att_index = SDfindattr(sd_id,att_name);
  status = SDattrinfo(sd_id,att_index,att_name,&data_type,&count);
  if (status == FAIL)
    error_exit(rout,"There was an error when calling SDattrinfo - add_offset.");

  if ((attr = (void *) malloc(count*DFKNTsize(data_type))) == NULL)
    error_allo(rout,att_name);
  status = SDreadattr(sd_id,att_index,attr);
  if (status == FAIL)
    error_exit(rout,"There was an error when calling SDreadattr - add_offset.");

  offset = (float32 *) attr;

  status = SDendaccess(sd_id);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot end access to specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }

  for (n=0; n<npts; n++) {
    buffer[n] = buf_i16[n]*scale_fac[0] + offset[0];
    if (buffer[n] < 0.0) buffer[n] = fill_emiss;
  }

  free(scale_fac);
  free(offset);
  free(buf_i16);

}

void cat_emiss(int nlat, int nlon1, int nlon2, float *temp1, float *temp2, float *emiss)
{

  int offset, ilat, ilon, index1, index2;

  for (ilat=0; ilat<nlat; ilat++) {
    offset = (nlon1 + nlon2) * ilat;
    for (ilon=0; ilon<nlon1; ilon++) {
      index1 = offset + ilon;
      index2 = (nlon1 * ilat) + ilon;
      emiss[index1] = temp1[index2];
    }
    for (ilon=0; ilon<nlon2; ilon++) {
      index1 = offset + nlon1 + ilon;
      index2 = (nlon2 * ilat) + ilon;
      emiss[index1] = temp2[index2];
    }
  }
}

void read_constant_emiss_imager(imagerL1 *imgr1, rtm_toa *rclr)
{

  char *rout = {"read_constant_emiss_imager"};
  int8 dateline_flg;
  int n, ilat1, ilat2, ilon1, ilon2, ilon1_2, ilon2_2;
  int32 istart[2], iedge[2],
        istart_2[2], iedge_2[2];

  if ((imgr1->index_rclr = (int *) malloc(imgr1->npts*sizeof(int))) == NULL)
    error_allo(rout,"imgr1->index_rclr");

  dateline_flg = 0;
  if (imgr1->bounds[1] > 180.0 || imgr1->bounds[3] > 180.0) dateline_flg = 1;

  if (dateline_flg == 0) {
    ilat1 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->bounds[2]-first_lat_sfcemiss)/dlat_sfcemiss)));
    ilat2 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->bounds[0]-first_lat_sfcemiss)/dlat_sfcemiss)));

    ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->bounds[1]-first_lon_sfcemiss)/dlon_sfcemiss)));
    ilon2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->bounds[3]-first_lon_sfcemiss)/dlon_sfcemiss)));

    if (ilat1 > ilat2) swap_int_values(&ilat1, &ilat2);
    if (ilon1 > ilon2) swap_int_values(&ilon1, &ilon2);

    istart[0] = ilat1;
    istart[1] = ilon1;

    iedge[0] = (ilat2-ilat1)+1;
    iedge[1] = (ilon2-ilon1)+1;

    rclr->npts = iedge[0]*iedge[1];

    constant_channels(imgr1->chflg, rclr);

    for (n=0; n<imgr1->npts; n++) {
      ilat1 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->lat[n]-first_lat_sfcemiss)/dlat_sfcemiss)));
      ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->lon[n]-first_lon_sfcemiss)/dlon_sfcemiss)));

      ilat2 = max(0,min((ilat1 - istart[0]),iedge[0]-1));
      ilon2 = max(0,min((ilon1 - istart[1]),iedge[1]-1));
      imgr1->index_rclr[n] = (iedge[1]*ilat2) + ilon2;
    }

  }
  else {

    ilat1 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->bounds[2]-first_lat_sfcemiss)/dlat_sfcemiss)));
    ilat2 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->bounds[0]-first_lat_sfcemiss)/dlat_sfcemiss)));

    ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->bounds[1]-first_lon_sfcemiss)/dlon_sfcemiss)));
    ilon2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(180.0-first_lon_sfcemiss)/dlon_sfcemiss)));

    ilon1_2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(-180.0-first_lon_sfcemiss)/dlon_sfcemiss)));
    ilon2_2 = max(0,min(nlon_sfcemiss-1, (int32) (fabs((imgr1->bounds[3]-360.0)-first_lon_sfcemiss)/dlon_sfcemiss)));

    if (ilat1 > ilat2) swap_int_values(&ilat1, &ilat2);
    if (ilon1 > ilon2) swap_int_values(&ilon1, &ilon2);
    if (ilon1_2 > ilon2_2) swap_int_values(&ilon1_2, &ilon2_2);

    istart[0] = ilat1;
    istart[1] = ilon1;

    iedge[0] = (ilat2-ilat1)+1;
    iedge[1] = (ilon2-ilon1)+1;

    istart_2[0] = ilat1;
    istart_2[1] = ilon1_2;

    iedge_2[0] = (ilat2-ilat1)+1;
    iedge_2[1] = (ilon2_2-ilon1_2)+1;

    rclr->npts = iedge[0]*iedge[1] + iedge_2[0]*iedge_2[1];

    constant_channels(imgr1->chflg, rclr);

    for (n=0; n<imgr1->npts; n++) {
      ilat1 = max(0,min(nlat_sfcemiss-1, (int32) (fabs(imgr1->lat[n]-first_lat_sfcemiss)/dlat_sfcemiss)));
      ilon1 = max(0,min(nlon_sfcemiss-1, (int32) (fabs(imgr1->lon[n]-first_lon_sfcemiss)/dlon_sfcemiss)));

      if (imgr1->lon[n] >= 0.0) {
        ilat2 = max(0,min((ilat1 - istart[0]),iedge[0]-1));
        ilon2 = max(0,min((ilon1 - istart[1]),iedge[1]-1));
        imgr1->index_rclr[n] = iedge_2[1]*ilat2 + ((iedge[1]*ilat2) + ilon2);
      }
      else {
        ilat2 = max(0,min((ilat1 - istart_2[0]),iedge_2[0]-1));
        ilon2 = max(0,min((ilon1 - istart_2[1]),iedge_2[1]-1));
        imgr1->index_rclr[n] = iedge[1]*ilat2 + ((iedge_2[1]*ilat2) + ilon2);
      }
    }
  }

  if ((rclr->flag = (unsigned char *) malloc(rclr->npts*sizeof(unsigned char))) == NULL)
    error_allo(rout,"rclr->flag");
  for (n=0; n<rclr->npts; n++) rclr->flag[n] = NO;

}

void constant_channels(int8 *chflg, rtm_toa *rclr)
{

  char *rout = {"constant_channels"};

  if (chflg[19]) {
    if ((rclr->sfc_emiss20 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss20");
    set_emiss_arr(rclr->sfc_emiss20, fill_emiss[19], rclr->npts);
  }
  if (chflg[20]) {
    if ((rclr->sfc_emiss21 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss21");
    set_emiss_arr(rclr->sfc_emiss21, fill_emiss[20], rclr->npts);
  }
  if (chflg[21]) {
    if ((rclr->sfc_emiss22 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss22");
    set_emiss_arr(rclr->sfc_emiss22, fill_emiss[21], rclr->npts);
  }
  if (chflg[22]) {
    if ((rclr->sfc_emiss23 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss23");
    set_emiss_arr(rclr->sfc_emiss23, fill_emiss[22], rclr->npts);
  }
  if (chflg[23]) {
    if ((rclr->sfc_emiss24 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss24");
    set_emiss_arr(rclr->sfc_emiss24, fill_emiss[23], rclr->npts);
  }
  if (chflg[24]) {
    if ((rclr->sfc_emiss25 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss25");
    set_emiss_arr(rclr->sfc_emiss25, fill_emiss[24], rclr->npts);
  }
  if (chflg[26]) {
    if ((rclr->sfc_emiss27 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss27");
    set_emiss_arr(rclr->sfc_emiss27, fill_emiss[26], rclr->npts);
  }
  if (chflg[27]) {
    if ((rclr->sfc_emiss28 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss28");
    set_emiss_arr(rclr->sfc_emiss28, fill_emiss[27], rclr->npts);
  }
  if (chflg[28]) {
    if ((rclr->sfc_emiss29 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss29");
    set_emiss_arr(rclr->sfc_emiss29, fill_emiss[28], rclr->npts);
  }
  if (chflg[29]) {
    if ((rclr->sfc_emiss30 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss30");
    set_emiss_arr(rclr->sfc_emiss30, fill_emiss[29], rclr->npts);
  }
  if (chflg[30]) {
    if ((rclr->sfc_emiss31 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss31");
    set_emiss_arr(rclr->sfc_emiss31, fill_emiss[30], rclr->npts);
  }
  if (chflg[31]) {
    if ((rclr->sfc_emiss32 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss32");
    set_emiss_arr(rclr->sfc_emiss32, fill_emiss[31], rclr->npts);
  }
  if (chflg[32]) {
    if ((rclr->sfc_emiss33 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss33");
    set_emiss_arr(rclr->sfc_emiss33, fill_emiss[32], rclr->npts);
  }
  if (chflg[33]) {
    if ((rclr->sfc_emiss34 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss34");
    set_emiss_arr(rclr->sfc_emiss34, fill_emiss[33], rclr->npts);
  }
  if (chflg[34]) {
    if ((rclr->sfc_emiss35 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss35");
    set_emiss_arr(rclr->sfc_emiss35, fill_emiss[34], rclr->npts);
  }
  if (chflg[35]) {
    if ((rclr->sfc_emiss36 = (float *) malloc(rclr->npts*sizeof(float))) == NULL)
      error_allo(rout,"rclr->sfc_emiss36");
    set_emiss_arr(rclr->sfc_emiss36, fill_emiss[35], rclr->npts);
  }

}

void set_emiss_arr(float *emiss, float const_emiss, long npts)
{

  int n;

  for (n=0; n<npts; n++) {
    emiss[n] = const_emiss;
  }

}
