/*$Id: utils.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"
#include "imagerL1_leocat.h"
#include "rtm_leocat.h"
#include "nwp_leocat.h"
#include "utils.h"

unsigned char ** allocate_2d_uchar8_ptr(char *var, int nrow, int ncol)

{

  char r_name[256];
  int irow=0;
  unsigned char **ptr;

  strcpy(r_name,"allocate_2d_uchar8_ptr");

  if ((ptr = (unsigned char **) malloc(sizeof(unsigned char *) * nrow)) == NULL)
    error_allo(r_name,var);

  for (irow=0; irow < nrow; irow++) {
    if ((ptr[irow] = (unsigned char *) malloc(sizeof(unsigned char) * ncol)) == NULL)
     error_allo(r_name,var);
  }

  return(ptr);

}

float ** allocate_2d_float_ptr(char *var, int nrow, int ncol)

{

  char r_name[256];
  int irow=0;
  float **ptr;

  strcpy(r_name,"allocate_2d_float_ptr");

  if ((ptr = (float **) malloc(sizeof(float *) * nrow)) == NULL)
    error_allo(r_name,var);

  for (irow=0; irow < nrow; irow++) {
    if ((ptr[irow] = (float *) malloc(sizeof(float) * ncol)) == NULL)
     error_allo(r_name,var);
  }

  return(ptr);

}

double ** allocate_2d_double_ptr(char *var, int nrow, int ncol)

{

  char r_name[256];
  int irow=0;
  double **ptr;

  strcpy(r_name,"allocate_2d_double_ptr");

  if ((ptr = (double **) malloc(sizeof(double *) * nrow)) == NULL)
    error_allo(r_name,var);

  for (irow=0; irow < nrow; irow++) {
    if ((ptr[irow] = (double *) malloc(sizeof(double) * ncol)) == NULL)
     error_allo(r_name,var);
  }

  return(ptr);

}

float *** allocate_3d_float_ptr(char *var, int nver, int nrow, int ncol)

{

  char r_name[256];
  int iver=0,irow=0;
  float ***ptr;

  strcpy(r_name,"allocate_3d_float_ptr");

  if ((ptr = (float ***) malloc(sizeof(float **) * nver)) == NULL)
    error_allo(r_name,var);

  for (iver=0; iver < nver; iver++) {
    if ((ptr[iver] = (float **) malloc(sizeof(float *) * nrow)) == NULL)
     error_allo(r_name,var);
    for (irow=0; irow < nrow; irow++) {
      if ((ptr[iver][irow] = (float *) malloc(sizeof(float) * ncol)) == NULL)
       error_allo(r_name,var);
    }
  }

  return(ptr);

}

unsigned char ** calloc_2d_uchar8_ptr(char *var, int nrow, int ncol)

{

  char r_name[256];
  int irow=0;
  unsigned char **ptr;

  strcpy(r_name,"calloc_2d_uchar8_ptr");

  if ((ptr = (unsigned char **) calloc(nrow,sizeof(unsigned char *))) == NULL)
    error_allo(r_name,var);

  for (irow=0; irow < nrow; irow++) {
    if ((ptr[irow] = (unsigned char *) calloc(ncol,sizeof(unsigned char))) == NULL)
     error_allo(r_name,var);
  }

  return(ptr);

}

float ** calloc_2d_float_ptr(char *var, int nrow, int ncol)

{

  char r_name[256];
  int irow=0;
  float **ptr;

  strcpy(r_name,"calloc_2d_float_ptr");

  if ((ptr = (float **) calloc(nrow,sizeof(float *))) == NULL)
    error_allo(r_name,var);

  for (irow=0; irow < nrow; irow++) {
    if ((ptr[irow] = (float *) calloc(ncol,sizeof(float))) == NULL)
     error_allo(r_name,var);
  }

  return(ptr);

}

double ** calloc_2d_double_ptr(char *var, int nrow, int ncol)

{

  char r_name[256];
  int irow=0;
  double **ptr;

  strcpy(r_name,"calloc_2d_double_ptr");

  if ((ptr = (double **) calloc(nrow,sizeof(double *))) == NULL)
    error_allo(r_name,var);

  for (irow=0; irow < nrow; irow++) {
    if ((ptr[irow] = (double *) calloc(ncol,sizeof(double))) == NULL)
     error_allo(r_name,var);
  }

  return(ptr);

}

void destroy_2d_uchar8_ptr(int nrow, unsigned char **ptr)

{

  int irow=0;

  for (irow=0; irow < nrow; irow++) free(ptr[irow]);
  free(ptr);

}

void destroy_2d_float_ptr(int nrow, float **ptr)

{

  int irow=0;

  for (irow=0; irow < nrow; irow++) free(ptr[irow]);
  free(ptr);

}

void destroy_2d_double_ptr(int nrow, double **ptr)

{

  int irow=0;

  for (irow=0; irow < nrow; irow++) free(ptr[irow]);
  free(ptr);

}

int leap_year_check (int year)
{
  int ileap = 0;

  if (((year % 4) == 0 && (year % 100) != 0) ||
      (year % 400) == 0) ileap = 1;

  return (ileap);
}

int jday2month (int jday, int ileap)

{
  int month;

  if (jday < 32)
    month = 1;
  else if (jday < 60+ileap)
    month = 2;
  else if (jday < 91+ileap)
    month = 3;
  else if (jday < 121+ileap)
    month = 4;
  else if (jday < 152+ileap)
    month = 5;
  else if (jday < 182+ileap)
    month = 6;
  else if (jday < 213+ileap)
    month = 7;
  else if (jday < 244+ileap)
    month = 8;
  else if (jday < 274+ileap)
    month = 9;
  else if (jday < 305+ileap)
    month = 10;
  else if (jday < 335+ileap)
    month = 11;
  else
    month = 12;

  return(month);
}

int jday2day (int jday, int ileap)

{
  int day;

  if (jday < 32)
    day = jday;
  else if (jday < 60)
    day = jday - 31;
  else if ((jday == 60) && (ileap == 1))
    day = jday - 31;
  else if ((jday == 60) && (ileap == 0))
    day = jday - 59;
  else if (jday < 91+ileap)
    day = jday - (59 + ileap);
  else if (jday < 121+ileap)
    day = jday - (90 + ileap);
  else if (jday < 152+ileap)
    day = jday - (120 + ileap);
  else if (jday < 182+ileap)
    day = jday - (151 + ileap);
  else if (jday < 213+ileap)
    day = jday - (181 + ileap);
  else if (jday < 244+ileap)
    day = jday - (212 + ileap);
  else if (jday < 274+ileap)
    day = jday - (243 + ileap);
  else if (jday < 305+ileap)
    day = jday - (273 + ileap);
  else if (jday < 335+ileap)
    day = jday - (304 + ileap);
  else
    day = jday - (334 + ileap);

  return(day);
}

/*******************************************************************************
/ Returns saturation vapor pressure over water, in hPa, given temperature in K.
/
*******************************************************************************/
double svpwat(double t) {

     double a0 =  .999996876e0;
     double a1 = -.9082695004e-2;
     double a2 =  .7873616869e-4;
     double a3 = -.6111795727e-6;
     double a4 =  .4388418740e-8;
     double a5 = -.2988388486e-10;
     double a6 =  .2187442495e-12;
     double a7 = -.1789232111e-14;
     double a8 =  .1111201803e-16;
     double a9 = -.3099457145e-19;

     double b = .61078e+1;

     t -= 273.16;

     return b / pow(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*a9)))))))), 8);
}



/*******************************************************************************
/ Returns saturation vapor pressure over ice, in hPa, given temperature in K.
/ The Goff-Gratch equation (Smithsonian Met. Tables,  5th ed., pp. 350, 1984)
*******************************************************************************/
double svpice(double t) {

     double a;

     a = 273.16 / t;

     return pow(10, -9.09718 * (a - 1.) -
                     3.56654 * log10(a) +
                     0.876793 * (1. - 1./a) +
                     log10(6.1071));
}

/*******************************************************************************
/ Returns saturation mixing ratio over water, in g/kg, given pressure in hPa and
/ temperature in K.
*******************************************************************************/
double satmixwat(double p, double t) {

     double es;

     es = svpwat(t);
/*
     return 622. * es / (p - es);
*/
     return 622. * es / p;
}



/*******************************************************************************
/ Returns saturation mixing ratio over ice, in g/kg, given pressure in hPa and
/ temperature in K.
*******************************************************************************/
double satmixice(double p, double t) {

     double es;

     es = svpice(t);
/*
     return 622. * es / (p - es);
*/
     return 622. * es / p;
}




/*******************************************************************************
/ Returns saturation mixing ratio in g/kg, given pressure in hPa and
/ temperature in K.
*******************************************************************************/
double satmix(double p, double t) {

     if (t > 253.)
          return satmixwat(p, t);
     else
          return satmixice(p, t);
}



/*******************************************************************************
/
*******************************************************************************/
double rh_to_mr_wat(double rh, double p, double t) {

     return rh * 0.01 * satmixwat(p, t);
}



double mr_to_rh_wat(double mr, double p, double t) {

     return mr * 100. / satmixwat(p, t);
}

void array_minmax_double (long npts, double *array, double *minval, double *maxval)

{

 long
  n;

  *minval = 9999.0;
  *maxval = -9999.0;
  for (n=0; n<npts; n++) {
    if (array[n] != MISSING_FLOAT) {
      *minval = min(array[n],*minval);
      *maxval = max(array[n],*maxval);
    }
  }

}

void array_minmax_float (long npts, float *array, float *minval, float *maxval)

{

 long
  n;

  *minval = 9999.0;
  *maxval = -9999.0;
  for (n=0; n<npts; n++) {
    if (array[n] != MISSING_FLOAT) {
      *minval = min(array[n],*minval);
      *maxval = max(array[n],*maxval);
    }
  }

}

void array_minmax_int (long npts, int *array, int *minval, int *maxval)

{

 long
  n;

  *minval = 9999;
  *maxval = -9999;
  for (n=0; n<npts; n++) {
    if (array[n] != MISSING_INT) {
      *minval = min(array[n],*minval);
      *maxval = max(array[n],*maxval);
    }
  }

}

void array_minmax_sub_float (int i1, int i2, int j1, int j2, int ncol, float *array, float *minval, float *maxval)

{

 long
  i, j, index;

  *minval = 9999.0;
  *maxval = -9999.0;
  for (j=j1; j<=j2; j++) {
    for (i=i1; i<=i2; i++) {
      index = (j*ncol) + i;
      if (array[index] != MISSING_FLOAT) {
        *minval = min(array[index],*minval);
        *maxval = max(array[index],*maxval);
      }
    }
  }
}

void array_minmax_sub_double (int i1, int i2, int j1, int j2, int ncol, double *array, double *minval, double *maxval)

{

 long
  i, j, index;

  *minval = 9999.0;
  *maxval = -9999.0;
  for (j=j1; j<=j2; j++) {
    for (i=i1; i<=i2; i++) {
      index = (j*ncol) + i;
      if (array[index] != MISSING_FLOAT) {
        *minval = min(array[index],*minval);
        *maxval = max(array[index],*maxval);
      }
    }
  }
}

int locate(float *xx, int n, float x)
{
  int j;
  int i,jl,jm,ju;
  jl = 0;
  ju = n + 1;
  for (i=0; i<(2*n); i++) {
   if (ju-jl <= 1)
     break;
   jm = (ju + jl)/2;
   if ((xx[n-1] >= xx[0]) == (x >= xx[jm-1]))
     jl=jm;
   else
     ju=jm;
  }

  if (x == xx[0])
    j = 1;
  else if (x == xx[n-1])
    j = n - 1;
  else
    j = jl;

  return(j-1);
}

/*******************************************************************************
/
*******************************************************************************/
double lin_int(double x1, double x2, double x, double y1, double y2) {

  return (x - x1) / (x2 - x1) * (y2 - y1);
}

/*******************************************************************************
/
*******************************************************************************/
int lin_interp_sorted(double *x, double xx,
                      double *y, int n, int *i0, double *f0, double *yy) {

     int i;
     int ii;

     int last;

     last = n - 1;

     if (x[0] < x[1]) {
          if (xx <= x[0]) {
               *i0 = 0;
               *f0 = 0.;
               *yy = y[0];
               return -1;
          }
          else
          if (xx >= x[last]) {
               *i0 = last;
               *f0 = 0.;
               *yy = y[last];
               return  1;
          }
          else {
               for (i = 1; 1; ++i) {
                    if (xx < x[i]) {
                         ii = i-1;
/*
                         *yy = y[ii] + lin_int(x[ii], x[i], xx, y[ii], y[i]);
*/
                         *i0 = ii;
                         *f0 = (xx - x[ii]) / (x[i] - x[ii]);
                         *yy = y[ii] + *f0 * (y[i] - y[ii]);
                         return 0;
                    }
               }
          }
     }
     else {
          if (xx >= x[0]) {
               *i0 = 0;
               *f0 = 0.;
               *yy = y[0];
               return -1;
          }
          else
          if (xx <= x[last]) {
               *i0 = last;
               *f0 = 0.;
               *yy = y[last];
               return  1;
          }
          else {
               for (i = 1; 1; ++i) {
                    if (xx > x[i]) {
                         ii = i-1;
/*
                         *yy = y[ii] + lin_int(x[ii], x[i], xx, y[ii], y[i]);
*/
                         *i0 = ii;
                         *f0 = (xx - x[ii]) / (x[i] - x[ii]);
                         *yy = y[ii] + *f0 * (y[i] - y[ii]);
                         return 0;
                    }
               }
          }
     }
}

void time_elapsed(double total_sec, int *dd, int *hh, int *mm, int *ss)

{
  *dd = (int)(total_sec/86400);
  *hh = (int)(((long)total_sec % 86400)/3600);
  *mm = (int)(((long)total_sec % 3600)/60);
  *ss = (int)((long)total_sec % 60);
}

void print_time_elapsed(int dd, int hh, int mm, int ss, int option)

{

  if (option == 0) {
    fprintf(stdout,"%sTotal processing time:"
            "dd=%d, hh=%d, mm=%d, ss=%d\n",
	    EXE_PROMPT,dd,hh,mm,ss);
  }
  else {
    fprintf(stdout,"%sTotal processing time:\n"
            "\t Days = %d\n"
	    "\t Hours = %d\n"
	    "\t Minutes = %d\n"
	    "\t Seconds = %d\n",EXE_PROMPT,
            dd,hh,mm,ss);
  }

}

int eqangle_index(float lat, float lon, float dlat, float dlon,
                  float first_lat, float first_lon,
		  float last_lat, float last_lon, int *j, int *i)
{

  int factor = 1;
  int nlon;

  float x, x0, dx;
  float y, y0, dy;

  y = min(max(lat,  -89.99),  89.99);
  y0 = first_lat;
  dy = dlat / factor;
  if (first_lat > last_lat) dy *= (-1.0);
  *j = (int) ((y - y0 + 0.5*dy) / dy);

  x = min(max(lon, 0.0), 359.99);
  x0 = first_lon;
  dx = dlon / factor;
  *i = (int) ((x - x0 + 0.5*dx) / dx);
  if (*i == 360)
    *i = 0;

  nlon = 360.0/dx;

  return((nlon*(*j)) + *i);

}

void i2xy(int index, int ncols, int *x, int *y)
{

  *x = index % ncols;
  *y = index/ncols;

}

float convert_lon(int reverse, float lon)

{

  if (reverse == 0) {
    if (lon < 0.0) lon += 360.0;
  }
  else {
    if (lon > 180.0) lon -= 360.0;
  }

  return(lon);
}

void swap_int_values(int *val1, int *val2)
{
  int temp;

  temp = *val1;
  *val1 = *val2;
  *val2 = temp;
}

void eqangle_index2latlon(float first_lon, float dlon, float first_lat, float dlat,
                          int x, int y, float *lon, float *lat)

{

  char *r_name = {"eqangle_index2latlon"};

  *lat = first_lat + (dlat*y);
  if (*lat > 90.0) *lat = 180.0 - (*lat);

  *lon = first_lon + (dlon*x);
  if (*lon >= 360.0) *lon -= 360.0;

  if (*lat < -90.0 || *lat > 90.0)
    error_exit(r_name,"Equal angle grid latitudes are out of bounds-cannot recover.");
  if (*lon < 0.0 || *lon > 360.0)
    error_exit(r_name,"Equal angle longitudes are out of bounds-cannot recover.");

}

/******************************************************************************/
/******************************************************************************/

void read_umd_lmask_map(unsigned char *lmask_map,char *sfc_dir_name,char *sub_dir_name)
{
  FILE *fp;

  char fileWithPath[200];
  sprintf(fileWithPath,"%s/%s/%s",sfc_dir_name,sub_dir_name,lmask_file_name);
//printf("lmask_file_name = %s\n",lmask_file_name);
  if ((fp = fopen(fileWithPath,"r")) == NULL) {
    fprintf(stderr,"%sCannot open the land mask file %s\nUSER:check file location and format - aborting\n",EXE_PROMPT,fileWithPath);
    exit(EXIT_FAILURE);
  }
  fread(lmask_map,sizeof(unsigned char),nlon_lmask*nlat_lmask,fp);
  fclose(fp);

}

/******************************************************************************/
/******************************************************************************/

void read_IGBP_map(unsigned char *lmask_map, char *sfc_dir_name, char *sub_dir_name)
{

   char fileWithPath[200];
   int IGBP_fileID;
   int h4status;
   int IGBP_dsetIDx;
   int IGBP_dsetID;
   char * IGBP_dsetName = "IGBP_Land_Cover_Type";
   int start[2];
   int edge[2];

   sprintf(fileWithPath,"%s/%s",sfc_dir_name,IGBP_file_name);
   printf("\nIGBP_file_name = %s\n", IGBP_file_name);
   printf("fileWithPath = %s\n", fileWithPath);

// Open the file and get the file id
   IGBP_fileID = -1;
   IGBP_fileID = SDstart(fileWithPath, DFACC_READ);
   if (IGBP_fileID == FAIL) {
     fprintf(stderr,"%sInvalid HDF file, %s\n",EXE_PROMPT,fileWithPath);
     exit (EXIT_FAILURE);
   }
// printf("IGBP_fileID = %d\n",IGBP_fileID);

// Get the index of the dataset in the file from the dataset name.
   IGBP_dsetIDx = -1;
   if (IGBP_fileID != -1){
     IGBP_dsetIDx = SDnametoindex(IGBP_fileID,IGBP_dsetName);
//   printf("IGBP_dsetIDx = %d\n",IGBP_dsetIDx);
   }

// Connect to dataset and get the dataset id
   IGBP_dsetID = -1;
   if (IGBP_dsetIDx != -1){
//   printf("Connecting to dataset in IGBP file %s\n",fileWithPath);
     IGBP_dsetID = SDselect(IGBP_fileID,IGBP_dsetIDx);
     if (IGBP_dsetID == FAIL) {
       fprintf(stderr,"%s-cannot select specified SDS, =%s\n",EXE_PROMPT,IGBP_dsetName);
       exit(EXIT_FAILURE);
     }
//   printf("IGBP_dsetID = %d\n",IGBP_dsetID);
   }

// Read the IGBP dataset
// printf("\nDataset ID for file is %s is %d\n",fileWithPath,IGBP_dsetID);
   start[0] = 0;
   start[1] = 0;
   edge[0] = nlat_IGBP;
   edge[1] = nlon_IGBP;
   h4status = SDreaddata(IGBP_dsetID,start,NULL,edge, (VOIDP) lmask_map);
   if (h4status == FAIL) {
     fprintf(stderr,"%s-cannot read specified SDS %s, freeing memory.\n",EXE_PROMPT,IGBP_dsetName);
     free(lmask_map);
     exit(EXIT_FAILURE);
   }

// End access to dataset
// printf("Ending access to dataset in IGBP file %s\n",fileWithPath);
   h4status = SDendaccess(IGBP_dsetID);
   if (h4status == FAIL) {
     fprintf(stderr,"%s-cannot end access to specified SDS, =%s\n",EXE_PROMPT,IGBP_dsetName);
     free(lmask_map);
     exit(EXIT_FAILURE);
   }

// Close file
// printf("Closing IGBP file %s\n",fileWithPath);
   h4status = SDend(IGBP_fileID);
   if (h4status == FAIL) {
     fprintf(stderr,"%s-cannot close %s\n",EXE_PROMPT,fileWithPath);
     free(lmask_map); 
     exit(EXIT_FAILURE); 
   }

}

/******************************************************************************/
/******************************************************************************/

void read_elev_map(int16 *elev_map,char *sfc_dir_name,char *sub_dir_name)
{
  FILE *fp;
  char fileWithPath[200];
  sprintf(fileWithPath,"%s/%s/%s",sfc_dir_name,sub_dir_name,elev_file_name);
//printf("elev_file_name = %s\n",elev_file_name);
  if ((fp = fopen(fileWithPath,"r")) == NULL) {
    fprintf(stderr,"%sCannot open the surface elevation file %s\nUSER:check file location and format - aborting\n",EXE_PROMPT,fileWithPath);
    exit(EXIT_FAILURE);
  }
  fread(elev_map,sizeof(int16),nlon_elev*nlat_elev,fp);
  fclose(fp);

}

/******************************************************************************/
/******************************************************************************/

void read_sst_map(float *sst_map,char *sfc_dir_name,char *sub_dir_name)
{
  char *rout = {"read_sst_map"};
  FILE *fp;
  int8 *temp;
  int n, count;
  char fileWithPath[200];
  sprintf(fileWithPath,"%s/%s/%s",sfc_dir_name,sub_dir_name,sst_file_name);
//printf("sst_file_name = %s\n",sst_file_name);

  count = nlon_sst * nlat_sst * 12;

  if ((temp = (int8 *) malloc(sizeof(int8) * count)) == NULL)
    error_allo(rout,"temp");

  if ((fp = fopen(fileWithPath,"r")) == NULL) {
    fprintf(stderr,"%sCannot open the climatological sst file %s\nUSER:check file location and format - aborting\n",EXE_PROMPT,fileWithPath);
    exit(EXIT_FAILURE);
  }
  fread(temp,sizeof(int8),count,fp);
  fclose(fp);

  for (n=0; n<count; n++) {
    sst_map[n] = min_sst_clim + (max_sst_clim - min_sst_clim)*(temp[n] + 128)/255.0;
    if (sst_map[n] == min_sst_clim) {
      sst_map[n] = MISSING_FLOAT;
    }
    else {
      sst_map[n] = sst_map[n] + 273.15;
    }
  }

  free(temp);

}

/******************************************************************************/
/******************************************************************************/

unsigned char * umd_lmask_grid2pixel(float * plat, float * plon, long npts, unsigned char *lmask_map)
{
  char *rout = {"umd_lmask_grid2pixel"};
  unsigned char *pixel_stype;
  long i, ilat, ilon, index;
  float dlon_lmask, dlat_lmask;

  if ((pixel_stype = (unsigned char *) malloc(npts*sizeof(unsigned char))) == NULL)
    error_allo(rout,"pixel_stype");

  dlat_lmask = fabs(last_lat_lmask - first_lat_lmask)/nlat_lmask;
  dlon_lmask = fabs(last_lon_lmask - first_lon_lmask)/nlon_lmask;

  for (i=0; i<npts; i++) {
    ilat = max(0,min(nlat_lmask-1, (long) (fabs(plat[i]-first_lat_lmask)/dlat_lmask)));
    if (plon[i] < first_lon_lmask) {
      ilon = max(0,min(nlon_lmask-1, (long) ((360.0 - fabs(plon[i]-first_lon_lmask))/dlon_lmask)));
    }
    else {
      ilon = max(0,min(nlon_lmask-1, (long) ((plon[i]-first_lon_lmask)/dlon_lmask)));
    }
    index = (nlon_lmask*ilat) + ilon;
    pixel_stype[i] = lmask_map[index];
  }

  return(pixel_stype);

}

/******************************************************************************/
/******************************************************************************/

unsigned char * IGBP_grid2pixel(float *plat, float *plon, long npts, unsigned char *lmask_map)
{
  char *rout = {"IGBP_grid2pixel"};
  unsigned char *pixel_stype;
  long i, ilat, ilon, index;
  float dlon_IGBP, dlat_IGBP;
  
  if ((pixel_stype = (unsigned char *) malloc(npts*sizeof(unsigned char))) == NULL)
    error_allo(rout,"pixel_stype");
  
  dlat_IGBP = fabs(last_lat_IGBP - first_lat_IGBP)/nlat_IGBP;
  dlon_IGBP = fabs(last_lon_IGBP - first_lon_IGBP)/nlon_IGBP;
    
  for (i=0; i<npts; i++) {
    ilat = max(0,min(nlat_IGBP-1, (long) (fabs(plat[i]-first_lat_IGBP)/dlat_IGBP)));
    if (plon[i] < first_lon_IGBP) {
      ilon = max(0,min(nlon_IGBP-1, (long) ((360.0 - fabs(plon[i]-first_lon_IGBP))/dlon_IGBP)));
    }
    else {
      ilon = max(0,min(nlon_IGBP-1, (long) ((plon[i]-first_lon_IGBP)/dlon_IGBP)));
    }
    index = (nlon_IGBP*ilat) + ilon;
    pixel_stype[i] = lmask_map[index];
  }
  
  return(pixel_stype);
  
} 

/******************************************************************************/
/******************************************************************************/

float * elev_grid2pixel(float * plat, float * plon, long npts, int16 *elev_map)
{
  char *rout = {"elev_grid2pixel"};
  float *pixel_elev;
  long i, ilat, ilon, index;
  float dlon_elev, dlat_elev, xlon;

  if ((pixel_elev = (float *) malloc(npts*sizeof(float))) == NULL)
    error_allo(rout,"pixel_elev");

  dlat_elev = fabs(last_lat_elev - first_lat_elev)/nlat_elev;
  dlon_elev = fabs(last_lon_elev - first_lon_elev)/nlon_elev;

  for (i=0; i<npts; i++) {
    if (plon[i] < 0.0) {
      xlon = plon[i] + 360.0;
    }
    else {
      xlon = plon[i];
    }

    ilat = max(0,min(nlat_elev-1, (long) (fabs(plat[i]-first_lat_elev)/dlat_elev)));
    ilon = max(0,min(nlon_elev-1, (long) ((xlon-first_lon_elev)/dlon_elev)));

    index = (nlon_elev*ilat) + ilon;
    pixel_elev[i] = elev_map[index];
  }

  return(pixel_elev);

}

/******************************************************************************/
/******************************************************************************/

float * sst_grid2pixel(float * plat, float * plon, long npts, int month, float *sst_map)
{
  char *rout = {"sst_grid2pixel"};
  float *pixel_sst;
  long i, ilat, ilon, index;
  float dlon_sst, dlat_sst, xlon;

  if ((pixel_sst = (float *) malloc(npts*sizeof(float))) == NULL)
    error_allo(rout,"pixel_sst");

  dlat_sst = fabs(last_lat_sst - first_lat_sst)/nlat_sst;
  dlon_sst = fabs(last_lon_sst - first_lon_sst)/nlon_sst;

  for (i=0; i<npts; i++) {
    if (plon[i] < 0.0) {
      xlon = plon[i] + 360.0;
    }
    else {
      xlon = plon[i];
    }

    ilat = max(0,min(nlat_sst-1, (long) (fabs(plat[i]-first_lat_sst)/dlat_sst)));
    ilon = max(0,min(nlon_sst-1, (long) ((xlon-first_lon_sst)/dlon_sst)));

    index = ((nlon_sst*ilat) + ilon) + ((month - 1)*nlon_sst*nlat_sst);
    pixel_sst[i] = max(270.0,sst_map[index]);
  }

  return(pixel_sst);

}

/******************************************************************************/
/******************************************************************************/

void nwp_grid2pixel(nwp_params *nwp, imagerL1 *imgr1)
{

//      printf("Inside nwp_grid2pixel()...\n");

	char *rout = {"nwp_grid2pixel"};
	int i, j;
	int factor = 1;
	long int n;

	float x, x0, dx;
	float y, y0, dy;

	//RAF Mods
	float lat1, lat2, lon1, lon2, lat00, lat10, lon00, lon10, dyi, pyi, dxi, pxi;
	profile_params *nwp2;
	nwp2 = nwp->map;
	int k, i1_temp, i2_temp;

//      printf("Allocating memory for imgr1 quantities...\n");

	if ((imgr1->index_nwp = (int *) malloc(imgr1->npts*sizeof(int))) == NULL)
		error_allo(rout,"imgr1->index_nwp");
	if ((imgr1->index_vza = (int *) malloc(imgr1->npts*sizeof(int))) == NULL)
		error_allo(rout,"imgr1->index_vza");
	if ((imgr1->j1 = (int *) malloc(imgr1->npts*sizeof(int))) == NULL)
		error_allo(rout,"imgr1->j1");
	if ((imgr1->j2 = (int *) malloc(imgr1->npts*sizeof(int))) == NULL)
		error_allo(rout,"imgr1->j2");
	if ((imgr1->i1 = (int *) malloc(imgr1->npts*sizeof(int))) == NULL)
		error_allo(rout,"imgr1->i1");
	if ((imgr1->i2 = (int *) malloc(imgr1->npts*sizeof(int))) == NULL)
		error_allo(rout,"imgr1->i2");
	if ((imgr1->py = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
		error_allo(rout,"imgr1->py");
	if ((imgr1->px = (float *) malloc(imgr1->npts*sizeof(float))) == NULL)
		 error_allo(rout,"imgr1->px");

//      printf("Finished allocating memory for imgr1 quantities...\n");

	for (n=0; n<imgr1->npts; n++) {

          if( (imgr1->lat[n] >= -90.0 && imgr1->lat[n] <= 90.0) &&
              (imgr1->lon[n] >= -180.0 && imgr1->lon[n] <= 180.0) ) {

		y = min(max(imgr1->lat[n],  -89.99),  89.99);
		y0 = nwp->first_lat;
		dy = nwp->dlat / factor;
		if (nwp->first_lat > nwp->last_lat) dy *= (-1.0);
		j = (int) ((y - y0 + 0.5*dy) / dy);

		x = min(max(imgr1->lon[n], -179.99), 179.99);
		//  RAF added test for bounds = -180
		if (imgr1->bounds[1] > 179.99 || imgr1->bounds[3] > 179.99 ||
			imgr1->bounds[1] <= -179.99 || imgr1->bounds[3] <= -179.99 ||
			(nwp->first_lon == 0.0 && nwp->last_lon == 359.0) ) {
		  if (x < 0.0) x += 360.0;
		}

		x0 = nwp->first_lon;
		//  RAF added check for special condition.
		if(imgr1->bounds[1] <= 0.0 && imgr1->bounds[3] <= 0.0 && x0 == 180.0) x0 = -180.0;
		if (nwp->first_lon < 180.0 && nwp->last_lon > 180.0 && x < 0.0) x += 360.0;
		dx = nwp->dlon / factor;
		i = (int) ((x - x0 + 0.5*dx) / dx);
		if (i == 360)
			i = 0;

		//  Added by RAF for spatial interpolation in leocatBridge

		k = min(max( (nwp->nlon * j) + i, 0),nwp->ncells - 1);

		lat1 = y0 + ((j-1) * dy);
		lat2 = y0 + ((j+1) * dy);
		if( fabs(y-lat1) <= fabs(y-lat2) ) {
			imgr1->j1[n] = j-1;
			imgr1->j2[n] = j;
		}else {
			imgr1->j1[n] = j;
			imgr1->j2[n] = j+1;
		}
		lat00 = y0 + (imgr1->j1[n] * dy);
		lat10 = y0 + (imgr1->j2[n] * dy);
		dyi = lat00 - lat10;
		pyi = lat00 - y;
		imgr1->py[n] = pyi / dyi;
		if(imgr1->j1[n] < 0 || fabs(dyi) < 0.000001 || imgr1->py[n] > 1.00 || imgr1->py[n] < 0.00) {
//                      printf("problem with j index \n");
			imgr1->py[n] = 9.99;
		}

		lon1 = x0 + ((i-1) * nwp->dlon);
		lon2 = x0 + ((i+1) * nwp->dlon);
		if(lon1 < 0.0 && x0 == 0.0) lon1 = lon1 + 360.0;
		if( fabs(x-lon1) < fabs(x-lon2) ) {
			imgr1->i1[n] = i-1;
			imgr1->i2[n] = i;
			i1_temp = i-1;
			i2_temp = i;
		}else {
			imgr1->i1[n] = i;
			imgr1->i2[n] = i+1;
			i1_temp = i;
			i2_temp = i+1;
		}

		if(i1_temp < 0) {
			i1_temp = i1_temp + 360;
			i2_temp = i2_temp + 360;
		}

		if(imgr1->i1[n] <= -1) imgr1->i1[n] = imgr1->i1[n] + 360;
		if(imgr1->i2[n] <= -1) imgr1->i2[n] = imgr1->i2[n] + 360;
		if(imgr1->i1[n] >= 360) imgr1->i1[n] = imgr1->i1[n] - 360;
		if(imgr1->i2[n] >= 360) imgr1->i2[n] = imgr1->i2[n] - 360;

		lon00 = x0 + (i1_temp * nwp->dlon);
		lon10 = x0 + (i2_temp * nwp->dlon);
		dxi = lon00 - lon10;
		pxi = lon00 - x;
		imgr1->px[n] = pxi/dxi;
                if(fabs(dxi) < 0.000001 || imgr1->px[n] > 1.00 || imgr1->px[n] < 0.00) {
//  	          printf("problem with i index \n");
//   	          printf("%ld %f %f %f %f %f %f %d %f %f %f %f %d %d %f %f %f \n", n, imgr1->lon[n], dx, imgr1->bounds[1],
//   			imgr1->bounds[3], x0, x, i, lon1, lon2, lon00, lon10, imgr1->i1[n], imgr1->i2[n], dxi, pxi, imgr1->px[n]);
			imgr1->px[n] = 9.99;
                }

//                  if(n<100) {
//                    printf("interp: %ld %f %f %f %f %d %f %f %d %d %f %f %f \n", n, nwp->first_lat, nwp->last_lat, y, y0, j,
//    		             lat1, lat2, imgr1->j1[n], imgr1->j2[n], dyi, pyi, imgr1->py[n]);
//                    printf("%ld %f %f %f %f %d %f %f %f %f %d %d %f %f %f \n", n, imgr1->lon[n], x0, nwp->last_lon, x, i, lon1,
//       		     lon2,lon00, lon10, imgr1->i1[n], imgr1->i2[n], dxi, pxi, imgr1->px[n]);
//                    printf("interpc: %ld %f %f %d %d %f %d %d %f %d %d \n", n, y, x, imgr1->j1[n], imgr1->j2[n], imgr1->py[n],
//    		              imgr1->i1[n], imgr1->i2[n], imgr1->px[n], j, i);
//                  }


		//  End of RAF additions

		imgr1->index_nwp[n] = min(max(nwp->nlon*j + i,0),nwp->ncells-1);
		imgr1->index_vza[n] = vza_bin(imgr1->cos_satzen[n]);
		imgr1->index_vza[n] = min(max(imgr1->index_vza[n],0),nwp->rtm_nvzen-1);

          }
          else {

	    /*Added by G.Britzolakis on 01/13/2014 to check for fill values in the index_nwp*/
	    imgr1->index_nwp[n] = -999;
            imgr1->px[n] = 9.99;
            imgr1->py[n] = 9.99;

          }

	}
//      printf("Leaving nwp_grid2pixel()...\n");

}

/*******************************************************************************
/
*******************************************************************************/
int vza_bin(float cos_vza) {

     return (int) (cos_vza / RTM_VZA_BINSIZE);
}



float vza_mid(float cos_vza) {

     int bin;
     float temp;

     bin = vza_bin(cos_vza);
     temp = bin * RTM_VZA_BINSIZE + RTM_VZA_BINSIZE / 2.;
     if (temp > 1.0) temp = 1.0;

     return acos(temp)/DTOR;
}
