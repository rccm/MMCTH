#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"
#include "sounderL1_leocat.h"

void read_hdf_all(int32, char8 *, int32 dimen[], void **);
void read_ascii_ch_prop(char *, int, sounderL1 *);

void airs_level1b_read (char *filename, sounderL1 *sndr1)
{

  char *rout = {"airs_level1b_read"};
  int32 id;
  int32 status;
  int32 dimen[MAX_VAR_DIMS];
  int nn;
  int n;
  int index;
  float32 *buffer_f32;
  float64 *buffer_f64;
  
  id = SDstart(filename, DFACC_READ);
  if (id == FAIL) {
    fprintf(stderr,"%sInvalid HDF file, %s\n",EXE_PROMPT,filename);
    exit(EXIT_FAILURE);
  }
  
  fprintf(stdout,"%sReading in AIRS Latitude.\n",EXE_PROMPT);
  read_hdf_all(id, "Latitude", dimen, (void **) &buffer_f64);
  if ((sndr1->lat = (float32 *) malloc(dimen[0]*dimen[1]*sizeof(float32))) == NULL)
      error_allo(rout,"sndr1->lat");
  for (n=0; n<dimen[0]*dimen[1]; n++) sndr1->lat[n] = (float32) buffer_f64[n];
  free(buffer_f64);
    
  fprintf(stdout,"%sReading in AIRS Longitude.\n",EXE_PROMPT);
  read_hdf_all(id, "Longitude", dimen, (void **) &buffer_f64);
  if ((sndr1->lon = (float32 *) malloc(dimen[0]*dimen[1]*sizeof(float32))) == NULL)
      error_allo(rout,"sndr1->lon");
  for (n=0; n<dimen[0]*dimen[1]; n++) sndr1->lon[n] = (float32) buffer_f64[n];
  free(buffer_f64);
  
  fprintf(stdout,"%sReading in AIRS Satellite Zenith Angle.\n",EXE_PROMPT);
  read_hdf_all(id, "satzen", dimen, (void **) &sndr1->satzen);
  
  fprintf(stdout,"%sReading in AIRS Solar Zenith Angle.\n",EXE_PROMPT);
  read_hdf_all(id, "solzen", dimen, (void **) &sndr1->solzen);
  
  fprintf(stdout,"%sReading in AIRS Relative Azimuth.\n",EXE_PROMPT);
  read_hdf_all(id, "satazi", dimen, (void **) &buffer_f32);
  read_hdf_all(id, "solazi", dimen, (void **) &sndr1->relaz);
  for (n=0; n<dimen[0]*dimen[1]; n++) sndr1->relaz[n] = 180.0 - (sndr1->relaz[n] - buffer_f32[n]);
  free(buffer_f32);
  
  fprintf(stdout,"%sReading in AIRS Surface Elevation.\n",EXE_PROMPT);
  read_hdf_all(id, "topog", dimen, (void **) &sndr1->zsfc);
  
  fprintf(stdout,"%sReading in AIRS Land Fraction.\n",EXE_PROMPT);
  read_hdf_all(id, "landFrac", dimen, (void **) &sndr1->landfrac);
  
  fprintf(stdout,"%sReading in AIRS State.\n",EXE_PROMPT);
  read_hdf_all(id, "state", dimen, (void **) &sndr1->state);
  
  fprintf(stdout,"%sReading in AIRS CalFlag.\n",EXE_PROMPT);
  read_hdf_all(id, "CalFlag", dimen, (void **) &sndr1->calflg);
      
  fprintf(stdout,"%sReading in AIRS Radiances.\n",EXE_PROMPT);
  read_hdf_all(id, "radiances", dimen, (void **) &buffer_f32);
  if ((sndr1->chnum = (int *) malloc(dimen[2]*sizeof(int))) == NULL)
      error_allo(rout,"sndr1->chnum");
  if ((sndr1->wvnum = (float *) malloc(dimen[2]*sizeof(float))) == NULL)
      error_allo(rout,"sndr1->wvnum");
  if ((sndr1->nedt = (float *) malloc(dimen[2]*sizeof(float))) == NULL)
      error_allo(rout,"sndr1->nedt");
  if ((sndr1->fwhm = (float *) malloc(dimen[2]*sizeof(float))) == NULL)
      error_allo(rout,"sndr1->fwhm");
  if ((sndr1->ab_state = (int *) malloc(dimen[2]*sizeof(int))) == NULL)
      error_allo(rout,"sndr1->ab_state");
  if ((sndr1->rad_qc = (int *) malloc(dimen[2]*sizeof(int))) == NULL)
      error_allo(rout,"sndr1->rad_qc");
  if ((sndr1->use_flg = (int *) malloc(dimen[2]*sizeof(int))) == NULL)
      error_allo(rout,"sndr1->use_flg");
  sndr1->rad = calloc_2d_float_ptr("sndr1->rad", dimen[0]*dimen[1], dimen[2]);
  sndr1->bt = calloc_2d_float_ptr("sndr1->bt", dimen[0]*dimen[1], dimen[2]);
  read_ascii_ch_prop("./data/airs/detectors/L2.chan_prop.2003.11.19.v6.8.1.txt", dimen[2], sndr1);
  index = 0;
  for (nn=0; nn<dimen[0]*dimen[1]; nn++) {
    for (n=0; n<dimen[2]; n++) {
      sndr1->rad[nn][n] = buffer_f32[index];
      index++;
    }
  }
  free(buffer_f32);
  
  status = SDend(id);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot close %s\n",EXE_PROMPT,rout,filename);
    exit(EXIT_FAILURE);
  }
  
}

void read_hdf_all(int32 id, char8 *sds_name, int32 dimen[MAX_VAR_DIMS], void **data)
{

  char *rout = {"read_hdf_all"};
  char8 sds_name_check[MAX_NC_NAME];
  int n;
  int32 sd_id;
  int32 sds_index;
  int32 status;
  int32 data_type;
  int32 rank;
  int32 nattr;
  int32 *start;
  int32 *edge;
  int32 npts;
  void *buffer;
  
  sds_index = SDnametoindex(id,sds_name);
  
  sd_id = SDselect(id,sds_index);
  if (sd_id == FAIL) {
    fprintf(stderr,"%s%s-cannot select specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }
  
  status = SDgetinfo(sd_id,&sds_name_check[0],&rank,&dimen[0],&data_type,&nattr);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot get info about specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }
  
  if ((start = (int32 *) calloc(rank,sizeof(int32))) == NULL)
      error_allo(rout,"start");
      
  if ((edge = (int32 *) calloc(rank,sizeof(int32))) == NULL)
      error_allo(rout,"edge");  
  
  npts = 1;
  for (n=0; n<rank; n++) {
    start[n] = 0;
    edge[n] = dimen[n];
    npts *= dimen[n];
  }
  
  if ((buffer = (void *) malloc(npts*DFKNTsize(data_type))) == NULL)
      error_allo(rout,"buffer");
      
  *data = buffer;
      
  status = SDreaddata(sd_id,start,NULL,edge,buffer);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot read specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }
  
  status = SDendaccess(sd_id);
  if (status == FAIL) {
    fprintf(stderr,"%s%s-cannot end access to specified SDS, %s\n",EXE_PROMPT,rout,sds_name);
    exit(EXIT_FAILURE);
  }
  
  free(start);
  free(edge);
}

void read_ascii_ch_prop(char *filename, int nchan, sounderL1 *sndr1)
{
  const int nheader = 123;
  FILE *fptr;
  char header[500];
  int n;
  char array_name[5];
  int cal_index;
  float Cij;
  float centx;
  float centy;
  float rta;
  char comments[500];
  
  if ((fptr = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"%sCannot open AIRS channel property file %s - aborting\n",EXE_PROMPT,filename);
    exit(EXIT_FAILURE);
  }
  
  for (n=0; n<nheader; n++) {
    if (fgets(header, 500, fptr) == NULL) {
      fprintf(stderr,"%sError reading AIRS channel property file %s - aborting\n",EXE_PROMPT,filename);
      exit(EXIT_FAILURE);
    }
  }
  
  for (n=0; n<nchan; n++) {
    if (fscanf(fptr,"%d %f %s %d %f %f %f %f %f %f %d %d %d",
        &sndr1->chnum[n],&sndr1->wvnum[n],&array_name[0],&cal_index,&sndr1->nedt[n],&sndr1->fwhm[n],
	&Cij,&centx,&centy,&rta,&sndr1->ab_state[n],&sndr1->rad_qc[n],&sndr1->use_flg[n]) == EOF) {
      fprintf(stderr,"%sError reading AIRS channel property file %s - aborting\n",EXE_PROMPT,filename);
      exit(EXIT_FAILURE);
    }
    if (fgets(comments, 500, fptr) == NULL) {
      fprintf(stderr,"%sError reading AIRS channel property file %s - aborting\n",EXE_PROMPT,filename);
      exit(EXIT_FAILURE);
    }
  }
  
  fclose(fptr);
}
