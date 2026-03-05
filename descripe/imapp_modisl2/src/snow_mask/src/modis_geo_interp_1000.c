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
#include <stdio.h>
#include <stdlib.h>

float * modis_geo_interp_1000 (int ncols_5, int nrows_5, float *ptr5000)

{

  int
   nscan,
   scan,
   ncols_1000,
   i,j,ip,jp,
   j0,j1,
   k0,k1,
   m,n,
   DEBUG;
   
  long
   icol,irow,
   npts_1,
   index,
   index1,
   index2,
   index3,
   index4;
  
  float
   *ix,
   *jy,
   *temp,
   *ptr1000,
   dx,dx1,dy,dy1,
   slope,
   b;
  
  DEBUG = 9999;
   
  if (ncols_5 != 270 && ncols_5 != 271) {
    fprintf(stderr,"\nInput 5 km array must have 270 or 271 columns, not %d columns-exiting\n",ncols_5);
    exit(EXIT_FAILURE);
  }
  
  ncols_1000 = 1354;
  nscan = nrows_5/2;
  npts_1 = (long)ncols_1000*nscan*10;
  
  if (DEBUG < nscan) {
    fprintf(stderr,"In modis_geo_interp_1000\n");
    fprintf(stderr,"Interpolating %d earth scans.\n",nscan);
    fprintf(stderr,"Size of 5 km data: ncols=%d, nrows=%d.\n",ncols_5,nrows_5);
  }
  
  ptr1000 = (float *) malloc(npts_1*sizeof(float));
  if (ptr1000 == NULL) {
    fprintf(stderr,"Cannot allocate memory for ptr1000\n");
    perror("modis_geo_interp_1000");
    exit(EXIT_FAILURE);
  }
  
  ix = (float *) malloc(ncols_1000*sizeof(float));
  if (ix == NULL) {
    fprintf(stderr,"Cannot allocate memory for ix\n");
    perror("modis_geo_interp_1000");
    exit(EXIT_FAILURE);
  }
  
  jy = (float *) malloc(10*sizeof(float));
  if (jy == NULL) {
    fprintf(stderr,"Cannot allocate memory for jy\n");
    perror("modis_geo_interp_1000");
    exit(EXIT_FAILURE);
  }
  
  /*Set x and y interpolate locations.*/
  for (icol=0; icol < ncols_1000; icol++) ix[icol] = (float)icol*0.2 - 0.4;
  for (irow=0; irow < 10; irow++) jy[irow] = ((float)irow*0.2) - 0.4;
  
  /*Allocate memory for a temp array that will hold a single 5 km earth scan (10 lines/detectors)*/
  /*Memory for an extra line is allocated to avoid seg faults. This should have no bearing*/
  /*on the results since the last two scan lines are treated in a separate manner.*/
  temp = (float *) malloc((3*ncols_5)*sizeof(float));
  if (temp == NULL) {
    fprintf(stderr,"Cannot allocate memory for temp\n");
    perror("modis_geo_interp_1000");
    exit(EXIT_FAILURE);
  }
  
  for (scan=0; scan < nscan; scan++) {
  
    if (DEBUG < nscan) fprintf(stderr,"scan=[%d]\n",scan);
    
    j0 = 2*scan;
    j1 = 2*scan + 1;
    k0 = 10*scan;
    k1 = 10*scan + 9;
    
    /*Get region of small array that is to be expanded.*/
    index = 0;
    for (m=0; m<2; m++) {
      for (n=0; n<ncols_5; n++) {
	temp[index] = ptr5000[((j0+m)*ncols_5)+n];
	index++;
      }
    }
    
    /*Interpolate out to larger array.*/
    for (irow=k0; irow<=k1; irow++) {
      j = (int)jy[irow-k0];
      jp = j + 1;
      dy = jy[irow-k0] - (float)j;
      dy1 = 1.0 - (float)dy;
      for (icol=0; icol<ncols_1000; icol++) {
        if (scan == DEBUG) fprintf(stderr,"irow=[%d],icol=[%d]\n",irow,icol);
	i = (int)ix[icol];
	ip = i + 1;
	dx = ix[icol] - (float)i;
	dx1 = 1.0 - (float)dx;
	index1 = i + (j*(ncols_5));
        index2 = i + (jp*(ncols_5));
        index3 = ip + (j*(ncols_5));
        index4 = ip + (jp*(ncols_5));
	
	if (scan == DEBUG) {
	  fprintf(stderr,"i=%d  ip=%d  ix=%f  j=%d  jp=%d  jy=%f\n",i,ip, ix[icol],j,jp,jy[irow]);
	  fprintf(stderr,"dx=%f  dx1=%f  dy=%f  dy1=%f\n",dx,dx1,dy,dy1);
	  fprintf(stderr,"index1=%d  index2=%d  index3=%d  index4=%d\n",index1,index2,index3,index4);
	  fprintf(stderr,"temp1=%f  temp2=%f  temp3=%f  temp4=%f\n",temp[index1],temp[index2],temp[index3],temp[index4]);
	}
	
	/*Perform bilinear interplation.*/
	
	/*Standard bilinear interpolation*/
	if (ip < ncols_5) {
	  ptr1000[((irow)*ncols_1000)+icol] = temp[index1]*dx1*dy1 +
                                              temp[index2]*dx1*dy  +
		       	                      temp[index3]*dx*dy1  +
			                      temp[index4]*dx*dy;
        }
	/*Some special nearest neighbor treatment must be used when ip is out of bounds*/
	else {	  
	  ptr1000[((irow)*ncols_1000)+icol] = temp[(ncols_5-1) + (j*(ncols_5))]*dx1*dy1 +
                                              temp[(ncols_5-1) + (jp*(ncols_5))]*dx1*dy  +
		       	                      temp[(ncols_5-1)+(j*(ncols_5))]*dx*dy1  +
			                      temp[(ncols_5-1)+(jp*(ncols_5))]*dx*dy;
	}
      }
    }
    
    /*Treat scan line boundaries differently.*/
    for (icol=0; icol<(ncols_1000); icol++) {
      /*Replace first two and last two output pixels along track with 
        linearly extrapolated values*/
      slope = (ptr1000[((k0+7)*ncols_1000)+icol] - ptr1000[((k0+2)*ncols_1000)+icol])/
	      (jy[7] - jy[2]);
      b = ptr1000[((k0+7)*ncols_1000)+icol] - slope*jy[7];
      ptr1000[((k0+0)*ncols_1000)+icol] = slope*jy[0] + b;
      ptr1000[((k0+1)*ncols_1000)+icol] = slope*jy[1] + b;
      ptr1000[((k0+8)*ncols_1000)+icol] = slope*jy[8] + b;
      ptr1000[((k0+9)*ncols_1000)+icol] = slope*jy[9] + b;
    }
    
    for (irow=k0; irow<=k1; irow++) {
      /*Replace first two outputpixels across track with linearly extrapolated values.*/
      slope = (ptr1000[(irow*ncols_1000)+7] - ptr1000[(irow*ncols_1000)+2]) / (ix[7] - ix[2]);
      b = ptr1000[(irow*ncols_1000)+7] - slope*ix[7];
      ptr1000[(irow*ncols_1000)+0] = slope*ix[0] + b;
      ptr1000[(irow*ncols_1000)+1] = slope*ix[1] + b;
      
      /*Replace last six outputpixels across track with linearly extrapolated values.*/
      slope = (ptr1000[(irow*ncols_1000)+1347] - ptr1000[(irow*ncols_1000)+1342]) / (ix[1347] - ix[1342]);
      b = ptr1000[(irow*ncols_1000)+1347] - slope*ix[1347];
      ptr1000[(irow*ncols_1000)+1348] = slope*ix[1348] + b;
      ptr1000[(irow*ncols_1000)+1349] = slope*ix[1349] + b;
      ptr1000[(irow*ncols_1000)+1350] = slope*ix[1350] + b;
      ptr1000[(irow*ncols_1000)+1351] = slope*ix[1351] + b;
      ptr1000[(irow*ncols_1000)+1352] = slope*ix[1352] + b;
      ptr1000[(irow*ncols_1000)+1353] = slope*ix[1353] + b;
    }
  
  }
  
  free(ix);
  free(jy);
  free(temp);
  
  return(ptr1000);
  
}
  
