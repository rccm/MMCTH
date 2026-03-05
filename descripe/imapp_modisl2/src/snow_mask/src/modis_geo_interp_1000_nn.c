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
/* This version of the modis geo interpolate code gives the "nearest neighbor"
for each data point. Orignal code devloped by Michael Pavolonis. */

#include <stdio.h>
#include <stdlib.h>

float * modis_geo_interp_1000_nn (int ncols_5, int nrows_5, float *ptr5000)

{

  int irow, icol, i5, j5, iclose5, jclose5, irem, jrem;
  int ncols_1000, nrows_1000;
  int rems[5] = {0, 0, 0, 1, 1};
  long npts_1000, index1, index5;
  float *ptr1000;
  
  if (ncols_5 != 270 && ncols_5 != 271) {
    fprintf(stderr,"\nInput 5 km array must have 270 or 271 columns, not %d columns-exiting\n",ncols_5);
    exit(EXIT_FAILURE);
  }

  ncols_1000 = 1354;
  nrows_1000 = nrows_5 * 5;
  npts_1000 = (long) ncols_1000*nrows_1000;

  ptr1000 = (float *) malloc(npts_1000*sizeof(float));
  if (ptr1000 == NULL) {
    fprintf(stderr,"Cannot allocate memory for ptr1000\n");
    perror("modis_geo_interp_1000_az");
    exit(EXIT_FAILURE);
  }
  
  /* Loop through each row in the 1 km, assigning values from the "closest"
     value in the 5 km array. */

  for (irow=0; irow<nrows_1000; irow++) {
	  i5 = irow / 5;
	  irem = irow % 5;
	  iclose5 = i5 + rems[irem];
	  if (iclose5 > nrows_5-1) iclose5 = nrows_5 - 1;
	  for (icol=0; icol<ncols_1000; icol++) {
		  j5 = icol / 5;
		  jrem = icol % 5;
		  jclose5 = j5 + rems[jrem];
		  if (jclose5 > ncols_5-1) jclose5 = ncols_5 - 1;
		  index1 = irow * ncols_1000 + icol;
		  index5 = iclose5 * ncols_5 + jclose5;
		  ptr1000[index1] = ptr5000[index5];
//		  if (ptr5000[index5] < 29) printf("zero prt5000 at i,j: %d, %d\n", iclose5, jclose5);
	  }
  }

  return(ptr1000);
  
}
  
