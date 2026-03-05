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
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  This program is for calcuating the clear sky inversion depth and strength
  in either Antarctica or the artic from MODIS DB and DAAC. Currently the 
  high altitude equation is used for the Antarctic and low altitude equation.
  The equations that are used are from Liu and Key (2003). 

  EXTRA LIBRARIES REQUIRED: hdf.h, subpixel_cloudfrac.h

  REQUIRED INPUTS (taken in as command line arguments):

  hemisphere - 0= Antarctic, 1 = Arctic
  filee1 = Level 1b file name and location
  filee2 = MOD35 file and location

  OUTPUT:

  inverdepthbin = clear sky inversion depth in meters * 10, 
                  output to binary file 'invdepth.dat'

  inverdepthbin = clear sky inversion strength in degrees C * 10, 
                  output to binary file 'invstrength.dat'

 Program created by William Straka III 
 University of Wisconsin - Madison
 Space Science and Engineering Center (SSEC/CIMSS)
 VERSION 1.0 - 9/5/05

 Originally written in IDL by Yinghui Liu, UW Madison
 University of Wisconsin - Madison
 Space Science and Engineering Center (SSEC/CIMSS)

 Under the Guidance and funding of Dr. Jeffrey Key, NOAA/NESDIS/ASPB

EXAMPLE Execution : 
   inversions.exe [Hemisphere] [MODIS L1B 1km file] [MODIS cloud mask file]
   
Example for Antarctica:

	inversions.exe 0 a1.05200.0711.1000m.hdf a1.05200.0711.mod35.hdf

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#include <hdf.h>

/* Function prototypes */

void modis_level1b_read (char *filename, int option, float band, void **buf,
                         int *nlines, int *npixels, float *earth2sun,
			 int SUB_FLG);
void modis_level35_read (char *filename, int option, void **buf, int *nlines, 
                         int *npixels);
float * modis_geo_interp_1000 (int ncols_5, int nrows_5, float *ptr5000);

/* Main program */

int main(short argc, void *argv[])

{
	
system("date");      //Print the system time for the log

char *filee1, *filee2;
FILE *fp;
int hemis;
int i;

/* NEED ALLL VARIABLES HERE!!!!*/

int nlines1, npixels1, looper;
float earth2sun;
float *mod_tb11, *mod_tb12, *mod_tb7,*inverdepth, *inverstrength;
uchar8 *mod_cldmsk;
short *inverdepthbin, *inverstrengthbin;
float mind, maxd, mins, maxs, aved, aves, stot, dtot;
float bt7_11, bt11_12;
int igood;

/*get argument values*/

hemis = atoi(argv[1]);
filee1 = argv[2];
filee2 = argv[3];

/* Read in MODIS1B data. Developed by Michael Pavolonis*/
/* need 11, 12 and 7.2 micron brightness temps for depth and strenghts*/
/* Channel 31, 32, 28*/

modis_level1b_read(filee1, 5, 31, (void **) &mod_tb11,  &nlines1, &npixels1, &earth2sun, 0);
modis_level1b_read(filee1, 5, 32, (void **) &mod_tb12,  &nlines1, &npixels1, &earth2sun, 0);
modis_level1b_read(filee1, 5, 28, (void **) &mod_tb7,  &nlines1, &npixels1, &earth2sun, 0);

/* Read in MOD35 data. Developed by Michael Pavolonis */

modis_level35_read(filee2, 1, (void **) &mod_cldmsk, &nlines1,&npixels1);

/*loop over all pixels*/

looper = npixels1*nlines1;
igood = 0;
stot = 0.;
dtot = 0.;
mins = 100.;
maxs = 0.;
mind = 10000.;
maxd = 0.;

inverdepth = (float *) malloc(looper*sizeof(float));
inverdepthbin = (short *) malloc(looper*sizeof(short));
inverstrength = (float *) malloc(looper*sizeof(float));
inverstrengthbin = (short *) malloc(looper*sizeof(short));

for (i=0; i<looper; i++) {
inverstrength[i] = -1.0;
inverdepth[i] = -1.0;

if(mod_cldmsk[i] > 1) {      /* Cloud mask is "cloudy" when <= 1) */
  bt7_11 = mod_tb7[i]-mod_tb11[i];
  bt11_12 = mod_tb11[i] - mod_tb12[i];
  if (bt7_11 >= -20. && bt7_11 <= 100.) {
  /* Use high altitude equations for Antarctica */
  if (hemis == 0) {
	inverstrength[i] = 23.6 + (1.28*bt7_11) - (2.61*bt11_12) -
			(0.059*mod_tb11[i]) + (0.035*pow(bt7_11,2));

	inverdepth[i] = 1806.5 + (33.9*bt7_11) + (103.7*bt11_12) -
			(5.8*mod_tb11[i]) + (0.2*pow(bt7_11,2));
  }

  else {
	inverstrength[i] = 32.2 + (0.84*bt7_11) - (4.63*bt11_12) -
			(0.081*mod_tb11[i]) + (0.021*pow(bt7_11,2));

	inverdepth[i] = 720.3 + (44.1*bt7_11) - (133.5*bt11_12) -
			(0.45*mod_tb11[i]) + (1.27*pow(bt7_11,2));
  }
  if (inverstrength[i] > 0. && inverstrength[i] < 100. && inverdepth[i] > 0. && inverdepth[i] < 2000.) {
    if (inverstrength[i] < mins) mins = inverstrength[i];
    if (inverstrength[i] > maxs) maxs = inverstrength[i];
    if (inverdepth[i] < mind) mind = inverdepth[i];
    if (inverdepth[i] > maxd) maxd = inverdepth[i];
    stot += inverstrength[i];
    dtot += inverdepth[i];
    igood += 1;
  }
  else {
    inverstrength[i] = -1.0;
    inverdepth[i] = -1.0;
  }
  }
}

inverdepthbin[i] = (short)(inverdepth[i]*10);
inverstrengthbin[i] = (short)(inverstrength[i]*10);
}

printf("Number of good (retrievable) and total pixels: %d, %d\n", igood, looper);
if (igood > 0) {
  printf("Average, min, max strength of good pixels: %6.3f, %6.3f, %6.3f\n", stot/(float)igood, mins, maxs);
  printf("Average, min, max depth of good pixels: %6.3f, %6.3f, %6.3f\n", dtot/(float)igood, mind, maxd);
}

/* output data to file */

fp = fopen("invdepth.dat", "w"); 
fwrite(inverdepthbin,sizeof(short),looper,fp);
fclose(fp);

fp = fopen("invstrength.dat", "w"); 
fwrite(inverstrengthbin,sizeof(short),looper,fp);
fclose(fp);

free(inverdepth);
free(inverdepthbin);
free(inverstrength);
free(inverstrengthbin);

}
