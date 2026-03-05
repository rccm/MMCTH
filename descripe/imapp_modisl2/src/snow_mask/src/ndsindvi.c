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

  This program is for comparing the snow ice flags
  This will produce a binary snow mask as well as the NDVI dataset
  	NDVI results in -1 to 1, so we need to figure this out

  These values were chosen so that the binary reader in McIDAS can read them in
  
  Comparison is going to be done with MOD35 processing flag

  EXTRA LIBRARIES REQUIRED: hdf.h, anisdata.h, subpixel_cloudfrac.h

  REQUIRED INPUTS (taken in as command line arguments):

  filee1 = Level 1b file name and location
  filee2 = MOD35 file and location

  OUTPUT:
  snow ice binary map (short) in binary file output: ndsialb.dat, ndvistr.dat,


     In order to get it to display in McIDAS, we will do a similar measure as albedo

  	Output in same directory as executable. Displayed in McIDAS


Procedures in this file used to calculate ndsindvi. They are:
binsnow (binary snow/ice mask), based on NDSI ATBD
ndvi (NDVI)


  Others are included in compile statement

 Compilation:  ./compile

 Program adapted by William Straka III (SSEC/CIMSS) and Jeff Key (NOAA/ASPT)
 Original code from NGST


EXAMPLE ANTARCTIC RUN: 
	./snowice ./a1.05200.0711.1000m.hdf ./a1.05200.0711.mod35.hdf AQUA

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

void modis_level1b_read (char *filename, int option, float band, void **buf, int *nlines, int *npixels, float *earth2sun, int SUB_FLG);
void modis_level35_read (char *filename, int option, void **buf, int *nlines, int *npixels);
float ndsi (float modr1, float modr2, float modr4, float modr6, float bt31, float solzen, float landmsk, float lat);
float * modis_geo_interp_1000 (int ncols_5, int nrows_5, float *ptr5000);
float * modis_geo_interp_1000_nn (int ncols_5, int nrows_5, float *ptr5000);
float zen2scanner(float sza);

/* Main program */

int main(int argc, void *argv[])
{

/*input variables */
char *filee1, *filee2, *sater, *snowmap;

FILE *fp;
int i, satdiff;
int igood, ibad;
char sattest[5]="TERRA";
char satter[5];
char sattest1[5];

/* NEED ALLL VARIABLES HERE!!!!*/

int nlines1, nlines2, npixels1, errout, npixels2, looper;
float earth2sun, snang;
float *mod_r1, *mod_r2, *mod_r6, *mod_r4, *mod_bt31, *tsurf2, *tsurf1;
float *mod_lat5, *mod_lat;
uchar8 *mod_cldmsk, *mod_simsk, *mod_landmsk;


float *ndvi, *sibin;
short *snmapbin, *ndvibin;

float *mod_solzen5, *mod_satzen5, *mod_relaz5;
float *mod_solzen1,*mod_satzen1, *mod_relaz1;




/*set errout to 1 if you want to check for errouneous values*/

errout = 1;

/*get argument values*/

if ( argc !=5){
printf("You must input the correct number of command line arguments\n");
printf("Your command line should be as follows:\n");
printf("./snowndsi.exe SAT FIL1KM MASK OUTFILE\n");
printf("SAT is the satellite platform (Aqua or Terra)\n");
printf("FIL1KM is the MODIS L1B 1000 meter resolution HDF file\n");
printf("MASK is the MODIS L1B 1000 meter Cloud Mask  HDF file\n");
printf("OUTFILE is the Output 1000 meter NDSI Snow Map Binary file\n");
exit(0);
}

sater = argv[1];
filee1 = argv[2];
filee2 = argv[3];
snowmap = argv[4];

printf("Satellite name: %s \n", sater);
printf("Level 1 file name: %s \n", filee1);
printf("MOD35 file name: %s \n", filee2);
printf("Output Snow Map file name: %s \n", snowmap);


strcpy(satter, sater);
/* This will see if we are dealing with TERRA or AQUA, since they handle the NDSI differently */

satdiff = strcmp(satter, "Terra");

/* Prints the value of the comparison between input and TERRA. If TERRA it
should be 0. */

/* Read in MODIS1B data. Developed by Michael Pavolonis*/
/* Read in Channel 1, 2, 4, and 6 (or 7) Reflectances*/
modis_level1b_read(filee1, 4, 1, (void **) &mod_r1,  &nlines1, &npixels1, &earth2sun, 0);
modis_level1b_read(filee1, 4, 2, (void **) &mod_r2,  &nlines1, &npixels1, &earth2sun, 0);
modis_level1b_read(filee1, 4, 4, (void **) &mod_r4,  &nlines1, &npixels1, &earth2sun, 0);

/*Read in CH31 BT */
modis_level1b_read(filee1, 5, 31, (void **) &mod_bt31,  &nlines1, &npixels1, &earth2sun, 0);

/* Read in Latitude */
modis_level1b_read(filee1, 6, 31, (void **) &mod_lat5,  &nlines2, &npixels2, &earth2sun, 0);
mod_lat= modis_geo_interp_1000(npixels2,nlines2,mod_lat5);


if (satdiff < 0) {
printf("AQUA, reading in CH 7\n");
modis_level1b_read(filee1, 4, 7, (void **) &mod_r6,  &nlines1, &npixels1, &earth2sun, 0);
}

if (satdiff >= 0) {
printf("TERRA, reading in CH 6\n");
modis_level1b_read(filee1, 4, 6, (void **) &mod_r6,  &nlines1, &npixels1, &earth2sun, 0);
}

printf("Level 1b data read\n");

/* Read in MOD35 cloud mask data. Developed by Michael Pavolonis */

modis_level35_read(filee2, 1, (void **) &mod_cldmsk, &npixels1,&nlines1);


/* Read in MOD35 land mask data. Developed by Michael Pavolonis */
/*  0=Water, 1=Coastal, 2=Desert, 3=Land */

modis_level35_read(filee2, 2, (void **) &mod_landmsk, &npixels1,&nlines1);

printf("MOD35 data read\n");


/* solzen, satzen, relaz */

modis_level1b_read(filee1, 8, 31, (void **) &mod_solzen5,  &nlines2, &npixels2, &earth2sun, 02);
modis_level1b_read(filee1, 9, 31, (void **) &mod_satzen5,  &nlines2, &npixels2, &earth2sun, 2);
modis_level1b_read(filee1, 10, 31, (void **) &mod_relaz5,  &nlines2, &npixels2, &earth2sun, 2);

/* interpolation of solzen, satzen and relaz. */

mod_satzen1= modis_geo_interp_1000(npixels2,nlines2,mod_satzen5);
mod_solzen1= modis_geo_interp_1000(npixels2,nlines2,mod_solzen5);
mod_relaz1= modis_geo_interp_1000_nn(npixels2,nlines2,mod_relaz5);

printf("Angles computed and interpolated to 1km\n");

/*loop over all pixels*/

looper = npixels1*nlines1;

/*Snow/ice mask, NDVI calculation */
sibin = (float *) malloc(looper*sizeof(float));
ndvi = (float *) malloc(looper*sizeof(float));

/*Binary for output */
snmapbin = (short *) malloc(looper*sizeof(short));
ndvibin = (short *) malloc(looper*sizeof(short));

igood = 0;
ibad = 0;


for (i=1; i<looper; i++) {


	sibin[i]= -1. ;
	
/*TEST IN AREAS THAT ARE CLEAR AND VALID */
	if(mod_cldmsk[i] > 2 && mod_r1[i] != -999. && mod_bt31[i] != -999. && mod_r2[i] !=-999. && mod_r4[i] !=-999.) {

/*      printf("Band 6 values: %f \n", mod_r6[i]); */
//        printf("Band 6 values: %f \n", mod_r6[i]); 
	if(mod_r6[i] != -999.) {
	if(mod_solzen1[i] <= 85.) {
	if(mod_solzen1[i] >= 0.) {
	if(mod_satzen1[i] >= 0.) {

/* Calculate NDSI AND NDVI */
	sibin[i]=ndsi(mod_r1[i], mod_r2[i], mod_r4[i], mod_r6[i], mod_bt31[i], mod_solzen1[i],mod_landmsk[i], mod_lat[i]);
	
	ndvi[i] = (mod_r2[i] - mod_r1[i]) / (mod_r2[i] + mod_r1[i]);
	
	/*To get the NDVI in the range we want, we add 2.*/
	
	ndvi[i] = ndvi[i] +2;
	
	
	
	}}}}} /* End of test for good data */




	snmapbin[i]=(short)(sibin[i]*1000);
	/*NDVI is multiplied by 10 to get decimal place, then 10 again for 
		scaling in McIDAS */
	ndvibin[i]=(short)(ndvi[i] *100);

	if (snmapbin[i] > -1) {
			igood += 1;
	}
	else {
		ibad +=1;
	}
	
}


/* output data to file. Specified folder must exist*/


printf("Number of good (retrievable), bad (out of range), and total pixels: %d, %d, %d\n", igood, ibad, looper);


printf("calculation complete\n");


fp = fopen(snowmap, "w"); 
fwrite(snmapbin,sizeof(short),looper,fp);
fclose(fp);

printf("File stored\n");


free(sibin); 
free(snmapbin);

free(ndvi);
free(ndvibin);


exit (0);
}


/*********************** Snow ice binary map ********************/

float ndsi (float modr1, float modr2, float modr4, float modr6, float bt31, float solzen, float landmsk, float lat)
{



/*  Only pixels with valid */
  float ndsi_thre1 = .4;
  float ndsi_thre3 = 1.0;
  float ndsi_thre2 = .1;
  float sm_m02 = 0.25;
  float sm_ref2 = 0.11;
  float sm_ref2_low = 0.9;

// test ref2_low with .6

  float sm_mnir = 0.15;
  float sm_mnir_low = 0.9;
  float sm_ndsi [] = {0.40, 0.15, 0.70, 0.60};
  float sm_ndvi [] = {0.10, 0.05, -0.02};
  float sm_bt11 [] = {280.0, 275.0, 253.0, 243.0, 233.0, 220.0};



  float snmap, ndsi, ref, ndvi;
  
  /* fill value for binary map is no snow = 0.01, snow = 1. -1 indicates a non-NDSI valid pixel*/

snmap = 0.01;

/* normalize reflectance by solar zenith */
modr1= modr1 / cos(solzen*DEGRAD);
modr2= modr2 / cos(solzen*DEGRAD);
modr4= modr4 / cos(solzen*DEGRAD);
/* Channel 6 or 7 reflectance, depending on satellite*/
modr6= modr6 / cos(solzen*DEGRAD); 

/* bt31 = channel 31 brightness temp*/


ndsi = (modr4 - modr6) / (modr4 + modr6);
ndvi = (modr2 - modr1) / (modr2 + modr1);

//printf("sm_ref2, sm_bt11[1] , sm_ndsi[0]: %f, %f, %f\n", sm_ref2, sm_bt11[1], sm_ndsi[0]);
		
if (modr4 >= .1) {

	if (ndsi >= ndsi_thre1) {
		if (modr2 > .11){
			snmap = 1.;
		}
	}
	
}


/*simplified Klein et al. NDVI test, as performed by MODIS MOD35.

Determines if snow or ice is present on the Earth's surface
       Uses NDSI algorithm by Dorothy Hall, George Riggs and 
       Vince Salmonson.
       Other spectral algorithms provided by R. Frey 
       
       Based off of MERSI cloud mask NDSI calculations.
       
       Basically, these are reasonableness tests
       
       */
	
	/* LAND PIXEL TESTS */

	if ((landmsk != 0.) && (landmsk != 1.)) {
		if ( (ndsi > sm_ndsi[1]) &&  (modr4 < sm_m02) && (modr6 < sm_mnir) ) {
		snmap = 1.;
		}
		else if  ((ndsi > sm_ndsi[1]) && (ndvi > sm_ndvi[0])) {
		snmap = 1.;	
		}
		
		/* Reasonableness tests for land. Note, since we don't have elev. map,
			we use warmer BT31 threshold*/
		
		if ((lat > -60) && (lat <=20 )) {
			if (bt31 < sm_bt11[2]) {
				snmap = 0.01;
			}
		}
		
		if ((lat > 20) && (lat <=40 )) {
			if ((bt31 < sm_bt11[4])) {
				snmap = 0.01;
			}

/*ref1 test */
		if (modr1 < 0.2) {
			snmap = 0.01;
			}			

		}
		
		if ((lat > 40) && (lat <=60 )) {
			if (bt31 < sm_bt11[5]) {
				snmap = 0.01;
			}
		}

	}

	/* TEST OCEAN PIXELS */
	if (landmsk == 0.) {
		if ((modr2 > sm_ref2) && (bt31 <= sm_bt11[1]) && (ndsi > sm_ndsi[0])) {
			snmap = 1.;
		}
/* Deep Ocean test */
		if ((lat >= -55.0) && (lat <= 40.0)) {
			if (bt31 < sm_bt11[3]) {
				snmap = 0.01;
				}
/*ref1 test */
		if (modr1 < 0.2) {
			snmap = 0.01;
			}			
/* Mid latitude BT31 test */
		if (bt31 > 280.0 ) {
			snmap = 0.01;
			}

/*Open ocean test using NDVI. */
		if (ndvi < 0.05) {
			snmap = 0.01;
			}			

		}
	}

	/* TEST COAST PIXELS */
	if (landmsk == 1.){
	
		if ( (ndsi > sm_ndsi[1]) &&  (modr4 < sm_m02) && (modr6 < sm_mnir) ) {
		snmap = 1.;
		}
		else if ((modr2 > sm_ref2) && (bt31 <= sm_bt11[1]) && (ndsi > sm_ndsi[0])){
			snmap = 1.;	
		}
		
		/*ref1 test */
		if (modr1 < 0.2) {
			snmap = 0.01;
			}			

		
	}

return snmap;

}

