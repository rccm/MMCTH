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

  This program is for calcuating the clear sky surface tempuratur in either
  Antarctica or the artic.

  EXTRA LIBRARIES REQUIRED: hdf.h, anisdata.h, subpixel_cloudfrac.h

  REQUIRED INPUTS (taken in as command line arguments):

  hemisphere - 0= Antarctic, 1 = Arctic
  filee1 = Level 1b file name and location
  filee2 = MOD35 file and location

  OUTPUT:
  Surface tempurature (short) in binary file output: surftemp.dat
  	Output in same directory as executable. Displayed in McCIDAS

To compile: 

gcc -o bin/surftemp surftemp.c modis_level1b_read.c modis_level35_read.c modis_geo_interp_1000.c zen2scanner.c -I$INCLUDES -L$LIBS -lmfhdf -ldf -ljpeg -lz -lm

Procedures in this file used to calculate surftemp. They are:
surftemp_clear

 Program adapted by William Straka III (SSEC/CIMSS) and Jeff Key (NOAA/ASPT)
 Original code from Jeff Key (NOAA/ASPT)

 Version 1.2a 10/21/05 BETA VERSION!!!!

 Has snow/ice/land/ocean databits

344, 2005, 01:10:50 UTC - Cloud mask changed back to original code for ALBEDO ONLY
     based off of testing and observations of "true-color" images, surface albedo is fine for
     clear and probbly clear pixels. The error is in the surface temp algorithm.
     problem likely lies in sza = asin(ReHsat*sin(scanang*DEGRAD)/Rearth); line. NO CHANGE MADE YET - WCS3

347, 2005 - changing sza = asin(ReHsat*sin(scanang*DEGRAD)/Rearth) line to sza=satzen
348, 2005 - changing sza not the answer

EXAMPLE ANTARCTIC RUN: 
	./surftempwcs3 0 ./a1.05200.0711.1000m.hdf ./a1.05200.0711.mod35.hdf

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
float surftemp_clear (float modtb4, float modtb5, float modsatzen, int hemis, short simsk, short landmsk);
float * modis_geo_interp_1000 (int ncols_5, int nrows_5, float *ptr5000);
float zen2scanner(float sza);

/* Main program */

int main(int argc, void *argv[])
{

char *filee1, *filee2;
FILE *fp;
int hemis;
int i;

/* NEED ALLL VARIABLES HERE!!!!*/

int nlines1, nlines2, npixels1, errout, npixels2, looper;
float earth2sun, snang;
float *mod_tb4, *mod_tb5, *tsurf2;
float *mod_solzen5, *mod_satzen5, *mod_relaz5;
uchar8 *mod_cldmsk, *mod_simsk, *mod_landmsk;
//short *mod_cldmsk, *mod_simsk, *mod_landmsk;
float *mod_solzen1,*mod_satzen1, *mod_relaz1;
short *tempsurfbin;

/*set errout to 1 if you want to check for errouneous values*/

errout = 1;

/*get argument values*/

hemis = atoi(argv[1]);
filee1 = argv[2];
filee2 = argv[3];

/* Read in MODIS1B data. Developed by Michael Pavolonis*/
/* Read in Channel 31 and 32 (tb4, tb4 from AVHRR) */

modis_level1b_read(filee1, 5, 31, (void **) &mod_tb4,  &nlines1, &npixels1, &earth2sun, 0);
modis_level1b_read(filee1, 5, 32, (void **) &mod_tb5,  &nlines1, &npixels1, &earth2sun, 0);

/* solzen, satzen, relaz */

modis_level1b_read(filee1, 8, 31, (void **) &mod_solzen5,  &nlines2, &npixels2, &earth2sun, 02);
modis_level1b_read(filee1, 9, 31, (void **) &mod_satzen5,  &nlines2, &npixels2, &earth2sun, 2);
modis_level1b_read(filee1, 10, 31, (void **) &mod_relaz5,  &nlines2, &npixels2, &earth2sun, 2);

/* interpolation of solzen, satzen and relaz. Developed by Michael Pavolonis*/

mod_satzen1= modis_geo_interp_1000(npixels2,nlines2,mod_satzen5);
mod_solzen1= modis_geo_interp_1000(npixels2,nlines2,mod_solzen5);
mod_relaz1= modis_geo_interp_1000(npixels2,nlines2,mod_relaz5);

/* Read in MOD35 cloud mask data. Developed by Michael Pavolonis */

modis_level35_read(filee2, 1, (void **) &mod_cldmsk, &npixels1,&nlines1);

/* Read in MOD35 land mask data. Developed by Michael Pavolonis */

modis_level35_read(filee2, 2, (void **) &mod_landmsk, &npixels1,&nlines1);

/* Read in MOD35 snow/ice data. Developed by Michael Pavolonis */

modis_level35_read(filee2, 3, (void **) &mod_simsk, &npixels1,&nlines1);

/*loop over all pixels*/

looper = npixels1*nlines1;

tsurf2 = (float *) malloc(looper*sizeof(float));
tempsurfbin = (short *) malloc(looper*sizeof(short));

for (i=1; i<looper; i++) {

//IF STATEMENT CONTROLLING LOOP. Cloud mask number must be greater than 1 and data must exist

//CLOUD MASK CONTROL STATEMENT EDITED SO THAT ONLY CLEAR SKY (==3) IS DONE
//this is for testing on errors that have cropped up in surface temp data uncomment other
//if statement to include "probably clear" areas -WCS3 12/9

//EDITED so that data from erroneous values will be stored

       if (mod_cldmsk[i] > 1 && mod_tb4[i] != -999 && mod_tb5[i]!= -999 && mod_satzen1[i] != -999 && mod_tb4[i] < 390 && mod_tb5[i] <390)
	   tsurf2[i]= surftemp_clear(mod_tb4[i],mod_tb5[i],mod_satzen1[i],hemis, mod_simsk[i],mod_landmsk[i] );
       else
		tsurf2[i] = -1.0;

	tempsurfbin[i] = (short)tsurf2[i]*10;

}

/* output data to file. Specified folder must exist*/

printf("calculation complete\n");


fp = fopen("surftemp.dat", "w"); 
fwrite(tempsurfbin,sizeof(short),looper,fp);
fclose(fp);

free(tsurf2);
free(tempsurfbin);

}

/*********************** CLEAR SKY SURFACE TEMPERATURE ********************/

float surftemp_clear (float modtb4, float modtb5, float satzen, int hemis, short simsk, short landmsk)
{
/*--------------------------------------------------------------------------
  This function loops through the input arrays and calls 'surftemp_clear'
  (below) for each pixel.  

  Input (as pointers to the specified type):
    ncols, nrows   : Dimensions of 'tsurf' and the arrays pointed to by 'array' 
                     [short integer].  While these arrays are 2D, they are
                     treated here as 1D.  If scalar values are passed in,
                     both of these values must be 1.  If a 1D array is passed
                     in, then 'nrows' must be 1.
    tb4,tb5        : AVHRR channels 4 and 5 brightness temperatures (K) [float]
    scanang        : Sensor scan angle (degrees) [float]
    hemis          : Hemisphere (0=southern,1=northern) [short integer] 
    surftype       : Surface type as LAND, OCEAN, ICE, SNOW (see 'globals.pro' for 
                     values).  LAND means snow-free land (vegetation cover), 
                     OCEAN means open water, ICE means sea ice (with or without 
                     snow), and SNOW means snow-covered land.
                     However, corrections are currently done for only land
                     and ice surfaces.  Any other surface type uses the ice
                     procedure. [char]
    cmask          : Cloud mask array, cloudy=CLOUDY, clear=CLEAR ('globals.h') 
                     [char]
  Output:
    tsurf          : Surface temperature (K) [float]
                     If badpix=TRUE, MISSING values are returned.  If the 
                     pixel is CLOUDY then the values that this variable
                     has on input is unchanged.

  Additional Variables (from globals.h):
    TRUE, FALSE    : true and false flags to be used for returning error 
                    condition

  Calls: modis_level1b_read (reading the 1B data)

  This is a modification of Jeff Key's surftemp C program by William Straka
  Done for running stand alone in polar regions
  Revision data:6/22/05
----------------------------------------------------------------------------*/

/*  parameters. hemisphere is the only thing passed in. Need to figure this out*/

FILE *fp;

  float  Rearth = 6357.0;            /* Radius of the earth */
  float  Hsat = 850.0;               /* Nominal height of the satellite */
  float  ReHsat = 7207.0;            /* Rearth + Hsat */
  float  degrad = 0.017453292;        /*degrees to radian*/
  float scanang,tsurf,sza;

  short temprange;

//  printf("surftemp_clear: %f, %f, %f, %d", modtb4, modtb5, satzen, hemis);

/********************** CLEAR SKY SURFACE TEMPERATURE ***********************/


/*---------------------------------------------------------------------------- 
  Surface temperature retrieval for clear sky conditions.  
 
  This procedure operates on scalar values only.
  
  Input: 
     tb4       - Channel 4 brightness temperature (K) 
     tb5       - Channel 5 brightness temperature (K) 
     scanang   - scan angle (degrees) 
     hemis     - hemisphere (0=southern,1=northern) 
     satellite - satellite number (7,9,11,12,14,15(K),16(L)) 
     surftype  - A value representing the surface type (global variables
                 LAND, OCEAN, ICE, SNOW see 'globals.h' for values).
                 LAND means snow-free land (vegetation cover), OCEAN means 
                 open water, ICE means sea ice (with or without snow), and
                 SNOW means snow-covered land.
  
  Returned as value of function: 
     tsurf     - Clear sky surface temperature (K)
 
  Called by: call_surftemp_clear  
 
  Algorithm Summary:
 
  The clear sky temperature of ice and snow (and currently open ocean) 
  is estimated as a function of the split-window brightness temperature 
  difference, the channel 4 temperature, and the sensor scan angle.  
  The algorithm uses coefficients for three temperature ranges and is
  based on simulations using LOWTRAN with Arctic Ocean and Antarctic 
  radiosonde data.  The form of the correction equation is:
 
       tsurf = a + b*Tb4 + c*(Tb4-Tb5) + d*(Tb4-Tb5)*(sec(scanang)-1)
 
  where there is a unique set of coefficients a, b, c, and d for the 
  temperature ranges: < 240K, 240-260K, and > 260K.
 
  
---------------------------------------------------------------------------*/

  

/* SNOW-ICE COEFFICIENTS: 
   a,b,c,d[i,j]  where i=0,1 hemisphere: 0=Antarctic,1=Arctic 
                         j=0..2 temp range: 0 = <240, 1 = 240-260, 2 = >260 */

  float a[2][3] =
  {{-0.159480,
    -3.329456,
    -5.207360},

   {-1.571123,
    -2.372697,
    -4.295305}};

  float b[2][3] =
  {{0.999926,
    1.012946,
    1.019429},

   {1.005477,
    1.008604,
    1.015018}};

  float c[2][3] =
  {{1.390388,
    1.214573,
    1.510250},

   {1.853279,
    1.694824,
    1.949525}};

  float d[2][3] =
  {{-0.413575,
    0.131017,
    0.260355},

   {-0.790518,
    -0.205252,
    0.197133}};

//LAND COEFFICIENTS NOAA-15 DATA USED

  float a1[1][3] =
  {{26.0307,
    32.0778,	
    44.1398}};

  float b1[1][3] =
  {{3.9893,
    3.5843,
    3.6687}};

  float c1[1][3] =
  {{-2.9655,
    -2.5600,	
    -2.6554}};

  float d1[1][3] =
  {{-164.3357,
    -165.3660,	
    -182.2059}};

  float e1[1][3] =
  {{132.5978,
    127.3099,	
    134.4543}};

//OCEAN COEFFICENT, Used NOAA-17 coefficient

  float b1sst = 1.01015;

  float b2sst = 2.58150;

  float b3sst = 1.00054;

  float b4sst = 276.590;


	
/* Calculate Zenith to scan angle (adapted from JK pro file */

scanang=zen2scanner(satzen);

/*loop over pixel locations here */

/* Use the channel 4 temperature to get an index for the correct temperature range.  */

  if (modtb4 < 240.)
    temprange = 0;
  else if ((modtb4 >= 240.) && (modtb4 < 260.))
    temprange = 1;
  else /* (modtb4 >= 260)*/
    temprange = 2;

  /* Operate on the appropriate surface type. */
  
 /* SNOW/ICE ALGORITHM */

    if (simsk == 0) {
	    tsurf = a[hemis][temprange] + b[hemis][temprange] * modtb4 +
        	    c[hemis][temprange] * (modtb4 - modtb5) + 
        	    d[hemis][temprange] * (modtb4 - modtb5) * 
        	    ((1.0/cos(scanang*DEGRAD)) - 1.0);
	}
    else {
	    
	if (landmsk == 0) {
    /* OCEAN */
    
//	    sza = asin(ReHsat*sin(scanang*DEGRAD)/Rearth);
	    sza = satzen*DEGRAD;

	    tsurf = b1sst*modtb4 + b2sst*(modtb4 - modtb5) +
	            b3sst*(modtb4 - modtb5)*((1.0/cos(sza)) - 1.0) - 
	            b4sst + 273.16;

	if (tsurf > 285) {
//		fp = fopen ("/home/wstraka/errordata.txt","w");
/*		fp = fopen ("/home/oper/errordata.txt","w");
		fprintf(fp,"Ch31, CH32, Satzen: %f, %f, %f\n", modtb4, modtb5, sza);
		fprintf(fp,"Snow and land mask, temp: %u, %u, %f\n", simsk,landmsk, tsurf);
		fprintf(fp,"Scan angle: %f\n", scanang);
		fclose(fp);  */
         }
		

	}
	else {


    /* SNOW-FREE LAND */
   /* can't do right now, until we figure out how to get emissivities right - 9/12/06, WCS3*/

/* Land emissivities.  These are from the ASTER Spectral Library
   (http://speclib.jpl.nasa.gov), unless indicated otherwise.  */ 

	/*  float lemiss4 = 0.98389; */ /* Grass */
	/*  float lemiss5 = 0.98912; */
	
	/*  float lemiss4 = 0.91210; */  /* Dry grass */
	/*  float lemiss5 = 0.91518; */

	/*  float lemiss4 = 0.98891; */ /* Conifers */
	/*  float lemiss5 = 0.99095; */
	
//	  float lemiss4 = 0.985;    /* From John Collins.  Vegetation type? */
//	  float lemiss5 = 0.975; 


/*	    tsurf = a1[1][temprange] + (b1[1][temprange] * modtb4) + 
	            (c1[1][temprange] * modtb5) + (d1[1][temprange] * lemiss4) + 
	            (e1[1][temprange] * lemiss5); */

	    tsurf = -1.0;

	}
}

return tsurf;

}
