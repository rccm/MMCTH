#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include "modis_grib_readC.h"

int modis_grib_driver_(char* grib_LUN,
                       char *ESDT_name, char *errmsg,
                       int *year, int *month, int *day, int *hour,
                       int ESDT_name_size, int errmsg_size) {
/*
!C**********************************************************************
!Description:
  This subroutine calls modis_grib_setup to write grib data to a binary file.

!Input parameters:
 grib_LUN        Input grib data file logical unit number
 temp_pcf_LUN    Output grib data file logical unit number
 ESDT_name       ESDT name for grib file
 ESDT_name_size  size of ESDT name for grib file
 errmsg_size     size of error message string

!Output Parameters:
 year            year of grib file data (YYYY)
 month           month of grib file data (MM)
 day             day of grib file data (DD)
 hour            hour of grib file data (hh)
 errmsg          Error message string

!Return Value:
 grib_driver returns 0 if exited correctly


!Revision History:
  1997 Apr/May devine(neal.devine@gsfc.nasa.gov)  Original Development.

!Team-unique Header:

     This software is developed by the MODIS Science Data Support
     Team for the National Aeronautics and Space Administration,
     Goddard Space Flight Center, under contract NAS5-32373.

!Design Notes:

!END********************************************************************
*/

  int numofrecords, value;
  int numberofchars;
  grib_record_t *grib_records;
  int GDAS_0ZF, OZ_DAILY, SEA_ICE;
  char temp_pcf_LUN[500];

  /*  The records to read for the met ancillary data  */
  /* rec_num, name, data, OUT: nx,ny, date(year,month,day,hour) */

#ifdef USE_GDAS
	grib_record_t met_records[] = {
      0, "TMP:10 mb:kpds=11,100,10 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:20 mb:kpds=11,100,20 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:30 mb:kpds=11,100,30 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:50 mb:kpds=11,100,50 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:70 mb:kpds=11,100,70 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:100 mb:kpds=11,100,100 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:150 mb:kpds=11,100,150 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:200 mb:kpds=11,100,200 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:250 mb:kpds=11,100,250 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:300 mb:kpds=11,100,300 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:350 mb:kpds=11,100,350 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:400 mb:kpds=11,100,400 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:450 mb:kpds=11,100,450 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:500 mb:kpds=11,100,500 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:550 mb:kpds=11,100,550 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:600 mb:kpds=11,100,600 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:650 mb:kpds=11,100,650 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:700 mb:kpds=11,100,700 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:750 mb:kpds=11,100,750 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:800 mb:kpds=11,100,800 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:850 mb:kpds=11,100,850 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:900 mb:kpds=11,100,900 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:925 mb:kpds=11,100,925 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
  	  0, "TMP:950 mb:kpds=11,100,950 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
  	  0, "TMP:975 mb:kpds=11,100,975 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "TMP:1000 mb:kpds=11,100,1000 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  
	  0, "RH:100 mb:kpds=52,100,100 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:150 mb:kpds=52,100,150 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:200 mb:kpds=52,100,200 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:250 mb:kpds=52,100,250 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:300 mb:kpds=52,100,300 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:350 mb:kpds=52,100,350 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:400 mb:kpds=52,100,400 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:450 mb:kpds=52,100,450 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:500 mb:kpds=52,100,500 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:550 mb:kpds=52,100,550 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:600 mb:kpds=52,100,600 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:650 mb:kpds=52,100,650 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:700 mb:kpds=52,100,700 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:750 mb:kpds=52,100,750 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:800 mb:kpds=52,100,800 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:850 mb:kpds=52,100,850 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:900 mb:kpds=52,100,900 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:925 mb:kpds=52,100,925 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
  	  0, "RH:950 mb:kpds=52,100,950 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
  	  0, "RH:975 mb:kpds=52,100,975 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:1000 mb:kpds=52,100,1000 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  
	  0, "UGRD:10 m above gnd:kpds=33,105,10 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "VGRD:10 m above gnd:kpds=34,105,10 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  	
	  0, "PRES:sfc:kpds=1,1,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1,

 	  0, "TMP:sfc:kpds=11,1,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
	  0, "RH:2 m above gnd:kpds=52,105,2 grid=3", NULL,  -1,-1, -1,-1,-1,-1,

	  0, "PRMSL:MSL:kpds=2,102,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1, 
	  
	  0, "TOZNE:atmos col:kpds=10,200,0 grid=3", NULL, -1, -1, -1, -1, -1, -1,
#ifdef CT_CODE
	  0, "O3MR:10 mb:kpds=154,100,10 grid=3", NULL, -1, -1, -1, -1, -1, -1,
	  0, "O3MR:20 mb:kpds=154,100,20 grid=3", NULL, -1, -1, -1, -1, -1, -1,
	  0, "O3MR:30 mb:kpds=154,100,30 grid=3", NULL, -1, -1, -1, -1, -1, -1,
	  0, "O3MR:50 mb:kpds=154,100,50 grid=3", NULL, -1, -1, -1, -1, -1, -1,
	  0, "O3MR:70 mb:kpds=154,100,70 grid=3", NULL, -1, -1, -1, -1, -1, -1,
	  0, "O3MR:100 mb:kpds=154,100,100 grid=3", NULL, -1, -1, -1, -1, -1, -1
#endif		
  };

#endif

//	  0, "TMP:2 m above gnd:kpds=11,105,2 grid=3", NULL,  -1,-1, -1,-1,-1,-1,

	
#ifdef USE_ECMWF
	

	grib_record_t met_records[] = {
// temperature profile
		97, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		92, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		87, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		82, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		77, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		72, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		67, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		62, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		57, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		52, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		47, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		42, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		37, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		32, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		27, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		22, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		17, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		12, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		7, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		2, NULL, NULL,  -1,-1, -1,-1,-1,-1,
// RH profile		
		98, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		93, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		88, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		83, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		78, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		73, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		68, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		63, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		58, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		53, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		48, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		43, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		38, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		33, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		28, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		23, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		18, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		13, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		8, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		3, NULL, NULL,  -1,-1, -1,-1,-1,-1,
// U and V surface wind components		
		130, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		131, NULL, NULL,  -1,-1, -1,-1,-1,-1,
// surface pressure		
		127, NULL, NULL,  -1,-1, -1,-1,-1,-1,
// SKT and TD2M (not RH2M as in GDAS, so adjust accordingly
		135, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		133, NULL, NULL,  -1,-1, -1,-1,-1,-1,
// These are the O3 levels that are same as GDAS ones		
		205, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		221, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		233, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		253, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		269, NULL, NULL,  -1,-1, -1,-1,-1,-1,
		285, NULL, NULL,  -1,-1, -1,-1,-1,-1, 
// this is the LSM flag, we need it for cloud mask
		134, NULL, NULL,  -1,-1, -1,-1,-1,-1
		
		

// we have no total ozone and no pressure at mean sea level from ECMWF
// so we will use surface pressure and TOAST ozone if we use ECMWF
// there are 20 levels in ECMWF between 10 and 1000 mb
		
	};		
	
#endif
	
	
	
	
/*
 0, "TMP:2 m above gnd:kpds=11,105,2 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
 0, "TMP:sfc:kpds=11,1,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
*/
	
  /*  The record to read for the ozone ancillary data  */
  /* rec_num, name, data, OUT: nx,ny, date(year,month,day,hour) */
  grib_record_t ozn_records[] = {
  0, "TOZNE:atmos col:kpds=10,200,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1};


  /*  The record to read for the ice ancillary data  */
  /* rec_num, name, data, OUT: nx,ny, date(year,month,day,hour) */
  grib_record_t ice_records[] = {
  0, "ICEC:MSL:kpds=91,102,0 grid=255", NULL,  -1,-1, -1,-1,-1,-1};


  GDAS_0ZF = OZ_DAILY = SEA_ICE = 0;

  /*  Zero out the error message string  */
  strncpy((char *) errmsg, " ", errmsg_size);

  if (strncmp( ESDT_name,"GDAS_0ZF",8  ) == 0) {
    GDAS_0ZF = 1;
    grib_records = met_records;
    numofrecords = sizeof(met_records)/sizeof(grib_record_t);

  } else if (strncmp( ESDT_name,"OZ_DAILY",8 ) == 0) {
    OZ_DAILY = 1;
    grib_records = ozn_records;
    numofrecords = sizeof(ozn_records)/sizeof(grib_record_t);
  
  /* rhucek 06/05/03: change number of characters to compare from 8 to 7
     and removed the trailing blank in the literal "SEA_ICE "             
  */
  } else if (strncmp( ESDT_name,"SEA_ICE",7 ) == 0) {
    SEA_ICE = 1;
    grib_records = ice_records;
    numofrecords = sizeof(ice_records)/sizeof(grib_record_t);

  } else {
      sprintf(errmsg, "data product type not supported: %s\n", ESDT_name);
      sprintf(errmsg,"%sOPERATOR ACTION: Contact SDST.", errmsg);
      return (-1);
  } /* end if (strcmp(ESDT_name ...)) */

  /*
   * unpack grib file, write contents to binary format file.
   * log error, if unable to write to binary file.
  */
  
	sprintf(temp_pcf_LUN, "%s.bin", grib_LUN);
  
  
  value = modis_grib_setup (grib_LUN, temp_pcf_LUN, grib_records,
                            numofrecords, errmsg);

  if (value < 0) {
    if (GDAS_0ZF) {
      sprintf(errmsg, "%s - GRIB MET FILE.", errmsg);
    } else if (OZ_DAILY) {
      sprintf(errmsg, "%s - GRIB OZN FILE.", errmsg);
    } else if (SEA_ICE) {
      sprintf(errmsg, "%s - GRIB ICE FILE.", errmsg);
    }

    numberofchars = strlen(errmsg);
    errmsg[numberofchars] = '\0';

    return (-1);
  } /* end if (value < 0) */


  *year = grib_records->date.year;
  *month = grib_records->date.month;
  *day = grib_records->date.day;
  *hour = grib_records->date.hour;

  return (0);
}




int modis_grib_setup (char* inlun, char* outlun, grib_record_t *grib_records,
                int numofrecords, char *errmsg) {

/*
!C**********************************************************************
!Description:
  This subroutine checks the grib files, call modis_grib_read to get the
  grib data, and creates a temporary binary file to write the grib data.

!Input parameters:
 inlun          Input grib data file logical unit number
 outlun         Output grib data file logical unit number
 grib_records   Grib data structure to store data
 numofrecords   Number of records to be written to output file

!Output Parameters:
 errmsg         Error message string

!Return Value:
 modis_grib_setup returns 0 if exited correctly


!Revision History:
  1997 Apr/May devine(neal.devine@gsfc.nasa.gov)  Original Development.
  1997 Dec wolf(walter.wolf@ssec.wisc.edu)  Subroutine modifications.

!References and Credits:
     wgrib v1.5.0 (7-10-96) Wesley Ebisuzaki (NCEP/NCAR Reanalysis Project)

!Team-unique Header:

     This software is developed by the MODIS Science Data Support
     Team for the National Aeronautics and Space Administration,
     Goddard Space Flight Center, under contract NAS5-32373.

!Design Notes:

!END********************************************************************
*/

   int  version;
   int   ret_pgs;
   FILE  *output, *input;

   int i, value;
   int  nx, ny, maxnx, maxny;
   int  records;
   long nxny;
   long len_grib, pos = 0;
   float *array[100];
   unsigned char *buffer, *msg;
   unsigned char *pds, *gds, *bms, *bds;

  /*
   *  This section is used to find the maximum dimensions of
   *  the data arrays in the grib file.  This maximum dimension is
   *  used to allocate the space for each data array read out of the
   *  grib file.
   */

   buffer = (unsigned char *) malloc(1000000);
   version = 1;

	input = fopen(inlun, "r");
   if ( input == NULL ) {
      sprintf(errmsg,"PGS_IO_Gen_Open(lun=%ld,retpgs %d).\nOPERATOR ACTION: Check that file exists",
               inlun, ret_pgs);
      return (-1);
   }

   maxnx = 0;
   maxny = 0;
   pos = 0;
   records = 0;
   msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
   while (msg != NULL) {
      read_grib(input, pos, len_grib, buffer);
      ret_pgs = parse_grib_message( buffer, &pds, &gds, &bms, &bds,
                                   &nx, &ny, &nxny, errmsg);
      if ( ret_pgs != 0) {
         sprintf(errmsg,"modis_grib_setup(): Error: %s", errmsg);
         return (-1);
      }
      records += 1;
      pos += len_grib;
      if (nx > maxnx) maxnx = nx;
      if (ny > maxny) maxny = ny;

      msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
   }
   if (records == 0) {
      sprintf(errmsg,"ERROR: No Records in file.\nOPERATOR ACTION: Contact SDST.");
      return (-1);
   }

	fclose(input);

  /*  Open up the temporary file  */
   version = 1;
   output = fopen(outlun, "w");
   if ( output == NULL ) {
      sprintf(errmsg,"ERROR: PGS_IO_Gen_Temp_Open(lun=%ld,retpgs %d).\nOPERATOR ACTION: Check that file exists",
               outlun, ret_pgs);
      return (-1);
   }

  /*  Allocate space for each data array  */
   for (i = 0; i < numofrecords; i++) {
      if ((array[i] = (float *) malloc(sizeof(float) * maxnx * maxny)) == NULL){
         sprintf(errmsg,"memory allocation problems.\nOPERATOR ACTION: Contact SDST.");
         return(-1);
      }
      grib_records[i].data = array[i];
   }

  /*  Read the data from the grib file  */
   value = modis_grib_read(inlun,1, grib_records, numofrecords, errmsg );
   if (value < 0) return (-1);

  /* Write grib data to single output file. */
   for ( i = 0 ; i < numofrecords ; i++ ) {
      fwrite( (char *)grib_records[i].data, sizeof(float),
                grib_records[i].nx * grib_records[i].ny, output );
      if (ferror(output)) {
         sprintf(errmsg,"error writing %ld rec # %d.\nOPERATOR ACTION: Contact SDST.", outlun, i);
         return(-1);
      }
   }

  /*  Free up the allocated space  */
   free(buffer);
   for (i = 0; i < numofrecords; i++) free(array[i]);

  /*  Close the temporary file  */
	fclose(output);

   return (0);
}
