#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include "PGS_SMF.h"
#include "PGS_IO_Gen.h"
#include "modis_grib_read.h"

int modis_grib_driver_(int *inlun1, int *outlun1, int *inlun2, int *outlun2,
                       int *inlun3, int *outlun3, char *errmsg, 
                       int errmsg_size) {
/*
!C**********************************************************************
!Description:
  This subroutine calls modis_grib_setup to write grib data to a binary file.

!Input parameters:
 inlun1        Input Met grib data file logical unit number
 inlun2        Input Ozone grib data file logical unit number
 inlun3        Input Ice grib data file logical unit number
 outlun1       Output Met grib data file logical unit number
 outlun2       Output Ozone grib data file logical unit number
 outlun3       Output Ice grib data file logical unit number
 errmsg_size   size of error message string

!Output Parameters:
 errmsg        Error message string

!Return Value:
 grib_driver returns 0 if exited correctly


!Revision History:
  06/04/1999 fhliang corrected an error: changed 'outlun' to 'inlun'
  around L.234.

  1997 Dec wolf(walter.wolf@ssec.wisc.edu)  Subroutine modifications.

  1997 Apr/May devine(neal.devine@gsfc.nasa.gov)  Original Development.

!Team-unique Header:
 
!References and Credits:
     wgrib v1.5.0 (7-10-96) Wesley Ebisuzaki (NCEP/NCAR Reanalysis Project)

     External subroutines of importance:  none

     This software is developed by the MODIS Science Data Support
     Team for the National Aeronautics and Space Administration,
     Goddard Space Flight Center, under contract NAS5-32373.

!Design Notes:

PDL:

!END********************************************************************
*/

    int j, numofrecords, value;
    int i, numberofchars;

   /*  The records to read for the met ancillary data  */

    grib_record_t met_records[] = {
    /* rec_num, name, data, OUT: nx,ny, date(year,month,day,hour) */
    0, "TMP:1000 mb:kpds=11,100,1000 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:925 mb:kpds=11,100,925 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:850 mb:kpds=11,100,850 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:700 mb:kpds=11,100,700 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:500 mb:kpds=11,100,500 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:400 mb:kpds=11,100,400 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:300 mb:kpds=11,100,300 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:250 mb:kpds=11,100,250 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:200 mb:kpds=11,100,200 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:150 mb:kpds=11,100,150 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:100 mb:kpds=11,100,100 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:70 mb:kpds=11,100,70 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:50 mb:kpds=11,100,50 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:30 mb:kpds=11,100,30 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:20 mb:kpds=11,100,20 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:10 mb:kpds=11,100,10 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "RH:1000 mb:kpds=52,100,1000 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "RH:925 mb:kpds=52,100,925 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "RH:850 mb:kpds=52,100,850 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "RH:700 mb:kpds=52,100,700 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "RH:500 mb:kpds=52,100,500 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "RH:400 mb:kpds=52,100,400 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "RH:300 mb:kpds=52,100,300 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "LAND:sfc:kpds=81,1,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "TMP:sfc:kpds=11,1,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "PRMSL:MSL:kpds=2,102,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "PWAT:atmos col:kpds=54,200,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "UGRD:10 m above gnd:kpds=33,105,10 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "VGRD:10 m above gnd:kpds=34,105,10 grid=3", NULL,  -1,-1, -1,-1,-1,-1,
    0, "ICEC:sfc:kpds=91,1,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1};


   /*  The record to read for the ozone ancillary data  */

    grib_record_t ozn_records[] = {
    /* rec_num, name, data, OUT: nx,ny, date(year,month,day,hour) */
    0, "TOZNE:atmos col:kpds=10,200,0 grid=3", NULL,  -1,-1, -1,-1,-1,-1};


   /*  The record to read for the ice ancillary data  */

    grib_record_t ice_records[] = {
    /* rec_num, name, data, OUT: nx,ny, date(year,month,day,hour) */
    0, "ICEC:MSL:kpds=91,102,0 grid=255", NULL,  -1,-1, -1,-1,-1,-1};


  /*  Zero out the error message string  */

    strncpy((char *) errmsg, " ", 72);

  /*  This section is set up to read decide which grib file to read  */

   for (j = 1; j < 4; j++) {
      switch (j) {

         case 1: numofrecords = sizeof(met_records)/sizeof(grib_record_t);
                 value = modis_grib_setup (*inlun1, *outlun1, met_records, 
                                           numofrecords, errmsg);
                 if (value < 0) {
                    sprintf (errmsg, "%s: GRIB MET FILE", errmsg);
                    numberofchars = strlen(errmsg);
                    for (i = numberofchars; i < errmsg_size; i++) {
                      errmsg[i] = ' ';
                    }
                    return (-1);
                 }
                 break;

         case 2: numofrecords = sizeof(ozn_records)/sizeof(grib_record_t);
                 value = modis_grib_setup (*inlun2, *outlun2, ozn_records, 
                                           numofrecords, errmsg);
                 if (value < 0) {
                    sprintf (errmsg, "%s: GRIB OZN FILE", errmsg);
                    numberofchars = strlen(errmsg);
                    for (i = numberofchars; i < errmsg_size; i++) {
                      errmsg[i] = ' ';
                    }
                    return (-1);
                 }
                 break;

         case 3: numofrecords = sizeof(ice_records)/sizeof(grib_record_t);
                 value = modis_grib_setup (*inlun3, *outlun3, ice_records, 
                                           numofrecords, errmsg);
                 if (value < 0) {
                    sprintf (errmsg, "%s: GRIB ICE FILE", errmsg);
                    numberofchars = strlen(errmsg);
                    for (i = numberofchars; i < errmsg_size; i++) {
                      errmsg[i] = ' ';
                    }
                    return (-1);
                 }
                 break;

      }
   }
   return (0);
}

int modis_grib_setup (long inlun, long outlun, grib_record_t *grib_records, 
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
 errmsg        Error message string

!Return Value:
 modis_grib_setup returns 0 if exited correctly


!Revision History:
  1997 Apr/May devine(neal.devine@gsfc.nasa.gov)  Original Development.
  1997 Dec wolf(walter.wolf@ssec.wisc.edu)  Subroutine modifications.

!Team-unique Header:

!References and Credits:
     wgrib v1.5.0 (7-10-96) Wesley Ebisuzaki (NCEP/NCAR Reanalysis Project)

     External subroutines of importance:  ww_grib_read.f

     This software is developed by the MODIS Science Data Support
     Team for the National Aeronautics and Space Administration,
     Goddard Space Flight Center, under contract NAS5-32373.

!Design Notes:

PDL:

!END********************************************************************
*/

   PGSt_integer  version;
   PGSt_SMF_status   ret_pgs;
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
   ret_pgs = PGS_IO_Gen_Open( inlun, PGSd_IO_Gen_Read, &input, version);
   if ( ret_pgs != PGS_S_SUCCESS ) {
      sprintf (errmsg,"PGS_IO_Gen_Open(lun=%d,retpgs %d) [OPERATOR ACTION: Check that file exists]",
/* 06/04/1999 fhliang: changed following line
               outlun, ret_pgs);
*/
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
      if ( ret_pgs != PGS_S_SUCCESS ) {
         sprintf (errmsg,"modis_grib_setup: Error: %s", errmsg);
         return (-1);
      }
      records += 1;
      pos += len_grib;
      if (nx > maxnx) maxnx = nx;
      if (ny > maxny) maxny = ny;

      msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
   }
   if (records == 0) {
      sprintf (errmsg,"ERROR: No Records in file [OPERATOR ACTION: Contact SDST]");
      return (-1);
   }

  /*  Close the file  */

   PGS_IO_Gen_Close(input);


  /*  Open up the temporary file  */
  
   version = 1;
   ret_pgs = PGS_IO_Gen_Temp_Open( PGSd_IO_Gen_NoEndurance, outlun, 
                                      PGSd_IO_Gen_Write, &output);
   if ( ret_pgs != PGS_S_SUCCESS ) {
      sprintf (errmsg,"ERROR: PGS_IO_Gen_Temp_Open(lun=%d,retpgs %d) [OPERATOR ACTION: Check that file exists]",
               outlun, ret_pgs);
      return (-1);
   }

  /*  Allocate space for each data array  */

   for (i = 0; i < numofrecords; i++) {
      if ((array[i] = (float *) malloc(sizeof(float) * maxnx * maxny)) == NULL){
         sprintf (errmsg,"memory allocation problems [OPERATOR ACTION: Contact SDST]");
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
         sprintf (errmsg,"error writing %d rec # %d [OPERATOR ACTION: Contact SDST]", outlun, i);
         return(-1);
      }
   }

  /*  Free up the allocated space  */

   free(buffer);
   for (i = 0; i < numofrecords; i++) free(array[i]);

  /*  Close the temporary file  */

   ret_pgs = PGS_IO_Gen_Close(output);
   if ( ret_pgs != PGS_S_SUCCESS ) {
      sprintf (errmsg,"ERROR: PGS_IO_Gen_Close(lun=%d,retpgs %d) [OPERATOR ACTION: Contact SDST]",
               outlun, ret_pgs);
      return (-1);
   }

   return (0);
}
