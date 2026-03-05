
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include "PGS_SMF.h"
#include "PGS_IO_Gen.h"
#include "modis_grib_read.h"

int modis_grib_read(PGSt_Logical lun, PGSt_integer version,
		     grib_record_t *grib_records, int N_grib_records,
                     char *errmsg )

/*
!C**********************************************************************
!Description:
  This function opens a grib file with given lun & version and calls
  wgrib based code to decode the specified records along with date & hour
  and nx & ny for grid.

!Input Parameters: 
  PGSt_Logical  lun            -- LogicalUnitNumber of input grib file.
  PGSt_integer  version        -- version of input grib file.
  grib_record_t *grib_records  -- Struct specifying grib data to decode:
    int   rec_num;    -- if > 0 specifies grib rec # to decode.
    char  *name;      -- specifies inventory of grib record in 
                         format based on "wgrib -v" or -V:
                         "*kpds=%d,%d,%d grid=%d*" where the %d's in order are:
                         pds octet  9   : parameter & units
		         pds octet 10   : indicator of type of level or layer.
		         pds octet 11&12: Height, pressure, etc. of level/layer.
		         pds octet  7   : Grid Id. (from wgrib -V, 1st line)
      Example: "UGRD:1000 mb:kpds=33,100,1000 grid=3" 
                (can include informative info not used by code)
    float *data;       -- pointer to buffer for decoded data.  
                          If NULL, data set to malloced mem.
    int   nx;          -- OUT: nx for grid.
    int   ny;          -- OUT: ny for grid.
    grib_date_t date;  -- OUT struct grib_date_t: 
                          {  int   year,month,day,hour; }


  int           N_grib_records -- Number of grib records to decode from file.

!Output Parameters: 
  errmsg   Error message character string
  see data, nx, ny, date above.

!Return Value: 
  modis_grib_read returns 0 on Success

!Revision History:
  1997 Apr/May devine(neal.devine@gsfc.nasa.gov)  Original Development.
  1997 Dec wolf(walter.wolf@ssec.wisc.edu)       Subroutine update

!Team-unique Header:
     This software is developed by the MODIS Science Data Support
     Team for the National Aeronautics and Space Administration,
     Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
     wgrib v1.5.0 (7-10-96) Wesley Ebisuzaki (NCEP/NCAR Reanalysis Project)

     External subroutines of importance:  none

!Design Notes:

PDL: not today!

!END********************************************************************
*/
{
  PGSt_SMF_status   ret_pgs;
  unsigned char *buffer, *msg;
  double        temp;
  int           irec,nrec_read, /* prev_rec_num = 0, */ ierror, nrec_in_file;
  int           is_required_reading;
  long          len_grib, pos = 0, nxny, buffer_size;
  unsigned char *pds, *gds, *bms, *bds;
  FILE          *input;

 /*  Allocate temporary buffer space  */

  if ((buffer = (unsigned char *) malloc(BUFF_ALLOC0)) == NULL){
      sprintf (errmsg,"modis_grib_read:  memory allocation problems [OPERATOR ACTION: Contact SDST]");
      return(-1);
  }
  buffer_size = BUFF_ALLOC0;

 /*  Open the grib file for reading  */

  ret_pgs = PGS_IO_Gen_Open( lun, PGSd_IO_Gen_Read, &input, version);
  if ( ret_pgs != PGS_S_SUCCESS ) {
     sprintf (errmsg,
          "modis_grib_read:  ERROR: PGS_IO_Gen_Open(lun=%d,retpgs %d) [OPERATOR ACTION: Check that file exists]",
          lun, ret_pgs);
     return (-1);
  }

  /*
   * Find and read all grib records in file requested by user.
   */

  nrec_in_file = 0;
  for ( nrec_read = 0 ; nrec_read < N_grib_records ;) {

    /*
     * Seek the next grib record in file requested by user.
     */ 

    for ( is_required_reading = FALSE ; !is_required_reading ; ) {
      msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
      if (msg == NULL) {
         sprintf (errmsg,
          "modis_grib_read: ERROR: Reached end of file before data was found [OPERATOR ACTION: Contact SDST]");
         return (-1);
      }
      nrec_in_file++;

      /*
       *  Check this grib record for match against any requested, 
       *  set is_required_reading if so.
       */
      for ( irec = 0 ; irec < N_grib_records ; irec++ ) {

        /* 
         *  If a specified grib record number matches current record in 
         *  file, read this record later.
         *  If name info was specified, verify it against PDS. 
         */

        if ( grib_records[irec].rec_num == nrec_in_file ) {
          is_required_reading = TRUE;
          if ( grib_records[irec].name != NULL && 
               !Name_Matches_PDS(grib_records[irec].name,(char *)(msg+8)) ) {
             sprintf (errmsg,
            "modis_grib_read: ERROR: specified rec number & name inconsistent [OPERATOR ACTION: Contact SDST]");
             return (-1);
          }
        }

       /* 
        *  If a specified grib record number is invalid (<1) and 
        *  specified name info matches PDS, read this record later. 
        */

        if ( grib_records[irec].rec_num < 1 && 
             Name_Matches_PDS(grib_records[irec].name,(char *)(msg+8)) ) {
          grib_records[irec].rec_num = nrec_in_file;
          is_required_reading = TRUE;
        }
      }
      if ( !is_required_reading ) pos += len_grib;

    }

    /*
     * Read full grib record into buffer starting from start of grib record.  
     * realloc buffer if too small.
     */

     if ( len_grib+1 > buffer_size ) {
       buffer_size = len_grib + 1;
       free( buffer );
       buffer = (unsigned char *) malloc( buffer_size );
       if (buffer == NULL) {
          sprintf (errmsg,"modis_grib_read: memory allocation problems [OPERATOR ACTION: Contact SDST]");
          return(-1);
       }
     }

    /*  Read the data from the grib file  */

     read_grib(input, pos, len_grib, buffer);
     pos += len_grib;

    /*
     * For each requested record == current record in file:
     *   parse grib message, decode numeric data & date.
     */

     for ( irec = 0 ; irec < N_grib_records ; irec++ ) 
       if (grib_records[irec].rec_num == nrec_in_file) {

      /*
       * parse grib message
       */

      ierror = parse_grib_message( buffer, &pds, &gds, &bms, &bds, 
				   &grib_records[irec].nx, 
                                   &grib_records[irec].ny, &nxny, errmsg);
      if ( ierror != GRIB_SUCCESS ) {
         sprintf (errmsg,"modis_grib_read: %s [OPERATOR ACTION: Contact SDST]", errmsg);
         return(-1);
      }

      /* malloc space for data if no pointer provided. */

      if ( grib_records[irec].data == NULL ) {
	if ( (grib_records[irec].data = 
              (float *) malloc( nxny * sizeof(float) )) == NULL ) {
          sprintf (errmsg,"modis_grib_read: memory allocation problems [OPERATOR ACTION: Contact SDST]");
          return(-1);
        }
      }

      /*
       * decode numeric data
       */

      temp = int_power(10.0, - PDS_DecimalScale(pds));
      BDS_unpack(grib_records[irec].data, bds + 11, BMS_bitmap(bms), 
                 BDS_NumBits(bds), nxny, temp*BDS_RefValue(bds),
                 temp*int_power(2.0, BDS_BinScale(bds)));

      /*
       * Read date from pds.
       */

      grib_records[irec].date.year = PDS_Year4(pds);
      grib_records[irec].date.month = PDS_Month(pds);
      grib_records[irec].date.day  = PDS_Day(pds);
      grib_records[irec].date.hour = PDS_Hour(pds);

      nrec_read++;
    }

  }/* while nrec_read < N_grib_records */

 /*  Free up the allocated space  */

  free(buffer);

  ret_pgs = PGS_IO_Gen_Close(input);
  if ( ret_pgs != PGS_S_SUCCESS ) {
     sprintf (errmsg,"ERROR: PGS_IO_Gen_Close(lun=%d,retpgs %d) [OPERATOR ACTION: Contact SDST]",
              lun, ret_pgs);
     return (-1);
  }
  return( 0 );
}


int Name_Matches_PDS( char *name, char *pds ) 
/*
!C**********************************************************************

!Description:
  Determine whether info in name in format "...kpds=%d,%d,%d grid=%d..."
  is same as that in PDS in octets 9, 10, 11&12, and 7.  These define
  the parameter, units, height, and grid type which should be enough
  to fully specify the desired record.

!Input Parameters: 
  char *name  -- string in format "...kpds=%d,%d,%d grid=%d..."
  char *pds   -- pointer to pds in grib record.

!Output Parameters: N/A

!Return Value:
  (int) 1 if name matches pds, 0 if not.

!Revision History:
  97-10-15 N. Devine -- Original, using wgrib macros.

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and credits:
    N. Devine (devine@gsfc.nasa.gov)

!Design Notes:

!END*************************************************************
*/
{
  int pds7,pds9,pds10,pds11to12;
  char * kpds_in_name = strstr(name,"kpds=") + 5;
  char * grid_in_name = strstr(name,"grid=") + 5;

/*
 *  Due to problems with defined grid number defined, a change has been made.
 *  The new if statement checks the grid number of the record read from
 *  the grib file.  If the number is greater than 126 and less than 255
 *  (not between 201 and 216) then we do not check the grid number of
 *  the requested record verses the record read from the file.  We still
 *  check the other variables.  (Note that this if statement assumes that
 *  the maximum grid number is 255.)
 *
 *  Changed by Walter Wolf  (Univ of Wisconsin-Madison)  01/23/98)
 */

  if ((PDS_Grid(pds) > 126 && PDS_Grid(pds) < 201) || 
      (PDS_Grid(pds) > 216 && PDS_Grid(pds) <= 255)) {
     if ( kpds_in_name-5 != NULL && grid_in_name-5 != NULL  &&
          sscanf(kpds_in_name,"%d,%d,%d",&pds9,&pds10,&pds11to12) == 3  && 
          pds9 == PDS_PARAM(pds)  &&  pds10 == PDS_KPDS6(pds)  &&  
          pds11to12 == PDS_KPDS7(pds)  &&
          sscanf(grid_in_name,"%d",&pds7) == 1 ) {
       return ( TRUE );
     }else{
       return ( FALSE );
     }
  }
  else {
     if ( kpds_in_name-5 != NULL && grid_in_name-5 != NULL  &&
          sscanf(kpds_in_name,"%d,%d,%d",&pds9,&pds10,&pds11to12) == 3  && 
          pds9 == PDS_PARAM(pds)  &&  pds10 == PDS_KPDS6(pds)  &&  
          pds11to12 == PDS_KPDS7(pds)  &&
          sscanf(grid_in_name,"%d",&pds7) == 1  &&  pds7 == PDS_Grid(pds) ) {
       return ( TRUE );
     }else{
       return ( FALSE );
     }
  }
}



/* ********************************************************************** */
/* ********************************************************************** */
#ifndef min
   #define min(a,b)  ((a) < (b) ? (a) : (b))
#endif

#define NTRY 100
/* #define LEN_HEADER_PDS (28+42+100) */
#define LEN_HEADER_PDS (28+8)

unsigned char *seek_grib(FILE *file, 
			 long *pos, 
			 long *len_grib, 
			 unsigned char *buffer, 
			 unsigned int buf_len) 
/*
!C**********************************************************************
!Description:
  find next grib header.
!Input Parameters:
  FILE *file            -- input grib file pointer (what do you think?)
  long *pos             -- Position from file start to start reading from ( = 0 for 1st call).
  unsigned int buf_len  -- Size of buffer
  unsigned char *buffer -- Buffer used to read grib record into.
!Output Parameters:
  long *pos             -- Position of 'GRIB' in file, or pos + buf_len - 
                             LEN_HEADER_PDS if 'GRIB' not found in NTRY tries.
  long *len_grib        -- Length of grib record as recorded after 'GRIB'.
  unsigned char *buffer -- Buffer used to read grib record into.
!Return Value:
  (char *) to start of GRIB header+PDS, NULL if not found


!Revision History:
 * adapted from SKGB (Mark Iredell)
 *
 * v1.1 9/94 Wesley Ebisuzaki
 * v1.2 3/96 Wesley Ebisuzaki handles short records at end of file
 * v1.3 8/96 Wesley Ebisuzaki increase NTRY from 3 to 100 for the folks
 *      at Automation decided a 21 byte WMO bulletin header wasn't long 
 *      enough and decided to go to an 8K header.  

!Team-unique Header:

  wgrib is a portable program to read grib files that were created by
  the NCEP/NCAR Reanalysis Project.  Of course, the program is not 
  restricted to Reanalysis files but Eugenia Kalnay is happy whenever
  she sees the phrase "NCEP/NCAR Reanalysis".  

!References and credits:
  adapted from SKGB (Mark Iredell)

!Design Notes:

!END*************************************************************
*/
{

    int i, j, len;

    for (j = 0; j < NTRY; j++) {

        if (fseek(file, *pos, SEEK_SET) == -1) {
            *len_grib = 0;
            return (unsigned char *) NULL;
        }
     
        i = fread(buffer, sizeof (unsigned char), buf_len, file);
     
        len = i - LEN_HEADER_PDS;
     
        for (i = 0; i < len; i++) {
            if (buffer[i] == 'G' && buffer[i+1] == 'R' && buffer[i+2] == 'I'
                && buffer[i+3] == 'B' && buffer[i+7] == 1) {
                    *pos = i + *pos;
                    *len_grib = (buffer[i+4] << 16) + (buffer[i+5] << 8) +
                            buffer[i+6];
                    return (buffer+i);
            }
        }
	*pos = *pos + (long) (buf_len - LEN_HEADER_PDS);
    }

    *len_grib = 0;
    return (unsigned char *) NULL;
}



/* ********************************************************************** */
/* ********************************************************************** */
int read_grib( FILE          *file, 
		      long          pos, 
		      long          len_grib, 
		      unsigned char *buffer)
/*
!C**********************************************************************
!Description:
  Read full GRIB record of length len_grib into buffer.

!Input Parameters:
  FILE          *file
  long          pos      -- File position from start to begin reading from.
  long          len_grib -- Number of chars to read from file.

!Output Parameters:
  unsigned char *buffer  -- Buffer to read grib record into

!Return Value:
  (int) (i == len_grib) -- if false, error has occured.

!Revision History:
 * 1997/10 N. Devine -- Added prolog.
 * v1.0 9/94 Wesley Ebisuzaki

!Team-unique Header:

  wgrib is a portable program to read grib files that were created by
  the NCEP/NCAR Reanalysis Project.  Of course, the program is not 
  restricted to Reanalysis files but Eugenia Kalnay is happy whenever
  she sees the phrase "NCEP/NCAR Reanalysis".  

!References and credits:
    Wesley Ebisuzaki
    N. Devine (devine@gsfc.nasa.gov)

!Design Notes:

!END*************************************************************
*/
{

    int i;

    if (fseek(file, pos, SEEK_SET) == -1) {
	    return 0;
    }

    i = fread(buffer, sizeof (unsigned char), len_grib, file);
    return (i == len_grib);
}


/* ********************************************************************** */
int parse_grib_message( unsigned char *buffer,
			unsigned char **pds_p,
			unsigned char **gds_p,
			unsigned char **bms_p,
			unsigned char **bds_p, 
			int *nx, int *ny, long *nxny, char *errmsg)
/* parse grib message 
!C**********************************************************************
!Description:
  parse grib message for pointers to sections (pds,gds,bms,bds), and determine
  grid size and check for GRIB record end.

!Input Parameters:
  unsigned char *buffer -- Buffer containing grib record.

!Output Parameters:
  unsigned char *pds_p    -- Pointer to Product Definition Sector.
  unsigned char *gds_p    -- Pointer to Grid Definition Sector.
  unsigned char *bms_p    -- Pointer to Bit Map Sector (NULL if not present).
  unsigned char *bds_p    -- Pointer to Binary Data Sector.
  (Note these values are assigned at end, local char* values are otherwise
   used to avoid errors in macros which reference pds[], etc.)

!Return Value:
  int   -- GRIB_SUCCESS
        -- GRIB_MISSING_END_SECTION   , probably fatal.
        -- GRIB_MISSING_GDS_NUM_DATAPOINTS , not fatal?

!Revision History:
 * 11/94 - v1.0  Wesley Ebisuzaki
 * ...
 * 1997/10 N. Devine -- Removed this code from main(), removed error messages.

!Team-unique Header:

  wgrib is a portable program to read grib files that were created by
  the NCEP/NCAR Reanalysis Project.  Of course, the program is not 
  restricted to Reanalysis files but Eugenia Kalnay is happy whenever
  she sees the phrase "NCEP/NCAR Reanalysis".  

!References and credits:
    Wesley Ebisuzaki
    N. Devine (devine@gsfc.nasa.gov)

!Design Notes:

!END*************************************************************
*/
{
        unsigned char *msg, *pds, *gds, *bms, *bds, *pointer;

	msg = buffer;
        pds = (msg + 8);
        pointer = pds + PDS_LEN(pds);

        if (PDS_HAS_GDS(pds)) {
            gds = pointer;
            pointer += GDS_LEN(gds);
        }
        else {
            gds = NULL;
        }

        if (PDS_HAS_BMS(pds)) {
            bms = pointer;
            pointer += BMS_LEN(bms);
        }
        else {
            bms = NULL;
        }

        bds = pointer;
        pointer += BDS_LEN(bds);

        /* end section - "7777" in ascii */
        if (pointer[0] != 0x37 || pointer[1] != 0x37 ||
            pointer[2] != 0x37 || pointer[3] != 0x37) {
           sprintf (errmsg,"parse_grib_message: Missing End Section [OPERATOR ACTION: Contact SDST]");
           return (-1);
        }

	/* figure out size of array */
	if (gds != NULL) {
	    /* this doesn't work for spherical harmonics */
	    GDS_grid(gds, nx, ny, nxny);
	}
	else if (bms != NULL) {
	    *nxny = *nx = BMS_nxny(bms);
	    *ny = 1;
	}
	else {
	    if (BDS_NumBits(bds) == 0) {
                *nxny = *nx = 1;
               sprintf (errmsg,
                        "parse_grib_message: Missing GDS, constant record [OPERATOR ACTION: Contact SDST]");
               return (-1);
	    }
	    else {
	        *nxny = *nx = BDS_NValues(bds);
	    }
	    *ny = 1;
	}

	*pds_p = pds;
	*gds_p = gds;
	*bms_p = bms;
	*bds_p = bds;
	return( GRIB_SUCCESS );
}


/* ********************************************************************** */
int GDS_grid(unsigned char *gds, int *nx, int *ny, long int *nxny)
/*
!C**********************************************************************

!Description:
    Determine number of x & y grid points, and product, from grib GDS.

!Input Parameters: 
  unsigned char *gds -- pointer to GDS in grib record.

!Output Parameters: 
  int      *nx       -- number of x dim grid points
  int      *ny       -- number of x dim grid points
  long int *nxny     -- nx * ny

!Revision History:
  11/94 - v1.0  Wesley Ebisuzaki
  97-10-15 N. Devine -- Added prolog.

!Team-unique Header:

  wgrib is a portable program to read grib files that were created by
  the NCEP/NCAR Reanalysis Project.  Of course, the program is not 
  restricted to Reanalysis files but Eugenia Kalnay is happy whenever
  she sees the phrase "NCEP/NCAR Reanalysis".  

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    int i, ix, iy, pl;
    long int isum;

    *nx = ix = GDS_LatLon_nx(gds);
    *ny = iy = GDS_LatLon_ny(gds);
    *nxny = ix * iy;

    /* thin grid */
    if (GDS_Gaussian(gds) | GDS_LatLon(gds)) {
	if (ix == 65535) {
	    *nx = -1;
	    /* reduced grid */
	    isum = 0;
	    pl = GDS_PL(gds);
	    for (i = 0; i < iy; i++) {
		isum += gds[pl+i*2]*256 + gds[pl+i*2+1];
	    }
	    *nxny = isum;
	}
	return 0;
    }
    return 0;
}


/* ********************************************************************** */
double int_power(double x, int y)
/*
!C**********************************************************************

!Description:
    Determine number of x & y grid points, and product, from grib GDS.

!Input Parameters: 
  double x
  int    y

!Output Parameters: None

!Return Value:
  (double) x**y

!Revision History:
  11/94 - v1.0  Wesley Ebisuzaki
  97-10-15 N. Devine -- Added prolog.

!Team-unique Header:

  wgrib is a portable program to read grib files that were created by
  the NCEP/NCAR Reanalysis Project.  Of course, the program is not 
  restricted to Reanalysis files but Eugenia Kalnay is happy whenever
  she sees the phrase "NCEP/NCAR Reanalysis".  

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

	double value;

	if (y < 0) {
		y = -y;
		x = 1.0 / x;
	}
	value = 1.0;

	while (y) {
		if (y & 1) {
			value *= x;
		}
		x = x * x;
		y >>= 1;
	}
	return value;
}

/* ********************************************************************** */
/* ********************************************************************** */
static unsigned int mask[] = {0,1,3,7,15,31,63,127,255};

void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
		       int n_bits, int n, double ref, double scale)
/*
!C**********************************************************************
!Description: 
  Converts packed grid-data in grib BinaryDataSection to array of float.

!Input Parameters:
  unsigned char *bits   -- BDS octet 12.
  unsigned char *bitmap -- Bit Map (BMS octet 7), or NULL if BMS not defined.
  int           n_bits  -- Nbits/datum, BDS octet 11.
  int           n       -- Number of data to write to flt.
  double        ref     -- (10^-D)*R
  double        scale   -- (10^-D)*2^E
    R = reference value     , BDS octet 7-10
    E = binary scale factor , BDS octet 5-6
    D = decimal scale factor, PDS octet 27-28
    All octet numbers above are 1 based, not 0

!Output Parameters:
  float         *flt    -- decoded grib data.

!Return Value: N/A

!Revision History:
  1996	   wesley ebisuzaki -- Original
  1996/04 -- v1.1 faster
  97-10-15 N. Devine -- Added prolog.

!Team-unique Header:

  wgrib is a portable program to read grib files that were created by
  the NCEP/NCAR Reanalysis Project.  Of course, the program is not 
  restricted to Reanalysis files but Eugenia Kalnay is happy whenever
  she sees the phrase "NCEP/NCAR Reanalysis".  

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

    int t_bits, c_bits, j_bits;
    unsigned int i, j, map_mask, tbits, jmask, bbits;
    long int jj;

    tbits = bbits = map_mask = 0;

    /* assume integer has 32+ bits */
    if (n_bits <= 25) {
        jmask = (1 << n_bits) - 1;
        t_bits = 0;

        if (bitmap) {
            while (n-- > 0) {
	        if (bitmap) {
		    if (map_mask == 0) {
			map_mask = 128;
			bbits = *bitmap++;
		    }
	            if ((bbits & map_mask) == 0) {
		        *flt++ = UNDEFINED;
	                map_mask >>= 1;
		        continue;
	            }
	            map_mask >>= 1;
	        }
	        while (t_bits < n_bits) {
	            tbits = (tbits * 256) + *bits++;
	            t_bits += 8;
	        }
	        t_bits -= n_bits;
	        j = (tbits >> t_bits) & jmask;
	        *flt++ = ref + scale*j;
            }
        }
        else {
	    for (i = 0; i < n; i++) {
                while (t_bits < n_bits) {
                    tbits = (tbits * 256) + *bits++;
                    t_bits += 8;
                }
                t_bits -= n_bits;
	        j = (tbits >> t_bits) & jmask;
                flt[i] = (tbits >> t_bits) & jmask;
            }
	    /* at least this vectorizes :) */
	    for (i = 0; i < n; i++) {
		flt[i] = ref + scale*flt[i];
	    }
        }
    }
    else {
        c_bits = 8;
        while (n-- > 0) {
	    if (bitmap) {
	        j = (*bitmap & map_mask);
	        if ((map_mask >>= 1) == 0) {
		    map_mask = 128;
		    bitmap++;
	        }
	        if (j == 0) {
		    *flt++ = UNDEFINED;
		    continue;
	        }
	    }

	    jj = 0;
	    j_bits = n_bits;
	    while (c_bits <= j_bits) {
	        if (c_bits == 8) {
		    jj = (jj << 8) + *bits++;
		    j_bits -= 8;
	        }
	        else {
		    jj = (jj << c_bits) + (long)(*bits & mask[c_bits]);
		    bits++;
		    j_bits -= c_bits;
		    c_bits = 8;
	        }
	    }
	    if (j_bits) {
	        c_bits -= j_bits;
	        jj = (jj << j_bits) + (long)((*bits >> c_bits) & mask[j_bits]);
	    }
	    *flt++ = ref + scale*jj;
        }
    }
    return;
}

/* ********************************************************************** */
double ibm2flt(unsigned char *ibm)
/* 
!C**********************************************************************
!Description:
  Called by BDS_RefValue() macro to convert BDS octet 7-10 into a double.

!Input Parameters:
  unsigned char *ibm -- pointer to BDS octet 7-10 in grib record.

!Output Parameters: N/A

Return Value:
  (double) value stored in BDS octet 7-10.

!Revision History:
  1996	v1.1   wesley ebisuzaki -- Original
  97-10-15 N. Devine -- Added prolog.

!Team-unique Header:

  wgrib is a portable program to read grib files that were created by
  the NCEP/NCAR Reanalysis Project.  Of course, the program is not 
  restricted to Reanalysis files but Eugenia Kalnay is happy whenever
  she sees the phrase "NCEP/NCAR Reanalysis".  

!References and credits:

!Design Notes:

!END*************************************************************
*/
{

	int positive, power;
	unsigned int abspower;
	long int mant;
	double value, exp;

	positive = (ibm[0] & 0x80) == 0;
	mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
	power = (int) (ibm[0] & 0x7f) - 64;
	abspower = power > 0 ? power : -power;


	/* calc exp */
	exp = 16.0;
	value = 1.0;
	while (abspower) {
		if (abspower & 1) {
			value *= exp;
		}
		exp = exp * exp;
		abspower >>= 1;
	}

	if (power < 0) value = 1.0 / value;
	value = value * mant / 16777216.0;
	if (positive == 0) value = -value;
	return value;
}






