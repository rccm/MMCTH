/*
!C ********************************************************************
!Description:

  Integer function get_pxldat.c

  Collects all inputs neccessary to calculate the confidence of clear
  sky for a given 1-km MODIS pixel.

!Input arguments:
  none

  Input values stored in variables found in granule.h, pixel.h,
  mask_processing_constants.h

!Output arguments:
  none

  Output values stored in structure variable 'pxout' of type "pixel.out".

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>
#include "granule.h"
#include "pixel.h"
#include "mask_processing_constants.h"
#include "stats.h"


int get_pxldat(STATS * stats_out)


{

/*  Declarations */

    extern int check_line_edge(int, int, int);
    extern int check_elem_edge(int, int, int);
    extern int check_reg_uniformity();
    extern void snow_mask();

    unsigned char lsf;
    unsigned char nise_val;

    int return_code;
    int band;
    int New_Zealand;

    long int j;

    float radval;
    float cos_refang, cos_sctang;
    float vrat;
    float elev;

    static float pxlmin_sza=1000.0;
    static float pxlmax_sza=-1000.0;

/*  Initializations */

//  printf("\nExecuting get_pxldat \n");
//  (void)fflush(stdout);

    return_code = 0;
    pxout.nbands_used = 0;

    pxin.line_edge = 0;
    pxin.elem_edge = 0;
    pxin.uniform = 0;
    pxin.land = 0;
    pxin.water = 0;
    pxin.coast = 0;
    pxin.desert = 0;
    pxin.snow = 0;
    pxin.map_snow = 0;
    pxin.ndsi_snow = 0;
    pxin.ice = 0;
    pxin.day = 0;
    pxin.night = 0;
    pxin.visusd = 0;
    pxin.vrused = 1;
    pxin.sunglint = 0;
    pxin.polar = 0;
    pxin.sh_ocean = 0;
    pxin.sh_lake = 0;
    pxin.hi_elev = 0;
    pxin.Greenland = 0;
    pxin.Antarctic = 0;
    pxin.thinCirrusSolar = 0;
    pxin.thinCirrusIR = 0;
    pxin.nonCloudObstruction = 0;

//  Check for data border pixels.
    pxin.line_edge = check_line_edge(pxin.scan, nscans, lines_in_edge);
    pxin.elem_edge = check_elem_edge(pxin.elem, neles, elems_in_edge);

//  Granule index of current pixel.
    j = (pxin.scan * neles) + pxin.elem;

//  Set angular and position data.
    pxin.lat = *(g_lat+j);
    pxin.lon = *(g_lon+j);
    pxin.vza = *(g_vza+j);
    pxin.sza = *(g_sza+j);
    pxin.raz = *(g_raz+j);

//  Save minimum and maximum sza.
    if(rint(pxin.sza) != rint(badsza_data)) {
      if(pxin.sza < pxlmin_sza) pxlmin_sza = pxin.sza;
      stats_out->min_sza = pxlmin_sza;
      if(pxin.sza > pxlmax_sza) pxlmax_sza = pxin.sza;
      stats_out->max_sza = pxlmax_sza;
    }

//  Get ecosystem type (from Olson ecosystem file).
    pxin.eco_type = *(g_eco+j);

//  Calculate glint angle.
    cos_refang = min( (sin(pxin.vza*dtr) * sin(pxin.sza*dtr) * cos(pxin.raz*dtr) +
    		           cos(pxin.vza*dtr) * cos(pxin.sza*dtr) ) , 1.0);
    pxin.refang = acos(cos_refang) * rtd;

//  Calculate scattering angle.
    cos_sctang = -1.0 * (cos(pxin.sza*dtr) * cos(pxin.vza*dtr) -
                         sin(pxin.sza*dtr) * sin(pxin.vza*dtr) * cos(pxin.raz*dtr));
    pxin.sctang = acos(cos_sctang) * rtd;

//  Gather radiance data for current pixel.
    for ( band=0; band<36; ++band) {

      pxin.rad[band] = bad_data;

      if( (g_rad[band]) != NULL) {

        pxout.nbands_used = pxout.nbands_used + 1;

        radval = *((g_rad[band])+j);

        if(band == 25){
    	  if(radval <= 1.3 && radval > -99.0) pxin.rad[band] = radval;
        }
        else if(band == 1) {
          if(radval > 0.0 && radval <= 1.3) pxin.rad[band] = radval;
//        if(radval > 1.3 && pxin.rad[0] <= 2.0) pxin.rad[band] = radval;
          if(radval > 1.3 && pxin.rad[0] <= 2.0) pxin.rad[band] = 1.3;
        }
        else if(band < 19) {
          if(radval > 0.0 && radval <= 1.3)  pxin.rad[band] = radval;
        }
        else if(band == 21) {
          if(radval > 0.0 && radval < 330.0) pxin.rad[band] = radval;
        }
        else {
          if(radval > 0.0 && radval < 1000.0) pxin.rad[band] = radval;
        }
      }
//    printf("rads: %ld %d %f \n", j, band, pxin.rad[band]);
    }

//  Check geolocation. If missing, cloud mask will execute but results will
//  be zeroed out in create_cloud_mask.c.

    pxin.bad_geo = 0;
    if(rint(pxin.lat) == rint(bad_data) || rint(pxin.lon) == rint(bad_data) ||
       rint(pxin.sza) == rint(badsza_data)) pxin.bad_geo = 1;

//  Determine day or night processing.

    if(pxin.sza <= 85.0) {
      pxin.day = 1;
      pxin.visusd = 1;
    }
    else {
      pxin.night = 1;
    }

//  Determine geographic region and surface type processing flags.
    if(pxin.lat > 60.0 || pxin.lat < -60.0) pxin.polar = 1;

//  Determine geometric sun-glint.
    if(pxin.day && pxin.refang <= 36.0) pxin.sunglint = 1;

//  Set processing path flags.

//  Land/sea flag from geolocation file.
    lsf = *(g_lsf+j);

//  Get NDVI background values from "Moody" five-year mean data.
    pxin.ndvibk = *(g_ndvibk+j);

//  Get surface emissivity for bands 22, 31, and 32 (Seeman, Borbas data).
//  pxin.emissfc31 = *(g_emissfc31+j);
//  pxin.emissfc32 = *(g_emissfc32+j);
//  pxin.emissfc22 = *(g_emissfc22+j);

//  Get clear BTs for bands 22, 31, and 32.
//  pxin.bt22clr = *(g_bt22clr+j);
//  pxin.bt31clr = *(g_bt31clr+j);
//  pxin.bt32clr = *(g_bt32clr+j);


//  Force consistency between lsf and ecosystem type for water.
    if(lsf == 0 || (lsf >= 5 && lsf <= 7)) pxin.eco_type = 14;

    if(lsf != 255) {
      if (lsf == 1 || lsf == 4) {
        pxin.land = 1;
        if(pxin.eco_type == 14) {

//        Fix-up for missing ecosystem data in eastern Greenland and
//        north-eastern Siberia.  Without this, these regions become
//        completely "coast".

          if( (pxin.lat < 64.0) ||
              (pxin.lat >= 67.5 && ((pxin.lon < -40.0 && pxin.lon > -168.5) || pxin.lon > -12.5)) ||
              ((pxin.lat >= 64.0 && pxin.lat < 67.5) && ((pxin.lon < -40.0 && pxin.lon > -168.5) ||
              pxin.lon > -30.0)) ) {

              pxin.coast = 1;

          }

        }
      }
      else if (lsf == 2) {
        pxin.coast = 1;
        pxin.land = 1;
      }
      else if (lsf == 3) {
        pxin.land = 1;
        pxin.sh_lake = 1;

//      Need shallow lakes to be processed as "coast" for day, but not night.

        if (pxin.day) pxin.coast = 1;
      }
      else {
        pxin.water = 1;
        if(lsf == 0) pxin.sh_ocean = 1;
      }

//  If land/sea flag is missing, then calculate visible ratio to determine if land or water.

    }
    else if (rint(pxin.rad[0]) != rint(bad_data) &&
             rint(pxin.rad[1]) != rint(bad_data)) {

      vrat = pxin.rad[1] / pxin.rad[0];
      if (vrat > 0.9) {
        pxin.land = 1;
      }
      else {
        pxin.water = 1;
      }

//  If all else fails, call it land.

    }
    else {
      pxin.land = 1;
      pxin.water = 0;
    }

//  Check surface elevation.
//  First, define "Greenland".
    if(pxin.land) {
      if(pxin.lat >= 60.0 && pxin.lat < 67.0) {
        if(pxin.lon >= -60.0 && pxin.lon < -30.0) {
          pxin.Greenland = 1;
        }
      }
      else if(pxin.lat >= 67.0 && pxin.lat < 75.0) {
        if(pxin.lon >= -60.0 && pxin.lon < -10.0) {
          pxin.Greenland = 1;
        }
      }
      else if(pxin.lat >= 75.0) {
        if(pxin.lon >= -70.0 && pxin.lon < -10.0) {
          pxin.Greenland = 1;
        }
      }
    }
    if( (*(g_elev+j) > 2000.0) || (*(g_elev+j) > 200.0 && pxin.Greenland && pxin.land) ||
    	  (pxin.lat >= 75.7 && pxin.lat <= 79.0 && pxin.lon >= -73.0 && pxin.lon <= -50.0 && pxin.land) ) {
      pxin.hi_elev = 1;
    }

//  Define "Antarctic" flag.
    if(pxin.lat < -60.0) pxin.Antarctic = 1;

//  Get 11 um brightness temperature threshold adjustment for deserts.
    pxin.tbadj = ( *(g_elev+j) / 1000.0) * 5.0;

//  Get surface temperature from NWP and SST fields.
    if(pxin.land) {
//    Use GDAS.
      pxin.sfctmp = *(g_sfctmp+j);
    }
    else {
//    Use Reynolds SST.
      pxin.sfctmp = *(g_reynSST+j);
    }
//  printf("sfctmp: %ld %d %d %f %f\n", j,pxin.land,pxin.water,(*(g_reynSST+j)),pxin.sfctmp);

//  Get surface temperature flag that indicate coastal regions where SST test should no be performed.
    pxin.sfct_flag = *(g_sfct_flag+j);

//  Get total precipitable water from NWP fields. 
    pxin.tpw = *(g_tpw+j);

//  Define "desert".

//  Use background NDVI to define "desert".
    if( (pxin.land) && pxin.ndvibk < 0.3) pxin.desert = 1;
    if( pxin.land && pxin.lat < -70.0) pxin.desert = 1;

//  Global.
/*
    if(pxin.eco_type == 8 || pxin.eco_type == 46 ||
       pxin.eco_type == 50 || pxin.eco_type == 51 ||
       pxin.eco_type == 59 || pxin.eco_type == 71 ||
       pxin.eco_type == 11 || pxin.eco_type == 9 ||
       pxin.eco_type == 52) {
       pxin.desert = 1;
    }
//  High-elevation.
    else if( (*(g_elev+j) > 2000.0) && (pxin.eco_type == 42) ) {
      pxin.desert = 1;
//    Check for locations where there should be no high-elevation desert.
      if( (pxin.lat <= 10.0 && pxin.lat >= -10.0) && (pxin.lon >= 90.0) ) {
        pxin.desert = 0;
      }
      else if( (pxin.lat >= -30.0 && pxin.lat <= -10.0) &&
               (pxin.lon >= 160.0 && pxin.lon <= 180.0) ) {
        pxin.desert = 0;
      }
      else if( (pxin.lat >= 10.0 && pxin.lat <= 26.0) &&
               (pxin.lon >= 120.0 && pxin.lon <= 180.0) ) {
        pxin.desert = 0;
      }
    }
//  Africa.
    else if(pxin.lat <= 20.0 && pxin.lat >= -35.0 && pxin.lon <= 60.0 && pxin.lon >= -20.0) {
      if(pxin.eco_type ==  7 || pxin.eco_type == 41 || pxin.eco_type == 43 ||
    	 pxin.eco_type == 58 || pxin.eco_type == 36 || pxin.eco_type == 91 ||
         pxin.eco_type == 32 || pxin.eco_type == 29) {
         pxin.desert = 1;
      }
    }
//  Eurasia.
    else if(pxin.lat <= 70.0 && pxin.lat >= -60.0 && pxin.lon <= 180.0 && pxin.lon >= -20.0) {
      if(pxin.eco_type == 11 || pxin.eco_type == 2) {
        pxin.desert = 1;
      }
//    Exclude New Zealand.
      if(pxin.lat >= -50.0 && pxin.lat <= -30.0 && pxin.lon >= 160.0 && pxin.lon <= 180.0) {
        pxin.desert = 0;
      }
    }
//  Add in Australia.
    if(pxin.lat <= -11.0 && pxin.lat > -40.0 && pxin.lon <= 155.0 && pxin.lon > 110.0) {
      if(pxin.eco_type == 43 || pxin.eco_type == 41 || pxin.eco_type == 91) {
        pxin.desert = 1;
      }
    }
*/

//  Determine whether or not visible ratio test may be used over land surfaces.
    if(pxin.eco_type == 2 || pxin.eco_type == 8 || pxin.eco_type == 11 || pxin.eco_type == 40 ||
       pxin.eco_type == 41 || pxin.eco_type == 46 || pxin.eco_type == 51 || pxin.eco_type == 52 ||
       pxin.eco_type == 59 || pxin.eco_type == 71 || pxin.eco_type == 50) {
       pxin.vrused = 0;
    }

//  Get NISE snow cover.
    nise_val = *(g_nise+j);
    if (nise_val == 103 || nise_val == 104 || nise_val == 101 ||
       (nise_val > 25 && nise_val < 102 && pxin.water) ||
       (nise_val > 25 && nise_val < 101 && pxin.land) ||

       (pxin.day && nise_val == 200 && pxin.land && (pxin.lat >= 44.0 || pxin.lat <= -40.0)) ||

       (pxin.night && nise_val == 200 && pxin.land && (pxin.lat >= 60.0 || pxin.lat <= -40.0 ||
       (pxin.lat >= 44.0 && (pxin.lon > 30.0 || pxin.lon < -10.0))) ) ||

       (nise_val == 200 && pxin.water && pxin.lat < -60.0) ||

       (pxin.day && nise_val == 200 && pxin.water && pxin.lat >=44.0) ||

       (pxin.night && nise_val == 200 && pxin.water && (pxin.lat >= 60.0 || 
       (pxin.lat >= 44.0 && (pxin.lon > 30.0 || pxin.lon < -10.0))) ) ) {

      pxin.map_snow = 1;
    }

    if(pxin.day) {
//    Run quick version of D.Hall's snow detection algorithm.
      (void) snow_mask();

      if(pxin.water) {

        if(pxin.lat >= -60.0 && pxin.lat <= 25.0) {
          if(pxin.map_snow && pxin.ndsi_snow) {
            pxin.ice = 1;
          }
        }
        else if( (nise_val == 252 || nise_val == 200 || nise_val == 253) &&
        		 (pxin.lat >=44.0 || pxin.lat <= -40.0) ) {
          if(pxin.ndsi_snow) {
            pxin.ice = 1;
          }
        }
        else if( (lsf == 3 || lsf == 5) && pxin.ndsi_snow) {
          pxin.ice = 1;
        }
        else if(pxin.map_snow && pxin.ndsi_snow) {
          pxin.ice = 1;
        }

      }
      else if(pxin.land) {

        if(pxin.lat >= -60.0 && pxin.lat <= 25.0) {
//        Define New Zealand region which receives snow but snow map does not show it.
          if( (pxin.lat >= -48.0 && pxin.lat <= -34.0) && (pxin.lon >= 165.0 && pxin.lon <= 180.0) ) {
            New_Zealand = 1;
          }
          else {
            New_Zealand = 0;
          }
          if(pxin.map_snow && pxin.ndsi_snow) {
            pxin.snow = 1;
          }
          else if(New_Zealand && pxin.ndsi_snow) {
            pxin.snow = 1;
          }
        }
        else if(pxin.lat < -60.0) {
          pxin.snow = 1;
        }
        else {
          if(pxin.ndsi_snow) {
            pxin.snow = 1;
          }
        }

      }

    }
    else {

      pxin.ice = pxin.map_snow;
      pxin.snow = pxin.map_snow;
      elev = *(g_elev+j);

      if(pxin.map_snow && pxin.sfctmp > 280.00 && elev < 500.0) {
        pxin.ice = 0;
        pxin.snow = 0;
      }

      if(pxin.lat > 86.0) pxin.ice = 1;
      if(pxin.land && pxin.lat < -60.0) pxin.snow = 1;

    }

//  printf("\nlat/lon:       %d %d %ld %f %f \n", pxin.scan, pxin.elem, j, pxin.lat, pxin.lon );
//  printf("angles:          %f %f %f %f %f \n", pxin.vza, pxin.sza, pxin.raz, pxin.refang, pxin.sctang );
//  printf("Day/night:       %d %d %d \n", pxin.day, pxin.night, pxin.visusd);
//  printf("Polar/sun-glint: %d %d \n", pxin.polar, pxin.sunglint);
//  printf("Processing path: %d %d %d %d %d %d %d \n", lsf, pxin.eco_type,pxin.land,pxin.water,pxin.coast,pxin.sh_ocean,pxin.sh_lake);
//  printf("hi_elev, tbadj:  %d %d %d %f \n", pxin.Antarctic, pxin.Greenland, pxin.hi_elev, pxin.tbadj);
//  printf("desert, vrused, ndvi:  %d %d %d %f \n", pxin.eco_type, pxin.vrused, pxin.desert, pxin.ndvibk);
//  printf("surface temp:    %f %d %d %d %d %d \n", pxin.sfctmp, nise_val, pxin.map_snow, pxin.ndsi_snow, pxin.snow, pxin.ice);


//  Check regional uniformity
    pxin.uniform = check_reg_uniformity();
//  if(pxin.uniform == 0) {
//    printf("non-uniformity: %d %d %d %d %d\n", pxin.scan, pxin.elem, pxin.line_edge, pxin.elem_edge, pxin.uniform);
//  }


    return return_code;

}
