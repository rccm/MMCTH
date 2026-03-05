/*
!C ********************************************************************
!Description:

  Function snow_mask

  Subroutine which implements a quick and dirty version of the snow
  algorithm by Dorothy Hall, George Riggs and Vince Salmonson.

!Input arguments:
  none

  Input values stored in variables found in pixel.h, granule.h, 
  mask_processing _constants.h, and thresholds.h.

!Output arguments:
  none

  Output stored in 'pxin.ndsi_snow' found in pixel.h.

!Revision History:
  R. Frey           10/2007  Converted to C

!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>
#include "granule.h"
#include "pixel.h"
#include "mask_processing_constants.h"
#include "thresholds.h"

void snow_mask()

{ 

// Declarations
   int I_bad;
   
   float m04;
   float m02;
   float m26;
   float m20;
   float m29;
   float m31;
   float m35;
   float nirband;
   float ndsi;
   float diff;
   float diff2;
   float sth85_11;
   float sth37_11;

// Initializations
   m04 = pxin.rad[3];
   m02 = pxin.rad[1];
   m26 = pxin.rad[25];
   m20 = pxin.rad[19];
   m29 = pxin.rad[28];
   m31 = pxin.rad[30];
   m35 = pxin.rad[34];

   diff = 0.0;
   diff2 = 0.0;
   
// Select NIR channel for NDSI based on platform name.
// Band 6 for Terra, band 7 for Aqua.
   if(sid_terra) {
     nirband = pxin.rad[5];
   }     
   else if(sid_aqua) {
     nirband = pxin.rad[6];
   }  
   else {
     printf("Platform name not recognized in snow_mask \n");
   }

   I_bad = rint(bad_data);
   
// Begin snow/ice detection.
   
// First, check the 11 micron brightness temperature.
   if (rint(m31) != I_bad && m31 <= sm_bt11[0]) {

//   Get NDSI index
     if (rint(m04) != I_bad && rint(nirband) != I_bad) {
       ndsi  = (m04 - nirband) / (m04 + nirband);
//     printf("NDSI value is: %f %f %f %f %f %f \n", ndsi, sm_ndsi[0], m04, nirband, m02, sm_ref2[0]);

       if (ndsi > sm_ndsi[0] && m02 > sm_ref2[0]) { 
         pxin.ndsi_snow = 1;

//       Check for false snow detection.
//       Make sure NDSI is not flagging a thin cirrus cloud.
         if ( !pxin.Greenland && pxin.land) {
           if (rint(m26) != I_bad && rint(m35) != I_bad) {
             if (m26 > sm_ref3[0] && m35 < sm_co2[0]) {
               pxin.ndsi_snow = 0;
             }
           }
         }
         
//       If in sunglint region, disregard if between -60 and 50 degrees latitude.
         if(pxin.water && pxin.sunglint && pxin.lat <= 50.0 && pxin.lat >= -60.0) {
           pxin.ndsi_snow = 0;
         }
         
//       Check for ice clouds mis-identified as snow. Note modified BTD thresholds for Greenland.
         if(rint(m29) != I_bad && rint(m31) != I_bad) {
           diff = m29 - m31;
           sth85_11 = sm85_11[0];
           if(rint(m20) != I_bad && !(pxin.Greenland && pxin.land)) {
             diff2 = m20 - m31;
             sth37_11 = sm37_11[0];
             if(diff >= sth85_11 && diff2 >= sth37_11) {
               pxin.ndsi_snow = 0;
             }
           }  
           else {
             if(pxin.Greenland && pxin.land) {
               sth85_11 = sm85_11[0] + 1.5;
             }  
             else {
               sth85_11 = sm85_11[0];
             }
             if(diff >= sth85_11) {
               pxin.ndsi_snow = 0;
             }
           }
//         printf("ndsi check ice: %d %d %f %f %f %f %f \n", pxin.Greenland, pxin.land, diff, sth85_11, m20, diff2, sth37_11);
         }
           
//       Check for water clouds mis-identified as snow.
         if(m20 != I_bad && m31 != I_bad) {
           if( !(pxin.land && pxin.Greenland) ) {
             diff = m20 - m31;
             if(pxin.hi_elev) {
               sth37_11 = sm37_11hel[0];
             }  
             else {
               sth37_11 = sm37_11[0];
             }
             if(diff >= sth37_11) {
               pxin.ndsi_snow = 0;
             }
           }
//         printf("ndsi check water: %f %f \n", diff, sth37_11);
         }

//       Check value of NIR reflectance.         
         if( !pxin.hi_elev && !(pxin.Greenland && pxin.land) ) {
           if(nirband > sm_mnir[0]) {
             pxin.ndsi_snow = 0;
           }
//         printf("ndsi check NIR refl: %f %f \n", nirband, sm_mnir[0]);
         }
        

       } 
     }
   }
    

}
   
