/*
!C *********************************************************************
!Description:

  Routine for filling in snow and ice information in coastal
  regions.  For all coast pixels in the NISE data, search both
  latitudinally and longitudinally for non-coast pixels which
  indicate snow or ice.  Stop search when non-snow/ice, non-coast
  pixels are found.  For those coast pixels which have snow/ice
  pixels indicated on their boundaries, set flag to indicate 
  snow/ice.

!Input Parameters:
  h_index     Hemisphere index; 1=NH, 2=NH
  map_nise    Input/output snow/ice data grid (granule.h)
  xsize       Size of gris in horizontal direction
  ysize       Size of grid in vertical direction

!Output Parameters:
  map_nise    Input/output snow/ice data grid (granule.h)

!Revision History:
  Original subroutine:
      03-09-04  R. Frey
      10-23-07  R. Frey converted to C

!Team-Unique Header:

!References and Credits:
  See Cloud Mask ATBD-MOD-06.

!End ******************************************************************/

#include <stdio.h>
#include <math.h>
#include "granule.h"

void massage_snowice(int h_index, int xsize, int ysize)

{
	
/* Declarations */
	
   int k[4];
   int row;
   int col;
   int m;
   long int idx;	
   unsigned char nise_out;
   unsigned char *nise_in;
   unsigned char n;

// printf("Entering massage_snowice \n");
   
// Set pointer to appropriate data.   
   if(h_index == 1) {
     nise_in = g_nisenorth;
   }
   else {
	 nise_in = g_nisesouth;
   }
   
// Loop through map data.

   idx = 0;
   
   for ( row=0; row<ysize; ++row ) {
     for ( col=0; col<xsize; ++col ) {

       nise_out = 0;
     
//     printf("processing: %d %d %ld %d \n", row, col, idx, *(nise_in+idx));
       
//     Check for coast pixel (value = 252).
       if( *(nise_in+idx) == 252) {

//       Check if on map edges.  Cannot search on all four sides.
//       No data in corners.

         if(col == 0) {
//         k[0] = ichar(map_nise(j+1,i))
           k[0] = *(nise_in + (idx+1) );
//         k[1] = ichar(map_nise(j,i-1))
           k[1] = *(nise_in + (idx-xsize) );
//         k[2] = ichar(map_nise(j,i+1))
           k[2] = *(nise_in + (idx+xsize) );
           k[3] = 0;
         }
         else if(row == 0) {
//         k(1) = ichar(map_nise(j,i+1))
           k[0] = *(nise_in + (idx+xsize) );
//         k(2) = ichar(map_nise(j+1,i))
           k[1] = *(nise_in + (idx+1) );
//         k(3) = ichar(map_nise(j-1,i))
           k[2] = *(nise_in + (idx-1) );
           k[3] = 0;
         }  
         else if(col == (xsize-1) ) {
//          k(1) = ichar(map_nise(j-1,i))
        	k[0] = *(nise_in + (idx-1) );
//          k(2) = ichar(map_nise(j,i-1))
        	k[1] = *(nise_in + (idx-xsize) );
//          k(3) = ichar(map_nise(j,i+1))
        	k[2] = *(nise_in + (idx+xsize) );
            k[3] = 0;
         }   
         else if(row == (ysize-1) ) {
//         k(1) = ichar(map_nise(j,i-1))
           k[0] = *(nise_in + (idx-xsize) );
//         k(2) = ichar(map_nise(j+1,i))
           k[1] = *(nise_in + (idx+1) );
//         k(3) = ichar(map_nise(j-1,i))
           k[2] = *(nise_in + (idx-1) );
           k[3] = 0;
           }
           else {
//          Get all four adjacent values.
//          k(1) = ichar(map_nise(j,i-1))
        	k[0] = *(nise_in + (idx-xsize) );
//          k(2) = ichar(map_nise(j,i+1))
        	k[1] = *(nise_in + (idx+xsize) );
//          k(3) = ichar(map_nise(j-1,i))
        	k[2] = *(nise_in + (idx-1) );
//          k(4) = ichar(map_nise(j+1,i))
        	k[3] = *(nise_in + (idx+1) );
           }

//         Fill in "missing" snow/ice values (=200).
           for (m=0; m<4; ++m) {

             n = k[m];
             if(n != 0) {

               if( (n == 103 || n == 104 || n == 200) || (n > 25 && n < 102) ) {
                 nise_out = 200;
//               printf("processing: %d %d %ld %d %d \n", row, col, idx, *(nise_in+idx), nise_out);
               }

             }
           }

           if (nise_out == 200) *(nise_in+idx) = nise_out;

       }       

       idx++;
       
     }
   }



}
