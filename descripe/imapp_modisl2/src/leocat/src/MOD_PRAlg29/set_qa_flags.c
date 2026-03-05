/*
!C *********************************************************************
!Description:

  Function set_qa_flags.c
  Sets QA flags for a 1-km MODIS cloud mask pixel.

!Input arguments:
  int num_tests_performed   number spectral tests performed

!Output arguments:
  unsigned char *qabits    pointer to 10-element cloud mask QA byte array

!Revision History:
  R. Frey           02/2010

!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include "pixel.h"


void set_qa_flags(int num_tests_performed, unsigned char *qabits)


{

/*  Declarations */

    extern void set_bit(int, unsigned char[]);

//  Set cloud mask usefulness (bit 0) and confidence (bits 1-2) QA flags.

    if(num_tests_performed > 0) {

      (void) set_bit(0, pxout.qabits);

      if(num_tests_performed < 3) {
        (void) set_bit(3, pxout.qabits);
      }
      else if(num_tests_performed < 7 || pxin.sunglint) {
        (void) set_bit(2, pxout.qabits);
        (void) set_bit(3, pxout.qabits);
      }
      else {
        (void) set_bit(1, pxout.qabits);
        (void) set_bit(2, pxout.qabits);
        (void) set_bit(3, pxout.qabits);
      }

    }

//  Number spectral  used

    if(pxout.nbands_used > 14) {
      (void) set_bit(48, pxout.qabits);
      (void) set_bit(49, pxout.qabits);
    }
    else if (pxout.nbands_used > 7) {
      (void) set_bit(49, pxout.qabits);
    }
    else if (pxout.nbands_used > 0) {
      (void) set_bit(48, pxout.qabits);
    }

//  Number spectral tests used

    if(num_tests_performed > 6) {
      (void) set_bit(50, pxout.qabits);
      (void) set_bit(51, pxout.qabits);
    }
    else if (num_tests_performed > 3) {
      (void) set_bit(51, pxout.qabits);
    }
    else if (num_tests_performed > 0) {
      (void) set_bit(50, pxout.qabits);
    }

//  Ancillary data

//  Clear Radiance Origin
    (void) set_bit(56, pxout.qabits);
    (void) set_bit(57, pxout.qabits);

//  Surface Temperature Over Land
//  qa bits 58, 59 == 0 -> GDAS used

//  Surface Temperature Over Ocean
//  qa bits 60, 61 == 0 -> Reynolds used

//  Surface Winds (not used but indicate "other")
    (void) set_bit(62, pxout.qabits);
    (void) set_bit(63, pxout.qabits);

//  Ecosystem Map (indicate Olson)
    (void) set_bit(64, pxout.qabits);

//  Snow Mask (indicate "other")
    (void) set_bit(66, pxout.qabits);
    (void) set_bit(67, pxout.qabits);

//  Ice Cover (indicate "other")
    (void) set_bit(68, pxout.qabits);
    (void) set_bit(69, pxout.qabits);

//  Land/Sea Mask
//  qa bits 70, 71 == 0 -> USGS 1-km 6-level used

//  DEM
//  qa bit 72 == 0 -> EOS DEM

//  Precipitable Water
//  qa bits 73, 74 == 0 -> NCEP GDAS

    return;

}
