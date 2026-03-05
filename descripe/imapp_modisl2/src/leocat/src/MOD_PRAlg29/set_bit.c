/*
!C *********************************************************************
!Description:

  Function set_bit.c
  Sets a single bit within the 48-bit cloud mask.

!Input arguments:
  int bit_num       0-based bit number to set (0-47)

!Output arguments:
  char *testbits    pointer to 6-element cloud mask byte array

!Revision History:
  R. Frey           05/2007

!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>


void set_bit(int bit_num, unsigned char *testbits)


{

/*  Declarations */

    int byte_num;
    int pos;

/*  Determine which byte (0-5) of the array contains the bit of
    interest. */
    byte_num = bit_num / 8;

/*  Determine the position of the bit within the current byte. */
    pos = bit_num - (byte_num * 8);

/*  Set the bit. */
    *(testbits + byte_num) |= (0x01<<pos);

    return;

}
