/*
!C *********************************************************************
!Description:

  Function clear_bit.c
  Clears a single bit (=0) within the 48-bit cloud mask.

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


void clear_bit(int bit_num, char *testbits)


{

/*  Declarations */

    char cmask = 0;

    int byte_num;
    int pos;

/*  Determine which byte (0-5) of the array contains the bit of
    interest. */
    byte_num = bit_num / 8;

/*  Determine the position of the bit within the current byte. */
    pos = bit_num - (byte_num * 8);

/*  Clear the bit. */
    cmask |= (1<<pos);
    *(testbits + byte_num) = ~cmask & *(testbits + byte_num);

    return;

}
