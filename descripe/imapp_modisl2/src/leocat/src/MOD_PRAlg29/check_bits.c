/*
!C *********************************************************************
!Description:

  Function check_bits.c
  Function for checking if a bit has been set in a 6-byte array
    variable.

!Input arguments:
  int bit_num       0-based bit number to check (0-47)
  char *testbits    pointer to 6-element cloud mask byte array

!Output arguments:
  none

  Function returns 1 if bit is set, 0 if not set (through returned
    variable 'bit_test'.

!Revision History:
  R. Frey           05/2007

!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>


int check_bits(int bit_num, unsigned char *testbits)


{

/*  Declarations */

    unsigned char cmask;

    int byte_num;
    int bit_test;
    int pos;

/*  Initialize output. */
    bit_test = -1;
    cmask = 0;

/*  Determine which byte (0-5) of the array contains the bit of
    interest. */
    byte_num = bit_num / 8;

/*  Determine the position of the bit within the current byte. */
    pos = bit_num - (byte_num * 8);

    cmask |= (0x01<<pos);
    if( (cmask & *(testbits + byte_num)) > 0) {
      bit_test = 1; }
    else if( (cmask & *(testbits + byte_num)) == 0) {
      bit_test = 0; }
    else {
      printf("Problem in check_bits.c: bit_test = -1\n");
    }
    
//  printf("bit_test: %d %d %d %d %d %d \n", bit_num, byte_num, pos, cmask, *(testbits + byte_num), bit_test);

    return bit_test;

}
