/*
!C *********************************************************************
!Description:

  Integer function get_thr_values.c    
  Returns values of a particular cloud test threshold.
  Called from get_cm_thresholds.c 

!Input arguments:
  char name         threshold name as found in thresholds file
  int nvals         number of values associated with threshold

  Additional inputs (definitions, constants, etc.) in thresholds.h

!Output arguments:
  float *thr_ptr    pointer to threshold value(s)

  int return_code   successful completion is zero, otherwise non-zero
                    (returned through function call)

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "thresholds.h"

 
int get_thr_values(char name[], int nvals, float *thr_ptr)


{

/*  Declarations */
  
    char line[max_line_len];
    char first[max_thrname_len];
    const char delm[] = " :,";
    char *value;
    
    float thr;

    int jj;
    int kk;
    int return_code;
    int icmp;

/*  Initializations */

    jj = 0;
    return_code = 0;

/*  Get threshold value(s) for input threshold name. */

    while (jj < line_count) {

      strncpy(first,lines[jj],strlen(name));
      first[strlen(name)] = '\0';

      icmp = strcmp(first,name);
      if(icmp == 0) {

        strcpy(line,lines[jj]);
        value = strtok(line, delm);

        kk = 0;
        while(kk < nvals) {

          value = strtok(NULL, delm);
          if(value == NULL)
            break;

          sscanf(value, "%f", &thr);
          *(thr_ptr + kk) = thr;
          kk++;

        }
    
        break;
      }

      jj++;

    }

    if(jj == line_count) {
      printf("Problem getting threshold values for %s\n", name);
      return_code = -1;
    }

    return return_code;

}
