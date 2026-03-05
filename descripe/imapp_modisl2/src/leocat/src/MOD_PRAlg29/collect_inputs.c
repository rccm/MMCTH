/*
!C ********************************************************************
!Description:

  Integer function collect_inputs.c
  Collects all inputs necessary for executing the MODIS cloud mask.
  Called from modis_cm_main.c

!Input arguments:
  none

  Inputs from ...

!Output arguments:
  none

  int return_code   successful completion is zero, otherwise non-zero
                    (returned through function call)

  Output parameters in ...

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*int collect_inputs(char *sfc_dir_name)*/
int collect_inputs(char *algData_dir_name) // GPC


{

/*  Declarations */

    extern int get_cm_thresholds(char *);

    int ret_code;
    int return_code;

/*  Initializations */

    return_code = 0;

/*  Get cloud mask thresholds */

    /*ret_code = get_cm_thresholds(sfc_dir_name); // GPC*/
    ret_code = get_cm_thresholds(algData_dir_name);
    ret_code = 0;

    if(ret_code != 0) {
      printf("Problem collecting cloud mask input data\n");
      return_code = -1;
    }









    return return_code;

}
