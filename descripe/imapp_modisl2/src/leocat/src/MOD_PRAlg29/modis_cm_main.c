/*
!C *********************************************************************
!Description:

  MODIS Cloud Mask main program (modis_cm_main.c)
  Determines if a given pixel is obstructed by cloud or thick aerosol.
  Taken from MOD35 FORTRAN code.

!Input arguments:

!Output arguments:

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void modis_cm_main(char *algData_dir_name) // GPC

{

    extern int collect_inputs(char *);
    extern int create_cloud_mask();

    int irt;

    /*irt = collect_inputs(sfc_dir_name);*/
    irt = collect_inputs(algData_dir_name); // GPC

    irt = create_cloud_mask();

//  printf("Returning from modis_cm_main\n");
//  (void)fflush(stdout);

}
