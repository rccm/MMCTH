/*
!C *********************************************************************
!Description:
 
  Integer function check_elem_edge.c
  Determines if current element ('element') is at the border of       
    the input data.
  Called from get_pxldat.c
 
!Input arguments:
  int element       element number in current scan line (0-1353)   
  int n_eles        total number of elements in a scan line (1354)
  int einedge       number elements defined in "edge" region

!Output arguments:
  none

  int loc_elem_edge   value of 1 if on data edge, 0 otherwise
                      (returned through function call)
 
!Revision History:
  R. Frey           05/2007
  
 
!Team-unique Header:
 
!References and Credits:
 
!END ******************************************************************/

/* Includes */

#include <stdio.h>


int check_elem_edge(int element, int n_eles, int einedge)


{



    int loc_elem_edge;

    loc_elem_edge = 0;

    if((element+1) <= einedge || (element+1) > (n_eles - einedge)) {
      loc_elem_edge = 1;
    }

    return loc_elem_edge;

}
