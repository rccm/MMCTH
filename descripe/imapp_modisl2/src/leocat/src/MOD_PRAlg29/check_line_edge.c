/*
!C *********************************************************************
!Description:
 
  Integer function check_line_edge.c
  Determines if current scan line ('scan_line') is at the border of 
    the input data.
  Called from create_cloud_mask.c
 
!Input arguments:
  int scan_line     scan line number of the data (0-based)
  int n_scans       total number of scan lines in the data
  int linedge       number lines defined to be in "edge" region
 
!Output arguments:
  none

  int loc_line_edge     value of 1 if on data edge, 0 otherwise
                        (returned through function call)
 
!Revision History:
  R. Frey           05/2007

 
!Team-unique Header:
 
!References and Credits:
 
!END ******************************************************************/


/* Includes */

#include <stdio.h>


int check_line_edge(int scan_line, int n_scans, int linedge)

{

    int loc_line_edge;

    loc_line_edge = 0;

    if((scan_line+1) <= linedge || (scan_line+1) > (n_scans - linedge)) {
      loc_line_edge = 1;
    }

    return loc_line_edge;

}
