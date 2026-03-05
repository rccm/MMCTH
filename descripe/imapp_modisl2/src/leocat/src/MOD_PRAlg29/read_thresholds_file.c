/*
!C *********************************************************************
!Description:

  Integer function read_thresholds_file.c
  Reads thresholds file and checks file version.
  Called from get_cm_thresholds.c

!Input arguments:
  none

  Inputs (definitions, constants, etc.) in thresholds.h, including
    threshold file name and current version string

!Output arguments:
  none

  int return_code   successful completion is zero, otherwise non-zero
                    (returned through function call)
                    -1 if file not opened
                    -2 if incorrect version

  Output string array 'lines' contains each line of text found in
    thresholds file.  Defined in thresholds.h.

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "thresholds.h"
#include "granule.h"
#include "stats.h"


/*int read_thresholds_file(char *sfc_dir_name)*/
int read_thresholds_file(char *algData_dir_name) // GPC


{

/*  Declarations */

    FILE *file_ptr;

    char *string_ptr;
    char line[max_line_len];
    char first[id_len];
    char *thrfile_dir_name;
    char *thrfile;
    char *first_id;

    int  icmp;
    int  iver;
    int return_code;

/*  Initializations */

    file_ptr = NULL;
    string_ptr = NULL;
    iver = 0;
    return_code = 0;

/*  Open file */

//  printf("read_thresholds_file alg data dir: %s \n", algData_dir_name);
    thrfile_dir_name = strcat(algData_dir_name, "/");
//  printf("read_thresholds_file thr dir: %s \n", thrfile_dir_name);

    if(sid_terra == 1) {
      thrfile = strcat(thrfile_dir_name, THR_FILE_NAME_TERRA);
      (void) strcpy(g_threshFile, THR_FILE_NAME_TERRA);
      first_id = FIRST_ID_TERRA;
    }
    if(sid_aqua == 1) {
      thrfile = strcat(thrfile_dir_name, THR_FILE_NAME_AQUA);
      (void) strcpy(g_threshFile, THR_FILE_NAME_AQUA);
      first_id = FIRST_ID_AQUA;
    }
    printf("read_thresholds_file thr file name: %s \n", thrfile);

    file_ptr = fopen(thrfile, "r");
    if(file_ptr == NULL) {
      printf("Cannot open file %s\n", thrfile);
      return_code = -1;
    }

    else {

/*    Read file */

      line_count = 0;

      while (1) {

        string_ptr = fgets(lines[line_count], sizeof(lines[line_count]), file_ptr);

        if(string_ptr == NULL)
          break;

/*      Save each line of file in array 'lines' */
        strcpy(line, lines[line_count]);
        line[strlen(line)-1] = '\0';
        strcpy(lines[line_count], line);

/*      Get thresholds version number */
        strncpy(first,line,id_len);
        icmp = strcmp(first,first_id);
        if(icmp == 0) {
          printf("\nVersion of thresholds file ok\n");
          printf("%s\n\n", first);
          iver = 1;
        }

        line_count++;

      }

/*   Close file */
     (void) fclose(file_ptr);

     printf("Found %d lines in file %s\n", line_count, thrfile);

     if(iver == 0) {
       printf("Incorrect version number found in %s\n", thrfile);
       return_code = -2;
     }

   }

   return return_code;

}
