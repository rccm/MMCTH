#include <stdio.h>
#include <stdlib.h>
#include "mapic.h" 

int createMODISarray(MODFILE *file, char const *arrayname, 
                     char const *groupname, char const *data_type, 
                     long int rank, long int const dimsizes[])
/* 
!C****************************************************************************
*!Purpose:     Initializes an array in the MODIS-HDF file to store a multi-
*             dimensional array of data.
*
*!Description: Function createMODISarray is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
*	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      function which allow MODIS Science Team-supplied software
*	      to read in Level 1B radiance bands and write out output 
*	      products and metadata to HDF files.  The functionality of the
*	      M-API is defined in the MODIS API User's Guide, Version 1, 
*	      dated 6/16/95.
*
*	      createMODISarray creates an HDF SDS structure to store 
* 	      a new data array into a MODIS HDF file. It must be called before
*             the data may be written to the file using putMODISarray (PMAR)
*	      or the associated labels for each array dimension may (optionally)
*	      be stored using putMODISarlabel (PMARLB).
*
*	      The groupname string provides the facility to place the new array
*	      in an HDF 'Vgroup' data group. If a Vgroup with the name groupname
*	      does not exist, the array structure will not be created. The array
*	      may be placed in the file outside of any Vgroup by replacing
*	      groupname with NULL in C and a blank string (" ") in FORTRAN.
*
*	      If an array with the name arrayname is written outside of a Vgroup,
*	      it must not already exist in the file. This is to prevent the
*	      confusion caused by multiple data objects within the same file.
* 	      SDS's with the same name are permitted, however, if they are 
*	      placed in separate Vgroups.
*
* !Input Parameters:	MODFILE *file 	: MODISfile name
*			char *arrayname	: input array name
*			char *groupname : data group name in which to put array
*					  ( if NULL, don't put in any group )
*			char *data_type : data type of the array  (int8, uint8,
*                                         int16, uint16, int32, uint32,
*                                         float32, float64)
*			long int rank	: dimensions of the input array
*			long int dimsizes: the size of each array dimension
*
* !Output Parameters:	NONE
*
*Returns:               MAPIOK if successful or MFAIL if an error occurs.
*
*Externally references: MODFILE 				(mapi.h)
*			PGS_SMF_MAX_MSGBUF_SIZE         	(PGS_SMF.h)
*			MAX_VAR_DIMS				(netcdf.h)
*			NULLstr					(mapic.h)
*			MAPIERR					(mapic.h)
*			NULLMODFIL				(mapic.h)
*			DFACC_READ				(hdf.h)
*			SDnametoindex				(mfhdf.h)
*			NO_OBJECT				(mapic.h)
*
*			datatype_to_DFNT 			(mapic.h)
*			SDcreate				(mfhdf.h)
*			SDidtoref				(
*			searchMODISgroup			(mapi.h)
*                       SDendaccess				(mfhdf.h) 
*			addMODISgroup				(mapic.h)
*                       FAIL            			(hdf.h)
*			NULL					(stdio.h)
*			MFAIL					(mapi.h)
*			MAPIOK					(mapi.h)
* !Revision History:
*		Qi Huang	1996/01/24
*		Version 2.0
*
*$Log: createMODISarray.c,v $
*Revision 5.1  2005/04/04 18:17:57  vlin
*constant safe for pointer arguments.
*
*Revision 1.1  1998/02/06 22:26:06  fshaw
*Initial revision
*
 * Revision 1.9  1997/05/30  17:56:26  qhuang
 * Modified data_type in input parameters section in prolog.
 *
 * Revision 1.8  1997/02/05  17:39:47  qhuang
 * Added for loop check to dimsizes.
 *
 * Revision 1.6  1996/08/09  15:47:00  qhuang
 * Removed some ','s in error messages, casted some argument in HDF routine
 * call.
 *
 * Revision 1.5  1996/07/31  17:21:14  qhuang
 * Version 2.1.  Ring super structure and other changes make this version
 * much faster.
 *
 * Revision 1.4  1996/05/08  17:14:52  qhuang
 * Added NO_OBJECT ihcluded in external references, included programmer's
 * name, limited length of arrayname and groupname strings in ERROR message,
 * cast dimension size movement from dimsizes to sds_dimsizes, cast sd_id
 * in SD* function calls, retrieved sds_ref for SDS using SDidtoref, corrected
 * searchMODISgroup,SDendaccess and SDidtoref ERROR message.
 *
 * Revision 1.3  1996/01/29  22:17:53  qhuang
 * Added sds_id for SDS using SDidtoref. Used sds_id to store in Vgroup
 * as argument to addMODISgroup.
 *
 * Revision 1.2  1996/01/24  16:52:44  qhuang
 * Original development and testing for mapi version2.0.
 *
 * Revision 1.1  1996/01/24  16:51:24  qhuang
 * Initial revision
 *
*
* !Team-unique Header:
*             Portions developed at the National Center for Supercomputing
*             Applications at the Univ. of Illinois at Urbana-Champaign.
*
* !References and Credits:
*              This software is developed by the MODIS Science Data Support
*	      Team for the National Aeronautics and Space Administration,
*	      Goddard Space Flight Center, under contract NAS5-32373.
*
* !Design Notes:
*
!END*********************************************************************
*/
{
  /* buffer to hold the error/warning message */
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE]; 
  char *funcname="createMODISarray";    /* name of this routine */	
  int32   sds_dimsizes[MAX_VAR_DIMS];
  int32   ret_search, number_type, sds_id, sds_ref;
  int     i, status_code;

  status_code = MFAIL;

  /* Input checks */
  if ( NULLstr(arrayname) ){
    sprintf(buff, "ERROR: createMODISarray unable to make a new\n"
			"\tarray without an array name input\n");
    MAPIERR(buff, funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) ){
    sprintf(buff, "ERROR: createMODISarray unable to make a new %.*s\n"
			"\t array with a NULL file MODFILE structure\n",
			  MAX_NC_NAME,arrayname);			
    MAPIERR(buff, funcname);
    return(MFAIL);
  }

  if ( dimsizes == NULL ){
    sprintf(buff,"ERROR: createMODISarray unable to make a new %.*s\n"
			"\t array without array dimension input\n",
			  MAX_NC_NAME,arrayname);
    MAPIERR(buff, funcname);
    return(MFAIL);
  }

  for (i=0; i<rank; i++)
    if (dimsizes[i] <= 0)
      {
        sprintf(buff,"ERROR: createMODISarray unable to create a new\n"
                "\t %.*s array with %ld dimsize.\n",
                              MAX_NC_NAME,arrayname, dimsizes[i]);
        MAPIERR(buff, funcname);
        return(MFAIL);
      }

  if ( data_type == NULL || NULLstr(data_type)  ){
    sprintf(buff,"ERROR: createMODISarray unable to make a new %.*s\n"
			"\t array without an array data type input\n",
			  MAX_NC_NAME,arrayname);
    MAPIERR(buff, funcname);
    return(MFAIL);
  }

  if ( (rank < 1) || (rank > MAX_VAR_DIMS) ){
    sprintf(buff,"ERROR: createMODISarray unable to create a new\n"
			"\t %.*s array with %ld dimensions.\n",
			MAX_NC_NAME,arrayname,rank);
    MAPIERR(buff, funcname);
    return(MFAIL);
  }

  for (i=0; i< rank; i++)
    sds_dimsizes[i] = (int32)dimsizes[i];

  if ( file->access == DFACC_READ ){
    sprintf(buff,"ERROR: createMODISarray unable to make new\n"
		"\t%.*s arrayname in file opened for read only\n",
			MAX_NC_NAME,arrayname);
    MAPIERR(buff, funcname);
    return(MFAIL);
  }

  if ( NULLstr(groupname) ) 
     groupname = NULL;
  
  if ( groupname != NULL){
    if ( (ret_search = searchMODISgroup(file,groupname,NULL,arrayname,
	NULL,DFTAG_NDG)) == MFAIL ){
       sprintf(buff,"ERROR: createMODISarray detected an error in\n"
                    "\t searchMODISgroup while searching the\n"
                    "\t %.*s array.\n",MAX_NC_NAME,arrayname);
       MAPIERR(buff, funcname);
       return(MFAIL);
    }else if (ret_search != NO_OBJECT){
       sprintf(buff,"ERROR: createMODISarray found %.*s array\n"
                    "\t already exists in data group %.*s\n",
                    MAX_NC_NAME,arrayname,VGNAMELENMAX,groupname);
       MAPIERR(buff,funcname);
       return(MFAIL);
    }
  }else if ( SDnametoindex( (int32) file->sd_id,arrayname) != MFAIL ){
  /* then array already exists. */
      sprintf(buff,"ERROR: createMODISarray found %.*s array already\n"
                   "\t exists.\n",MAX_NC_NAME,arrayname);
      MAPIERR(buff,funcname);
      return(MFAIL);
  }

  if ((number_type = datatype_to_DFNT(data_type)) == MFAIL ){
    sprintf(buff,"ERROR: createMODISarray unable to create %.*s\n"
			"\t array of data type %s\n",MAX_NC_NAME,
                        arrayname,data_type);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }
 
  if ((sds_id =  SDcreate( (int32)file->sd_id,arrayname,number_type,(int32)rank,sds_dimsizes)) == MFAIL ){
    sprintf(buff,"ERROR: createMODISarray detected FAIL from HDF\n"
			"\t procedure SDcreate attempting to create\n"
			"\t %.*s array\n",MAX_NC_NAME,arrayname);
    MAPIERR(buff,funcname);
  }

  else{

      if ( groupname != NULL ){

         if ( (sds_ref = SDidtoref(sds_id))  == MFAIL ){
             sprintf(buff,"ERROR: createMODISarray detected FAIL from\n"
			"\t HDF procedure SDidtoref attempting\n"
			"\t to obtain the reference id of\n"
			"\t %.*s array\n",MAX_NC_NAME,arrayname);
	     MAPIERR(buff,funcname);
         }

         else if ( addMODISgroup(file,groupname,NULL,DFTAG_NDG,sds_ref) == MFAIL){
             sprintf(buff,"ERROR: createMODISarray unable to add\n"
			"\t %.*s array to data group\n"
			"\t %.*s because of failure in\n"
			"\t addMODISgroup.\n",
			MAX_NC_NAME,arrayname,VGNAMELENMAX,groupname);
	     MAPIERR(buff,funcname);
         }

         else 
             status_code = MAPIOK;
      }

      else 
          status_code = MAPIOK;
   
      if (status_code == MAPIOK) {

          if ( addid( file, arrayname, groupname, sds_id, DFTAG_NDG, NULL) == MFAIL){
              sprintf(buff,"ERROR: createMODISarray detected an error\n"
			"\t in M-API internal routine addid\n"
			"\t attempting to keep the array id for\n"
			"\t the %.*s array.",
			MAX_NC_NAME,arrayname);
	      MAPIERR(buff,funcname);
              status_code = MFAIL;    
          }
      }

      if (status_code == MFAIL){

          if ( SDendaccess(sds_id) == MFAIL){    
                sprintf(buff,"ERROR: createMODISarray detected an error\n"
                        "\t in HDF procedure attempting to close\n"
                        "\t the opened array which has errors.\n");
                MAPIERR(buff,funcname);
          }
      } 
  }
  return(status_code);
}
