#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"


 int putMODISfileinfo (MODFILE *file, char const *attribute, 
                       char const *data_type, long int n_elements, 
                       void const *value)
/*
!C****************************************************************************
*
*Routine:     int putMODISfileinfo (MODFILE *file, char *attribute, 
*				char *data_type, long int *n_elements,
*				void *value)
*					
*!Description: Subroutine putMODISfileinfo is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
* 	      Utility, abbreviated M-API. The M-API Utility consists of
*	      subroutines which allow MODIS Science Team-supplied software
*	      to read and write data from/to HDF files. The functionality of
*             M-API is defined in the MODIS API User's Guide.
*
*	      putMODISfileinfo (PMFIN) store an attribute = value(s) metadata
*             pair to the indicated MODIS-HDF file. If the attribute already
*             exists, the value(s) will be updated.
*
*	      NOTE :  An assumption is made that int32 arrays are equivalent 
*		      to long int arrays.
* 
* !Input Parameters:MODFILE *file	: Address of MODISfile structure
*					  that references the MODIS-HDF file
*                                         recieving the metadata.
*		    char *attribute	: Name to assign the attribute. Provided
*                                         macros for accepted MODIS file
*                                         metadata names are listed in Appendix
*                                         A, MODIS API Supplied Constants.
*                   char *data_type     : Data type of the value.
*               long int n_elements     : Number of metadata values to extract
*                                         from value.
*                        void value     : Address of data to store in the
*                                         attribute. If the attribute already
*                                         exists, the value will be updated.
*						
* !Output Parameters:NONE
*
*Returns:   0 if successful, -1 on error, MFAIL (-1)  MAPIOK (0)
*
*Externally defined:	NULLstr		(mapic.h)
*                       MAX_REC         (mapic.h)
*                       datatype_to_DFNT(mapic.h)
*			SDsetattr	(mfhdf.h)
*                       DFKNTsize       (hproto.h)
*			TXT	        (mapi.h)
*                       MFAIL           (mapi.h)
*                       MAPIOK          (mapi.h)
*			MODFILE		(mapi.h)
* 	 	        PGS_SMF_MAX_MSGBUF_SIZE(PGS_SMF.h)
* 	 	        MAX_NC_NAME(netcdf.h)
* 	 	        NULLMODFIL(mapic.h)
* 	 	        DATATYPELENMAX(mapic.h)
* 	 	        VOIDP(hdfi.h)
*			DFACC_READ(hdf.h
* 	 	        MAPIERR(mapic.h)
*
* !Revision History:
*
*$Log: putMODISfileinfo.c,v $
*Revision 5.1  2005/04/04 18:49:11  vlin
*constant safe for pointer arguments.
*
*Revision 1.1  1998/02/06 22:26:06  fshaw
*Initial revision
*
 * Revision 1.5  1996/03/25  15:00:22  qhuang
 * Added 'MAX_NC_NAME,attribute' in else if ((data_t = datatype_t-_DFNT ...
 * conditional block.
 *
 * Revision 1.4  1996/03/07  22:16:53  qhuang
 * version 2 design, separated input checking, new message, removed the
 * specifical case for TXT type which uses strlen to obtain the length.
 *
 * Revision 1.3  1995/11/09  20:39:50  qhuang
 * minor changes
 *
 * Revision 1.2  1995/11/01  16:17:03  qhuang
 * Added capability to pass status messages to log files.
 *
 * Revision 1.1  1995/11/01  16:16:02  qhuang
 * Initial revision
 *
*
*	1.    Mitchell Weiss/RDC		Jul 95
*	     
*	      Original Development / Testing
*
*       2.    Xiao-Yang Ding/RDC               08/09/95
*
*             Revision after walk through.
* !Team-unique Header:
*
*              This software is developed by the MODIS Science Data Support
*	       Team for the National Aeronautics and Space Administration,
*	       Goddard Space Flight Center, under contract NAS5-32373.
*
*              Portions developed at the National Center for Supercomputing
*              Applications at the Univ. of Illinois at Urbana-Champaign.
!END**********************************************************************
*/

   {
   char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning 
                                            message */
   char *funcname="putMODISfileinfo";    /* name of this module */
   int   status;                         
   int32 data_t = 0;                            /* integer (data_type)      */
   long int size = 0L;

   status=MFAIL;

   /* Input checks: */
   if ( NULLstr(attribute) )
   {
     sprintf(buff,"ERROR: putMODISfileinfo unable to write a file\n"
			"\t attribute without an attribute name input\n");
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( NULLstr(data_type) )
   {
     sprintf(buff,"ERROR: putMODISfileinfo unable to write the %.*s\n"
			"\t file attribute without data type information.\n",
			MAX_NC_NAME,attribute);
     return(MFAIL);
   }

   if ( value == NULL )
   {
     sprintf(buff,"ERROR: putMODISfileinfo unable to write the %.*s\n"
			"\t file attribute without the value buffer.\n",
			MAX_NC_NAME,attribute);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( NULLMODFIL(file) )
   {
     sprintf(buff,"ERROR: putMODISfileinfo unable to continue with an\n"
			"\t invalid MODIS file structure input.\n");
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( n_elements < 1 )
   {
     sprintf(buff,"ERROR: putMODISfileinfo unable to store %ld\n"
			"\t %.*s file attribute values.\n",
			n_elements,MAX_NC_NAME,attribute);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if (file->access == DFACC_READ)
   {
      sprintf(buff,"ERROR: putMODISfileinfo unable to write the %.*s\n"
			"\t array attribute in a file opened for read only.\n",
			MAX_NC_NAME,attribute);
      MAPIERR(buff,funcname);
   }                   
   else if ( ( data_t = datatype_to_DFNT(data_type) ) == MFAIL )
   {
     sprintf(buff,"ERROR: putMODISfileinfo unable to write the %.*s\n"
			"\t file attribute of data type %.*s.\n",
			MAX_NC_NAME,attribute,
			DATATYPELENMAX,data_type);
     MAPIERR(buff,funcname);
   }

   else if ( ( size =  n_elements * DFKNTsize(data_t) ) > MAX_REC )
   {
     sprintf(buff,"ERROR: putMODISfileinfo unable to write the\n"
			"\t %.*s file attribute with a %ld byte\n"
			"\t value which exceeds the MAX_REC,32k.\n",
			MAX_NC_NAME,attribute,size);
     MAPIERR(buff,funcname);
   }
			
   else if ( SDsetattr((int)file->sd_id,attribute,data_t,(int32)n_elements,
			(VOIDP)value) == FAIL )
   {
     sprintf(buff,"ERROR: putMODISfileinfo detected FAIL from HDF\n"
			"\t procedure SDsetattr attempting to write the\n"
			"\t %.*s file attribute.\n",
			MAX_NC_NAME,attribute);
     MAPIERR(buff,funcname);
   }

   else
     status = MAPIOK;

   return status;
   }                                           /* end putMODISfileinfo      */
