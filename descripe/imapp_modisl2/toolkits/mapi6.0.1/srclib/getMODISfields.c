#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mapic.h"

/* define macro MAX_V to return the maximum of two values */
#define MAX_V(a,b) (((a) > (b))?(a):(b))

int getMODISfields (MODFILE *file, char const *tablename, char const *groupname,  
                    long int *stringlen, long int *recno, long int *fieldno, 
                    char *fieldname, char *data_type, char *classname)
/*
!C*****************************************************************************
*
*!Description: Function getMODISfields is part of a larger software system
*             called the MODIS Applications Programming Interface (API)
*	      Utility, abbreviated M-API.  The M-API Utility consists of
*	      subroutines which allow MODIS Science Team-supplied software
*	      to read in Level 1B radiance bands and write out output 
*	      products and metadata to HDF files.  The functionality of the
*	      M-API is defined in the MODIS API User's Guide, Version 1.1, 
*	      dated 6/16/95.
*
*	      getMODISfields retrieves the essential characteristics of an 
*             HDF Vdata table structure contained in a MODIS-HDF file.  This
*             provides the information neede for properly reading data
*             from the table structure using getMODIStable or to write to it 
*             using putMODIStable.  If any of the outpu parameters are set
*             to NULL, then that data are not retrieved.  An error(MFAIL) will
*             be returned if 1) The output strings are not long enough to
*             contain the data type or field name strings for all the
*             Vdata's fields, 2) an unknown (e.g. not supported by the MODIS
*             API) number type is encountered or 3) an HDF routine FAILs.
*             The dataq type string (if requested) will be returned 
*             truncated to the point where the fault occurred.  stringlen,
*             the address of the length of the data_type and filename output
*             strings, is a required input if either of these strings is to
*             be retrieved.  It is normally revised to contain the actual
*             array length required to hold the larger of the two output
*             strings.  If the later two errors occur, however, *stringlen is
*             set to 0.
*
*             The groupname string provides the facility to select a table
*             structure existing in a particular HDF 'Vgroup' data group.
*             Alternatively, the entire file will be searched for a table
*             structure named tablename if groupname = NULL.
*
*             The essential characteristics returned by getMODISfields 
*             are the number of records, number of fields,  the field
*             names, data types, and the class of a specified MODIS-HDF
*             file table structure.
*
* !Input Parameters:
* file          IN:  Address of MODFILE structure that is used to reference
*               the MODIS-HDF file containing the target array      
*               structure.
* tablename     IN:  ASCII string name of the target table structure.
* groupname     IN:  ASCII string name of the dat group containing the target table
*               structure.  If set to NULL the entire file will be searched 
*               for the table structure name tablename. 
* stringlen     IN/OUT:  Address of minimum length of fieldname and data_type arrays.
*               Return the minimum array length actually required to hold the
*               longer of the two strings. It is set to 0 if a functional
*               error occurs.
*
* !Output Parameters:
* recno         OUT:  Number of records (rows) present in the table structure.
* fieldno       OUT:  Number of fields (columns) present in the table structure.
* fieldname     OUT:  Array of comma-delimited ASCII string table headers.
* data_type     OUT:  Array of comma-delimited data types for each table field.
*               Possible C data types:
*                      "int8"
*                      "uint8"
*                      "int16"
*                      "uint16"
*                      "int32"
*                      "uint32"
*                      "int64"
*                      "float32"
*                      "float64"     
* classname     OUT:  ASCII string for the class name of the table structure.
*               Provided array should be at least VSNAMELENMAX(64) bytes long.
* 
*Returns:       MAPIOK if successful
*               MFAIL if an error occurs
*
*Externally defined:
*                       emptyVdata(mapic.h)
*                       VFdatatypes(mapic.h)
*                       VSgetfields(vproto.h)
*                       VSgetclass(vproto.h)
*                       VSinquire(vproto.h)
*			VSQureryref(vg.h)
*
*                       NULLstr(mapic.h)
*                       MODFILE(mapi.h)
*                       MAPIOK(mapi.h)
*                       MFAIL(mapi.h)
*                       FAIL(hdf.h)
*                       VSFIELDMAX(hdf.h)
*                       FIELDNAMELENMAX(hdf.h)
*			PGS_SMF_MAX_MSGBUF_SIZE(PGS_SMF.h)
*			MAPIERR(mapic.h)
*			NULLMODFIL(mapic.h)
*			VSNAMELENMAX(hdf.h)
*			getMODIStableid(mapic.h)
*
* !Revision History:
*		Qi Huang	1996/10/08
*		Version 2.2
*		Ring super structure and other changes make this version
*		much faster.
*
*		Qi Huang 1996/02/13
*		Version 2.0 development and testing
*
*$Log: getMODISfields.c,v $
*Revision 5.1  2005/04/04 18:49:11  vlin
*constant safe for pointer arguments.
*
*Revision 1.1  1998/02/06 22:26:06  fshaw
*Initial revision
*
 * Revision 1.5.1.2  1996/03/13  20:56:04  qhuang
 * replace MFAIL in call MVSfind with FAIL.
 *
 * Revision 1.5.1.1  1996/02/29  18:06:58  qhuang
 * *** empty log message ***
 *
 * Revision 1.5  1996/02/13  18:31:11  qhuang
 * Version2.0.  Added group handling, new messages and input checking.
 *
 * Revision 1.4  1996/02/12  15:33:25  qhuang
 * *** empty log message ***
 *
 * Revision 1.4 1995/12/29  Xia W. Chang xia@modis-xl.gsfc.nasa.gov
 * Replace HDF VSfind routine to mapi MVSfind.
 * 
 * Revision 1.3  1995/11/09  20:02:19  qhuang
 * Expanded comments in the code to elaborate the information from the PDL.
 * Added comoments to first two char declarations.
 *
 * Revision 1.2  1995/11/02  21:21:48  qhuang
 * Added capability to pass status messages to log files.
 *
 * Revision 1.1  1995/11/02  21:21:05  qhuang
 * Initial revision
 *
*
*              Revision 01.00  1995/04/24   12:00 EDT
*	      David J. Lerda/RDC
*
*	      Create the initial template for MODIS API utility function 
*	      getMODISfields.
*
*              Frank Chen      1995/08/27
*
*              Update the code to comply with the PDL as of 1995/08/23.
*
*              Frank Chen      1995/09/07
*
*              Update the code according to the result of first walk-through
*              at 1995/09/06 to change the original coding style and eliminate
*              extra include files like mapi.h and hdf.h. Module test code is
*              moved to another file also.
*              
*              Liping Di, Hughes STX, Oct. 17, 1995
*              Using a local variable recno32 for the call to VSinqure for the portability.
*              Using NULL directly in VSinquire.
*              
*              Liping Di, Hughes STX, Oct 31, 1995
*              delete the local_class. Modifying the part of retrieving data_type code which 
*              does not agree with design. Modifying the part of return stringlen which does 
*              not agree with the design.
*              
* !Team-unique Header:
*
*              This software is developed by the MODIS Science Data Support 
*              Team for the National Aeronautics and Space Administration,
*	      Goddard Space Flight Center, under contract NAS5-32373.
*                
*              HDF portions developed at the National Center for Supercomputing
*              Applications at the University of Illinois at Urbana-Champaign.
!END*************************************************************************
*/


{
 char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
 char *funcname="getMODISfields";

 DATAID        *did;
 int status_code;        	/*  return value for function.  MAPIOK = successful */
 int32 number_of_fields = 0; 	/*  number of fields in the Vdata  */
 long int field_len = 0;     	/* string length of fieldname */
 long int dtype_len = 0;     	/* string length of data_type */

 /*char local_fields[VSFIELDMAX * (FIELDNAMELENMAX + 1)];*/
 char local_fields[32767];      /* use 32767 to replace VSFIELDMAX * (FIELDNAMELENMAX + 1) */

 int32 recno32; 

 status_code = MAPIOK;     /* Set status code to MAPIOK, success. */

 /* Input checks: */
 if ( NULLstr(tablename) )
 {
   sprintf(buff,"ERROR: getMODISfields unable to access a table without\n"
			"\t a table name input.\n");
   MAPIERR(buff,funcname);
     
   if ( stringlen != NULL )
     *stringlen = 0;	/* output length is indeterminant. */
   
   return(MFAIL);
 }

 if ( NULLMODFIL(file) )
 {
   sprintf(buff,"ERROR: getMODISfields unable to access the %.*s\n"
			"\t table with a NULL file MODFILE structure.\n",
			VSNAMELENMAX,tablename);
   MAPIERR(buff,funcname);
   
   if ( stringlen != NULL )
     *stringlen = 0;	/* output length is indeterminant. */

   return(MFAIL);
 }

 if ( (did = getMODIStableid(file,tablename,groupname,"r")) == NULL )
 {
   sprintf(buff,"ERROR: getMODISfields detected FAIL from procedure\n"
			"\t getMODIStableid attempting to access the\n"
			"\t %.*s table.\n",VSNAMELENMAX,tablename);
   MAPIERR(buff,funcname);
   return(MFAIL);
 }

 if ( (fieldno != NULL) || (fieldname != NULL) )
 {
   /* Call VSgetfields to retrieve the fieldnames and number of fields
   in the Vdata.  Both outputs are passed to locally declared, initialized
   variables: Retrieve the fieldname string to a VSFIELDMAX*(FIELDNAMELENMAX
   + 1) array.  (NOTE: This exceeds the ANSI minimum required array support
   size, so this may fail on some wimpy compilers). */

   if ( (number_of_fields = VSgetfields((int32)did->id,local_fields)) == FAIL )
   {
     sprintf(buff, "ERROR: getMODISfields detected FAIL from HDF\n"
                 	"\t procedure VSgetfields.\n");
     MAPIERR(buff,funcname);
     status_code = MFAIL;    /* Fail */
   }
   else /* VSgetfields successful */
   {
     if ( fieldno != NULL )
       *fieldno = number_of_fields;

     field_len = (long)strlen(local_fields);
     if (fieldname != NULL) 
     {
       if ( stringlen == NULL )
       {
	 sprintf(buff,"ERROR: getMODISfields unable to fit\n"
			"\t %.*s table's %ld\n"
			"\t byte field names\n"
			"\t string into output string of\n"
			"\t unknown length.\n",VSNAMELENMAX,tablename,
			(long)field_len);
	 MAPIERR(buff,funcname);
	 status_code = MFAIL;
       }

       else if ( field_len > *stringlen )
	    {
	      sprintf(buff,"ERROR: getMODISfields unable to fit\n"
				"\t %.*s table's %ld\n"
				"\t byte field names into\n"
				"\t %ld byte output string.\n",
				VSNAMELENMAX,tablename,
				(long)field_len,*stringlen);
	      MAPIERR(buff,funcname);
	    }

       else
	 strcpy(fieldname,local_fields);

     } /* end of if; if( fieldname != NULL ) */     

   } /* end of else; (VSgetfields successful) */   

 } /* end of if; if ((fieldno != NULL) || (fieldname != NULL)) */

 if ( data_type != NULL )
 {
   if ( stringlen == NULL )
   {
     sprintf(buff,"ERROR: getMODISfields unable to fit %.*s\n"
			"\t table's data types string into output\n"
			"\t string of unknown length.\n",
			VSNAMELENMAX,tablename);
     MAPIERR(buff,funcname);
     status_code = MFAIL;
   }

   else
   {
     dtype_len = *stringlen;
     if ( VFdatatypes((int32)did->id,&dtype_len, data_type) == MFAIL )
     {
       if ( dtype_len == 0 )  /* a functional error occurred. */
       {
	 sprintf(buff,"ERROR: getMODISfields detected FAIL\n"
			"\t retrieving the data type\n"
			"\t string for the %.*s table\n"
			"\t using VFdatatypes.\n",
			VSNAMELENMAX,tablename);
	 MAPIERR(buff,funcname);
	 status_code = MFAIL;
       }

       else  /* the datatype string is too short for VFdatatypes' output */
       {
         sprintf(buff,"ERROR: getMODISfields unable to fit\n"
			"\t the %.*s tables %ld\n"
			"\t byte data types into\n"
			"\t %ld	byte output string.\n",
			VSNAMELENMAX,tablename,dtype_len,*stringlen);
         MAPIERR(buff,funcname);
       }

     } /* end of if; if ( VFdatatypes(... */

   } /* end of else; */

 } /* end of if; if ( data_type != NULL ) */

 if ( classname != NULL )
 /* Retrieve the Vdata's class name into classname using VSgetclass. */
   VSgetclass((int32)did->id, classname);

 if ( recno != NULL ) 
 {
  /* Call VSinquire to retrieve the number of records in the Vdata. */
   if ( VSinquire((int32)did->id,&recno32, NULL, NULL, NULL, NULL) == FAIL )
   {
     sprintf(buff, "ERROR: getMODISfields detected FAIL from HDF\n"
                     "\t procedure VSinquire.\n");
     MAPIERR(buff,funcname);
     status_code = MFAIL;   /* Fail */
   }

   /* ELSE IF there is only one record in the Vdata AND emptyVdata
       indicates that it is only a dummy record. */
   else if ( ((*recno = (long)recno32) == 1) && (emptyVdata(file,VSQueryref((int32)did->id)) == 1) ) 
          *recno = 0;   /* no records in the Vdata */

 } /* end of if; recon != NULL */ 

 if (stringlen != NULL) 
 {
   if ( status_code == MFAIL )
   /* array length information may be unreliable */
     *stringlen = 0;
   else
   { 
     if( (*stringlen < field_len) || (*stringlen < dtype_len) )
       status_code = MFAIL;    /* Fail */

     *stringlen = MAX_V(field_len,dtype_len);
   }
 }

 return(status_code); 
}  
