#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"

int  getMODIStable(MODFILE *file, char const *tablename, char const *groupname,
                   char const *fieldname, long int start, long int recno,
                   long int *buffsize, unsigned char *data)

/*
!C************************************************************************
*!Purpose:  Retrieves data records from a MODIS-HDF file table structure.
*
*!Description: 
*
*    Function getMODIStable is part of a larger software system called 
*    the MODIS Application Program Interface (API) Utility, abbreviated 
*    M-API.  The M-API Utility consists of subroutines which allow MODIS
*    Science Team supplied software to read and write data and metadata 
*    from/to HDF files.  The functionality of the M-API is defined in 
*    the MODIS API Specification.
*
*    getMODIStable retrieves one or more fields of data from one or more
*    records of an HDF Vdata table structure contained in a MODIS-HDF 
*    file.  The data are placed in the data buffer in consecutive records
*    and in the order that the input fieldnames are listed.  The length 
*    of this buffer must be able to contain all the fields requested 
*    times the number of records requested.  If the buffsize input 
*    indicates that it is too small to contain the actual quantity of 
*    data requested, getMODIStable will fail, but it will return the 
*    actual buffer size required.  
*
*    The groupname string provides the facility to select a table 
*    structure placed in a particular HDF 'Vgroup' data group.  
*    Alternatively, the entire file will be searched for a table 
*    structure named tablename if groupname = NULL in C or grpnm = ' ' 
*    in FORTRAN.
*
* !Input Parameters:
* file	IN: 	Address of MODFILE structure that is used to 
* 		reference the MODIS-HDF file containing the target 
* 		table structure.
* tablename	IN:	ASCII string name of the target table 
* 		structure.
* groupname	IN:	ASCII string name of the data group 
*		containing the target table structure.  If set to NULL 
*		the entire file will be searched for the table 
*		structure named tablename.
* fieldname	IN:	Array of comma-delimited ASCII string table 
* 		headers.  The data from each field will appear in the 
* 		same order as the headers.
* start		IN: 	Zero-based record location to begin reading 
* 		the data from the table structure.
* recno		IN:	Number of records to retrieve from the table 
* 		structure.
* buffsize	IN/OUT: Address of the data buffer size on input, in 
*		bytes.  The buffer must be at least this size.  buffsize 
*		will normally return the number of bytes of data 
*		successfully retrieved.  If the buffer is too small, 
*		however, the routine returns MFAIL and buffsize 
*		will contain the size a buffer must be to contain the 
*		output data.  It is set to 0 if a functional error occurs 
*		and making this output size determination is 
*		unreliable.
* 
* !Output Parameters:
* data		OUT:	Address of the data buffer.
* 
* Returns:	MAPIOK if successful, MFAIL if an error occurs.
*
* External reference:
*          MODFILE                 (mapi.h)
*          PGS_SMF_MAX_MSGBUF_SIZE (mapic.h)
*	   DATAID		   (mapi.h)
*	   VDINFO		   (mapic.h)
*          MAPIOK                  (mapi.h)
*          MFAIL                   (mapi.h)
*	   MAPIERR		   (mapic.h)
*          NULLstr                 (mapic.h)
*	   NULLMODFIL		   (mapic.h)
*          FULL_INTERLACE          (hdf.h)
*	   getMODIStableid	   (mapic.h)
*          VSNAMELENMAX		   (hdf.h)
*	   VSQuerycount		   (hdf.h)
*          VSsetfields             (vproto.h)
*          VSsizeof                (vproto.h)
*          VSseek                  (vproto.h)
*          VSread                  (vproto.h)
*	   VSQueryref		   (vg.h)
*          emptyVdata              (mapic.h)
*	   MAPIWARN		   (mapic.h)
*
* !Revision History:
*		Qi Huang	1996/09/23
*		Version 2.2
*
*		Ring super structure and other changes make this version
*		much faster.
*
* $Log: getMODIStable.c,v $
* Revision 6.1  2010/07/13 19:28:30  kuyper
* Corrected format specifier and corresponding argument's type to match.
*
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.4  1996/05/23  17:19:03  qhuang
 * (Note: This information replaces 'minor changes' in revision 1.3.) Added 
 * comments to first two character declarations, deleted one of the two character
 * definitions of the function name, in 'save input buffsize' SMF error message,
 * indicated which parameter was empty, corrected redundant logic in if state
 * ment, used NULLstr(tablename) and delete (tablename == NULL).
 *
* Revision 1.3  1995/11/07  18:50:12  qhuang
* minor changes
*
* Revision 1.2  1995/10/31  18:02:19  qhuang
* Added capability to pass status messages to log files.
*
* Revision 1.1  1995/10/31  17:59:43  qhuang
* Initial revision
*
* !Team-unique Header:
*
*    This software is developed by the MODIS Science Data Support
*    Team for the National Aeronautics and Space Administration,
*    Goddard Space Flight Center, under contract NAS5-32373.
*
* !References and Credits:
*    Portions developed at the National Center for Supercomputing
*    Applications at the Univ. of Illinois at Urbana-Champaign.
*
*    Frank Chen 
*    General Science Corporation
*    frank@modis-xl.gsfc.nasa.gov  Sept. 3, 1995
*
* !Design Notes:
*
*!END**********************************************************************
*/
 {
   char buff[PGS_SMF_MAX_MSGBUF_SIZE]; /* buffer for error/warning message */
   char *funcname="getMODIStable";     /* name of this routine */
   int status_code;             /* return value for routine. 
                                   MAPIOK = successful, MFAIL = fail */
   DATAID	*did;
   VDINFO	*vdinfo;
   int32 	nrecords;             /* total number of record */
   int32 	rec_byte_size;         /* byte size of Vdata field(s) */
   long int 	buffsize_in;        /* input *buffsize value */
   int32	output_size;
 
   /* initialize return status to MFAIL */
   status_code = MFAIL;
 
   /* input checks */
   if (buffsize == NULL) 
      {
        sprintf(buff, "ERROR: getMODIStable unable to continue without\n" 
               		"\t buffer size information. \n");
        MAPIERR(buff,funcname);
        return(MFAIL);
      }

   /* Save the contents of *buffsize to a local variable buffsize_in.  Buffsize 
      will now be used for output only.  Initialize *buffsize to 0. */ 
   buffsize_in = *buffsize;
   *buffsize = 0L;
 
   if (NULLstr(tablename)) 
   {
     sprintf(buff, "ERROR: getMODIStable unable to read from a table\n"
               		"\t without a table name input. \n");
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if (NULLstr(fieldname))
   {
     sprintf(buff,"ERROR: getMODIStable unable to read from a table\n"
			"\t without a fieldname input.\n");
     MAPIERR(buff,funcname);
     return(MFAIL);
   }
 
   if (NULLMODFIL(file))
   {
     sprintf(buff, "ERROR: getMODIStable unable to read from the %.*s" 
       			"\t table with a NULL file MODFILE structure. \n",
			VSNAMELENMAX,tablename);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if (data == NULL)
   {
     sprintf(buff, "ERROR: getMODIStable unable to read from the %.*s\n"
       			"\t table without a data buffer. \n",
			VSNAMELENMAX, tablename);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( (did = getMODIStableid(file,tablename,groupname,"r")) == NULL )
   {
     sprintf(buff,"ERROR: getMODIStable detected FAIL from procedure\n"
			"\t getMODIStableid attempting to access the\n"
			"\t %.*s table.\n",VSNAMELENMAX, tablename);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   vdinfo = (VDINFO *)did->info;

   /* Check that the requested read is consistent with the content of
      the retrieved Vdata: */

   if ( VSQuerycount((int32)did->id,&nrecords) == FAIL )
   {
     sprintf(buff,"ERROR: getMODIStable detected FAIL from HDF procedure\n"
			"\t VSQuerycount attempting to read the %.*s\n"
			"\t table.\n",VSNAMELENMAX, tablename);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( (start < 0) || (recno < 1) || ((start + recno) > nrecords) )
   {
     sprintf(buff,"ERROR: getMODIStable unable to read data from the\n"
			"\t %.*s table from invalid table structure locations.\n",
			VSNAMELENMAX, tablename);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   /* define the read fields */
   if ( (vdinfo->read_fields == NULL) ||
	(strcmp(fieldname,vdinfo->read_fields) != 0) )
   {
     if ( VSsetfields((int32)did->id,fieldname) == FAIL )
     {
       sprintf(buff,"ERROR: getMODIStable detected FAIL from HDF\n"
			"\t procedure VSsetfields attempting to read\n"
			"\t the %.*s table.\n",VSNAMELENMAX, tablename);
       MAPIERR(buff,funcname);
       return(MFAIL);
     }

     if ( (rec_byte_size = VSsizeof((int32)did->id,(char *)fieldname)) == FAIL )
     {
       sprintf(buff,"ERROR: getMODIStable detected FAIL from HDF\n"
			"\t procedure VSsizeof attempting to read\n"
			"\t the %.*s table.\n",VSNAMELENMAX, tablename);

       MAPIERR(buff,funcname);
       return(MFAIL);
     }

     /* reset the vdinfo structure */
     if ( vdinfo->read_fields != NULL )
       free(vdinfo->read_fields);
     vdinfo->read_fields = malloc(strlen(fieldname) + 1);
     if ( vdinfo->read_fields == NULL )
     {
     /*vdinfo->read_fields = NULL;*/
       vdinfo->size = 0;
       sprintf(buff,"ERROR: getMODIStable is out of memory\n"
			"\t attempting to read the %.*s table.\n",
			VSNAMELENMAX, tablename);
       MAPIERR(buff,funcname);
       return(MFAIL);
     }

     /* Copy fieldname to vdinfo->read_fields */
     strcpy(vdinfo->read_fields,fieldname);
     vdinfo->size = rec_byte_size;
   } /* end of if */

   output_size = (int32)(vdinfo->size * recno);
   if ( output_size > buffsize_in ) 
   {
     sprintf(buff, "ERROR: getMODIStable unable to fit %ld bytes "
        "of %.*s table's data into a %ld output buffer. \n",
        (long)output_size,VSNAMELENMAX,tablename,(long)buffsize_in);
     MAPIERR(buff,funcname);
     /* Revise the contents of buffsize with output_size (This is
	the minimum array size capable of holding the output data). */
     *buffsize = output_size;
   }
   else /* requested read is consistent with the content of the 
	   retrieved Vdata */
   {
     if (VSseek((int32)did->id,(int32)start) == FAIL) 
     {
       sprintf(buff, "ERROR: getMODIStable detected FAIL from HDF\n"
     			"\t procedure VSseek attempting to read the %.*s table.\n",
         		VSNAMELENMAX, tablename);
       MAPIERR(buff,funcname); 
     }
     else if (VSread((int32)did->id,data,(int32)recno,FULL_INTERLACE)==FAIL) 
     {
       sprintf(buff, "ERROR: getMODIStable detected FAIL from HDF\n"
	     			"\t procedure VSread attempting to read the\n"
				"\t %.*s table.\n",VSNAMELENMAX, tablename);
       MAPIERR(buff,funcname); 
     }
     else 
     {
       status_code = MAPIOK; /* successful operation */
       *buffsize = vdinfo->size * recno;

       /* check to determine that read is not a dummy record */
       if ((nrecords == 1) && (emptyVdata(file,VSQueryref((int32)did->id)) == 1))
       {
         sprintf(buff, "WARNING: getMODIStable retrieved dummy\n"
                		"\t record from empty %.*s table.\n",
				VSNAMELENMAX,tablename);
         MAPIWARN(buff,funcname); 
       }
     }/* end of else */
   }/* end of else; else (requested read is ... */

   return(status_code);
 }
