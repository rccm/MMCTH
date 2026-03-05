#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "mapic.h"

#define RECORDS_TO_STORE 	1
#define MAXLEN_DATATYPES	300
#define MAXLEN_FIELDNAME	290

int createMODIStable (MODFILE *file, char const *tablename, 
                      char const *classname, char const *groupname, 
                      char const *fieldname, char const *data_type)
/*
!C*****************************************************************************
*
*!Description: Subroutine createMODIStable is part of a larger software
*	     system called the MODIS Applications Programming
*	     Interface (API) Utility, abbreviated M-API.  M-API
*	     consists of subroutines that allow MODIS Science
*	     Team-supplied software to read and write data from/to HDF
*	     files.  The functionality of M-API is defined in the
*	     MODIS API User's Guide.
*
*	     createMODIStable creates an HDF Vdata structure in a MODIS
*	     HDF file to store a new data table.  It must be called before
*	     the data may be written to the file using putMODIStable.
*
*!Input Parameters:
*
*	     MODFILE *file	:	MODISfile to receive new table
*	     char *tablename	:	Name of the new table
*	     char *classname	:	Class name of the new table (optional)
*            char *groupname	:	Data group to place the new table in.
*				        (If NULL, then don't put in any group)
*	     char *fieldname	:	Array of comma delimited ASCII string
*					table headers.  The headers should be 
*					in the same order that the data for 
*					each table row will be written in.
*            char *data_type	:	Array of comma delimited data types
*					for each table field.  The data type
*					strings should be in the same order 
*					that the data for each table row will
*					be in.
*
*!Output Parameters:	None
*
*Return values:
*             MAPIOK            successful
*	      MFAIL             otherwise
*
*Externally defined:
*			MODFILE				(mapi.h)
*			PGS_SMF_MAX_MSGBUF_SIZE		PGS_SMF.h
*			NULLstr              		(mapic.h)
*			NULLMODFIL			(mapic.h)
*			DFAAC_READ           		(hdf.h)
*			VSNAMELENMAX			(hdf.h)
*			VGNAMELENMAX			(hdf.h)
*			searchMODISgroup		(mapic.h)
*			NO_OBJECT			(mapic.h)
*			MVSfind              		(mapic.h)
*			MODISsizeof          		(mapi.h)
*			MAX_REC              		(mapic.h)
*			VSFIELDMAX           		(hlimits.h)
*			FULL_INTERLACE			(hdf.h)
*			VSQueryref			(vg.h)
*			parse_string         		(mapic.h)
*			VSattach             		(vproto.h)
*			VSdetach			(vproto.h)
*			VSsetname            		(vproto.h)
*			VSsetclass           		(vproto.h)
*			VSsetfields          		(vproto.h)
*			VSwrite              		(vproto.h)
*			VSsizeof             		(vproto.h)
*			addMODISgroup			(mapic.h)
*			DFTAG_VH			(hdf.h)
*			datatype_to_DFNT     		(mapic.h)
*			VSfdefine            		(vproto.h)
*                       MAPIOK               		(mapi.h)
*                       MFAIL                		(mapi.h)
*			MAPIERR				(mapic.h)
*			addid				(mapic.h)
*
* !Revision History:
* $Log: createMODIStable.c,v $
* Revision 5.1  2005/04/04 18:17:57  vlin
* constant safe for pointer arguments.
* calls to function parse_string() updated.
*
* Revision 4.1  2003/01/27 20:45:05  kuyper
* Corrected to zero-initialize the dummy record buffer.
*
* !Team-unique Header:
*
*         This software is developed by the MODIS Science Data Support 
*         Team for the National Aeronautics and Space Administration, 
*         Goddard Space Flight Center, under contract NAS5-32373.
*
*         HDF portions developed at the National Center for Supercomputing 
*         Applications at the University of Illinois at Urbana-Champaign.
*
*!END
***************************************************************************/

{
   char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
   char *funcname="createMODIStable"; 
   int32  vdata_ref=-1;                        /* vdata reference num.
						   Set for new Vdata        */
   int  vdata_id=0;                             /* vdata identifier         */
   int  length = 0, size=0;
   char *field_names[VSFIELDMAX-1];         	/* array of field pointers  */
   char *data_types[VSFIELDMAX]={0};        	/* array of data type ptrs. */
   unsigned char *dumb_buf=NULL;		/* dummy write buffer	    */
   register int i=0;				/* generic counter    	    */
   int status_code; 	                        /* success/failure          */
   char *separator = ",";
   char *field_tmp="";                        
   char *data_tmp="";                           /* with parse_string        */
   int data_num=0;                              /* number of data types in 
						   data_type string         */
   int field_num=0;                     /* number of field names in a Vdata */
   status_code = MAPIOK;		/* success */

   /* Input checks: */
   if ( NULLstr(tablename) )
   {
     sprintf(buff,"ERROR: createMODIStable unable to make a new table\n"
			"\t without a table name input.\n");
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

   if ( NULLMODFIL(file) )
   {
     sprintf(buff,"ERROR: createMODIStable unable to make a new %.*s\n"
			"\t table with a NULL file MODFILE structure.\n",
			VSNAMELENMAX,tablename);
     MAPIERR(buff,funcname);
     return(MFAIL);
   }

  if ( NULLstr(fieldname) )
  {
    sprintf(buff,"ERROR: createMODIStable unable to make a new %.*s\n"
			"\t table without field name input.\n",
			VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLstr(data_type) )
  {
    sprintf(buff,"ERROR: createMODIStable unable to make a new %.*s\n"
			"\t table without field data types input.\n",
			VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( file->access == DFACC_READ )
  {
    sprintf(buff,"ERROR: createMODIStable unable to make a new %.*s\n"
			"\t table in file opened for read only.\n",
			VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( !NULLstr(groupname) )
  {
    if ( ( vdata_ref = searchMODISgroup(file,groupname,NULL,
                    tablename,classname,DFTAG_VH) ) == MFAIL )
    {
      sprintf(buff,"ERROR: createMODIStable unable to find data\n"
			"\t group %.*s to place new %.*s\n"
			"\t table in.\n",VGNAMELENMAX,groupname,
			VSNAMELENMAX,tablename);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }

    else if ( vdata_ref != NO_OBJECT )
         {
	   sprintf(buff,"ERROR: createMODIStable found the %.*s\n"
				"\t table already exists in the %.*s data\n"
				"\t group.\n",VSNAMELENMAX,tablename,
				VGNAMELENMAX,groupname);
           MAPIERR(buff,funcname);
           return(MFAIL);
         }
  }

  else if ( ( vdata_ref = MVSfind( (int32)file->hdf_id,tablename) ) != FAIL )
       {
	 sprintf(buff,"ERROR: createMODIStable found the %.*s table\n"
			"\t already exists.\n",VSNAMELENMAX,tablename);
	 MAPIERR(buff,funcname);
	 return(MFAIL);
       }
 
  if( (size=(int)MODISsizeof(data_type)) >MAX_REC )  /* get and check data_type size */
  {
    sprintf(buff, "ERROR: createMODIStable unable to create %.*s \n"
			"\t table with %d byte records.\n", 
			VSNAMELENMAX,tablename, size);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }
 
  if( size==MFAIL )		/* Check for unidentified data type */
  {
    sprintf(buff, "ERROR: createMODIStable unable to create %.*s \n"
			"\t table with %.*s data types.\n", 
			VSNAMELENMAX,tablename,MAXLEN_DATATYPES,data_type);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }
 
  /* Parse up the fieldname and datatype and error check: */
  if ( (field_tmp=(char *)malloc((size_t)(strlen(fieldname) + 1))) == NULL )
  {
    sprintf(buff, "ERROR: createMODIStable unable to allocate memory for\n"
			"\t the field name temporary buffer used to create\n"
			"\t the %.*s table.\n",VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* Copy fieldname into the allocated string to parse up field elements. */
  strcpy(field_tmp,fieldname);

   length = VSFIELDMAX;
   field_num = parse_string(field_tmp,separator,length,field_names);
   if (field_num <= 0 ) {
      sprintf(buff, "ERROR: createMODIStable found %.*s table to have\n"
			"\t no fields in the fieldname string %.*s.\n", 
	                VSNAMELENMAX,tablename,MAXLEN_FIELDNAME,fieldname);
    MAPIERR(buff,funcname); 
    status_code=MFAIL;
   }
   else if( field_num > VSFIELDMAX ) {
     sprintf(buff, "ERROR: createMODIStable unable to support the creation\n"
	              "\t of %d fields in the field name string\n"
	       	      "\t %.*s for the %.*s table.\n", field_num,
                      MAXLEN_FIELDNAME,fieldname,VSNAMELENMAX,tablename);
     MAPIERR(buff,funcname);
     status_code=MFAIL;
   }
 
   if ( (data_tmp=(char *)malloc((size_t)(strlen(data_type) + 1))) == NULL )
   {
     sprintf(buff, "ERROR: createMODIStable unable to allocate memory for\n"
       		      "\t the data type temporary buffer used to create\n"
			"\t the %.*s table.\n",VSNAMELENMAX,tablename);
     MAPIERR(buff,funcname);
     free(field_tmp);
     return(MFAIL);
   } 

   /* Copy data_type into an allocated string to parse up field elements. */
   strcpy(data_tmp,data_type);

   length = VSFIELDMAX;
   data_num = parse_string(data_tmp,separator,length,data_types);
   if (data_num != field_num ) {
      sprintf(buff, "ERROR: createMODIStable found %.*s table to\n"
		"\t have %d data types in the data type\n"
		"\t string %.*s instead of %d.\n", VSNAMELENMAX,tablename,
                data_num,MAXLEN_DATATYPES,data_type,field_num); 
     MAPIERR(buff,funcname);
     status_code=MFAIL;
   }

   if( status_code == MFAIL )
   {
     free(data_tmp);
     free(field_tmp);
     return(MFAIL);
   }

   /* Attach to new Vdata (vdata_ref=-1) */
   vdata_ref = -1;
   if ( (vdata_id=VSattach((int32)file->hdf_id, (int32)vdata_ref, "w")) == FAIL )
   {
     sprintf(buff, "ERROR: createMODIStable detected FAIL from HDF\n"
			  "\t procedure VSattach while attempting to create\n"
			  "\t the %.*s table.\n",VSNAMELENMAX,tablename);
     MAPIERR(buff,funcname);

     free(data_tmp);
     free(field_tmp);
     return(MFAIL);
   }

   VSsetname (vdata_id, tablename);   /* Set the current table name */

   if(!NULLstr(classname))
     VSsetclass (vdata_id, classname);  /* Set class name  */

   /* define each field as 1 order */	    
   for (i=0; i<field_num; i++)
     if( VSfdefine(vdata_id,field_names[i],datatype_to_DFNT(data_types[i]),1) == FAIL )
   {
      sprintf(buff, "ERROR: createMODIStable detected FAIL from HDF\n"
		 	  "\t procedure VSfdefine for %.*s and %.*s\n"
			  "\t of the %.*s table.\n",FIELDNAMELENMAX,field_names[i],
			  DATATYPELENMAX,data_types[i],VSNAMELENMAX,tablename);
      MAPIERR(buff,funcname);
      status_code=MFAIL;
      break;
   }

   if( status_code == MAPIOK )
   {
     if( (VSsetfields(vdata_id,fieldname)) == FAIL )
     {
       sprintf(buff, "ERROR: createMODIStable detected FAIL from\n"
	                     "\t HDF procedure VSsetfields creating\n"
			     "the %.*s table.\n",VSNAMELENMAX,tablename);
       MAPIERR(buff,funcname);
       status_code = MFAIL;
     }
   }
 
   /* create dummy buffer */
   if ( status_code == MAPIOK )
   {
     if ((dumb_buf = (unsigned char *)
	 calloc((size_t)VSsizeof(vdata_id, (char *)fieldname), 1)) == NULL )
     {
       sprintf(buff, "ERROR: createMODIStable unable to allocate\n"
			    "\t memory for dummy field buffer used to\n"
			    "\t create the %.*s table.\n",VSNAMELENMAX,tablename);
       MAPIERR(buff,funcname);
       status_code=MFAIL;
     }
   }

   
   if ( status_code == MAPIOK ) 
   {
     if ( (VSwrite (vdata_id, dumb_buf, RECORDS_TO_STORE, FULL_INTERLACE)) == FAIL )
     {
       sprintf(buff, "ERROR: createMODIStable detected FAIL from HDF\n"
		          "\t procedure VSwrite creating the %.*s\n" 
		          "\t table.\n", VSNAMELENMAX,tablename);
       MAPIERR(buff,funcname);
       status_code=MFAIL;       
     }
     free(dumb_buf);  /* frees dumb_buf even if it was written incorrectly */
   }
   
   if( ( status_code == MAPIOK ) && !NULLstr(groupname) )
   {
     if ( ( vdata_ref = VSQueryref(vdata_id) ) == FAIL )
     {
       sprintf(buff,"ERROR: createMODIStable detected FAIL from HDF\n"
			"\t procedure VSQueryref while createing the\n"
			"\t %.*s table.\n",VSNAMELENMAX,tablename);
       MAPIERR(buff,funcname);
       status_code = MFAIL;
     }
   }

   free(data_tmp);
   free(field_tmp);

   if ( (status_code == MAPIOK) && !NULLstr(groupname) )
     if ( addMODISgroup(file, groupname,NULL,DFTAG_VH, vdata_ref) == MFAIL )
     {
       sprintf(buff,"ERROR: createMODIStable unable to insert the\n"
			"\t %.*s table into the %.*s data\n"
			"\t group.\n",VSNAMELENMAX,tablename,VGNAMELENMAX,groupname);
       MAPIERR(buff,funcname);
       status_code = MFAIL;
     }

   if ( status_code == MAPIOK )
     if ( addid(file,tablename,groupname,vdata_id,DFTAG_VH,"w") == MFAIL )
     {
       sprintf(buff,"ERROR: createMODIStable detected an error in M-API\n"
			"\t internal routine addid attempting to keep the\n"
			"\t table id for the %.*s table.\n",
			VSNAMELENMAX,tablename);
       MAPIERR(buff,funcname);
       status_code = MFAIL;
     }
  
   if ( status_code == MFAIL )
     if ( VSdetach(vdata_id) == FAIL )
     {
       sprintf(buff,"ERROR: createMODIStable detected an error in HDF\n"
			"\t procedure VSdetach attempting to close the opend\n"
			"\t table which has errors.\n");
       MAPIERR(buff,funcname);
     }

   return(status_code);
}
