#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"

void nccrmtbl(intf *modfil, _fcd tbname, intf *ltb, _fcd class, intf
              *lcl, _fcd group, intf *lgr, _fcd field, intf *lfi,
	      _fcd dtype, intf *ldt, intf *ret) 
/*
!C**********************************************************************
* 
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
*               createMODIStable. This C function is only called by FORTRAN function 
* 		CRMTBL. This function is a M-API internal routine.
* 
*!Description:  Function ccrmtbl is part of a larger software system called the 
*		MODIS Applications Programming Interface (API) Utility, 
*		abbreviated M-API. The M-API Utility consists of subroutines 
*		which allow MODIS Science Team-supplied software to read  and 
*		write data and metadata from/to HDF files. The functionality of 
*		the M-API is defined in the MODIS Application Program Interface 
*		(API) Specification.
* 
*		ccrmtbl is a C function which is callable from FORTRAN. This 
*		function will call createMODIStable to form a table structure. In 
*		M-API, ccrmtbl is low-level routine which is called only by 
*		CRMTBL. 
*		
*		In order to be callable from the FORTRAN in different platforms 
*		using function name ccrmtbl, this function is called nccrmtbl in 
*		the actual C code. nccrmtbl is redefined in mapic.h according to 
*		compiler's FORTRAN naming conventions/conversion of each 
*		platform, so that the object name of nccrmtbl will always the 
*		object name of a FORTRAN function named ccrmtbl.
* 
nput parameters:
* !Input Parameters:
* 		modfil	IN: 	Array  of FORTRAN integer array that is used 
* 	                 	to reference the MODIS-HDF file receiving the
*			        new table.
*	        tbname	IN:	FORTRAN character string that will be the 
*	      		        name of the table.
* 	        ltb	IN:	FORTRAN integer address of the memory size 
*			        of tbname.
* 		class	IN:	FORTRAN character string that will be the 
* 				class name of the table.
* 		lcl	IN:	FORTRAN integer address of the memory size 
* 				of class.
* 		group	IN:	ASCII string name of the data group to place 
* 				the new table in.
*		lgr	IN:	FORTRAN integer address of the memory size 
*				of group.
* 		field	IN:	FORTRAN comma-delimited character string 
* 				table headers.
*		lfi	IN:	FORTRAN integer address of the memory size 
*			        of field.
*		dtype	IN:	FORTRAN character string of comma-
* 				delimited data types. for each table field.
*	        ldt	IN:	FORTRAN integer address of the memory size 
*			        of dtype.
* 
* !Output Parameters:
*	        ret	OUT:	FORTRAN integer address of the *ret(MFAIL, 
* 				MAPIOK) 
* 
* Returns:	none
* 
* External references:
*			   MODFILE			 (mapi.h)
*			   _fcd				 (hdfi.h)
*			   intf				 (hdfi.h)
*                          MFAIL              		 (mapi.h)
*                          MAPIOK             		 (mapi.h)
*			   PGS_SMF_MAX_MSGBUF_SIZE	 (PGS_SMF.h)
*			   HDf2cstring			 (hproto.h)
*			   DATATYPELENMAX		 (mapic.h)
*			   MTYPEf2c			 (mapic.h)
*			   MAPIERR			 (mapic.h)
*			   createMODIStable		 (mapi.h)
*			   HDfreespace			 (hproto.h)
*			   VOIDP			 (hdfi.h)	
*
* !Revision History:
* $Log: ccrmtbl.c,v $
* Revision 5.1  2005/04/04 18:03:52  vlin
* constant safe for pointer arguments.
*
*		Qi Huang	1996/10/30
*		Version 2.2
*		Ring super structure and other changes make this version
*		much faster.
*
* !Team-unique Header:
*       This software is developed by the MODIS Science Data Support
*       Team for the National Aeronautics and Space Administration,
*       Goddard Space Flight Center, under contract NAS5-32373.
* 
*      HDF portions developed at the National Center for Supercomputing
*      Applications at the University of Illinois at Urbana-Champaign.
*
!END*********************************************************************
*/

{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char *funcname="nccrmtbl";       /* name of this module */
  MODFILE *file;      		/* pointer to structure for C HDF file */
  char const *tablename;   	/* pointer to C string for table name */
  char const *classname;   	/* pointer to C string for class name */
  char const *groupname;   	/* pointer to C string for group name */
  char const *fieldname;   	/* pointer to C string for field name */
  char *f_data_type = "";  /* pointer to C string for FORTRAN data types */
  char *c_data_type = "";     /* pointer to c_data_type allocated */
  long int c_length;      	/* length of c_data_type */
  int fieldno; 			/* number of fields */
  int i;             		/* local variable for controling FOR loop */

  *ret = MAPIOK;

  /* Convert the FORTRAN character strings to C character strings by
     using HDf2cstring. */
  groupname = HDf2cstring(group, (intn)*lgr);
  tablename = HDf2cstring(tbname, (intn)*ltb);
  classname = HDf2cstring(class, (intn)*lcl);
  fieldname = HDf2cstring(field, (intn)*lfi);
  f_data_type = HDf2cstring(dtype, (intn)*ldt);

  /* Set file by memcpy */
  memcpy(&file,&modfil[3],sizeof(MODFILE *));

  if (f_data_type != NULL)
  {
      /* Count number of fields by counting number of commas in the dtype.*/
      fieldno = 1;
      for (i=0; i< (int)strlen(f_data_type); i++) 
        if (f_data_type[i] == ',') 
          fieldno ++; 

      /* Allocate a memory with the size of( fieldno *
         DATATYPELENMAX + 1) for storing C data type. */

      if ( (c_data_type = (char *)malloc((size_t)(fieldno * DATATYPELENMAX + 1))) )
      {
	c_length = fieldno * DATATYPELENMAX + 1;
        *ret = MTYPEf2c(f_data_type,c_data_type,&c_length);

        if (*ret == MFAIL)
        { 
          sprintf(buff,"ERROR: nccrmtbl failed at data_type conversion\n");
          MAPIERR(buff,funcname);
        }
      }
      else
      {
        *ret = MFAIL;
        sprintf(buff,"ERROR: nccrmtbl out of memory\n");
        MAPIERR(buff,funcname);
      }
  }/* end of if */

  if ( *ret != MFAIL )
    *ret = createMODIStable(file,tablename,classname,groupname,fieldname,
           c_data_type);

  /* Free all allocated memory. */
  free((char *)c_data_type);
  if (tablename)
    HDfreespace((VOIDP)tablename);
  if (classname)
    HDfreespace((VOIDP)classname);
  if (fieldname)
    HDfreespace((VOIDP)fieldname);
  if (f_data_type)
    HDfreespace((VOIDP)f_data_type);
  if (groupname)
    HDfreespace((VOIDP)groupname);

  return;
}
