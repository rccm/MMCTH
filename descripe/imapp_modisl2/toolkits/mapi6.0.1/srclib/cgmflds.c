#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mapic.h"
#include "hproto.h"

void ncgmflds(intf *modfil, _fcd tbname, intf *ltb, _fcd group, intf *lgr,
              intf *strln, intf *recno, intf *fldno, _fcd fldnm, intf *lfl,
              _fcd dtype, intf *ldt, _fcd clss, intf *lcl, intf *ret)

/*
!C*****************************************************************************
* 
*!Description:	
*      A wrapping function interfacing between C and FORTRAN for 
* getMODISfields. This C function is only called by FORTRAN function 
* GMFLDS. This function is a M-API internal routine.
*      Function cgmflds is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read  and 
* write data and metadata from/to HDF files. The functionality of the M-API
* is defined in the MODIS Application Program Interface (API) Specification.
* 
* cgmflds is a C function callable from FORTRAN. This function will 
* call getMODISfields to get a table structure. In M-API, cgmflds is 
* low-level routine called only by FORTRAN function GMFLDS. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cgmflds, this function is called ncgmflds in 
* the actual C code. ncgmflds is redefined in mapic.h according to 
* FORTRAN compiler's naming conventions/conversion of each 
* platform, so that the object name of ncgmflds will always the 
* object name of a FORTRAN function named cgmflds.
* 
*!Input Parameters:
* 
* modfil	IN: 	Array  of FORTRAN integer array that is used 
* 		to reference the MODIS-HDF file.
* tbname	IN:	FORTRAN character string of the name of the vdata (table).
* ltb	IN:	FORTRAN integer address of the memory size 
* 		of tbname.
* group	IN:	ASCII string name of the data group to which 
* 		the table belongs.
* lgr	IN:	FORTRAN integer address of the memory size of group.
* lfl	IN:	FORTRAN integer address of the memory size of fldnm.
* ldt	IN:	FORTRAN integer address of the memory size of dtype.
* lcl	IN:	FORTRAN integer address of the memory size of class.
* 
*!Output Parameters:
* 
* strln	OUT:	The minimum length required to store either 
* 		fldnm or dtype. 0 if a function error occures.
* recno	OUT:	FORTRAN integer address of number of 
* 		records in the table.
* fldno	OUT:	FORTRAN integer address of number of fields in the table.
* fldnm	OUT:	FORTRAN comma-delimited character string table headers.
* dtype	OUT:	FORTRAN character string of comma-
* 		delimited data types for each table field.
* clss	OUT:	FORTRAN character string of table's class name.
* ret		OUT:	FORTRAN interger address of the status(MFAIL, MAPIOK) 
*
* Return values:	none
* 
* Externally defined:
*
*             PGS_SMF_MAX_MSGBUF_SIZE    (mapic.h)
*	      MODFILE			 (mapi.h)
*             getMODISfields             (mapi.h)
*             MAPIERR                    (mapic.h)
*             MTYPEC2f                   (mapic.h)
*             MFAIL                      (mapi.h)
*	      MAPIOK			 (mapi.h)
*             DATATYPELENMAX             (mapic.h)
*             HDf2cstring                (hproto.h)
*             HDc2fstr                   (hproto.h)
*             HDfreespace                (hproto.h)
*             _fcdtocp                   (hdfi.h)
*
*!Revision History:
* $Log: cgmflds.c,v $
* Revision 6.1  2010/07/13 18:58:05  kuyper
* Clarified code by adding curly brackets.
*
* Revision 5.1  2005/04/04 18:03:52  vlin
* Header file "hproto.h" included.
*
*
*             Xia W. Chang, xia@ltpmail.gsfc.nasa.gov  , 12/12/95
*             Version 1.4 original development from prototype by Liping Di.
*
*!Team-unique Header:
*
*             This software is developed by the MODIS Science Data Support
*             Team for the National Aeronautics and Space Administration,
*             Goddard Space Flight Center, under contract NAS5-32373.
*
*             Portions developed at the National Center for Supercomputing
*             Applications at the Univ. of Illinois at Urbana-Champaign.
*
!END*******************************************************************
*/
{
  char buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char *funcname = "cgmflds";

  /* Declare local character pointers */
  char *tablename, *data_type, *classname, *fieldname, *groupname;
  long int stringlen, c_recno, fieldno;
  long int f_length;
  MODFILE *file;

  stringlen = 0l;
  *strln = 0;

  /* Convert FORTRAN input strings to C string  */
  tablename = HDf2cstring(tbname, (intn)*ltb);
  groupname = HDf2cstring(group, (intn)*lgr);

  /* Set file by memcpy */
  memcpy(&file,&modfil[3],sizeof(MODFILE *));

  /* get the number of fields so the data_type array can be allocated. */
  if ((*ret = getMODISfields(file, tablename, groupname, NULL,
		NULL, &fieldno, NULL, NULL, NULL)) == MFAIL)
  {
    if(groupname) HDfreespace((VOIDP)groupname);
    if(tablename) HDfreespace((VOIDP)tablename);
    return;
  }

  /* allocate the memory for data_type  */
  stringlen = fieldno * DATATYPELENMAX;
  if ((data_type = (char*) malloc((size_t)stringlen)) == NULL ) 
  {
    sprintf(buff, "ERROR: ncgmflds out of memory \n");
    MAPIERR(buff, funcname);
    *ret = MFAIL;
    HDfreespace((VOIDP)groupname);
    HDfreespace((VOIDP)tablename); 
    return;
  } 

  /*  directly convert the FORTRAN string memory */
  fieldname = _fcdtocp(fldnm);
  classname = _fcdtocp(clss);

  stringlen = *lfl;

  /* get the field information  */	
  *ret = getMODISfields(file, tablename, groupname, &stringlen, 
		   &c_recno,&fieldno, fieldname, data_type, classname);

  /*  free groupname and tablename  */
  if (groupname) HDfreespace((VOIDP)groupname);
  if (tablename) HDfreespace((VOIDP)tablename);

  *recno = (intf)c_recno;
  *fldno = (intf)fieldno;

  /* convert C data type to FORTRAN data type */
  *strln = 0;
  if ( (*ret != MFAIL) || ( (*ret == MFAIL) && (stringlen > 0) ) )
  {
    *strln = (intf)stringlen;
    if (*ret == MAPIOK) 
	*strln = (intf)(strlen(fieldname) + 1);
    HDc2fstr(classname, (intn)(*lcl));
    HDc2fstr(fieldname, (intn)(*lfl));
    f_length= *ldt;

    if ( MTYPEc2f(data_type, _fcdtocp(dtype), &f_length) == MAPIOK )
    {
      /* add trailing blanks to dtype using HDc2fstr  */
      HDc2fstr(_fcdtocp(dtype), (intn)(*ldt));
      if (*ret != MFAIL) 
	*strln = *strln -1 ;
	if (*strln < f_length) 
	*strln = (intf)f_length;
    } 
    else 
    {
      *ret = MFAIL;
      if ( f_length > *ldt )
      {
        if ( f_length > *strln ) 
          *strln = (intf)f_length; 
      }
      else 
	*strln = 0;
    }
  }

  /* free data_type   */
  if ( data_type )
    free(data_type);

  return;
}
