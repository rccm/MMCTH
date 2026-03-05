#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mapic.h"

void ncgmecin(intf modfil[MODFILLEN], _fcd pvlname, intf *lpv,
		_fcd pname, intf *lpn, _fcd dtype, intf *ldt, intf *nms,
		VOIDP pvalue, intf *ret)
/*
!C**********************************************************************
* 
*!Purpose:	A wrapping function interfacing between C and FORTRAN for 
* getMODISECSinfo. This C function is only called by FORTRAN function 
* GMECIN. This function is a M-API internal routine.
* 
*!Description: Function cgmecin is part of a larger software system called the 
* MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines 
* which allow MODIS Science Team-supplied software to read  and 
* write data and metadata from/to HDF files. The functionality of 
* the M-API is defined in the MODIS Application Program Interface 
* (API) Specification.
* 
* cgmecin is a C function which is callable from FORTRAN. This 
* function will call getMODISECSinfo to form a table structure. In 
* M-API, cgmecin is low-level routine which is called only by 
* GMEATTR. 
* 
* In order to be callable from the FORTRAN in different platforms 
* using function name cgmecin, this function is called ncgmecin 
* in the actual C code. ncgmecin is redefined in mapic.h according 
* to compiler's FORTRAN naming conventions/conversion of each 
* platform, so that the object name of ncgmecin will always the 
* object name of a FORTRAN function named cgmecin.
* 
*!Input parameters:
* modfil	IN: 	MODFILE file array that is used to reference 
*		the MODIS-HDF file containing the target PVL 
*		attribute.
* pvlname	IN:	FORTRAN ASCII string name of the HDF 
*		 attribute which contains the PVL text block.
* lpv		IN:	Number of bytes in pvlname.
* pname		IN:	FORTRAN ASCII string name of a parameter 
*		whose value will be retrieved.
* lpn		IN:	Number of bytes in pname
* ldt		IN:	number of bytes in dtype.
* 
* nms		IN/OUT:	number of values in pvalue
* dtype		IN/OUT:	Type of data in pvalue. There are only 3 
*		data types:
*				'CHARACTER*(*) '
*				'INTEGER32'
*				'FLOAT64'
*		Users may also use 'FLOAT32', for which the function 
*		will convert the value from float64 to float32.
* 
*!Output parameters:
* pvalue	OUT: 	buffer for the value. User should allocate 
* 		enough memory for this buffer.
* ret		OUT: 	Return status of the call to C routine. 
*
* Returns:	none
* 
* External References:
*			   MODFILLEN			 (mapic.h)
*			   MODFILE			 (mapi.h)
*			   HDf2cstring			 (hproto.h)
*			   P_ADDR  			 (mapic.h)
*			   DATATYPELENMAX		 (mapic.h)
*			   MTYPEf2c			 (mapic.h)
*			   MTYPEc2f			 (mapic.h)
*			   HDfreespace			 (hproto.h)
*			   VOIDP			 (hdfi.h)	
*			   TXT				 (mapi.h)
*			   getMODISECSinfo		 (mapi.h)
*			   MAPIOK			 (mapi.h)
*
*!Revision history:
*		Qi Huang	1996/09/04
*		Version 2.1
*
*		Ring super structure and other changes make this
*		version much faster.
*
* 		Qi Huang	1996/05/07
*		Version 2.0
*		Original development and testing
*
* $Log: cgmecin.c,v $
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.3  1996/09/18  16:29:46  qhuang
 * Version 2.1
 *
 * Revision 1.2  1996/06/03  20:47:04  qhuang
 * Added NCSA acknowledgement in prolog, added HDc2fstring to External
 * references, added comments to all local variables.
 *
 * Revision 1.1  1996/05/07  17:59:14  qhuang
 * Initial revision
 *
*
*!Team-unique header:
* This software is developed by the MODIS Science Data SupportTeam for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
*!References and Credits
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
*!Design Notes
*
!END********************************************************************
*/
{
  MODFILE	*file;		/* MODFILE structure */
  char		*cpvlname;	/* pointer to converted C string for pvlname */
  char		*cpname;	/* pointer to converted C string for pname   */
  char		*fdtype;	/* pointer to converted C string for dtype   */
  char		*cvalue;	/* local pointer variable to point to the
				   same place where pvalue points to */
  char		cdtype[DATATYPELENMAX]; /* character array to hold retrieved
					C data type */
  long int	lnms;		/* local variable to hold *nms passed in */
  long int	c_length = DATATYPELENMAX; /* local variable to hold length
						of C data type retrieved */
  long int	f_length = *ldt;            /* local variable to hold length
					of FORTRAN data type retrieved */
  int		i;		/* local variable to control FOR loop */
  

  /* Set file by memcpy */
  memcpy(&file,&modfil[P_ADDR],sizeof(MODFILE *));

  /* Convert pvlname to cpvlname,pname to cpname, dtype to fdtype by
	 using HDf2cstring. */ 
  cpvlname = HDf2cstring(pvlname, (intn)*lpv);
  cpname = HDf2cstring(pname, (intn)*lpn);
  fdtype = HDf2cstring(dtype, (intn)*ldt);
  lnms = *nms;

  *cdtype = '\0';
  MTYPEf2c(fdtype,cdtype,&c_length); 

  if (fdtype) HDfreespace((VOIDP)fdtype);


  if (strlen(cpvlname) == 0)
  {
    HDfreespace((VOIDP)cpvlname);
    cpvlname = NULL;
  }

  if (strlen(cpname) == 0)
  {
    HDfreespace((VOIDP)cpname);
    cpname = NULL;
  }
  /* Call getMODISECSinfo. */
  *ret = getMODISECSinfo (file, cpvlname,cpname, cdtype, &lnms,pvalue);

  if ( (*ret == MAPIOK) || (lnms != 0) )
  {
    fdtype = _fcdtocp(dtype);
    MTYPEc2f(cdtype,fdtype,&f_length);
    HDc2fstr(fdtype,(intn)*ldt);
  }

  if ( (*ret == MAPIOK) && (strcmp(cdtype,TXT) == 0) )
  {
    cvalue = (char *)pvalue;
    for (i=(int)(lnms%65536 - 1); i< *nms; i++)
      cvalue[i] = ' ';
    *nms = (intf)lnms - 1;
  }
  else
    *nms = (intf)lnms;

  if (cpvlname)
    HDfreespace((VOIDP)cpvlname);

  if (cpname)
    HDfreespace((VOIDP)cpname);

  return;
}
