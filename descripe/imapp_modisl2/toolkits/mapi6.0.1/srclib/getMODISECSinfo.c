#include <stdio.h> 
#include <stdlib.h>
#include "mapic.h"

int getMODISECSinfo(MODFILE *file, char const *PVLAttrName, 
                    char const *parmName, char *data_type, 
                    long int *n_elements, void *value)
/*
!C**********************************************************************
*!Description:
*   To retrieve  the value of a parameter from a HDF PVL text block.
*   Function getMODISECSinfo is part of a larger software system 
* called the MODIS Applications Programming Interface (API) 
* Utility, abbreviated M-API. The M-API Utility consists of 
* subroutines which allow MODIS Science Team-supplied software 
* to read  and write data and metadata from/to HDF files. The 
* functionality of the M-API is defined in the MODIS Application 
* Program Interface (API) Specification.
* 	
* In HDF-EOS, parameters are collected together to form a text block 
* using PVL. Then the text block is stored in HDF as a single 
* attribute. getMODISECSinfo (GMEATTR) retrieve the value of a 
* parameter from the PVL text block.
* 
* In order to obtain value of a parameter inside a PVL text block, 
* the function reads the PVL text block specified by PVLAttrName 
* from the MODIS file, creates the internal ODL tree structure from 
* the PVL text block, and search the tree structure to retrieve the 
* value of a parameter. The tree structure is then saved internally 
* for consecutive searches in the same PVL text block for code 
* efficiency. If multiple parameters will be retrieved from the 
* same PVL block, just set PVLAttrName to the HDF PVL attribute 
* name in the first call and set to NULL in C and ' ' in FORTRAN in 
* the consecutive calls. If the next call is to retrieve the value of a 
* parameter in a different PVL text block, set the PVLAttrName to 
* the new PVL attribute name. The saved old tree structure will be 
* deleted automatically and a new ODL tree will be created and 
* saved. If you will no longer call getMODISECSinfo in your 
* program and want to release the memory occupied by the saved 
* tree, just set both PVLAttrName and parmName to NULL in C or ' ' 
* in FORTRAN.
* 
*!Input parameters:
* file		IN: 	Address of MODFILE structure that is used to 
* 		reference the MODIS-HDF file containing the target 
* 		PVL attribute.
* PVLAttrName	IN:	ASCII string name of the HDF attribute 
* 		which contains the PVL text block. Set PVLAttrName 
* 		to NULL in C or ' ' in FORTRAN while parmName is 
* 		not NULL in C or not ' ' in FORTRAN will result in 
* 		searching the last PVL text block for the value of 
* 		parmName parameter.
* parmName	IN:	ASCII string name of a parameter whose 
*		value will be retrieved.  Set both PVLAttrName and 
*		parmName to NULL in C or ' ' in FORTRAN will 
*		release the memory occupied by the internal ODL 
*		tree. The parmName could parameter name only or 
*		combination of name and class represented as 
*		"name.class".
* 
* n_elements	IN/OUT: Address of the number of memory elements 
*		as data_type available in the value array. The 
*		attribute's value will not be retrieved unless 
*		*n_elements indicates that there is sufficient space 
*		available in value. getMODISECSinfo replaces this 
*		input with the number of elements required to 
*		contain the metadata. If the parameter cannot be 
*		found, *n_element will be left unchanged, or set to 0 
*		if a function error occurs. This argument must not 
*		be the address of a constant.
*		
*		SPECIAL CASE for multiple strings:
*		If there are multiple character strings for the 
*		parameters, strings will be packed together and 
*		returned in value . The separator between strings is 
*		'\0'. The low 16 bit of n_elements will return the 
*		total bytes in the values, including the '\0's between 
*		the strings and the '\0' at the end of last string. The 
*		part above the low 16 bits will return number of 
*		strings packed - 1. To obtain how many string 
*		retrieved, do the calculation:
*		n_strings  =  *n_elemets/65536  + 1
*		n_bytes  =  *n_elements%65536
*		Therefore, if *n_elements is less than 65536, there is 
*		only one strings in value and *n_elements is the 
*		number of bytes (characters) in the string, 
*		including the last '\0'.
* 
* data_type	IN/OUT: Data type of value. Output replaces with the 
*		data type of the retrieved  metadata. There are only 3 
*		data types in PVL:
*				"char *"
*				"int32"
*				"float64"
*		But you might use "float32" in input. If the value of 
*		the parameter is in "float64" type, this function will 
*		return in "float32" if users set the input value of 
*		dtype to "float32". This argument must be a 
*		variable.The memory size for dtype should be at 
*		least DATATYPELENMAX bytes long
* 
*!Output parameters:
* 
* value		OUT: 	buffer for the value. User should allocate 
*		enough memory for this buffer. If there are 
*		multiple data values in character type, the value 
*		will be placed consecutively. If the data value type is 
*		"char *", string will be separated by '\0'.
* 
* Return values:	
*           MAPIOK            successful
*           MFAIL             otherwise
*
* Externally defined:
* 		PGS_SMF_MAX_MSGBUF_SIZE		(mapic.h)
*		AGGREGATE			(odldef.h)
*		PARAMETER			(odldef.h)
*		PGSd_MET_CLASS_L		(PGS_MET.h)
*		RemoveAggregate			(odldef.h)
*		MAPIERR				(mapic.h)
*		MFAIL				(mapi.h)
*		NULLstr				(mapic.h)
*		MAX_NC_NAME			(netcdf.h)
*		MPVL2ODL			(mapic.h)
*		PGS_MET_NameAndClass		(PGS_MET.h)
*		FindObject			(odldef.h)
*		FindParameter			(odldef.h)
*		FirstValue			(odl+D149def.h)
*		datatype_to_DFNT		(mapic.h)
*		TV_STRING			(odldef.h)
*		TV_SYMBOL			(odldef.h)
*		TV_INTEGER			(odldef.h)
*		TV_REAL				(odldef.h)
*		DFNT_CHAR8			(hdf.h)
*		DFNT_CHAR			(hdf.h)
*		MAPIOK				(mapi.h)
*		DFNT_INT32			(hdf.h)
*		DFNT_FLOAT32			(hdf.h)
*		DFNT_FLOAT64			(hdf.h)
*		NextValue			(odldef.h)
*		DFNT_to_datatype		(mapic.h)
*		FAIL				(hdf.h)
*
*!Revision History:
* $Log: getMODISECSinfo.c,v $
* Revision 6.1  2010/07/13 19:30:38  kuyper
* Clarified code by adding extra curly brackets.
*
* Revision 5.1  2005/04/04 18:47:54  vlin
* constant safe for pointer arguments.
*
* 		Qi Huang	1996/05/02
*		Version 2.0
*		Original development and testing
*
*!Team-unique header:
* This software is developed by the MODIS Science Data SupportTeam for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
!END********************************************************************
*/
{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning message */
  char *funcname="getMODISECSinfo";    /* name of this routine */
  AGGREGATE		objectNode;
  PARAMETER		parmNode;
  VALUE			valueNode;
  int32			type;
  int32			n_elements_in;
  int32			nbytes;
  static AGGREGATE	aggNode = NULL;
  char			class[PGSd_MET_CLASS_L] = "\0";
  char			name[PGSd_MET_NAME_L] = "\0";
  int			status;

  status = MAPIOK;

  /* free the saved tree */
  if ( (PVLAttrName == NULL) && (parmName == NULL) )
  {
    if (aggNode != NULL)
      RemoveAggregate(aggNode);
    
    aggNode = NULL;
    return(MAPIOK);
  }

  /* Input checks: */
  if ( n_elements == NULL )
  {
    sprintf(buff,"ERROR: getMODISECSinfo can not continue without the\n"
			"\t n_elements input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  n_elements_in = (int32)*n_elements;
  *n_elements = 0;

  if ( NULLstr(parmName) )
  {
    sprintf(buff,"ERROR: getMODISECSinfo unable to access an ECS\n"
			"\t metadata without the parameter name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLstr(PVLAttrName) && (aggNode == NULL) )
  {
    sprintf(buff,"ERROR: getMODISECSinfo unable to access the %s\n"
			"\t metadata without the name of the global\n"
			"\t attribute it is stored within.\n",
			parmName);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( data_type == NULL )
  {
    sprintf(buff,"ERROR: getMODISECSinfo unable to access the %s\n"
			"\t metadata from ECS global attribute %.*s\n"
			"\t without the data type input.\n",
			parmName,MAX_NC_NAME,PVLAttrName);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* Create new ODL tree */
  if ( !NULLstr(PVLAttrName) )
  {
    if(aggNode != NULL)RemoveAggregate(aggNode);
    aggNode = MPVL2ODL(file,PVLAttrName,"locAgg");
    if ( aggNode == NULL )
    {
      sprintf(buff,"ERROR: getMODISECSinfo detected fails in\n"
			"\t procedure MPVL2ODL while attempting to\n"
			"\t retrieve parameter %s from ECS\n"
			"\t global attribute %.*s.\n",
			parmName,MAX_NC_NAME,PVLAttrName);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }  
  }

  PGS_MET_NameAndClass((char *)parmName,name,class);
  objectNode = FindObject(aggNode,name,class);
  if ( objectNode == NULL )
  {
    /* the parameter might be in a flat PVL form */
    parmNode = FindParameter(aggNode,name);
  }

  else
    parmNode = FindParameter(objectNode,PGSd_MET_ATTR_VALUE_STR);   

  if ( parmNode == NULL )
  {
    *n_elements = n_elements_in;
    sprintf(buff,"ERROR: getMODISECSinfo can not find the %s\n"
			"\t metadata.\n",parmName);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* get the values */
  valueNode = FirstValue(parmNode);

  if ( valueNode == NULL )
  {
    sprintf(buff,"ERROR: getMODISECSinfo found the value for parameter\n"
			"\t %s is undefined.\n",parmName);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( (type = datatype_to_DFNT(data_type)) == MFAIL )
  {
    status = MFAIL;
    type = 0;
  }

  nbytes = 0;
  if ( value == NULL )
  {
    sprintf(buff,"ERROR: getMODISECSinfo unable to access the %s\n"
			"\t metadata without the output data buff.\n",
			parmName);
    MAPIERR(buff,funcname);
    status = MFAIL;
  }

  while ( valueNode != NULL )
  {
    if ( (valueNode->item.type == TV_STRING) 
	|| (valueNode->item.type == TV_SYMBOL) )
    {
      if ( (type != DFNT_CHAR8) && (type != DFNT_CHAR) )
      {
	status = MFAIL;
	type = DFNT_CHAR8;
      }
      
      nbytes = nbytes + (int32)strlen(valueNode->item.value.string) + 1;

      if ( nbytes > n_elements_in )
	status = MFAIL;

      if ( status == MAPIOK )
      {
	strcpy(value,valueNode->item.value.string);
	value = (char *)value + strlen(valueNode->item.value.string) + 1;
      }
    }

    else if ( valueNode->item.type == TV_INTEGER )
    {
      if ( type != DFNT_INT32 )
      {
	status = MFAIL;
	type = DFNT_INT32;
      }

      if ( (*n_elements + 1) > n_elements_in )
	status = MFAIL;

      if ( status == MAPIOK )
      {
        *(int32 *)value = (int32) valueNode->item.value.integer.number;
	value = (int32 *)value + 1;
      }
    }

    else if ( valueNode->item.type == TV_REAL )
    {
      if ( (type != DFNT_FLOAT32) && (type != DFNT_FLOAT64) )
      {
	status = MFAIL;
	type = DFNT_FLOAT64;
      }

      if ( (*n_elements + 1) > n_elements_in )
	status = MFAIL;

      if ( status == MAPIOK )
      {
	if ( type == DFNT_FLOAT64 )
	{
	  *(float64 *)value = (float64) valueNode->item.value.real.number;
	  value = (float64 *)value + 1;
	}

	else
	{
	  *(float32 *)value = (float32) valueNode->item.value.real.number;
	  value = (float32 *)value + 1;
	}
      }
    }

    else
    {
      sprintf(buff,"ERROR: getMODISECSinfo found unknown ODL value\n"
			"\t type %ld for parameter %s.\n",
			(long)valueNode->item.type,parmName);
      MAPIERR(buff,funcname);
      return(MFAIL);
    }

    *n_elements= *n_elements + 1;
    valueNode = NextValue(valueNode);

  } /* end of while */

  /* Set the M-API data type */
  DFNT_to_datatype(type,data_type);

  /* Special for strings */
  if ( (type == DFNT_CHAR8) || (type == DFNT_CHAR) )
    *n_elements = (*n_elements - 1)*65536 + nbytes; /* *65536 is 16 
							bit shift */

  return(status);
}
