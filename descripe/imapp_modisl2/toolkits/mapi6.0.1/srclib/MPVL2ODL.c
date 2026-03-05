#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"
#include <PGS_IO.h>
#include <PGS_MET.h>

AGGREGATE MPVL2ODL(MODFILE *file, char const *PVLAttrName, char const *aggName)
/*
!C**********************************************************************
* 
*!Description: Function MPVL2ODL is part of a larger software system
* called the MODIS Applications Programming Interface (API) Utility, 
* abbreviated M-API. The M-API Utility consists of subroutines which
* allow MODIS Science Team-supplied software to read  and write data
* and metadata from/to HDF files. The functionality of the M-API is 
* defined in the MODIS Application Program Interface (API) 
* Specification.
* 	
* In HDF-EOS, attributes are collected togather to form a text block 
* using PVL. Then the text block are stored in HDF as a family of
* attributes. MPVL2ODL reads in a MODIS PVL (Parameter = Value) 
* attribute family in a HDF file and converts it to ODL stored in a tree 
* structure for searching individual attribute in PVL. This is M-API 
* internal routine called by getMODISECSinfo.
* 
* The function seaches the MODIS HDF global attributes first. If the 
* PVLAttrname attributes can not be found in the global attributes, 
* the function will then search the attributes attached to individual 
* SDSs. Once the PVL attribute(s) is(are) found, the function will retrieve 
* it(them) and convert to ODL.
* 
*!Input parameters:
* file		Address of MODFILE structure that is used to reference
*		the MODIS-HDF file containing the target PVL attribute.
* PVLAttrName	ASCII string name of the HDF attribute which contains
* 		the PVL text block.
* aggName	ASCII string name used to name the AGGREGATE structure.
* 
*!Output parameters:		
*		None.
* 
* Returns:	AGGREGATE if successful, NULL if an error occurs.
*		(Note AGGREGATE is a pointer to a tree structure defined in 
*		odldef.h)
*
* Externally Defined:
*               AGGREGATE                       (odldef.h)
*               errno                           (errno.h)
*               FAIL                            (hdf.h)
*               MAX_NC_NAME			(netcdf.h)
*               MAX_ORDER                       (hlimits.h)
*		PGS_SMF_MAX_MSGBUF_SIZE		(mapic.h)
*		PGS_S_SUCCESS			(PGS_SMF.h)
*               PGS_IO_Gen_Open                 (PGS_IO_Gen_Wrap.h)
*               PGSt_PC_Logical                 (PGS_PC.h)
*               PGSd_MET_TEMP_ATTR_FILE         (PGS_MET.h)
*               SUCCEED                         (hdf.h)
* 
* Called By:
*               getMODISECSinfo                 M-API Routine to get ECS 
*                                               attribute values.
*
* Routines Called:
*		
*		DFKNTsize			(hproto.h)
*		MAPIERR				(mapic.h)
*		NULLstr				(mapic.h)
*		NULLMODFIL			(mapic.h)
*		PGS_IO_Gen_Close		(PGS_IO_Gen_Wrap.h)
*               PGS_IO_Gen_Open                 (PGS_IO_Gen_Wrap.h)
*		PGS_MET_ConvertToOdl		(PGS_MET.h)
*               PGS_SMF_SetUNIXMsg              (PGS_SMF.h)
*		SDendaccess			(mfhdf.h)
*		SDfindattr			(mfhdf.h)
*		SDfileinfo			(mfhdf.h)
*		SDattrinfo			(mfhdf.h)
*		SDreadattr			(mfhdf.h)
*               
*
*!Revision History:
* $Log: MPVL2ODL.c,v $
* Revision 5.1  2005/04/04 16:48:04  vlin
* constant safe for pointer arguments.
*
*
*!Team-unique Header:
* This software is developed by the MODIS Science Data SupportTeam for 
* the National Aeronautics and Space Administration, Goddard Space 
* Flight Center, under contract NAS5-32373.
*
* Portions developed at the National Center for Supercomputing
* Applications at the Univ. of Illinois at Urbana-Champaign.
*
* Design Notes
* see PGS_MET_HDFToODL() from the SDP Toolkit. The key difference is
* that the file is already open when MPVL2ODL is called.
* 
!END********************************************************************
*/
{
  char	     buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char	     *funcname = "MPVL2ODL";
  int32	     s_id;
  int32	     attr_index; 
  int32	     loop_count; 
  int32	     num_SDS;
  int32	     nglobal_attr;
  char	     *buffer;
  char	     attr_name[MAX_NC_NAME];
  int32	     num_type; 
  int32	     count;  /*number of elements being read*/
  intn	     ret;
  size_t     buff_size;
  int        MAXILEN = 100;
  int        attrCount; /*count for family of attribute*/
  char       *suffix; /*added to handle family of attribute*/

  /* Added file handle to open and write scratch metadata file 6/17/97 */
  PGSt_IO_Gen_FileHandle *fileHandle = NULL;

  AGGREGATE 		aggNode = NULL;		/* ODL tree node */
 
  /* Input checks: */
  if ( NULLstr(PVLAttrName) )
  {
    sprintf(buff,"ERROR: MPVL2ODL unable to access an HDF PVL attribute\n"
	    "\t without the PVL attribute name input.\n");
    MAPIERR(buff,funcname);
    return(NULL);
  }
  
  if ( NULLstr(aggName) )
  {
    sprintf(buff,"ERROR: MPVL2ODL can not continue without naming the\n"
	    "\t returning tree structure.\n");
    MAPIERR(buff,funcname);
    return(NULL);
  }
  
  if ( NULLMODFIL(file) )
  {
    sprintf(buff,"ERROR: MPVL2ODL unable to access the array\n"
	    "\t with a NULL file MODFILE structure.\n");
    MAPIERR(buff,funcname);
    return(NULL);
  }
  
  /* Following added suffix and for loop to handle the family attribute.
     Almost a rewrite. Nov-28-2000. */
  
  strncpy(attr_name, PVLAttrName, MAX_NC_NAME);
  suffix = &attr_name[strlen(attr_name)];
  
  /* get the number of SDS in hdf file */
  ret = SDfileinfo((int32)file->sd_id, &num_SDS, &nglobal_attr);
  if (ret != SUCCEED)
  {
    sprintf(buff, "ERROR: MPVL2ODL unable to retrieve the number\n"
	    "\t of datasets in file %.*s\n", MAXILEN, file->filename);
    MAPIERR(buff, funcname);
    return(NULL);
  }
  
  /* Open a temporary file for scratching */
  ret = PGS_IO_Gen_Open((PGSt_PC_Logical) PGSd_MET_TEMP_ATTR_FILE, 
			PGSd_IO_Gen_Write, &fileHandle, 1);
  if (ret != PGS_S_SUCCESS)
  {
    sprintf(buff, "ERROR: MPVL2ODL detected error from SDP toolkit\n"
	    "\t procedure PGS_IO_Gen_Open while attempting to open\n"
	    "\t a scratch file.\n");
    MAPIERR(buff, funcname);
    return(NULL);
  }
  
  /* set attrCount to 0, use for suffix counting */
  attrCount = 0;
  
  /* use count, the number of the elements in attribute to control
     the loop for searching the multiple attributes. */
  count = MAX_ORDER;
  while (count == MAX_ORDER)
  {
    /* Add suffix to attribute name. no suffix for first round
       attribute search */

    if (attrCount != 0)
      sprintf(suffix, ".%d", attrCount);
    
    /* Search the global attribute first, then search the local 
       attribute */

    attr_index = SDfindattr((int32)file->sd_id, attr_name);
    s_id = (int32)file->sd_id;
    for (loop_count=0; loop_count<=(num_SDS-1) && attr_index == FAIL; 
	 loop_count++)  
    {
      if (s_id != file->sd_id)
	SDendaccess(s_id);
      
      s_id = SDselect((int32)file->sd_id, loop_count);
      if (s_id == FAIL)
	s_id = (int32)file->sd_id;
      else
	attr_index = SDfindattr(s_id, attr_name);
    } /* end of for (loop_count ...*/
    
    if (attr_index == FAIL)
    {
      sprintf(buff, "ERROR: MPVL2ODL unable to find the PVL\n"
	      "\t attribute %.*s in file %.*s\n", MAX_NC_NAME,
	      attr_name, MAXILEN, file->filename);
      MAPIERR(buff, funcname);
      count = FAIL;
    }
    else /* found the attr_index */
    {
      /* read the attribute info from hdf file */
      ret = SDattrinfo(s_id, attr_index, attr_name, &num_type, 
		       &count);
      if (ret != SUCCEED)
      {
	sprintf(buff, "ERROR: MPVL2ODL detected FAIL from HDF\n"
		"\t procedure SDattrinfo while attempting to retrieve\n"
		"\t the PVL attribute %.*s from HDF file %.*s\n", 
		MAX_NC_NAME, attr_name, MAXILEN, file->filename);
	MAPIERR(buff, funcname);
	count = FAIL;
      }
      else /* Got the attribute info. */
      {
	/* Locate the buffer for reading the attribute into it */

	buff_size = DFKNTsize(num_type)*count;
	buffer = (char *)malloc(buff_size);
	if (buffer == NULL)
	  count = FAIL;
	else
	{
	  ret = SDreadattr(s_id, attr_index, buffer);
	  if (ret != SUCCEED)
	    count = FAIL;
	  else
	  {

	    /* write the buffer into fileHandle */
	    ret = (int)fwrite(buffer, sizeof(char), count, fileHandle);
	    if (ret != count)
	    {
	      PGS_SMF_SetUNIXMsg(errno, "Temporary Metadata Attribute File",
				 funcname);
	      count = FAIL;
	    } /* end of if (fwrite ret!= count...) */
	    
	  } /* end of if (SDreadattr) */

	  free(buffer); /* free buffer for next round of reading */

	}/* end of if (buffer == NULL) */

      } /* end of if (SDattrinfo ret != SUCCEED) */

    } /* end of if (attr_index == FAIL) */

    if (s_id != file->sd_id)
      SDendaccess(s_id);
    
    attrCount ++;
  } /* end of while */
  
  /* Close the scratch file. */
  PGS_IO_Gen_Close(fileHandle);
  
  if (count == FAIL)
    return(NULL);
  
  /* Convert to ODL */
  if (PGS_MET_ConvertToOdl((char *)aggName, NULL, &aggNode) != PGS_S_SUCCESS )
  {
    sprintf(buff, "ERROR:MPVL2ODL detected FAIL from PGS procedure\n"
	    "\t PGS_MET_ConvertToOdl while attempting to\n"
	    "\t retrieve the PVL attribute %.*s from HDF\n"
	    "\t file %.*s.\n", MAX_NC_NAME, PVLAttrName, MAXILEN,
	    file->filename);
    MAPIERR(buff, funcname);
  }

  return(aggNode);
} 

