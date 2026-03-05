#include <stdlib.h>  
#include <stdio.h>
#include <string.h>
#include "mapic.h"
#define  ATTRIBUTE  "Attr0.0"
#include "local_nc.h"

int32 searchMODISgroup(MODFILE *file,  char const *groupname, 
                       char const *classname, char const *objectname,
                       char const *objectclass, int32 objecttype)

/*****************************************************************************
!C
 
!Description:    
        To search a Vgroup to find if an HDF object is in the specified Vgroup.

	     Function searchMODISgroup is part of a larger software system 
             called the MODIS Applications Programming Interface (API) 
             Utility, abbreviated M-API. The M-API Utility consists of 
             subroutines which allow MODIS Science Team-supplied software 
             to read  and write data and metadata from/to HDF files. The 
             functionality of the M-API is defined in the MODIS Application
             Program Interface (API) Specification.

             searchMODISgroup searches an HDF Vgroup structure in a MODIS 
             HDF file to find if an HDF object is in the Vgroup. Both the 
             group and the object are specified by their name and class 
             name. However, the classname is an optional feature. If class
             names are set to NULL, only name comparison is performed. 
             Because SDS (array) has no class name, the objectclass for an 
             SDS is always ignored. If the specified object exists, the 
             function will return the reference id for Vdata and Vgroup, 
             and index for SDS. If the object does not exist, the function
             will return NO_OBJECT, which is defined in mapic.h as -2.

!Input Parameters:
    file         MODFILE structure that is used to reference the HDF file
    groupname    ASCII string of the name of the data group.
    classname    (Optional) ASCII string of the class name of the data group. 
                 Set to NULL for not comparing the group classname.
    objectname   ASCII string of the object name to be searched.
    objectclass  (Optional)ASCII string of the class name of the data object. 
                 Set to NULL  for not comparing the object class.
    objecttype   Type of the object; The valid objects are:
		 DFTAG_NDG (for SDS)
                 DFTAG_VH  (for Vdata, or attribute if the object class 
                            is set to Attr0.0)
		 DFTAG_VG  (for Vgroup). 

!Output Parameters:
          None

Return Values: 
          NO_OBJECT if the object is not in the specified Vgroup 
          MFAIL if the function is failed
          Reference id if the object of Vdata or Vgroup object is found
          SDS index if the function successfully find the SDS object

Externally Defined:
             ATTRIBUTE                    (local_nc.h)
             CDF                          (local_nc.h)
             DFTAG_NDG                    (hdf.h)
             DFTAG_VH                     (hdf.h)
             DFTAG_VG                     (hdf.h)
             DIMENSION                    (local_nc.h)
             DIM_VALS01                   (local_nc.h)
             MAX_VAR_DIMS                 (netcdf.h)
             MAX_NC_NAME                  (netcdf.h)
             NO_OBJECT                    (mapic.h)
             NULLMODFIL                   (mapic.h)
             NULLstr                      (mapic.h)
             PGS_SMF_MAX_MSGBUF_SIZE      (PGS_SMF.h)
             SDreftoindex                 (mfhdf.h)
             SDselect                     (mfhdf.h)
             SDgetinfo                    (mfhdf.h)
             SDendaccess                  (mfhdf.h)
             VGNAMELENMAX                 (hdf.h)
             VARIABLE                     (local_nc.h)
             VSNAMELENMAX                 (hdf.h)
             Vgetid                       (vproto.h)
             Vattach                      (vproto.h)
             Vgetclass                    (vproto.h)
             Vgetname                     (vproto.h)
             Vdetach                      (vproto.h)
             Vntagrefs                    (vproto.h)
             Vgettagref                   (vproto.h)
             VSattach                     (vproto.h)
             VSgetclass                   (vproto.h)
             VSgetname                    (vproto.h)
             VSdetach                     (vproto.h)

!Revision History:
 $Log: searchMODISgroup.c,v $
 Revision 5.1  2005/04/04 18:49:11  vlin
 1. constant safe for pointer arguments.
 2. Replaced "DIM_VALS" with "DIM_VALS01".

       Xia W. Chang ,  Jan 4, 1996 , xia@modis-xl.gsfc.nasa.gov
       Version 2.0 orininal development.
	
!Team-unique Header:

       This software is developed by the MODIS Science Data Support
       Team for the National Aeronautics and Space Administration, 
       Goddard Space Flight Center, under contract NAS5-32373.

       HDF portions developed at the National Center for Supercomputing
       Applications at the Univ. of Illinois at Urbana-Champaign.

!END
*****************************************************************************/
{
  char   buff[PGS_SMF_MAX_MSGBUF_SIZE];
  char   *funcname = "searchMODISgroup";
  char   vgname[VGNAMELENMAX+1];
  char   vgclass[VGNAMELENMAX+1];
  int32  vgroup_ref;
  int32  vgroup_id, vg_id;
  int32  tag, ref_id;
  int32  index;
  int32  sds_index,sds_id,sds_rank,sds_number_type,sds_nattrs;
  int32  sds_dimsize[MAX_VAR_DIMS];
  char   sds_name[MAX_NC_NAME+1];
  int32  vdata_id;
  char   vdata_class[VSNAMELENMAX+1], vdata_name[VSNAMELENMAX+1];
  int    ret = NO_OBJECT;
  int    MAXILEN = VGNAMELENMAX;
  int    status = -1;
 
  if (NULLMODFIL(file)) 
    {
      sprintf(buff, "ERROR: searchMODISgroup unable to continue with an"
	      "empty file pointer\n");
      MAPIERR(buff, funcname);
      return(MFAIL);
    }

  if (NULLstr(groupname))
    {
      sprintf(buff, "ERROR: searchMODISgroup unable to continue with an"
	      "empty groupname\n");
      MAPIERR(buff, funcname);
      return(MFAIL);
    }
  
  if (NULLstr(objectname))
    {
      sprintf(buff, "ERROR: searchMODISgroup unable to search group %.*s "
	      "\n\t with an empty object name\n", VGNAMELENMAX, groupname);
      MAPIERR(buff, funcname);
      return(MFAIL);
    }
  
  if ((objecttype != DFTAG_NDG) && (objecttype != DFTAG_VH) &&
      (objecttype != DFTAG_VG))
    {
      sprintf(buff, "ERROR: searchMODISgroup unable to search object %.*s in"
	      "\n\t group %.*s with an unrecognized object type %ld \n", 
	      MAXILEN, objectname, MAXILEN, groupname, (long)objecttype);
      MAPIERR(buff, funcname);
      return(MFAIL);
    }

  /* search for the Vgroup first   */
  vgroup_ref = -1;
  /* Loop calling Vgetid for each Vgroup reference id. */
  while((vgroup_ref = Vgetid( (int32)file->hdf_id, vgroup_ref)) != FAIL)
    {
      /* get the Vgroup ID    */
      if ((vgroup_id = Vattach( (int32)file->hdf_id, vgroup_ref, "r")) == FAIL)
	{
	  sprintf(buff, "ERROR: searchMODISgroup fails to search object "
		  "%.*s \n\t in group %.*s becaure Vattach fails\n",
		  MAXILEN, objectname, MAXILEN, groupname);
	  MAPIERR(buff, funcname);
	  return(MFAIL);
	}

      Vgetclass(vgroup_id, vgclass);       /* get Vgroup classname */

      if ((classname != NULL) || !( (strcmp(vgclass,VARIABLE) == 0) 
			        || (strcmp(vgclass,DIMENSION) == 0)
				|| (strcmp(vgclass,CDF) == 0) ) )
	  {
	    Vgetname(vgroup_id, vgname);
	    if (strcmp(vgname, groupname) == 0)
	      {
		if ((classname == NULL) || (strcmp(vgclass,classname) == 0))
		  break;
	      }
	  }
      Vdetach(vgroup_id);
    }                              /* end of while loop */


  if (vgroup_ref == FAIL)
    {
      sprintf(buff, "ERROR: searchMODISgroup unable to find the specified"
	      "\n\t Vgroup group %.*s", MAXILEN, groupname);
      MAPIERR(buff, funcname);
      return(MFAIL);
    }

  for(index=0; index<=(Vntagrefs(vgroup_id)-1); index++)
    { 
      /* get the tag and reference id */
      status = Vgettagref(vgroup_id, index, &tag, &ref_id);
      if (status == FAIL)
	{
	  sprintf(buff, "ERROR: searchMODISgroup fails to obtain %.*s's"
		  " tag and reference number\n", MAXILEN, objectname);
	  MAPIERR(buff, funcname);
	  Vdetach(vgroup_id);
	  return(MFAIL);
	}

      if (tag == objecttype)
	{
	  switch (objecttype) {
	  case DFTAG_NDG:
	    if ((sds_index = SDreftoindex((int32)file->sd_id, ref_id)) != FAIL)
	      {
		sds_id = SDselect((int32)file->sd_id, sds_index);
		status = SDgetinfo(sds_id,sds_name, &sds_rank, sds_dimsize,
				   &sds_number_type, &sds_nattrs);
		/* release the internal data structure  */
		status = SDendaccess(sds_id);
		if (strcmp(sds_name, objectname) == 0)
		  ret = sds_index;
	      }
	    else
	      {
		sprintf(buff, "Warning: Vgroup \"%.*s\" may contain non-"
			"existing SDS object with reference id %ld .\n", 
			MAXILEN, groupname, (long int)ref_id);
		MAPIWARN(buff, funcname);
	      }
	    break;
	  case DFTAG_VH:
	    if ((vdata_id = VSattach( (int32)file->hdf_id,ref_id, "r")) != FAIL)
	      {
		VSgetclass(vdata_id, vdata_class);  /* get the class name */
		if ((objectclass!=NULL) || ((strcmp(vdata_class,ATTRIBUTE)!=0)
		    && (strcmp(vdata_class, DIM_VALS01) != 0)))
		  {
		    VSgetname(vdata_id, vdata_name); /* get the vdata name */
		    if (strcmp(vdata_name, objectname) == 0)
		      if ((objectclass == NULL) || 
			  (strcmp(vdata_class,objectclass) == 0))
			ret = ref_id;
		  }
		VSdetach(vdata_id);
	      }
	    else
	      {
		sprintf(buff, "Warning: Vgroup \"%.*s\" contains non-exist"
			" Vdata object with reference id %ld .\n", 
			MAXILEN, groupname, (long int)ref_id);
		MAPIWARN(buff, funcname);
	      }
	    break;
	  case DFTAG_VG:
	    if ((vg_id = Vattach( (int32)file->hdf_id,ref_id,"r")) != FAIL)
	      {
		Vgetclass(vg_id, vgclass);     /* get the classname */
		if ((objectclass != NULL) || 
		    ((strcmp(vgclass, VARIABLE) != 0) &&
		     (strcmp(vgclass, DIMENSION) != 0) &&
		     (strcmp(vgclass, CDF) != 0)))
		  {
		    Vgetname(vg_id, vgname);
		    if (strcmp(vgname, objectname) == 0)
		      if ((objectclass == NULL) ||
			  (strcmp(vgclass, objectclass) == 0))
			ret = ref_id;
		  }
		Vdetach(vg_id);
	      }
	    else
	      {
		sprintf(buff, "Warning: Vgroup \"%.*s\" contains non-exist"
			" Vgroup object with reference id %ld .\n", 
			MAXILEN, groupname, (long int)ref_id);
		MAPIERR(buff, funcname);
	      }
	    break;
	  }              /* end of switch       */
	  
	}                /* end of if tag ==    */
	  
      if (ret != NO_OBJECT) break;
    }                   /* end of for loop     */

  Vdetach(vgroup_id);
  return(ret);
}
