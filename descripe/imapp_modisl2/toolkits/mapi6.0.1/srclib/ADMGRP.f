      FUNCTION ADMGRP(modfil,grpnm,clsnm,tag,ref)
C
C-----------------------------------------------------------------------
C!F77
C
C!Purpose:	To add an object to an HDF Vgroup.
C 
C!Description:	Function ADMGRP is part of a larger software system called the 
C MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API. The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read  and 
C write data and metadata from/to HDF files. The functionality of 
C the M-API is defined in the MODIS Application Program Interface 
C (API) Specification.
C 
C ADMGRP is an MAPI internal function which adds an HDF object into 
C a Vgroup. The group is specified by its name and class name. 
C However, the classname is an optional feature. If class names are set 
C to NULL, the function will find the Vgroup with the same name as 
C groupname without considering the classname of the group to insert 
C the object into the group. The object is specified by its tag and 
C reference number. The tag for SDS is DFTAG_NDG and for Vdata is 
C DFTAG_VH. To obtain the reference id, using sfidtoref for SDS and 
C VSQueryref for Vdata. No check for the tag or reference is 
C performed.
C 
C!Input Parameters:
C 
C modfil	IN: 	Array that is used to reference the MODIS-HDF 
C			 file.
C grpnm		IN:	ASCII string of the name of the data group.
C clsnm		IN:	(Optional) ASCII string of the class name of 
C 			the data group. Set to NULL for not comparing the 
C 			group classname.
C tag		IN:	The tag of the object. Valid values 
C 			include:DFTAG_NDG for SDS, DFTAG_VH for Vdata, 
C 			and DFTAG_VG for Vgroup. 
C ref		IN:	The reference number of the object.
C 
C!Output Parameters:.:		none
C 
C Returns:	MAPIOK if successful, otherwise, MFAIL.
C 
C External references:
C			MODFILLEN		(mapi.inc)
C			cadmgrp			(mapic.h)
C!Revision History:
C $Log: ADMGRP.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.1  1996/01/29  19:07:36  qhuang
c Initial revision
C 
C!Team-unique header:
C       This software is developed by the MODIS Science Data Support Team for 
C       the National Aeronautics and Space Administration, Goddard Space 
C       Flight Center, under contract NAS5-32373.
C
C !References and Credits:
C
C !Design Notes:
C
C----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER          modfil(MODFILLEN), tag,ref
      CHARACTER*(*)    grpnm, clsnm
      
      INTEGER  ret

      CALL cadmgrp(modfil,grpnm,len(grpnm),clsnm,len(clsnm),tag,ref,ret)
      ADMGRP = ret
      RETURN
      END
