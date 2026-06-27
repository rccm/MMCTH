      FUNCTION CPMFIL(modfil, mdHandles, HDFattrnms, NumHandles)

C--------------------------------------------------------------------------
C!F77 
C
C!Purpose:  Closes a new file created and accessed by the MODIS API
C           routines.
C 
C!Description:	Function CPMFIL is part of a larger software system 
C               called the MODIS Applications Programming Interface
C               (API) Utility, abbreviated M-API. The M-API Utility
C               consists of subroutines which allow MODIS Science
C               Team-supplied software to read and  write data and
C               metadata from/to HDF files. The functionality of the
C               M-API is defined in the MODIS Application Program
C               Interface (API) Specification.
C 
C               CPMFIL terminates the access of MODIS API routines
C               to a MODIS-HDF file created using openMODIS (OPMFIL). 
C
C           
C
C 
C!Input Parameters:
C         
C   modfil: Array that is used to reference a MODIS-HDF 
C           file in all other MODIS API routines.
C 
C   mdHandles: An array of character strings.  The memory size of the
C              array is [PGSd_MET_NUM_OF_GROUPS][PGSd_MET_GROUP_NAME_L],
C              where PGSd_MET_NUM_OF_GROUPS is 20 and PGSd_MET_GROUP_NAME_L
C              is 50.  This array is typedef-ed as PGSt_MET_all_handles.
C              Each row in the array stores a handles to an internal ODL
C              tree structure which will be written out a an ECS PVL
C              attribute.  Each handles, which is a string, should be
C              less than 50 characters and occupy, one row in the array.
C              Therefore, the maximum number of handles should be 20.
C
C   HDFattrnms: A character array with size of
C              [PGSd_MET_NUM_OF_GROUPS][MAX_ECS_NAME_L],
C              where PGSd_MET_NUM_OF_GROUPS is 20 and MAX_ECS_NAME_L
C              is 256.  Each string in this array is a character
C              string used as a global attribute name for storing
C              an ECS PVL text block which has a handle in the
C              corresponding row in mdHandles array. Each name,
C              which is a string, should be less that MAX_ECS_NAME_L
C              characters and occupies one row in the array.
C   
C   NumHandles: Specifies the number of actual handles contained in
C               mdHandles.  This may be set from 0 to
C               PGSd_MET_NUM_OF_GROUPS.
C
C!Output Parameters:	None.
C 
C Returns:	MAPIOK if successful, MFAIL if an error occurs. 
C 
C External References:
C           MODFILLEN                    (mapi.inc)
C           cclmfil		         (mapic.h)
C           PGSd_MET_GROUP_NAME_L        (PGS_MET.f)
C           PGSd_MET_NUM_OF_GROUPS       (PGS_MET.f)
C           MAX_ECS_NAME_L               (mapi.inc)
C
C!Revision History:
C $Log: CPMFIL.f,v $
C Revision 6.1  2010/07/13 15:45:49  kuyper
C Removed unused variables.
C
C Revision 1.2  1999/04/26 16:54:13  jayshree
C reorganizing mapi RCS
C
c Revision 1.8  1999/04/09  20:06:43  jayshree
c added len(mdHandles), len(HDFattrnms) - call to cpmfil.c
c
c Revision 1.7  1997/03/13  18:45:18  fshaw
c *** empty log message ***
c
c Revision 2.5  1996/07/23  13:29:19  fshaw
c removed printf's and changed if statemnet to test k instead of j
c
c Revision 2.4  1996/07/23  13:16:31  fshaw
c added len checking greater than 0
c
c Revision 2.3  1996/07/23  12:50:53  fshaw
c *** empty log message ***
c
c Revision 2.2  1996/07/22  21:07:58  fshaw
c added -1 to declarations to mdHandles and HDFattrnms
c
c Revision 1.1  1996/07/12  14:43:48  fshaw
c Initial revision
c
C 		
C!Team-unique Header:
C  Portions developed at the National Center for Supercomputing
C  Applications at the Univ. of Illinois at Urbana-Champaign.
C
C!References and Credits:
C  This software is developed by the MODIS Science Data Support
C  Team for the National Aeronautics and Space Administration,
C  Goddard Space Flight Center, under contract NAS5-32373.
C
C!DESIGN NOTES:
C
C!END-------------------------------------------------------------------
C
 
      IMPLICIT NONE
      INCLUDE  'mapic.inc'

      INTEGER          modfil(MODFILLEN)
      CHARACTER*(*) mdHandles(*)
      CHARACTER*(*) HDFattrnms(*)
      INTEGER  NumHandles
      INTEGER  ret

C     call ccpmfil to complete the MODIS file

C     replace first space found with NULL character


      CALL ccpmfil(modfil, mdHandles, len(mdHandles(1)), HDFattrnms, 
     &             len(HDFattrnms(1)), NumHandles, ret)
      CPMFIL = ret
      RETURN
      END
