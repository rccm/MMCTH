      FUNCTION PMDNAM(modfil, arrnm, grpnm, moddim, dnm) 
C-----------------------------------------------------------------------
C!F77
C!Purpose:	Assigns a name to a specific dimension of an array structure.
C 
C!Description: 	Function PMDNAM is part of a larger software system called the 
C MODIS Applications Programming Interface (API) Utility, 
C abbreviated M-API.  The M-API Utility consists of subroutines 
C which allow MODIS Science Team-supplied software to read  and 
C write data and metadata from/to HDF files. The functionality of 
C the M-API is defined in the MODIS Application Program Interface 
C (API) Specification.
C 
C PMDNAM (PMDNAM) associates an HDF dimension name with a 
C specified SDS array structure dimension. The SDS array must be 
C created (using createMODISarray) before it is possible to name 
C any of its dimensions. This routine does not create a "long_name" 
C dimension attribute. PMDNAM can produce such a dimension 
C label, however.
C 
C PMDNAM does more than apply an appellation to a dimension. An 
C HDF dimension name is an independent data object. It may be 
C shared by several array structure dimensions, but they all must 
C be of the same size. Any dimension attribute that is associated 
C with any one of these dimensions is immediately associated with 
C all the dimensions having that name. Likewise, updating a 
C dimension attribute for one dimension updates it for all 
C dimensions having the same name (they could only have one 
C "long name" dimension shared between them).
C 
C Naming an SDS dimension will also cause any dimension 
C attributes currently associated with that dimension to be lost. 
C Therefore it is most practical to name an array's dimensions, if 
C necessary, immediately after the array structure's creation and 
C before creating dimension attributes for it.
C 
C The grpnm string provides the facility to select an array 
C structure placed in a particular HDF 'Vgroup' data group. 
C Alternatively, the entire file will be searched for an array 
C structure named arrnm if the argument grpnm = NULL in C or 
C grpnm is a blank string ( " ") in FORTRAN.
C 
C!Input Parameters:
C modfil	IN: 	Array  that is used to reference the MODIS 
C 		HDF file receiving the dimension name.
C arrnm	IN:	ASCII string name of the target array 
C 		structure.
C grpnm	IN:	ASCII string name of the data group 
C		containing the target array structure.  If grpnm = ' ' 
C		the entire file will be searched for the array 
C		structure named arrnm.
C moddim	IN:	The dimension to be named (0-based). The 0 
C		dimension of an HDF SDS array structure is 
C		associated with the most rapidly varying array 
C		index.
C dnm	IN:	ASCII string name to gove to the dimension.
C 
C!Output Parameters:		NONE
C 
C Returns:	MAPIOK if successful, MFAIL if an error occurs.
C
C External reference:
C		MODFILLEN	(mapi.inc)
C 		cpmdnam		(mapic.h)
C
C!Revision History:
C		Qi Huang	1996/0819
C		Version 2.1
C		Orginal development and testing
C
C		Ring super structure and other changes make this
C		version much faster.
C $Log: PMDNAM.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.2  1996/08/26  21:08:58  qhuang
c Version 2.1
c
c Revision 1.1  1996/08/26  21:07:52  qhuang
c Initial revision
c
C 
C!Team-unique Header:
C This software is developed by the MODIS Science Data Support Team for 
C the National Aeronautics and Space Administration, Goddard Space 
C Flight Center, under contract NAS5-32373.
C 
C !References and Credits:
C		Qi Huang	1996/0819
C		Version 2.1
C		Orginal development and testing
C
C !Design Notes:
C
C----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER          modfil(MODFILLEN), moddim
      CHARACTER*(*) arrnm, grpnm, dnm
      
      INTEGER  ret

      CALL cpmdnam(modfil,arrnm,len(arrnm),grpnm,len(grpnm),
     +             moddim, dnm,len(dnm),ret)
      PMDNAM = ret
      RETURN
      END
