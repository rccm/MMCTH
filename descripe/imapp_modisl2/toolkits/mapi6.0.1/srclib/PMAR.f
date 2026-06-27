      FUNCTION PMAR(modfil, arrnm, grpnm, start, dims, ddata) 

C----------------------------------------------------------------------
C!F77
C
C!Purpose: Write an array or subarray of data to MODIS-HDF file array
C          structure.
C
C!Description:
C
C   Function PMAR is part of a larger software system called the MODIS
C   Applications Programming Interface (API) Utility, abbreviated M-API.
C   The M-API Utility consists of subroutines which allow MODIS Science
C   Team supplied software to read and write data and metadata from/to
C   HDF files.  The functionality of the M-API is defined in the MODIS
C   Application Program Interface Specification.
C
C   PMAR places a multi-dimensional array of data into an HDF SDS array
C   structure previously created using CRMAR.  The data in the array
C   must be of the data type the target array structure was created for.
C   In addition, the dimensions and placement of the input array in the
C   array structure must be consistent with the structure's rank and
C   dimensions.  If a PMAR error message occurs, the data insertion will
C   not be performed.  This routine may be called multiple times to fill
C   the array structure.  Data previously in the array structure will
C   be overwritten.
C
C   The groupname string provides the facility to select an array
C   structure placed in a particular HDF Vgroup data group.  The entire
C   file will be searched for an array structure named arrayname if
C   groupname = NULL in C or a blank string (' ') in FORTRAN.
C
C !Input Parameters:
C   modfil   IN:  Array used to reference the MODIS-HDF file
C   arrnm    IN:  The name of the array.
C   grpnm    IN:  Name of the data group to place the array,  If grpnm = ' ',
C            the entire file will be searched for the array structure
C            named arrayname.
C   start    IN:  Array containing the array structure location to begin
C            placing data into the array structure.  Start must have the
C            same number of elements as the target array has dimensions.
C   dims     IN:  Array describing the size of the array being inserted into
C            the array structure.  dims must have the same number of
C            elements as the target array structure has dimensions and
C            the product of the array dimensions must equal the number of
C            elements in data.
C   ddata     IN:  Multi-dimensional data array.
C
C !Output Parameters:NONE
C
C Returns:	MAPIOK if successful, MFAIL if an error occurs.
C
C Externals referece: 
C		MODFILLEN                (mapi.inc)
C		cpmar			 (mapic.h)
C
C !Revision History:
C		Qi Huang	1996/08/19
C		Version 2.1
C		Original development and testing
C		
C		Ring super structrue and other changes make this
C		version much faster.
C   $Log: PMAR.f,v $
C   Revision 1.1  1999/04/15 19:26:51  jayshree
C   Initial revision
C
c Revision 1.2  1996/08/26  21:20:09  qhuang
c Version 2.1
c
c Revision 1.1  1996/08/26  21:19:50  qhuang
c Initial revision
c
C
C !Team-unique Header:
C
C   This software is developed by the MODIS Science Data Support
C   Team for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C   Portions developed at the National Center for Supercomputing
C   Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !References and Credits:
C		Qi Huang	1996/08/19
C		Version 2.1
C		Original development and testing
C
C!Design Notes:
C
C!END
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER       modfil(MODFILLEN), start(*), dims(*)
      CHARACTER*(*) arrnm, grpnm
      BYTE          ddata(*)
      INTEGER       ret
    
      CALL cpmar(modfil,arrnm,len(arrnm),grpnm,len(grpnm),
     +           start,dims,ddata,ret)
   
      PMAR = ret
      RETURN
      END
