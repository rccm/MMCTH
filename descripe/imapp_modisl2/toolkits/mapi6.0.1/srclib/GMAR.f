        FUNCTION GMAR (modfil, arrnm, grpnm, start, dims, ddata)
C
C----------------------------------------------------------------------
C !F77
C
C !Description: 
C
C   Function GMAR is part of a larger software system called the MODIS 
C   Applications Programming Interface (API) Utility, abbreviated M-API.
C   The M-API Utility consists of subroutines which allow MODIS Science 
C   Team supplied software to read and write data and metadata from/to
C   HDF files.  The functionality of the M-API is defined in the MODIS 
C   Application Program Interface Specification.
C
C   GMAR returns a multi-dimensional array of data from an HDF SDS
C   array structure contained in a MODIS-HDF file.  The data array must 
C   be of the same data type as data in the target array structure.  
C   In addition, the dimensions and array region requested from the 
C   array structure must be consistent with the structure's rank and 
C   dimensions.  (The array structure's data type, rank, and dimensions 
C   may be retrieved using GMARIN).  If a GMAR error message occurs the
C   data retrieval will not be performed. 
C
C   The groupname  string provides the facility to select an array 
C   structure placed in a particular HDF Vgroup data group.  
C   Alternatively , the entire file will be searched for an array 
C   structure named arrayname  if groupname = NULL in C or a blank 
C   strin (' ') n FORTRAN.
C
C !Input Parameters:
C   modfil    Array used to reference the MODIS-HDF file 
C   arrnm     ASCII string that will be the name of the array.
C   grpnm     ASCII string name of the data group containing the 
C             target array structure.  If grpnm = ' ' the entire file 
C             will be searched for the array structure named arrnm.
C   start     Array containing the array structure  location to begin 
C             reading the data from the array structure.  start must 
C             have the same number of elements as the target array has 
C             dimensions.
C   dims      Array describing the size of the array being retrieved 
C             from ther array structure.  dims must have the same number
C             of elements as the target array structure has dimensions 
C             and the product of the array dimensions must equal the 
C             number of elements in data.
C
C !Output Parameters:
C   ddata      Multi-dimensional data array.
C
C Returns: MAPIOK if successful, MFAIL if an error occurs.
C
C Externals:
C     MODFILLEN                (mapi.inc)
C     cgmar                    (mapic.inc)
C
C !Revision History:
C   $Log: GMAR.f,v $
C   Revision 1.1  1999/04/15 19:26:51  jayshree
C   Initial revision
C 
C                Frederick J. Shaw Sept 9, 1996
C                Orginal development and Testing
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
C
C !Design Notes:
C
C !END
C----------------------------------------------------------------------
C
        IMPLICIT NONE
        INCLUDE 'mapic.inc'
C Declaration of variables

      byte          ddata(*)
      character*(*) arrnm, grpnm
      integer       modfil(MODFILLEN), start(*), dims(*)  
      integer       ret

      CALL cgmar(modfil, arrnm, len(arrnm), grpnm, len(grpnm),
     &        start, dims, ddata, ret)

      GMAR = ret
      return
      end
