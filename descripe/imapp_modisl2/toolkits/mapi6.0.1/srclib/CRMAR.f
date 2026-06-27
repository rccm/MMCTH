      FUNCTION CRMAR (modfil, arrnm, grpnm, dtype, rank, dims)
C
C--------------------------------------------------------------------
C!F77
C
C!Description:    
C
C     Function CRMAR is part of a larger software system called the 
C     MODIS Applications Programming Interface (API) Utility, 
C     abbreviated M-API. The M-API Utility consists of subroutines 
C     which allow MODIS Science Team-supplied software to read and 
C     write data and metadata from/to HDF files. The functionality of 
C     the M-API is defined in the MODIS API Specification.
C
C     CRMAR creates an HDF SDS structure to store a new data array into
C     a MODIS HDF file.  It must be called before the data may be 
C     written to the file using PMAR or the associated attributes for 
C     whole array or each array dimension may (optionally) be stored 
C     using PMARIN or PMDMIN.
C
C     The grpnm string provides the facility to place the new array in 
C     an HD' Vgroup' data group.  If a Vgroup with the name grpnm 
C     does not exist, the array structure will not be created.  The 
C     array may be placed in the file outside of any Vgroup by replacing 
C     grpnm with NULL in C and a blank string (' ') in FORTRAN.
C
C     If an array with the name arrnm is written outside of a Vgroup, 
C     it must not already exist in the file.  This is to prevent the 
C     confusion caused by multiple data objects with the same name.  
C     Arrays with the same name may be stored in the same file, however,
C     if they are placed in separate Vgroups.
C
C!Input Parameters:
C  integer       modfil(3)  MODIS-HDF file structure 
C  character*(*) arrnm      Array name
C  character*(*) grpnm      Group name, if set to blank, the array 
C                           will not be placed in a data group.
C  character*(*) dtype      Data type of the array.
C                           Recommended FORTRAN data types:
C                                'CHARACTER*(*)'
C                                'INTEGER*1'
C				 'UINTEGER*1'
C				 'INTEGER*2'
C				 'UINTEGER*2'
C				 'INTEGER*4'
C				 'UINTEGER*4'
C				 'INTEGER*8'
C				 'UINTEGER*8'
C				 'REAL*4'
C				 'REAL*8'
C  integer       rank       Number of dimensions in the array.
C  integer       dim        Size of each array dimension.
C
C!Output Parameters:NONE.
C
C Returns:       MAPIOK if successful, MFAIL if an error occurs.
C
C!Revision History:
C $Log: CRMAR.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
C
C!Team-unique Header:
C
C   This software is developed by the MODIS Science Data Support Team 
C   for the National Aeronautics and Space Administration, Goddard 
C   Space Flight Center, under contract NAS5-32373.
C
C!References and Credits:
C
C    Written by   Vicky Lin        01/31/96
C    Research and Data Systems Corporation
C    SAIC/GSC MODIS SCIENCE DATA SUPPORT OFFICE
C    7501 FORBES BLVD, SEABROOK MD 20706
C
C    vlin@ltpmail.gsfc.nasa.gov
C
C!Design Notes:
C     Externals:
C               MODFILLEN               (mapi.inc)
C               ccrmar                  (mapic.inc)
C
C!END
C----------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'mapic.inc'

      CHARACTER*(*) arrnm, grpnm, dtype
      INTEGER       rank, ret, modfil(MODFILLEN), dims(*)

      call ccrmar(modfil, arrnm, len(arrnm), grpnm, len(grpnm), dtype,
     &            len(dtype), rank, dims, ret)

      CRMAR = ret
      return
      end
