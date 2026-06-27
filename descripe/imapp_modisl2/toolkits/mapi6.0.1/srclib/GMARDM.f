      FUNCTION GMARDM(modfil,arrnm,grpnm,dtype,rank,dims)
C
C--------------------------------------------------------------------
C!F77
C
C!Description:    
C
C   GMARDM retrieves the essential characteristics (the rank, dimensions
C   and data type) of an HDF SDS array structure contained in a MODIS-
C   HDF file.  This provides the information needed for properly reading
C   data from the array structure using GMAR.
C
C   The groupname  string provides the facility to select an array 
C   structure placed in a particular HDF Vgroup data group.  
C   Alternatively, the entire file will be searched for an array 
C   structure named arrayname if groupname is a blank string (' ').
C
C   Proper dimensioning of dimsizes to provide sufficient elements for
C   the dimensions of the array structure may at first appear to require
C   precognition.  The easiest solution is to provide a generous 
C   dimsizes array.  Another approach is to use the rank variable as an
C   input containing the number of elements in dimsizes.  If dimsizes is
C   inadequate for the multi-dimensional array structure in question, 
C   getMODISardims will fail gracefully but will return the rank of the
C   array structure, allowing for the dimension information to be 
C   retrieved on the second call.
C
C!Input Parameters:
C   modfil   Array that is used to reference the MODIS-HDF file 
C            containing the target array structure.
C   arrnm    ASCII string name of the target array structure.
C   grpnm    ASCII string name of the data group containing the target 
C            array structure.  If grpnm = ' ', the entire file will be
C            searched for the array structure named arrayname.
C
C!Output Parameters:
C   dtype         Data type of the array.  
C                 Possible FORTRAN data types are:
C                 'INTEGER*1', 'UINTEGER*1', 'INTEGER*2', 'UINTEGER*2',
C                 'INTEGER*4', 'UINTEGER*4', 'INTEGER*8', 'UINTEGER*8',
C                 'REAL*4', and 'REAL*8'
C   rank  IN/OUT  The number of elements in the array dimensions on 
C                 input.  This is replaced with the number of dimensions 
C                 in the array structure for output.  It is set to 0 if
C                 a functional error occures.  No dimensional 
C                 information will be provided if rank = 0.
C   dimsizes      Array describing the size of each dimension of the the
C                 array structure.  The dimensions will not be provided 
C                 unless dimsizes contains sufficient elements for the 
C                 rank of the array.  If the dimsizes is too small, the
C                 array rank information will be return in the rank 
C                 argument and the function will return MFAIL.
C
C!Revision History:
C
C $Log: GMARDM.f,v $
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
C   Written by   Vicky Lin        02/06/96
C   Research and Data Systems Corporation
C   SAIC/GSC MODIS SCIENCE DATA SUPPORT OFFICE
C   7501 FORBES BLVD, SEABROOK MD 20706
C
C   vlin@ltpmail.gsfc.nasa.gov
C
C!Design Notes:
C
C   Externals:  MODFILLEN               (mapi.inc)
C               GMARDM                  (mapi.inc)
C               cgmardm                 (mapic.inc)
C
C   Returns:  MAPIOK if successful, MFAIL if an error occurs. 
C
C----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      include 'mapic.inc'
      character*(*) arrnm, grpnm, dtype
      integer       ret, rank, modfil(MODFILLEN), dims(*)

      call cgmardm(modfil,arrnm,len(arrnm),grpnm,len(grpnm),dtype,
     &     len(dtype),rank,dims,ret)
      GMARDM = ret
      return
      end
