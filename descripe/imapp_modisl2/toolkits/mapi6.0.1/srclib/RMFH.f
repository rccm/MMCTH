      FUNCTION RMFH(modfil)
C
C---------------------------------------------------------------------------
C !F77
C
C !Purpose:  To release a MAPI file id handle, writes accessed  objects to a
C           file and  frees  memory allocated by MAPI.
C
C !Description:	Function RMFH is part of a larger software system called the
C               MODIS Applications Programming Interface (M-API) utility,
C               The M-API utility consists of subroutines
C               which allow MODIS Science Team-supplied software to read and
C               write data and metadata from/to HDF files. The functionality
C               of the M-API is defined in the MODIS Application Program 
C               Interface (M-API) Specification.
C
C               Function RMFH releases the MODFIL structure created by the
C               routine  CMFH.  This routine also automatically closes all the 
C               objects openned with M-API routines.  The M-API file id handle 
C               created by calling CMFH should not be passed to the MAPI routines
C               CLMFIL or CPMFIL.
C
C               Users may open an HDF-EOS file using the HDF-EOS open file routine;
C               call CMFH to create the M-API file id handle.  Use M-API routines to
C               access data object(s) in the opened  file; after finishing object 
C               access call RMFH to release this file handle created by calling this
C               routine RMFH automatically closes all the objects opened with M-API
C               routines); finally call the  HDF-EOS closing routine to close the
C               file.
C
C !Input  parameters:
C       modfil	IN:  Array that is used to reference a MODIS-HDF file in all other
C                     MODIS API routines.  
C
C !Output   parameters:
C       None
C
C !Returns:	MAPIOK if successful, or MFAIL if an error occurs.
C
C !Revision History:
C $Log: RMFH.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c Revision 1.4  1998/02/27  15:01:16  fshaw
c fixed problem found by FORCHECK
c
c Revision 1.3  1997/12/23  22:58:15  fshaw
c correctd typos
c
c Revision 1.2  1997/12/22  21:53:57  fshaw
c corrected 12/17 W/T defects
c
C Revision 1.1  1997/11/24  15:58:46  fshaw
C Initial revision
C
C !Team_unique header:
C    This software is developed by the MODIS Science Data Support Team for
C    the National Aeronautics and Space Administration, Goddard Space
C    Flight Center, under contract NAS5-32373.
C 
C !Design Notes:
C   none.
C
C !END

      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER  modfil(MODFILLEN)

      INTEGER  ret

      CALL crmfh(modfil, ret)
      RMFH = ret
      RETURN
      END
