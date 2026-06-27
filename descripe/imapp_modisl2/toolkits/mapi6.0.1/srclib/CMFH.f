       FUNCTION CMFH (swfid, modfil)
C
C  !F77
C
C  !Purpose: To create a M-API file id handle from a HDF-EOS file id handle.
C
C  !Description: Function CMFH is part of a larger software system called the
C      MODIS Applications Programming Interface (API) Utility, abbreviated M-API.
C      The M-API Utility consists of subroutines which allow MODIS Science
C      Team-supplied software to read and write data and metadata from/to HDF files.
C      The functionality of the M-API is defined in the MODIS Application Program
C      Interface (API) Specification.
C
C      Function CMFH provides a M-API file handle for a previously opened HDF-EOS
C      file (swath, grid, or point).  The routine RMFH must be called for each call
C      to this  routine. NOTE: This file id handle should not be passed to M-API
C      routine CLMFIL or CPMFIL.  The file must be opened and closed using either
C      HDF-EOS or M-API.  These can not be intermixed.
C
C      Users may open an HDF-EOS file using the HDF-EOS open file routine.  Call
C      CMFH to create the M-API file handle.  Use M-API routines to access data
C      object(s) in the opened file.  Once one is finished  accessing the file, 
C      call RMFH to release this file handle created by calling this routine. 
C      (It also writes all the objects opened with the M-API routines to disk
C      and closes them.)  The final step before exiting the program is to call
C      the HDF-EOS close file routine to close the  HDF-EOS file.
C
C
C !Input parameter:
C	swfid    IN: HDF-EOS file id handle.
C 
C !Output parameter:
C       modfil  OUT: Array that is used to reference the file in all other MODIS API
C       routines.  The array will  contain all zeroes if an error occurs.
C
C !Returns:  MAPIOK if successful, MFAIL if an error occurs. 
C
C !Revision History: 
C      Frederick J. Shaw  (fshaw@ltpmail.gsfc.nasa.gov)
C      General Sciences Corp.
C    
C   $Log: CMFH.f,v $
C   Revision 1.2  2001/01/24 19:05:38  pliu
C   Added !F77 in prolog.
C   Added RCS log keyword.
C
C Revision 1.1  1999/04/15  19:26:51  jayshree
C Initial revision 
C
C !Team-unique header:
C      This software is developed by the MODIS Science Data Support Team for the
C      National Aeronautics and Space Administration, Goddard Space Flight Center,
C      under contract NAS5-32373.
C 
C !Design Notes:
C      swfid is declared as an INTEGER.  However, in the HDF-EOS documentation it
C      is declared as INTEGER*4.  The routine ccmfh is expecting a pointer to an
C      INTEGER.  
C 
C !END

      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER  swfid 
      INTEGER  modfil(MODFILLEN)
      INTEGER  ret

      CALL ccmfh(swfid, modfil, ret) 
      CMFH = ret
      RETURN
      END
