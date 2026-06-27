        FUNCTION GMDNAM(modfil, arrnm, grpnm, moddim, dnam)
        IMPLICIT NONE
        INCLUDE 'mapic.inc'
C----------------------------------------------------------------------
C !F77
C
C !Purpose:     Obtains the name of a specific dimension in an array structure.
C
C !Description: Function GMDNAM is part of a larger software system called the
C               MODIS Applications Programming Interface (API) Utility,
C               abbreviated M-API.  The M-API Utility consists of subroutines
C               which allow MODIS Science Team-supplied software to read and
C               write data and metadata from/to HDF files. The functionality
C               of the M-API is defined in the MODIS Application Program
C               Interface (API) Specification.
C
C               GMDNAM  retrieves the name of an HDF dimension associated
C               with an array structure given the array's name and the
C               dimensions number. If the dimension name cannot be found,
C               the routine will return MFAIL (-1). This routine does not
C               retrieve a " long_name" dimension attribute. getMODISdiminfo
C               (GMDMIN) can retrieve such a dimension label (if it exists),
C               however.
C
C               The grpnm string provides the facility to select an array
C               structure placed in a particular HDF 'Vgroup' data group.
C               Alternatively, the entire file will be searched for an array
C               structure named arrnm if the argument is a blank string (' ')
C               in FORTRAN.
C
C
C  !Input parameters:
C    
C     modfil  IN: Array that is used to reference the MODIS HDF file containing
C               the dimension name.
C
C     arrnm   IN: ASCII string name of the target array structure.
C
C     grpnm   IN: ASCII string name of the data group containing the target
C                 array
C               structure. If grpnm = ' ' the entire file will be searched for
C               the array structure named arrnm.
C     moddim     IN: The dimension number which the dimension name is attached to
C               (0-based). The 0 dimension of an HDF SDS array structure is
C               associated with the most rapidly varying array index.
C
C   !Output parameters:
C
C     dname   OUT: ASCII string for the dimension name. Provided array should be
C                at least 256 bytes long.
C
C   !Returns:        MAPIOK if successful, MFAIL if an error occurs.
C
C   Externals:
C
C     MODFILLEN                (mapi.inc)
C     MODFILE                  (mapi.inc)
C     cgmdnam                  (mapic.inc)
C
C !Revision History:
C   $Log: GMDNAM.f,v $
C   Revision 1.1  1999/04/15 19:26:51  jayshree
C   Initial revision
C
c Revision 1.1  1996/08/21  17:09:23  fshaw
c Initial revision
c 
C
C !TEAM-UNIQUE HEADER:
C
C   This software is developed by the MODIS Science Data Support
C   Team for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C   Portions developed at the National Center for Supercomputing
C   Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !REFERENCES AND CREDITS
C
C !DESIGN NOTES:
C
C !END
C----------------------------------------------------------------------

C Declaration of variables

      INTEGER modfil(MODFILLEN), moddim, ret
      CHARACTER*(*) arrnm, grpnm, dnam

      CALL cgmdnam(modfil, arrnm, len(arrnm), grpnm, len(grpnm),
     &        moddim, dnam, len(dnam), ret)
      GMDNAM = ret
      return
      end
