      FUNCTION CEIL(X)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C     Return the smallest integer not less than X
C
C!INPUT PARAMETERS:
C     X        Floating point scalar
C
C!OUTPUT PARAMETERS:
C     CEIL     Smallest integer not less than X
C
C!REVISION HISTORY:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C!DESIGN NOTES:
C     CEIL(3.7) has the value 4
C     CEIL(-3.7) has the value -3
C
C!END
C-----------------------------------------------------------------------

      implicit none

      integer ceil
      real x

      ceil = int(x)
      if (real(ceil) .lt. x) ceil = ceil + 1

      END
