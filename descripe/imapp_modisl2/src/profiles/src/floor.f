      FUNCTION FLOOR(X)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C     Return the largest integer not greater than X
C
C!INPUT PARAMETERS:
C     X        Floating point scalar
C
C!OUTPUT PARAMETERS:
C     FLOOR    Largest integer not greater than X
C
C!REVISION HISTORY:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C!DESIGN NOTES:
C     FLOOR(3.7) has the value 3
C     FLOOR(-3.7) has the value -4
C
C!END
C-----------------------------------------------------------------------

      implicit none

      integer floor
      real x

      floor = int(x)
      if (real(floor) .gt. x) floor = floor - 1

      END
