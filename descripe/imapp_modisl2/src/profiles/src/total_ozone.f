      SUBROUTINE TOTAL_OZONE(P, O3PPMV, NL, LS, TOTAL)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Compute total column ozone given a top-down profile of
C     ozone partial pressures.
C
C !INPUT PARAMETERS:
C     P        Pressure profile (millibars)
C     O3PPMV   Ozone profile (parts per million by volume)
C     NL       Number of levels
C     LS       Surface level in profile
C
C !OUTPUT PARAMETERS:
C     TOTAL    Total column ozone (Dobsons)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Parameters
      real oct
      parameter (oct = 0.78961)

c ... Scalar arguments
      real total
      integer ls, nl

c ... Array arguments
      real o3ppmv(nl), p(nl)

c ... Local scalars
      real phd1, phd2, tem1, tem2
      integer ip

c ... Integrate the profile
      total = 0.0
      phd1 = o3ppmv(1)
      tem1 = p(1)
      do ip = 2, ls
          phd2 = o3ppmv(ip)
          tem2 = p(ip)
          total = total + 0.5 * (phd1 + phd2) * (tem2 - tem1)
          phd1 = phd2
          tem1 = tem2
      end do
      total = total * oct

      END
