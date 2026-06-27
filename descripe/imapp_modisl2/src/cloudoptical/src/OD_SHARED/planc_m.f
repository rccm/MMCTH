      REAL FUNCTION PLANC_M(V, T)

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:
c    Compute monochromatic Planck radiance given brightness temperature
c    (Radiance units: milliWatts per square meter per steradian per
c    inverse centimeter)
c
c!INPUT PARAMETERS:
c    V (REAL)          Wavenumber (inverse centimeters)
c    T (REAL)          Brightness temperature (Kelvin)
c
c!OUTPUT PARAMETERS:
c    PLANC_M (REAL)    Monochromatic Planck radiance (milliWatts per
c                      square meter per steradian per
c                      inverse centimeter)
c
c!REVISION HISTORY:
c
c!TEAM-UNIQUE HEADER:
c    Liam.Gumley@ssec.wisc.edu
c
c!END
c-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Include files
      include 'fundamental_constants.inc'

c ... Arguments
      real v, t

c ... Local variables
      double precision vs

c ... Set default return value
      planc_m = -1.0
      
c ... Check input parameters and return if they are bad
      if (v .le. 0.0 .or. t .le. 0.0) return
                  
c ... Convert wavenumber to inverse meters
      vs = 1.0d+2 * dble(v)
      
c ... Compute Planck radiance
      planc_m = sngl(1.0d+5 * (c1 * vs**3) /
     &  (exp(c2 * vs / dble(t)) - 1.0d+0))
            
      END
