      REAL FUNCTION PLANCK_M(W, T)

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:
c    Compute monochromatic Planck radiance given brightness temperature
c    (Radiance units: Watts per square meter per steradian per micron)
c
c!INPUT PARAMETERS:
c    W (REAL)           Wavelength (microns)
c    T (REAL)           Brightness temperature (Kelvin)
c
c!OUTPUT PARAMETERS:
c    PLANCK_M (REAL)    Monochromatic Planck radiance (Watts per
c                       square meter per steradian per micron)
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
      real w, t

c ... Local variables
      double precision ws

c ... Set default return value
      planck_m = -1.0
      
c ... Check input parameters and return if they are bad
      if (w .le. 0.0 .or. t .le. 0.0) return
                  
c ... Convert wavelength to meters
      ws = 1.0d-6 * dble(w)
      
c ... Compute Planck radiance
      planck_m = sngl(1.0d-6 * (c1 / ws**5) /
     &  (exp(c2 / (ws * dble(t))) - 1.0d+0))

      END
