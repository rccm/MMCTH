      REAL FUNCTION BRIGHT_M(W, R)

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:
c    Compute brightness temperature given monochromatic Planck radiance
c    (Radiance units: Watts per square meter per steradian per micron)
c
c!INPUT PARAMETERS:
c    W (REAL)           Wavelength (microns)
c    R (REAL)           Monochromatic Planck radiance (Watts per
c                       square meter per steradian per micron)
c
c!OUTPUT PARAMETERS:
c    BRIGHT_M (REAL)    Brightness temperature (Kelvin)
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
      real w, r

c ... Local variables
      double precision ws

c ... Set default return value
      bright_m = -1.0
      
c ... Check input parameters and return if they are bad
      if (w .le. 0.0 .or. r .le. 0.0) return
                  
c ... Convert wavelength to meters
      ws = 1.0d-6 * dble(w)
      
c ... Compute brightness temperature
      bright_m = sngl(c2 /
     &  (ws * log(c1 / (1.0d+6 * dble(r) * ws**5) + 1.0d+0)))

      END
