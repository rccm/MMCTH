module platform
  implicit none
  character*5 platform_name
  parameter (platform_name = 'TERRA')

end module platform


module planck_bright
 use platform

 implicit none
 private
 public modis_bright, bright_m, modis_planck, &
        brite_m, planc_m, planck_m
! ... Planck constant (Joule second)
  double precision h
  parameter (h = 6.62606876d-34)

! ... Speed of light in vacuum (meters per second)
  double precision c
  parameter (c = 2.99792458d+08)

! ... Boltzmann constant (Joules per Kelvin)
  double precision k
  parameter (k = 1.3806503d-23)

! ... Derived constants
  double precision c1, c2
  parameter (c1 = 2.0d+0 * h * c * c)
  parameter (c2 = (h * c) / k)

! ... Platform name hardcoded to Terra.
  ! character*5 platform_name
  ! parameter (platform_name = 'TERRA')

  real cwn_terra(16), tcs_terra(16), tci_terra(16)

! c     TERRA MODIS DETECTOR-AVERAGED SPECTRAL RESPONSE
! c     (LIAM GUMLEY 2003/06/05)
!
! c     BAND 20 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
! c     BAND 21 TEMPERATURE RANGE WAS  180.00 K TO  400.00 K
! c     BAND 22 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
! c     BAND 23 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
! c     BAND 24 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
! c     BAND 25 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
! c     BAND 27 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
! c     BAND 28 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
! c     BAND 29 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
! c     BAND 30 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
! c     BAND 31 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
! c     BAND 32 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
! c     BAND 33 TEMPERATURE RANGE WAS  180.00 K TO  330.00 K
! c     BAND 34 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
! c     BAND 35 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K
! c     BAND 36 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K

! c     BANDS
! c      20,  21,  22,  23,
! c      24,  25,  27,  28,
! c      29,  30,  31,  32,
! c      33,  34,  35,  36,

! ... Effective central wavenumbers (inverse centimeters)
  data cwn_terra/ &
   2.641767E+03, 2.505274E+03, 2.518031E+03, 2.465422E+03,&
   2.235812E+03, 2.200345E+03, 1.478026E+03, 1.362741E+03,&
   1.173198E+03, 1.027703E+03, 9.081998E+02, 8.315149E+02,&
   7.483224E+02, 7.309089E+02, 7.188677E+02, 7.045309E+02/

! ! ... Temperature correction slopes (no units)
  data tcs_terra/&
   9.993487E-01, 9.998699E-01, 9.998604E-01, 9.998701E-01,&
   9.998825E-01, 9.998849E-01, 9.994942E-01, 9.994937E-01,&
   9.995643E-01, 9.997499E-01, 9.995880E-01, 9.997388E-01,&
   9.999192E-01, 9.999171E-01, 9.999174E-01, 9.999264E-01/

! !c ... Temperature correction intercepts (Kelvin)
  data tci_terra/&
   4.744530E-01, 9.091094E-02, 9.694298E-02, 8.856134E-02,&
   7.287017E-02, 7.037161E-02, 2.177889E-01, 2.037728E-01,&
   1.559624E-01, 7.989879E-02, 1.176660E-01, 6.856633E-02,&
   1.903625E-02, 1.902709E-02, 1.859296E-02, 1.619453E-02/

contains

  real function modis_bright(rad, band, units)
    !Original description.
    ! C
    ! C!DESCRIPTION:
    ! C    Compute brightness temperature for a MODIS infrared band
    ! C    on Terra or Aqua.
    ! C
    ! C    Spectral responses for each IR detector were obtained from MCST:
    ! C    ftp://ftp.mcst.ssai.biz/pub/permanent/MCST in directories
    ! C    PFM_L1B_LUT_4-30-99 (Terra) and FM1_RSR_LUT_07-10-01 (Aqua).
    ! C
    ! C    An average spectral response for all detectors in each band was
    ! C    computed. The detector-averaged spectral response data were used
    ! C    to compute the effective central wavenumbers and temperature
    ! C    correction coefficients included in this module.
    ! C
    implicit none

!    c ... Arguments
    real rad
    integer band, units

!   Local variables
    real cwn, tcs, tci
    integer index

!c ... Set default return value
    modis_bright = -1.0

!c ... Check input parameters and return if they are bad
    if (rad .le. 0.0 .or.  &
       band .lt. 20 .or.  &
       band .gt. 36 .or.  &
       band .eq. 26) return

!c ... Get index into coefficient arrays
    if (band .le. 25) then
      index = band - 19
    else
      index = band - 20
    endif

!c ... Get the coefficients for Terra or Aqua
   !  if (platform_name(1:5) .eq. 'Terra' .or.
   ! &    platform_name(1:5) .eq. 'terra' .or.
   ! &    platform_name(1:5) .eq. 'TERRA') then
    cwn = cwn_terra(index)
    tcs = tcs_terra(index)
    tci = tci_terra(index)
   !  else if (platform_name(1:4) .eq. 'Aqua' .or.
   ! &         platform_name(1:4) .eq. 'aqua' .or.
   ! &         platform_name(1:4) .eq. 'AQUA') then
   !    cwn = cwn_aqua(index)
   !    tcs = tcs_aqua(index)
   !    tci = tci_aqua(index)
!    else

!        call message('modis_bright.f',
!    &    'Platform name not recognized ' //
!     &    '[OPERATOR ACTION: Contact SDST]', 0, 2)
!    endif

!c ... Compute brightness temperature
    if (units .eq. 1) then

!c ...   Radiance units are
!c ...   Watts per square meter per steradian per micron
      modis_bright = (bright_m(1.0e+4 / cwn, rad) - tci) / tcs

    else

!c ...   Radiance units are
!c ...   milliWatts per square meter per steradian per wavenumber
      modis_bright = (brite_m(cwn, rad) - tci) / tcs

    endif
    END FUNCTION

  real function bright_m(w, r)

! c-----------------------------------------------------------------------
! c!
! c
! c!DESCRIPTION:
! c    Compute brightness temperature given monochromatic Planck radiance
! c    (Radiance units: Watts per square meter per steradian per micron)
! c
! c!INPUT PARAMETERS:
! c    W (REAL)           Wavelength (microns)
! c    R (REAL)           Monochromatic Planck radiance (Watts per
! c                       square meter per steradian per micron)
! c
! c!OUTPUT PARAMETERS:
! c    BRIGHT_M (REAL)    Brightness temperature (Kelvin)
! c
! c!REVISION HISTORY:
! c
! c!TEAM-UNIQUE HEADER:
! c    Liam.Gumley@ssec.wisc.edu
! c
! c!END
! c-----------------------------------------------------------------------

    implicit none

!c ... Include files
!    include 'fundamental_constants.inc'

!c ... Arguments
    real w, r

!c ... Local variables
    double precision ws

!c ... Set default return value
    bright_m = -1.0

!c ... Check input parameters and return if they are bad
    if (w .le. 0.0 .or. r .le. 0.0) return

!c ... Convert wavelength to meters
    ws = 1.0d-6 * dble(w)

!c ... Compute brightness temperature
    bright_m = sngl(c2 /    &
     (ws * log(c1 / (1.0d+6 * dble(r) * ws**5) + 1.0d+0)))

    END FUNCTION

  real function brite_m(v, r)

! c-----------------------------------------------------------------------
! c!F77
! c
! c!DESCRIPTION:
! c    Compute brightness temperature given monochromatic Planck radiance
! c    (Radiance units: milliWatts per square meter per steradian per
! c    inverse centimeter)
! c
! c!INPUT PARAMETERS:
! c    V (REAL)          Wavenumber (inverse centimeters)
! c    R (REAL)          Monochromatic Planck radiance (milliWatts per
! c                      square meter per steradian per
! c                      inverse centimeter)
! c
! c!OUTPUT PARAMETERS:
! c    BRITE_M (REAL)    Brightness temperature (Kelvin)
! c
! c!REVISION HISTORY:
! c
! c!TEAM-UNIQUE HEADER:
! c    Liam.Gumley@ssec.wisc.edu
! c
! c!END
! c-----------------------------------------------------------------------
    implicit none


! c ... Include files
!     include 'fundamental_constants.inc'

! c ... Arguments
    real v, r

! c ... Local variables
    double precision vs

! c ... Set default return value
    brite_m = -1.0

! c ... Check input parameters and return if they are bad
    if (v .le. 0.0 .or. r .le. 0.0) return

! c ... Convert wavenumber to inverse meters
    vs = 1.0d+2 * dble(v)

! c ... Compute brightness temperature
    brite_m = sngl(c2 *    &
     vs / log(c1 * vs**3 / (1.0d-5 * dble(r)) + 1.0d+0))

   END FUNCTION

   real function modis_planck(temp, band, units)

! C-----------------------------------------------------------------------
! C!F77
! C
! C!DESCRIPTION:
! C    Compute Planck radiance for a MODIS infrared band on Terra or Aqua.
! C
! C    Spectral responses for each IR detector were obtained from MCST:
! C    ftp://ftp.mcst.ssai.biz/pub/permanent/MCST in directories
! C    PFM_L1B_LUT_4-30-99 (Terra) and FM1_RSR_LUT_07-10-01 (Aqua).
! C
! C    An average spectral response for all detectors in each band was
! C    computed. The detector-averaged spectral response data were used
! C    to compute the effective central wavenumbers and temperature
! C    correction coefficients included in this module.
! C
! C    NOTE: The plaform name ('Terra' or 'Aqua') is passed to this
! C    function via the common block defined in 'platform_name.inc'.
! C
! C!INPUT PARAMETERS:
! C    TEMP (REAL)     Brightness temperature (Kelvin)
! C    BAND (LONG)     MODIS IR band number (20-25, 27-36)
! C    UNITS (LONG)    Flag defining radiance units
! C                    0 => milliWatts per square meter per
! C                         steradian per inverse centimeter
! C                    1 => Watts per square meter per
! C                         steradian per micron
! C
! C!OUTPUT PARAMETERS:
! C    MODIS_PLANCK  Planck radiance (units are determined by UNITS)
! C                  Note that a value of -1.0 is returned if
! C                  TEMP .LE. 0.0, or BAND is not in range 20-25, 27-36.
! C
! C!REVISION HISTORY:
! C    Liam.Gumley@ssec.wisc.edu
! C
! C!TEAM-UNIQUE HEADER:
! C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
! C
! C!END
! C-----------------------------------------------------------------------
    implicit none

!c ... Arguments
    real temp
    integer band, units
!c ... Local variables
    real cwn, tcs, tci
    integer index

!c ... Set default return value
    modis_planck = -1.0

! c ... Check input parameters and return if they are bad
    if (temp .le. 0.0 .or.  &
       band .lt. 20 .or.   &
       band .gt. 36 .or.   &
       band .eq. 26) return

! c ... Get index into coefficient arrays
    if (band .le. 25) then
      index = band - 19
    else
      index = band - 20
    endif

! c ... Get the coefficients for Terra or Aqua
!     if (platform_name(1:5) .eq. 'Terra' .or.
!    &    platform_name(1:5) .eq. 'terra' .or.
!    &    platform_name(1:5) .eq. 'TERRA') then
    cwn = cwn_terra(index)
    tcs = tcs_terra(index)
    tci = tci_terra(index)
!     else if (platform_name(1:4) .eq. 'Aqua' .or.
!    &         platform_name(1:4) .eq. 'aqua' .or.
!    &         platform_name(1:4) .eq. 'AQUA') then
!       cwn = cwn_aqua(index)
!       tcs = tcs_aqua(index)
!       tci = tci_aqua(index)
!     else
! CCC
!       print *, ' modis_planck.f: platform_name not recognized.'
! CCC        call message('modis_planck.f',
! CCC     &    'Platform name not recognized ' //
! CCC     &    '[OPERATOR ACTION: Contact SDST]', 0, 2)
!     endif

!c ... Compute Planck radiance
    if (units .eq. 1) then

!c ...   Radiance units are
!c ...   Watts per square meter per steradian per micron
      modis_planck = planck_m(1.0e+4 / cwn, temp * tcs + tci)

    else

!c ...   Radiance units are
!c ...   milliWatts per square meter per steradian per wavenumber
      modis_planck = planc_m(cwn, temp * tcs + tci)

    endif

  END FUNCTION


  real function planc_m(v, t)

! c-----------------------------------------------------------------------
! c!F77
! c
! c!DESCRIPTION:
! c    Compute monochromatic Planck radiance given brightness temperature
! c    (Radiance units: milliWatts per square meter per steradian per
! c    inverse centimeter)
! c
! c!INPUT PARAMETERS:
! c    V (REAL)          Wavenumber (inverse centimeters)
! c    T (REAL)          Brightness temperature (Kelvin)
! c
! c!OUTPUT PARAMETERS:
! c    PLANC_M (REAL)    Monochromatic Planck radiance (milliWatts per
! c                      square meter per steradian per
! c                      inverse centimeter)
! c
! c!REVISION HISTORY:
! c
! c!TEAM-UNIQUE HEADER:
! c    Liam.Gumley@ssec.wisc.edu
! c
! c!END
! c-----------------------------------------------------------------------

    implicit none

! c ... Include files
!     include 'fundamental_constants.inc'

!c ... Arguments
    real v, t

!c ... Local variables
    double precision vs

!c ... Set default return value
    planc_m = -1.0

!c ... Check input parameters and return if they are bad
    if (v .le. 0.0 .or. t .le. 0.0) return

!c ... Convert wavenumber to inverse meters
    vs = 1.0d+2 * dble(v)

!c ... Compute Planck radiance
    planc_m = sngl(1.0d+5 * (c1 * vs**3) /  &
     (exp(c2 * vs / dble(t)) - 1.0d+0))

  END FUNCTION

  real function planck_m(w, t)
! c-----------------------------------------------------------------------
! c!F77
! c
! c!DESCRIPTION:
! c    Compute monochromatic Planck radiance given brightness temperature
! c    (Radiance units: Watts per square meter per steradian per micron)
! c
! c!INPUT PARAMETERS:
! c    W (REAL)           Wavelength (microns)
! c    T (REAL)           Brightness temperature (Kelvin)
! c
! c!OUTPUT PARAMETERS:
! c    PLANCK_M (REAL)    Monochromatic Planck radiance (Watts per
! c                       square meter per steradian per micron)
! c
! c!REVISION HISTORY:
! c
! c!TEAM-UNIQUE HEADER:
! c    Liam.Gumley@ssec.wisc.edu
! c
! c!END
! c-----------------------------------------------------------------------

    implicit none

! c ... Include files
    ! include 'fundamental_constants.inc'

! c ... Arguments
    real w, t

! c ... Local variables
    double precision ws

! c ... Set default return value
    planck_m = -1.0

! c ... Check input parameters and return if they are bad
    if (w .le. 0.0 .or. t .le. 0.0) return

! c ... Convert wavelength to meters
    ws = 1.0d-6 * dble(w)

! c ... Compute Planck radiance
    planck_m = sngl(1.0d-6 * (c1 / ws**5) /  &
     (exp(c2 / (ws * dble(t))) - 1.0d+0))

  END FUNCTION

end module planck_bright
