module emission_rte
!  use planck_bright
  implicit none
contains
  real function fm_modrad_emis( tskn, esfc, psfc, pres, &
                              temp, tau, band, lsfc )

! c-----------------------------------------------------------------------
! c
! c!F77
! c
! c!Description:
! c    Compute MODIS IR band radiance at top of atmosphere.
! c
! c!Input Parameters:
! c    TSKN      Surface skin temperature (K)
! c    TEMP      Temperature profile (K) at NL levels
! c    TAU       Transmittance profile (no units) at NL levels
! c    BAND      MODIS IR band number (20-25,27-36)
! c    NL        Number of levels in TEMP and TAU profiles
! c
! c SWS aded: PRES Pressure profile (hPa)
! c           ESFC surface emissivity
! c           PSFC surface pressure
! c   changed NL to LSFC level of surface in pressure array
! c
! c!Output Parameters:
! c    MODRAD    Radiance at top of atmosphere
! c              (milliWatts per square meter per steradian per
! c              inverse centimeter)
! c
! c!Revision History:
! c SWS 10 October 2003: modified fm_modrad by accounting for surface
! c   emissivity and increment between surface pressure and closest
! c   level (based on tbbe_modis.f)
! c
! c!Team-unique Header:
! c
! c!End
! c
! c-----------------------------------------------------------------------
  use planck_bright
  implicit none

!c ... Arguments

  real tskn, esfc, psfc, temp(*), tau(*), pres(*)
  integer band, lsfc

!c ... Local variables

  real b1, b2, rad, t1, t2, tau1, tau2
  real dpre, dtau, taus, taub, drad, refl, taur, bs
  integer i

!c ... Parameters

  integer units
  parameter ( units = 0 )

!c ... External function

!  real modis_planck_shift, modis_planck
!  external modis_planck_shift, modis_planck

!c ... Initialize

  tau1 = tau( 1 )
  t1 = temp(  1)
!  b1 = modis_planck_shift( t1, band, units )
  b1 = modis_planck( t1, band, units )
  rad = 0.0
!c SWS added refl
  refl = 0.0

!c ... SWS added: compute differences
  dpre = alog(psfc/pres(lsfc))/alog(pres(lsfc)/pres(lsfc-1))
  dtau = tau(lsfc-1)-tau(lsfc)
  taus = tau(lsfc)-dtau*dpre

!c ... Compute RTE

  do i = 2, lsfc
    tau2 = tau( i )
    t2 = temp( i )
!    b2 = modis_planck_shift( t2, band, units )
     b2 = modis_planck( t2, band, units )
    drad = 0.5 * (b1 + b2) * (tau1 - tau2)
    rad = rad + drad

    if(taus.gt.0.1.and.esfc.lt.1.00) then
!c * Do not add reflected component for last level unless psfc .gt. 1000.
       if((i.ne.lsfc).or.(i.eq.lsfc.and.psfc.gt.1000.)) then
          taur=taus/(0.5*(tau1+tau2))
          refl=refl+taur*taur*drad
       endif
    endif

    tau1 = tau2
    b1 = b2
  end do

!c * Add (subtract) increment of atmospheric radiance to reach surface.
!c    dpre will be negative if psfc < 1000 mb
!c    drad falls out as the delta radiance of layer

  rad=rad+drad*dpre

!c * Add increment of reflected radiance for layer down to surface.
  if(taus.gt.0.1 .and. esfc.lt.1.00) then
     if(psfc.lt.1000.) then
        taub=0.5*(tau(lsfc-1)+taus)
        dpre=1.0+dpre
     else
        taub=0.5*(tau(lsfc)+taus)
     endif
     taur=taus/taub
     refl=refl+taur*taur*drad*dpre
  endif

  rad=rad+(1.-esfc)*refl
!  bs = modis_planck_shift(tskn, band, units)
  bs = modis_planck(tskn, band, units)
  rad=rad+esfc*bs*taus
  rad=amax1(rad,.001)

!c ... Set return value

  fm_modrad_emis = rad

  end function

  real function fm_modrad_emis_twolayer( tskn, esfc, psfc, pres, &
                              temp, tau_in, band, lsfc, &
                              upper_cld_ctp, upper_emissivity)

! c-----------------------------------------------------------------------
! c
! c!F77
! c
! c!Description:
! c    Compute MODIS IR band radiance at top of atmosphere.
! c
! c!Input Parameters:
! c    TSKN      Surface skin temperature (K)
! c    TEMP      Temperature profile (K) at NL levels
! c    TAU       Transmittance profile (no units) at NL levels
! c    BAND      MODIS IR band number (20-25,27-36)
! c    NL        Number of levels in TEMP and TAU profiles
! c
! c SWS aded: PRES Pressure profile (hPa)
! c           ESFC surface emissivity
! c           PSFC surface pressure
! c   changed NL to LSFC level of surface in pressure array
! c
! c!Output Parameters:
! c    MODRAD    Radiance at top of atmosphere
! c              (milliWatts per square meter per steradian per
! c              inverse centimeter)
! c
! c!Revision History:
! c SWS 10 October 2003: modified fm_modrad by accounting for surface
! c   emissivity and increment between surface pressure and closest
! c   level (based on tbbe_modis.f)
! c
! c!Team-unique Header:
! c
! c!End
! c
! c-----------------------------------------------------------------------
  use planck_bright
  implicit none

!c ... Arguments

  real tskn, esfc, psfc, temp(*), tau_in(*), pres(*), upper_emissivity
  real upper_cld_ctp
  integer band, lsfc

!c ... Local variables
  real b1, b2, rad, t1, t2, tau1, tau2, tau(lsfc)
  real dpre, dtau, taus, taub, drad, refl, taur, bs
  integer i, upper_cld_i

!c ... Parameters

  integer units
  parameter ( units = 0 )

!c ... External function

!  real modis_planck_shift, modis_planck
!  external modis_planck_shift, modis_planck

!c ... Initialize
!find upper_cld_ctp
  do i=1,lsfc
    if (pres(i) .le. upper_cld_ctp) then
      upper_cld_i = i
    end if
  enddo
  
!update transmission profile to include gray cloud at upper_cld_i.
  do i=1,lsfc
    if (i .lt. upper_cld_i) then
      tau(i) = tau_in(i)
    else !modify all taus below the upper cloud layer.
      tau(i) = tau_in(i)*(1.0 - upper_emissivity)
    endif
  enddo

  tau1 = tau( 1 )
  t1 = temp(  1)
!  b1 = modis_planck_shift( t1, band, units )
  b1 = modis_planck( t1, band, units )
  rad = 0.0
!c SWS added refl
  refl = 0.0

!c ... SWS added: compute differences
  dpre = alog(psfc/pres(lsfc))/alog(pres(lsfc)/pres(lsfc-1))
  dtau = tau(lsfc-1)-tau(lsfc)
  taus = tau(lsfc)-dtau*dpre

!c ... Compute RTE

  do i = 2, lsfc
    tau2 = tau( i )
    t2 = temp( i )
!    b2 = modis_planck_shift( t2, band, units )
     b2 = modis_planck( t2, band, units )
    drad = 0.5 * (b1 + b2) * (tau1 - tau2)
    rad = rad + drad

!     if(taus.gt.0.1.and.esfc.lt.1.00) then
! !c * Do not add reflected component for last level unless psfc .gt. 1000.
!        if((i.ne.lsfc).or.(i.eq.lsfc.and.psfc.gt.1000.)) then
!           taur=taus/(0.5*(tau1+tau2))
!           refl=refl+taur*taur*drad
!        endif
!     endif

    tau1 = tau2
    b1 = b2
  end do

!c * Add (subtract) increment of atmospheric radiance to reach surface.
!c    dpre will be negative if psfc < 1000 mb
!c    drad falls out as the delta radiance of layer

  rad=rad+drad*dpre

! lower cloud is always black so this is not needed.
!c * Add increment of reflected radiance for layer down to surface.
!   if(taus.gt.0.1 .and. esfc.lt.1.00) then
!      if(psfc.lt.1000.) then
!         taub=0.5*(tau(lsfc-1)+taus)
!         dpre=1.0+dpre
!      else
!         taub=0.5*(tau(lsfc)+taus)
!      endif
!      taur=taus/taub
!      refl=refl+taur*taur*drad*dpre
!   endif
!
!   rad=rad+(1.-esfc)*refl
!  bs = modis_planck_shift(tskn, band, units)
  bs = modis_planck(tskn, band, units)
  rad=rad+esfc*bs*taus
  rad=amax1(rad,.001)

!c ... Set return value

  fm_modrad_emis_twolayer = rad

  end function


end module emission_rte
