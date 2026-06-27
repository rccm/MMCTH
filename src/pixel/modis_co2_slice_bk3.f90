!!!!This is the model before working on debugging why no retrievals happen. 
module co2
  integer nb_wavelen
  parameter (nb_wavelen = 7)
  integer plev
  parameter (plev = 101)
  integer max_solutns
  parameter (max_solutns = 5 )
  integer tau_lev
  parameter (tau_lev = 101)
  integer nsct
  parameter (nsct = 4)
  integer nbct
  parameter (nbct = 5)
  integer ntbct
  parameter (ntbct = 7)

! Fill Value for all cases
  real badvalue
  parameter(badvalue = -9999.0)

  real pp(plev)
  integer kch(nsct, 2)
  integer mbnds(nb_wavelen)
  real freq(nb_wavelen)
  real rmin(nbct)
  integer std_doy(12), leap_doy(12)
  integer jdet
  real emisrw, emadj(nsct)

! LEGACY: Ice cloud emissivity adjustments for each channel combination.
  data emadj / 1.000, 1.000, 1.000, 1.000 /

  data jdet /0/

  data emisrw /1.0/

  data std_doy  /0,31,59,90,120,151,181,212,243,273,304,334/
  data leap_doy /0,31,60,91,121,152,182,213,244,274,305,335/

! LEGACY: Emissivity adjustment using window channel removed (29, 32).
  data mbnds/36,35,34,33,31,29,32/
  data kch /1,2,2,3,2,3,4,4/

! Define wavenumbers of CO2 bands (36, 35, 34, 33, 31), Legacy 29,32.
  data freq /7.055309E+02, 7.196677E+02, 7.317089E+02,&
                  7.483224E+02, 9.081998E+02, 1.173198E+03,&
                  8.315149E+02/

! Min radiance difference to not be counted as channel noise.
  data rmin /-1.0, -1.0, -100.0, -1.0, -0.5/

! 101 MODIS pressure levels for ancillary inputs and CO2-slicing CTP solutions.
  data pp   / 0.0050,    0.0161,    0.0384,    0.0769,    0.1370,&
     0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972,&
     1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204,&
     5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492,&
    14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829,&
    29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882,&
    51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396,&
    83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775,&
   125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784,&
   180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338,&
   247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369,&
   328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738,&
   424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200,&
   535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398,&
   661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897,&
   802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236,&
   958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000/


end module co2

subroutine co2cld_onepixel_misr(wprof, tprof, hprof, psfc, pmsl, surftmp, view, &
                                trad, met_date, rlat, rlon, landsea, misr_ctp, &
                                Cloud_Top_Pressure, &
                                Cloud_Top_Height, &
                                Cloud_Effective_Emissivity, &
                                Cloud_Optical_Depth, &
                                Processing_Mask)
!-------------------------------------------------------------------------------
! Name: co2cld_onepixel_misr
!
! Purpose:
!   Retrieves upper cloud-top pressure (CTP), height (CTH), temperature, and
!   effective emissivity for a two-layer cloud system using CO2 slicing, with
!   MISR CTP as the lower cloud pressure (MODIS MOD06 ATBD, section 3.3.4,
!   equations 11-15). Processes a single pixel with valid MISR CTP.
!
! Inputs:
!   wprof             : Water vapor mixing ratio profile (g/kg), 101 levels: Top-Down
!   tprof             : Temperature profile (K), 101 levels:  Top-Down 
!   psfc              : Surface pressure (hPa)
!   pmsl              : Mean sea level pressure (hPa)
!   surftmp           : Surface temperature (K)
!   view              : Satellite zenith angle (degrees)
!   tcold             : Brightness temperatures (K) for MODIS bands 36,35,34,33,31
!   met_date          : Array of year, month, day, hour
!   rlat, rlon        : Latitude, longitude (degrees)
!   landsea           : Land/sea flag (0=water, 1=land, 2=coast)
!   misr_ctp          : MISR lower cloud CTP (hPa)
!
! Outputs:
!   Cloud_Top_Pressure        : Upper cloud CTP (hPa)
!   Cloud_Top_Temperature     : Upper cloud temperature (K)
!   Cloud_Effective_Emissivity: Upper cloud effective emissivity (dimensionless)
!   Cloud_Optical_Depth       : Upper cloud optial depth (dimensionless)                        
!   Processing_Mask           : Quality flag (2=success, 1=emissivity invalid, 0=failed)
!
! Notes:
!   - Uses MISR CTP as lower cloud pressure (P_l), searching for upper cloud (P_u).
!   - Validates inputs to prevent crashes or incorrect retrievals.
!   - Eliminates 12, 8.6 μm radiance profiles, as MISR CTP defines lower boundary.
!   - 11 μm radiance is only used in emissivity retrievals directly.
!   - Assumes clozo101 and tran_modisd101 are functional for ozone and transmittance.
!   - Variable names align with MODIS PGE (MOD06) conventions.
!   - Debug messages aid troubleshooting (e.g., invalid inputs, retrieval steps).
!   - Updated the logic of the Processing Flag (Nov 18, 2025):

!   - Processing Flag Logic: For each call to co2cld_onepixel_misr:
!
! 000
! no CTP, no emissivity
! 010
! 2-layer CTP from 36/35, no emissivity
! 012
! 2-layer CTP from 36/35, emissivity from 2-layer call
! 020
! 2-layer CTP from 35/33, no emissivity
! 022
! 2-layer CTP from 35/33, emissivity from 2-layer call
! 100
! 1-layer CTP from 36/35, no emissivity
! 101
! 1-layer CTP from 36/35, emissivity from 1-layer call
! 200
! 1-layer CTP from 35/33, no emissivity
! 201
! 1-layer CTP from 35/33, emissivity from 1-layer call

!    1st digit (A) — 1-layer CTP (this call only if misr_ctp ≈ psfc)
!
!    0: no 1-layer CTP (or this call is 2-layer)
!    1: band 36/35
!    2: band 35/33
!
!    2nd digit (B) — 2-layer CTP (this call only if misr_ctp < psfc)
!
!    0: no 2-layer CTP (or this call is 1-layer)
!    1: band 36/35
!    2: band 35/33
!
!    3rd digit (C) — emissivity source in this call
!
!    0: emissivity = –9999 (no retrieval)
!    1: emissivity from a 1-layer CTP solution (this call has misr_ctp = psfc)
!    2: emissivity from a 2-layer CTP solution (this call has misr_ctp < psfc)
!
! Author: Arka Mita (rewritten by Guangyu Zhao; edited further by Arka Mitra)
! Date: April 2x, 2025; July 2025; Oct 2025 -...
!-------------------------------------------------------------------------------

  use co2, only: nbct, plev, pp, nb_wavelen, ntbct, mbnds, jdet, freq, &
                 kch, nsct, rmin, emadj, badvalue, emisrw
  use, intrinsic :: ieee_arithmetic
  use transmission, only: tran_modisd101, wstd,tstd
  use surfemis, only: getiremis, assign_eco_emis
  use emission_rte, only: fm_modrad_emis
  use planck_bright, only: modis_bright, modis_planck, modis_planck_shift, &
                           modis_bright_shift
  implicit none
! Inputs
  real, intent(in) :: wprof(plev)               ! Water vapor mixing ratio (g/kg)
  real, intent(in) :: tprof(plev)               ! Temperature profile (K)
  real, intent(in) :: hprof(plev)               ! Height profile (km)
  real, intent(in) :: psfc                      ! Surface pressure (hPa)
  real, intent(in) :: pmsl                      ! Mean sea level pressure (hPa)
  real, intent(in) :: surftmp                   ! Surface temperature (K)
  real, intent(in) :: view                      ! Satellite zenith angle (degrees)
  real, intent(in) :: trad(nbct)               ! Brightness temperatures (K)
  integer, intent(in) :: met_date(4)            ! Year, month, day, hour
  real, intent(in) :: rlat, rlon                ! Latitude, longitude
  integer, intent(in) :: landsea                ! Land/sea flag
  real, intent(in) :: misr_ctp                  ! MISR lower cloud CTP (hPa)

! Outputs
  real, intent(out) :: Cloud_Top_Pressure       ! Upper cloud CTP (hPa)
  real, intent(out) :: Cloud_Top_Height       ! Upper cloud CTH (km)
  real, intent(out) :: Cloud_Effective_Emissivity ! Upper cloud effective emissivity
  real, intent(out) :: Cloud_Optical_Depth
  integer, intent(out) :: Processing_Mask       ! Quality flag

! Local variables
  real :: tcold(nbct)  
  integer :: ltrp, isp, imisr, imslp, ll, k, id, lmin, lco2,jday
  integer ::lwin, iw1 
  integer :: date2doy   
  integer :: nl                     ! Level where water/temperature profiles = 0
  integer ::krto(nsct)
  real :: t11(plev)
  real :: tmin, ptrp,rwarm,Iobs
  real :: rclr(nbct), robs(nbct), delr(nbct), ra(plev,nbct)
  real :: sum, db, fm, fm1, fm2
  integer :: is1, k1, k2
  real :: rwcld, bot, top, ratio
  real :: z(plev), zs                           ! Geopotential height (km)
  integer :: ctp_flag(nsct)                     ! Per-band-pair flag
  real :: ctp_pres(nsct)                        ! Temporary CTP storage
  real :: sfc_emis(ntbct)                       ! Surface emissivity
  real :: ozprof(plev)                          ! Ozone profile (ppmv)
  real :: taup(plev,ntbct)                      ! Optical depth
  real :: freqemis(nb_wavelen), emis(nb_wavelen), rho(nb_wavelen)
  real :: emis_out, rho_out
  real :: lev(nsct),fmsav, tmisr
  integer :: met_year, met_month, met_day
  real :: tpad(plev), wpad(plev), ozpad(plev)
  integer            :: newunit, CO2_Slicing_Flag, kwc
  character(len=40)  :: fname
  logical ::ok, neg, start
  real ::  den, num, deltaI
  integer :: ipco2              ! Index of selected band pair
  real :: pfco2, ecaco2         ! Final CTP and emissivity
  real :: ecawin                ! Window-based emissivity (assumed 1.0 for single pixel)
  logical :: found_solution     ! Flag for valid CTP
  real, parameter :: xi_vis = 2.13          ! Minnis et al. 1990, ice clouds
  real :: tau_vis 
  integer :: digit1, digit2, digit3 ! three digits of the Processing Flag
  logical :: is_one_layer_call ! Logical flag to figure out if MISR or SFC Pressure used

! Initialize Processing Flag Pipeline
  digit1 = 0
  digit2 = 0
  digit3 = 0

  Cloud_Top_Height = -9999.0

  is_one_layer_call = (abs(misr_ctp - psfc) <= 1.0)

  found_solution   = .false.
  Processing_Mask  = 0


! Initialize outputs
  Cloud_Top_Pressure = badvalue
  Cloud_Effective_Emissivity = badvalue
  Cloud_Optical_Depth = badvalue
  CO2_Slicing_Flag = 0
!  Processing_Mask = 3
  ctp_flag = 0
  ctp_pres = badvalue

  do k = 1, nbct
     tcold(k) = modis_bright_shift(trad(k),mbnds(k), 1)
  end do 
!-------------------------------------------------------------------------------
! Validate inputs
!-------------------------------------------------------------------------------
  if (misr_ctp <= 0.0 .or. misr_ctp > psfc .or. misr_ctp < pp(1)) then
    ! write(*,'(A,F10.2,F10.2,F10.2,A)') 'ERROR: Invalid misr_ctp = ', misr_ctp, psfc, pp(1), &
    !                        ' hPa, returning badvalue'
    return
  endif
  if (psfc <= pp(1) .or. psfc > 1100.0) then
    ! write(*,'(A,F10.2,A)') 'ERROR: Invalid psfc = ', psfc, ' hPa, returning badvalue'
    return
  endif
  if (any(tcold <= 0.0)) then
    ! write(*,'(A)') 'ERROR: Invalid tcold values, returning badvalue'
    return
  endif

! Extract meteorological date
  met_year = met_date(1)
  met_month = met_date(2)
  met_day = met_date(3)

!-------------------------------------------------------------------------------
! Find pressure levels
!-------------------------------------------------------------------------------
! Surface level
  isp = 0
  do ll = 1, plev
    if (psfc <= pp(ll)) then
      isp = ll
      exit
    end if
  enddo
  if (isp == 0) then
    ! write(*,'(A,F10.2,A)') 'ERROR: Cannot find isp for psfc = ', psfc, &
    !                        ' hPa, returning badvalue'
    return
  endif
! MISR CTP level
  imisr = 0
  do ll = 1, plev
    if (misr_ctp <= pp(ll)) then
      imisr = ll - 1
      exit
    end if
  enddo
  if (imisr == 0) then
    ! write(*,'(A,F10.2,A)') 'ERROR: Cannot find imisr for misr_ctp = ', misr_ctp, &
    !                        ' hPa, returning badvalue'
    return
  endif
  is1  = imisr - 1

! Mean sea level pressure level (for height calculation): is this nessary?
  imslp = 0
  do ll = 1, plev
    if (pmsl <= pp(ll)) then
      imslp = ll
      exit
    end if
  enddo
  if (imslp == 0) then
    imslp = isp
    ! write(*,'(A,F10.2,A)') 'WARNING: Using isp for imslp, pmsl = ', pmsl, ' hPa'
  endif
! ! First level with water/temperature profiles = 0
!   nl = 0
!   do ll = 1, plev
!     if (wprof(ll) < 0 .or. tprof(ll) <0 ) then
!       nl = ll
!       exit
!     end if
!   enddo
! Tropopause level (coldest temperature or inflection point)
  tmin = 350.0
  lmin = isp
  do ll = 44, isp
    if (tprof(ll) <= 0.0) cycle
    if (tprof(ll) <= tmin) then
      tmin = tprof(ll)
      lmin = ll
    end if
  enddo
  ltrp = -99
  do ll = lmin, min(71, isp-1)
    if (tprof(ll-1) >= tprof(ll) .and. tprof(ll+1) > tprof(ll)) then
      ltrp = ll
      exit
    end if
  enddo
  if (ltrp == -99) ltrp = lmin
  if (ltrp < 2) ltrp = 2  
  ptrp = pp(ltrp)
  if (ltrp >= imisr) then
    ! write(*,'(A,I3,A,F10.2,A)') 'ERROR: Tropopause level ', ltrp, &
    !                             ' at ', ptrp, ' hPa below misr_ctp, returning badvalue'
    return
  endif
!  write(*,*) imisr, isp, pp(imisr), pp(isp)
!-------------------------------------------------------------------------------
! Compute auxiliary profiles
!-------------------------------------------------------------------------------
! Geopotential height
  zs = 0.0
  z = 0.0
  call height(pp, tprof, wprof, zs, imslp, z)

! Ozone profile
  call clozo101(rlat, met_month, ozprof)
  if (all(ozprof == 0.0)) then
    write(*,'(A)') 'WARNING: clozo101 returned zero ozone profile'
  endif

! Surface emissivity
  call assign_eco_emis(landsea, emis, rho, freqemis)
  do k = 1, ntbct
    call getiremis(nb_wavelen, freqemis, freq(k), emis, rho, emis_out, rho_out)
    sfc_emis(k) = emis_out
  enddo

  tpad   = tprof            ! copy original top → surface
  wpad   = wprof
  ozpad  = ozprof
  tmisr  = tpad(imisr-1)      ! storing CTT for later usage 
  !if (tmisr < 0) then
  !  tmisr = surftmp
  !endif

  do ll = imisr-1, plev       ! layers that are physically below the terrain
    tpad(ll)  = tmisr        ! hold at surface T
    wpad(ll)  = 0.0        ! zero vapour → negligible extra τ
    ozpad(ll) = 0.0
  end do 
   
  jday = date2doy( met_date(1), met_date(2), met_date(3) )

   
! Transmittance profiles
  do k = 1, nbct
    call tran_modisd101(met_year, jday, tpad, wpad, ozpad, view, mbnds(k), jdet, taup(:,k), plev )               ! nl = 101
    !call tran_modisd101(met_year, jday, tstd, wstd, ozpad, view, mbnds(k), jdet, taup(:,k), plev )               ! nl = 101
    if (taup(1,k) <= 0.0 .or. taup(1,k) > 1.0) then
        write(*,*) 'ERORR: Transmittance for ',  mbnds(k), ': ', taup(1,k)
    endif
  enddo

  ! 2.  Clear / cloudy radiances and reference terms  -----------------------------
  do k = 1, nbct - 1
     robs(k) = modis_planck_shift( tcold(k), mbnds(k), 0 )
  
     rclr(k) = fm_modrad_emis( tmisr, 1.0, pp(imisr), pp, tpad, &
                               taup(1,k), mbnds(k), imisr )
  
     delr(k) = robs(k) - rclr(k)
  
     ! --- accumulate RHS integral from Pl up to TOA -----------------------------
     sum = 0.0
     do ll = imisr-2, 1, -1
        db  = modis_planck_shift( tpad(ll+1), mbnds(k), 0 ) - &
              modis_planck_shift( tpad(ll  ), mbnds(k), 0 )
        sum = sum + 0.5 * ( taup(ll+1,k) + taup(ll,k) ) * db
        ra(ll,k) = sum
     end do
  end do
  
  !write(*,'("Pressure levels: trop=",i3,"  Pl=",i3,"  ptrop=",f7.2,"  Pmisr=",f7.2)') &
  !        ltrp, imisr, pp(ltrp), pp(imisr)
  
  ! 3.  Zero‑crossing search for each CO₂‑band pair  ------------------------------
  lev  = badvalue
  krto = 0
  
  pair_loop: do id = 1, nsct
     k1 = kch(id,1)
     k2 = kch(id,2)
  
     ! Skip pair if clear–cloud signals are below noise threshold
     if (delr(k1) > rmin(k1)+101 .or. delr(k2) > rmin(k2)+101) cycle
  
    !  write(*,'("Working on Pair:",i2)') id
     ok    = .false.
     start = .false.
  
     do ll = ltrp, imisr              ! search downward until the MISR layer
        if (delr(k2) == 0.0  .or. ra(ll,k2) == 0.0) cycle     ! avoid /0
  
        fm1 = delr(k1)              / delr(k2)
        fm2 = ra(ll,k1) * emadj(id) / ra(ll,k2)
  
        if (ieee_is_nan(fm1) .or. ieee_is_nan(fm2)) cycle
  
        fm  = fm1 - fm2
        neg = fm < 0.0
  
        if (ll == ltrp) start = neg
        if (neg .neqv. start) then      ! sign change ⇒ zero‑crossing bracket
           ok      = .true.
           lev(id) = ll - 1
           exit
        end if
     end do
  
     if (ok) then
        krto(id) = 1
        lev(id)  = max( 1, min( plev, int(lev(id)) ) )
        ctp_pres(id) = pp(lev(id)-1)
     else
        krto(id)     = 0
        lev(id)      = badvalue
        ctp_pres(id) = badvalue
     end if
  end do pair_loop
  
!-------------------------------------------------------------------------------
! Select best CTP (top-down: 36/35, 35/34, 35/33, 33/33)
!-------------------------------------------------------------------------------
  do id = 1, nsct
    if (krto(id) == 1 .and. ctp_pres(id) /= badvalue) then
      if (id == 1 .and. ctp_pres(id) < 450.0 .and. ctp_pres(id) < misr_ctp) then
        ! Band pair 36/35
        ipco2 = id
        lco2 = nint(lev(id)-1)
        pfco2 = ctp_pres(id)
        found_solution = .true.
        Cloud_Top_Pressure = pp(lco2) 
        Cloud_Top_Height = hprof(lco2)
!        Processing_Mask = 5
      end if
      if (id == 3 .and. ctp_pres(id) < 650.0 .and. ctp_pres(id) < ctp_pres(1)) then
        ! Band pair 35/33
        ipco2 = id
        lco2  = nint(lev(id)-1)
        pfco2 = ctp_pres(id)
        found_solution = .true.
        Cloud_Top_Pressure = pp(lco2)
        Cloud_Top_Height = hprof(lco2)
!        Processing_Mask = 5
      end if
    end if
  end do

! Update the Processing Flag for CTP retrieval
  if (found_solution) then
     if (ipco2 == 1) then          ! Band pair 36/35
        if (is_one_layer_call) then
           digit1 = 1              ! 1st digit (1-layer CTP)
        else
           digit2 = 1              ! 2nd digit (2-layer CTP)
        end if
     else if (ipco2 == 3) then     ! Band pair 35/33
        if (is_one_layer_call) then
           digit1 = 2
        else
           digit2 = 2
        end if
     end if
  end if

! If no valid solution found, check fallback

!===============================================================================
!  EFFECTIVE CLOUD EMISSIVITY  (11 µm window)  –  UNIT‑CONSISTENT VERSION
!==============================================================================

  if (found_solution) then
    Cloud_Top_Pressure  = pfco2 
    kwc = 5                     ! 11‑µm window band (band‑31)
    ll  = lco2                  ! level index of selected upper cloud
    
    robs(5) = modis_planck( tcold(kwc), mbnds(kwc), 0 )
    rclr(5) = fm_modrad_emis( tmisr, 1.0, pp(imisr), pp, tpad, &
                               taup(1,kwc), mbnds(kwc), imisr )
  
    delr(5) = robs(5) - rclr(5)
    rwcld = fm_modrad_emis( tpad(lco2), emisrw, pp(lco2),  pp, tpad, taup(1,kwc), &
                            mbnds(kwc), lco2 )     ! milliWatts per square meter per steradian per cm-1
 
    ! -- 4. Numerator & denominator of Eq. 14 ----------------------------------
    num = delr(5)          ! should be ≤ 0
    den = rwcld - rclr(5)          ! should be ≤ 
    !write(*,*) id, pp(lev(1)), pp(lev(3)), pp(lco2)
    ! -- 6. Emissivity solution -------------------------------------------------
    if (abs(den) > 0.1) then
       ratio = num / den                       ! εA  (0…1)
       if (ratio >= 0.01 .and. ratio <= 1.0) then
          Cloud_Effective_Emissivity = ratio * 1
          CO2_Slicing_Flag          = 1
       else
          Cloud_Effective_Emissivity = badvalue
           Cloud_Optical_Depth = badvalue
       end if
    else
       Cloud_Effective_Emissivity = badvalue
       Cloud_Optical_Depth = badvalue
    end if
    
! Update the Processing Flag for cloud emissivity retrieval    
  if (Cloud_Effective_Emissivity == badvalue) then
     digit3 = 0
  else
     if (is_one_layer_call) then
        digit3 = 1      ! Emissivity from 1-layer CTP
     else
        digit3 = 2      ! Emissivity from 2-layer CTP
     end if
  end if

 
 !===============================================================================  
  !  εA is ratio (window emissivity) * ecawin  (ecawin = 1.0 here)
    if (Cloud_Effective_Emissivity >= 0.9999) then
      tau_vis = 99.9            ! saturation flag for optically thick cloud
    else if (Cloud_Effective_Emissivity <= 0.0) then
      tau_vis = badvalue             ! fully transparent (should not occur)
    else
      tau_vis = -log( 1.0 - Cloud_Effective_Emissivity ) * xi_vis
    end if
    Cloud_Optical_Depth = tau_vis
end if

! Finalize Processing Mask (read logic at the top)
 Processing_Mask = 100*digit1 + 10*digit2 + digit3


end subroutine co2cld_onepixel_misr

subroutine get_clearsky_bias(wprof, tprof, psfc, surftmp, view, rlat, landsea,  &
                             met_date, clr_obs, npix,                            &
                             bias_w_um, n_used)
!-------------------------------------------------------------------------------
! Mean clear-sky bias per channel, meant as a diagnostic test to ascertain the impact
! of not implementing a bias-correction for the modeled clear-sky radiances vis-a-vis
! the observed clear-sky radiances for a given pixel. In the operational MODIS algo,
! the bias-correction methodology employs a 8-day running-mean bias for 1-degree lat-long
! bins. Instead here, the mean speactral biases are calculated over cs-pixels for each 
! granule. We note that for the initial test granules, this approach shows promise and 
! saves us from having to calculate long-term running-means over multiple satellite 
! overpasses over a given location. To enfore speed, we go with max 20k cs pixels/granule.
! The max 20k cs pixels/granule condition can be removed in operational deployment.
!
! Inputs:
!   wprof(plev,npix)   : water vapor [g/kg], top->down
!   tprof(plev,npix)   : temperature [K],    top->down
!   psfc(npix)         : surface pressure [hPa]
!   surftmp(npix)      : skin temperature [K]
!   view(npix)         : satellite zenith [deg]
!   rlat(npix)         : latitude [deg] (for ozone clim)
!   landsea(npix)      : 0 water, 1 land, 2 coast (for emissivity)
!   met_date(4)        : [year, month, day, hour]
!   clr_obs(nbct,npix) : observed clear-sky radiance [W m^-2 sr^-1 µm^-1]
!   npix               : number of pixels
!
! Outputs:
!   bias_w_um(nbct)    : mean(clr_mod - clr_obs) per channel [W m^-2 sr^-1 µm^-1]
!   n_used(nbct)       : number of valid pixels used per channel
!
! Procedure (per pixel/channel), exactly like your Python:
!   - taup = tran_modisd101(...)
!   - rad_wn = fm_modrad_emis(..., taup(1), ...)       ! wavenumber units
!   - tb_mod = modis_bright(rad_wn, band, 0)           ! -> brightness temp
!   - rad_w_um = modis_planck(tb_mod, band, 1)         ! -> W m^-2 sr^-1 µm^-1
!   - accumulate (rad_w_um - clr_obs)
!-------------------------------------------------------------------------------
  use co2,          only: plev, pp, nbct, mbnds, freq, ntbct
  use transmission, only: tran_modisd101
  use surfemis,     only: assign_eco_emis, getiremis
  use emission_rte, only: fm_modrad_emis
  use planck_bright,only: modis_bright, modis_planck
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Inputs
  integer, intent(in) :: met_date(4), npix
  real,    intent(in) :: wprof(plev, npix), tprof(plev, npix)
  real,    intent(in) :: psfc(npix), surftmp(npix), view(npix), rlat(npix)
  integer, intent(in) :: landsea(npix)
  real,    intent(in) :: clr_obs(nbct, npix)

  ! Outputs
  real,    intent(out) :: bias_w_um(nbct)
  integer, intent(out) :: n_used(nbct)

  ! Locals
  integer :: pix, kban, band, jday, mon, yr, lsfc, k, n_pix
  real    :: tbuf(plev), wbuf(plev), oz(plev), taup(plev)
  real    :: freqemis(7), emis_ir(7), rho_ir(7)
  real    :: emiss_ch, rho_ch
  real    :: rad_wn, tb_mod, rad_w_um, rad_obs
  real    :: sum_bias(nbct)
  real    :: dmin, dcur
  integer, external :: date2doy

  sum_bias = 0.0
  n_used   = 0

  yr  = met_date(1)
  mon = met_date(2)
  jday = date2doy(met_date(1), met_date(2), met_date(3))
  if (npix > 20000) then
      n_pix = 20000
  else
      n_pix = npix
  end if

  do pix = 1, npix

    ! copy profiles (top -> down)
    tbuf(:) = tprof(:,pix)
    wbuf(:) = wprof(:,pix)

    ! ozone clim for this pixel/month
    call clozo101(rlat(pix), mon, oz)

    ! emissivity lookup tables for this pixel
    call assign_eco_emis(landsea(pix), emis_ir, rho_ir, freqemis)

    ! find surface index for RTE
    dmin = 1.0e30
    lsfc = 1
    do k = 1, plev
      dcur = abs(pp(k) - psfc(pix))
      if (dcur < dmin) then
        dmin = dcur
        lsfc = k
      end if
    end do

    ! loop channels (expected: 36,35,34,33,31 with nbct=5)
    ! skip correcting the 11 um channel here; that's another story
    ! look to official MODIS code for guidance
    do kban = 1, nbct - 1
       band = mbnds(kban)

       ! channel emissivity at nominal center wavenumber
       call getiremis(nemis=7, freqemis=freqemis, freq=freq(kban),  &
                      emisir=emis_ir, rhoir=rho_ir,                 &
                      emiss=emiss_ch, rho=rho_ch)

       ! transmittance (top->down). Positional args to match your code.
       call tran_modisd101( yr, jday, tbuf, wbuf, oz, view(pix), band, 0, taup, plev )

       ! forward-model clear-sky radiance at TOA (wavenumber units)
       rad_wn = fm_modrad_emis( surftmp(pix), emiss_ch, psfc(pix), pp, tbuf, taup(1), band, lsfc )
       if (ieee_is_nan(rad_wn) .or. rad_wn <= 0.0) cycle

       ! -> Tb (wavenumber units) then back to W m^-2 sr^-1 µm^-1
       tb_mod = modis_bright(rad_wn, band, 0)
       if (ieee_is_nan(tb_mod) .or. tb_mod <= 0.0) cycle
       rad_w_um = modis_planck(tb_mod, band, 1)

       ! observed clear radiance for this channel/pixel
       rad_obs = clr_obs(kban, pix)
       if (ieee_is_nan(rad_obs) .or. rad_obs <= 0.0) cycle

       sum_bias(kban) = sum_bias(kban) + (rad_w_um - rad_obs)
       n_used(kban)   = n_used(kban) + 1
    end do
  end do

  do kban = 1, nbct
     if (n_used(kban) > 0) then
        bias_w_um(kban) = sum_bias(kban) / real(n_used(kban))
     else
        bias_w_um(kban) = 0
     end if
  end do

end subroutine get_clearsky_bias

subroutine compute_multilayer_score(mod_cth_final, misr_cth_final, &
                                    badvalue_in, z_single, ml_flag)
  use, intrinsic :: ieee_arithmetic
  implicit none

  real,    intent(in)  :: mod_cth_final
  real,    intent(in)  :: misr_cth_final
  real,    intent(in)  :: badvalue_in
  real,    intent(out) :: z_single
  integer, intent(out) :: ml_flag

  real :: d_cth
  real :: mu_single, sigma_single
  real, parameter :: max_valid_cth = 20000.0

  z_single = badvalue_in
  ml_flag  = 0

  ! Only evaluate pixels where both MODIS and MISR CTH are valid.
  if (.not. ieee_is_finite(misr_cth_final)) return
  if (.not. ieee_is_finite(mod_cth_final)) return
  if (misr_cth_final <= -500.0 .or. misr_cth_final > max_valid_cth) return
  if (mod_cth_final <= 0.0 .or. mod_cth_final > max_valid_cth) return

  ! 3.5*MAD-filtered calibration constants, units = meters.
  ! D = MODIS_CTH - MISR_CTH.
  if (misr_cth_final > 4000.0) then
    mu_single    = -349.0
    sigma_single = 1971.1
  else
    mu_single    = -241.0
    sigma_single = 676.1
  end if

  if (sigma_single <= 0.0) return

  d_cth = mod_cth_final - misr_cth_final
  z_single = (d_cth - mu_single) / sigma_single

  ! One-sided multilayer evidence flag:
  ! 0 = not evaluated
  ! 1 = weak_or_no_evidence, confidence < 90%, Z < 1.2816
  ! 2 = moderate_evidence,  90% <= confidence < 95%, 1.2816 <= Z < 1.6449
  ! 3 = strong_evidence,    95% <= confidence < 99%, 1.6449 <= Z < 2.3263
  ! 4 = very_strong_evidence, confidence >= 99%, Z >= 2.3263

  if (z_single < 1.2816) then
    ml_flag = 1
  else if (z_single < 1.6449) then
    ml_flag = 2
  else if (z_single < 2.3263) then
    ml_flag = 3
  else
    ml_flag = 4
  end if

end subroutine compute_multilayer_score

subroutine compute_multilayer_score_mode(mod_cth_final, misr_cth_final, processing_mask, &
                                    badvalue_in, z_single, ml_flag)
  use, intrinsic :: ieee_arithmetic
  implicit none

  real,    intent(in)  :: mod_cth_final
  real,    intent(in)  :: misr_cth_final
  integer, intent(in)  :: processing_mask
  real,    intent(in)  :: badvalue_in
  real,    intent(out) :: z_single
  integer, intent(out) :: ml_flag

  real :: d_cth
  real :: mu_single, sigma_single
  real :: mod_bias, misr_bias
  real :: mod_precision, misr_precision
  real, parameter :: max_valid_cth = 20000.0

  z_single = badvalue_in
  ml_flag  = 0

  ! Only evaluate pixels where both MISR and final MM CTH are valid.
  if (.not. ieee_is_finite(misr_cth_final)) return
  if (.not. ieee_is_finite(mod_cth_final)) return
  if (misr_cth_final <= -500.0 .or. misr_cth_final > max_valid_cth) return
  if (mod_cth_final <= 0.0 .or. mod_cth_final > max_valid_cth) return
  ! Skip no-retrieval, MODIS-only, and MISR-only cases.
  if (processing_mask == 0 .or. processing_mask == 1 .or. processing_mask == 2) return
  ! -----------------------------------------------------------
  ! Temporary calibration constants, units = meters.
  ! Replace later with final CATS single-layer calibration table.
  !
  ! -----------------------------------------------------------
  if (misr_cth_final > 4000.0) then
    ! high cloud
    mod_bias        = -1125.0
    mod_precision   = 106
    misr_bias      = -325
    misr_precision = 403
  else
    ! low cloud
    mod_bias        = -675
    mod_precision   = 510
    misr_bias      = -325
    misr_precision = 276
  end if

  mu_single = mod_bias - misr_bias
  sigma_single = sqrt(mod_precision * mod_precision + misr_precision * misr_precision)


  d_cth = mod_cth_final - misr_cth_final
  z_single = (d_cth - mu_single) / sigma_single

  if (sigma_single <= 0.0) return
  ! One-sided multilayer evidence flag:
  ! 0 = not evaluated
  ! 1 = weak_or_no_evidence, confidence < 90%, Z < 1.2816
  ! 2 = moderate_evidence,  90% <= confidence < 95%, 1.2816 <= Z < 1.6449
  ! 3 = strong_evidence,    95% <= confidence < 99%, 1.6449 <= Z < 2.3263
  ! 4 = very_strong_evidence, confidence >= 99%, Z >= 2.3263

  if (z_single < 1.2816) then
    ml_flag = 1
  else if (z_single < 1.6449) then
    ml_flag = 2
  else if (z_single < 2.3263) then
    ml_flag = 3
  else
    ml_flag = 4
  end if
end subroutine compute_multilayer_score_mode

subroutine process_selected_pixels(wprof, tprof, hprof, psfc, pmsl, surftmp, &
                                  view, trad, rlat, rlon, landsea, misr_ctp, &
                                  misr_cth,mod_cth,mod_ctp,mod_method, mod_opt,mod_emi, &
                                  met_date, npix, &
                                  Cloud_Top_Pressure, &
                                  Cloud_Top_Height, &
                                  Cloud_Effective_Emissivity, &
                                  Cloud_Optical_Depth, &
                                  Processing_Mask, &
                                  Multilayer_ZScore, &
                                  Multilayer_flag, &
                                  MM_RetrievalStatus, &
                                  ! ------------ NEW OPTIONAL ARGS -------------
                                  clr_obs, cs_idx)
  use co2, only: nbct, plev, pp, nb_wavelen, ntbct, mbnds, jdet, freq, &
                 kch, nsct, rmin, emadj, badvalue, emisrw
  use, intrinsic :: ieee_arithmetic
  implicit none
  ! ---------- Required inputs (unchanged) ----------
  integer, intent(in) :: npix


  real,    intent(in) :: wprof(plev, npix)
  real,    intent(in) :: tprof(plev, npix)
  real,    intent(in) :: hprof(plev, npix)
  real,    intent(in) :: psfc(npix), pmsl(npix)
  real,    intent(in) :: surftmp(npix)
  real,    intent(in) :: view(npix)
  real,    intent(in) :: trad(nbct, npix)      ! (as before)
  real,    intent(in) :: rlat(npix), rlon(npix)
  integer, intent(in) :: landsea(npix)
  real,    intent(in) :: misr_ctp(npix)
  real,    intent(in) :: misr_cth(npix)
  real,    intent(in) :: mod_cth(npix)
  real,    intent(in) :: mod_ctp(npix)
  real,    intent(in) :: mod_opt(npix)
  real,    intent(in) :: mod_emi(npix)
  integer, intent(in) :: mod_method(npix)

  integer, intent(in) :: met_date(4)

  ! ---------- Outputs (unchanged) ----------
  real,    intent(out) :: Cloud_Top_Pressure(npix)
  real,    intent(out) :: Cloud_Top_Height(npix)
  real,    intent(out) :: Cloud_Effective_Emissivity(npix)
  real,    intent(out) :: Cloud_Optical_Depth(npix)
  integer, intent(out) :: Processing_Mask(npix)
  real,    intent(out) :: Multilayer_ZScore(npix)
  integer, intent(out) :: Multilayer_flag(npix)
  integer, intent(out) :: MM_RetrievalStatus(npix)

  ! ---------- NEW optional inputs ----------
  real,    intent(in), optional :: clr_obs(nbct, npix)
  integer, intent(in), optional :: cs_idx(:)

  ! ---------- Locals ----------
  real :: wprof_pix(plev), tprof_pix(plev), hprof_pix(plev)
  real :: psfc_pix, pmsl_pix, surftmp_pix, view_pix
  real :: trad_pix(nbct), trad_adj(nbct)
  integer :: met_date_pix(4), landsea_pix
  real :: rlat_pix, rlon_pix, misr_ctp_pix, misr_cth_pix, mod_cth_pix,mod_ctp_pix
  integer ::  mod_method_pix
  real :: Cloud_Top_Pressure_pix, Cloud_Top_Height_pix, Cloud_Optical_Depth_pix
  real :: Cloud_Effective_Emissivity_pix
  integer :: Processing_Mask_pix, flag_digit


  integer :: pix, i, k
  real    :: bias_w_um(nbct)
  integer :: n_used(nbct)
  logical :: use_bias
  logical :: mm_failed
  logical :: misr_cth_valid, mod_cth_valid
  integer :: digit1, digit2
  real, parameter :: max_valid_cth = 20000.0

  ! For bias computation subsets
  integer :: ncs
  real, allocatable    :: w_cs(:,:), t_cs(:,:), psfc_cs(:), surftmp_cs(:), view_cs(:), rlat_cs(:)
  integer, allocatable :: landsea_cs(:), cs_list(:)
  real, allocatable    :: clr_obs_cs(:,:)
  ! --------- keep the old file outputs (unchanged) ----------
  ! open(unit=98, file='output_1layer.txt', status='unknown', action='write')
  ! open(unit=99, file='output_2layer.txt',  status='unknown', action='write')

  ! --------- Initialize outputs (unchanged) ----------
  Cloud_Top_Pressure         = badvalue
  Cloud_Top_Height           = badvalue
  Cloud_Effective_Emissivity = badvalue
  Cloud_Optical_Depth        = badvalue
  Processing_Mask            = 1
  Multilayer_flag   = 0
  Multilayer_ZScore = badvalue
  MM_RetrievalStatus = 0

  ! ===========================================================
  ! 1) Compute scene clear-sky bias if data are provided
  ! ===========================================================
  use_bias  = present(clr_obs) .and. present(cs_idx)
  bias_w_um = 0.0
  n_used    = 0

  if (use_bias) then
    ncs = size(cs_idx)
    if (ncs > 0) then
      allocate( w_cs(plev,ncs), t_cs(plev,ncs) )
      allocate( psfc_cs(ncs), surftmp_cs(ncs), view_cs(ncs), rlat_cs(ncs) )
      allocate( landsea_cs(ncs), cs_list(ncs) )
      allocate( clr_obs_cs(nbct,ncs) )

      do i = 1, ncs
        cs_list(i)    = cs_idx(i)
        w_cs(:,i)     = wprof(:, cs_list(i))
        t_cs(:,i)     = tprof(:, cs_list(i))
        psfc_cs(i)    = psfc( cs_list(i) )
        surftmp_cs(i) = surftmp( cs_list(i) )
        view_cs(i)    = view( cs_list(i) )
        rlat_cs(i)    = rlat( cs_list(i) )
        landsea_cs(i) = landsea( cs_list(i) )
        clr_obs_cs(:,i) = clr_obs(:, cs_list(i))
      end do

      call get_clearsky_bias( w_cs, t_cs, psfc_cs, surftmp_cs, view_cs, rlat_cs, landsea_cs, &
                              met_date, clr_obs_cs, ncs, bias_w_um, n_used )

      deallocate( w_cs, t_cs, psfc_cs, surftmp_cs, view_cs, rlat_cs, landsea_cs, cs_list, clr_obs_cs )
    end if
  end if

  ! ===========================================================
  ! 2) Choose pixel list to process:
  !    - if multi_idx present -> process only those
  !    - else process all (old behavior)
  ! ===========================================================
  do pix = 1, npix
      ! extract current pixel (unchanged)
      wprof_pix   = wprof(:, pix)
      tprof_pix   = tprof(:, pix)
      hprof_pix   = hprof(:, pix)
      psfc_pix    = psfc(pix)
      pmsl_pix    = pmsl(pix)
      surftmp_pix = surftmp(pix)
      view_pix    = view(pix)
      trad_pix    = trad(:, pix)
      rlat_pix    = rlat(pix)
      rlon_pix    = rlon(pix)
      landsea_pix = landsea(pix)
      misr_ctp_pix= misr_ctp(pix)
      misr_cth_pix= misr_cth(pix)
      mod_cth_pix= mod_cth(pix)
      mod_ctp_pix= mod_ctp(pix)
      mod_method_pix= mod_method(pix)
      met_date_pix= met_date

      misr_cth_valid = ieee_is_finite(misr_cth_pix) .and. &
                       misr_cth_pix > -500.0 .and. misr_cth_pix <= max_valid_cth
      mod_cth_valid = ieee_is_finite(mod_cth_pix) .and. &
                      mod_cth_pix > 0.0 .and. mod_cth_pix <= max_valid_cth

      if (.not. misr_cth_valid .and. .not. mod_cth_valid) then
        Cloud_Top_Pressure(pix)         = -9999.0
        Cloud_Top_Height(pix)           = -9999.0
        Cloud_Effective_Emissivity(pix) = -9999.0
        Cloud_Optical_Depth(pix)        = -9999.0
        Processing_Mask(pix)            = 0
        Multilayer_flag(pix)          =  0
        cycle
      end if
      Multilayer_flag(pix)          =  0
      if (.not. misr_cth_valid) then
        Cloud_Top_Pressure(pix)         = mod_ctp(pix)
        Cloud_Top_Height(pix)           = mod_cth(pix)
        Cloud_Effective_Emissivity(pix) = mod_emi(pix)
        Cloud_Optical_Depth(pix)        = mod_opt(pix)
        Processing_Mask(pix)            = 1
        cycle
      end if
      if (.not. mod_cth_valid) then
        Cloud_Top_Pressure(pix)         = misr_ctp_pix
        Cloud_Top_Height(pix)           = misr_cth_pix
        Cloud_Effective_Emissivity(pix) = mod_emi(pix)
        Cloud_Optical_Depth(pix)        = mod_opt(pix)
        Processing_Mask(pix)            = 2
        cycle
      end if

!# Turn off misr cloud top height <=700 on Mar 17, 2026
!#      if ((mod_cth_pix - misr_cth_pix) <= 1000.0 .or. mod_method_pix >= 5 .or. misr_ctp_pix <= 700.0) then
!# Turn off everything for testing on Mar 26, 2026, Should turn back on with fixes

      call compute_multilayer_score(mod_cth_pix, misr_cth_pix, &
                              badvalue, Multilayer_ZScore(pix), Multilayer_flag(pix))

      if (mod_cth_pix < misr_cth_pix) then
          Cloud_Top_Pressure(pix)         = misr_ctp_pix
          Cloud_Top_Height(pix)           = misr_cth_pix
          Cloud_Effective_Emissivity(pix) = mod_emi(pix)
          Cloud_Optical_Depth(pix)        = mod_opt(pix)
          Processing_Mask(pix)            = 3
          cycle
      end if
      if (mod_method_pix >= 5) then
         if (Multilayer_flag(pix) >= 2 .and. Multilayer_flag(pix) <= 4) then
          Cloud_Top_Pressure(pix)         = mod_ctp_pix
          Cloud_Top_Height(pix)           = mod_cth_pix
          Cloud_Effective_Emissivity(pix) = mod_emi(pix)
          Cloud_Optical_Depth(pix)        = mod_opt(pix)
          Processing_Mask(pix)            = 5
          cycle
         else
          Cloud_Top_Pressure(pix)         = misr_ctp_pix
          Cloud_Top_Height(pix)           = misr_cth_pix
          Cloud_Effective_Emissivity(pix) = mod_emi(pix)
          Cloud_Optical_Depth(pix)        = mod_opt(pix)
          Processing_Mask(pix)            = 4
          cycle
        end if
      end if

      ! --------- NEW: robust MISR-CTP fallback (no cycle/skip) ----------
      if (.not. ieee_is_finite(misr_ctp_pix) .or. misr_ctp_pix <= 0.0 .or. &
          misr_ctp_pix < pp(1) .or. misr_ctp_pix > psfc_pix) then
         misr_ctp_pix = psfc_pix
      end if

      ! --- apply bias (if available) to all nbct channels -------------
      trad_adj = trad_pix
      if (use_bias) then
        do k = 1, nbct
          trad_adj(k) = trad_adj(k) + bias_w_um(k)
        end do
      end if

      call co2cld_onepixel_misr(wprof_pix, tprof_pix, hprof_pix, psfc_pix, pmsl_pix, &
                                surftmp_pix, view_pix, trad_adj, met_date_pix, &
                                rlat_pix, rlon_pix, landsea_pix, psfc_pix, &
                                Cloud_Top_Pressure_pix, Cloud_Top_Height_pix, &
                                Cloud_Effective_Emissivity_pix, Cloud_Optical_Depth_pix, &
                                Processing_Mask_pix)

      Cloud_Top_Pressure(pix)         = Cloud_Top_Pressure_pix
      Cloud_Top_Height(pix)           = Cloud_Top_Height_pix
      Cloud_Effective_Emissivity(pix) = Cloud_Effective_Emissivity_pix
      Cloud_Optical_Depth(pix)        = Cloud_Optical_Depth_pix
      Processing_Mask(pix)            = Processing_Mask_pix
      ! write(98, '(I8, 6F16.6, I8)') pix, psfc_pix, misr_ctp_pix, Cloud_Top_Pressure(pix), &
      !                          Cloud_Top_Height(pix), Cloud_Effective_Emissivity(pix), &
      !                          Cloud_Optical_Depth(pix), Processing_Mask(pix)

      call co2cld_onepixel_misr(wprof_pix, tprof_pix, hprof_pix, psfc_pix, pmsl_pix, &
                                surftmp_pix, view_pix, trad_adj, met_date_pix, &
                                rlat_pix, rlon_pix, landsea_pix, misr_ctp_pix, &
                                Cloud_Top_Pressure_pix, Cloud_Top_Height_pix, &
                                Cloud_Effective_Emissivity_pix, Cloud_Optical_Depth_pix, &
                                Processing_Mask_pix)

      mm_failed = .false.
      MM_RetrievalStatus(pix) = 3
      if (Cloud_Top_Pressure_pix <= Cloud_Top_Pressure(pix)) then
        Cloud_Top_Pressure(pix)         = Cloud_Top_Pressure_pix
        Cloud_Top_Height(pix)           = Cloud_Top_Height_pix * 1000.0
        Cloud_Effective_Emissivity(pix) = Cloud_Effective_Emissivity_pix
        Cloud_Optical_Depth(pix)        = Cloud_Optical_Depth_pix
        digit1 = Processing_Mask_pix / 100
        digit2 = mod(Processing_Mask_pix, 100) / 10
        if (digit1 == 1 .or. digit2 == 1) then
          Processing_Mask(pix) = 8
          MM_RetrievalStatus(pix) = 1
        else if (digit1 == 2 .or. digit2 == 2) then
          Processing_Mask(pix) = 9
          MM_RetrievalStatus(pix) = 1
        else
          MM_RetrievalStatus(pix) = 3
          mm_failed = .true.
        end if
      else
        mm_failed = .true.
        MM_RetrievalStatus(pix) = 2
      end if


      ! If the MM retrieval failed or is not accepted, fall back to MISR or MODIS
      ! depending on the MODIS-MISR multilayer evidence flag.
      if (mm_failed) then
        if (Multilayer_flag(pix) >= 2 .and. Multilayer_flag(pix) <= 4) then
          Cloud_Top_Pressure(pix)         = mod_ctp_pix
          Cloud_Top_Height(pix)           = mod_cth_pix
          Cloud_Effective_Emissivity(pix) = mod_emi(pix)
          Cloud_Optical_Depth(pix)        = mod_opt(pix)
          Processing_Mask(pix)            = 7
        else
          Cloud_Top_Pressure(pix)         = misr_ctp_pix
          Cloud_Top_Height(pix)           = misr_cth_pix
          Cloud_Effective_Emissivity(pix) = mod_emi(pix)
          Cloud_Optical_Depth(pix)        = mod_opt(pix)
          Processing_Mask(pix)            = 6
        end if
      end if
      ! write(99, '(I8, 6F16.6, I8)') pix, psfc_pix, misr_ctp_pix, Cloud_Top_Pressure(pix), &
      !                          Cloud_Top_Height(pix), Cloud_Effective_Emissivity(pix), &
      !                          Cloud_Optical_Depth(pix), Processing_Mask(pix)
  end do
  ! close(98)
  ! close(99)
end subroutine process_selected_pixels

subroutine height_to_log_pressure(geo_height, misr_cth, rows, cols, pressure_at_cth)
    implicit none
    integer, intent(in) :: rows, cols
    real(4), intent(in) :: geo_height(101, rows, cols)
    real(4), intent(in) :: misr_cth(rows, cols)
    real(4), intent(out) :: pressure_at_cth(rows, cols)

    integer :: i, j, k
    real(8) :: interp_log_p, h_below, h_above, p_below, p_above
    integer, parameter :: num_levels = 101
    real(8), parameter :: int_press_o(101) = (/4.9999999e-03, 1.6100001e-02, 3.8400002e-02, 7.6899998e-02, 1.3699999e-01, &
                  2.2440000e-01, 3.4540001e-01, 5.0639999e-01, 7.1399999e-01, 9.7530001e-01, &
                  1.2972000e+00, 1.6872000e+00, 2.1526000e+00, 2.7009001e+00, 3.3397999e+00, &
                  4.0770001e+00, 4.9204001e+00, 5.8776002e+00, 6.9566998e+00, 8.1654997e+00, &
                  9.5118999e+00, 1.1003800e+01, 1.2649200e+01, 1.4455900e+01, 1.6431801e+01, &
                  1.8584700e+01, 2.0922400e+01, 2.3452600e+01, 2.6182899e+01, 2.9121000e+01, &
                  3.2274399e+01, 3.5650501e+01, 3.9256599e+01, 4.3100101e+01, 4.7188202e+01, &
                  5.1527802e+01, 5.6125999e+01, 6.0989498e+01, 6.6125298e+01, 7.1539803e+01, &
                  7.7239601e+01, 8.3231003e+01, 8.9520401e+01, 9.6113800e+01, 1.0301720e+02, &
                  1.1023660e+02, 1.1777750e+02, 1.2564560e+02, 1.3384621e+02, 1.4238480e+02, &
                  1.5126640e+02, 1.6049590e+02, 1.7007840e+02, 1.8001830e+02, 1.9032030e+02, &
                  2.0098869e+02, 2.1202769e+02, 2.2344150e+02, 2.3523380e+02, 2.4740849e+02, &
                  2.5996909e+02, 2.7291910e+02, 2.8626169e+02, 3.0000000e+02, 3.1413690e+02, &
                  3.2867529e+02, 3.4361761e+02, 3.5896649e+02, 3.7472409e+02, 3.9089261e+02, &
                  4.0747379e+02, 4.2446979e+02, 4.4188190e+02, 4.5971179e+02, 4.7796069e+02, &
                  4.9662979e+02, 5.1571997e+02, 5.3523218e+02, 5.5516687e+02, 5.7552478e+02, &
                  5.9630621e+02, 6.1751123e+02, 6.3913977e+02, 6.6119202e+02, 6.8366730e+02, &
                  7.0656543e+02, 7.2988568e+02, 7.5362750e+02, 7.7778967e+02, 8.0237140e+02, &
                  8.2737128e+02, 8.5278802e+02, 8.7862012e+02, 9.0486591e+02, 9.3152362e+02, &
                  9.5859113e+02, 9.8606659e+02, 1.0139476e+03, 1.0422319e+03, 1.0709170e+03, &
                  1.1000000e+03 /)

    real(8) :: log_press(101)

    ! Compute log_press inside the subroutine
    do k = 1, num_levels
        log_press(k) = log(int_press_o(k))
    end do

    ! Initialize pressure_at_cth with NaN or -9999 as placeholder
    pressure_at_cth = -9999.0

    do i = 1, rows
        do j = 1, cols
            ! Skip invalid misr_cth values
            if (misr_cth(i, j) < 0) cycle

            ! Find the interpolation levels for each pixel
            do k = 1, num_levels-1
                if (misr_cth(i, j) <= geo_height(k, i, j) .and. misr_cth(i, j) >= geo_height(k+1, i, j)) then
                    h_below = geo_height(k, i, j)
                    h_above = geo_height(k+1, i, j)
                    p_below = log_press(k)
                    p_above = log_press(k+1)

                    ! Perform linear interpolation
                    interp_log_p = p_below + (misr_cth(i, j) - h_below) * (p_above - p_below) / (h_above - h_below)
                    pressure_at_cth(i, j) = exp(interp_log_p)
                end if
            end do

        end do
    end do                    
end subroutine height_to_log_pressure

integer function date2doy(iyear, imonth, iday)
  implicit none
  integer, intent(in) :: iyear, imonth, iday
  integer :: m
  integer, dimension(12) :: dmon = (/ 31, 28, 31, 30, 31, 30, &
                                     31, 31, 30, 31, 30, 31 /)

  ! --- adjust February length for leap-years --------------------------
  if ( (mod(iyear,4)   == 0 .and. mod(iyear,100) /= 0) .or. &
       (mod(iyear,400) == 0) ) dmon(2) = 29

  ! --- accumulate days for preceding months ---------------------------
  date2doy = iday
  do m = 1, imonth-1
     date2doy = date2doy + dmon(m)
  end do
  end function date2doy
