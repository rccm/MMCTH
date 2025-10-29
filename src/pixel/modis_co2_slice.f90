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

!  c     Ice cloud emissivity adjustments for each channel combination.
  data emadj / 1.000, 1.000, 1.000, 1.000 /

  data jdet /0/

  data emisrw /1.0/

  data std_doy  /0,31,59,90,120,151,181,212,243,273,304,334/
  data leap_doy /0,31,60,91,121,152,182,213,244,274,305,335/

  data mbnds/36,35,34,33,31,29,32/
  data kch /1,2,2,3,2,3,4,4/

!     Define wavenumbers of CO2 bands used (36, 35, 34, 33, 31).
  data freq /7.055309E+02, 7.196677E+02, 7.317089E+02,&
                  7.483224E+02, 9.081998E+02, 1.173198E+03,&
                  8.315149E+02/

  data rmin /-1.0, -1.0, -100.0, -1.0, -0.5/

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


subroutine co2cld_onepixel_misr(wprof, tprof, psfc, pmsl, surftmp, view, trad, &
                                met_date, rlat, rlon, landsea, misr_ctp, &
                                Cloud_Top_Pressure, &
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
!   - Eliminates 11, 12, 8.6 μm radiance profiles, as MISR CTP defines lower boundary.
!   - Assumes clozo101 and tran_modisd101 are functional for ozone and transmittance.
!   - Variable names align with MODIS PGE (MOD06) conventions.
!   - Debug messages aid troubleshooting (e.g., invalid inputs, retrieval steps).
!
! Author: Arka Mita (rewritten by Guangyu Zhao)
! Date: April 2x, 2025
!-------------------------------------------------------------------------------

  use co2, only: nbct, plev, pp, nb_wavelen, ntbct, mbnds, jdet, freq, &
                 kch, nsct, rmin, emadj, badvalue, emisrw
  use, intrinsic :: ieee_arithmetic
  use transmission, only: tran_modisd101, wstd,tstd
  use surfemis, only: getiremis, assign_eco_emis
  use emission_rte, only: fm_modrad_emis
  use planck_bright, only: modis_bright, modis_planck
  implicit none
! Inputs
  real, intent(in) :: wprof(plev)               ! Water vapor mixing ratio (g/kg)
  real, intent(in) :: tprof(plev)               ! Temperature profile (K)
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
  real :: taup(plev,ntbct)                     ! Optical depth
  real :: freqemis(nb_wavelen), emis(nb_wavelen), rho(nb_wavelen)
  real :: emis_out, rho_out
  real :: lev(nsct),fmsav
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
  real, parameter :: xi_vis = 2.56          ! Minnis et al. 1990, ice clouds
  real :: tau_vis 
  
  
! Initialize outputs
  Cloud_Top_Pressure = badvalue
  Cloud_Effective_Emissivity = badvalue
  Cloud_Optical_Depth = badvalue
  CO2_Slicing_Flag = 0
  Processing_Mask = 3
  ctp_flag = 0
  ctp_pres = badvalue
   
  do k = 1, nbct
     tcold(k) = modis_bright(trad(k),mbnds(k), 1)
  end do 
!-------------------------------------------------------------------------------
! Validate inputs
!-------------------------------------------------------------------------------
  if (misr_ctp <= 0.0 .or. misr_ctp > psfc .or. misr_ctp < pp(1)) then
    write(*,'(A,F10.2,F10.2,F10.2,A)') 'ERROR: Invalid misr_ctp = ', misr_ctp, psfc, pp(1), &
                           ' hPa, returning badvalue'
    return
  endif
  if (psfc <= pp(1) .or. psfc > 1100.0) then
    write(*,'(A,F10.2,A)') 'ERROR: Invalid psfc = ', psfc, ' hPa, returning badvalue'
    return
  endif
  if (any(tcold <= 0.0)) then
    write(*,'(A)') 'ERROR: Invalid tcold values, returning badvalue'
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
    write(*,'(A,F10.2,A)') 'ERROR: Cannot find isp for psfc = ', psfc, &
                           ' hPa, returning badvalue'
    return
  endif
! MISR CTP level
  imisr = 0
  do ll = 1, plev
    if (misr_ctp <= pp(ll)) then
      imisr = ll
      exit
    end if
  enddo
  if (imisr == 0) then
    write(*,'(A,F10.2,A)') 'ERROR: Cannot find imisr for misr_ctp = ', misr_ctp, &
                           ' hPa, returning badvalue'
    return
  endif

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
    write(*,'(A,F10.2,A)') 'WARNING: Using isp for imslp, pmsl = ', pmsl, ' hPa'
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
    write(*,'(A,I3,A,F10.2,A)') 'ERROR: Tropopause level ', ltrp, &
                                ' at ', ptrp, ' hPa below misr_ctp, returning badvalue'
    return
  endif

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

  do ll = isp, plev       ! layers that are physically below the terrain
    tpad(ll)  = tprof(isp-1) ! hold at surface T
    wpad(ll)  = 0.0        ! zero vapour → negligible extra τ
    ozpad(ll) = ozprof(isp-1)
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
  
  ! Write taup profile into a file (should comment it off during processing) 
  ! write(fname,'("taup_",I4.4,"_",I3.3,".bin")') met_year, jday 
  ! open(newunit=newunit, file=fname, status='replace', access='stream', &
  !    form='unformatted', action='write')
  ! write(newunit) taup(1:plev,1:nbct)  ! taup is REAL(KIND=8) by default in your module
  ! close(newunit)

!-------------------------------------------------------------------------------
!  ⬛  CO₂‑slicing search  (ATBD Eq. 11–13)  ‑‑ with NaN/Inf safeguards ⬛
!-------------------------------------------------------------------------------
  
  ! 1.  Index of MISR low‑cloud pressure (Pl)
  imisr = 0
  do ll = 1, plev
     if (misr_ctp <= pp(ll)) then
        imisr = ll
        exit
     end if
  end do
  if (imisr == 0) return          ! safety: MISR CTP above model top
  is1  = imisr - 1
  
  ! 2.  Clear / cloudy radiances and reference terms  -----------------------------
  do k = 1, nbct
     robs(k) = modis_planck( tcold(k), mbnds(k), 0 )
  
     rclr(k) = fm_modrad_emis( tpad(imisr), 1.0, misr_ctp, pp, tpad, &
                               taup(:,k), mbnds(k), imisr )
  
     delr(k) = robs(k) - rclr(k)
  
     ! --- accumulate RHS integral from Pl up to TOA -----------------------------
     sum = 0.0
     do ll = is1, 1, -1
        db  = modis_planck( tpad(ll+1), mbnds(k), 0 ) - &
              modis_planck( tpad(ll  ), mbnds(k), 0 )
        sum = sum - 0.5 * ( taup(ll+1,k) + taup(ll,k) ) * db
        ra(ll,k) = sum
     end do
  end do
  
  write(*,'("Pressure levels: trop=",i3,"  Pl=",i3,"  ptrop=",f7.2,"  Pmisr=",f7.2)') &
          ltrp, imisr, pp(ltrp), pp(imisr)
  
  ! 3.  Zero‑crossing search for each CO₂‑band pair  ------------------------------
  lev  = badvalue
  krto = 0
  
  pair_loop: do id = 1, nsct
     k1 = kch(id,1)
     k2 = kch(id,2)
  
     ! Skip pair if clear–cloud signals are below noise threshold
     if (delr(k1) > rmin(k1)+100 .or. delr(k2) > rmin(k2)+100) cycle
  
    !  write(*,'("Working on Pair:",i2)') id
     ok    = .false.
     start = .false.
  
     do ll = ltrp, imisr              ! search downward until the MISR layer
        if (delr(k2) == 0.0  .or. ra(ll,k2) == 0.0) cycle     ! avoid /0
  
        fm1 = delr(k1)            / delr(k2)
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
        ctp_pres(id) = pp(lev(id))
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
        lco2 = nint(lev(id))
        pfco2 = ctp_pres(id)
        found_solution = .true.
        Cloud_Top_Pressure = pfco2 
        Processing_Mask = 5
        ! write(*,'(A,I2,A,F10.2)') 'INFO: Selected band pair 36/35 (id=', id, &
        !                           ') CTP: ', pfco2
        exit
      else if (id == 3 .and. ctp_pres(id) < 650.0 .and. ctp_pres(id) < misr_ctp) then
        ! Band pair 35/33
        ipco2 = id
        lco2 = nint(lev(id))
        pfco2 = ctp_pres(id)
        found_solution = .true.
        Cloud_Top_Pressure = pfco2 
        Processing_Mask = 5
        ! write(*,'(A,I2,A,F10.2)') 'INFO: Selected band pair 35/33 (id=', id, &
        !                           ') CTP: ', pfco2
        exit
      end if
    end if
  end do 
  
! If no valid solution found, check fallback

!===============================================================================
!  EFFECTIVE CLOUD EMISSIVITY  (11 µm window)  –  UNIT‑CONSISTENT VERSION
!===============================================================================
  if (found_solution .and. imisr - lco2 < 1 ) then
    Cloud_Top_Pressure        = pfco2 
    kwc = 5                     ! 11‑µm window band (band‑31)
    ll  = lco2                  ! level index of selected upper cloud
 
    ! -- 0. Recompute OBSERVED radiance at 11 µm in W m‑2 sr‑1 µm‑1 ------------
    !     tcold(kwc) is the brightness temperature already derived earlier
    Iobs = modis_planck( tcold(kwc), mbnds(kwc), 0 )
 
    ! -- 1. Build a true TRANSMITTANCE array for 11 µm --------------------------
    do k = 1, plev
       T11(k) = taup(k,kwc)                ! already transmittance
       if (T11(k) > 1.0) T11(k) = 1.0      ! safety clamp
       if (T11(k) < 0.0) T11(k) = 0.0
    end do
    ! write(*,'("tau(Pl)=",f6.3," tau(Pu)=",f6.3," tau(TOA)=",f6.3)') &
    !     T11(imisr), T11(lco2), T11(1)
    ! -- 2. Simulated cloudy and warm‑reference radiances ----------------------
    write(*,*) ll, imisr
    ! write(*,*) T11(1),imisr,ll,pp(ll),isp, ltrp
    ! write(*,'("DEBUG levels  Pl=",i3,"  Psfc=",i3,"  Trop=",i3)') &
          ! imisr, isp, ltrp
    rwcld = fm_modrad_emis( tpad(ll), emisrw, pp(ll),   pp, tpad, T11,   &
                            mbnds(kwc), ll )            ! Ic(Pu)
    rwarm = fm_modrad_emis( tpad(imisr), emisrw, pp(imisr), pp, tpad, T11,&
                            mbnds(kwc), imisr )          ! Ic(Pl)

    if (imisr < 1 .or. imisr >= plev) then
      write(*,'("SKIP pixel: imisr=",i4," outside 1..",i4)') imisr, plev
      return
    end if
                     
    ! -- 3. ΔI integral from Pl down to Ps (Eq. 12) ----------------------------
    deltaI = 0.0
    do k = imisr, isp-1
      if (k+1 > plev) exit   
      db = modis_planck( tpad(k+1), mbnds(kwc), 0 ) -  &
           modis_planck( tpad(k  ), mbnds(kwc), 0 )
      deltaI = deltaI + 0.5*(T11(k)+T11(k+1)) * db
    end do
 
    ! -- 4. Numerator & denominator of Eq. 14 ----------------------------------
    num = (Iobs  - rwarm) - deltaI          ! should be ≤ 0
    den = (rwcld - rwarm) - deltaI          ! should be ≤ 0
 
    ! -- 5. DEBUG PRINTS --------------------------------------------------------
    ! write(*,'("----- EMISSIVITY DEBUG -----")')
    ! write(*,'("Pu lvl =",i4,"  Pl lvl =",i4)') ll, imisr
    ! write(*,'("Iobs =",f9.3,"  rwarm =",f9.3,"  rwcld =",f9.3)')  &
    !       Iobs, rwarm, rwcld
    ! write(*,'("deltaI=",f9.3,"  num =",f9.3,"  den =",f9.3)')    &
    !       deltaI, num, den
    ! write(*,'("Tb_obs =",f7.2,"  Tb_warm =",f7.2,"  Tb_cld =",f7.2)') &
          ! modis_bright(Iobs , mbnds(kwc), 0),                    &
          ! modis_bright(rwarm, mbnds(kwc), 0),                    &
          ! modis_bright(rwcld, mbnds(kwc), 0)
 
    ! -- 6. Emissivity solution -------------------------------------------------
    if (abs(den) > 0.1 .and. num <= 0.0 .and. den < 0.0) then
       ratio = num / den                       ! εA  (0…1)
       if (ratio >= 0.01 .and. ratio <= 1.0) then
          Cloud_Effective_Emissivity = ratio * 1
          CO2_Slicing_Flag          = 1
          ! Processing_Mask           = 2
          ! write(*,'("SUCCESS  CTP=",f8.2," hPa   ECA=",f6.3)')   &
          !       Cloud_Top_Pressure, Cloud_Effective_Emissivity
       else
          ! write(*,'("WARNING  εA out of range: ",f6.3)') ratio
          ! Processing_Mask = 1
          return
       end if
    else
      !  write(*,'("WARNING  emissivity sign / magnitude problem")')
      !  Processing_Mask = 1
       return
    end if
 
 !===============================================================================  
  !  εA is ratio (window emissivity) * ecawin  (ecawin = 1.0 here)
    if (Cloud_Effective_Emissivity >= 0.9999) then
      tau_vis = 99.9            ! saturation flag for optically thick cloud
    else if (Cloud_Effective_Emissivity <= 0.0) then
      tau_vis = badvalue             ! fully transparent (should not occur)
    else
      tau_vis = -log( 1.0 - Cloud_Effective_Emissivity ) / xi_vis
    end if
    Cloud_Optical_Depth = tau_vis
end if

end subroutine co2cld_onepixel_misr


subroutine process_selected_pixels(wprof, tprof, psfc, pmsl, surftmp, view, trad, &
                                  rlat, rlon, landsea, misr_ctp, met_date, npix, &
                                  Cloud_Top_Pressure, &
                                  Cloud_Effective_Emissivity, &
                                  Cloud_Optical_Depth, &
                                  Processing_Mask)
  use co2, only: nbct, plev, pp, nb_wavelen, ntbct, mbnds, jdet, freq, &
                 kch, nsct, rmin, emadj, badvalue, emisrw
  implicit none
! Input arrays for selected pixels
  integer, intent(in) :: npix                         ! Number of selected pixels
  real, intent(in) :: wprof(plev, npix)              ! Water vapor profiles
  real, intent(in) :: tprof(plev, npix)              ! Temperature profiles
  real, intent(in) :: psfc(npix), pmsl(npix)         ! Surface, MSL pressure
  real, intent(in) :: surftmp(npix)                  ! Surface temperature
  real, intent(in) :: view(npix)                     ! View angle
  real, intent(in) :: trad(nbct, npix)              ! Brightness temperatures
  real, intent(in) :: rlat(npix), rlon(npix)         ! Latitude, longitude
  integer, intent(in) :: landsea(npix)               ! Land/sea flag
  real, intent(in) :: misr_ctp(npix)                 ! MISR CTP
  integer, intent(in) :: met_date(4)                 ! Year, month, day, hour

! Output arrays
  real, intent(out) :: Cloud_Top_Pressure(npix)      ! Upper cloud CTP
  real, intent(out) :: Cloud_Effective_Emissivity(npix) ! Effective emissivity
  real, intent(out) :: Cloud_Optical_Depth(npix)  ! Cloud Optical Depth
  integer, intent(out) :: Processing_Mask(npix)      ! Quality flag

! Local variables for a single pixel
  real :: wprof_pix(plev), tprof_pix(plev)
  real :: psfc_pix, pmsl_pix, surftmp_pix, view_pix
  real :: trad_pix(nbct)
  integer :: met_date_pix(4), landsea_pix
  real :: rlat_pix, rlon_pix, misr_ctp_pix
  real :: Cloud_Top_Pressure_pix, Cloud_Top_Height_pix, Cloud_Optical_Depth_pix
  real :: Cloud_Effective_Emissivity_pix
  integer :: CO2_Slicing_Flag_pix, Processing_Mask_pix
  integer :: pix

! Initialize outputs
  Cloud_Top_Pressure = badvalue
  Cloud_Effective_Emissivity = badvalue
  Cloud_Optical_Depth = badvalue
  Processing_Mask = 1

! Loop over selected pixels
  do pix = 1, npix
    ! Extract data for the current pixel
    wprof_pix = wprof(:, pix)
    tprof_pix = tprof(:, pix)
    psfc_pix = psfc(pix)
    pmsl_pix = pmsl(pix)
    surftmp_pix = surftmp(pix)
    view_pix = view(pix)
    trad_pix = trad(:, pix)
    rlat_pix = rlat(pix)
    rlon_pix = rlon(pix)
    landsea_pix = landsea(pix)
    misr_ctp_pix = misr_ctp(pix)
    met_date_pix = met_date

    ! Validate MISR CTP
    if (misr_ctp_pix <= 0.0 .or. misr_ctp_pix > psfc_pix .or. &
        misr_ctp_pix < pp(1)) cycle

    ! Call the subroutine for the current pixel
    
    call co2cld_onepixel_misr(wprof_pix, tprof_pix, psfc_pix, pmsl_pix, &
                              surftmp_pix, view_pix, trad_pix, met_date_pix, &
                              rlat_pix, rlon_pix, landsea_pix, misr_ctp_pix, &
                              Cloud_Top_Pressure_pix, &
                              Cloud_Effective_Emissivity_pix, &
                              Cloud_Optical_Depth_pix, &
                              Processing_Mask_pix)

    ! Store results
    Cloud_Top_Pressure(pix) = Cloud_Top_Pressure_pix
    Cloud_Effective_Emissivity(pix) = Cloud_Effective_Emissivity_pix
    Cloud_Optical_Depth(pix) = Cloud_Optical_Depth_pix
    Processing_Mask(pix) = Processing_Mask_pix
  end do

end subroutine process_selected_pixels

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