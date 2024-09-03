module granule

  implicit none
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

contains

subroutine co2cld_onepixel_misr_grn(wprof, tprof, psfc,pmsl, surftmp, view, tcold, ctp, &
              eca, met_date, rlat, rlon, landsea, misr_ctp)

  use transmission, only: tran_modisd101
  use surfemis, only: getiremis, assign_eco_emis
  use emission_rte, only: fm_modrad_emis
  use planck_bright, only: modis_bright, modis_planck
  use co2, only: nbct, plev, leap_doy, std_doy, pp, nb_wavelen, &
                ntbct, mbnds, jdet, freq, kch, nsct, &
                rmin, emadj, badvalue, emisrw
  implicit none
  save

  real tcold(nbct), wprof(plev), tprof(plev), surftmp
!f2py intent(in) tcold, wprof, tprof, surftmp
  real view, ctp(nsct), eca(2), rlat, rlon, pmsl, psfc
!f2py intent(in) view,rlat, rlon, pmsl, psfc
  integer met_date(4), landsea
!f2py intent(in) met_date, landsea
!f2py intent(out) ctp, eca
  real misr_ctp
!f2py intent(in) misr_ctp

  real z(plev), ozone(plev), zs, taup(plev, nb_wavelen), tmin, ptrp
  integer met_year, met_day, met_month, jday, ltrp, lmin

  real freqemis(nb_wavelen), emis(nb_wavelen), rho(nb_wavelen)
  real emis_out, rho_out, sfc_emis(ntbct)

  real emiswc, emis12, emis86, ppsfc, tss, rclr_s(plev)
  real rclr_s2(plev), rclr_s3(plev), ttpp(plev)

  real tw
  integer lwin, iw1, lco2

  integer is1
  real rclr(nbct), sum, db, ra(plev, nbct), robs(nbct), delr(nbct)

  integer krto(nsct), k1,k2, id
  real fmsav, rwcld(nsct), amo(nsct), bot, top, ratio(nsct)
  logical ok, neg, start
  real fm1, fm2, fm, lev(nsct)

  integer ll, isp, imslp, nl, k, kban, iisp,ngch
  integer imisr

  ngch = 0

  met_month = met_date(2)
  met_year = met_date(1)
  met_day = met_date(3)

! c     Calculate day-of-year.
  if( mod(met_year,4) .eq. 0 ) then
    jday = leap_doy(met_month) + met_day
  else
    jday = std_doy(met_month) + met_day
  end if

  call clozo101(rlat, met_month, ozone)

!     Find level of surface pressure.
  do ll = 1,plev
    if(psfc .le. pp(ll)) then
      isp = ll
      go to 100
    end if
  enddo
100  continue

! Find level of misr_ctp note this is only approximate.
! there won't be a perfect matching level.
  do ll = 1,plev
    if (misr_ctp .le. pp(ll)) then
      imisr = ll
      go to 900
    end if
  enddo
900 continue

!c     Find level of mean sea level pressure.
  do ll = 1,plev
    if(pmsl .le. pp(ll)) then
      imslp = ll
      go to 200
    end if
  enddo
200  continue

!c     Get geopotential height profile (km).
  zs = 0.0
  z = 0.0
  nl = imslp
  call height(pp,tprof,wprof,zs,nl,z)

! get surface emissivity
  call assign_eco_emis(landsea,emis,rho,freqemis)

  do k = 1, ntbct
    call getiremis(nb_wavelen,freqemis,freq(k),emis,rho,emis_out, &
                  rho_out)
    sfc_emis(k) = emis_out
  enddo
!c     Calculate transmittance profiles (101 level fast model).
  do k = 1, ntbct
    kban = mbnds(k)
    call tran_modisd101(met_year, jday, tprof, wprof, ozone, view, &
                       kban, jdet, taup(1,k), plev)
  enddo


!c----------------------------------------------------------------------

!c     Find tropopause level.

!c----------------------------------------------------------------------

!c     Find level of coldest temperature between 100 mb and the surface.
  tmin = 350.0
  do ll = 44,isp
    if(tprof(ll) .le. tmin) then
      tmin = tprof(ll)
      lmin = ll
    end if
  enddo

! c     Look for inflection point in temperature profile.
! c     If found do not change tropopause level.
! c     In isothermal conditions, choose the bottom level.
! c     If temperature increases monotonically below 100 mb, do not
! c     adjust tropopause level.

  ltrp = -99
  do ll = lmin, 71
    if(tprof(ll-1) .ge. tprof(ll) .and. tprof(ll+1) .gt. tprof(ll)) then
      ltrp = ll
      go to 300
    end if
  enddo

300  continue
  if(ltrp .eq. -99) ltrp = lmin
  ptrp = pp(ltrp)

! c----------------------------------------------------------------------
! c     Get IRW (11 micron) radiance profile, then convert to BT.
   do k = ltrp,isp
    if(k .eq. isp) then
      emiswc = sfc_emis(5)
      emis12 = sfc_emis(7)
      emis86 = sfc_emis(6)
      ppsfc = psfc
      tss = surftmp
    else
      emiswc = 1.0
      emis12 = 1.0
      emis86 = 1.0
      ppsfc = pp(k)
      tss = tprof(k)
    end if
! c       tss = tp(k)
    iisp = k
! 11 micron.
    rclr_s(k) = fm_modrad_emis(tss,emiswc,ppsfc,pp,tprof,taup(1,5), &
                            mbnds(5),iisp)
    ttpp(k) = modis_bright(rclr_s(k),mbnds(5),0)
! c       Also get 12 and 8.6 micron radiance profile.
    rclr_s2(k) = fm_modrad_emis(tss,emis12,ppsfc,pp,tprof,taup(1,7), &
                              mbnds(7),iisp)
    rclr_s3(k) = fm_modrad_emis(tss,emis86,ppsfc,pp,tprof,taup(1,6), &
                              mbnds(6),iisp)
  enddo


! c----------------------------------------------------------------------
!
! c     Compute the "IR window" cloud height.  Compare 11 micron BT to
! c     atmospherically-corrected temperatures from radiance profile.
!
! c----------------------------------------------------------------------

  tw = tcold(5)

  lwin = ltrp
  do ll = ltrp,isp
    if(ttpp(ll) .lt. tw) then
      lwin = ll
    end if
  enddo

! c     Find which level tw is closest to and define beginning level at
! c     which to begin comparison of LHS and RHS of co2-slicing equation.
  if(lwin .eq. isp) then
    iw1 = lwin - 1
  else if(lwin .eq. ltrp) then
    iw1 = lwin
  else if( abs(tw - ttpp(lwin)) .gt. abs(tw - ttpp(lwin + 1)) ) then
    lwin = lwin + 1
    iw1 = lwin - 1
  else
    iw1 = lwin
  end if
! CCC iw1 appears to be the beginning level of co2-slicing equation.
  lco2 = lwin


! c----------------------------------------------------------------------
!
! c       Perform radiative transfer calculations for co2-slicing method.
!
! c----------------------------------------------------------------------
!
! c       Define index of first level of integration in computing RHS of
! c       co2-slicing equation.

! now we modify things to use MISR ctp as the bottom surface instead of actual surface.
!  is1 = isp - 1
  is1 = imisr - 1
  do k = 1,nbct
    robs(k) = modis_planck(tcold(k), mbnds(k), 0)
! c         Compute TOA clear radiance from input profiles.
    rclr(k) = fm_modrad_emis(tprof(imisr),1.0,misr_ctp,pp,tprof,taup(1,k),&
                             mbnds(k),imisr) !clear is now 'reference' which includes misr cloud contribution
    delr(k) = robs(k) - rclr(k)
! c         Compute components of RHS
    sum = 0.0
    do ll = is1,1,-1 !this is modified to calculate only from misr_ctp upwards rather than from surface.
      db = modis_planck(tprof(ll + 1),mbnds(k),0)  &
               - modis_planck(tprof(ll),mbnds(k),0)
      sum = sum - 0.5 * (taup(ll + 1,k) + taup(ll,k)) * db
      ra(ll,k) = sum
    enddo
  enddo

! c       Compute cloud height for each channel combination (when possible).

  do id = 1,nsct

    krto(id) = 0
    fmsav = 1000.0
    k1 = kch(id,1)
    k2 = kch(id,2)

! c         Check if values of cold minus warm are within instrument noise.
    if(delr(k1) .le. rmin(k1)) then
      if(delr(k2) .le. rmin(k2)) then

! c               Find minimum difference between LHS and RHS of co2-slicing
! c               equation. Apply emissivity correction.

          ok = .false.
!          do ll = ltrp, iw1
          do ll = ltrp, imisr !now the maximum pressure to search for solutions is misr_ctp.
            fm1 = delr(k1) / delr(k2)
            fm2 = (ra(ll,k1) * emadj(id)) / ra(ll,k2)
            write(*,*) pp(ll), fm1, fm2
            fm = fm1 - fm2
            if(fm .lt. 0.0) then
              neg = .true.
            else
              neg = .false.
            end if
            if(ll .eq. ltrp) start = neg
            if(neg .neqv. start) then

              ok = .true.
              lev(id) = ll - 1
              fmsav = abs(fm)
              go to 2000
            end if

!  c                 write(*,'(4i5,3f12.5,2l5)') id,k1,k2,ll,fm1,fm2,fm,neg,start

          enddo
2000           continue

          if(fmsav .lt. 1000.0) then
            krto(id) = 1
          end if

          if( (.not. ok) ) krto(id) = 0
          if( (.not. ok) .and. neg) then
            krto(id) = 1
            lev(id) = iw1
          end if

      end if
    end if

  enddo

!
! ! c       Compute effective cloud emissivities for every channel combination
! ! c       for which there was a successful cloud height retrieval.

        eca(1) = rclr(5)
        eca(2) = delr(5)



  do id=1,nsct
    ctp(id) = pp(lev(id))

  end do
  return
end subroutine

end module granule

subroutine mm_cth(dims, wpr, tpr, ps, pm, sft, vza, tcld, ctp, &
              eca, met_date, lat, lon, surftype, misr_ctp, &
              mm_flag)
              
  use co2, only: nbct, plev, nsct
  use granule, only: co2cld_onepixel_misr_grn
  
  implicit none
  save

  integer dims
!f2py intent(in) dims
  real tcld(dims, nbct), wpr(dims, plev), tpr(dims, plev)
!f2py intent(in) tcld, wpr, tpr
  real vza(dims), ctp(dims, nsct), eca(2), lat(dims), lon(dims)
!f2py intent(in) vza, lat, lon
  real sft(dims), pm(dims), ps(dims)
!f2py intent(in) sft, lon, pm, ps
  integer met_date(4), surftype(dims), mm_flag(dims)
!f2py intent(in) met_date, surftype, mm_flag
!f2py intent(out) ctp, eca
  real misr_ctp(dims)
!f2py intent(in) misr_ctp

  integer i
  ctp(:,:) = 0
  eca      = 0

  do i = 1, dims
      if(mm_flag(i).eq.-1) then
!  outside the MISR swath
          ctp(i,:) = -9999
      else if(mm_flag(i).eq.0) then
!  not eligible to run MM_CTH
          ctp(i,:) = -999
      else
!  eligible to run MM_CTH
          call co2cld_onepixel_misr_grn(wpr(i,:), tpr(i,:), ps(i), &
                      pm(i), sft(i), vza(i), tcld(i,:), ctp(i,:), &
                   eca, met_date, lat(i), lon(i), surftype(i), &
                   misr_ctp(i))

        end if
    end do    

end subroutine