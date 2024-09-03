module modrad_2lay

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

subroutine rad_multi(wprof, tprof, psfc,pmsl, surftmp, view, rad, &
              met_date, rlat, rlon, landsea, imisr, iup, amup)

  use transmission, only: tran_modisd101
  use surfemis, only: getiremis, assign_eco_emis
  use emission_rte, only: fm_modrad_emis
  use planck_bright, only: modis_bright, modis_planck
  use co2, only: nbct, plev, leap_doy, std_doy, pp, nb_wavelen, &
                ntbct, mbnds, jdet, freq, kch, nsct, &
                rmin, emadj, badvalue
  implicit none
  save

  real wprof(plev), tprof(plev), surftmp
!f2py intent(in) wprof, tprof, surftmp
  real view, rad(nbct), eca, rlat, rlon, pmsl, psfc, amup
!f2py intent(in) view,rlat, rlon, pmsl, psfc, amup
  integer met_date(4), landsea, imisr, iup
!f2py intent(in) met_date, landsea, imisr, iup
!f2py intent(out) rad

  real z(plev), ozone(plev), zs, taup(plev, nb_wavelen), tmin, ptrp
  integer met_year, met_day, met_month, jday, ltrp, lmin

  real freqemis(nb_wavelen), emis(nb_wavelen), rho(nb_wavelen)
  real emis_out, rho_out, sfc_emis(ntbct), delmisr(nbct), delup(nbct)

  real emiswc, emis12, emis86, ppsfc, tss, rclr_s(plev)
  real rclr_s2(plev), rclr_s3(plev), ttpp(plev)

  real tw
  integer lwin, iw1, lco2

  integer is1
  real rclr(nbct), sum, db, ra(plev, nbct), robs(nbct), delr(nbct)

  integer krto(nsct), k1,k2, id, lev(nsct)
  real fmsav, rwcld(nsct), amo(nsct), bot, top, ratio
  logical ok, neg, start
  real fm1, fm2, fm

  integer ll, isp, imslp, nl, k, kban, iisp,ngch

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



! c----------------------------------------------------------------------
!
! c       Perform radiative transfer calculations for co2-slicing method.
!
! c----------------------------------------------------------------------
!
! c       Define index of first level of integration in computing RHS of
! c       co2-slicing equation.

  is1 = isp - 1
  do k = 1,nbct
!    robs(k) = modis_planck(tcold(k), mbnds(k), 0)
! c         Compute TOA clear radiance from input profiles.
    rclr(k) = fm_modrad_emis(surftmp,sfc_emis(k),psfc,pp,tprof,taup(1,k),&
                             mbnds(k),isp)
                             
    
! c         Compute components of RHS (but from MISR CTP, not surface)
    sum = 0.0
    do ll = is1,iup+1,-1
!      if (ll>imisr) then
!          sum = 0.
!      else
      db = modis_planck(tprof(ll + 1),mbnds(k),0)  &
               - modis_planck(tprof(ll),mbnds(k),0)
      sum = sum - 0.5 * (taup(ll + 1,k) + taup(ll,k)) * db
!      ra(ll,k) = sum
!      endif
    enddo
    delup(k) = sum
    
! c        Compute correction term on RHS from MISR CTP
    sum = 0.0
    do ll = is1,imisr+1,-1
      db = modis_planck(tprof(ll + 1),mbnds(k),0)  &
               - modis_planck(tprof(ll),mbnds(k),0)
      sum = sum - 0.5 * (taup(ll + 1,k) + taup(ll,k)) * db
    enddo
    delmisr(k) = sum
    
    rad(k) = rclr(k) + (1. - amup) * delmisr(k) + amup * delup(k)
    
  enddo
return

end subroutine

end module modrad_2lay