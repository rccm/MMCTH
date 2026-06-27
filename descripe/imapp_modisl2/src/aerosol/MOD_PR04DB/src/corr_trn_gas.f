!c    program corr_trn_gas
      subroutine corr_trn_gas(sza_in,vza_in,pwv_in,to3_in,   
     1                       srb1,srb2,srb4,srb9,srb10)  
!c
!c
!c  Written by Myeong-Jae Jeong (MJ)
!c  Ver. 1.1 Aug 4, 2011
!c
!c     This code removes the contributions of gas absorption
!c     from MODIS atmospheric bands (B1 through B7) reflectances using 
!c     PW and ozone data (from NCEP Reanalysis dataset)
!c.    This code inherited from "reform_aot500_5wav_allsites_v5.pro"
!c
!c     Routines based on "Algorithm for remote sensing of tropospheric aerosol
!c     from MODIS: Collection 005: Revision 2; Feb 2009" by Remer et al.
!c     MODIS ATBD for MOD04/MYD04.
!c
!c... Inputs
!c...    sza_in, vza_in: solar and sensor zenith angles (deg)
!c...    pwv_in: precipitable water (cm or g/cm^2)
!c...    to3_in: total ozone (Dobson Unit)
!c...    srb1,srb2,srb4,srb9,srb10: TOA reflectances (before corr.)
!c...
!c... Output
!c...    srb1,srb2,srb4,srb9,srb10: TOA reflectances (After corr.)
!c...
!c...
      real*4 sza_in, vza_in, wvl_ib
      real*4 pwv_in, to3_in, trn_h2o, trn_o3, trn_co2
      real*8 tg_b1,tg_b2,tg_b3,tg_b4,tg_b5,tg_b6,tg_b7
      real*8 srb1,srb2,srb4,srb9,srb10
      integer iband, iop_clim_h2o, iop_clim_o3

!c===================================
!c  Input parameters to compute transmission factor
!c-----------------------------------
!c  set to 1 to use simple climatology instead of obs/reanalysis data
      iop_clim_h2o=0 
      iop_clim_o3=0
!c===================================
c. examples of input data format
c.    sza_in=30.0          ! solar zenith angle in deg.
c.    vza_in=45.0          ! satellite viewing zenith angle in deg.
c.    pwv_in=2.0           ! precipitable water in cm unit(or g/cm^2) (<--NCEP data)
c.    to3_in=300.0         ! assume constant total ozone amount (300DU)

!c... Initialization of transmittance factors
      tg_b1=1.0
      tg_b3=1.0
      tg_b7=1.0
      tg_b2=1.0
      tg_b5=1.0

!c==========================================================================
!c   call subroutines to compute transmission factor
!c--------------------------------------------------------------------------
!     iband=1  !; MODIS B01
      call get_h2o_trn(1, iop_clim_h2o, sza_in, vza_in, pwv_in, 
     1                 trn_h2o, wvl_ib)
      call get_o3_trn(1, iop_clim_o3, sza_in, vza_in, to3_in, 
     1                trn_o3, wvl_ib)
      call get_co2_trn(1, sza_in, vza_in, trn_co2, wvl_ib)
      tg_b1=trn_h2o*trn_o3*trn_co2  !; compute transmission factor due to H2O, O3, and CO2
!--------------------------
!     iband=3  !; MODIS B03
      call get_h2o_trn(3, iop_clim_h2o, sza_in, vza_in, pwv_in, 
     1                 trn_h2o, wvl_ib)
      call get_o3_trn(3, iop_clim_o3, sza_in, vza_in, to3_in, 
     1                trn_o3, wvl_ib)
      call get_co2_trn(3, sza_in, vza_in, trn_co2, wvl_ib)
      tg_b3=trn_h2o*trn_o3*trn_co2  !; compute transmission factor due to H2O, O3, and CO2
!--------------------------
!     iband=7  !; MODIS B07
      call get_h2o_trn(7, iop_clim_h2o, sza_in, vza_in, pwv_in, 
     1                 trn_h2o, wvl_ib)
      call get_o3_trn(7, iop_clim_o3, sza_in, vza_in, to3_in, 
     1                trn_o3, wvl_ib)
      call get_co2_trn(7, sza_in, vza_in, trn_co2, wvl_ib)
      tg_b7=trn_h2o*trn_o3*trn_co2  !; compute transmission factor due to H2O, O3, and CO2
!--------------------------
!     iband=2  !; MODIS B02 ; added on Sep 17, 2009
      call get_h2o_trn(2, iop_clim_h2o, sza_in, vza_in, pwv_in, 
     1                 trn_h2o, wvl_ib)
      call get_o3_trn(2, iop_clim_o3, sza_in, vza_in, to3_in, 
     1                trn_o3, wvl_ib)
      call get_co2_trn(2, sza_in, vza_in, trn_co2, wvl_ib)
      tg_b2=trn_h2o*trn_o3*trn_co2  !; compute transmission factor due to H2O, O3, and CO2
!--------------------------
!     iband=5  !; MODIS B05 ; added on Sep 17, 2009
      call get_h2o_trn(5, iop_clim_h2o, sza_in, vza_in, pwv_in, 
     1                 trn_h2o, wvl_ib)
      call get_o3_trn(5, iop_clim_o3, sza_in, vza_in, to3_in, 
     1                trn_o3, wvl_ib)
      call get_co2_trn(5, sza_in, vza_in, trn_co2, wvl_ib)
      tg_b5=trn_h2o*trn_o3*trn_co2  !; compute transmission factor due to H2O, O3, and CO2
!--------------------------

!c    ************************************
!c      usage:  refl_crr=refl_obs*trn_gas
!c    ************************************
!cc...**********************************************
!cc...***     end of gas correction routines     *** ; --> v3
!cc...**********************************************
!c==========================================================================

!      print *, 'tg_b1, tg_b3, tg_b7, tg_b2, tg_b5'
!      print *, tg_b1, tg_b3, tg_b7, tg_b2, tg_b5

!c Final reflectance corrections for gas absroption effect
      srb1=tg_b1*srb1    ! MODIS B01
      srb2=tg_b3*srb2    ! MODIS B03
      srb4=tg_b7*srb4    ! MODIS B07
      srb9=tg_b2*srb9    ! MODIS B02
      srb10=tg_b5*srb10  ! MODIS B05
c...        myd_rb01out=tg_b1*myd_rb01[i]
c...        myd_rb03out=tg_b3*myd_rb03[i]
c...        myd_rb07out=tg_b7*myd_rb07[i]
c...        myd_rb02out=tg_b2*myd_rb02[i]  
c...        myd_rb05out=tg_b5*myd_rb05[i]  

!c--------------------------------------------------
!c--------------------------------------------------
      return
      end
!c--------------------------------------------------


!c--------------------------------------------------
      subroutine get_h2o_trn(iband, iop_clim_h2o, sza_in, 
     1           vza_in, pwv_in, trn_h2o, wvl_ib)

!c  Written by Myeong-Jae Jeong (MJ)
!c  Ver. 1.1 Aug 4, 2011
!c
      real*4 wvlarr(7), k1_h2o(7), k2_h2o(7), k3_h2o(7), tau_h2o(7)
      real*4 am_sat, pi, xdtor, trn_h2o, wvl_ib, sza_in, vza_in
      real*4 pwv_in
      integer ib, iband, iop_clim_h2o

      pi=3.141592
      xdtor=pi/180.0
!c    ib=iband-1  !; IDL
      ib=iband    !; FORTRAN
      
      data wvlarr /0.65,0.86,0.47,0.55,1.24,1.64,2.12/
      data k1_h2o /-5.73888,  -5.32960,  0.00000, 0.00000, 
     +             -6.39296,  -7.76288,  -4.05388/
      data k2_h2o /0.925534,  0.824260, 0.00000, 0.00000,  
     +             0.942186,  0.979707,  0.872951/
      data k3_h2o /-0.0188365,-0.0277443,0.00000, 0.00000, 
     +             -0.0131901, 0.007784, -0.0268464/
      data tau_h2o /1.543e-2,  1.947e-2, 0.0000,  0.0000,   
     +              1.184e-2,  9.367e-3,  5.705e-2/

      am_sat=1.0/cos(xdtor*sza_in) + 1.0/cos(xdtor*vza_in)
      wvl_ib=wvlarr(ib)

      if(iop_clim_h2o.eq.1) then 
         trn_h2o=exp(am_sat*tau_h2o(ib))
      else
         trn_h2o=exp(exp(k1_h2o(ib) + k2_h2o(ib)*alog(am_sat*pwv_in)  
     1               + k3_h2o(ib)*(alog(am_sat*pwv_in))**2))
      endif
      if(iband.eq.3.or.iband.eq.4) trn_h2o=1.0

!c--------------------------------------------------
!c--------------------------------------------------
      return
      end
!c--------------------------------------------------


!c--------------------------------------------------
      subroutine get_o3_trn(iband, iop_clim_o3, sza_in, vza_in, 
     1           to3_in, trn_o3, wvl_ib)

!c  Written by Myeong-Jae Jeong (MJ)
!c  Ver. 1.1 Aug 4, 2011
!c
      real*4 wvlarr(7), k_o3(7), tau_o3(7)
      real*4 am_sat, pi, xdtor, trn_co2, wvl_ib, sza_in, vza_in
      integer ib, iband, iop_clim_o3

      pi=3.141592
      xdtor=pi/180.0
!c    ib=iband-1  !; IDL
      ib=iband    !; FORTRAN
      
      data wvlarr /0.65,0.86,0.47,0.55,1.24,1.64,2.12/
      data k_o3 /5.09e-5, 0.00000, 4.26e-6, 1.05e-4, 
     +           0.00000, 0.00000, 0.00000/
      data tau_o3 /2.478e-2,0.00000, 2.432e-3,2.957e-2,
     +           0.00000, 0.00000, 0.00000/

      am_sat=1.0/cos(xdtor*sza_in) + 1.0/cos(xdtor*vza_in)
      wvl_ib=wvlarr(ib)
      if(iop_clim_o3.eq.1) then 
         trn_o3=exp(am_sat*tau_o3(ib))
      else 
         trn_o3=exp(am_sat*k_o3(ib)*to3_in)
      endif  

      if(iband.eq.2.or.iband.ge.5) trn_o3=1.0  
!c--------------------------------------------------
!c--------------------------------------------------
      return
      end
!c--------------------------------------------------


!c--------------------------------------------------
      subroutine get_co2_trn(iband,sza_in,vza_in,trn_co2,wvl_ib)

!c  Written by Myeong-Jae Jeong (MJ)
!c  Ver. 1.1 Aug 4, 2011
!c
      real*4 wvlarr(7), tau_co2(7)
      real*4 am_sat, pi, xdtor, trn_co2, wvl_ib, sza_in, vza_in
      integer ib, iband

      pi=3.141592
      xdtor=pi/180.0
!c    ib=iband-1  !; IDL
      ib=iband    !; FORTRAN

      data wvlarr /0.65,0.86,0.47,0.55,1.24,1.64,2.12/
      data tau_co2 /0.00000, 0.00000, 0.00000, 0.00000
     +            , 4.196e-4, 8.260e-3, 2.164e-2/

      am_sat=1.0/cos(xdtor*sza_in) + 1.0/cos(xdtor*vza_in)
      wvl_ib=wvlarr(ib)
      trn_co2=exp(am_sat*tau_co2(ib))

      if(iband.le.4) trn_co2=1.0   
!c--------------------------------------------------
!c--------------------------------------------------
      return
      end
!c--------------------------------------------------



