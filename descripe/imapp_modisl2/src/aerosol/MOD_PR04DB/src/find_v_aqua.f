! @TODO: search() can't handle data that is at the last node of the table
! TODO: Convert all of this to F95
      subroutine find_v_aqua(realbuf, outbuf, tmpvg, qa_flag)
      
      use landcover, only: get_landcover
      use calendars, only: gdatetime, season_from_doy, gregorian_from_doy
      use modis_surface, only: get_brdfcorr_sr,
     1                         get_aot500,
     1                         terrain_flag, 
     1                         terrain_flag_new,
     1                         get_LER412, 
     1                         get_LER470,
     1                         get_LER650, 
     1                         get_LER865
      use seawifs_surface_pressure, only: get_elevation
      
      real realbuf(26), outbuf(20)
c
c     Taking into account the terrain effect on Rayleigh
c
      include 'aottbl.inc'
      include 'newaottbl.inc'

c     ---input parameters
      real pdiff(3), tmp(14), xnvalm6(6), Dstar1, tmpvg(6)
c     ---intermediate parameters
      integer xlatp, xlonp, index_ii, index_ia, doy, ilat, ilon
      real nval(10,46,30), yy(10), yyw(8)
      real nnvalx(4,4,2,10), nnvalxw(4,4,2,8)
      real nvalxw(10,46,30,8), yy2(8), aa550
      real xnvalm(8)
      real terrain_flag_new5
      real trflg, pteran, x1, x2, xr470
      real aot_mod(6),r470new, xday, r470ss2
      real sza
      real  outbufvg(20)
      real  model_frac
      
      real            ::  toa_ndvi, tmp412, tmp470, ndvi_thold
      real						::	px_elev
      integer         ::  lc, gzflg, season, mod_sfc
      type(gdatetime) ::  gdt1
      integer         ::  status
      
      real            ::  r412_tbl, r470_tbl, r650_tbl
      real            ::  r412_135, r470_135, r650_135
      
      real    ::  tau_x412, tau_x412ss, tau_x412ss2, tau_x412_91
      real    ::  tau_x412ss_995, tau_x412ss2_995, tau_x412ss_96
      real    ::  tau_x412ss_97, tau_x412ss_94, tau_x412ss_95
      real    ::  tau_x412ss_98, tau_x412ss_91, tau_smoke, tau_x412ss2_98
      
      real    ::  tau_x412_new_91, tau_x412_new_93, tau_x412_new_94
      real    ::  tau_x412_new_96, tau_x412_new_995
      real    ::  tau_x470_new_91, tau_x470_new_92, tau_x470_new_93
      real    ::  tau_x470_new_94, tau_x470_new_95, tau_x470_new_96, tau_x470_new_995
      
      integer ::  tau_x412_flag, tau_x412ss_flag, tau_x412ss2_flag, tau_x412_91_flag
    
      real    ::  tau_x470, tau_x470ss, tau_x470ss2, tau_x470_new
      integer ::  tau_x470_flag2, tau_X470ss_flag, tau_x470_new_flag, tau_x470_new_91_flag
      real    ::  tau_x650
      integer ::  tau_x650_flag
      integer ::  tau_x412_flag2, tau_x470_flag, tau_x412_flag_91
      
      logical ::  abs_aero_flag
      
      logical ::  sr_fail_flag, do_veg
      integer ::  alg_flag, brdf_flag
      
      logical dflag, dflag2
      logical ::  debug, use_alternate_brdf
      integer ::  lprint
      
c     ---output parameters
      real xtau(3), alpha, ssa(3), tau550, sfc_typ
      integer qa_flag(4)
c     ---common parameter      
      common  /xday/ doy
			
			integer handle_lut_out_of_bounds
	    
c-----------------------------------------------------------------------
c      Initialization
c-----------------------------------------------------------------------

      mm = 10     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 30     ! rel. azimuth
      ma = 10     ! tau
      mw =  8     ! ssa

c-----------------------------------------------------------------------
c      Start processing the data
c-----------------------------------------------------------------------
      
      use_alternate_brdf = .false.
      
c     Load the input data into local storage

      xlat         = realbuf(1)
      xlong        = realbuf(2)
      sza          = realbuf(3)
      xthet        = realbuf(4)
      xphi         = realbuf(5)
      ref650       = realbuf(6)
      abs_aero_flag= .false.
      
      mod_sfc      = -999
      
      lprint = -999
      debug = .false.
c      if (abs(xlat-19.02752) < 0.001 .AND. abs(xlong+9.995494) < 0.001) then  
c        debug = .true.
c        lprint = 1
c      end if
      
      itmp = xlat
      xtmp = xlat - itmp
      xlatp = 10*(itmp+90) + 10*xtmp + 1
      itmp = xlong
      xtmp = xlong - itmp
      xlonp = 10*(itmp+180) + 10*xtmp + 1
      if (xlatp.gt.1800) xlatp = 1800
      if (xlonp.gt.3600) xlonp = 3600
      
!   -- convert geolocation into array indices.
c    	ilat = floor(xlat*10.0) + 900 + 1
c    	ilon = floor(xlong*10.0) + 1800 + 1
c    	if (ilat > 1800) ilat = 1800
c    	if (ilon > 3600) ilon = 3600
c			if (xlatp /= ilat .OR. xlonp /= ilon) print *, 'diff: ', xlat, xlong, xlatp, xlonp, ilat, ilon
			
c    	xlatp = ilat
c    	xlonp = ilon
    	
      trflg = terrain_flag(xlonp,xlatp)
      terrain_flag_new5 = terrain_flag_new(xlonp,xlatp)
			gzflg = terrain_flag_new5
			
      do i = 1, 6
         xnvalm6(i) = realbuf(i+6)
         ss = -1.*xnvalm6(i)/100.
         xnvalm6(i) = 10.**ss   
      enddo

      do i = 1, 14
         tmp(i) = realbuf(i+12)
      enddo
      
      stdv        = realbuf(13)
      toa_ndvi    = realbuf(14)

      Dstar1   = tmp(4)
      pdiff(1) = tmp(6)
      pdiff(2) = tmp(7)
      pdiff(3) = tmp(8)
      pteran   = tmp(9)
      xr470    = tmp(13)

!      if (Dstar1 > 1.06) go to 10
!      if (Dstar1 > 1.1) go to 10

      do i = 1, 3
      if (pdiff(i).lt.0.0.and.pdiff(i).gt.-1.E-4) 
     1    pdiff(i) = 0.0
!      if (trflg.lt.0.0) pdiff(i) = 0.0
      enddo
      if (xr470.lt.0.) pdiff(2) = 0.
      
!     -- get elevation of current pixel.
      px_elev = 0.0
      px_elev = get_elevation(xlat, xlong, status)
      if (status /= 0) then
        print *, "ERROR: Failed to get elevation for pixel, ", i, j, xlat, xlong, status
        return
      end if
!      print *, 'elev: ', xlat, xlong, px_elev


!     -- disable pressure corrections
!      pdiff(:) = 0.0
      
      x1 = sza
      x2 = xthet
      refl1= xnvalm6(5)                    ! 412 nm
      refl3= xnvalm6(3)                    ! 470 nm
      refl6= xnvalm6(1)                    ! 650 nm
      
!      refl1= xnvalm6(5) + pdiff(1)/1.30    ! 412 nm
!      refl3= xnvalm6(3) + pdiff(2)/1.30    ! 470 nm
!      refl6= xnvalm6(1)                    ! 650 nm
      rr412_mod = tmp(10)
      rr470_mod = tmp(11)
      band26    = tmp(5) 
      
!     -- apply pressure correction to TOA reflectance over Morocco and West Sudan.
      if (px_elev > 750.0 .AND. (xlat > 28.0 .AND. xlat < 37.0 .AND. xlong > -12.0 .AND. 
     1    xlong < 10.0)) then 
        refl3 = xnvalm6(3) + pdiff(2)/1.30    ! 470 nm
      end if
      
      if (px_elev > 900.0 .AND. (xlat > 10.5 .AND. xlat < 19.5 .AND. xlong > 20.5 .AND. 
     1    xlong < 29.0)) then
        refl3 = xnvalm6(3) + pdiff(2)/1.50    ! 470 nm
       end if
               
      if (refl6 > 0.2) go to 10
 
      rat = refl6 / refl1
      rat_650_470 = ref650 / rr470_mod
      rat_650_412 = ref650 / rr412_mod
      rat1 = rr470_mod / rr412_mod

      iflag = 1
			
			if (debug) then
        print *, 'find_v, in: ', xlat, xlong, sza, xthet, xphi, Dstar1
        print *, 'ler, 412, 470, 650: ', realbuf(22), realbuf(23), realbuf(6)
      endif
    
 5    continue
c---------------------------------------------------------------
c     Initialization
c---------------------------------------------------------------
c     --intermediate parameters
      dflag      = .false.
      dflag2     = .false.
      w0_x       = -999.
      w0_int     = -999.
      w0_int_470 = -999.
      w0_x412    = -999.
      w0_x470    = -999.
      
      tau_x412          = -999.
      tau_x412_91       = -999.
      tau_x412ss        = -999.
      tau_x412ss_995    = -999.0
      tau_x412ss2_995   = -999.0
      tau_x412ss2_98   = -999.0
      tau_x412ss_96     = -999.0
      tau_x412ss_97     = -999.0
      tau_x412ss_94     = -999.0
      tau_x412ss_95     = -999.0
      tau_x412ss_98     = -999.0
      tau_x412ss_91     = -999.0
      tau_smoke         = -999.0
      tau_x412_new_91   = -999.0
      tau_x412_new_93   = -999.0
      tau_x412_new_94   = -999.0
      tau_x412_new_96   = -999.0
      tau_x412_new_995  = -999.0
      tau_x470_new      = -999.0
      tau_x470_new_91   = -999.0
      tau_x470_new_92   = -999.0
      tau_x470_new_93   = -999.0
      tau_x470_new_94   = -999.0
      tau_x470_new_95   = -999.0
      tau_x470_new_96   = -999.0
      tau_x470_new_995  = -999.0
      
      tau_x470          = -999.
      tau_x470ss        = -999.
      tau_x470ss2       = -999.
      tau_x650          = -999.
      tau_x412_flag     = -999
      tau_x470_flag     = -999
      tau_x650_flag     = -999
      tau_x412_flag2    = -999
      tau_x470_flag2    = -999
      tau_x412_flag_91  = -999
      xxrat             = -999.
      xxrat2            = -999.
      aot               = -999.
      r412              = -999.0
      r412new           = -999. 
      r412ss            = -999.0
      r412ss2           = -999.0
      r470              = -999.0
      r470ss            = -999.0
      r470ss2           = -999.0
      r470new           = -999.
      r470_sav          = -999.0
      r650              = -999.0
      tau_x412_new      = -999.
      tau_x470_new      = -999.
      
      tau_x412ss2       = -999.
      tau_x470_new_91   = -999.
      sr_fail_flag      = .false.
      
c     -- output parameters
      tau550     = -999.
      alpha      = -999.
      sfc_typ    = -999.
      alg_flag   = -999
      
      do i = 1, 3
         xtau(i) = -999.
         ssa(i) = -999.
      enddo

      do i = 1, 4
         qa_flag(i) = 0
      enddo

      do i = 1, 6
         aot_mod(i) = -999.
      enddo

      do i = 1,20
         outbuf(i) = -999.
      enddo
      
      if (sza.gt.72.0) go to 10
      
      xday = real(doy)

c---------------------------------------------------------------
c     Screen for pixels outside reasonable ranges of reflectance
c---------------------------------------------------------------
c      if (refl1.gt.0.0.and.refl1.lt.0.09.and.
c     1    refl6.gt.0.0.and.refl6.lt.0.14) go to 11
c      if (refl1.gt.0.09.and.refl1.lt.0.50.and.
c     1    res.gt.6.0) go to 11
c      go to 10

11    continue
      if (xphi.gt.179.99) go to 10
      if (xphi.lt.6.0) xphi = 6.
      
c     -- sun glint mask
c
      cc     = 3.14159/180.
      psi    = acos(cos(sza*cc)*cos(xthet*cc) +
     1         sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      glint_ang = psi/cc

c      if (abs(psi/cc).lt.35.0) go to 10
c
c     -- scattering angle (scat_ang)
c
      cc     = 3.14159/180.
      psi    = acos(cos(sza*cc)*cos(xthet*cc) -
     1         sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      scat_ang = 180. - psi/cc
 
c      if (scat_ang .gt. 175.) go to 10
      
c--------------------------------------------------------
c   Load Surface Reflectance
c--------------------------------------------------------
c     -- get base surf. reflc. values from surf. coeff. tables and save.
      r412_tbl = get_LER412(xlatp, xlonp, toa_ndvi, scat_ang, xphi)
      r470_tbl = get_LER470(xlatp, xlonp, toa_ndvi, scat_ang, xphi)
      r650_tbl = get_LER650(xlatp, xlonp, toa_ndvi, scat_ang, xphi)
      r865_tbl = get_LER865(xlatp, xlonp)
      
      r412_135 = get_LER412(xlatp, xlonp, 0.1, 135.0, 95.0)
      r470_135 = get_LER470(xlatp, xlonp, 0.1, 135.0, 95.0)
      r650_135 = get_LER650(xlatp, xlonp, 0.1, 135.0, 95.0)
      
c     -- set to surface reflectance values to default table values.
      r412 = r412_tbl
      r470 = r470_tbl
      r650 = r650_tbl
      r865 = r865_tbl

      if (r865.lt.12.0.and.glint_ang.lt.30.0) go to 10    ! sun glint mask

      if (debug) print *, 'glint, scat. ang, r865: ', glint_ang, scat_ang, r865
      if (debug) print *, 'xlatp, xlonp, toa_ndvi, xphi: ', xlatp, xlonp, toa_ndvi, xphi
      if (debug) print *, "r412, r470, r650: ", r412, r470, r650

c     -- out of scope, skip.      
c      if (r412.gt.30.0) go to 10
c      if (r470.gt.50.0) go to 10
c      if (r650.gt.60.0) go to 10
      
!     -- initialize other surface variables
      r470ss    = r470_tbl
      r470new   = r470_tbl
      
      r412ss    = r412_tbl
      r412ss2   = r412_tbl
      r412new   = r412_tbl
      
c-- 470 nm
c     Backward direction
c      if (xphi.ge.90.) then
c        r412 = xsfc412bwd(xlonp,xlatp)
c        r470 = xsfc470bwd(xlonp,xlatp)
c        r650 = xsfc650bwd(xlonp,xlatp)
c      endif
c      
cc     Forward direction
c      if (xphi.lt.90.0) then
c        r412 = xsfc412fwd(xlonp,xlatp)
c        r470 = xsfc470fwd(xlonp,xlatp)
c        r650 = xsfc650fwd(xlonp,xlatp)
c      end if
      
!---------------------------------------------------------------------------------------------------
! START NEW AERONET-DERIVED SURFACE REFLECTIVITY
!---------------------------------------------------------------------------------------------------
!     -- get landcover, then surface reflectivity for pixel..  
      lc = get_landcover(xlat, xlong, status)
      if (status .ne. 0) then
        print *, "ERROR: Failed to get land cover value: ", i, j, xlat, xlong, status
        return
      end if
      
!     -- turn off tropical Sahel transition zone (zone 27) if it's not winter.
      if (gzflg == 27 .AND. (xday >= 60 .AND. xday < 335)) gzflg = 5

      brdf_flag = -1
      tmp412=-999.0 ; tmp470=-999.0 ; tmp670=-999.0
      season = season_from_doy(2005, doy)   ! year isn't currently available in here.
      gdt1 = gregorian_from_doy(2005,doy)   ! just use 2005 for now.
      brdf_flag = get_brdfcorr_sr(xlat, xlong, xphi, scat_ang, px_elev, gdt1%month, 
     &						toa_ndvi, stdv, gzflg, lc, tmp412, tmp470, tmp670, use_alternate_brdf=use_alternate_brdf, 
     &            debug=debug)
      if (brdf_flag == 0) then
        r412 = tmp412
        r470 = tmp470
      end if
      
!     -- insert special function to calculate surface reflectivity over N.Africa coast 
!     --  e.g. Blida, Saada -- geozone 2.  Overriding values from BRDF in modis_surface.f95.
!      if (gzflg == 2) then
!       mod_sfc = 8
!       call rsfc470(dflag2,mod_sfc,xday,xthet,xphi,scat_ang,r470)
!       if (dflag2) go to 10
!      endif

!     -- special function surface reflectivity for Hamim and all of Arabian Peninsula.
      if (gzflg >= 6 .AND. gzflg <= 11) then 
        if (r650_135 > 32.0 .AND. (r650_135/r412_135) > 3.7) then
          
!         -- r412_135+0.5 == value from table was a bit too low here. Adjust upwards to
!         -- bring the AOT back down.
          mod_sfc = 9
          call newsfc412_arab(dflag2,mod_sfc,xday,xlatp,xlonp,xthet,xphi,
     &                        scat_ang,terrain_flag_new5,r412_135,r412new)
!          if (xday >= 152 .AND. xday <= 243) then     ! summer
            if ((r412_135 < 7.25) .OR. r650_135 > 36.0 .OR. use_alternate_brdf) then
              continue
            else
              r412 = r412new
            end if
!          else
!            r412 = r412new
!          end if 
        end if
      end if

      
!     -- use surface table value as last resort.  Skip if undefined.
      if (r412 < -900.0) r412 = r412_tbl
      if (r470 < -900.0) r470 = r470_tbl
      if (r650 < -900.0) r650 = r650_tbl
      
      if (r412 < -900.0 .OR. r470 < -900.0 .OR. r650 < -900.0) then
        if (toa_ndvi >= 0.2) then
          sr_fail_flag = .true.
          go to 637               ! go to vege. retrieval, no need for SR database values.
        else
          go to 10
        end if
      end if
      
      r650_sav = r650
      r470_sav = r470
      r412_sav = r412
      
c     --- DIRECT COPY/PASTE FROM SEAWIFS CODE
!-- 470 nm

      r470new = r470_sav

!     -- start tweaking r412ss, r412ss2, r470ss---------------------------------------
!     -- new bi-directional surface 2
      ddx = 0.0
      if (xday > 150.0 .and. xday < 258.0)  ddx = 1.0   ! Jun, Jul

      ddx2 = 0.0
      if (xday > 150.0 .and. xday < 258.0)  ddx2 =  0.5 ! Jun, Jul

      r470ss  = r470 - ddx     
      r470ss2 = r470      

      dda1 = r470_tbl*0.1 +0.5
      if (xday > 150.0 .and. xday < 258.0)  dda1 = r470_tbl*0.1  ! Jun, Jul
      dda2 = r470_tbl*0.0
      if (r470_tbl > 12.0)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl+dda1+dda2*(scat_ang-175.)/5.
      r470ss = r470ss  - ddx - ddx2
      endif

!     Libya - Egypt 1
      if (xlat >21. .and. xlong > 12.0) then
      dda1 = r470_tbl*0.08 + 0.3
      if (r412_tbl > 12.0) dda1 = r470_tbl*0.08 +1.8
      dda2 = r470_tbl*0.0
      if (r470_tbl > 12.0)  then
       if (scat_ang >= 130.0 .and. scat_ang < 175.0)  
     1   r470ss = r470_tbl + dda1*(scat_ang-130.)/45.
       if (scat_ang >= 175.0)  
     1   r470ss = r470_tbl+dda1+dda2*(scat_ang-175.)/5.
      endif
      endif
 
!     Libya - Egypt 2
      if (xlat >15. .and. xlong > 22.0) then
      dda1 = r470_tbl*0.08 + 0.3
      if (r412_tbl > 12.0) dda1 = r470_tbl*0.08 +1.8
      dda2 = r470_tbl*0.0
      if (r470_tbl > 12.0)  then
       if (scat_ang >= 130.0 .and. scat_ang < 175.0)  
     1   r470ss = r470_tbl + dda1*(scat_ang-130.)/45.
       if (scat_ang >= 175.0)  
     1   r470ss = r470_tbl+dda1+dda2*(scat_ang-175.)/5.
      endif
      endif

!     N. Algeria 1
      if (xlat >31.5 .and. xlong > 4. .and. xlong < 10.) then
      dda1 = r470_tbl*0.08 + 0.5
      dda2 = r470_tbl*0.05 + 0.5
      if (r470_tbl > 15.0)  then
       if (scat_ang >= 100.0 .and. scat_ang < 170.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/70.
       if (scat_ang >= 170.0)
     1   r470ss = r470_tbl+dda1+dda2*(scat_ang-170.)/10.
      r470ss = r470ss  - ddx - ddx2
      endif
      endif

!     NE Mauritania 1
      if (xlat>22.5 .and. xlat<30.0.and. xlong>-11.0.and. xlong< -5.) then
      if (r412_tbl > 9.)  then
        r470ss = r470ss + 0.9
      endif
      endif
 
!     NE Mauritania 2
      if (xlat>22.5 .and. xlat<26.0.and. xlong>-13.0.and. xlong< -11.001) then
      if (r412_tbl > 9.)  then
        r470ss = r470ss + 0.9
      endif
      endif

!     NE Mauritania 3
      if (xlat>20.0 .and. xlat<30.0.and. xlong>-20.0.and. xlong< -12.0) then
      if (r412_tbl > 9.)  then
        r470ss = r470ss + 0.9
      endif
      endif

      if (gzflg >= 6 .and. gzflg <= 11) then

      ddx = 0.0
      if (xday > 90.0 .and. xday < 121.0)   ddx = 0.3  ! Apr
      if (xday > 120.0 .and. xday < 151.0)   ddx = 0.5  ! May
      if (xday > 150.0 .and. xday < 258.0)  ddx =  1.0 ! Jun, Jul, Aug
 
      ddx2 = 0.0
      if (xday > 150.0 .and. xday < 258.0)  ddx2 =  0.5 ! Jun, Jul, Aug

      r470ss  = r470 - ddx
      r470ss2 = r470
 
      dda1 = r470_tbl*0.1 +0.5
      if (xday > 150.0 .and. xday < 258.0)  dda1 = r470_tbl*0.1  ! Jun, Jul
      dda2 = r470_tbl*0.0
      if (r470_tbl > 12.0)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl+dda1+dda2*(scat_ang-175.)/5.
      r470ss = r470ss  - ddx - ddx2
      endif

      dda1 = r470_tbl*0.1 +1.0
      if (xday > 90.0 .and. xday < 121.0)  dda1 = r470_tbl*0.1  ! Apr
      dda2 = r470_tbl*0.0
      if (r412_tbl > 11.0) dda1 = r470_tbl*0.1 +3.0
      if (r412_tbl > 12.0) dda1 = r470_tbl*0.1 +4.0
      if (xday > 120.0 .and. xday < 152.0 .and. r412_tbl > 11.0)  dda1 = dda1 - 0.7  ! May
      if (xday > 151.0 .and. xday < 182.0 .and. r412_tbl > 11.0)  dda1 = dda1 + 0.4  ! Jun
       if (r412_tbl > 9.0)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl+dda1+dda2*(scat_ang-175.)/5.
      r470ss = r470ss  - ddx - ddx2
       endif

      if (xday > 150.0 .and. xday < 244.0) then       ! Jun,Jul, Aug
          if (r412_tbl > 13.0 .and. rat1 > 1.43) r470ss = r470_tbl
          if (Dstar1 > 1.1)  r470ss = r470_tbl
      endif

!     Southern Arabian Pen
      if (xday > 152.0 .and. xday < 244.0) then   ! Jun, Jul, Aug

      ddx  = 0.0
      if (xday > 152.0 .and. xday < 182.) ddx  = 0.5    ! Jun
      if (xday > 181.0 .and. xday < 244.0) ddx  = 1.5   ! Jul, Aug
      dda1 = 1.5
      dda2 = 0.0
      if (r412_tbl > 11.0) dda1 = 1.8
      if (r412_tbl > 12.0) dda1 = 2.0

      if (xlat < 18.0.and. xlong <= 49.0) then
      if (r412_tbl > 9.)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl-ddx + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl-ddx +dda1+dda2*(scat_ang-175.)/5.
      endif
      endif

      if (xlat < 18.5 .and. xlong > 49.0.and. xlong <= 52.5) then
      if (r412_tbl > 9.)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl-ddx + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl-ddx +dda1+dda2*(scat_ang-175.)/5.
      endif
      endif

      if (xlat < 23.0 .and. xlong > 52.5) then
      if (r412_tbl > 9.)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl-ddx + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl-ddx +dda1+dda2*(scat_ang-175.)/5.
      endif
      endif

      endif      ! end of southern Arabian Pen

      endif

      if (r470 >= 24.0)   r470 = 23.9
      if (r470 < 1.0)     r470 = 1.0
      if (r470ss >= 24.0) r470ss = 23.9
      if (r470ss < 1.0 .AND. (gzflg <= 11 .OR. gzflg == 27)) go to 10 !r470ss = 1.0

      if (lprint > 0) print *,'final r470,r470ss =', r470,r470ss

!-- 412 nm

      ddx = 0.0
      if (xday > 32.0 .and. xday < 60.0)  ddx = 0.5            ! Feb
      if (xday > 59.0 .and. xday < 152.0)   ddx = 1.0         ! Mar,Apr, May
      if (xday > 151.0 .and. xday < 213.0)  ddx = 1.0 + 0.8    ! Jun, Jul
      if (xday > 212.0 .and. xday < 244.0)  ddx = 0.5          ! Aug
      if (xday > 243.0 .and. xday < 335.0)  ddx = 0.           ! Sep, Oct, Nov

      ddx1 = 0.0
      if (xday < 32.0)   ddx1 = 0.0                         ! Jan
      if (xday > 59.0 .and. xday < 121.0)   ddx1 = 0.0      ! Mar,Apr
      if (xday > 120.0 .and. xday < 244.0)  ddx1 = 0.0      ! May, Jun,Jul,Aug
      if (xday > 243.0 .and. xday < 305.0)  ddx1 = 0.0      ! Sep, Oct
      if (xday > 304.0)  ddx1 = 0.0                         ! Nov, Dec

      ddx2 = 0.0
      if (xday < 32.0)   ddx2 = 0.0                         ! Jan
      if (xday > 59.0 .and. xday < 152.0)   ddx2 = 0.5      ! Mar,Apr,May
      if (xday > 243.0 .and. xday < 305.0)  ddx2 = 0.0      ! Sep, Oct
      if (xday > 304.0)  ddx2 = 0.0                         ! Nov, Dec

      ddx3 = 0.3
      if (xday > 0.0 .and. xday < 152.0)   ddx3 = 1.2       ! Jan,Feb,Mar,Apr,May
      ddx33 = 0.4
      if (xday > 32.0 .and. xday < 60.0)    ddx33 = 0.7     ! Feb
      if (xday > 59.0 .and. xday < 152.0)   ddx33 = 2.5     ! Mar,Apr, May
      if (xday > 151.0 .and. xday < 244.0)  ddx33 = 0.7     ! Jun,Jul,Aug
      if (xday > 243.0 .and. xday < 305.0)  ddx33 = 0.4     ! Sep, Oct
      if (xday > 304.0) ddx33 = 0.6                         ! Nov, Dec

      r412ss2 = r412_tbl - ddx

      if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1   r412ss2 = r412_tbl- ddx + ddx1*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1   r412ss2 = r412_tbl- ddx+ddx1+0.2*(scat_ang-170.)/10.

      if (r412_tbl > 9.0)  then
       if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1   r412ss2 = r412_tbl + ddx2*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1     r412ss2 = r412_tbl+ddx2+0.0*(scat_ang-170.)/10.
      endif

      if (r412_tbl > 11.0)  then
       if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1   r412ss2 = r412_tbl + ddx3*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1   r412ss2 = r412_tbl+ddx3+0.0*(scat_ang-170.)/10.
      endif

      if (r412_tbl > 12.0)  then
       if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1    r412ss2 = r412_tbl + ddx33*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1    r412ss2 = r412_tbl+ddx33+0.0*(scat_ang-170.)/10.
      endif
            
!     -- new bi-directional surface 2

      ddx = 0.0
      if (xday > 0.0 .and. xday < 32.0)   ddx = 0.5 + 0.7    ! Jan
      if (xday > 31.0 .and. xday < 60.0)  ddx = 0.5 + 0.4    ! Feb
      if (xday > 59.0 .and. xday < 152.0)  ddx = 0.5 + 0.4    ! Mar, Apr, May
      if (xday > 151.0 .and. xday < 244.0)  ddx = 0.5 + 0.4    ! Jun, Jul, Aug
      if (xday > 243.0 .and. xday < 305.0)  ddx = 0.5 + 0.7     ! Sep, Oct
      if (xday > 304.0 .and. xday < 335.0)  ddx = 0.5 + 0.4     ! Nov
      if (xday > 334.0)                     ddx = 0.5 + 0.4     ! Dec
      ddx2 = 0.0
      if (xday > 0.0 .and. xday < 122.0)  ddx2 = 1.5    ! Mar, Apr
      if (xday > 121.0 .and. xday < 213.0)  ddx2 = 1.5    ! May,Jun
      if (xday > 212.0 .and. xday < 305.0)  ddx2 = 1.5    ! Aug, Sep, Oct
      if (xday > 304.0 .and. xday < 335.0)  ddx2 = 0.     ! Nov
      ddx3 = 1.5
      if (xday < 32.0)  ddx3 = 1.    
 
      r412ss = r412_tbl - ddx

      if (xday < 152.0 .or. xday > 244.0)  then   ! angular adjustment 

      if (scat_ang >= 100.0 .and. scat_ang < 175.0)       
     1   r412ss = r412_tbl- ddx + 1.0*(scat_ang-100.)/75. 
      if (scat_ang >= 175.0)  
     1   r412ss = r412_tbl- ddx+1.0+0.3*(scat_ang-175.)/5.

      if (xday > 59.0 .and. xday < 152.0)  then   ! spring
       if (xlat > 25.) then
       if (scat_ang >= 140.0 .and. scat_ang < 170.0) 
     1     r412ss = r412_tbl- ddx + 1.0*(scat_ang-140.)/30.
       if (scat_ang >= 170.0)  
     1     r412ss = r412_tbl- ddx+1.0+0.3*(scat_ang-170.)/10.
       endif
      endif

      if (r412_tbl > 9.5)  then
       if (scat_ang >= 150.0 .and. scat_ang < 175.0) 
     1    r412ss = r412_tbl- ddx + ddx3*(scat_ang-150.)/25.
       if (scat_ang >= 175.0)  
     1    r412ss = r412_tbl- ddx+ddx3+0.5*(scat_ang-175.)/5.
      endif

      endif     ! end of angular adjustment

!     E. Central Algeria, long narrow desert
      dda1 = 2.5                                       ! winter, spring,fall
      dda2 = 9.5
      dda3 = 120.0
      dda4 = 0.5
      dda5 = 55.
       if (xday > 151.0 .and. xday < 244.0)  then       ! summer
        dda1 = 1.1
        dda2 = 9.5
        dda3 = 100.0
        dda4 = 0.0
        dda5 = 75.
       endif
       if (xlat>27. .and. xlat<36. .and. xlong>2.5 .and. xlong<11.5) then
       if (r412_tbl > dda2)  then
        if (scat_ang >= dda3 .and. scat_ang < 175.0) 
     1    r412ss = r412_tbl- ddx + dda1*(scat_ang-dda3)/dda5
        if (scat_ang >= 175.0)  
     1    r412ss = r412_tbl- ddx+dda1+dda4*(scat_ang-175.)/5.
       endif
       endif

!     Libya - Egypt 1
      if (xlat >20. .and. xlong > 14.9) then
      dda1 = r412_tbl*0.08 +0.5
      if (r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 0.8
      if (xday > 59.0 .and. xday < 91.0 .and.r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 1.5   ! Mar
      if (xday > 151.0 .and. xday < 244.0 .and.r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 0.4   ! summer
      if (xday > 151.0 .and. xday < 244.0)  dda1 =  dda1 - 0.5      ! summer
      dda2 = r412_tbl*0.0
      if (r412_tbl > 9.6)  then
       if (scat_ang >= 120.0 .and. scat_ang < 175.0)  
     1   r412ss = r412_tbl- ddx + dda1*(scat_ang-120.)/55.
       if (scat_ang >= 175.0)  
     1   r412ss = r412_tbl- ddx+dda1+dda2*(scat_ang-175.)/5.
      endif
      endif

!     Libya - Egypt 2
      if (xlat >15. .and. xlong > 22.0) then
      dda1 = r412_tbl*0.08 +0.5
      if (r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 0.8
      if (xday > 59.0 .and. xday < 91.0 .and.r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 1.5   ! Mar
      if (xday > 151.0 .and. xday < 244.0 .and.r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 0.4   ! summer
      if (xday > 151.0 .and. xday < 244.0)  dda1 =  dda1 - 0.5      ! summer
      dda2 = r412_tbl*0.0
      if (r412_tbl > 9.6)  then
       if (scat_ang >= 120.0 .and. scat_ang < 175.0)  
     1   r412ss = r412_tbl- ddx + dda1*(scat_ang-120.)/55.
       if (scat_ang >= 175.0)  
     1   r412ss = r412_tbl- ddx+dda1+dda2*(scat_ang-175.)/5.
      endif
      endif
 
      if (xday > 151.0 .and. xday < 335.0)  go to 33   ! winter, spring
!     N. Algeria 1
      if (xlat >31.5 .and. xlong > 4. .and. xlong < 10.) then
      dda1 = r412_tbl*0.08 + 0.5
      dda2 = r412_tbl*0.05 + 0.5
      if (r412_tbl > 9.4)  then
       if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1   r412ss = r412_tbl + dda1*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1   r412ss = r412_tbl+dda1+dda2*(scat_ang-170.)/10.
      endif
      endif
33    continue

!     N. Algeria 2
      if (xlat >30.0 .and. xlat <32.0 .and. xlong > 4.8 .and. xlong < 7.) then
      dda1 = r412_tbl*0.08
      dda2 = r412_tbl*0.05
      if (r412_tbl > 9.0)  then
        if (scat_ang >= 150.0 .and. scat_ang < 175.0) 
     1    r412ss = r412_tbl + dda1*(scat_ang-150.)/25.
      if (scat_ang >= 175.0)  
     1    r412ss = r412_tbl+dda1+dda2*(scat_ang-175.)/5.
      endif
      endif

!     Chad - Libya border
      if (xday > 31.0 .and. xday < 60.0)  go to 36  ! use for all months except for Feb
      if (xlat >20. .and. xlat <25. .and. xlong >15.0.and. xlong <17.5) then
      if (r412_tbl > 12.0)  then
        r412ss = r412ss + 0.5
      endif
      endif
36    continue

!     NE Mauritania 1
      if (xlat>22.5 .and. xlat<30.0.and. xlong>-11.0.and. xlong< -5.) then
      if (r412_tbl > 9.)  then
        r412ss = r412ss + 0.5
      endif
      endif

!     NE Mauritania 2
      if (xlat>22.5 .and. xlat<26.0.and. xlong>-13.0.and. xlong< -11.001) then
      if (r412_tbl > 9.)  then
        r412ss = r412ss + 0.5
      endif
      endif

!     Morocco 1
      if (xlat>30.0 .and. xlat<36.0 .and. xlong<= -7.5) then
      if (r412_tbl > 5.8)  then
        if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1    r412ss = r412_tbl + 4.0*(scat_ang-100.)/75.
      if (scat_ang >= 175.0)
     1    r412ss = r412_tbl+4.0+0.0*(scat_ang-175.)/5.
      endif
      endif

!     Morocco 2
      if (xlat>31.5 .and. xlat<36.0 .and. xlong<= -5.0.and. xlong> -7.5) then
      if (r412_tbl > 5.8)  then
        if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1    r412ss = r412_tbl + 4.0*(scat_ang-100.)/75.
      if (scat_ang >= 175.0)
     1    r412ss = r412_tbl+4.0+0.0*(scat_ang-175.)/5.
      endif
      endif


      if (xday > 151.0 .and. xday < 335.0)  go to 37   ! winter, spring

!     Central Algeria
      if (xlat >22.5 .and. xlat <30. .and. xlong > -4.0.and. xlong < 1.) then
      if (r412_tbl > 9.)  then
        r412ss = r412ss + 0.5
      endif
      endif

!     Mali - Algeria border
      if (xlat>20. .and. xlat<22.501.and. xlong>-2.0.and. xlong< 2.) then
      if (r412_tbl > 9.)  then
        r412ss = r412ss + 0.5
      endif
      endif

37    continue

      if (lprint > 0) print *,'scat_ang,b,f r412ss =',r412_tbl,scat_ang,r412ss

      if (r412ss2 > 20.0) r412ss2 = 19.9
      if (r412ss2 < 1.0) r412ss2 = 1.0      
      if (r412 > 20.0 .and. r412 < 40.0) r412 = 19.9
      if (r412 > 40.0) go to 10
      if (r412 < 1.0)  r412 = 1.0
      if (r412ss >= 20.0 .and. r412ss < 40.0) r412ss = 19.9
      if (r412ss < 1.0)  r412ss = 1.0
      
      if (r470 >= 24.0)   r470 = 23.9
      if (r470 < 1.0)     r470 = 1.0
      if (r470ss >= 24.0) r470ss = 23.9
      if (r470ss < 1.0)   r470ss = 1.0
      
      if (lprint > 0) print *,'final r412,r412ss =',r412,r412ss
!      write (888,'(3(F10.5,1X))') xlat, xlong, r412ss
      
      if (r650 >= 47.0) r650 = 46.9
      if (r650 < 1.0)   r650 = 1.0

C     -- END DIRECT COPY SEAWIFS CODE

      if (r412.gt.12.0.and.r412.lt.80.) qa_flag(4) = 2
            
c--------------------------------------------------------
c       Cloud Screening
c--------------------------------------------------------
c-- band26 values should be fill value (-999) so most of these
c-- test are non-functional. Confirmed w/ Dr. Hsu 2013-03-01 --CB.
      if (debug) print *,'--- start cloud screening ----'
      if (debug) print *, 'band26, xnvalm6(1), xvnalm6(5): ', band26, xnvalm6(1), xnvalm6(5)
      if (band26.lt.0.0.and.xnvalm6(5).lt.0.0.and.
     1    xnvalm6(1).gt.0.0) go to 620
      if (band26.gt.0.0.and.xnvalm6(5).lt.0.0.and.
     1    xnvalm6(1).gt.0.0) go to 10
      rat = ref650/ rr412_mod
      rat1 = rr470_mod / rr412_mod
      if (debug) print *, 'rat, rat1, ref650: ', rat, rat1, ref650
      if (ref650.gt.45.0.and.rat.lt.1.4) go to 10
      if (ref650.gt.56.0.and.rat.lt.1.3) go to 10
      
      if (debug) print *, 'trflg, rr412_mod, rat1, r412: ', trflg, rr412_mod, rat1, r412
      if (trflg.gt.0.0.and.rr412_mod.gt.10.0.and.
     1    rat1.gt.1.25) go to 50

      if (r412.gt.6.0) then
        if (band26.gt.0.0.and.rr412_mod.gt.6.0) go to 10
      else
        if (band26.gt.0.0.and.rr412_mod.gt.10.0.and.rat1.lt.1.25)
     1      go to 10
      endif
      if (debug) print *, '--- end cloud screening ---'
      
50    continue

c--------------------------------------------------------
c   Special Case: For thick strong-absorbing dust plume, 
c                 go to 3-channels
c--------------------------------------------------------

      dd = r650 / r412
      if (r650.lt.30.0.and.r650.gt.9.0.and.dd.gt.3.2) go to 80
      if (trflg.lt.0.0) go to 80
      if (rat1.gt.1.40.and.ref650.gt.30.0) then
        if (r650.gt.27.0.and.refl6.lt.0.09) go to 80
        if (r650.gt.25.0) go to 80
        xxrat = 0.8

c       if (rat1 > 1.5) go to 10
c       if (rat1 > 1.1) go to 10
      
        abs_aero_flag = .true.    ! -- set flag for use below.
!      print *, 'going to 620: ', xlat, xlong, r650, dd 
        go to 620
      endif
 80   continue
      
!     -- vegetated area in zone 1 next to DMN_Maine_Soroa site in summer. Tweak.
      if (gzflg == 1 .AND. (gdt1%month >= 6 .AND. gdt1%month <= 8)) then
        if (r650_135 < 16.0) then
          if (toa_ndvi < 0.18) then
            r412 = min(r412 + 2.0, 20.0)
            r470new = min(r470new + 2.0, 24.0)
            r650 = min(r650 + 2.0, 47.0)
          else
            r412 = max(r412 - 2.0, 1.0)
            r470new = max(r470new - 2.0, 1.0)
            r650 = max(r650 - 2.0, 1.0)
          end if
        end if
      end if
      
c--------------------------------------------------------
c     Preliminary Retrieval on AOT
c--------------------------------------------------------

c
c     For Moderate AOT, Use 412 - 470 nm Pair
c
c     Retrieving 470 nm AOT
c
      if (debug) print *, '--- starting aot470 --- '
      refl = refl3
      x3 = xphi

      tau_x470_flag = -999
      tau_x470_flag2 = -999
      
!     -- retrieve using most aerosol models -- in preparation for call to get_aot500.
      tau_x470_new_91 = -999.0 ; tau_x470_new_92 = -999.0 ; tau_x470_new_93 = -999.0
      tau_x470_new_94 = -999.0 ; tau_x470_new_95 = -999.0 ; tau_x470_new_96 = -999.0 
      tau_x470_new_995 = -999.0
      
      imod = 3                                ! w0 = 0.96
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1             imod,r470,tau_x470,tau_x470_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag, tau_x470)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
            
!      if (gzflg < 0 .or. gzflg > 11) go to 81    ! gzflg > 11 everywhere except N. Africa/SAP.

      imod = 3                                ! w0 = 0.96
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1              imod,r470ss,tau_x470ss,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470ss)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      if (tau_x470.lt.0.0601.and.Dstar1.gt.1.1.and.rat1.gt.1.6)
     1    then
      imod = 3                                ! w0 = 0.96
      r470ss = r470ss - 1.0
      if (r470ss .lt. 1.0) r470ss = 1.0
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1              imod,r470ss,tau_x470ss,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      go to 313
      endif
 
      if (tau_x470.lt.0.5.and.Dstar1.gt.0.98.and.rat1.gt.1.46)
     1    then
      imod = 3                                ! w0 = 0.96
      r470ss = r470ss - 1.0
      if (Dstar1.gt.1.04) r470ss = r470ss - 0.5
      if (r470ss .lt. 1.0) r470ss = 1.0
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1              imod,r470ss,tau_x470ss2,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      endif
313   continue
      
      if (r470new >= 24.0) go to 81
			if (r470new .lt. 1.0) r470new = 1.0
			
      imod = 3                                ! w0 = 0.96
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1              imod,r470new,tau_x470_new,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      tau_x470_new_96 = tau_x470_new
      
      imod = 4                                ! w0=0.995
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_995,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_995)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      
      imod = 2      
      model_frac = 0.5                        ! w0=0.95
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_95,tau_x470_flag2,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_95)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 2                                ! w0=0.94
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_94,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_94)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 1          
      model_frac = 2.0/3.0                    ! w0=0.93
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_93,tau_x470_flag2,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_93)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 1
      model_frac = 1.0/3.0                    ! w0=0.92
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_92,tau_x470_flag2,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_92)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod =  1                               ! w0 = 0.91
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1             imod,r470new,tau_x470_new_91,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_91)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
 
      if (r470 < 12.0 .and. rr470_mod > 11.0) then
        rat_470_412 = rr470_mod / rr412_mod
        if (rat_470_412 > 1.85) then
          imod = 1                              ! w0 = 0.91
          call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,                 
     1                  imod,r470new,tau_x470_new_91,tau_x470_flag2,trflg,0.0,debug) 
          if (dflag) go to 10                                                 
          status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_91)
          if (status /= 0) then
            print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
            return
          end if
        endif
      endif
 
      if (xthet > 62.0) then
        rat_470_412 = rr470_mod / rr412_mod
        if (rat_470_412 > 1.7) then
          imod = 1                              ! w0 = 0.91
          call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,             
     1                  imod,r470new,tau_x470_new_91,tau_x470_flag2,trflg,0.0,debug)
          if (dflag) go to 10
          status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_91) 
          if (status /= 0) then
            print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
            return
          end if
        end if
      end if
      if (debug) print *, '--- end aot470 ---'
81    continue

C      if (lprint > 0) print *,'tau_x470,tau_x470ss,tau_x470_new=',
C     1    tau_x470,tau_x470ss,tau_x470_new
     

c--------------------------------------------------------
c          Retrieving 412 nm AOT
c--------------------------------------------------------
 
      refl = refl1
      x3 = xphi

      tau_x412_flag = -999
      tau_x412_flag2 = -999
      
      r412new = r412                         ! transitional zone
      
      imod = 5                                ! w0 = 0.94
      call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,   
     1              imod,r412,tau_x412,tau_x412_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag2, tau_x412)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
!---
      if (r412 > 11.0 .and. tau_x412 < 0.4) then
        tau_x412_91 = tau_x412 * 2.
        go to 630
      else if (r412 > 12.0) then
        tau_x412_91 = tau_x412 * 2.
        go to 630
      endif

      refl = refl1
      x3 = xphi
 
      tau_x412_flag_91 = -999

      imod = 4                               ! w0 = 0.91
      call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1             imod,r412,tau_x412_91,tau_x412_flag_91,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag_91, tau_x412_91)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
630   continue
      
      tau_x412_flag = -999
      tau_x412_new_91 = -999.0 ; tau_x412_new_93 = -999.0 ; tau_x412_new_94 = -999.0
      tau_x412_new_96 = -999.0 ; tau_x412_new_995 = -999.0
      
      tau_x412ss_995 = -999.0 ; tau_x412ss2_995 = -999.0 ; tau_x412ss_96 = -999.0
      tau_x412ss_97 = -999.0 ; tau_x412ss_94 = -999.0 ; tau_x412ss_95 = -999.0
      tau_x412ss_98 = -999.0
    
      if (r412new < 1.0) go to 631 
      if (r412new >= 20.0) go to 631
      
      imod = 8      ! w0=0.995
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_995,tau_x412_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_995)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
!      print *
!      print *, 'aero_412, tau_x412_new_995: ', refl, sza, xthet, xphi, imod, r412new, tau_x412_new_995, tau_x412_flag, trflg, dflag
!      print *
      
      imod = 6      ! w0=0.96
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_96,tau_x412_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_96)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 5      ! w0=0.94
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_94,tau_x412_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_94)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 4      ! w0=0.93
      model_frac = 2.0/3.0
      if (dflag) print *,'calling aero_412, model_frac: imod, model_frac: ', imod, model_frac
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_93,tau_x412_flag,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_93)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 4      ! w0=0.91
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_91,tau_x412_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_91)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
            
!     -- arabian peninsula
      if (gzflg >= 6 .and. gzflg <= 11) then
        tau_x412_new = tau_x412_new_94
      end if

      if (gzflg == 6) then
        if (xday >= 60.0 .and. xday < 274.0) go to 602
        tau_x470_new = tau_x412_new
      end if
      
602   continue

      if (gzflg == 8) then
        if (xday >= 182.0 .and. xday < 274.0) go to 603
        tau_x470_new = tau_x412_new
      end if
603   continue

      if (gzflg == 7) tau_x470_new = tau_x412_new

      if (gzflg == 9) tau_x470_new = tau_x412_new

      if (gzflg == 10) then
        if (xday >= 60.0 .and. xday < 274.0) go to 605
        tau_x470_new = tau_x412_new
      end if 
605   continue

631   continue

      if ((gzflg <= 5 .and. gzflg > 0) .OR. gzflg == 27) then

        tau_x412_flag = -999
!        if (r412ss > 19.89) print *, 'out, r412ss: ', r412ss,xlat,xlong

        imod = 4                                ! w0 = 0.91
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss_91,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_91)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
        
        imod = 5                                ! w0 = 0.94
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
  
        tau_x412ss_94 = tau_x412ss
        imod = 8                                ! w0=0.995
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss_995,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_995)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
        
        imod =  7     ! w0=0.98
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss_98,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_98)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
        

        imod = 6      ! w0=0.96
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss_96,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_96)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
        
        imod = 5      ! w0=0.95
        model_frac = 0.5
        call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma, 
     1      imod,r412ss,tau_x412ss_95,tau_x412_flag,trflg,model_frac,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_95)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if

!        if (Dstar1 .lt. 1.03 .and. scat_ang.gt.160.0
!     1      .and. tau_x412ss.gt.0.5) 
!     1       tau_x412ss = tau_x412ss_98
 
!        if (Dstar1 .lt. 0.98 .and. scat_ang.gt.150.0
!     1      .and.scat_ang.le.160.0
!     1      .and. tau_x412ss.gt.0.5) 
!     1       tau_x412ss = tau_x412ss_96

!       General surfaces,  tau_x412ss = default (0.94)
        if (xday < 60.0 .or. xday > 334.0) then                 ! Dec, Jan, Feb
          if (tau_x412ss_97 <0.6) tau_x412ss = tau_x412ss_96     
        endif
        if (xday > 243.0 .and. xday < 274.0) then               ! Sep
          if (tau_x412ss_97 <0.6) tau_x412ss = tau_x412ss_96     
        endif
        if (xday > 273.0 .and. xday < 335.0) then               ! Oct, Nov
          if (tau_x412ss < 0.5)   tau_x412ss = tau_x412ss_96
        endif
        if (xday > 59.0 .and. xday < 121.0) then               ! Mar,Apr
          if (tau_x412ss < 0.5)   tau_x412ss = tau_x412ss_96
        endif

!      Bright Surfaces        
        if (r412_tbl > 12.0) then                
          if (tau_x412ss_97 < 0.6) tau_x412ss = tau_x412ss_96 
          if (xday > 243.0 .or. xday < 60.0) then             ! Sep,Nov,Dec,Jan,Feb
            if (tau_x412ss_995 < 0.5) tau_x412ss = tau_x412ss_995 
            if (xday > 273.0 .and. xday < 305.0) then             ! Oct
              if (tau_x412ss_97 < 0.6) tau_x412ss = tau_x412ss_96 
            endif   
            if (xlat >10.0 .and. xlat <21.0 .and. xlong > 10.0 .and. xlong < 20.0)   
     1          tau_x412ss = tau_x412ss_96                        
          endif   
          if (xday > 59.0 .and. xday < 121.0) then             ! Mar,Apr
            if (tau_x412ss_995 < 0.5) tau_x412ss = tau_x412ss_995
            if (xlat >10.0 .and. xlat <21.0 .and. xlong > 10.0 .and. xlong < 20.0)
     1          tau_x412ss = tau_x412ss_96
            if (xday > 90.0 .and. xday < 121.0) then             ! Apr
            if (xlat >10.0 .and. xlat <21.0 .and. xlong > 10.0 .and. xlong < 20.0)
     1          tau_x412ss = tau_x412ss_98
            endif   
          endif   
        endif   

!      Deserts near Libya and Egypt         
        if (xlat >20. .and. xlong > 12.9) then               
          if (r412_tbl > 8.5) then
              if (Dstar1 .lt. 1.01) then
              dd = (xlong - 11.9) /3.0
              if (xlong > 14.9) dd = 1.
              tau_x412ss = tau_x412ss-(tau_x412ss-tau_x412ss_995)*dd
              endif   
          endif   
        endif   
        if (xlat >15. .and. xlong > 22.0) then                
          if (r412_tbl > 8.5) then
              if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
          endif
        endif   
        if (xlat >19.0 .and. xlat <22.0 .and. xlong > 30.0) then
          if (r412_tbl > 8.5) then
              if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
          endif
        endif   

!     NE Mauritania 1
      if (xlat>20.0 .and. xlat<30.0 .and. xlong< -12.5) then
      if (r412_tbl > 9.)  then
          if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
      endif
      endif
 
!     NE Mauritania 2
      if (xlat>26.0 .and. xlat<30.0.and. xlong>=-12.5.and. xlong< -11.001) then
      if (r412_tbl > 9.)  then
          if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
      endif
      endif

!     Morocco 1
      if (xlat>30.0 .and. xlat<36.0 .and. xlong<= -7.5) then
      if (r412_tbl > 5.8)  then
          if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
      endif
      endif   

!     Morocco 2
      if (xlat>31.5 .and. xlat<36.0 .and. xlong<= -5.0.and. xlong> -7.5) then
      if (r412_tbl > 5.8)  then
          if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
      endif
      endif  

        go to 635
633     continue
        if (tau_x412ss_97 < 0.6) tau_x412ss = tau_x412ss_96
        if (xday > 181.0 .and. xday < 274.0) then               ! Jul, Aug, Sep
          if (tau_x412ss_995 <0.5) tau_x412ss = tau_x412ss_995
        endif
        if (xday > 273.0 .and. xday < 335.0) then               ! Oct, Nov
          if (tau_x412ss_995 <0.5) tau_x412ss = tau_x412ss_995
        endif   
        if (xday > 59.0 .and. xday < 121.0) then               ! Mar,Apr
          if (tau_x412ss_995 <0.5) tau_x412ss = tau_x412ss_995
        endif

635   continue
      endif   
       
      if (gzflg >= 6 .and. gzflg <= 11) then

        tau_x412_flag = -999
        if (r412ss2 > 20.0) print *, 'out, r412ss2: ', r412ss2
        if (r412ss2 < 1.0) print *, 'less, r412ss2: ', r412ss2

        imod = 5                                ! w0 = 0.94
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss2,tau_x412ss2,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss2)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if

        imod =  7                               ! w0=0.98
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1               imod,r412ss2,tau_x412ss2_98,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss2_98)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if

        imod =  8                               ! w0=0.995
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1               imod,r412ss2,tau_x412ss2_995,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss2_98)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
      end if
 
c  If surface is vegetated, derive the surface reflectivity at 490 nm and 670 nm
c    from 865 nm surface LER and NDVI.  Specifically exclude Sahara and Arabia via gzflg
c    flag.  Exclude areas with AERONET-based BRDF info via brdf_flag == 0.
c--------------------------------------------------------------------------------------------------
 637  continue
      ndvi_thold = 0.1
      if (lc == 6) then
        ndvi_thold = 0.2
      else
        ndvi_thold = 0.1
      end if

!     -- increase NDVI threshold over high elevation area in the United States (13)
      if (gzflg == 13) then   
        if (px_elev >= 500.0 .OR. lc == 6) then
            ndvi_thold = 0.3
        end if
      endif
      
!     -- increase NDVI threshold over high elevation areas in North and South America.
!     -- this covers Central America which does not currently have a region assigned to it.
      if (xlong < -30.0 .AND. px_elev >= 500) then
          ndvi_thold = 0.3
      end if
      
      do_veg = toa_ndvi >= ndvi_thold .AND. (gzflg < -900 .OR. gzflg > 11) .AND. brdf_flag /= 0

c     -- exclude vegetated algorithm from India regions and Asia desert areas.
      if (gzflg == 15 .OR. gzflg == 20 .OR. gzflg == 19 .OR. gzflg == 21 .OR. gzflg == 23) do_veg = .false.
      if (gzflg == 24 .OR. gzflg == 27) do_veg = .false.
			if (debug) print *, "ndvi, gzflg, brdf_flag, sr_fail_flag, do_veg: ", toa_ndvi, gzflg, brdf_flag, do_veg

c     -- skip pixel if no surf. reflc. value and not suitable for vege. retrieval.
      if ((do_veg .eqv. .false.) .AND. (sr_fail_flag .eqv. .true.)) go to 10

      if (do_veg) then ! .OR. oob_samer) then
       
c      if (ndvi670 >= 0.1 .AND. (gzflg < -900 .OR. gzflg > 11) .AND. brdf_flag /= 0) then
        alg_flag = 1
c        print *, 'doing veg retrieval: ', itmp, jtmp, alg_flag
c       Get current season and tweak for swf_aero_veg input.
c        iopss = season_from_doy(yr, doy)  
c        iopss = iopss - 1
c        if (iopss == 0) iopss = 4

c... *******************************************************************
c... *******************************************************************
c... *******************************************************************
c-----------------------------------------------------------------------
c... do retrieval over vegetated surfaces (NDVI'=>0.1)
c-----------------------------------------------------------------------
c...  imod=2  ! ssa_490=0.94 --> You may open this line and remove
c...                             do loop  "do 2500 imod=1,4", or use
c...  --> "call swf_aero_veg(nvalx470,nvalx650,iopss,2,sza,xthet,"
c...  ------------------------------------------------------
c        gdt1 = gregorian_from_doy(yr,doy)
c        ioprg=0   ! region index initialization @@new@@
c        idx=int((xlong-(-180.0))/0.10)+1    ! @@new@@
c        idy=int((xlat-(-90.0))/0.10)+1     ! @@new@@
c        if(idx.ge.1.and.idx.le.3600.and.idy.ge.1.and.idy.le.1800) then   ! @@new@@
c           ioprg=veg_regions(idx,idy)   ! @@new@@
c        else          ! @@new@@
c           ioprg=0    ! @@new@@
c        endif         ! @@new@@
c        imod = 2
c        select case (ioprg)
c          case (6)                        ! S. Africa
c            select case (gdt1%month)  
c              case (6:11)       ! Jun-Nov use more absorbing model (0.89).
c                imod = 1
c              case default
c                imod = 2
c            end select
c          case default
c            imod = 2
c        end select
        
c        tau_x470_flag = -999
c        tau_x650_flag = -999
        
        call find_v_veg_test(itmp,jtmp,realbuf,tmpvg,outbufvg)
        
c       translate outbufvg back to local variables.
        do i=1,3
          xtau(i) = outbufvg(i)
          ssa(i) = outbufvg(i+3)
          outbuf(i+3) = ssa(i)
        enddo
        tau550 = outbufvg(7)
        alpha = outbufvg(8)
        r412 = -999.0
        r470 = outbufvg(11)
        r650 = outbufvg(12)
        tau_x650_flag = outbufvg(10)
        xthet = outbufvg(13)
        scat_ang = outbufvg(14)
        sfc_typ = outbufvg(15)
        if (debug) print *, "veg final: aot550, aot: ", tau550, xtau(1), xtau(2), xtau(3)
c     Return fill for 412 nm, no dark target retrieval yet.
c        xtau(1) = -999.0
c        tau_x412_flag2 = -999
c        if (lprint > 0) print *, 'after veg: ', i, j, xtau(2), xtau(4), tau550, & 
c        & alpha, tau_x470_flag, tau_x650_flag
c        
c!       -- reset AOT to zone-specific value if we ran off the LUT.
c        status = handle_lut_out_of_bounds(gzflg, tau_x470_flag, xtau(2))
c        if (status /= 0) then
c          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
c          return
c        end if
c        status = handle_lut_out_of_bounds(gzflg, tau_x650_flag, xtau(4))
c        if (status /= 0) then
c          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
c          return
c        end if
c        status = handle_lut_out_of_bounds(gzflg, tau_x470_flag, tau550)
c        if (status /= 0) then
c          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
c          return
c        end if
        
      
c        write (500,'(7(F10.5,1X))') xlat, xlong, r470, r650, xtau(2), xtau(4), tau550
        
c        if (lprint > 0) print *, 'after LUT check, veg: ', i, j, xlat, xlong, ndvi670, &
c        &   r470, r650, xtau(2), xtau(4), tau550, alpha
        
c!     Check consistency
        if (xtau(2) < -900.0) then
          xtau(2) = -999.0
          xtau(3) = -999.0
          tau550  = -999.0
          alpha   = -999.0
          tau_x470_flag = -999
          tau_x650_flag = -999
        end if
        if (debug) print *, 'final: ', xlat, xlong, xtau(2), xtau(3), tau550, alpha
c        write(888,'(7(F12.5), I3)') xlat, xlong, xtau(2), xtau(4), r470, r650, z_ndvi, proc_flag(i,j)
c-----------------------------------------------------------------------
c... end of retrieval over vegetated surfaces
c-----------------------------------------------------------------------
c... *******************************************************************
c... *******************************************************************
c... *******************************************************************

c         r470 = get_ssr_490(yr, doy, ler650(i,j)*100.0, ler865(i,j)*100.0, status)
c        if (status < 0) then
c          print *, "ERROR: Unable to derive vegetated surface "//
c     1      "reflectivity for 490 nm: ", status
c          stop
c          cycle
c        end if
c       
c        r650 = get_ssr_670(yr, doy, ler650(i,j)*100.0, ler865(i,j)*100.0, status)
c        if (status < 0) then
c          print *, "ERROR: Unable to derive vegetated surface "//
c     1      "reflectivity for 670 nm: ", status
c          stop
c          cycle
c        end if
        go to 865
      end if
      
      alg_flag = 0
c      if (lprint > 0) print *,'tau_x412,tau_x412_new,tau_x412ss,tau_x412ss2=',  
c     1    tau_x412,tau_x412_new,tau_x412ss,tau_x412ss2
c--------------------------------------------------------
c          Retrieving 412 nm SSA (412 - 470 nm)
c--------------------------------------------------------
      if (Dstar1.gt.1.05.and.rat1.gt.1.6) go to 620
      if (tau_x470.lt.0.2.and.tau_x470.gt.0.0) go to 805
      if (rr412_mod.gt.20.0.and.tau_x470.lt.0.0) go to 620
      if (tau_x470.lt.0.0) go to 620
      
      refl = refl1
      x3 = xphi
      tau_x = tau_x470
      w0_x = -999.
     
      if (tau_x  >=  tau(10)) tau_x = tau(10)-0.0001     ! search2 dies if tau_x >= tau(10), 3.5
      call search2(dflag,tau_x,tau,10,index_ia,frac_ia) ! quick fix to keep this going... 
      if (dflag) go to 10                               ! @TODO fix these search things so they stay in the tables!

      call aero_412_abs(dflag,refl,x1,x2,x3,mm,nn,ll,
     1    r412,index_ia,frac_ia,w0_x)
      if (dflag) go to 10

      w0_int = w0_x

      if (w0_x.lt.0.0) go to 10
 
      w0_int_470 = w0_x +(0.976 -w0_x)*(470.-412.)/(650.-412.)
 
      if (w0_int_470.lt.0.0) w0_int_470 = -999.
c--------------------------------------------------------
c          Retrieving 650 nm AOT
c--------------------------------------------------------
      if (xthet.gt.60.0.and.xphi.lt.90.0.and.tau_x470.gt.0.5)
     1    go to 620
      if (xlong < -30.0) then	! mostly N.America
      	if (r650.lt.30.0.and.tau_x470.lt.0.7.and.tau_x470.gt.0.0)
     1    go to 805
      else
      	if (r650.lt.30.0.and.r650.gt.15.0.and.tau_x470.lt.0.7.and.tau_x470.gt.0.0)
     1 		go to 805
     	end if
c      if (r650.lt.30.0.and.tau_x470.lt.0.7.and.tau_x470.gt.0.0)
c     1    go to 805
      if (r650.ge.30.0.and.tau_x470.lt.1.0.and.tau_x470.gt.0.0)
     1    go to 805
      if (scat_ang.gt.165.0.and.tau_x470.lt.1.8) go to 805
      
620   continue
      refl = refl6
      x3 = xphi

      tau_x650_flag = -999
      call aero_650(dflag,refl,x1,x2,x3,mm,nn,ll,ma,    
     1    r650,tau_x650,tau_x650_flag,tau_x470_flag2,   
     1    tau_x412, tau_x470,tau_x412_flag_91,trflg)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x650_flag, tau_x650)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
!--------------------------------------------------------
!          Retrieving 412 and 470 nm SSA
!--------------------------------------------------------
      if (tau_x650 < 0.0) go to 805
c      if (tau_x650 > 3.5) go to 805

!--- 412 nm

      refl = refl1
      x3 = xphi
      tau_x = tau_x650
      if (rat1.lt.1.8.and.tau_x470ss.gt.0.8) tau_x = tau_x470ss
      w0_x = -999.
      if (tau_x650 > 3.5) tau_x = 3.5

      if (tau_x  >  tau(10)) then
         goto 805
      endif
      if (tau_x  >=  tau(10)) tau_x = tau(10)-0.0001 
      call search2(dflag,tau_x,tau,10,index_ia,frac_ia)
      if (dflag) go to 10

      call aero_412_abs(dflag,refl,x1,x2,x3,mm,nn,ll,   
     1                 r412,index_ia,frac_ia,w0_x)
      if (dflag) go to 10

      w0_x412 = w0_x

      if (w0_x < 0.0) go to 10

!--- 470 nm 
 
      refl = refl3
      x3 = xphi
      tau_x = tau_x650
      w0_x = -999.
      if (tau_x650 > 3.5) tau_x = 3.5
      if (tau_x  >=  tau(10)) tau_x = tau(10)-0.0001 
      call search2(dflag2,tau_x,tau,10,index_ia,frac_ia)
      if (dflag2) go to 10
      
! ADDED BY COREY
      if (r470  >  24.0) go to 805
! END ADDED BY COREY
      call aero_470_abs(dflag2,refl,x1,x2,x3,mm,nn,ll,  
     1                 r470,index_ia,frac_ia,w0_x470)
      if (dflag2) go to 10
 
      if (w0_x470 < 0.0) w0_x470 = -999.

 805  continue

!
!     -- Selecting Models
!

      call aero_mod (tau_x412,tau_x470,tau_x650,tau_x412_91,aot_mod)

      if (tau_x412 > 0.0 .and. tau_x470 > 0.0) xxrat = tau_x412 / tau_x470
      if (tau_x650 > 0.0 .and. tau_x470 > 0.0) xxrat2 =  tau_x470 / tau_x650

      if (xxrat < 0.0) then
         alpha   = -999.
         go to 806
      end if

      dd      = alog(412./470.)
      alpha   = alog(xxrat)
      alpha   = -1.*alpha/dd
      
 806  continue
      dd    = ref650 / rr470_mod
      dd1   = ref650 / rr412_mod
      dd2   = rr470_mod / rr412_mod
      sfcdd = r650 / r412
      view  = xthet

      if (r650 < 30.0 .and. r650 > 9.0 .and. sfcdd < 2.4)  go to 850
      if (r650 < 30.0 .and. r650 > 9.0 .and. sfcdd > 2.4 .and. 
     1   sfcdd < 3.2 .and. trflg > 0.)  go to 850

      aot = aot_mod(3)
      if (tau_x412_91 < 0.0 .and. tau_x650 > 0.0 .and. r650 < 30.0) aot = aot_mod(1)
      if (tau_x412_91 < 0.0 .and. r650 >= 30.0 .and. tau_x412 > 0.0) aot = aot_mod(4)
      if (aot < 0.3 .and. aot > 0.0) go to 860
      if (w0_int < 0.92 .and. aot > 0.0) aot = aot_mod(2)
      if (r650 < 30.0 .and. tau_x412_91 > 1.2                       
     1      .and. tau_x412_91 > tau_x650 .and. tau_x650 > tau_x470  
     1      .and. w0_int >=  0.939) aot = aot_mod(1)
      if (tau_x412_91 < 0.0 .and. tau_x470 < 0.0 .and. tau_x650 > 0.0) aot = aot_mod(1)
      if (aot < 0.0 .and. tau_x650 > 0.0) aot = aot_mod(1)
      if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 > 1.35) alpha = -0.4
      if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 <= 1.35) alpha = 1.8
      if (tau_x412 < 0.0 .and. tau_x470 > 0.0 .and. dd1 <= 1.4) alpha = 1.8
      if (aot < 0.0) go to 10
      go to 860

850   continue
      
      aot = aot_mod(2)
      if (tau_x650 > aot) aot = aot_mod(1)
      if (aot < 0.0) go to 10
      if (tau_x650 > 1.0 .and. dd2 < 1.1 .and. r650 > 15.) go to 10
      if (tau_x412 > 0.0 .and. tau_x650 > 0.4 .and. dd2 < 1.15   
     1     .and. r650 > 15.) go to 10
      if (tau_x412 > 0.0 .and. tau_x650 > 1.9 .and. dd2 < 1.2) go to 10
      if (tau_x412 > 0.0 .and. tau_x470 > 1.9 .and. tau_x650 > 0.4    
     1     .and. dd2 < 1.2) go to 10
      if (tau_x412 > 0.0 .and. aot > 1.0 .and. dd2 < 1.3 .and.    
     1    tau_x650 > 0.4 .and. r650 > 13.0 .and. r650 <= 15.) go to 10

      if (tau_x412 > 0.0 .and. r412 > 18.0 .and. tau_x650 > 0.4   
     1     .and. dd < 1.4 .and. r650 <= 15.) go to 10

      if (tau_x412 > 1.0 .and. xxrat > 1.2 .and. r650 > 15.0 .and. tau_x650 > 0.4) go to 10
      if (tau_x650 > 1.5 .and. xxrat > 1.05 .and. tau_x412 < tau_x650 .and. tau_x412 > 0.0) go to 10
      if (tau_x650 > 1.5 .and. tau_x412 < 0.0 .and. tau_x470 < 0.0 .and.  
     1   dd2 < 1.6 .and. w0_x412 > 0.96) go to 10
      if (tau_x650 > 1.2 .and. tau_x412 < 0.0 .and. xxrat2 > 1.2 .and. w0_x412 > 0.97) go to 10
      

      if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 > 1.35)   alpha = -0.4
      if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 <= 1.35)  alpha = 1.8
      if (tau_x412 < 0.0 .and. tau_x470 > 0.0 .and. dd1 <= 1.4)   alpha = 1.8

      if (tau_x650 > 2.0 .and. tau_x412 < 0.0 .and. tau_x470 > 2.0      
     1     .and. tau_x470 < 3.0 .and. xxrat2 < 1.2 .and. xxrat2 > 1.0) 
     1    go to 10 
      
      if (tau_x650 > 2.0 .and. tau_x412 < 0.0 .and. tau_x470 > 3.0  
     1     .and. xxrat2 < 1.45 .and. xxrat2 > 1.0) go to 10

      if (tau_x412 < 0.0 .and. tau_x470 > 0.0 .and. xxrat2 > 2.)    
     1    aot = aot_mod(5)
      if (tau_x412 > 1.5 .and. tau_x470 > 0.0 .and. tau_x650 < 0.3) 
     1    aot = aot_mod(5)
      if (alpha > 1.0 .and. tau_x470 > 0.2) aot = aot_mod(5)
      if (alpha > 1.0 .and. tau_x470 <= 0.2) aot = aot_mod(5)*0.75
      if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 <= 1.35)  aot = aot_mod(6)

860   continue
      tau550 = aot
c     --- Additional Cloud Screening
      if ((gzflg <= 6 .OR. gzflg == 27) .and. (gzflg > 0 .and. gzflg /= 2)) then
      if (tau_x412ss > 0.6.and.rat1 <1.4.and.Dstar1 < 1.0.and.r412_tbl < 10.) go to 10
      endif
      
      !     -- skip all of the stuff below, just use aot.
      if (abs_aero_flag .eqv. .true.) then
        goto 864
      end if
      ! -- force Tinga_Tingana scheme for AOT @ 500nm over barren surfaces even when gzflag is undefined.
      if (lc == 6 .AND. (gzflg /= 16 .AND. gzflg /= 2) .AND. gzflg /= 14 .AND. gzflg /= 21) then  
        tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 12, lc, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     alpha, status, (lprint > 0))
        if (status /= 0) then     
          tau550 = aot
        endif
      endif
       
      if (gzflg > 0) then
!       -- get AOT @ 500nm based on AERONET regions and case studies 
        tau550 = get_aot500(xlat, xlong, px_elev, scat_ang, season, toa_ndvi, gzflg, lc, stdv, 
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,        
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,     
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,      
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,     
     1                     alpha, status, (lprint > 0))
        if (status /= 0) then! revert to manual assignment       
          tau550 = aot
        endif

!       -- force Kanpur (gzflg=15, lc=3) AOT models over gzflg 20 (Thar Desert)
        if (gzflg == 20) then
          tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 15, 3, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     alpha, status, (lprint > 0))
          if (status /= 0) then     
            tau550 = aot
          endif
        endif
        
!       -- force Beijing (gzflg=16, low elev) AOT models over gzflg 21 (Asian desert)
        if (gzflg == 21) then
          tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 16, 4, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     alpha, status, (lprint > 0))
          if (status /= 0) then     
            tau550 = aot
          endif
        endif
        
!       -- force SW Asia (Pakistan, Iraq, etc) (gzflg=23) to use Tinga_Tingana
        if (gzflg == 23) then
          tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 12, 6, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     alpha, status, (lprint > 0))
          if (status /= 0) then     
            tau550 = aot
          endif
        endif
        
!       -- force Tinga_Tingana scheme for AOT @ 500nm over barren surfaces
        if (lc == 6 .AND. (gzflg /= 16 .AND. gzflg /= 2) .AND. gzflg /= 14 .AND. gzflg /= 20 .AND. 
     1      gzflg /= 21 .AND. gzflg /= 23 .AND. gzflg /= 24) then
          tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 12, lc, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     alpha, status, (lprint > 0))
          if (status /= 0) then     
            tau550 = aot
          endif
          
          if ((gzflg <= 6 .OR. gzflg == 27) .and. (gzflg > 0 .and. gzflg /= 2)) then
          tau550 = tau_x412ss    
          if (Dstar1 > 1.06 .and. tau_x470ss > 0.6) tau550 = tau_x470ss
          if (xday > 31.0.and.xday < 60.0) then   ! Feb
            if (Dstar1 > 1.06 .and. tau_x470ss > 0.4) tau550 = tau_x470ss
          endif
          if (tau550 < 0.1) tau550 = tau550 + 0.05    ! check for geo zone

          if (xday > 59.0.and.xday < 121.0) then  ! Mar, Apr
          tau550 = tau_x412ss                     

          if (xlat <20.0 .and. xlong > 7.5 .and. xlong < 20.0) go to 820
          if (xlat >=20.0 .and. xlat <22.5 .and. xlong > 7.5 .and. xlong < 14.8) go to 820
          if (r412_tbl > 10.0) tau550 = tau_x412ss_995
 820      continue

          if (Dstar1 > 1.06 .and. tau_x470ss > 0.6 .and. r412_tbl < 12.0) tau550 = tau_x470ss
          if (r412_tbl >= 12.0) then
          if (Dstar1 > 1.09 .and. tau_x470ss > 0.8 .and. xlong > 2.5) tau550 = tau_x470ss
          if (Dstar1 > 1.06 .and. tau_x470ss > 0.6 .and. xlong <= 2.5) tau550 = tau_x470ss
          endif
          endif      

          if (xday > 120.0.and.xday < 152.0) then     ! May
          tau550 = tau_x470ss
          if (Dstar1 > 1.06 .and. tau_x412ss > 0.6) tau550 = tau_x412ss
          if (Dstar1 < 1.03 .and. tau550 > 0.6) then
            if (tau_x412ss < tau_x470ss.and.r412_tbl >= 10.0) tau550 = tau_x412ss_98
          endif
            if (xlat >10.0 .and. xlat <21.0 .and. xlong > 10.0 .and. xlong < 20.0) then
            if (r412_tbl >= 11.0.and.tau550 < 0.6)
     1          tau550 = tau_x412ss_94
          endif
          endif

          if (xlat < 18. .and. xlong > 25.0) then
            if (Dstar1 > 1.09.and.rat1 > 1.6) then
             if (tau_x650 > 0.0.and.tau_x650 > tau550) tau550 = tau_x650
             if (tau_x650 < 0.0) tau550 = tau_x470ss *2.
            endif
          endif  

          if (xday > 151.0.and.xday < 258.0.and.rat1 > 1.6.and.
     1        Dstar1 > 1.2.and.tau_x650 < 0.0.and.tau550 < 0.8)
     1        tau550 = tau_x470ss * 2.
 
          if (Dstar1 > 1.2.and.rat1 > 1.6.and.tau_x412ss < 0.5 .and.tau550< 0.7)
     1        tau550 = Dstar1
 

          end if
 
          if (gzflg >= 6 .and. gzflg <= 11) then
               tau550 = tau_x470ss
               if (tau550 < 0.45) tau550 = (tau_x470ss+tau_x412ss2)/2.
               if (xday > 120.0 .and. xday < 151.0) then                   ! May
                tau550 = tau_x470ss
                if (tau550 < 0.3.and.r412_tbl < 11.0) tau550 = tau_x412ss2 
               end if
               if (xday > 90.0 .and. xday < 121.0) then                    ! Apr
                   tau550 = tau_x412ss2
                   if (tau550 < 0.5 .and. tau_x470ss > tau550) tau550 = tau_x470ss
               end if 
          if (xlat > 29.5 .and. Dstar1 < 0.97) then
             if (tau550 > 0.6) tau550 = tau_x412ss2_995
          end if 

          if (xday < 182.0) then
          if (xlat > 20.0.and. xlat < 28.0.and.xlong > 42.0.and. xlong < 50.0) then
             if (Dstar1 < 1.06.and.r412_tbl > 12.0.and.tau550 > 0.6) go to 10
          end if
          end if

          if (Dstar1 > 1.1 .and. tau_x470ss <0.45) tau550 = Dstar1

          if (xlat > 20.5.and. xlat < 25.5.and.xlong > 47.5.and. xlong < 50.0) then
             if (Dstar1 < 1.06.and.r412_tbl < 9.5.and.tau550 > 0.9) go to 10
          end if

          end if 
           
        end if    ! end of barren
 
      end if  
      
!     -- for Arabian peninsula, just use AOT412 w/ w0=0.96. Override all of that above.
      if (gzflg >= 6 .and. gzflg <= 11) then
        if ((r650_135 > 32.0) .AND. (r650_135/r412_135) > 3.7) then
          if ((r412_135 < 7.25) .OR. r650_135 > 36.0 .OR. use_alternate_brdf) then
            continue
           
          else if (xday >= 152 .AND. xday <= 243) then     ! summer
            tau550 = (tau_x412_new_94+tau_x412_new_96)/2.0    ! orig = 0.96

          else if (xday >= 244 .AND. xday < 335) then
            tau550 = (tau_x412_new_94+tau_x412_new_96)/2.0      ! orig = 0.995

          else
            tau550 = tau_x412_new_94
          end if 
          
        end if
      end if
      
      if (lprint > 0) print *, 'after get_aot500: xlat, xlong, tau550: ', gzflg, xlat, xlong, tau550
!     Weakly absorbing dust
!      if (Dstar1 < 1.01.and.tau550 > 0.8.and.rat1 < 1.6.and.
!     1    gzflg <6 .and. gzflg >0)! .and. lc .eq. 6) 
!     1    tau550 = tau_x412ss_98
!      if (lprint > 0) print *, 'after absorbing dust check: tau550, Dstar1, rat1: ',  
!     1  tau550, Dstar1, rat1, rr470_mod, rr412_mod
      

!     Smoke pixels
      if (((gzflg <6 .and. gzflg >0) .OR. gzflg == 27) .and. xlat <20.0 .and. lc <5) then
        if (xday > 335.0 .or. xday < 32.0) then
          if (tau550 < 0.065 .or. tau_x470 <0.05) then
            dd412 = rr412_mod - r412_tbl
            dd470 = rr470_mod - r470_tbl
            dd650 = ref650    - r650_tbl
          
            if (dd412 < 0.0 .AND. dd470 < 0.0 .AND. dd650 < 0.0) then
              call smoke_mod(xthet,scat_ang,dd412,dd470,dd650,  
     1                       tau_x412ss,tau_x412ss_91,tau_smoke)
              if (debug) print *, 'smoke_mod, dd412, dd470, dd650,t412ss, t412_91, tsmoke: ',  
     1                              dd412, dd470, dd650, tau_x412ss, tau_x412ss_91, tau_smoke
              tau550 = tau_smoke
              if (tau_smoke < 0.25) then
                tau550 = tau_smoke *2.0
                if (debug) print *, 'smoke detected, aot550 reset: ', tau550
              end if
            end if        
          end if
        end if
      end if
      if (lprint > 0) print *, 'after smoke check: tau550: ', tau550

!     -- adjust southern Sahel in Winter
!      if ((gzflg == 5 .AND. r650_135 < 9.0) .AND. (xday < 60.0 .OR. xday > 334)) then
!        if (tau550 < 1.5) tau550 = tau550 + 0.4
!      end if
      
!      if (xlat >28. .and.xlat <31. .and. xlong >9.0 .and. xlong <12.0 &
!       &  .and. r412_tbl > 9. .and. tau_x412ss <0.2 .and. &
!       &  tau_x470ss > 0.3) tau550 = tau_x470ss
!      if (xlat >27. .and.xlat <30. .and. xlong >3.0 .and. xlong <10.0 &
!       &  .and. r412_tbl > 9. .and. tau_x412ss <0.2 .and. &
!       &  tau_x470ss > 0.3) tau550 = tau_x470ss
      
      if (lprint > 0) print *,                                  
     1   'gzflg,lc, tau_x412ss,tau_x470ss,tau550 =',   
     1    gzflg,lc,tau_x412ss,tau_x470ss,tau550

864   if (lprint > 0) then
        print *, 'lat,lon,gzflg,r650,tau_x650,tau_x470,tau550: ', xlat, xlong, gzflg, r650, 
     1   tau_x650, tau_x470, tau550
        print *, 'abs_aero_flag: ', abs_aero_flag
      end if
      
c     -- Try to detect and remove smoke over clouds. For example, see RGB:
c     --  MYD021KM.A2007081.0640*.hdf
      if (debug) print *, 'smoke detection: ', xlat, xlong, r650, tau_x650, refl6, w0_x470
      if (r650 <8.0 .and. tau_x650 > 3.49 .and. refl6 > 0.1 .and.w0_x470 > 0.999) go to 10
      if (xday > 0.0 .and. xday < 121.0) then
 				if (gzflg==16 .and. (w0_x470 > 0.96 .and. alpha < 1.0) .and. Dstar1 <1.0 .and. tau550 >1.0) then
 					go to 10
 				end if
      endif
      
!     -- set tau550 to tau_x470ss over high elevation to minimize the effects of 
!     -- Rayleigh scattering. 
!     -- Morocco region
      if (px_elev > 750.0 .AND. (xlat > 28.0 .AND. xlat < 37.0 .AND. xlong > -12.0 .AND. 
     1    xlong < 10.0)) then
        tau550 = tau_x470
      endif
      
!     -- West Sudan region
      if (px_elev > 750.0 .AND. (xlat > 10.5 .AND. xlat < 19.5 .AND. xlong > 20.5 .AND. 
     1    xlong < 29.0)) then
        tau550 = tau_x470
      endif
      
c     -- check D* for dust plumes globally except for N.Africa (1-5, 26, 27), 
c     -- Arabian Peninsula (6-11), Taklimakan Desert (24), and Beijing (East China, 16) 
c     -- zones. Limit to locations at less than 4500 m elevation. Otherwise,
c			-- abnormally high AOT's are seen over Tibetan Plateau.
      if ((gzflg < 0 .OR. gzflg > 11) .AND. (gzflg /= 24 .AND. gzflg /= 26 .AND. 
     & gzflg /= 27 .AND. gzflg /= 16)) then
     		if (px_elev < 4750.0) then
        	if (Dstar1 >1.05) tau550 = tau_x650
        	if (Dstar1 >1.05 .and. tau_x470 > tau550) tau550 = tau_x470
        end if
      end if
      
      xtau(1) = tau_x412
      xtau(2) = tau_x470
      xtau(3) = tau_x650
      
      if (w0_x412 > 0.0) w0_x = w0_x412
      if (w0_x470 > 0.0) w0_int_470 = w0_x470

      if (alpha < 1.0)  ssa(1) = w0_x
      if (alpha < 1.0)  ssa(2) = w0_int_470
      if (ssa(1) > 0.0) ssa(3) = 0.976
c-- Set Surface Type

      sfc_typ = -999.

      if (terrain_flag_new5.gt.0.0) sfc_typ = 7.
      if (terrain_flag_new5.lt.5.5.and.terrain_flag_new5.gt.0.0)
     1    sfc_typ = 10.

C      dd = xsfc650(xlonp,xlatp) - xsfc650_bk2(xlonp,xlatp)
C      dd1 = xsfc470b_bk(xlonp,xlatp) - xsfc470_bk(xlonp,xlatp)
 
      if (terrain_flag_new5.lt.5.5.and.terrain_flag_new5.gt.0.0.and.
     1    xlat.le.15.0.and.xlong.lt.26.0) then
          sfc_typ = 1.
C          if (xlong.gt.10.0.and.abs(dd1).ge.0.2) then
C          if (xsfc470_bk(xlonp,xlatp).gt.10.0) sfc_typ = 2.
C          if (xsfc470_bk(xlonp,xlatp).le.10.0.and.
C     1        xsfc470_bk(xlonp,xlatp).gt.9.0.and.
C     1        xsfc650_bk(xlonp,xlatp).lt.22.0) sfc_typ = 2.
C          endif
C          if (xlong.gt.10.0) then
C          if (xsfc470_bk(xlonp,xlatp).gt.10.0.and.
C     1        xsfc650_bk(xlonp,xlatp).lt.23.0) sfc_typ = 2.
C          endif
      endif
C      if (terrain_flag_new5.lt.5.5.and.terrain_flag_new5.gt.0.0.and.
C     1    xlat.gt.15.0.and.xlat.le.20.0.and.dd.ge.1.6) then
C          sfc_typ = 1.
C          if (xlong.gt.10.0.and.abs(dd1).ge.0.2) then
C          if (xsfc470_bk(xlonp,xlatp).gt.10.0) sfc_typ = 2.
C          if (xsfc470_bk(xlonp,xlatp).gt.10.0.and.
C     1        xsfc650_bk(xlonp,xlatp).lt.23.0) sfc_typ = 2.
C          if (xsfc470_bk(xlonp,xlatp).le.10.0.and.
C     1        xsfc470_bk(xlonp,xlatp).gt.9.0.and.
C     1        xsfc650_bk(xlonp,xlatp).lt.22.0) sfc_typ = 2.
C          endif
C      endif
C      if (terrain_flag_new5.lt.5.5.and.terrain_flag_new5.gt.0.0.and.
C     1    xlong.ge.26.0.and.xlat.le.20.0.and.dd.ge.1.6) then
C          sfc_typ = 1.
C          if (xlong.gt.10.0.and.abs(dd1).ge.0.2) then
C          if (xsfc470_bk(xlonp,xlatp).gt.10.0) sfc_typ = 2.
C          if (xsfc470_bk(xlonp,xlatp).gt.10.0.and.
C     1        xsfc650_bk(xlonp,xlatp).lt.23.0) sfc_typ = 2.
C          if (xsfc470_bk(xlonp,xlatp).le.10.0.and.
C     1        xsfc470_bk(xlonp,xlatp).gt.9.0.and.
C     1        xsfc650_bk(xlonp,xlatp).lt.22.0) sfc_typ = 2.
C          endif
C      endif

      if (terrain_flag_new5.gt.2.5.and.terrain_flag_new5.lt.3.5)
     1    sfc_typ = 3.

      if (terrain_flag_new5.gt.1.5.and.terrain_flag_new5.lt.2.5) then
          sfc_typ = 5.
C          if (xsfc470_bk(xlonp,xlatp).gt.6.7.and.xlong.lt.-7.9)
C     1        sfc_typ = 4.
      endif

C      if (terrain_flag_new5.gt.4.5.and.terrain_flag_new5.lt.5.5.and.
C     1    xlat.gt.31.3.and.xlong.lt.-7.9.and.xsfc470_bk(xlonp,xlatp)
C     1    .gt.6.7.and.xsfc470_bk(xlonp,xlatp).lt.9.5)
C     1    sfc_typ = 4.

      if (terrain_flag_new5.gt.3.5.and.terrain_flag_new5.lt.4.5)
     1    sfc_typ = 6.

      if (terrain_flag_new5.gt.8.5.and.terrain_flag_new5.lt.9.5)
     1    sfc_typ = 8.

      if (terrain_flag_new5.gt.9.5.and.terrain_flag_new5.lt.10.5)
     1    sfc_typ = 9.

      if (terrain_flag_new5.gt.5.5.and.terrain_flag_new5.lt.6.5.and.
     1    xlat.gt.24.6.and.xlat.lt.25.2.and.xlong.gt.46.1.and.
     1    xlong.lt.46.8.and.r412_tbl.ge.7.0.and.
     1    r412_tbl.le.9.9)
     1    sfc_typ = 9.

c     -- Set Flags
      if (tau_x650_flag.gt.0.or.tau_x470_flag.gt.0) qa_flag(4)= 3
      if (tau_x412.lt.0.0.and.tau_x470.lt.0.0.and.
     1    tau_x650.lt.0.0) then
      qa_flag(1)= 0
      qa_flag(2)= 0
      qa_flag(3)= 0
      qa_flag(4)= 2
      endif
      if (tau_x412.gt.0.0.or.tau_x470.gt.0.0.or.
     1    tau_x650.gt.0.0) qa_flag(1)= 1
      if (tau_x412.lt.0.05.or.tau_x470.gt.0.05)
     1    qa_flag(2)= 1
      if (tau_x412.ge.0.05.and.tau_x412.lt.0.3)
     1    qa_flag(2)= 2
      if (tau_x412.ge.0.3)
     1    qa_flag(2)= 3
      if (alpha.lt.0.5.and.alpha.gt.0.0)
     1    qa_flag(3)= 1
      if (alpha.ge.0.5.and.alpha.lt.1.4)
     1    qa_flag(3)= 2

 865  continue

!     -- in zone 1, in summer, sometimes the BRDF for low NDVI pixels is too high 
!     -- resulting in an abnormally low AOT.  In these cases, perform a second retrieval
!     -- using a constant fit (an alternate fit) for 470nm that will bring up 
!     -- those low areas.
      if (tau550 < 0.1 .AND. (gdt1%month >= 6 .AND. gdt1%month <= 8) .AND. gzflg == 1) then
        if (.NOT. use_alternate_brdf) then
          use_alternate_brdf = .true.
          goto 5
        end if
      end if
      
!     -- over Arabian Peninsula, trouble with some very low AOT's in certain areas. Redo 
!     -- the retrieval if in summer.
      if (.NOT. use_alternate_brdf .AND. (tau550 < 0.25 .AND.  
     &    (gdt1%month >= 6 .AND. gdt1%month< = 8) .AND. gzflg >= 6 .AND. gzflg <= 11)) then
        if (r650_135 > 32.0 .AND. (r650_135/r412_135) > 3.7) then
          use_alternate_brdf = .true.
          goto 5
        end if
      end if
      
c-------------------------------------------------------------
c Set output buffer
c-------------------------------------------------------------
      if (debug) print *, 'final spectra aot: ', xtau(1), xtau(2), xtau(3)
      if (debug) print *, 'final tau550, ae: ', tau550, alpha
      do i=1,3
         outbuf(i) = xtau(i)
         outbuf(i+3) = ssa(i)
      enddo
      outbuf(7) = tau550
      outbuf(8) = alpha
      outbuf(9) = r412
      outbuf(10) = 1.0*tau_x650_flag
      outbuf(11) = r470
      outbuf(12) = r650
      outbuf(13) = xthet
      outbuf(14) = scat_ang
      outbuf(15) = sfc_typ
      outbuf(18) = float(alg_flag)
10    continue
      return
      end
c
c-------------------------------------------------------------
c
      subroutine new_intep(x1a,x2a,x3a,ya,m,n,l,ia,x1,x2,x3,y,dy,
     1                     mbeg,nbeg,frac)
      dimension x1a(m),x2a(n),x3a(l),ya(4,4,2,10)
      dimension xx2a(4), xx1a(4)
      dimension yntmp(4),ymtmp(4),yltmp(2)
 
      do 12 j=1,4
        do 11 k=1,4
          yltmp(1)=ya(j,k,1,ia)
          yltmp(2)=ya(j,k,2,ia)
          yntmp(k) = yltmp(1)*(1.-frac) + yltmp(2)*frac
          xx2a(k) = x2a(k+nbeg)
11      continue
        call polint(xx2a,yntmp,4,x2,ymtmp(j),dy)
        xx1a(j) = x1a(j+mbeg)
12    continue
      call polint(xx1a,ymtmp,4,x1,y,dy)
      return
      end
c
c-------------------------------------------------------------
c
      subroutine new_intepw(x1a,x2a,x3a,ya,m,n,l,ia,x1,x2,x3,y,dy,
     1                     mbeg,nbeg,frac)
      dimension x1a(m),x2a(n),x3a(l),ya(4,4,2,8)
      dimension xx2a(4), xx1a(4) 
      dimension yntmp(4),ymtmp(4),yltmp(2)
       
      do 12 j=1,4 
        do 11 k=1,4 
          yltmp(1)=ya(j,k,1,ia)
          yltmp(2)=ya(j,k,2,ia)
          yntmp(k) = yltmp(1)*(1.-frac) + yltmp(2)*frac
          xx2a(k) = x2a(k+nbeg)
11      continue
        call polint(xx2a,yntmp,4,x2,ymtmp(j),dy)
        xx1a(j) = x1a(j+mbeg)
12    continue 
      call polint(xx1a,ymtmp,4,x1,y,dy) 
      return
      end
c
c-------------------------------------------------------------
c
      subroutine polint(xa,ya,n,x,y,dy)
      parameter (nmax=50) 
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n 
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
c          if(den.eq.0.) print *,'den=',den
c          if(den.eq.0.)pause
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
c
c     --------------
c
      subroutine search(dflag,xbar,x,n,i)
c
c     purpose
c       locate position in table of point at which interpolation is
c       required
c
c     usage
c       call  search (xbar, x, n, i)
c
c     description of parameters
c       xbar   - point at which interpolation is required
c       x      - array containing independent variable
c       n      - number of points in x array
c       i      - index specifying segment containing xbar
c
      logical dflag
      dimension x(n)
      data b/.69314718/
      icnt = 0
      if (n.lt.2) go to 15
      if(x(1).gt.x(2)) go to 17
      m = int((log(float(n)))/b)
      i=2**m
      if (i .ge. n) i = n-1
      k=i
   10 k=k/2
      if (k .eq. 0) icnt = icnt + 1
      if (icnt .ge. 2) goto 27
      if (xbar.ge.x(i).and.xbar.lt.x(i+1))return
      if (xbar.gt.x(i)) go to 12
      i = i-k
      go to 10
   12 i = i+k
      if (i.lt.n) go to 10
      i=n-1
      go to 10
   15 print *, "Search n is less than 2."
      return
   17 print *, "Search table is not in increasing order."
      return

   27 continue
      do 22 i=1,n-1
         if (xbar.ge.x(i).and.xbar.le.x(i+1)) return
   22 continue           
      write(6,*) "setting dflag = true"
      dflag = .true.
      return

      end

      subroutine search2(dflag,xbar,x,n,i,fx)
c
c       call  search (xbar, x, n, i)
c
c     description of parameters
c       xbar   - point at which interpolation is required
c       x      - array containing independent variable
c       n      - number of points in x array
c       i      - index specifying segment containing xbar
c       fx     - fraction from i (between i and i+1)
c
      dimension x(n)
      logical dflag

      call search(dflag,xbar,x,n,i)

      fx   = (xbar-x(i))/ (x(i+1)-x(i))

      end

c--------------------------------------------------------
      subroutine aero_412_abs(dflag,refl,x1,x2,x3,mm,nn,ll,
     1    r412,index_ia,frac_ia,w0_x)

      include 'aottbl.inc'
      include 'newaottbl.inc'

      real nvalxw(10,46,30,8), nnvalxw(4,4,2,8), yyw(8)
      logical dflag
      data pi   /3.14159/

      dflag = .false.
      index_ii = r412
      frac = (r412-sfc_ref412(index_ii))/
     1       (sfc_ref412(index_ii+1)-sfc_ref412(index_ii))
 
      if (index_ii.lt.1.or.index_ii.gt.20) 
     1    print *,'index_ii = ', index_ii
      if (frac.lt.0.0.or.frac.gt.1.0) 
     1    print *,'aero_412_abs frac on sfc=', frac
 
      call search(dflag,x3,phi,ll,ii)
      if (dflag) return
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))
 
      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif
 
      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif

      do iw = 1, 8
       do i = 1, 4
        do j = 1, 4
        dd1 = nvalx412(mbeg+i,nbeg+j,ii,index_ia,iw,index_ii)*
     1  (1.-frac) + 
     1  nvalx412(mbeg+i,nbeg+j,ii,index_ia,iw,index_ii+1)*frac
        dd2= nvalx412(mbeg+i,nbeg+j,ii,index_ia+1,iw,index_ii)*
     1  (1.-frac)+
     1   nvalx412(mbeg+i,nbeg+j,ii,index_ia+1,iw,index_ii+1)
     2   *frac   
 
         nnvalxw(i,j,1,iw) = dd1* (1.-frac_ia) + dd2*frac_ia   
 
        dd1 = nvalx412(mbeg+i,nbeg+j,ii+1,index_ia,iw,index_ii)*
     1  (1.-frac) +
     1   nvalx412(mbeg+i,nbeg+j,ii+1,index_ia,iw,index_ii+1)*frac 
        dd2= nvalx412(mbeg+i,nbeg+j,ii+1,index_ia+1,iw,index_ii)* 
     1  (1.-frac)+
     1   nvalx412(mbeg+i,nbeg+j,ii+1,index_ia+1,iw,index_ii+1)
     2   *frac  
 
        nnvalxw(i,j,2,iw) = dd1* (1.-frac_ia) + dd2*frac_ia
 
        enddo
       enddo
      enddo
 
c---     interpolating W0 tables
 
      do 800 iw = 1, 8
 
      call new_intepw(theta0, theta, phi, nnvalxw, mm, nn, ll, iw,
     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
      yyw(iw) = y/pi
 800  continue
 
      if (refl.le.yyw(1)) then
        w0_x = 0.82
        return
      endif
 
      if (refl.ge.yyw(8)) then
        w0_x = 1.0
        return
      endif
 
      call search2(dflag,refl,yyw,8,index_ii,frac)
      if (dflag) return
      w0_x = frac*w0(index_ii+1) + (1.-frac)*w0(index_ii)

      return
      end

c--------------------------------------------------------
      subroutine aero_470_abs(dflag2,refl,x1,x2,x3,mm,nn,
     1    ll,r470,index_ia,frac_ia,w0_x)
 
      include 'aottbl.inc'
      include 'newaottbl.inc'
 
      real nvalxw(10,46,30,4), nnvalxw(4,4,2,4), yyw(4)
      logical dflag2
      data pi   /3.14159/
 
      dflag2 = .false.
      index_ii = r470
 
      frac = (r470-sfc_ref470(index_ii))/
     1       (sfc_ref470(index_ii+1)-sfc_ref470(index_ii))
 
      if (index_ii.lt.1.or.index_ii.gt.24) 
     1    print *,'index_ii r470abs = ', index_ii
      if (frac.lt.0.0.or.frac.gt.1.0) 
     1    print *,'frac on sfc=', frac
 
      call search(dflag2,x3,phi,ll,ii)
      if (dflag2) return
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))
 
      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif
 
      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif
 
      do iw = 1, 4
       do i = 1, 4
        do j = 1, 4
        dd1 = nvalx470(mbeg+i,nbeg+j,ii,index_ia,iw,index_ii)*
     1  (1.-frac) + 
     1  nvalx470(mbeg+i,nbeg+j,ii,index_ia,iw,index_ii+1)*frac
        dd2= nvalx470(mbeg+i,nbeg+j,ii,index_ia+1,iw,index_ii)*
     1  (1.-frac)+
     1   nvalx470(mbeg+i,nbeg+j,ii,index_ia+1,iw,index_ii+1)
     2   *frac   
 
         nnvalxw(i,j,1,iw) = dd1* (1.-frac_ia) + dd2*frac_ia   
 
        dd1 = nvalx470(mbeg+i,nbeg+j,ii+1,index_ia,iw,index_ii)*
     1  (1.-frac) +
     1   nvalx470(mbeg+i,nbeg+j,ii+1,index_ia,iw,index_ii+1)*frac 
        dd2= nvalx470(mbeg+i,nbeg+j,ii+1,index_ia+1,iw,index_ii)* 
     1  (1.-frac)+
     1   nvalx470(mbeg+i,nbeg+j,ii+1,index_ia+1,iw,index_ii+1)
     2   *frac  
 
        nnvalxw(i,j,2,iw) = dd1* (1.-frac_ia) + dd2*frac_ia
 
        enddo
       enddo
      enddo
 
c---     interpolating W0 tables
 
      do 800 iw = 1, 4
 
      call new_intepw(theta0, theta, phi, nnvalxw, mm, nn, ll, iw,
     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
      yyw(iw) = y/pi
 800  continue
 
      if (refl.le.yyw(1)) then
        w0_x = -999.
        return
      endif
 
      if (refl.ge.yyw(4)) then
        w0_x = 1.0
        return
      endif
 
      call search2(dflag2,refl,yyw,4,index_ii,frac)
      if (dflag2) return
      w0_x = frac*w0_470(index_ii+1) + (1.-frac)*w0_470(index_ii)
 
      return
      end

c--------------------------------------------------------
      subroutine aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,  
     1    imod,r470,tau_x470,tau_x470_flag,trflg,model_frac,debug)
     
        implicit none
        
        include 'aottbl.inc'
        include 'newaottbl.inc'
        
        logical, intent(inout)      ::  dflag
        real, intent(in)            ::  refl
        real, intent(in)            ::  x1
        real, intent(in)            ::  x2
        real, intent(in)            ::  x3
        integer, intent(in)         ::  mm
        integer, intent(in)         ::  nn
        integer, intent(in)         ::  ll
        integer, intent(in)         ::  ma
        integer, intent(in)         ::  imod
        real, intent(in)            ::  r470
        real, intent(inout)         ::  tau_x470
        integer, intent(inout)      ::  tau_x470_flag
        real, intent(inout)         ::  trflg
        real, intent(in)            ::  model_frac
        logical, intent(in)         ::  debug
        
        
        
!       -- NOTE: if model_frac = 0.0, the aerosol model = imod
!       -- if model_frac > 0.0, the aerosol model will be interpolated between
!       --  imod and imod + 1, using the value of model_frac

        real, dimension(4,4,2,10) ::  nnvalx, nnvalx1, nnvalx2 
        real, dimension(10)       ::  yy
        real, dimension(8)        ::  yy2
      
        real      ::  frac_ia, w0_x
        integer   ::  index_ia
        real      ::   tau_x
      
        integer ::  ia, i, ii, index_ii, nsm, nsn, mbeg, nbeg, iw, j
        real    ::  dif, dift, frac, xfrac, dd1, dd2, y, dy
        real    ::  xxrat
      
        real, parameter   ::  pi = 3.14159
        
        tau_x470_flag = -999  
        
        if (debug) print *, 'aero_470, in: ', refl, x1, x2, x3, r470, imod, model_frac

        index_ii = r470
 
        frac = (r470-sfc_ref470(index_ii)) / (sfc_ref470(index_ii+1)-sfc_ref470(index_ii))
        if (debug) print *, 'aero_470, sfc indx: ', r470, index_ii
 
        if (index_ii < 1 .or. index_ii > 24) print *,'index_ii = ', index_ii, r470
        if (frac < 0.0 .or. frac > 1.0) print *,'aero_470 frac on sfc=', frac
 
        call search(dflag,x3,phi,ll,ii)
        if (dflag) return
        xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))
        if (debug) print *, 'aero_470, raa: ', x3, ii, phi(ii), phi(ii+1)
        
        nsm = 1
        dif = x1 - theta0(1)
        do i = 1, mm
          dift = x1 - theta0(i)
          if (dift > 0.0 .and. dift < dif) then
            nsm = i
            dif = dift
          end if
        end do
        mbeg = nsm - 2
        if (mbeg <= 0) then
          mbeg = 0
        else if (mbeg > mm-4) then
          mbeg = mm-4
        end if
 
        nsn = 1
        dif = x2 - theta(1)
        do i = 1, nn
          dift = x2 - theta(i)
          if (dift > 0. .and. dift < dif) then
            nsn = i
            dif = dift
          end if
        end do
        nbeg = nsn - 2
        if (nbeg <= 0) then
          nbeg = 0
        else if (nbeg > nn-4) then
          nbeg = nn-4
        end if
        if (debug) print *, 'aero_470, mbeg, nbeg: ', mbeg, nbeg
        
!       -- interpolate between models if requested.
 !       if (present(model_frac)) then
!       -- temporary check to avoid out-of-bounds condition with imod=4
!       -- change to above when this file is converted to a proper module
        if (imod < 4) then
          nnvalx1(:,:,:,:)  = -999.0
          nnvalx2(:,:,:,:)  = -999.0
          nnvalx(:,:,:,:)   = -999.0
          do ia = 1, 10
            do i = 1, 4
              do j = 1, 4
                nnvalx1(i,j,1,ia) = nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*
     1             (1.0-frac) + nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac
                nnvalx1(i,j,2,ia) = nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*
     1             (1.0-frac) + nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac
  
                nnvalx2(i,j,1,ia) = nvalx470(mbeg+i,nbeg+j,ii,ia,imod+1,index_ii)*
     1             (1.0-frac) + nvalx470(mbeg+i,nbeg+j,ii,ia,imod+1,index_ii+1)*frac
                nnvalx2(i,j,2,ia) = nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod+1,index_ii)*
     1             (1.0-frac) + nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod+1,index_ii+1)*frac
  
                nnvalx(i,j,1,ia) = (1.0-model_frac) * nnvalx1(i,j,1,ia) +   
     1                                 model_frac * nnvalx2(i,j,1,ia)
                nnvalx(i,j,2,ia) = (1.0-model_frac) * nnvalx1(i,j,2,ia) +  
     1                              model_frac * nnvalx2(i,j,2,ia)
              end do
            end do
          end do
          
!       -- just use imod, no interpolation
        else      
          do ia = 1, 10
            do i = 1, 4
              do j = 1, 4
                nnvalx(i,j,1,ia) = nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*
     1             (1.0-frac) + nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac
                nnvalx(i,j,2,ia) = nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*
     1             (1.0-frac) + nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac
              end do
            end do
          end do
        end if
       
!---     interpolating AOT tables
        yy(:) = -999.0
        do ia = 1, 10
 
          call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,
     1               x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
          yy(ia) = y/pi
!        print *,'tau, i/f=', tau(ia), y/pi, dy
        end do
        
        if (refl <= yy(1)) then
          tau_x470 = 0.06 
          if (trflg > 0.0) tau_x470 = 0.02
          tau_x470_flag = -10
          if (debug) print *, 'aero_470, hit low bound: ', refl, yy(1)
          return
        end if
 
! Reflc off the charts!  Set AOT to max and set flag.
        if (refl >= yy(10)) then
          xxrat = 0.8
          tau_x470 = 3.5
          tau_x470_flag = 1
          if (debug) print *, 'aero_412, hit hi bound: ', refl, yy(10)
          return
        endif

!
!     Check if the reflectance increase with AOT
!
 
        if (yy(1) < yy(2)) go to 650
 
        if (refl < yy(4)) return
      
        yy2(:) = -999.0
        do i = 1, 7
          yy2(i) = yy(i+3)
        end do
 
        if (yy2(2) < yy2(1)) return
        call search2(dflag,refl,yy2,7,index_ii,frac)
        if (dflag) return
        tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
        tau_x470 = tau_x
        tau_x470_flag = 0
        if (debug) print *, 'aero_470, exit 2358, aot: ', tau_x470
        return
 
650     continue

!
!     Pass the monotonic order check
!
        call search2(dflag,refl,yy,10,index_ii,frac)
        if (dflag) return

!      print *,'after 2nd search'
        tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii) 
        tau_x470 = tau_x
        tau_x470_flag = 0
        if (debug) print *, 'aero_470, exit 2371, aot: ', tau_x470
!      print *,'tau_x470 =', tau_x470
 
        return
      end subroutine aero_470
      
C      subroutine aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
C     1    imod,r470,tau_x470,tau_x470_flag,trflg,model_frac)
C 
C      include 'aottbl.inc'
C      include 'newaottbl.inc'
C      
C      logical, intent(in), optional ::  model_frac
C      
C      real nnvalx(4,4,2,10), yy(10), yy2(8)
C      integer tau_x470_flag, imod
C      logical dflag
C      data pi   /3.14159/
C
C      index_ii = r470
C 
C      frac = (r470-sfc_ref470(index_ii))/
C     1       (sfc_ref470(index_ii+1)-sfc_ref470(index_ii))
C 
C      if (index_ii.lt.1.or.index_ii.gt.24)
C     1    print *,'index_ii = ', index_ii,xlat,xlong
C      if (frac.lt.0.0.or.frac.gt.1.0)
C     1    print *,'aero_470 frac on sfc=', frac, r470, imod
C 
C      call search(dflag,x3,phi,ll,ii)
C      if (dflag) return
C      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))
C 
C      nsm=1
C      dif=x1-theta0(1)
C      do i=1,mm
C        dift=x1-theta0(i)
C        if (dift.gt.0. .and. dift.lt.dif) then
C          nsm=i
C          dif=dift
C        endif
C      enddo
C      mbeg = nsm - 2
C      if (mbeg.le.0) then
C        mbeg = 0
C      else if (mbeg.gt.mm) then
C        mbeg = mm-4
C      endif
C 
C      nsn=1
C      dif=x2-theta(1)
C      do i=1,nn
C        dift=x2-theta(i)
C        if (dift.gt.0. .and. dift.lt.dif) then
C          nsn=i
C          dif=dift
C        endif
C      enddo
C      nbeg = nsn - 2
C      if (nbeg.le.0) then
C        nbeg = 0
C      else if (nbeg.gt.nn) then
C        nbeg = nn-4
C      endif
C 
C      do ia = 1, 10
C       do i = 1, 4
C        do j = 1, 4
C       nnvalx(i,j,1,ia) = nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*
C     1  (1.-frac) + nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac
C       nnvalx(i,j,2,ia) = nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*
C     1  (1.-frac) + nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac
C        enddo
C       enddo
C      enddo
C 
Cc---     interpolating AOT tables
C
C      do 105 ia = 1, 10
C 
C      call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,
C     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
C
C      yy(ia) = y/pi
Cc      print *,'tau, i/f=', tau(ia), y/pi, dy
C 105  continue
C      if (refl.le.yy(1)) then
C       tau_x470 = 0.06
C       if (trflg.gt.0.0) tau_x470 = 0.02
C       tau_x470_flag = -10
C       return
C      endif
C      if (refl.ge.yy(10)) then
C       xxrat = 0.8
C       tau_x470_flag = 1
C       return
C      endif
Cc
Cc
Cc     Check if the reflectance increase with AOT
Cc
C 
C      if (yy(1).lt.yy(2)) go to 650
C 
C      if (refl.lt.yy(4)) return
C 
C      do i = 1, 7
C      yy2(i) = yy(i+3)
C      enddo
C 
C      if (yy2(2).lt.yy2(1)) return
C      call search2(dflag,refl,yy2,7,index_ii,frac)
C      if (dflag) return
C      tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
C      tau_x470 = tau_x
C      return
C 
C650   continue
Cc
Cc     Pass the monotonic order check
Cc
C      call search2(dflag,refl,yy,10,index_ii,frac)
C      if (dflag) return
Cc      print *,'after 2nd search'
C      tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)
C 
C      tau_x470 = tau_x
C
Cc      print *,'tau_x470 =', tau_x470
C      return
C      end

c--------------------------------------------------------
! TODO: check/use/add status flag      
      subroutine aero_650(dflag,refl,x1,x2,x3,mm,nn,ll,ma,r650,tau_x650,  
     1                    tau_x650_flag,tau_x470_flag,tau_x412, tau_x470, 
     1                    tau_x412_flag_91,trflg)
        implicit none
        
        include 'aottbl.inc'
        include 'newaottbl.inc'
        
        logical, intent(inout)        ::  dflag
        real, intent(in)              ::  refl
        real, intent(in)              ::  x1
        real, intent(in)              ::  x2
        real, intent(in)              ::  x3
        integer, intent(in)           ::  mm
        integer, intent(in)           ::  nn
        integer, intent(in)           ::  ll
        integer, intent(in)           ::  ma
        real, intent(in)              ::  r650
        real, intent(inout)           ::  tau_x650
        integer, intent(inout)        ::  tau_x650_flag
        integer, intent(in)           ::  tau_x470_flag
        real, intent(in)              ::  tau_x412
        real, intent(in)              ::  tau_x470
        integer, intent(in)           ::  tau_x412_flag_91
        real, intent(in)              ::  trflg
        
        real, dimension(4,4,2,10)     ::  nnvalx
        real, dimension(10)           ::  yy
        real, dimension(8)            ::  yy2
        real, dimension(4)            ::  yy3
        real, dimension(6)            ::  yy5
     
        real, parameter               ::  pi = 3.14159

        real                          ::  frac, xfrac, dif, dift
        real                          ::  y, dy, w0_x, tau_x
        integer                       ::  index_ii, ii, nsm, nsn, mbeg, nbeg
        integer                       ::  ia,i, j
        
        tau_x650_flag = -999
             
        index_ii = (r650 + 1.0) / 2.0
 
        frac = (r650-sfc_ref650(index_ii)) / (sfc_ref650(index_ii+1)-sfc_ref650(index_ii))
 
        if (index_ii < 1.or.index_ii > 24) print *,'index_ii = ', index_ii
        if (frac < 0.0.or.frac > 1.0) print *,'aero_650 frac on sfc=', frac
 
        if (index_ii < 1) then
          index_ii = 1
          frac = 0.0
        end if
 
        call search(dflag,x3,phi,ll,ii)
        if (dflag) return
        xfrac = (x3 - phi(ii)) / (phi(ii + 1) - phi(ii))
 
        nsm = 1
        dif = x1 - theta0(1)
        do i = 1, mm
          dift = x1 - theta0(i)
          if (dift > 0.0 .and. dift < dif) then
            nsm = i
            dif = dift
          end if
        end do
        mbeg = nsm - 2
        if (mbeg <= 0) then
          mbeg = 0
        else if (mbeg > mm-4) then
          mbeg = mm - 4
        end if
 
        nsn = 1
        dif = x2 - theta(1)
        do i = 1, nn
          dift = x2 - theta(i)
          if (dift > 0. .and. dift < dif) then
            nsn=i
            dif=dift
          end if
        end do
        nbeg = nsn - 2
        if (nbeg <= 0) then
          nbeg = 0
        else if (nbeg > nn-4) then
          nbeg = nn-4
        end if
 
        do ia = 1, 10
          do i = 1, 4
            do j = 1, 4
              nnvalx(i,j,1,ia) = nvalx650(mbeg+i,nbeg+j,ii,ia,index_ii)*    
     1          (1.-frac) + nvalx650(mbeg+i,nbeg+j,ii,ia,index_ii+1)*frac   
              nnvalx(i,j,2,ia) = nvalx650(mbeg+i,nbeg+j,ii+1,ia,index_ii)*  
     1          (1.-frac) + nvalx650(mbeg+i,nbeg+j,ii+1,ia,index_ii+1)*frac   
            end do
          end do
        end do
 
!---     interpolating AOT tables
 
      do ia = 1, 10
        call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,  
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
        yy(ia) = y/pi
      end do
      
      if (refl <= yy(1) .and. yy(1) < yy(2)) then
        tau_x650 = 0.002
        if (trflg > 0.0) tau_x650 = 0.02
        tau_x650_flag = -10
        return
      end if
 
! Reflc off the charts!  Set AOT to max and set flag.
      if (refl >= yy(10)) then
        tau_x650  = 3.5
        w0_x      = -999.
        tau_x650_flag = 1 
        return
      end if

!     -- check for absorbing dust over bright surface
      if (refl >= yy(7)) then
         do i = 1, 4
         yy3(i) = yy(i+6)
         enddo
 
         if (yy3(2) < yy3(1)) return
         call search2(dflag,refl,yy3,4,index_ii,frac)
         if (dflag) return
         tau_x = frac*tau(index_ii+1+6) + (1.-frac)*tau(index_ii+6)
         tau_x650 = tau_x
         tau_x650_flag = 0
         return
      endif
      
!     -- check if 470 or 412 retrievals were off the charts (refl > yy(10))
!      if (tau_x470_flag > 0) go to 670
!      if (tau_x412_flag_91 > 0) go to 680

!
!     Check if the reflectance increase with AOT
!
      if (yy(1) < yy(2)) go to 650
 
      if (refl < yy(4)) return
 
      do i = 1, 7
        yy2(i) = yy(i+3)
      end do
 
      if (yy2(2) < yy2(1)) return
      call search2(dflag,refl,yy2,7,index_ii,frac)
      if (dflag) return
      tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
      tau_x650 = tau_x
      tau_x650_flag = 0
      return

670   continue
c      if (refl < yy(8)) return
c
c      do i = 1, 3
c        yy3(i) = yy(i+7)
c      end do
c
c      if (yy3(2) < yy3(1)) return
c      call search2(dflag,refl,yy3,3,index_ii,frac)
c      if (dflag) return
c      tau_x = frac*tau(index_ii+1+7) + (1.-frac)*tau(index_ii+7)
c      tau_x650 = tau_x
c      tau_x650_flag = 0
c      return

680   continue
      if (refl < yy(5)) return

      do i = 1, 6
        yy5(i) = yy(i+4)
      end do

      if (yy5(2) < yy5(1)) return
      call search2(dflag,refl,yy5,6,index_ii,frac)
      if (dflag) return
!      print *,'after 1st search'
!      if (yy3(1) > yy3(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+4) + (1.-frac)*tau(index_ii+4)
      tau_x650 = tau_x
      tau_x650_flag = 0
      return
 
650   continue
!
!     Pass the monotonic order check
!
      call search2(dflag,refl,yy,10,index_ii,frac)
      if (dflag) return
!      print *,'after 2nd search'
      tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)
 
      tau_x650 = tau_x
      tau_x650_flag = 0
      
      return
      
      end subroutine aero_650


C      subroutine aero_650(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
C     1    r650,tau_x650,tau_x650_flag,tau_x470_flag,
C     2    tau_x412, tau_x470,tau_x412_flag_91,trflg)
C 
C      include 'aottbl.inc'
C      include 'newaottbl.inc'
C
C      real nnvalx(4,4,2,10), yy(10), yy2(8), yy3(3), yy5(6)
C      real tau_x650, tau_x412, tau_x470
C      real refl,x1,x2,x3,r650
C
C      integer tau_x650_flag, tau_x470_flag, tau_x412_flag_91
C      integer mm,nn,ll,ma
C      logical dflag
C      data pi   /3.14159/
C
C      index_ii = (r650+1.)/2.
C 
C      frac = (r650-sfc_ref650(index_ii))/
C     1       (sfc_ref650(index_ii+1)-sfc_ref650(index_ii))
C 
C      if (index_ii.lt.1.or.index_ii.gt.24) 
C     1    print *,'index_ii = ', index_ii,xlat,xlong
C      if (frac.lt.0.0.or.frac.gt.1.0) 
C     1    print *,'aero_650 frac on sfc=', frac
C 
C      if (index_ii.lt.1) then
C        index_ii = 1
C        frac = 0.
C      endif
C 
C      call search(dflag,x3,phi,ll,ii)
C      if (dflag) return
C      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))
C 
C      nsm=1
C      dif=x1-theta0(1)
C      do i=1,mm
C        dift=x1-theta0(i)
C        if (dift.gt.0. .and. dift.lt.dif) then
C          nsm=i
C          dif=dift
C        endif
C      enddo
C      mbeg = nsm - 2
C      if (mbeg.le.0) then
C        mbeg = 0
C      else if (mbeg.gt.mm) then
C        mbeg = mm-4
C      endif
C 
C      nsn=1
C      dif=x2-theta(1)
C      do i=1,nn
C        dift=x2-theta(i)
C        if (dift.gt.0. .and. dift.lt.dif) then
C          nsn=i
C          dif=dift
C        endif
C      enddo
C      nbeg = nsn - 2
C      if (nbeg.le.0) then
C        nbeg = 0
C      else if (nbeg.gt.nn) then
C        nbeg = nn-4
C      endif
C 
C      do ia = 1, 10
C       do i = 1, 4
C        do j = 1, 4
C          nnvalx(i,j,1,ia) = nvalx650(mbeg+i,nbeg+j,ii,ia,index_ii)*
C     1       (1.-frac) + nvalx650(mbeg+i,nbeg+j,ii,ia,index_ii+1)*frac   
C          nnvalx(i,j,2,ia) = nvalx650(mbeg+i,nbeg+j,ii+1,ia,index_ii)*
C     1       (1.-frac) + nvalx650(mbeg+i,nbeg+j,ii+1,ia,index_ii+1)*frac   
C        enddo
C       enddo
C      enddo
C 
Cc---     interpolating AOT tables
C 
C      do 600 ia = 1, 10
C 
C      call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,
C     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
C 
C      yy(ia) = y/pi
C 600  continue
C 
C      if (refl.le.yy(1).and.yy(1).lt.yy(2)) then
C       tau_x650 = 0.06
C       if (trflg.gt.0.0) tau_x650 = 0.02
C       return
C      endif
C 
C      if (refl.ge.yy(10)) then
C         tau_x650 = 4.0
C         w0_x = -999.
C         tau_x650_flag = 1 
C         return
C      endif
C
C      if (tau_x470_flag.gt.0) go to 670
C      if (tau_x412_flag_91.gt.0) go to 680
C
Cc
Cc     Check if the reflectance increase with AOT
Cc
C      if (yy(1).lt.yy(2)) go to 650
C 
C      if (refl.lt.yy(4)) return
C 
C      do i = 1, 7
C      yy2(i) = yy(i+3)
C      enddo
C 
C      if (yy2(2).lt.yy2(1)) return
C      call search2(dflag,refl,yy2,7,index_ii,frac)
C      if (dflag) return
C      tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
C      tau_x650 = tau_x
C      return
C
C670   continue
C      if (refl.lt.yy(8)) return
C
C      do i = 1, 3
C      yy3(i) = yy(i+7)
C      enddo
C
C      if (yy3(2).lt.yy3(1)) return
C      call search2(dflag,refl,yy3,3,index_ii,frac)
C      if (dflag) return
C      tau_x = frac*tau(index_ii+1+7) + (1.-frac)*tau(index_ii+7)
C      tau_x650 = tau_x
C      return
C
C680   continue
C      if (refl.lt.yy(5)) return
C
C      do i = 1, 6
C      yy5(i) = yy(i+4)
C      enddo
C
C      if (yy5(2).lt.yy5(1)) return
C      call search2(dflag,refl,yy5,6,index_ii,frac)
C      if (dflag) return
Cc      print *,'after 1st search'
Cc      if (yy3(1).gt.yy3(2)) print *,'yy=',refl,(yy(i),I=1,10)
C      tau_x = frac*tau(index_ii+1+4) + (1.-frac)*tau(index_ii+4)
C      tau_x650 = tau_x
C      return
C 
C650   continue
Cc
Cc     Pass the monotonic order check
Cc
C      call search2(dflag,refl,yy,10,index_ii,frac)
C      if (dflag) return
Cc      print *,'after 2nd search'
C      tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)
C 
C      tau_x650 = tau_x
C
C      return
C      end

c--------------------------------------------------------
      subroutine aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,imod,r412, 
     1                    tau_x412,tau_x412_flag,trflg,model_frac,debug)
        
        implicit none
        
        include 'aottbl.inc'
        include 'newaottbl.inc'

        logical, intent(inout)      ::  dflag
        real, intent(in)            ::  refl
        real, intent(in)            ::  x1
        real, intent(in)            ::  x2
        real, intent(in)            ::  x3
        integer, intent(in)         ::  mm
        integer, intent(in)         ::  nn
        integer, intent(in)         ::  ll
        integer, intent(in)         ::  ma
        integer, intent(in)         ::  imod
        real, intent(in)            ::  r412
        real, intent(inout)         ::  tau_x412
        integer, intent(inout)      ::  tau_x412_flag
        real, intent(in)            ::  trflg
        real, intent(in)            ::  model_frac
        logical, intent(in)         ::  debug
        
!       -- NOTE: if model_frac = 0.0, the aerosol model = imod
!       -- if model_frac > 0.0, the aerosol model will be interpolated between
!       --  imod and imod + 1, using the value of model_frac

        real, dimension(4,4,2,10) ::  nnvalx, nnvalx1, nnvalx2 
        real, dimension(10)       ::  yy
        real, dimension(8)        ::  yy2
        
        real, parameter             ::  pi = 3.14159
      
        integer                     ::  index_ii, ii, nsm, nsn, mbeg, nbeg
        integer                     ::  ia, j, i
        real                        ::  frac, xfrac, dif, dift, y, dy
        real                        ::  xxrat, tau_x
                
        if (debug) print *, 'aero_412, in: ', refl, x1, x2, x3, r412, imod, model_frac

        tau_x412 = -999.0
        tau_x412_flag = -999
                
        index_ii = r412
        if (r412 < 0.0) return
        frac = (r412-sfc_ref412(index_ii)) / (sfc_ref412(index_ii+1)-sfc_ref412(index_ii))
        if (debug) print *, 'aero_412, sfc indx: ', r412, index_ii
        if (index_ii < 1 .or. index_ii > 20) then
          print *,'index_ii = ', index_ii
          dflag = .false.
          return
        end if
        
        if (frac < 0.0 .or. frac > 1.0) then
          print *,'aero_412 frac on sfc=', frac, r412, sfc_ref412(index_ii)
          dflag = .false.
          return
        end if
      
        call search(dflag,x3,phi,ll,ii)
        if (dflag) return
        xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))
        if (debug) print *, 'aero_412, raa: ', x3, ii, phi(ii), phi(ii+1)

        nsm = 1
        mbeg = -999
        dif = x1 - theta0(1)
        do i = 1, mm
          dift = x1 - theta0(i)
          if (dift > 0.0 .and. dift < dif) then
            nsm = i
            dif = dift
          end if
        end do
        mbeg = nsm - 2
        if (mbeg <= 0) then
          mbeg = 0
        else if (mbeg > mm-4) then
          mbeg = mm-4
        end if
 
        nsn = 1
        nbeg = -999
        dif = x2 - theta(1)
        do i = 1, nn
          dift = x2 - theta(i)
          if (dift > 0.0 .and. dift < dif) then
            nsn = i
            dif = dift
          end if
        end do
        nbeg = nsn - 2
        if (nbeg <= 0) then
          nbeg = 0
        else if (nbeg > nn-4) then
          nbeg = nn-4
        end if
        if (debug) print *, 'aero_412, mbeg, nbeg: ', mbeg, nbeg

!       -- interpolate between models if requested.
!        if (present(model_frac)) then
        if (imod < 8) then              ! temp check to avoid out-of-bounds issue with imod=8
          nnvalx1(:,:,:,:) = -999.0
          nnvalx2(:,:,:,:) = -999.0
          do ia = 1, 10
            do i = 1, 4
              do j = 1, 4
                nnvalx1(i,j,1,ia) = nvalx412(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*
     1            (1.0-frac) + nvalx412(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac
                nnvalx1(i,j,2,ia) = nvalx412(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*
     1             (1.0-frac) + nvalx412(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac
  
                nnvalx2(i,j,1,ia) = nvalx412(mbeg+i,nbeg+j,ii,ia,imod+1,index_ii)*
     1             (1.0-frac) + nvalx412(mbeg+i,nbeg+j,ii,ia,imod+1,index_ii+1)*frac
                nnvalx2(i,j,2,ia) = nvalx412(mbeg+i,nbeg+j,ii+1,ia,imod+1,index_ii)*
     1             (1.0-frac) + nvalx412(mbeg+i,nbeg+j,ii+1,ia,imod+1,index_ii+1)*frac
  
                nnvalx(i,j,1,ia) = (1.0-model_frac) * nnvalx1(i,j,1,ia) +   
     1                                 model_frac * nnvalx2(i,j,1,ia)
                nnvalx(i,j,2,ia) = (1.0-model_frac) * nnvalx1(i,j,2,ia) +   
     1                              model_frac * nnvalx2(i,j,2,ia)
              end do
            end do
          end do
          
!       -- just use imod, no interpolation
        else      
          nnvalx(:,:,:,:) = -999.0
          do ia = 1, 10
            do i = 1, 4
              do j = 1, 4
                nnvalx(i,j,1,ia) = nvalx412(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*     
     1              (1.-frac) + nvalx412(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac   
                nnvalx(i,j,2,ia) = nvalx412(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*   
     1              (1.-frac) + nvalx412(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac   
              end do
            end do
          end do
        end if
        
!---     interpolating AOT tables
        yy(:) = -999.0
        do ia = 1, 10
!          if (debug) print *, mm, nn, ll, ia, x1, x2, x3, mbeg, nbeg, xfrac
!          if (debug) print *, nnvalx
          call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,    
     1                   x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
          yy(ia) = y/pi
!          if (debug) print *,'tau, i/f=', tau(ia), y/pi, dy
        end do
!        print *,'refl1=', refl
      
        if (refl <= yy(1)) then
          tau_x412 = 0.06 
          if (trflg > 0.0) tau_x412 = 0.02
          tau_x412_flag = -10       
          if (debug) print *, 'aero_412, hit low bound: ', refl, yy(1)
          return
        end if

! Reflc off the charts!  Set AOT to max and set flag.
        if (refl >= yy(10)) then
          xxrat = 0.8
          tau_x412 = 3.5
          tau_x412_flag = 1
          if (debug) print *, 'aero_412, hit hi bound: ', refl, yy(10)
          return
        end if
!
!     Check if the reflectance increase with AOT
!
        if (yy(1) < yy(2)) go to 650
 
        if (refl < yy(4)) return
        
        yy2(:) = -999.0
        do i = 1, 7
          yy2(i) = yy(i+3)
        end do
 
        if (yy2(2) < yy2(1)) return
        call search2(dflag,refl,yy2,7,index_ii,frac)
        if (dflag) return
        tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
        tau_x412 = tau_x
        tau_x412_flag = 0
        if (debug) print *, 'aero_412, exit 2355, aot: ', tau_x412
        return
 
650     continue
!
!     Pass the monotonic order check
!
        call search2(dflag,refl,yy,10,index_ii,frac)
        if (dflag) return
        tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)
 
        tau_x412 = tau_x
        tau_x412_flag = 0
        if (debug) print *, 'aero_412, exit 2367, aot: ', tau_x412
        return
 
      end subroutine aero_412
      
C      subroutine aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
C     1    imod,r412,tau_x412,tau_x412_flag,trflg,debug)
C 
C      include 'aottbl.inc'
C      include 'newaottbl.inc'
C      
Cc      logical, intent(in), optional ::  model_frac
C      logical, intent(in) ::  debug
C      
C      real nnvalx(4,4,2,10), yy(10), yy2(8)
C      integer tau_x412_flag
C      logical dflag
C      data pi   /3.14159/
C 
C      if (debug) print *, 'aero_412_old, in: ', refl, x1, x2, x3, r412, w0_x
C
C      index_ii = r412
C 
C      frac = (r412-sfc_ref412(index_ii))/
C     1       (sfc_ref412(index_ii+1)-sfc_ref412(index_ii))
C 
C      if (index_ii.lt.1.or.index_ii.gt.20)
C     1    print *,'index_ii = ', index_ii
C      if (frac.lt.0.0.or.frac.gt.1.0)
C     1    print *,'aero_412 frac on sfc=', frac
C      
C      if (debug) print *, 'aero_412_old, sfc indx: ', r412, index_ii
C      
C      call search(dflag,x3,phi,ll,ii)
C      if (dflag) return
C      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))
C      if (debug) print *, 'aero_412_old, raa: ', x3, ii, phi(ii), phi(ii+1)
C      
C      nsm=1
C      dif=x1-theta0(1)
C      do i=1,mm
C        dift=x1-theta0(i)
C        if (dift.gt.0. .and. dift.lt.dif) then
C          nsm=i
C          dif=dift
C        endif
C      enddo
C      mbeg = nsm - 2
C      if (mbeg.le.0) then
C        mbeg = 0
C      else if (mbeg.gt.mm-4) then
C        mbeg = mm-4
C      endif
C 
C      nsn=1
C      dif=x2-theta(1)
C      do i=1,nn
C        dift=x2-theta(i)
C        if (dift.gt.0. .and. dift.lt.dif) then
C          nsn=i
C          dif=dift
C        endif
C      enddo
C      nbeg = nsn - 2
C      if (nbeg.le.0) then
C        nbeg = 0
C      else if (nbeg.gt.nn-4) then
C        nbeg = nn-4
C      endif
C      if (debug) print *, 'aero_412_old, mbeg, nbeg: ', mbeg, nbeg
C      
C      do ia = 1, 10
C       do i = 1, 4
C        do j = 1, 4
C        nnvalx(i,j,1,ia) = nvalx412(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*
C     1   (1.-frac) + nvalx412(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac   
C       nnvalx(i,j,2,ia) = nvalx412(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*
C     1  (1.-frac) + nvalx412(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac   
C        enddo
C       enddo
C      enddo
C 
Cc---     interpolating AOT tables
C 
C      do 100 ia = 1, 10
C      if (debug) print *, mm, nn, ll, ia, x1, x2, x3, mbeg, nbeg, xfrac
C      call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,
C     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
C      if (debug) print *,'tau, i/f=', tau(ia), y/pi, dy
C      yy(ia) = y/pi
C 100  continue
C 
C      if (refl.le.yy(1)) then
C       tau_x412 = 0.06
C       if (trflg.gt.0.0) tau_x412 = 0.02
C       tau_x412_flag = -10
C       if (debug) print *, 'aero_412_old, hit low bound: ', refl, yy(1)
C       return
C      endif
C 
C      if (refl.ge.yy(10)) then
C       xxrat = 0.8
C       tau_x412_flag = 1
C       if (debug) print *, 'aero_412_old, hit hi bound: ', refl, yy(10)
C       return
C      endif
Cc
Cc     Check if the reflectance increase with AOT
Cc
C      if (yy(1).lt.yy(2)) go to 650
C 
C      if (refl.lt.yy(4)) return
C 
C      do i = 1, 7
C      yy2(i) = yy(i+3)
C      enddo
C 
C      if (yy2(2).lt.yy2(1)) return
C      call search2(dflag,refl,yy2,7,index_ii,frac)
C      if (dflag) return
C      tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
C      tau_x412 = tau_x
C      if (debug) print *, 'aero_412_old, exit 2482, aot: ', tau_x412
C      return
C 
C650   continue
Cc
Cc     Pass the monotonic order check
Cc
C      call search2(dflag,refl,yy,10,index_ii,frac)
C      if (dflag) return
C      tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)
C 
C      tau_x412 = tau_x
C      if (debug) print *, 'aero_412_old, exit 2494, aot: ', tau_x412
C      return
C      end

c------------------------------------------------------------
      subroutine aero_mod (tau_x412,tau_x470,tau_x650,
     1           tau_x412_91,aot_mod)

      real aot_mod(6)

      aot_mod(1) = tau_x650
      aot_mod(2) = tau_x470
      aot_mod(3) = tau_x412_91
      aot_mod(4) = tau_x412*1.9
      aot_mod(5) = tau_x470*2.
      aot_mod(6) = tau_x650*2.8
      return
      end
      
			integer function handle_lut_out_of_bounds(geo_zone, retr_flag, aot) result(status)
				implicit none
				
				integer, intent(in)	:: geo_zone
				integer, intent(in)	:: retr_flag
				real, intent(inout)	:: aot
				
				status = 0
												
				if (retr_flag == -10) then      ! -10 = return value from aero_412, etc indicating
					select case (geo_zone)        ! refl < yy(1) scenario.
						case (5:11, 15, 16, 19, 20, 21, 23, 24, 26, 27)
							aot = 0.06
						case (1:4, 12:14, 17, 18, 22, 25)
							aot = 0.02
						case (-999)
							aot = 0.02
						case default
							print *, "ERROR: Invalid geographical zone in handle_lut_out_of_bounds: ", geo_zone
							status = -1
							return
					end select
				end if
				
				return
				
			end function handle_lut_out_of_bounds
			
			
			subroutine smoke_mod(view,scat_ang,dd412,dd470,dd650,  
     1                     tau_x412ss,tau_x412ss_91,tau_smoke)
        implicit none
        
        real, intent(in)        ::  view
        real, intent(in)        ::  scat_ang
        real, intent(in)        ::  dd412
        real, intent(in)        ::  dd470
        real, intent(in)        ::  dd650
        real, intent(in)        ::  tau_x412ss, tau_x412ss_91
        real                    ::  tau_smoke
     
        if (dd412 <1.0 .and. dd470 <1.0) then
            tau_smoke = tau_x412ss * 2.85  
          if (dd470 > dd412 .and. dd650 > dd412)   
     1      tau_smoke = tau_x412ss_91 * 1.77
          if (abs(dd412-dd470) < 0.25) then
            if (dd412 <0.35) tau_smoke = tau_x412ss * 4.3
            if (scat_ang <158.0 .and. dd412 <0.8) then
            tau_smoke = tau_x412ss * 3.5
            if (view >45.0 .and. dd470 <0.2 .and. abs(dd412-dd470) < 0.24)  
     1        tau_smoke = tau_x412ss * 6.0
            endif
          endif
          if (abs(dd412-dd650) < 0.1 .and. dd412 >0.7)    
     1        tau_smoke = tau_x412ss_91
          if (dd650 > 1.0 .and. abs(dd412-dd470) < 0.25) then
              tau_smoke = tau_x412ss * 4.5
          if (scat_ang <158.0 .and. dd412 >0.8)    
     1        tau_smoke = tau_x412ss * 3.5
          endif
        endif
     
        if (dd412 >1.0 .and. dd470 <1.0) then
             tau_smoke = tau_x412ss * 2.7
          if (abs(dd412-dd470) < 0.25 .and. dd650 <0.)   
     1      tau_smoke = tau_x412ss * 2.4
        endif
     
        if (dd412 <1.0 .and. dd470 >1.0) then
             tau_smoke = tau_x412ss * 2.7
          if (abs(dd412-dd470) < 0.25 .and. dd650 <0.)    
     1        tau_smoke = tau_x412ss * 2.4
          if (dd650 - dd470 > -0.2) then
             tau_smoke = tau_x412ss * 2.1
          if (view < 45.0) tau_smoke = tau_x412ss_91
          endif
        endif
     
        if (dd412 >1.0 .and. dd470 >1.0) then
          if (dd412 > dd470) tau_smoke = tau_x412ss_91 * 1.15
          if (dd470 > dd650 .and. dd650 > dd412 .and. dd650 >1.0) then
            if (view >= 45.0) tau_smoke = tau_x412ss * 2.1
            if (view <  45.0) tau_smoke = tau_x412ss_91
          endif
          if (dd412 <2.0 .and. dd470 >2.0)   
     1      tau_smoke = tau_x412ss * 1.4
          if (dd650 > dd470 .and. dd470 > dd412) then
            if (view <= 25.0) tau_smoke = (tau_x412ss+tau_x412ss_91)/2.
            if (view > 25.0 .and. abs(dd412-dd470) < 0.2)  
     1        tau_smoke = (tau_x412ss+tau_x412ss_91)/2.
          endif
          if (dd650 > 3.5) tau_smoke = tau_x412ss
        endif
     
        if (dd650 <-1.0 .and. abs(dd412-dd650)>1.9 .and. view >45.) then
          if (dd412 >= 0.88) tau_smoke = tau_x412ss * 4.0
          if (dd412 <0.88 .and. dd412 >= 0.6) tau_smoke = tau_x412ss *6.0
          if (dd412 <0.6 .and. dd412 >= 0.35) tau_smoke = tau_x412ss *8.0
          if (dd412 <0.35) tau_smoke = tau_x412ss_91 *8.0
          if (abs(dd412-dd470) < 0.14 .and. dd650 > -1.7 .and. dd412 <0.4)  
     1        tau_smoke = tau_x412ss * 3.5
          if (dd412 <0.) tau_smoke = tau_x412ss_91 *9.0
        endif
     
        return
 
      end subroutine smoke_mod
			
