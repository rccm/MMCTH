c-- rsfc470 calculates a surface reflectance values including the effects of the 
c-- BRDF of the surface.
c-- mod_sfc is an index value designating a particular region of the world.
c------------------------------------------------------------
      subroutine rsfc470(dflag2,mod_sfc,xday,xthet,xphi,
     1                   scat_ang,r470)
c
      integer mod_sfc, nc_wav(6), index_iw
      logical dflag2
      real xday,r470p,xsfc470p,frac_iw
      real xthet,xphi,scat_ang,r470,xthet5
      real wav1(16),sfc1(16),wav2(16),sfc2(16)
      real wav3(16),sfc3(16),wav4(15),sfc4(15)
      real wav5(15),sfc5(15),wav6(15),sfc6(15)
      character *60 xname_wav(6)
      data xname_wav /'Harmim','Mezaira','Al_Khaznah',
     1  'Sai_Salam','SMART','Bright'/
      data nc_wav /16,16,16,15,15,15/
      data wav1 /-62.2694,-57.1621,-51.0789,-43.4577,
     1 -34.2301,-23.0985,-10.1670,-9.66400,3.82708,17.4518,
     1  29.4005, 39.5023, 47.6424, 54.2590,59.7253,64.2073/
      data sfc1 /12.8768,13.9229,13.7663,13.6870,13.2900,
     1   13.3372,13.1316,13.1165,12.3269,12.3000,12.1731,
     1   12.3783,12.0233,12.2451,11.9212,11.5149/
      data wav2 /-65.2144,-60.7404,-55.3780,-48.8869,-40.6035,
     1 -30.8647,-19.0751,-5.59166,8.29930,21.4201,32.7246,
     1  42.2137, 49.7960, 56.0005,61.1576,65.3304/
      data sfc2 /14.1221,14.3341,14.4184,14.8820,14.3835,
     1   14.0886,14.1124,13.7649,12.7891,12.9838,12.8853,
     1   12.9273,12.9842,12.9033,12.2712,11.6638/
      data wav3 /-65.0665,-60.5814,-55.2427,-48.7424,-40.7624,
     1 -30.8889,-19.1950,-5.84152,7.86039,21.0191,32.2845,
     1  41.6669, 49.4183, 55.6076,60.8122,65.1024/
      data sfc3 /14.5394,15.4657,16.0964,16.1359,15.7386,
     1   16.2421,15.8848,14.4071,14.5017,14.0462,14.2213,
     1   13.5992,13.6395,13.2211,12.9375,12.6472/
      data wav4 /-61.5411,-56.4792,-50.1862,-42.6219,-33.1951,
     1 -21.9909, -9.07847, 4.54637, 17.8186, 29.6754, 39.4847,
     1  47.5092,  54.1248, 59.5025, 63.9649/
      data sfc4 /17.9199,18.6786,18.6900,18.2124,18.1539,
     1   17.6234,16.5130,16.6544,16.2219,16.0444,16.2370,
     1   15.9920,15.7644,15.4508,16.0762/
      data wav5 /-62.1827,-57.1542,-51.0592,-43.6501,-34.5908,
     1 -23.3993, -10.6227, 3.20328, 16.5169, 28.6361, 38.6162,
     1  46.9379,  53.6876, 59.1098, 63.7301/
      data sfc5 /11.9229,12.3774,12.0428,11.1172,10.4741,
     1   10.5246,10.1789,9.89490,9.38978,9.49380,9.74338,
     1   9.96068,10.0407,9.62602,9.57829/
      data wav6 /-61.5411,-56.4792,-50.1862,-42.6219,-33.1951,
     1 -21.9909, -9.07847, 4.54637, 17.8186, 29.6754, 39.4847,
     1  47.5092,  54.1248, 59.5025, 63.9649/
      data sfc6 /19.9199,20.6786,20.6900,20.2124,20.1539,
     1   17.6234,15.5130,15.6544,15.2219,15.0444,16.2370,
     1   15.9920,15.7644,15.4508,16.0762/

      dflag2 = .false.
       
      if (mod_sfc.eq.8) then 
c---     Saada
c      --MAM
      r470p = 7.69
      if (xphi.le.90.0) then
       if (scat_ang.gt.150.0)
     1     r470 = 5.5
       if (scat_ang.gt.141.0.and.scat_ang.le.150.)
     1     r470 = 4.7 +0.8*(scat_ang-141.)/9.
       if (scat_ang.gt.131.0.and.scat_ang.le.141.)
     1     r470 = 4.7
       if (scat_ang.gt.111.0.and.scat_ang.le.131.)
     1     r470 = 3.2 +1.5*(scat_ang-111.)/20.
       if (scat_ang.gt.101.0.and.scat_ang.le.111.)
     1     r470 = 3.2 +0.0*(scat_ang-101.)/10.
       if (scat_ang.gt.99.0.and.scat_ang.le.101.)
     1     r470 = 2. +1.2*(scat_ang-99.)/2.
       if (scat_ang.le.99.0) r470 = 2.0
      endif
      r470p = 7.69
      if (xphi.gt.90.0) then
       if (xthet.le.30.0) then
       if (scat_ang.gt.168.0) r470 = 6.7
       if (scat_ang.gt.130.0.and.scat_ang.le.168.)
     1     r470 = 4.0 +2.7*(scat_ang-130.)/38.
      if (scat_ang.le.130.0) r470 = 4.0
       endif
       if (xthet.gt.30.0) then
       if (scat_ang.gt.168.0) r470 = 7.7
       if (scat_ang.gt.164.0.and.scat_ang.le.168.)
     1     r470 = 6.9 +0.8*(scat_ang-164.)/4.
       if (scat_ang.gt.150.0.and.scat_ang.le.164.)
     1     r470 = 5.6 +1.3*(scat_ang-150.)/14.
       if (scat_ang.gt.143.0.and.scat_ang.le.150.)
     1     r470 = 5.4 +0.2*(scat_ang-143.)/7.
       if (scat_ang.gt.140.0.and.scat_ang.le.143.)
     1     r470 = 3.5 +1.9*(scat_ang-140.)/3.
       if (scat_ang.le.140.0) r470 = 3.5
       endif
       if (xthet.gt.50.0.and.xthet.lt.59.0) then
       if (scat_ang.gt.143.0.and.scat_ang.le.150.)
     1     r470 = 6.0 -0.4*(scat_ang-143.)/7.
       if (scat_ang.gt.140.0.and.scat_ang.le.143.)
     1     r470 = 3.5 +2.5*(scat_ang-140.)/3.
       endif
      endif
 
c      --JAS
      if (xday.ge.145.0.and.xday.lt.274.0) then
      if (xphi.le.90.0) then
       if (scat_ang.gt.150.0)
     1     r470 = 6.8
       if (scat_ang.gt.141.0.and.scat_ang.le.150.)
     1     r470 = 6.2 +0.6*(scat_ang-141.)/9.
       if (scat_ang.gt.131.0.and.scat_ang.le.141.)
     1     r470 = 5.8 +0.4*(scat_ang-131.)/10.
       if (scat_ang.gt.114.0.and.scat_ang.le.131.)
     1     r470 = 5.2 +0.6*(scat_ang-114.)/17.
       if (scat_ang.gt.111.0.and.scat_ang.le.114.)
     1     r470 = 5.6 -0.4*(scat_ang-111.)/3.
       if (scat_ang.gt.100.0.and.scat_ang.le.111.)
     1     r470 = 4.5 +1.1*(scat_ang-100.)/11.
       if (scat_ang.le.100.0) r470 = 4.5
      endif
      if (xphi.gt.90.0) then
       if (scat_ang.gt.168.0) r470 = 8.2
       if (scat_ang.gt.164.0.and.scat_ang.le.168.)
     1     r470 = 8.1 +0.1*(scat_ang-164.)/4.
       if (scat_ang.gt.159.0.and.scat_ang.le.164.)
     1     r470 = 7.9 +0.2*(scat_ang-159.)/5.
       if (scat_ang.gt.155.0.and.scat_ang.le.159.)
     1     r470 = 8.2 -0.3*(scat_ang-155.)/4.
       if (scat_ang.gt.153.0.and.scat_ang.le.155.)
     1     r470 = 7.2 +1.0*(scat_ang-153.)/2.
       if (scat_ang.gt.150.0.and.scat_ang.le.153.)
     1     r470 = 7.2 +0.0*(scat_ang-150.)/3.
       if (scat_ang.gt.146.0.and.scat_ang.le.150.)
     1     r470 = 6.85 +0.35*(scat_ang-146.)/4.
       if (scat_ang.gt.143.0.and.scat_ang.le.146.)
     1     r470 = 6.0 +0.85*(scat_ang-143.)/3.
       if (scat_ang.gt.135.0.and.scat_ang.le.143.)
     1     r470 = 5.8 +0.2*(scat_ang-135.)/8.
       if (scat_ang.gt.130.0.and.scat_ang.le.135.)
     1     r470 = 6.5 -0.7*(scat_ang-130.)/5.
       if (scat_ang.le.130.0) r470 = 6.5
 
       if (xthet.le.35.0) then
       if (scat_ang.gt.168.0) r470 = 8.5
       if (scat_ang.gt.166.0.and.scat_ang.le.168.)
     1     r470 = 8.65 -0.15*(scat_ang-166.)/2.
       if (scat_ang.gt.164.0.and.scat_ang.le.166.)
     1     r470 = 8.5 +0.15*(scat_ang-164.)/2.
       if (scat_ang.gt.155.0.and.scat_ang.le.164.)
     1     r470 = 7.6 +0.9*(scat_ang-155.)/9.
       if (scat_ang.gt.153.0.and.scat_ang.le.155.)
     1     r470 = 7.7 -0.1*(scat_ang-153.)/2.
       if (scat_ang.gt.150.0.and.scat_ang.le.153.)
     1     r470 = 7.5 +0.2*(scat_ang-150.)/3.
       if (scat_ang.gt.143.0.and.scat_ang.le.150.)
     1     r470 = 6.0 +1.5*(scat_ang-143.)/7.
       if (scat_ang.gt.140.0.and.scat_ang.le.143.)
     1     r470 = 5.8 +0.2*(scat_ang-140.)/3.
       if (scat_ang.le.140.0) r470 = 5.8
      endif
      endif
      endif
 
c      --NDJ
      if (xday.ge.274.0.and.xday.lt.360.0) then
      if (xphi.le.90.0) then
       if (scat_ang.gt.150.0)
     1     r470 = 6.8
       if (scat_ang.gt.141.0.and.scat_ang.le.150.)
     1     r470 = 6.7 +0.1*(scat_ang-141.)/9.
       if (scat_ang.gt.131.0.and.scat_ang.le.141.)
     1     r470 = 6.7 -0.0*(scat_ang-131.)/10.
       if (scat_ang.gt.114.0.and.scat_ang.le.131.)
     1     r470 = 5.2 +1.5*(scat_ang-114.)/17.
       if (scat_ang.gt.111.0.and.scat_ang.le.114.)
     1     r470 = 5.6 -0.4*(scat_ang-111.)/3.
       if (scat_ang.gt.100.0.and.scat_ang.le.111.)
     1     r470 = 4.0 +1.6*(scat_ang-100.)/11.
       if (scat_ang.gt.96.0.and.scat_ang.le.100.)
     1     r470 = 4.5 -0.5*(scat_ang-96.)/4.
       if (scat_ang.le.96.0) r470 = 4.5
      endif
      if (xphi.gt.90.0) then
       if (scat_ang.gt.168.0) r470 = 8.2
       if (scat_ang.gt.164.0.and.scat_ang.le.168.)
     1     r470 = 8.1 +0.1*(scat_ang-164.)/4.
       if (scat_ang.gt.159.0.and.scat_ang.le.164.)
     1     r470 = 7.9 +0.2*(scat_ang-159.)/5.
       if (scat_ang.gt.155.0.and.scat_ang.le.159.)
     1     r470 = 8.2 -0.3*(scat_ang-155.)/4.
       if (scat_ang.gt.153.0.and.scat_ang.le.155.)
     1     r470 = 7.4 +0.8*(scat_ang-153.)/2.
       if (scat_ang.gt.150.0.and.scat_ang.le.153.)
     1     r470 = 7.2 +0.2*(scat_ang-150.)/3.
       if (scat_ang.gt.144.0.and.scat_ang.le.150.)
     1     r470 = 6.6 +0.6*(scat_ang-144.)/6.
       if (scat_ang.gt.142.0.and.scat_ang.le.144.)
     1     r470 = 6.9 -0.3*(scat_ang-142.)/2.
       if (scat_ang.gt.138.0.and.scat_ang.le.142.)
     1     r470 = 6.3 +0.6*(scat_ang-138.)/4.
       if (scat_ang.gt.135.0.and.scat_ang.le.138.)
     1     r470 = 5.6 +0.7*(scat_ang-135.)/3.
       if (scat_ang.gt.133.0.and.scat_ang.le.135.)
     1     r470 = 6.0 -0.4*(scat_ang-133.)/2.
       if (scat_ang.gt.132.0.and.scat_ang.le.133.)
     1     r470 = 4.0 +2.0*(scat_ang-132.)
       if (scat_ang.le.132.0) r470 = 4.0
 
       if (xthet.le.30.0) then
       if (scat_ang.gt.168.0) r470 = 8.3
       if (scat_ang.gt.166.0.and.scat_ang.le.168.)
     1     r470 = 8.65 -0.35*(scat_ang-166.)/2.
       if (scat_ang.gt.164.0.and.scat_ang.le.166.)
     1     r470 = 8.5 +0.15*(scat_ang-164.)/2.
       if (scat_ang.gt.155.0.and.scat_ang.le.164.)
     1     r470 = 7.6 +0.9*(scat_ang-155.)/9.
       if (scat_ang.gt.153.0.and.scat_ang.le.155.)
     1     r470 = 7.7 -0.1*(scat_ang-153.)/2.
       if (scat_ang.gt.150.0.and.scat_ang.le.153.)
     1     r470 = 7.5 +0.2*(scat_ang-150.)/3.
       if (scat_ang.gt.143.0.and.scat_ang.le.150.)
     1     r470 = 6.1 +1.4*(scat_ang-143.)/7.
       if (scat_ang.gt.138.0.and.scat_ang.le.143.)
     1     r470 = 6.1 +0.0*(scat_ang-138.)/5.
       if (scat_ang.gt.134.0.and.scat_ang.le.138.)
     1     r470 = 5.3 +0.8*(scat_ang-134.)/4.
       if (scat_ang.le.134.0) r470 = 5.3
       endif
      endif
      endif
      endif

      if (mod_sfc.eq.9) then
      xthet5 = xthet 
      if (xphi.gt.90.0) xthet5 = -1. * xthet
       if (xthet5.le.wav1(1)) then
          r470 = sfc1(1)
          return
       endif
       if (xthet5.ge.wav1(nc_wav(1))) then
          r470 = sfc1(nc_wav(1))
          return
       endif
       call search2(dflag2,xthet5,wav1,nc_wav(1),index_iw,frac_iw)
      if (dflag2) return
       r470 = frac_iw*sfc1(index_iw+1) + (1.-frac_iw)*sfc1(index_iw)
      endif

      if (mod_sfc.eq.10) then
      xthet5 = xthet
      if (xphi.gt.90.0) xthet5 = -1. * xthet
       if (xthet5.le.wav2(1)) then
          r470 = sfc2(1)
          return
       endif
       if (xthet5.ge.wav2(nc_wav(2))) then
          r470 = sfc2(nc_wav(2))
          return
       endif
       call search2(dflag2,xthet5,wav2,nc_wav(2),index_iw,frac_iw)
      if (dflag2) return
       r470 = frac_iw*sfc2(index_iw+1) + (1.-frac_iw)*sfc2(index_iw)
      endif

      if (mod_sfc.eq.11) then
      xthet5 = xthet
      if (xphi.gt.90.0) xthet5 = -1. * xthet
       if (xthet5.le.wav3(1)) then
          r470 = sfc3(1)
          return
       endif
       if (xthet5.ge.wav3(nc_wav(3))) then
          r470 = sfc3(nc_wav(3))
          return
       endif
       call search2(dflag2,xthet5,wav3,nc_wav(3),index_iw,frac_iw)
      if (dflag2) return
       r470 = frac_iw*sfc3(index_iw+1) + (1.-frac_iw)*sfc3(index_iw)
      endif

      if (mod_sfc.eq.12) then
      xthet5 = xthet
      if (xphi.gt.90.0) xthet5 = -1. * xthet
       if (xthet5.le.wav4(1)) then
          r470 = sfc4(1)
          return
       endif
       if (xthet5.ge.wav4(nc_wav(4))) then
          r470 = sfc4(nc_wav(4))
          return
       endif
       call search2(dflag2,xthet5,wav4,nc_wav(4),index_iw,frac_iw)
      if (dflag2) return
       r470 = frac_iw*sfc4(index_iw+1) + (1.-frac_iw)*sfc4(index_iw)
      endif

      if (mod_sfc.eq.14) then
      xthet5 = xthet
      if (xphi.gt.90.0) xthet5 = -1. * xthet
       if (xthet5.le.wav6(1)) then
          r470 = sfc6(1)
          return
       endif
       if (xthet5.ge.wav6(nc_wav(6))) then
          r470 = sfc6(nc_wav(6))
          return
       endif
       call search2(dflag2,xthet5,wav6,nc_wav(6),index_iw,frac_iw)
      if (dflag2) return
       r470 = frac_iw*sfc6(index_iw+1) + (1.-frac_iw)*sfc6(index_iw)
      endif

      if (mod_sfc.eq.15) then
      xthet5 = xthet
      if (xphi.gt.90.0) xthet5 = -1. * xthet
       if (xthet.le.wav5(1)) then
          r470 = sfc5(1)
          return
       endif
       if (xthet.ge.wav5(nc_wav(5))) then
          r470 = sfc5(nc_wav(5))
          return
       endif
       call search2(dflag2,xthet,wav5,nc_wav(5),index_iw,frac_iw)
      if (dflag2) return
       r470 = frac_iw*sfc5(index_iw+1) + (1.-frac_iw)*sfc5(index_iw)
      endif

      return
      end
