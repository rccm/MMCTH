       subroutine assign_eco_emis(landsea, emis, rho, freqemis) 

      implicit none

c**************************************************************************	
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c       Fortran code of assign_eco_emis.m (Matlab code from Suzanne Seemann) 
c       November 24, 2003  Eva Borbas 
c
c!Input Parameters:
c        input:  landsea (flag: 1 for land, 0 for ocean)
c
c!Output Parameters:
c        output: emis
c                rho
c                freqemis
c
c!Revision History:
c       25 nov 03: sws changed from wavelength to freq (wavenumber) and 
c                  reversed the order of emissivities to go in order of 
c                  ascending wavenumber
c 
c       15 jan 04: sws added input parameter landsea (1 land, 0 ocean 
c                  from modis mod03 flag).  
c
c       21 jan 04: sws moved igbp class 18 (tundra) from ice/snow to land
c
c       11 feb 04: sws simplified for land/ocean emissivity only - for mod06
c       16 feb 04: sws removed calls to anoise - no random noise added to emis
c       05 mar 04: raf removed lat, lon from calling arguments
c       24 May 04: G. Fireman Renamed from "assign_eco_emis_landSeaOnly_noNoise.f"
c                  to "assign_eco_emis.f" to avoid function name length errors.
c
c!Team-unique Header:
c
c!End
c
c***************************************************************************

        integer landsea

        integer nb_wavelen,i
        parameter(nb_wavelen=7)

        real emis(nb_wavelen),rho(nb_wavelen)

        real freq(nb_wavelen), freqemis(nb_wavelen)
        real seawater_emis(nb_wavelen)
        real land_emis(nb_wavelen)
        real pi

        data pi/3.14159265358979/

        data freq/700.0, 926.0, 1075.0, 1205.0, 1316.0, 
     &               2000.0, 2326.0/

        data seawater_emis/0.97, 0.99297, 0.98648, 0.9807, 0.97996, 
     &                     0.97999, 0.97646/

c        data land_emis/0.98, 0.97, 0.95, 0.95, 0.98, 0.95, 0.9/
c  new values for the first two land_emis - based on mod11 band 31 & 32  
c   averages of igbp 1-14 for jan, april, july, and october 2003
        data land_emis/0.9766, 0.9626, 0.95, 0.95, 0.98, 0.95, 0.9/

        do i=1,nb_wavelen
                emis(i)=0.
                rho(i)=0.
                freqemis(i) = freq(i)
        enddo

        do i=1,nb_wavelen

c     OCEAN: 
        if(landsea.eq.0) then 
           emis(i) = seawater_emis(i)

c     LAND:
        elseif(landsea.eq.1) then
              emis(i) = land_emis(i)
        endif

c     check for high emissivity	

           if(emis(i).gt.0.995) then
              emis(i) = 0.995 - (emis(i)-0.995)
           endif

c     compute rho profiles (1-emis)/pi

           rho(i) = (1 - emis(i))/pi

        enddo

        return
        end
