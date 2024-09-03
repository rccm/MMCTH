
module surfemis
contains

  subroutine assign_eco_emis(landsea, emis, rho, freqemis)

  implicit none

  ! c**************************************************************************
  ! c-----------------------------------------------------------------------
  ! c
  ! c!F77
  ! c
  ! c!Description:
  ! c       Fortran code of assign_eco_emis.m (Matlab code from Suzanne Seemann)
  ! c       November 24, 2003  Eva Borbas
  ! c
  ! c!Input Parameters:
  ! c        input:  landsea (flag: 1 for land, 0 for ocean)
  ! c
  ! c!Output Parameters:
  ! c        output: emis
  ! c                rho
  ! c                freqemis
  ! c
  ! c!Revision History:
  ! c       25 nov 03: sws changed from wavelength to freq (wavenumber) and
  ! c                  reversed the order of emissivities to go in order of
  ! c                  ascending wavenumber
  ! c
  ! c       15 jan 04: sws added input parameter landsea (1 land, 0 ocean
  ! c                  from modis mod03 flag).
  ! c
  ! c       21 jan 04: sws moved igbp class 18 (tundra) from ice/snow to land
  ! c
  ! c       11 feb 04: sws simplified for land/ocean emissivity only - for mod06
  ! c       16 feb 04: sws removed calls to anoise - no random noise added to emis
  ! c       05 mar 04: raf removed lat, lon from calling arguments
  ! c       24 May 04: G. Fireman Renamed from "assign_eco_emis_landSeaOnly_noNoise.f"
  ! c                  to "assign_eco_emis.f" to avoid function name length errors.
  ! c
  ! c!Team-unique Header:
  ! c
  ! c!End
  ! c
  ! c***************************************************************************

   integer landsea
!f2py intent(in) landsea
   integer nb_wavelen,i
   parameter(nb_wavelen=7)

   real emis(nb_wavelen),rho(nb_wavelen)
!f2py intent(out) emis, rho, freqemis
   real freq(nb_wavelen), freqemis(nb_wavelen)
   real seawater_emis(nb_wavelen)
   real land_emis(nb_wavelen)
   real pi

   data pi/3.14159265358979/

   data freq/700.0, 926.0, 1075.0, 1205.0, 1316.0, &
                 2000.0, 2326.0/

   data seawater_emis/0.97, 0.99297, 0.98648, 0.9807, 0.97996, &
                       0.97999, 0.97646/

  ! c        data land_emis/0.98, 0.97, 0.95, 0.95, 0.98, 0.95, 0.9/
  ! c  new values for the first two land_emis - based on mod11 band 31 & 32
  ! c   averages of igbp 1-14 for jan, april, july, and october 2003
   data land_emis/0.9766, 0.9626, 0.95, 0.95, 0.98, 0.95, 0.9/

   do i=1,nb_wavelen
           emis(i)=0.
           rho(i)=0.
           freqemis(i) = freq(i)
   enddo

   do i=1,nb_wavelen

  ! c     OCEAN:
   if(landsea.eq.0) then
      emis(i) = seawater_emis(i)

  ! c     LAND:
   elseif(landsea.eq.1) then
         emis(i) = land_emis(i)
   endif

  ! c     check for high emissivity

      if(emis(i).gt.0.995) then
         emis(i) = 0.995 - (emis(i)-0.995)
      endif

  ! c     compute rho profiles (1-emis)/pi

      rho(i) = (1 - emis(i))/pi

   enddo

   return
   end subroutine


   subroutine getiremis ( nemis, freqemis, freq, &
                       emisir, rhoir, emiss, rho)
   implicit none
  ! c***********************************************************************
  ! c!F77
  ! C
  ! c!Description:
  ! C ROUTINE: GETEMIS
  ! C
  ! C PURPOSE: Compute emissivities for IR and MW.
  ! C
  ! C  History of modifications
  ! C          This routine assumes that the emissivity is a linear function of
  ! C          frequency ( in microns ).
  ! c!Input Parameters:
  ! C   name       type      units       description
  ! C   ----       ----      -----       -----------
  ! C  nemis      integer                number of frequency limits
  ! C   freq       real    wave number   frequency of channel
  ! C freqemis     real                  frequency limit for emissivity calc
  ! C   emisir     real                  emissivities
  ! C   rhoir      real                  reflectivities
  ! C
  ! C!Output Parameters:
  ! C  output variables:
  ! C    name       type      units       description
  ! C    ----       ----      -----       -----------
  ! C    emiss      real                  emissivity
  ! C    rho        real                  reflectivity
  ! C
  ! C!Revision History:
  ! C
  ! C NOTE from Suzanne: It appears that freqemis and corresponding emisir
  ! C      must be in order of increasing wavenumber
  ! c
  ! c!Team-unique Header:
  ! c
  ! c!End
  ! c-----------------------------------------------------------------------

  ! C
  ! c***********************************************************************
  ! c  input variables
  ! c  ---------------
   integer nemis
   real    freqemis(nemis), emisir(nemis), rhoir(nemis), freq
! f2py intent(in) freqemis, emisir, rhoir, freq
  ! c  output variables
  ! c  ----------------
   real    emiss, rho
!f2py intent(out) emiss, rho

  ! c  local variables
  ! c  ---------------
   integer k

   real waveno, wv1, dwv

  ! c     **************************************************************
  ! c     ************************** Infrared **************************
  ! c     **************************************************************

   if ( freq .gt. 500.0 ) then

      if ( freq .le. freqemis(1) ) then
         emiss = emisir(1)
         rho   = rhoir(1)
      else if ( freq .ge. freqemis(nemis) ) then
         emiss = emisir(nemis)
         rho   = rhoir(nemis)
      else
         do k = 1, nemis-1
            if ( freq .lt. freqemis(k+1) ) go to 2100
         end do
         k = nemis - 1
  2100       waveno = 10000.0/freq
         wv1 = 10000.0/freqemis(k)
         dwv = (waveno - wv1) / (10000.0/freqemis(k+1) - wv1)

         emiss     = emisir(k) + dwv * ( emisir(k+1) - emisir(k) )
         rho       = rhoir(k)  + dwv * ( rhoir(k+1)  - rhoir(k) )

      end if

   else
     write(*,*), 'getirmis: Frequency exceeded the IR region'

   endif

   return
   end subroutine

end module surfemis
