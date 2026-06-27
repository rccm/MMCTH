      subroutine getiremis ( nemis, freqemis, freq,
     $                     emisir, rhoir, emiss, rho)
      implicit none
c***********************************************************************
c!F77
C
c!Description:
C ROUTINE: GETEMIS
C
C PURPOSE: Compute emissivities for IR and MW.
C
C  History of modifications
C          This routine assumes that the emissivity is a linear function of
C          frequency ( in microns ).
c!Input Parameters:
C   name       type      units       description
C   ----       ----      -----       -----------
C  nemis      integer                number of frequency limits
C   freq       real    wave number   frequency of channel
C freqemis     real                  frequency limit for emissivity calc
C   emisir     real                  emissivities
C   rhoir      real                  reflectivities
C
C!Output Parameters:
C  output variables:
C    name       type      units       description
C    ----       ----      -----       -----------
C    emiss      real                  emissivity
C    rho        real                  reflectivity
C
C!Revision History:
C
C NOTE from Suzanne: It appears that freqemis and corresponding emisir 
C      must be in order of increasing wavenumber
c
c!Team-unique Header:
c
c!End
c-----------------------------------------------------------------------

C
c***********************************************************************
c  input variables
c  ---------------
      integer nemis
      real    freqemis(nemis), emisir(nemis), rhoir(nemis), freq

c  output variables
c  ----------------
      real    emiss, rho

c  local variables
c  ---------------
      integer k

      real waveno, wv1, dwv

c     **************************************************************
c     ************************** Infrared **************************
c     **************************************************************

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
        call message( 'getiremis',
     &    'Frequency exceeded the IR region ' //
     &    ' [OPERATOR ACTION: Contact SDST]',
     &    0, 2 )
      endif

      return
      end

