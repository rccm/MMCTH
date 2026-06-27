      SUBROUTINE TPRECW(P, W, U, NL, NS, NP, PS)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Compute total  precipitable water  and integrated layer water vapour 
C             given a top-down profile of water vapor mixing ratios.
C
C !INPUT PARAMETERS:
C     P     Pressure profile (millibars)
C     W     Water vapor mixing ratio profile (grams/kilogram)
C     NL    Number of levels (assumed to be 101 levels in this version)
C     NS    Starting level in integration
C     NP    Surface level  or end level in profile integration
C     PS    surface pressure or pressure of the end level
C
C !OUTPUT PARAMETERS:
C     U     Total precipitable water vapor (centimeters)
C
C !REVISION HISTORY:
C
C  Nov 2002 SWS changed: compute increment using average of two 
C              levels around ps instead of only the one given
C  Sept 2003 SWS changed: fixed signs after 89th line 

C   Nov 2009 Eva Borbas (CIMSS/SSEC)
C              tprecw_new2.f has been modified to tprecwv.f for 
C               calculation TPW (10hPa-sfc)
C                           TPW_LO (700hPa-sfc)
C                           TPW_Hi (10hPa-300hPa)
C 
C     index 21 corresponds to 9.512 hPa (~10 hPa)
C     index 45 corresponds to 103.02 hPa (~100 hPa)
C     index 64 corresponds to 300.00 hPa (~300 hPa)
C     index 76 corresponds to 496.63 hPa (~500 hPa)
C     index 86 corresponds to 706.565 hPa (~700 hPa)

C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Scalar arguments
      real u, ps
      integer nl, np, np2, ns

c ... Array arguments
      real p(nl), w(nl)

c ... Local scalars
      real dp, f, p1, p2, s, w1, w2
      integer i

c ... Data statements
      data f/1961.33/

c ... Integrate the profile
c EvaB Nov 2009: added ns variable

      w1 = w(ns)
      p1 = p(ns)
      s = 0.0

C Nov 2009, EvaB added the following if condition
      if (ns+1.le.np) then

      do i = ns+1, np
          w2 = w(i)
          p2 = p(i)
          dp = abs(p2 - p1)
          s = s + (w1 + w2) * (dp / f)
          w1 = w2
          p1 = p2
      end do

c * Add (subtract) increment of atmospheric tpw to reach surface.
c    DP will be negative if PS < P(NP)

      DP= PS-P(NP)

C SWS changed 20 nov02: compute increment using average of two 
C   levels around ps instead of only the one given
c      S = S + 2*W(NP)*DP/F

C SWS changed 12 September 2003: fixed signs here
C NOTE this requires p values be ascending (TOA to surface)
C Nov 2009, EvaB added DP.ne.0 if condition

      if (DP.ne.0) then

      if(p(np).gt.ps) then
          np2 = np - 1
      else if(p(np).lt.ps .and. np+1.le.nl) then
          np2 = np + 1
      else
          np2 = np
      
      end if  

      s = s + (w(np2)+w(np))*dp/f
      
      endif
       
      endif
      
      U=S 

      END

