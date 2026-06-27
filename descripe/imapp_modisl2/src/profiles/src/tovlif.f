      SUBROUTINE TOVLIF( NLEV, P, T, TD, T500 )

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Compute temperature of parcel lifted dry-adiabatically to 500 mb.
C
C !INPUT PARAMETERS:
C     NLEV        Number of input levels in profile
C     P           Pressure levels in profile (mb)
C     T           Temperature profile (K)
C     TD          Dewpoint temperature profile (K)
C
C !OUTPUT PARAMETERS:
C     T500        Temperature of parcel lifted to 500 mb (K)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Arguments

      INTEGER nlev
      REAL p(*), t(*), td(*), t500

c ... Local variables

      REAL dpl, dpp, dsfc, pbndy, pl, ppl, s, ss, tbr, tl, tsfc, wbpt,
     &  wbr, x
      INTEGER k, kk, nlvl, iflag

c ... External functions

      REAL powt, wlift5
      EXTERNAL powt, wlift5

c ... Data statements

      DATA pbndy/100.0/

      nlvl = nlev
      t500 = 999999.
      if (nlvl.lt.2) return

c ... Find data values at top of boundary layer.

      pl = p(1) - pbndy
      wbr = 0.0
      tbr = 0.0
      s = 0.0

      iflag = 0

      do kk = 2,nlvl

         if ( p(kk) .lt. pl .and. iflag .eq. 0 ) then
            iflag = 1
            k = kk
         end if

         if ( iflag .eq. 0 ) then
            x = p(kk-1) - p(kk)
            wbr = wbr + (td(kk-1)+td(kk))*x/2.0
            tbr = tbr + (t(kk-1)+t(kk))*x/2.0
            s = s + x
         end if

      end do

      if ( iflag .eq. 0 ) return

      ss = p(k-1) - pl
      dpp = p(k-1) - p(k)
      tsfc = t(k-1) - ss* (t(k-1)-t(k))/dpp
      dsfc = td(k-1) - ss* (td(k-1)-td(k))/dpp
      wbr = wbr + (td(k-1)+dsfc)*ss/2.0
      tbr = tbr + (t(k-1)+tsfc)*ss/2.0
      ss = ss + s
      dpl = wbr/ss - 273.16
      tl = tbr/ss - 273.16
      ppl = p(1) - pbndy/2.0

c ... Compute the wet bulb potential temperature.

      wbpt = powt( tl, ppl, dpl )

c ... Compute temperature at the point where the moist adiabat
c ... crosses 500 mb. ideally, this is the parcel temp at 500 mb.

      t500 = wlift5( wbpt ) + 273.16

      END

      REAL FUNCTION POWT( T, P, TD )

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Compute wet-bulb potential temperature
C
C !INPUT PARAMETERS:
C     T       Temperature (celsius)
C     P       Pressure (millibars)
C     TD      Dew point (celsius)
C
C !OUTPUT PARAMETERS:
C     POWT    Wet-bulb potential temperature (celsius)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------
      
      IMPLICIT NONE

c ... Arguments

      REAL t, p, td

c ... Local variables

      REAL akap, atc, pt, tc

c ... External functions

      REAL tcon, wobf
      EXTERNAL tcon, wobf

c ... Data statements

c ... (Gas constant for dry air) /
c ... (Specific heat at constant pressure for dry air)

      DATA akap / 0.28541 /

c ... Freezing point of water (kelvin)

      DATA atc / 273.16 /

c ... Compute the potential temperature (celsius)

      pt = (t+atc)* (1000./p)**akap - atc

c ... Compute the lifting condensation level (lcl).

      tc = tcon(t,td)

c ... For the origin of the following approximation, see
c ... the documentation for the wobus function 'wobf'.

      powt = pt - wobf(pt) + wobf(tc)

      END

      REAL FUNCTION TCON( T, D )
      
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Compute temperature at lifting condensation level.
C
C !INPUT PARAMETERS:
C     T       Temperature (celsius)
C     D       Dew point temperature (celsius)
C
C !OUTPUT PARAMETERS:
C     TCON    Temperature (celsius) at lifting condensation level
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------
      
      IMPLICIT NONE

c ... Arguments

      REAL d, t

c ... Local variables

      REAL dlt, s

c ... Compute the dew point depression

      s = t - d

c ... The approximation below, a third order polynomial in s and t,
c ... is due to herman wobus. The source of data for fitting the
c ... polynomial is unknown.

      dlt = s*(1.2185+1.278e-03*t+s*(-2.19e-03+1.173e-05*s-5.2e-06*t))
      tcon = t - dlt

      END

      REAL FUNCTION WLIFT5( THETW )
      
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Compute temperature at 500mb for input theta-w lifted along moist 
C     adiabat.
C
C !INPUT PARAMETERS:
C     THETW     Theta-w
C
C !OUTPUT PARAMETERS:
C     WLIFT5    Temperature at 500mb (Celsius)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------
      
      IMPLICIT NONE

c ... Arguments

      REAL thetw

      wlift5 = -41.536 + thetw* (1.36083317+
     +         thetw* (1.91780552e-2+thetw* (1.3333332e-4+
     +         thetw* (-1.66611135e-5+thetw* (-2.46666673e-7+
     +         thetw* (8.805555540e-9))))))

      END

      REAL FUNCTION WOBF( T )

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C    Compute difference between wet-bulb potential temp for saturated air
C    and wet-bulb potential temp for completely dry air at given temp.
C
C !INPUT PARAMETERS:
C    T       Temperature (celsius)
C
C !OUTPUT PARAMETERS:
C    WOBF    Difference between wet-bulb potential temp for saturated
C            air and wet-bulb potential temp for completely dry air
C            at given temp (celsius)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !DESIGN NOTES:
C        Let wbpts = wet-bulb potential temperature for saturated
C   air at temperature t (celsius). Let wbptd = wet-bulb potential
C   temperature for completely dry air at the same temperature t.
C   The Wobus function wobf (in degrees celsius) is defined by
C                      wobf(t) = wbpts-wbptd.
C   Although wbpts and wbptd are functions of both pressure and
C   temperature, their difference is a function of temperature only.
C   To understand why, consider a parcel of dry air at temperature t and
C   pressure p. The thermodynamic state of the parcel is represented by
C   a point on a pseudoadiabatic chart. The wet-bulb potential tempera-
C   ture curve (moist adiabat) passing through this point is wbpts.
C   Now t is the equivalent temperature for another parcel saturated at
C   some lower temperature tw, but at the same pressure p.  To find tw,
C   ascend along the dry adiabat through (t,p).  At a great height, the
C   dry adiabat and some moist adiabat will nearly coincide.  Descend
C   descend along this moist adiabat back to p.  The parcel temperature
C   is now tw. The wet-bulb potential temperature curve (moist adiabat)
C   through (tw,p) is wbptd.  The difference (wbpts-wbptd) is propor-
C   tional to the heat imparted to a parcel saturated at temperature tw
C   if all its water vapor were condensed.  Since the amount of water
C   vapor a parcel can hold depends on temperature alone, (wbptd-wbpts)
C   must depend on temperature alone.
C
C        The Wobus function is useful for evaluating several thermo-
C   dynamic quantities.  By definition:
C                   wobf(t) = wbpts-wbptd.               (1)
C   If t is at 1000 mb, then t is a potential temperature pt and
C   wbpts = pt. Thus
C                   wobf(pt) = pt-wbptd.                 (2)
C   If t is at the condensation level, then t is the condensation
C   temperature tc and wbpts is the wet-bulb potential temperature
C   wbpt. Thus
C                   wobf(tc) = wbpt-wbptd.               (3)
C   If wbptd is eliminated from (2) and (3), there results
C                   wbpt = pt-wobf(pt)+wobf(tc).
C   If wbptd is eliminated from (1) and (2), there results
C                   wbpts = pt-wobf(pt)+wobf(t).
C
C        If t is an equivalent potential temperature ept (implying
C   that the air at 1000 mb is completely dry), then wbpts = ept
C   and wbptd = wbpt. Thus
C                   wobf(ept) = ept-wbpt.
C   This form is the basis for a polynomial approximation to wobf.
C   In Table 78 on pp.319-322 of the Smithsonian Meteorological
C   Tables by Roland List (6th revised edition), one finds wet-bulb
C   potential temperatures and the corresponding equivalent potential
C   temperatures listed together.  Herman Wobus, a mathematician for-
C   merly at the Navy Weather Research Facility, Norfolk, Virginia,
C   and now retired, computed the coefficients for the polynomial
C   approximation from that tabulated data.
C
C                                    Notes by T.W. Schlatter
C                                    NOAA/ERL/PROFS Program Office
C                                    August 1981
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Arguments

      REAL t

c ... Local variables

      REAL x, pol

      x = t - 20.0
      
      if ( x.gt.0.0 ) then
         pol=1.+x*(3.6182989e-3+x*(-1.3603273e-5+x*(4.9618922e-7)))
         wobf=29.930/pol**4+0.96*x-14.8
      else
         pol=1.+x*(-8.8416605e-3+x*(1.4714143e-4+x*(-9.6719890e-7)))
         wobf=15.130/pol**4
      endif

      END
