      REAL FUNCTION SVPWAT( TEMP )

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C Compute saturation vapor pressure over water for a given temperature.
C
C !INPUT PARAMETERS:
C REAL            TEMP        Temperature (K)
C
C !OUTPUT PARAMETERS:
C REAL            SVPWAT      Saturation vapor pressure over water (mb)
C
C !REVISION HISTORY:
C
c!Team-Unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

      implicit none

c ... arguments

      real temp
      
c ... local variables

      double precision b,s,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9

c ... data statements

      data a0/.999996876d0/,a1/-.9082695004d-2/,a2/.7873616869d-4/,
     *  a3/-.6111795727d-6/,a4/.4388418740d-8/,a5/-.2988388486d-10/,
     *  a6/.2187442495d-12/,a7/-.1789232111d-14/,a8/.1111201803d-16/,
     *  a9/-.3099457145d-19/,b/.61078d+1/
     
      t = dble( temp - 273.16 )
      s = a0 + 
     &  t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*a9))))))))
      s = b / s**8
      svpwat = real( s )

      END
