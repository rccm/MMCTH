      REAL FUNCTION DEWTEM( P, T, W )

C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     DETERMINES DEW POINT TEMPERATURE
C
C !INPUT PARAMETERS:
C     P = REFERENCE PRESSURE ( MB )
C     T = REFERENCE TEMPERATURE ( DEG K )
C     W = REFERENCE MIXING RATIO ( G/KG )
C
C !OUTPUT PARAMETERS:
C     DEWTEM = DEW POINT TEMPERATURE ( DEG K )
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !DESIGN NOTES:
C
C !END
C-----------------------------------------------------------------------

c ... Scalar Arguments

      REAL P,T,W

c ... Local Scalars

      REAL A,B,CENKEL,CON,EPS,WS,WW,Y

c ... External Functions

      REAL SATMIX
      EXTERNAL SATMIX

c ... Data statements

      DATA EPS,CON,CENKEL/0.621970585,6.1078,273.15/
      DATA A,B/17.2693882,35.86/

      DEWTEM = T
      WS = SATMIX(P,T)
      IF (W.LE.0.0) RETURN
      IF (WS.LE.0.0) RETURN
      IF (W.LT.WS) WW = W/1000.0
      IF (W.GE.WS) WW = WS/1000.0
      Y = ALOG(P/CON*WW/(EPS+WW))/A
      DEWTEM = (CENKEL-B*Y)/(1.0-Y)

      END
