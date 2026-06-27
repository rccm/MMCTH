      REAL FUNCTION SATMIX( P, T )

C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     DETERMINES SATURATION MIXING RATIO
C
C !INPUT PARAMETERS:
C     P = REFERENCE PRESSURE ( MB )
C     T = REFERENCE TEMPERATURE ( DEG K )
C
C !OUTPUT PARAMETERS:
C     SATMIX = SATURATION MIXING RATIO ( G/KG )
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

      REAL P,T

c ... Local Scalars

      REAL EPS,ES,WS

c ... External Functions

      REAL ESWAT
      EXTERNAL ESWAT

c ... Data statements

      DATA EPS/0.621970585/

      ES = ESWAT(1,2,T)
      IF (P.LT.50.0) WS = 0.0
      IF (P.GE.50.0) WS = EPS*ES/(P-ES)
      SATMIX = 1000.0*WS

      END
