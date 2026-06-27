      SUBROUTINE INDICES( P, T, W, PSFC, TSFC, WSFC, LSFC, 
     &  TLIFT, TT, KI )

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Compute atmospheric stability indices given temperature and water
C     vapor profiles.
C
C !INPUT PARAMETERS:
C     P           Pressure (mb) at the 40 standard levels
C     T           Temperature profile (K)
C     W           Water vapor mixing ratio profile (g/kg)
C     PSFC        Surface pressure (mb)
C     TSFC        Surface temperature (K)
C     WSFC        Surface water vapor mixing ratio (g/kg)
C     LSFC        Index of surface level in profile
C
C !OUTPUT PARAMETERS:
C     TLIFT       Lifted index (K)
C     TT          Total-totals index (K)
C     KI          K index (K)
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

      REAL p(*), t(*), w(*), psfc, tsfc, wsfc, tlift, tt, ki
      INTEGER lsfc
      
c ... Local variables

      REAL pst( 15 ), tst( 15 ), tsd( 15 )
      INTEGER i, nb

c ... External functions

      REAL dewtem
      EXTERNAL dewtem

c ... Compute lifted index
      
      nb = lsfc
      if( psfc .lt. p( nb ) ) nb = nb - 1
      do i = 2, 15
        pst( i ) = p( nb )
        tst( i ) = t( nb )
        tsd( i ) = dewtem( p( nb ), t( nb ), w( nb ) )
        nb = nb - 1
      end do
      pst( 1 ) = psfc
      tst( 1 ) = tsfc
      tsd( 1 ) = dewtem( psfc, tsfc, wsfc )
      i = 15
      call tovlif( i, pst, tst, tsd, tlift )
      tlift = t( 31 ) - tlift

c ... Compute total-totals (t850 + td850 - 2*t500)

      tt = t( 37 ) + dewtem( p( 37 ), t( 37 ), w( 37 ) ) -
     &  2.0 * t( 31 )
      
c ... Compute k-index (t850 + td850 - (t700 - td700) - t500)

      ki = t( 37 ) + dewtem( p( 37 ), t( 37 ), w( 37 ) ) -
     &  ( t( 35 ) - dewtem( p( 35 ), t( 35 ), w( 35 ) ) ) - t( 31 ) 

      END      
