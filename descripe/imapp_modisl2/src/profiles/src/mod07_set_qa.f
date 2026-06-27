      SUBROUTINE MOD07_SET_QA( LINE, PIXEL, NGOOD, NMISSING, NCLOUDY,
     &  FLAG )

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Set the appropriate MOD07 QA bits for a box of pixels.
C
C !INPUT PARAMETERS:
C     LINE        Line number at corner of this box
C     PIXEL       Pixel number at corner of this box
C     NGOOD       Number of good (clear) pixels in this box
C     NMISSING    Number of missing (bad radiance) pixels in this box
C     NCLOUDY     Number of cloudy pixels in this box
C     FLAG        Retrieval success flag (1=successful retrieval)
C
C !OUTPUT PARAMETERS:
C     None
C
C     The following arrays in COMMON /MOD07_DATA/ are filled:
C     PRODUCT_QA        MOD07 product QA bit array
C     WATER_VAPOR_QA    IR water vapor product QA bit array
C
C !REVISION HISTORY:
C
c    January, 2009, Eva Borbas: bug fixed between line 111-114


C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C


C !END
C-----------------------------------------------------------------------
      
      IMPLICIT NONE
      
      INCLUDE 'mod07_data.inc'

c ... Arguments

      INTEGER line, pixel, ngood, nmissing, ncloudy, flag

c ... Local variables

      INTEGER iout, jout
      
c ... Pixel and line coordinates in output 5x5 sampled array
      
      iout = ( pixel / isamp ) + 1
      jout = ( line / isamp ) + 1

c ... If retrieval was good, set product run time QA flags      

      if ( flag .eq. 0 ) then
      
c ...   Retrieved Temperature Profile

        call set_qa_bit( product_qa( 1, iout, jout ), 0 )
        call set_qa_bit( product_qa( 1, iout, jout ), 1 )
        
c ...   Retrieved Moisture Profile

        call set_qa_bit( product_qa( 1, iout, jout ), 4 )
        call set_qa_bit( product_qa( 1, iout, jout ), 5 )

c ...   Total Ozone

        call set_qa_bit( product_qa( 1, iout, jout ), 8 )
        call set_qa_bit( product_qa( 1, iout, jout ), 9 )

c ...   Lifted Index

        call set_qa_bit( product_qa( 1, iout, jout ), 12 )
        call set_qa_bit( product_qa( 1, iout, jout ), 13 )

c ...   K Index

        call set_qa_bit( product_qa( 1, iout, jout ), 16 )
        call set_qa_bit( product_qa( 1, iout, jout ), 17 )

c ...   Total Totals Index

        call set_qa_bit( product_qa( 1, iout, jout ), 20 )
        call set_qa_bit( product_qa( 1, iout, jout ), 21 )

      endif

c ... Ancillary data QA

c ... Guess Moisture Profile - not used

      call set_qa_bit( product_qa( 1, iout, jout ), 56 )
      call set_qa_bit( product_qa( 1, iout, jout ), 57 )

c ... Guess Temperature Profile - not used

      call set_qa_bit( product_qa( 1, iout, jout ), 58 )
      call set_qa_bit( product_qa( 1, iout, jout ), 59 )

c ... Reynolds blended SST - not used
      call set_qa_bit( product_qa( 1, iout, jout ), 62 )
      call set_qa_bit( product_qa( 1, iout, jout ), 63 )

c ... TOMS ozone profile - not used

      call set_qa_bit( product_qa( 1, iout, jout ), 66 )
      call set_qa_bit( product_qa( 1, iout, jout ), 67 )

c ... Surface Temperature - not used
c     EvaB, January, 2009

cEB      call set_qa_bit( product_qa( 1, iout, jout ), 62 )
cEB      call set_qa_bit( product_qa( 1, iout, jout ), 63 )
      call set_qa_bit( product_qa( 1, iout, jout ), 60 )
      call set_qa_bit( product_qa( 1, iout, jout ), 61 )

c ... Surface Pressure - not used

      call set_qa_bit( product_qa( 1, iout, jout ), 64 )
      call set_qa_bit( product_qa( 1, iout, jout ), 65 )

c ... Set optional run time QA flags

      product_qa( 4, iout, jout ) = ncloudy
      product_qa( 5, iout, jout ) = ngood
      product_qa( 6, iout, jout ) = nmissing
      call set_qa_bit( product_qa( 1, iout, jout ), 51 )

c ... Set IR water vapor QA flags

      if ( flag .eq. 0 ) then
        call set_qa_bit( water_vapor_qa( 1, iout, jout ), 0 )
        call set_qa_bit( water_vapor_qa( 1, iout, jout ), 1 )
      endif
      water_vapor_qa( 2, iout, jout ) = ncloudy
      water_vapor_qa( 3, iout, jout ) = ngood
      water_vapor_qa( 4, iout, jout ) = nmissing
      call set_qa_bit( water_vapor_qa( 1, iout, jout ), 32 )
      
      END
