      SUBROUTINE COPY_CLDMASK_META( CLDMSK_LRN, PCN4S, PCNND,
     &  PCNNT, PCNNG, PCNNI, PCNNL, PCNNW, PCNNS, PCNNV, PCNNR,
     &  PCNNC, MAX_SOL, MIN_SOL )

C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C    To extract Cloud Mask Product Specific Attributes from the
C    the cloud mask product HDF file.
C
C !INPUT PARAMETERS:
C    MCF_CMLRN     Variable containing the PCF logical
C                  reference number (LRN) to the cloud
C                  mask metadata configuration file (MCF)
C
C !OUTPUT PARAMETERS:
C    PCN4S         % of pixels in category 4 per granule
C    PCNND         % of day processed pixels per granule
C    PCNNT         % of night processed pixels per granule
C    PCNNG         % of sunglint found pixles per granule
C    PCNNI         % of snow/ice processed pixels per granule
C    PCNNL         % of land processed pixels per granule
C    PCNNW         % of water processed pixels per granule
C    PCNNS         % of shadow pixels found per granule
C    PCNNV         % of thin cirrus (vis) pixels per granule
C    PCNNR         % of thin cirrus (ir) pixels per granule
C    PCNNC         % of non-cloud obstruction pixels per granule
C    MAX_SOL       Maximum solar zenith angle for this granule
C    MIN_SOL       Minimum solar zenith angle for this granule
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'debug.inc'
      INCLUDE 'PGS_SMF.f'

c ... Arguments 

      INTEGER cldmsk_lrn
      REAL pcn4s, pcnnd, pcnnt, pcnng, pcnni, pcnnl,
     &  pcnnw, pcnns, pcnnv, pcnnr, pcnnc, max_sol, min_sol

c ... Local variables

      REAL meta_misg 
      PARAMETER (meta_misg = -99999.0)

      INTEGER i, rtn, version

      CHARACTER*36  attr, field_name(13)
      CHARACTER*70 buf_char(13)
      
c ... External functions

      INTEGER    pgs_met_getpcattr_s
      EXTERNAL   pgs_met_getpcattr_s

c ... Initialize holding arrays

      do i = 1, 13
         buf_char(i) = '  '
      end do

c ... Initialize output variables

      pcn4s   = meta_misg
      pcnnd   = meta_misg
      pcnnt   = meta_misg
      pcnng   = meta_misg
      pcnni   = meta_misg
      pcnnl   = meta_misg
      pcnnw   = meta_misg
      pcnns   = meta_misg
      pcnnv   = meta_misg
      pcnnr   = meta_misg
      pcnnc   = meta_misg
      max_sol = meta_misg
      min_sol = meta_misg

c ... Get CoreMetadata.0 attribute fields that have values retrieved 
c ... from mask product.  Note that the use of M-API parameters are 
c ... defined in "mapi.inc".

      Field_Name(1)  = 'PARAMETERVALUE.5'
      Field_Name(2)  = 'PARAMETERVALUE.8'
      Field_Name(3)  = 'PARAMETERVALUE.9'
      Field_Name(4)  = 'PARAMETERVALUE.10'
      Field_Name(5)  = 'PARAMETERVALUE.11'
      Field_Name(6)  = 'PARAMETERVALUE.12'
      Field_Name(7)  = 'PARAMETERVALUE.13'
      Field_Name(8)  = 'PARAMETERVALUE.14'
      Field_Name(9)  = 'PARAMETERVALUE.15'
      Field_Name(10) = 'PARAMETERVALUE.16'
      Field_Name(11) = 'PARAMETERVALUE.17'
      Field_Name(12) = 'PARAMETERVALUE.18'
      Field_Name(13) = 'PARAMETERVALUE.19'

      do i = 1, 13
      
        if ( i.ne.8 ) then
          attr = Field_Name(i)
          version = 1
          rtn = pgs_met_getpcattr_s( cldmsk_lrn, version,
     &    'CoreMetadata.0', attr, buf_char(i) )

          if ( rtn .ne. PGS_S_SUCCESS ) then
            call message( 'Copy_cldmask_meta',
     &      'Failed to extract PSA from Cloud Mask File ' //
     &      ' [OPERATOR ACTION: Contact SDST]', 0, 2 )
            buf_char(i) = '-9999900'
          endif
       endif

      end do

c ... Set return values

      read( buf_char(1), '(f8.2)' ) pcn4s
      read( buf_char(2), '(f8.2)' ) pcnnd
      read( buf_char(3), '(f8.2)' ) pcnnt
      read( buf_char(4), '(f8.2)' ) pcnng
      read( buf_char(5), '(f8.2)' ) pcnni
      read( buf_char(6), '(f8.2)' ) pcnnl
      read( buf_char(7), '(f8.2)' ) pcnnw
c      read( buf_char(8), '(f8.2)' ) pcnns
      read( buf_char(9), '(f8.2)' ) pcnnv
      read( buf_char(10),'(f8.2)' ) pcnnr
      read( buf_char(11),'(f8.2)' ) pcnnc
      read( buf_char(12),'(f8.2)' ) max_sol
      read( buf_char(13),'(f8.2)' ) min_sol

c ... Print debug information if required

      if ( debug .gt. 0 ) then

      write (h_output, '(/,72(''-''),/)' )
      write (h_output, '(''Copy_cldmask_meta debug info'',/)' )
      write (h_output,fmt='(1x,''Metadata Cloud Mask stats: '')')
      write (h_output,fmt='(1x,''% with confidence <  1% '',f9.2)')pcn4s
      write (h_output,fmt='(1x,''percent NCO found '',f9.2)')pcnnc
      write (h_output,fmt='(1x,''percent shadow found '',f9.2)')pcnns
      write (h_output,fmt='(1x,''percent day processed '',f9.2)')pcnnd
      write (h_output,fmt='(1x,''percent night processed '',f9.2)')pcnnt
      write (h_output,fmt='(1x,''percent sunglint found '',f9.2)')pcnng
      write (h_output,fmt='(1x,''percent snow/ice processed '',f9.2)')
     +                           pcnni
      write (h_output,fmt='(1x,''percent land processed '',f9.2)')pcnnl
      write (h_output,fmt='(1x,''percent water processed '',f9.2)')pcnnw
      write (h_output,fmt='(1x,''percent thin cirrus solar found '',
     +                    f9.2)')pcnnr
      write (h_output,fmt='(1x,''percent thin cirrus ir found '',f9.2)')
     +                           pcnnv
      write (h_output,fmt='(1x,''max solar zenith angle '',f9.2)')
     +                           max_sol
      write (h_output,fmt='(1x,''min solar zenith angle '',f9.2)')
     +                           min_sol

      endif

      END
