      SUBROUTINE MOD07_CHECK_BOX( LINE, PIXEL, NGOOD, NMISSING,
     &  NCLOUDY, NLAND )

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Check a 5x5 sampling box to see if there is enough good data
C     to perform a MOD07 retrieval.
C
C !INPUT PARAMETERS:
C     LINE        Line number within a 1km scan (1-10)
C     PIXEL       Pixel number within a 1km scan (1-1500)
C
C !OUTPUT PARAMETERS:
C     NGOOD       Number of good pixels found within box (clear sky)
C     NMISSING    Number of missing pixels within box (bad radiances)
C     NCLOUDY     Number of cloudy pixels within box
C     NLAND       Number of good land pixels found within box
C
C     The following array in COMMON /MOD07_DATA/ is filled:
C     BOX_FLAG    Valid data flag array for box (0=bad, 1=good)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !End
C-----------------------------------------------------------------------

      IMPLICIT NONE
      
      INCLUDE 'debug.inc'
      INCLUDE 'mod07_data.inc'
      
c ... Arguments

      INTEGER line, pixel, ngood, nmissing, ncloudy, nland

c ... Local variables

      INTEGER i, j, k, ibox, jbox
      BYTE mask_byte, mask_test, mask_value

c ... Set valid data flags in this pixel box (0=bad, 1=good)

      do jbox = 1, isamp
        do ibox = 1, isamp
          box_flag( ibox, jbox ) = 1
        end do
      end do

c ... Check all pixels in this box for valid radiance data

      nmissing = 0

      jbox = 1

      do j = line, line + isamp - 1

        ibox = 1

        do i = pixel, pixel + isamp - 1

          do k = 1, nbands
            if ( radiance1( i, j, k ) .le. 0.0 )
     &        box_flag( ibox, jbox ) = 0
          end do
          if ( box_flag( ibox, jbox ) .eq. 0 ) nmissing = nmissing + 1

          ibox = ibox + 1

        end do

        jbox = jbox + 1

      end do

c ... Check cloud mask for all pixels in this box

      ncloudy = 0
      
      nland = 0
      
      jbox = 1
 
      do j = line, line + isamp - 1

        ibox = 1
 
        do i = pixel, pixel + isamp - 1

c ...     Get first byte of cloud mask for this pixel

          mask_byte = mask1( 1, i, j )
 
c ...     Check that cloud mask was determined for this pixel

          mask_test = 1
          mask_value = iand( mask_byte, mask_test )
          if ( mask_value .eq. 0 ) box_flag( ibox, jbox ) = 0

c ...     Check cloud mask (99% and 95% prob. clear pixels are good)

          mask_test = 6
          mask_value = iand( mask_byte, mask_test )

          if ( mask_value .lt. 4 ) then
            box_flag( ibox, jbox ) = 0
            ncloudy = ncloudy + 1
          endif
 
c ...     Check surface pressure (must be in range 100 to 1100 hPa)

          if ( pre1( i, j ) .lt. 100.0 .or. pre1( i, j ) .gt. 1100.0 )
     &      box_flag( ibox, jbox ) = 0
 
c ...     Check precipitable water (must be in range 0 to 110 kg/m^2)
          if ( pwat1( i, j ) .lt. 0.0 .or. pwat1( i, j ) .gt. 110.0 )
     &      box_flag( ibox, jbox ) = 0

c ...     If pixel is good and is classified as land,
c ...     increment land pixel counter

          if (box_flag( ibox, jbox ) .eq. 1) then
            if ( int(land1( i, j )) .eq. 1 .or.
     &           int(land1( i, j )) .eq. 2 .or.
     &           int(land1( i, j )) .eq. 4 ) then
              nland = nland + 1
            endif
          endif
          
          ibox = ibox + 1

        end do

        jbox = jbox + 1

      end do

c ... Count good pixels in this box

      ngood = 0
      do jbox = 1, isamp
        do ibox = 1, isamp
          if (box_flag( ibox, jbox ) .eq. 1) ngood = ngood + 1
        end do
      end do
     
c ... Write debug info for nadir box

      if ( debug .eq. 1 ) then

        if ( line .eq. 6 .and. pixel .eq. 676 ) then

          write( h_output, '(/,a,/)' ) 'MOD07_CHECK_BOX DEBUG INFO'

          write( h_output,
     &      '(''Found'',i4,'' Missing Pixels for Nadir Box'')') nmissing

          write( h_output,
     &      '(''Found'',i4,'' Cloudy Pixels for Nadir Box'')') ncloudy
     
          write( h_output,
     &      '(''Found'',i4,'' Land Pixels for Nadir Box'')') nland

          write( h_output,
     &      '(''Found'',i4,'' Valid Data Flags for Nadir Box; '',i4,' //
     &      ' '' are required for retrieval'')') ngood, min_ngood

          do jbox = 1, isamp
            write( h_output, '(10i4)' )
     &        ( box_flag( ibox, jbox ), ibox = 1, isamp )
          end do

        endif

      endif

      END
