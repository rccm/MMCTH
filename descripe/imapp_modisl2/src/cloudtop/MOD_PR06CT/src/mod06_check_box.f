      subroutine mod06_check_box( line, pixel, nmissing, ncloud,
     &      ncloud_co2, ngood_msk, ngood_co2, ngood_IRphase, ngood_dayIRphase )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Check a 5x5 sampling box to see if there is enough good data
c    to perform a MOD06 retrieval.
c
c!Input Parameters:
c    LINE        Line number within a 1km scan (1-10)
c    PIXEL       Pixel number within a 1km scan (1-1500)
c
c!Output Parameters:
c    NGOOD       Number of good pixels found within box
c                (where BOX_FLAG=1)
c    NMISSING    Number of missing pixels within box (bad radiances)
c    NCLOUD      Number of cloudy pixels within box
c                   (from valid cloud mask pixels)
c    NCLOUD_CO2  Number of cloudy pixels within box
c                   (must also be valid pixel for retrieving CO2 cloud top pressure)
c    NGOOD_MSK   Number of good pixels as defined by bit 1 of cloud mask
c    NGOOD_CO2   Number of good pixels as defined by bit 1 of cloud mask
c                   for CO2-slicing algorithm
c    NGOOD_IRPHASE Number of good pixels as defined by bit 1 of cloud mask
c                   for IR Phasse algorithm
c    NGOOD_dayIRPHASE Number of good pixels as defined by bit 1 of cloud mask
c                   for daytime IR Phasse algorithm
c    The following array in COMMON /MOD07_DATA/ is filled:
c    BOX_FLAG    Valid data flag array for box
c                (0 = Bad, 1 = Good)
c                Good data implies that:
c                - Radiances are > 0.0 in all MODIS IR bands,
c                - Geolocation data are valid,
c                - Ancillary data are valid,
c                - Cloud mask was determined.
c    BOX_CLOUD   Clear/cloud flag array for box
c                (-1 = undefined, 0 = Clear, 1 = Cloud)
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none

      
      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'
      
c ... arguments

      integer line, pixel, nmissing, ncloud, ngood_msk, ncloud_co2
      integer ngood_IRphase, ngood_dayIRphase, ngood_co2
      integer nday, nglint

c ... local variables

      integer i, j, k, ibox, jbox
      byte mask_byte, mask_test, mask_value

c whuang 05/09/01: Added save statement for common / mod06_data / and
c                  / mod06_debug /
      save / mod06_data /, / mod06_debug /

c ... Set default return values

      nmissing    = 0
      ncloud      = 0
      ncloud_co2  = 0
      ngood_msk   = 0
      !** some additions by Baum 8-7-2000
      ngood_IRphase    = 0
      ngood_dayIRphase = 0
      ngood_co2        = 0
      
      nday   = 0  ! # pixels in box in daylight; based on dn_flag
      nglint = 0  !# pixels in box identified as containing sunglint
      
      do jbox = 1, isamp
        do ibox = 1, isamp
          box_flag( ibox, jbox )       = 1
          box_cloud( ibox, jbox )      = -1
	  box_flag_IR( ibox, jbox )    = 1
	  box_flag_co2( ibox, jbox )   = 1
	  box_flag_dayIR( ibox, jbox ) = 1
        end do
      end do

c ... Start loop over all pixels in this box

      jbox = 1
      do j = line, line + isamp - 1
        ibox = 1
        do i = pixel, pixel + isamp - 1

c ...     Check that cloud mask was determined for this pixel, and if so
c ...     check cloud mask
c ...     (99% and 95% confidence clear pixels are classified CLEAR)
c ...     (Low confidence clear and cloudy pixels are classified CLOUD)

          mask_byte = mask1( 1, i, j )
          mask_test = 1
          mask_value = iand( mask_byte, mask_test )
          if ( mask_value .eq. 1 ) then
            ngood_msk = ngood_msk + 1
            mask_test = 6
            mask_value = iand( mask_byte, mask_test )
            if ( mask_value .ge. 4 ) then
              box_cloud( ibox, jbox ) = 0
            else
              box_cloud( ibox, jbox ) = 1
            endif
          endif
	  
! now check sunglint flag in cloud mask
          mask_byte = mask1( 1, i, j )
          mask_test = 1
          mask_value = iand( mask_byte, mask_test )
          if ( mask_value .eq. 1 ) then
            mask_test = 16
            mask_value = iand( mask_byte, mask_test )
            if ( mask_value .lt. 16 ) nglint = nglint + 1
          endif

! ...     Check for valid radiance data for IR-only phase algorithm

          do k = 1, irphase_band
            if ( radiance1( i, j, ir_aband(k) ) .le. 0.0 ) then
              box_flag_IR( ibox, jbox ) = 0
            endif
          end do
          if (box_flag_IR(ibox,jbox) .eq. 0) nmissing = nmissing + 1


! ...     Check for valid radiance data for CO2 slicing algorithm

          do k = 1, co2_band
            if ( radiance1( i, j, ir_co2bandptr(k) ) .le. 0.0 ) then
              box_flag_co2( ibox, jbox ) = 0
            endif
          end do


! ...     Check for valid radiance data for daytime cloud phase algorithm

          do k = 1, dayphase_band
            if ( radiance1( i, j, day_bandptr(k) ) .le. 0.0 ) then
              box_flag_dayIR( ibox, jbox ) = 0
            endif
          end do
 
c ...     Check ancillary data

          if ( tem1( i, j ) .le. 0.0 ) box_flag( ibox, jbox ) = 0
c          if ( wat1( i, j ) .le. 0.0 ) box_flag( ibox, jbox ) = 0
          if ( pre1( i, j ) .le. 0.0 ) box_flag( ibox, jbox ) = 0
          do k = 1, 101
            if ( tprof1( k, i, j ) .le. 0.0 ) box_flag( ibox, jbox ) = 0
            if ( wprof1( k, i, j ) .le. 0.0 ) box_flag( ibox, jbox ) = 0
          end do

c ...     End loop over all pixels in this box

          ibox = ibox + 1
        end do
        jbox = jbox + 1
      end do

c ... Count Good/Bad and Clear/Cloudy pixels in this box

        do jbox = 1, isamp
          do ibox = 1, isamp
	    if (box_flag_co2(ibox,jbox) .eq. 1 .and. 
     +             box_flag(ibox,jbox)  .eq. 1) then
                ngood_co2 = ngood_co2 + 1
		! this means that both radiance data and ancillary data are good
                if(box_cloud(ibox,jbox) .ne. -1 ) ncloud_co2 = ncloud_co2 + box_cloud( ibox, jbox )
	    endif
	    ngood_IRphase = ngood_IRphase + box_flag_IR( ibox, jbox )
	    ngood_dayIRphase = ngood_dayIRphase + 
     +           box_flag_dayIR(ibox, jbox)
            if (box_cloud(ibox,jbox) .ne. -1 ) then
              ncloud = ncloud + box_cloud( ibox, jbox )
            end if
	    ! is this a "daytime" or "nighttime" box of data?
	    if (dn_flag(ibox,jbox) .eq. 1) nday=nday+1
          end do
        end do
	
c ... Write debug info for nadir box

      if ( debug .gt. 0 ) then
c       if ( line .eq. 6 .and. pixel .eq. 676 ) then

          write( h_output, '(/,a,/)' ) 'mod06_check_box debug'

          write( h_output,
     &      '(''Found'',i4,'' Missing Pixels for Box'')') nmissing

          write( h_output,
     &      '(''Found'',2i4,'' Cloudy Pixels for Box'')') ncloud, ncloud_co2
     
          write( h_output,
     &      '(''Found'',2i4,'' Good Data Pixels for Box'')') ngood_co2,
     &          ngood_IRphase

c       endif
      endif

      end
