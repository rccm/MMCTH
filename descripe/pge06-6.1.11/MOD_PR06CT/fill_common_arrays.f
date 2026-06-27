      subroutine fill_common_arrays( line, pixel )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Compute MOD07 retrieval products.
c
c!Input Parameters:
c    LINE                  Line number within a 1km scan (1-10)
c    PIXEL                 Pixel number within a 1km scan (1-1500)
c
c    The following arrays in COMMON /MOD07_DATA/ are used:
c    RADIANCE1             Radiances for IR bands
c    TEM1                  Surface temperature
c    MASK1                 Cloudmask
c    eco                   1km ecosystem type
c
c!Output Parameters:
c    FLAG                  Success flag (0=Success, -1=Failure)
c
c    The following arrays in COMMON /MOD07_DATA/ are filled:
c    BRIGHTNESS_TEMP       Brightness temperature
c    SFC_TEMP              Surface temperature
c    SFC_PRES              Surface pressure
c
c!Revision History:
c$Id: fill_common_arrays.f,v 1.6 1999/06/11 22:34:27 kis Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
      
      implicit none
      save
      
      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'
      
c ... arguments

      integer line, pixel
      
c ... local variables

      integer i, j, k, n, iout, jout, units
      real rad, rsum, ravg

c ... local arrays
      integer o_bands(n_output),m_bands(n_output)
      real tbrt(n_output)

c ... external functions

      real modis_bright
      external modis_bright

c ... data statments
      data o_bands / 8, 10, 11, 12, 13, 14, 15 /
      data m_bands / 29, 31, 32, 33, 34, 35, 36 /

c ... Initialize brightness temperature array for this box

      do i = 1, n_output
        tbrt( i ) = bad_value
      end do

c ... Loop over IR bands

      do k = 1 , n_output

c ...   Sum the good radiances within this box

        rsum = 0.0
        n = 0
        do j = 1, isamp
          do i = 1, isamp
            rad = radiance1( pixel + i - 1, line + j - 1, o_bands(k))
            if ( rad .gt. 0.0) then
              rsum = rsum + rad
              n = n + 1
            endif
          end do
        end do

c ...   If sufficient good pixels were found,
c ...   compute average radiance and convert to brightness temp
        
        if ( n .ge. min_ngood .and. n .gt. 0 ) then
          ravg = rsum / real( n )
          units = 1
          tbrt(k) = modis_bright( ravg, m_bands(k), units )
        endif

      end do

c-----------------------------------------------------------------------
c ... Store results in product arrays
c-----------------------------------------------------------------------

c ... Pixel and line coordinates in input 1km array (center of box)

      i = pixel + isamp / 2
      j = line + isamp / 2

c ... Pixel and line coordinates in output 5x5 sampled array
      
      iout = ( pixel / isamp ) + 1
      jout = ( line / isamp ) + 1

c ... Copy the brightness temperatures
      
      do k = 1 , n_output
        brightness_temp( iout, jout, k ) = tbrt( k )
      end do

c ... Parameters copied from input ancillary data
      
      sfc_temp( iout, jout ) = tem1( i, j )
      sfc_pres( iout, jout ) = pre1( i, j )
      day_night_flag( iout, jout ) = dn_flag( i, j )

c-----------------------------------------------------------------------
c ... Write debug information
c-----------------------------------------------------------------------

      if ( debug .gt. 0 ) then

        if ( line .eq. 6 .and. pixel .eq. 676 ) then

          write( h_output, '(/,a,/)' ) 'mod07_fill_common_arrays debug'

          write( h_output, '(''Input Data for Nadir Box'')' )
          write( h_output, '(''IR Bands: '',7i7)' )
     &      ( m_bands( i ), i = 1, n_output )
          write( h_output, '(''Avg Temp: '',7f7.2)' )
     &      ( tbrt( i  ), i = 1, n_output )

        endif

      endif

      return
      end
