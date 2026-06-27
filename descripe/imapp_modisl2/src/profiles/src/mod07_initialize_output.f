      SUBROUTINE MOD07_INITIALIZE_OUTPUT()
      
c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Initialize output product arrays for MOD07 at 5x5 sampling.
c
c !INPUT PARAMETERS:
c
c !OUTPUT PARAMETERS:
c     The following arrays in COMMON /MOD07_DATA/ are initialized:
c     BRIGHTNESS_TEMP       Brightness temperature
c     GUESS_TEMP_PROFILE    Guess temperature profile
c     GUESS_WVMR_PROFILE    Guess water vapour mixing ratio profile
c     RETR_TEMP_PROFILE     Retrieved temperature profile
c     RETR_DEWP_PROFILE     Retrieved dewpoint profile
c     RETR_WVMR_PROFILE     Retrieved water vapour mixing ratio profile
c     RETR_HITE_PROFILE     Retrieved geopotential height profile
c     RETR_OZONE_PROFILE     Retrieved ozone profile
c     SFC_TEMP              Surface temperature
c     SFC_PRES              Surface pressure
c     SFC_ELEV              Surface elevation
c     HEIGHT_TROPOPAUSE     Height of the tropopause
c     WATER_VAPOR           Water vapor for entire column, by integrating regression profiles
c     WATER_VAPOR_DIRECT    Water vapor for entire column, by direct regression
c     WATER_VAPOR_LOW       Water vapor for surface to 900 hPa
c     WATER_VAPOR_HIGH      Water vapor for 700 hPa to 300 hPa
c     TOTAL_OZONE           Total column ozone
c     TOTAL_TOTALS          Total totals index
c     LIFTED_INDEX          Lifted index
c     K_INDEX               K index
c     PRODUCT_QA            Product run time QA
c     WATER_VAPOR_QA        IR water vapor product QA
c
c !REVISION HISTORY:
c    May 2006: Eva Borbas included retr_ozone_proifle
c    26 February 2002: SWS included Water_Vapor_Direct
c     Jan 2009: dewpoint temperature and mixing ratio profiles are in the output now 

c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE
          
      INCLUDE 'mod07_data.inc'
      
c ... local variables

      INTEGER i, j, k

c ... Initialize all the output arrays

      do k = 1, nbands
        do j = 1, max_samp_line
          do i = 1, max_samp_pixel
            brightness_temp( i, j, k ) = bad_value
          end do
        end do
      end do

      do k = 1, outlevels
        do j = 1, max_samp_line
          do i = 1, max_samp_pixel
            guess_temp_profile( i, j, k ) = bad_value
            guess_wvmr_profile( i, j, k ) = bad_value
            retr_temp_profile( i, j, k ) = bad_value
            retr_dewp_profile( i, j, k ) = bad_value
            retr_wvmr_profile( i, j, k ) = bad_value
            retr_hite_profile( i, j, k ) = bad_value
c            retr_ozone_profile( i, j, k ) = bad_value
          end do
        end do
      end do
            
      do j = 1, max_samp_line
        do i = 1, max_samp_pixel
          sfc_temp( i, j ) = bad_value
          sfc_pres( i, j ) = bad_value
          sfc_elev( i, j ) = bad_value
          height_tropopause( i, j ) = bad_value
          water_vapor( i, j ) = bad_value
          water_vapor_direct( i, j ) = bad_value
          water_vapor_low( i, j ) = bad_value
          water_vapor_high( i, j ) = bad_value
          total_ozone( i, j ) = bad_value
          total_totals( i, j ) = bad_value
          lifted_index( i, j ) = bad_value
          k_index( i, j ) = bad_value
        end do
      end do                      

      do j = 1, max_samp_line
        do i = 1, max_samp_pixel
          do k = 1, nproduct_qa
            product_qa( k, i, j ) = 0
          end do
        end do
      end do

      do j = 1, max_samp_line
        do i = 1, max_samp_pixel
          do k = 1, nwater_vapor_qa
            water_vapor_qa( k, i, j ) = 0
          end do
        end do
      end do

      END
