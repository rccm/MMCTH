      subroutine mod06CT_initialize_output( )
      
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Initialize output product arrays for MOD06CT 5x5 sampling.
c
c!Input Parameters:
c   
c    None.
c
c!Output Parameters:
c
c    None.
c
c    The following arrays in COMMON /MOD06_DATA/ are initialized:
c    height_tropopause           tropopause height
c    cloud_h_method              index of method used to generate cloud
c                                  top height
c    cloudtop_pre_ir             cloud height from "window" channel only
c    spec_cloud_forcing          cloud - clear radiance difference
c    cloudtop_pres               cloud top pressure (hPa)
c    cloudtop_height             cloud top height (meters)
c    cloudtop_temp               cloud top temperature (K)
c    cloudtop_eff_emi            effective cloud emissivity
c    cloudtop_pres_from_ratios   cloud top pressure from all ratios(hPa)
c    cloud_fraction              cloud fraction from cloud mask
c
c!Revision History:
c $Id: mod06CT_initialize_output.f,v 1.4 1999/04/16 22:55:22 kis Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      save
          
      include 'mod06uw_data.inc'
      
c ... arguments - none

c ... local variables

      integer i, j, k

c ... Initialize all the output arrays

      do k = 1, max_solutns
        do j = 1, max_samp_line
          do i = 1, max_samp_pixel
            spec_cloud_forcing( i, j, k ) = bad_value
            cloudtop_pres_from_ratios( i, j, k ) = bad_value
          end do
        end do
      end do

      do j = 1, max_samp_line
        do i = 1, max_samp_pixel
          height_tropopause( i, j ) = bad_value
          cloud_h_method( i, j ) = bad_value
          cloudtop_pre_ir( i, j ) = bad_value
          cloudtop_pres( i, j ) = bad_value
          cloudtop_pres_nearnad( i, j ) = bad_value
          cloudtop_height( i, j ) = bad_value
          cloudtop_height_nearnad( i, j ) = bad_value
          cloudtop_temp( i, j ) = bad_value
          cloudtop_temp_nearnad( i, j ) = bad_value
          cloudtop_eff_emi( i, j ) = bad_value
          cloudtop_eff_emi_nearnad( i, j ) = bad_value
          cloud_fraction( i, j ) = bad_value
          cloud_fraction_nearnad( i, j ) = bad_value

        end do
      end do                      

      return
      end
