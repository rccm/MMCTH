module set_arrays

	use science_parameters, only: grid_xsize, grid_ysize, model_levels
	use mod06_run_settings, only: set_number_of_bands
	use global_model_grids
	use ct_core_arrays
	use core_arrays
	use nonscience_parameters, only: fillvalue_real
	use specific_other

	implicit none
	
contains

	subroutine allocate_arrays(data_wid, data_ht)


		integer, intent(in):: data_wid, data_ht
		integer :: i, j, model_layers

		allocate(model_info(grid_xsize, grid_ysize))
   
		model_layers = model_levels
   
		do i=1, grid_xsize
			do j=1, grid_ysize
		
				allocate(model_info(i,j)%mixr_profile(model_layers))
				allocate(model_info(i,j)%temp_profile(model_layers))
				allocate(model_info(i,j)%height_profile(model_layers))
				allocate(model_info(i,j)%pressure_profile(model_layers))
				allocate(model_info(i,j)%o3_profile(model_layers))
		
			end do
		end do

#ifndef SEVIRI_INST
		allocate(clear_sky_btemp(data_wid, data_ht,set_number_of_bands))
		allocate(clear_sky_rad(data_wid, data_ht,set_number_of_bands))
#endif

		allocate(surface_temperature(data_wid, data_ht))
		allocate(sensor_zenith_angle(data_wid, data_ht))
		
		call allocate_extra(data_wid, data_ht)

		allocate(latitude(data_wid, data_ht))
		allocate(longitude(data_wid, data_ht))
		allocate(cloudmask(data_wid, data_ht))
		
		allocate(band_measurements(data_wid, data_ht, set_number_of_bands))
		allocate(surface_emissivity_land(data_wid, data_ht, set_number_of_bands))
		
		allocate(cloud_top_pressure(data_wid, data_ht))
		allocate(cloud_top_temperature(data_wid, data_ht))
		allocate(cloud_top_height(data_wid, data_ht))
		allocate(cloud_height_method(data_wid, data_ht))
		allocate(ir_cloudphase(data_wid, data_ht))
		
		
		surface_temperature = fillvalue_real
		cloud_top_pressure = fillvalue_real
		cloud_top_temperature = fillvalue_real
		cloud_height_method = 127
		cloud_top_height = fillvalue_real		
		
#ifndef SEVIRI_INST
		clear_sky_rad = fillvalue_real
		clear_sky_btemp = fillvalue_real
#endif	
	
	end subroutine allocate_arrays


	subroutine deallocate_arrays
	
		integer :: i, j
		
		do i=1, grid_xsize
			do j=1, grid_ysize
		
				deallocate(model_info(i,j)%mixr_profile)
				deallocate(model_info(i,j)%temp_profile)
				deallocate(model_info(i,j)%height_profile)
				deallocate(model_info(i,j)%pressure_profile)
				deallocate(model_info(i,j)%o3_profile)
		
			end do
		end do

		deallocate(model_info)
	
#ifndef SEVIRI_INST	
		deallocate(clear_sky_rad)
		deallocate(clear_sky_btemp)
#endif		
		
		deallocate(cloud_height_method)
		deallocate(cloud_top_temperature)
		deallocate(cloud_top_pressure)
		deallocate(cloud_top_height)
		
		deallocate(band_measurements)
		deallocate(sensor_zenith_angle)
		deallocate(latitude)
		deallocate(longitude)
		deallocate(cloudmask)
		deallocate(surface_temperature)
		deallocate(surface_emissivity_land)
		deallocate(ir_cloudphase)
		
		
		call deallocate_extra
	
	
	end subroutine deallocate_arrays




end module set_arrays
