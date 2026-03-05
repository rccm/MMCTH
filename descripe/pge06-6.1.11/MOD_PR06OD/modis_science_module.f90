module modis_science_module

implicit none

private

public:: scienceinterface



contains

subroutine scienceinterface(threshold_solar_zenith,      &
                            threshold_sensor_zenith,     &
                            threshold_relative_azimuth,  &
                            debug,                       &
                            status)

   use modis_cloudstructure, only:cloudphase
   use core_arrays
   use libraryarrays
   use libraryinterpolates
   use science_parameters
   use corescience_module
   use multi_layer_clouds
   use clear_sky_restoral
   use mod06_run_settings
   use nonscience_parameters
   use interpolate_libraries
   use global_model_grids, only:get_model_idx
   use ancillary_module, only: get_above_cloud_properties, given_P_get_T
	use get_retrieval_uncertainty
	use retrieval_irw
	use atmospheric_correction_module
	use retrieval_prep_logic
	use rtm_support
	use cloud_phase
	use specific_other
	use modis_numerical_module, only: linearinterpolation, bisectionsearch
	use names, only: MY_UNIT_LUN
	
   implicit none

   logical, intent(in) :: debug
   real,    intent(in) :: threshold_solar_zenith,      &
                          threshold_sensor_zenith,     &
                          threshold_relative_azimuth

   integer, intent(out) :: status


   integer             :: xdimension, ydimension,i,j, jj,count_interpolations, &
                          retrievalcount
   integer             :: retrieval_failcount, library_failcount, cloud_failcount, atmoscorr_failcount
   integer             :: cloudstatus, retrievalstatus, librarystatus, atmoscorrstatus,cloudiness_degree_250m
   real                :: diff_solar_zenith,           &
                          diff_sensor_zenith,          &
                          diff_relative_azimuth
  
 	real :: scattering_angle, dscat, diff_scat_angle, diff_scat_angle_ss
	real :: dsol, dsen, drel, diff_wind_speed
	real :: cur_wind_speed
	real :: cloud_top_temperature_water, cloud_top_temperature_ice, ctt_1km
 
	integer :: istart, iend, jstart, jend, myrad
  
   logical :: ldebug, sunglint_dust_test, lowvariability_confidence_test, put_back_cloud
   logical :: ir_cloudphase_1km_watercloud, ir_cloudphase_1km_icecloud
	integer :: model_i, model_j
	integer :: uncertain_start

   integer :: na_band_used, R1R2wavelengthIdx(2), absorbingband_index, uncertainty_nonabsorbing_1621
   real:: corr_meas(set_number_of_bands), temp_meas(set_number_of_bands), alb_meas(set_albedo_bands)
   character(10) :: phase
   real    :: cloud_reflectance(2),delta_reflectance(2), albedo_holder(2) , &
              delta_transmittance(2), tauRegimeThreshold(2), unc_reflectance(2)

   integer :: start_time, end_time, cmax, crate

	integer :: cnt_sza, cnt_vza, cnt_raz, cnt_cm_switch, cnt_wspeed, cnt_scat, cnt_wspeed_only
	logical :: wind_speed_only, interp_MS, interp_SS
	real :: irw_pressure
	real :: Tc_liquid, Tc_ice, Pc_liquid, Pc_ice, Pw_liquid, Pw_ice
	real :: rad_clr(2), bt_clr(2)
	integer :: ice_near(4), liq_near(4), nearest_used(4)
	real :: RSS_ice(4), RSS_liq(4), RSS_final(4)
	
	real :: unc_tau_real, unc_tau_1621_real, unc_tau16_real, unc_tau37_real
	real :: unc_re21_real, unc_re16_real, unc_re37_real, unc_re1621_real
	real :: unc_lwp21_real, unc_lwp16_real, unc_lwp37_real, unc_lwp1621_real
	
	real :: emission_pw(20), emission_Tc(20), sigma_R37_pw(20)

	real :: alt_ray_liq, alt_ray_ice, temp_pres
	real :: aod550, irw_dummy
	
	integer :: re_idx_low, re_idx_hi	

    logical :: vis1km_test		! KGM 3-4-13
	type(cloudphase) :: ir_cloudphase

  real MODIS_PLANCK
  external MODIS_PLANCK

   status = 0
 
	xdimension = size(optical_thickness_final, 1)
	ydimension = size(optical_thickness_final, 2)
 

   lastinterp_solar_zenith     = fillvalue_real
   lastinterp_sensor_zenith    = fillvalue_real
   lastinterp_relative_azimuth = fillvalue_real
   lastinterp_scat_angle = fillvalue_real
   lastinterp_wind_speed = fillvalue_real
   lastinterp_scat_angle_ss = fillvalue_real
   count_interpolations   = 0
   retrievalcount = 0
   retrieval_failcount = 0
   library_failcount = 0
   cloud_failcount  = 0
   atmoscorr_failcount = 0
 
 
	cnt_sza = 0
	cnt_vza = 0
	cnt_raz = 0
	cnt_cm_switch = 0
	cnt_wspeed = 0
	cnt_scat = 0
	cnt_wspeed_only = 0
 
   call init_science_arrays

	call init_rtm_vars
	call init_half_radii

	COX_MUNK = .false.
	last_COX_MUNK = .false. 
	wind_speed_only = .false. 
	interp_MS = .false. 


	color_first_time = .true.
	last_2way_angle = fillvalue_real
	

	if (min_solar_zenith > 0.98 .or. min_sensor_zenith > 0.98) then 
		drel = 30.
	else if (min_solar_zenith > 0.85 .or. min_sensor_zenith > 0.85) then 
		drel = 10.
	else
		drel = threshold_relative_azimuth
	endif
	
   call init_retrieval(library_taus)

   do i = 1, xdimension


! do not process lines 1 and 1354 of the data because there are issues with cloud mask, particularly for the last line
! of the data
		if (iterationX == 1 .and. i==1 .or. iterationX == number_of_iterationsX .and. i==xdimension ) then 
			cloudsummary(i,:)%cloudobserved = .false.
			do j=1, ydimension
				call assign_retrieval_error(i,j)
				retrieval_failcount = retrieval_failcount + 1
			end do
			cycle
		endif

	 do j=1, ydimension

	
		pixX = i
		pixY = j

		cloudsummary(i,j)%cloudmask_determined = .false.
		cloudsummary(i,j)%cloudobserved = .false.
		cloudsummary(i,j)%watercloud = .false.
		cloudsummary(i,j)%icecloud = .false.
 		cloudsummary(i,j)%unknowncloud = .false.
		CSR_flag_array(i,j) = 0

		liq_near = 0
		ice_near = 0

		emission_uncertainty_pw_ice = fillvalue_real
		emission_uncertainty_pw_liq = fillvalue_real
		emission_uncertainty_Tc_ice = fillvalue_real
		emission_uncertainty_Tc_liq = fillvalue_real
		sigma_R37_PW_ice = fillvalue_real
		sigma_R37_PW_liq = fillvalue_real

! can't retrieve, sun too low or no cloud top properties
       if (solar_zenith_angle(i,j) > solar_zenith_threshold .or. solar_zenith_angle(i,j) < 0. ) then

			retrievalstatus = 1
			call assign_retrieval_error(i,j)
			retrieval_failcount = retrieval_failcount + 1
			cycle 
        endif



		scattering_angle = ScatAngle(solar_zenith_angle(i,j),sensor_zenith_angle(i,j),relative_azimuth_angle(i,j))
		lowvariability_confidence_test = .false.
		na_band_used = 0


! we have to set cloudy/not cloudy and surface type outside the cloud phase call so the retrieval actually works. 
		if (cloudmask(i,j)%cloudmask_determined) cloudsummary(i,j)%cloudmask_determined = .true.

		if (cloudmask(i,j)%confident_cloudy .or. cloudmask(i,j)%probablyclear_66) &
			cloudsummary(i,j)%cloudobserved = .true.

		if (cloudmask(i,j)%snowice_surface) then
			if (cloudmask(i,j)%land_surface) cloudsummary(i,j)%land_surface = .true.
			if (cloudmask(i,j)%water_surface) cloudsummary(i,j)%ocean_surface = .true.
			if (cloudmask(i,j)%desert_surface) cloudsummary(i,j)%desert_surface = .true.
			if (cloudmask(i,j)%coastal_surface) cloudsummary(i,j)%coastal_surface = .true.
		endif
		if (cloudmask(i,j)%land_surface .or. cloudmask(i,j)%coastal_surface .or. cloudmask(i,j)%desert_surface) then
			if (cloudsummary(i,j)%ocean_surface) then
				cloudsummary(i,j)%coastal_surface= .true.
				cloudsummary(i,j)%ocean_surface = .false.
			endif
		endif


                if (cloudsummary(i,j)%cloudmask_determined .and. cloudsummary(i,j)%cloudobserved .and. &
                                cloud_top_pressure(i,j) < 0. ) then 
                        retrievalstatus = 1
                        call assign_retrieval_error(i,j)
                        retrieval_failcount = retrieval_failcount + 1
                        cycle
                endif



                     
		corr_meas = fillvalue_real



! We have a cloud and a valid cloud top pressure, now we can attempt retrieval

		if (cloudsummary(i,j)%cloudobserved ) then


			if (iterationX == 1) then
				temp_meas = band_measurements(i, :, j)
				corr_meas = band_measurements(i, :, j)
			else
				temp_meas = band_measurements(i+1, :, j)
				corr_meas = band_measurements(i+1, :, j)
			endif

! DO_COX_MUNK is a model control flag. Set/unset this in mod06_run_settings.f90
			if (cloudsummary(i,j)%ocean_surface .and. .not. cloudsummary(i,j)%snowice_surface &
							.and. temp_meas(2) > 0. .and. DO_COX_MUNK) then 
				COX_MUNK = .true.
			else
				COX_MUNK = .false.
			endif


			const_C = pi / ( cos(solar_zenith_angle(i,j)*d2r) * solar_constant_37)
			
			dsol = threshold_solar_zenith
			dsen = threshold_sensor_zenith

			diff_scat_angle = abs(scattering_angle - lastinterp_scat_angle)
			diff_scat_angle_ss = abs(scattering_angle - lastinterp_scat_angle_ss)
		
			if (scattering_angle < 60. .or. (scattering_angle > 120. .and. scattering_angle <= 160.)) then 
				dscat = 2.0
				dsen = 1.0
			else  if (scattering_angle >= 60 .and. scattering_angle <= 120.) then 
				dscat = dscat3
			else ! more than 160 degrees
				dscat = 1.0
				dsen = 0.5 
			endif
			
		
			diff_solar_zenith = abs(solar_zenith_angle(i,j) - lastinterp_solar_zenith)
			diff_sensor_zenith= abs(sensor_zenith_angle(i,j)- lastinterp_sensor_zenith)
			diff_relative_azimuth= abs(relative_azimuth_angle(i,j)-lastinterp_relative_azimuth)
		 
			call get_model_idx(latitude(i,j), longitude(i,j), model_i, model_j)
			cur_wind_speed = model_info(model_i, model_j)%wind_speed
		
			diff_wind_speed = abs(cur_wind_speed - lastinterp_wind_speed)
		 
		 
			if (COX_MUNK .and. diff_wind_speed > threshold_wind_speed .and. &
				.not. (diff_solar_zenith > dsol) .and. &
				.not. (diff_sensor_zenith > dsen) .and. &
				.not. (diff_relative_azimuth > drel) .and. &
				.not. (diff_scat_angle > dscat) .and. &
				.not. (COX_MUNK .neqv. last_COX_MUNK)) then 
					cnt_wspeed_only = cnt_wspeed_only + 1
					wind_speed_only = .true.
			else
					wind_speed_only = .false. 
			endif

! interpolate the libraries
			interp_MS = .false.

			if ( diff_solar_zenith     > dsol  .or. &
				diff_sensor_zenith    > dsen .or. &
				diff_relative_azimuth > drel .or. &
				diff_scat_angle > dscat .or. &
				(COX_MUNK .neqv. last_COX_MUNK) .or. &
				(COX_MUNK .and. diff_wind_speed > threshold_wind_speed)) then 
				
					interp_MS = .true.

				if (librarystatus == 0 ) count_interpolations = count_interpolations+1
				if (librarystatus /= 0)  library_failcount = library_failcount+1

				if (diff_scat_angle > dscat) &
					lastinterp_scat_angle = scattering_angle
			
				if (diff_solar_zenith > dsol) &
					lastinterp_solar_zenith = solar_zenith_angle(i,j)
				if (diff_sensor_zenith > dsen) &
					lastinterp_sensor_zenith = sensor_zenith_angle(i,j)
				if (diff_relative_azimuth > drel) &
					lastinterp_relative_azimuth = relative_azimuth_angle(i,j)
			
				if (COX_MUNK .and. diff_wind_speed > threshold_wind_speed) &
					lastinterp_wind_speed = cur_wind_speed

				if (COX_MUNK .neqv. last_COX_MUNK) &
					last_COX_MUNK = COX_MUNK


			endif

			interp_SS = .false.

			if ((diff_scat_angle_ss > 0.1) .or. (COX_MUNK .neqv. last_COX_MUNK) .or. &
               (COX_MUNK .and. diff_wind_speed > threshold_wind_speed)) then				
               	interp_SS = .true.
				lastinterp_scat_angle_ss = scattering_angle
			end if

#ifdef RETRIEVE
				call libraryinterpolate(solar_zenith_angle(i,j),    &
                                  sensor_zenith_angle(i,j),   &
                                  relative_azimuth_angle(i,j), &
								  scattering_angle, &
								  cur_wind_speed, &
								  wind_speed_only, interp_MS, interp_SS, &
								  debug, &
                                  librarystatus, i, j)
#endif

            ctt_1km = cloud_top_temperature(i,j)	!(save 1km ctt for to pass to clouddecision)

			
! do the IRW retrieval
			irw_temperature(i,j) = fillvalue_real

			if (COX_MUNK) then 
				call get_rtm_parameters(platform_name, int_surface_emissivity_water(1,:,1), &
							sensor_zenith_angle(i,j), solar_zenith_angle(i,j),  model_i, model_j, i, j)
			else
				call get_rtm_parameters(platform_name, surface_emissivity_land(i,j,:),  &
							sensor_zenith_angle(i,j), solar_zenith_angle(i,j), model_i, model_j, i, j)			
			endif



			if (cloud_height_method(i,j) == 6) then 
			! retrieve regular temperature
				call retrieve_irw_temp(i,j, temp_meas(band_1100), &
						model_i, model_j, rtm_rad_atm_clr, rtm_trans_atm_clr, rtm_cloud_prof, &
						irw_temperature(i,j), irw_pressure, irw_dummy)		
				cloud_top_temperature(i,j) = irw_temperature(i,j)
				cloud_top_pressure(i,j) = irw_pressure
			! now do the stuff for uncertainty
				call retrieve_irw_temp(i,j, temp_meas(band_1100), &
						model_i, model_j, rtm_rad_atm_clr_low, rtm_trans_atm_clr_low, rtm_cloud_prof_low , &
						Tc_low_for_delta, temp_pres, irw_dummy )
				call retrieve_irw_temp(i,j, temp_meas(band_1100), &
						model_i, model_j, rtm_rad_atm_clr_high, rtm_trans_atm_clr_high, rtm_cloud_prof_high, &
						Tc_high_for_delta, temp_pres, irw_dummy )
				
			else if (cloud_height_method(i,j) > 0 .and. cloud_height_method(i,j) < 6) then 
				irw_temperature(i,j) = fillvalue_real
			! for uncertainty we need to find the temperature that fits a delta_P of 50 mb
				call given_P_get_T(cloud_top_pressure(i,j)-delta_Pc, model_info(model_i, model_j), Tc_low_for_delta)
				if (cloud_top_pressure(i,j) + delta_Pc > model_info(model_i, model_j)%Ps) then 
					Tc_high_for_delta = surface_temperature(i,j)
				else 
					call given_P_get_T(cloud_top_pressure(i,j)+delta_Pc, model_info(model_i, model_j), Tc_high_for_delta)
				endif
			endif

			if (cloud_top_pressure(i,j) < 0.) then
				cloud_top_pressure(i,j) = model_info(model_i, model_j)%Ps
				cloud_top_temperature(i,j) = surface_temperature(i,j)
			endif

! get above-cloud water vapor				
			call get_above_cloud_properties(model_info(model_i,model_j)%pressure_profile(:),&   
											model_info(model_i,model_j)%mixr_profile(:),  &   
											model_info(model_i,model_j)%surface_level, &
											cloud_top_pressure(i,j),   &
											abovecloud_watervapor(i,j), &
											status )



			Tc_liquid = fillvalue_real
			Tc_ice = fillvalue_real




! do atmospheric correction
			
			
			
			call atmospheric_correction(i,j, iterationX, corr_meas, model_info(model_i,model_j), &
								debug, atmoscorrstatus)

! now we need to compute derivatives that we'll hang on to in the uncertainty calculations
			Bprime_Tc = ( modis_planck(platform_name, Tc_high_for_delta, channel_37um, 1) - &
							modis_planck(platform_name, Tc_low_for_delta, channel_37um, 1)) / &
								(abovecloud_watervapor(i,j)*(2.*watervapor_error)) 

			Bprime_Ts = ( modis_planck(platform_name, surface_temperature(i,j) + delta_Ts, channel_37um, 1) - &
							 modis_planck(platform_name, surface_temperature(i,j) - delta_Ts, channel_37um, 1)) / &
								(2.0*delta_Ts) 
									
			Transprime_1way = ( thermal_correction_oneway_high(1) - thermal_correction_oneway_low(1)) / &
									(abovecloud_watervapor(i,j)*(2.*watervapor_error))
			Transprime_2way = ( thermal_correction_twoway_high(1) - thermal_correction_twoway_low(1)) / &
									(abovecloud_watervapor(i,j)*(2.*watervapor_error))



			albedo_real4 = surface_albedo(i,j,:)*albedo_fac
			if (COX_MUNK) albedo_real4 = 0.
			
			RSS_ice = -999.
			RSS_liq = -999.
			rayleigh_liq = fillvalue_real
			rayleigh_ice = fillvalue_real
			alt_ray_liq = fillvalue_real
			alt_ray_ice = fillvalue_real
			liq_near = 0
			ice_near = 0

			optical_thickness_liquid = fillvalue_real
			optical_thickness_16_liquid = fillvalue_real
			optical_thickness_37_liquid = fillvalue_real
		    optical_thickness_1621_liquid = fillvalue_real
		    effective_radius_16_liquid = fillvalue_real
		    effective_radius_21_liquid = fillvalue_real
  	    	effective_radius_37_liquid = fillvalue_real
 		 	effective_radius_1621_liquid = fillvalue_real

			optical_thickness_ice = fillvalue_real
			optical_thickness_16_ice = fillvalue_real
			optical_thickness_37_ice = fillvalue_real
  			optical_thickness_1621_ice = fillvalue_real
  			effective_radius_16_ice = fillvalue_real
  			effective_radius_21_ice = fillvalue_real
  			effective_radius_37_ice = fillvalue_real
  			effective_radius_1621_ice = fillvalue_real
  			
  			nearest_used = 0
  			RSS_final = fillvalue_real

#ifdef RETRIEVE
			
			call corescience (i, j, cloudsummary(i,j), corr_meas, &
					Tc_liquid, Tc_ice, &
					debug, na_band_used, liq_near, ice_near, &
					RSS_liq, RSS_ice, alt_ray_liq, alt_ray_ice, &
					retrievalstatus)

				
#endif


			if (COX_MUNK) &
				call set_cox_munk_albedo (albedo_real4(:), int_reflectance_water(1,:,1))


	
			if (retrievalstatus == 0 ) retrievalcount = retrievalcount+1
			if (retrieval_failcount /= 0) retrieval_failcount = retrieval_failcount +1


! there was no cloud
		else
         ! failure before retrieval 
			retrievalstatus = 1
			call assign_retrieval_error(i,j)
			retrieval_failcount = retrieval_failcount + 1
		endif




! now we do cloud phase
		cloudsummary(i,j)%cloudmask_determined = .true.
		cloudsummary(i,j)%cloudobserved = .false.
		cloudsummary(i,j)%watercloud = .false.
		cloudsummary(i,j)%icecloud = .false.
		cloudsummary(i,j)%unknowncloud = .false.

!      set the baum phase according to the "Cloud_Phase_Infrared_1km" SDS  (to pass to cloud phase and multi-layer alg.)

        ir_cloudphase%icecloud       = 0
        ir_cloudphase%watercloud     = 0
        ir_cloudphase%unknowncloud   = 1
		if (cloud_phase_infrared(i,j) == 1) ir_cloudphase%watercloud = 1
		if (cloud_phase_infrared(i,j) == 2) ir_cloudphase%icecloud = 1
		if (cloud_phase_infrared(i,j) == 1 .or. cloud_phase_infrared(i,j) == 2) ir_cloudphase%unknowncloud = 0

		call clouddecision(platform_name,                       &
							cloudmask(i,j),                     &
							corr_meas,                          &
							RSS_liq,                            &
							RSS_ice,                            &
                            optical_thickness_liquid,           & 
							optical_thickness_ice,              & 
							effective_radius_16_liquid,     	& 
							effective_radius_21_liquid,     	& 
							effective_radius_37_liquid,     	& 
							effective_radius_16_ice,     	    & 
							effective_radius_21_ice,     	    & 
							effective_radius_37_ice,     	    & 								
 							ctt_1km,                            &
 							cloud_mask_SPI(2,i,j)*0.01,              &
 							cloud_height_method(i,j),           &
 							ir_cloudphase,                      &
                            processing_information(i,j)%band_used_for_optical_thickness, &
                            cloudsummary(i,j),                  &
                            ldebug,                             &
							cloudstatus, i,j) 


! force phase to ice
		if (FORCE_ICE .and. cloudsummary(i,j)%cloudobserved) then 
			cloudsummary(i,j)%watercloud = .false.
			cloudsummary(i,j)%icecloud = .false.
			cloudsummary(i,j)%unknowncloud = .false.
			cloudsummary(i,j)%icecloud = .true.
		endif
! force phase to water
		if (FORCE_WATER .and. cloudsummary(i,j)%cloudobserved) then 
			cloudsummary(i,j)%watercloud = .false.
			cloudsummary(i,j)%icecloud = .false.
			cloudsummary(i,j)%unknowncloud = .false.
			cloudsummary(i,j)%watercloud = .true.
		endif


		if (cloudsummary(i,j)%cloudobserved) then 

! the channels get set regardless of phase, however
! 0.65um gets overwritten if rayleigh is applied later
				atm_corr_refl(band_0065,i,j) = corr_meas(band_0065) 
				atm_corr_refl(band_0086,i,j) = corr_meas(band_0086) 
				atm_corr_refl(band_0124,i,j) = corr_meas(band_0124)
				atm_corr_refl(band_0163,i,j) = corr_meas(band_0163)
				atm_corr_refl(band_0213,i,j) = corr_meas(band_0213)
				atm_corr_refl(band_0370-1,i,j) = fillvalue_real



! set the answers based on final phase decision and do the remaining science here
! we need to compute water path, multilayer and uncertainty and set the tau_out_of_bounds QA bit


			if (cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud) then 
				
! set the liquid water answers
				optical_thickness_final(i,j) = optical_thickness_liquid
				optical_thickness_1621_final(i,j) = optical_thickness_1621_liquid
				optical_thickness_16_final(i,j) = optical_thickness_16_liquid
				optical_thickness_37_final(i,j) = optical_thickness_37_liquid
				effective_radius_16_final(i,j) = effective_radius_16_liquid
				effective_radius_21_final(i,j) = effective_radius_21_liquid
				effective_radius_37_final(i,j) =  effective_radius_37_liquid
				effective_radius_1621_final(i,j) = effective_radius_1621_liquid
				if (Tc_liquid > 0.) then 
					irw_temperature(i,j) = Tc_liquid
					cloud_top_temperature(i,j) = Tc_liquid
				endif

				nearest_used = liq_near
				RSS_final = RSS_liq
				emission_pw = emission_uncertainty_pw_liq
				emission_Tc = emission_uncertainty_Tc_liq	
				sigma_R37_pw = sigma_R37_PW_liq

! set the rayleigh refl here
				if (.not. COX_MUNK) then
					if (alt_ray_liq > 0.) then 
						atm_corr_refl(band_0065, i,j) = alt_ray_liq
					else
						call bisectionsearch(water_radii, effective_radius_21_liquid, &
										re_idx_low,re_idx_hi)
						if (rayleigh_liq(re_idx_low) > 0. .and. rayleigh_liq(re_idx_hi) > 0.) then 
							atm_corr_refl(band_0065, i,j) = &
								linearinterpolation( (/water_radii(re_idx_low), water_radii(re_idx_hi) /), &
										(/rayleigh_liq(re_idx_low), rayleigh_liq(re_idx_hi)/), &
													effective_radius_21_liquid)	
						endif
					endif				
				endif
! set the tau out of bounds bit

				if (optical_thickness_final(i,j) > 150.) then!{
					optical_thickness_final(i,j) = 150.
					processing_information(i,j)%optical_thickness_outofbounds = 2
				else!}{
					processing_information(i,j)%optical_thickness_outofbounds = 0
				endif!}


				if (optical_thickness_16_final(i,j) > 150.) &
						optical_thickness_16_final(i,j) = 150.

				if (optical_thickness_37_final(i,j) > 150.) &
						optical_thickness_37_final(i,j) = 150.

				if (optical_thickness_1621_final(i,j) > 150.) &
						optical_thickness_1621_final(i,j) = 150.


! compute water path for regular and 1.6-2.1 retrievals			
				if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j)  > 0.) then 
					call compute_water_path(optical_thickness_final(i,j), &
										effective_radius_21_final(i,j), &
										liquid_water_density, &
										water_radii, &
										extinction_water(1, :), &
										liquid_water_path(i,j))
				else
					liquid_water_path(i,j) = fillvalue_real
				endif

				if (optical_thickness_16_final(i,j) > 0. .and. effective_radius_16_final(i,j)  > 0.) then 
					call compute_water_path(optical_thickness_16_final(i,j), &
										effective_radius_16_final(i,j), &
										liquid_water_density, &
										water_radii, &
										extinction_water(1, :), &
										liquid_water_path_16(i,j))
				else
					liquid_water_path_16(i,j) = fillvalue_real
				endif


				if (optical_thickness_37_final(i,j) > 0. .and. effective_radius_37_final(i,j)  > 0.) then 
					call compute_water_path(optical_thickness_37_final(i,j), &
										effective_radius_37_final(i,j), &
										liquid_water_density, &
										water_radii, &
										extinction_water(1, :), &
										liquid_water_path_37(i,j))
				else
					liquid_water_path_37(i,j) = fillvalue_real
				endif
			
				if (optical_thickness_1621_final(i,j) > 0. .and. effective_radius_1621_final(i,j)  > 0.) then 
					call compute_water_path(optical_thickness_1621_final(i,j), &
										effective_radius_1621_final(i,j), &
										liquid_water_density, &
										water_radii, &
										extinction_water(1, :), &
										liquid_water_path_1621(i,j))
				else
					liquid_water_path_1621(i,j) = fillvalue_real
				endif

			else if(cloudsummary(i,j)%icecloud) then 

! set the ice cloud answers		
				optical_thickness_final(i,j) = optical_thickness_ice
				optical_thickness_1621_final(i,j) = optical_thickness_1621_ice
				optical_thickness_16_final(i,j) = optical_thickness_16_ice
				optical_thickness_37_final(i,j) = optical_thickness_37_ice
				effective_radius_16_final(i,j) = effective_radius_16_ice
				effective_radius_21_final(i,j) = effective_radius_21_ice
				effective_radius_37_final(i,j) =  effective_radius_37_ice
				effective_radius_1621_final(i,j) = effective_radius_1621_ice
				if (Tc_ice > 0.) then 
					irw_temperature(i,j) = Tc_ice
					cloud_top_temperature(i,j) = Tc_ice
				endif

				nearest_used = ice_near
				RSS_final = RSS_ice
				emission_pw = emission_uncertainty_pw_ice
				emission_Tc = emission_uncertainty_Tc_ice	
				sigma_R37_pw = sigma_R37_PW_ice
				
				if (.not. COX_MUNK) then
					if (alt_ray_ice > 0.) then 
						atm_corr_refl(band_0065, i,j) = alt_ray_ice
					else
						call bisectionsearch(ice_radii,effective_radius_21_ice, &
									re_idx_low,re_idx_hi)
						if (rayleigh_ice(re_idx_low) > 0. .and. rayleigh_ice(re_idx_hi) > 0.) then 
							atm_corr_refl(band_0065, i,j) = &
								linearinterpolation( (/ice_radii(re_idx_low), ice_radii(re_idx_hi) /), &
									(/rayleigh_ice(re_idx_low), rayleigh_ice(re_idx_hi)/), &
													effective_radius_21_ice)	
						endif
					endif				
				endif

! set the tau out of bounds bit
! the new setting indicates flagging of tau only if it's more than 150. All others are considered perfectly valid

				if (optical_thickness_final(i,j) > 150.) then!{
					optical_thickness_final(i,j) = 150.
					processing_information(i,j)%optical_thickness_outofbounds = 2
				else!}{
					processing_information(i,j)%optical_thickness_outofbounds = 0
				endif!}


				if (optical_thickness_16_final(i,j) > 150.) &
						optical_thickness_16_final(i,j) = 150.

				if (optical_thickness_37_final(i,j) > 150.) &
						optical_thickness_37_final(i,j) = 150.

				if (optical_thickness_1621_final(i,j) > 150.) &
						optical_thickness_1621_final(i,j) = 150.


! computer the ice water path for regular and 1.6-2.1 retrievals
				if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j)  > 0.) then 
					call compute_water_path(optical_thickness_final(i,j), &
										effective_radius_21_final(i,j), &
										ice_water_density, &
										ice_radii, &
										extinction_ice(1, :), &
										liquid_water_path(i,j))
				else
					liquid_water_path(i,j) = fillvalue_real
				endif

				if (optical_thickness_16_final(i,j) > 0. .and. effective_radius_16_final(i,j)  > 0.) then 
					call compute_water_path(optical_thickness_16_final(i,j), &
										effective_radius_16_final(i,j), &
										ice_water_density, &
										ice_radii, &
										extinction_ice(1, :), &
										liquid_water_path_16(i,j))
				else
					liquid_water_path_16(i,j) = fillvalue_real
				endif

				if (optical_thickness_37_final(i,j) > 0. .and. effective_radius_37_final(i,j)  > 0.) then 
					call compute_water_path(optical_thickness_37_final(i,j), &
										effective_radius_37_final(i,j), &
										ice_water_density, &
										ice_radii, &
										extinction_ice(1, :), &
										liquid_water_path_37(i,j))
				else
					liquid_water_path_37(i,j) = fillvalue_real
				endif

			
				if (optical_thickness_1621_final(i,j) > 0. .and. effective_radius_1621_final(i,j)  > 0.) then 
					call compute_water_path(optical_thickness_1621_final(i,j), &
										effective_radius_1621_final(i,j), &
										ice_water_density, &
										ice_radii, &
										extinction_ice(1, :), &
										liquid_water_path_1621(i,j))
				else
					liquid_water_path_1621(i,j) = fillvalue_real
				endif


			
			endif

			if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j) > 0.) then
				IM_successful_retrieval_count = IM_successful_retrieval_count + 1
			endif

! assign the failure metric here: 

! nearest_used and RSS_final

		if (nearest_used(re16) == 1 .and. (.not. (optical_thickness_16_final(i,j) == MAX_TAU_RTRIEVED &
										    .and. effective_radius_16_final(i,j) /= fillvalue_real))) then 
			failure_metric_16(i,j)%tau = nint(optical_thickness_16_final(i,j)/retr_scale_factor)
			if (optical_thickness_16_final(i,j) < 0.) failure_metric_16(i,j)%tau = fillvalue_int2
			failure_metric_16(i,j)%re = nint(effective_radius_16_final(i,j)/retr_scale_factor)
			if (effective_radius_16_final(i,j) < 0.) failure_metric_16(i,j)%re = fillvalue_int2
			failure_metric_16(i,j)%RSS = nint(RSS_final(re16)*100./retr_scale_factor)
			optical_thickness_16_final(i,j) = fillvalue_real
			effective_radius_16_final(i,j) = fillvalue_real
			liquid_water_path_16(i,j) = fillvalue_real			
		endif

		if (nearest_used(re37) == 1 .and. (.not. (optical_thickness_37_final(i,j) == MAX_TAU_RTRIEVED &
										    .and. effective_radius_37_final(i,j) /= fillvalue_real))) then 
			failure_metric_37(i,j)%tau = nint(optical_thickness_37_final(i,j)/retr_scale_factor)
			if (optical_thickness_37_final(i,j) < 0) failure_metric_37(i,j)%tau = fillvalue_int2
			failure_metric_37(i,j)%re = nint(effective_radius_37_final(i,j)/retr_scale_factor)
			if (effective_radius_37_final(i,j) < 0.) failure_metric_37(i,j)%re = fillvalue_int2
			failure_metric_37(i,j)%RSS = nint(RSS_final(re37)*100./retr_scale_factor)
			optical_thickness_37_final(i,j) = fillvalue_real
			effective_radius_37_final(i,j) = fillvalue_real
			liquid_water_path_37(i,j) = fillvalue_real			
		endif

		if (nearest_used(re1621) == 1) then 
			failure_metric_1621(i,j)%tau = nint(optical_thickness_1621_final(i,j)/retr_scale_factor)
			if (optical_thickness_1621_final(i,j) < 0) failure_metric_1621(i,j)%tau = fillvalue_int2
			failure_metric_1621(i,j)%re = nint(effective_radius_1621_final(i,j)/retr_scale_factor)
			if (effective_radius_1621_final(i,j) < 0.) failure_metric_1621(i,j)%re = fillvalue_int2
			failure_metric_1621(i,j)%RSS = nint(RSS_final(re1621)*100./retr_scale_factor)
			optical_thickness_1621_final(i,j) = fillvalue_real
			effective_radius_1621_final(i,j) = fillvalue_real
			liquid_water_path_1621(i,j) = fillvalue_real			
		endif

		if (nearest_used(re21) == 1 .and. (.not. (optical_thickness_final(i,j) == MAX_TAU_RTRIEVED &
										    .and. effective_radius_21_final(i,j) /= fillvalue_real))) then 
			failure_metric(i,j)%tau = nint(optical_thickness_final(i,j)/retr_scale_factor)
			if (optical_thickness_final(i,j) < 0) failure_metric(i,j)%tau = fillvalue_int2
			failure_metric(i,j)%re = nint(effective_radius_21_final(i,j)/retr_scale_factor)
			if (effective_radius_21_final(i,j) < 0.) failure_metric(i,j)%re = fillvalue_int2
			failure_metric(i,j)%RSS = nint(RSS_final(re21)*100./retr_scale_factor)
			effective_radius_21_final(i,j) = fillvalue_real
			optical_thickness_final(i,j) = fillvalue_real
			liquid_water_path(i,j) = fillvalue_real			
		endif





! compute multilayer
			if (optical_thickness_final(i,j) > 0.) then 
				
					call Compute_Multilayer_Map(platform_name,  &
									transmit_correction_table,             &   
									temp_meas,   &   
									cloudsummary(i,j),           &   
									ir_cloudphase,          &   
									model_info(model_i,model_j)%pressure_profile,&   
									model_info(model_i,model_j)%mixr_profile(:),  &   
									model_info(model_i,model_j)%temp_profile(:),   &
									model_info(model_i,model_j)%surface_level, &
									cloud_top_pressure(i,j),     &   
									abovecloud_watervapor(i,j),  &   
									sensor_zenith_angle(i,j),                   &   
									solar_zenith_angle(i,j),                    &  
									relative_azimuth_angle(i,j), &
									optical_thickness_final(i,j),&
									optical_thickness_1621_final(i,j), &
									i, j,           &   
									cloud_layer_flag(i,j), ml_test_flag(i,j))              
			
			else
				cloud_layer_flag(i,j) = 0
				ml_test_flag(i,j) = 0
			endif





! *** ATTENTION **** 
! to do the 3.7um uncertainty, it is not enough to replace the re and tau with the 3.7um values. It is also 
! necessary to set the absorbingband_index to be 3.7um to feed the libraries in. In addition to that
! the albedo_holder and R1R2wavelengthIdx arrays MUST be fed with absorbingband_index-1 !!
! If you fail to do so, you will end up with a royal mess and not know why the numbers don't make sense. 
! -- G. Wind 7.5.2006


!  get retrieval uncertainty estimate
			if ((nearest_used(re21) == 0 .or. (nearest_used(re21) == 1 .and. optical_thickness_final(i,j) == MAX_TAU_RTRIEVED ))&
			 .and. (optical_thickness_final(i,j) .ge. 0.01) .and. (effective_radius_21_final(i,j) .ge. 0.01) .and. &
				(cloudsummary(i,j)%icecloud .or. cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud)) then!{
 
 
				absorbingband_index = band_0213
 
				if ( cloudsummary(i,j)%icecloud ) then!{
					phase = 'ice'
				else!}{
					phase = 'water'
				endif!}

				albedo_holder =  (/albedo_real4(na_band_used), &
							albedo_real4(absorbingband_index)/)
				cloud_reflectance = (/corr_meas(na_band_used), &
							      corr_meas(absorbingband_index)/) 
				R1R2wavelengthIdx = (/na_band_used, absorbingband_index/)
			
				if (iterationX == 1) then 
					uncertain_start = i
				else 
					uncertain_start	= i+1
				endif

				unc_reflectance(1) = spec_uncertain(na_band_used) * &
								exp (band_uncertainty(uncertain_start,na_band_used, j)*1.0 / uncertain_sf(na_band_used)) * 0.01
				unc_reflectance(2) = spec_uncertain(absorbingband_index) * &
								exp (band_uncertainty(uncertain_start,absorbingband_index, j)*1.0 / uncertain_sf(absorbingband_index)) * 0.01


				if (set_bands(na_band_used) <= 4) unc_reflectance(1) = max (VNIR_error, unc_reflectance(1))
				if (set_bands(na_band_used) >= 5) unc_reflectance(1) = max (SWIR_error, unc_reflectance(1))
				if (set_bands(absorbingband_index) >= 5) unc_reflectance(2) = max(SWIR_error, unc_reflectance(2))

! FIVE PERCENT			
!				unc_reflectance = 0.05
! UNC_REFL + 2%
!				unc_reflectance = sqrt ( unc_reflectance**2 + 0.02**2 )

				call getuncertainties(	optical_thickness_final(i,j), &
									effective_radius_21_final(i,j),         &
									liquid_water_path(i,j), &
									phase,                     & 
									R1R2wavelengthIdx,         &
									unc_reflectance, &
									albedo_holder,             &
									transmittance_twoway(na_band_used), &
									transmittance_twoway(absorbingband_index), &
									meandelta_trans(na_band_used), &
									meandelta_trans(absorbingband_index), &
									transmittance_stddev(na_band_used), &
									transmittance_stddev(absorbingband_index), &
									emission_pw, emission_Tc, sigma_R37_pw, &
									unc_tau_real ,    &
									unc_re21_real,   &
									unc_lwp21_real, i, j)

				optical_thickness_error(i, j) = nint(unc_tau_real / unc_scale_factor)
				effective_radius_21_error(i, j) = nint(unc_re21_real / unc_scale_factor)
				liquid_water_path_error(i,j)      = nint(unc_lwp21_real / unc_scale_factor)
							

			else!}{
				optical_thickness_error(i, j) = fillvalue_int2
				effective_radius_21_error(i, j) = fillvalue_int2
				liquid_water_path_error(i,j)      = fillvalue_int2
			endif!}

			if ( unc_tau_real .lt. epsilon(unc_tau_real) .or.  &
				unc_re21_real .lt. epsilon(unc_re21_real)  .or. &
				unc_lwp21_real .lt. epsilon(unc_lwp21_real)) then!{ 
			
				optical_thickness_error(i, j) = fillvalue_int2
				effective_radius_21_error(i, j) = fillvalue_int2
				liquid_water_path_error(i,j)       = fillvalue_int2
			endif!}




! get uncertainty estimate for 1.6um retrieval
			if ((nearest_used(re16) == 0 .or. (nearest_used(re16) == 1 .and. optical_thickness_16_final(i,j) == MAX_TAU_RTRIEVED ))&
				.and. (optical_thickness_16_final(i,j) .ge. 0.01) .and. (effective_radius_16_final(i,j) .ge. 0.01) .and. &
				(cloudsummary(i,j)%icecloud .or. cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud)) then!{
 

 
				absorbingband_index = band_0163
 
				if ( cloudsummary(i,j)%icecloud ) then!{
					phase = 'ice'
				else!}{
					phase = 'water'
				endif!}

				albedo_holder =  (/albedo_real4(na_band_used), &
							albedo_real4(absorbingband_index)/)
				cloud_reflectance = (/corr_meas(na_band_used), &
							      corr_meas(absorbingband_index)/) 
				R1R2wavelengthIdx = (/na_band_used, absorbingband_index/)
			
				if (iterationX == 1) then 
					uncertain_start = i
				else 
					uncertain_start	= i+1
				endif

				unc_reflectance(1) = spec_uncertain(na_band_used) * &
								exp (band_uncertainty(uncertain_start,na_band_used, j)*1.0 / uncertain_sf(na_band_used)) * 0.01
				unc_reflectance(2) = spec_uncertain(absorbingband_index) * &
								exp (band_uncertainty(uncertain_start,absorbingband_index, j)*1.0 / uncertain_sf(absorbingband_index)) * 0.01
	
				if (set_bands(na_band_used) <= 4) unc_reflectance(1) = max (VNIR_error, unc_reflectance(1))
				if (set_bands(na_band_used) >= 5) unc_reflectance(1) = max (SWIR_error, unc_reflectance(1))
				if (set_bands(absorbingband_index) >= 5) unc_reflectance(2) = max(SWIR_error, unc_reflectance(2))

				call getuncertainties(	optical_thickness_16_final(i,j), &
									effective_radius_16_final(i,j),         &
									liquid_water_path_16(i,j), &
									phase,                     & 
									R1R2wavelengthIdx,         &
									unc_reflectance, &
									albedo_holder,             &
									transmittance_twoway(na_band_used), &
									transmittance_twoway(absorbingband_index), &
									meandelta_trans(na_band_used), &
									meandelta_trans(absorbingband_index), &
									transmittance_stddev(na_band_used), &
									transmittance_stddev(absorbingband_index), &
									emission_pw, emission_Tc, sigma_R37_pw, &									
									unc_tau16_real ,    &
									unc_re16_real,   &
									unc_lwp16_real, i, j)
				
				optical_thickness_16_error(i, j) = nint(unc_tau16_real / unc_scale_factor)
				effective_radius_16_error(i, j) = nint(unc_re16_real / unc_scale_factor)
				liquid_water_path_16_error(i, j) = nint(unc_lwp16_real / unc_scale_factor)

			else!}{
				optical_thickness_16_error(i, j) = fillvalue_int2
				effective_radius_16_error(i, j) = fillvalue_int2
				liquid_water_path_16_error(i, j) = fillvalue_int2
			endif!}

			if ( unc_tau16_real .lt. epsilon(unc_tau16_real) .or.  &
				unc_re16_real .lt. epsilon(unc_re16_real) .or. &
				unc_lwp16_real .lt. epsilon(unc_lwp16_real) ) then!{ 
				optical_thickness_16_error(i, j) = fillvalue_int2
				effective_radius_16_error(i, j) = fillvalue_int2
				liquid_water_path_16_error(i, j) = fillvalue_int2
			endif!}



! get uncertainty estimate for 3.7um retrieval
! this is the initial part, without any emission uncertainty
! we are not using the transmittance table to do the atmospheric correction here, so 
! for the moment uncertainty due to PW table for 3.7um is set to be 0.0 
			if ((nearest_used(re37) == 0 .or. (nearest_used(re37) == 1 .and. optical_thickness_37_final(i,j) == MAX_TAU_RTRIEVED ))&
				.and. (optical_thickness_37_final(i,j) .ge. 0.01) .and. (effective_radius_37_final(i,j) .ge. 0.01) .and. &
				(cloudsummary(i,j)%icecloud .or. cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud)) then!{
				

				absorbingband_index = band_0370
 
				if ( cloudsummary(i,j)%icecloud ) then!{
					phase = 'ice'
				else!}{
					phase = 'water'
				endif!}

				albedo_holder =  (/albedo_real4(na_band_used), &
							albedo_real4(absorbingband_index-1)/)
				cloud_reflectance = (/corr_meas(na_band_used), &
							      corr_meas(absorbingband_index)/) 
				R1R2wavelengthIdx = (/na_band_used, absorbingband_index-1/)
			
				if (iterationX == 1) then 
					uncertain_start = i
				else 
					uncertain_start	= i+1
				endif
			
				unc_reflectance(1) = spec_uncertain(na_band_used) * &
								exp (band_uncertainty(uncertain_start,na_band_used, j)*1.0 / uncertain_sf(na_band_used)) * 0.01
				unc_reflectance(2) = spec_uncertain(absorbingband_index-1) * &
								exp (band_uncertainty(uncertain_start,absorbingband_index-1, j)*1.0 / &
									uncertain_sf(absorbingband_index-1)) * 0.01

				if (set_bands(na_band_used) <= 4) unc_reflectance(1) = max (VNIR_error, unc_reflectance(1))
				if (set_bands(na_band_used) >= 5) unc_reflectance(1) = max (SWIR_error, unc_reflectance(1))
				if (set_bands(absorbingband_index) >= 5 .and. set_bands(absorbingband_index) <= 7) &
															unc_reflectance(2) = max(SWIR_error, unc_reflectance(2))

				call getuncertainties(	optical_thickness_37_final(i,j), &
									effective_radius_37_final(i,j),         &
									liquid_water_path_37(i,j), &
									phase,                     & 
									R1R2wavelengthIdx,         &
									unc_reflectance, &
									albedo_holder,             &
									transmittance_twoway(na_band_used), &
									transmittance_twoway(absorbingband_index), &
									meandelta_trans(na_band_used), &
									meandelta_trans(absorbingband_index), &
									transmittance_stddev(na_band_used), &
									transmittance_stddev(absorbingband_index), &
									emission_pw, emission_Tc, sigma_R37_pw,&
									unc_tau37_real ,    &
									unc_re37_real,   &
									unc_lwp37_real, i, j)

				optical_thickness_37_error(i, j) = nint(unc_tau37_real / unc_scale_factor)
				effective_radius_37_error(i, j) = nint(unc_re37_real / unc_scale_factor)
				liquid_water_path_37_error(i, j) = nint(unc_lwp37_real / unc_scale_factor)

			else!}{
				optical_thickness_37_error(i, j) = fillvalue_int2
				effective_radius_37_error(i, j) = fillvalue_int2
				liquid_water_path_37_error(i, j) = fillvalue_int2
			endif!}

			if ( unc_tau37_real .lt. epsilon(unc_tau37_real) .or.  &
			    unc_re37_real .lt. epsilon(unc_re37_real) .or. &
				unc_lwp37_real .lt. epsilon(unc_lwp37_real) ) then!{ 
				optical_thickness_37_error(i, j) = fillvalue_int2
				effective_radius_37_error(i, j) = fillvalue_int2
				liquid_water_path_37_error(i, j) = fillvalue_int2
			endif!}



!  get 1621 retrieval uncertainty estimate
			if ((nearest_used(re1621) == 0 .or. (nearest_used(re1621) == 1 .and. optical_thickness_1621_final(i,j) == MAX_TAU_RTRIEVED ))&
				.and.  (optical_thickness_1621_final(i,j) .ge. 0.01) .and. (effective_radius_1621_final(i,j) .ge. 0.01) .and. &
				(cloudsummary(i,j)%icecloud .or. cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud)) then!{
 		

				uncertainty_nonabsorbing_1621 = band_0163
				absorbingband_index = band_0213

				if ( cloudsummary(i,j)%icecloud ) then!{
					phase = 'ice'
				else!}{
					phase = 'water'
				endif!}

				albedo_holder =  (/albedo_real4(uncertainty_nonabsorbing_1621), &
								albedo_real4(absorbingband_index)/)
				cloud_reflectance = (/corr_meas(uncertainty_nonabsorbing_1621), &
							corr_meas(absorbingband_index)/) 
				R1R2wavelengthIdx = (/uncertainty_nonabsorbing_1621, absorbingband_index/)
			
			
				if (iterationX == 1) then 
					uncertain_start = i
				else 
					uncertain_start	= i+1
				endif

				unc_reflectance(1) = spec_uncertain(uncertainty_nonabsorbing_1621) * &
								exp (band_uncertainty(uncertain_start,uncertainty_nonabsorbing_1621, j)*1.0 / &
								uncertain_sf(uncertainty_nonabsorbing_1621)) * 0.01
				unc_reflectance(2) = spec_uncertain(absorbingband_index) * &
								exp (band_uncertainty(uncertain_start,absorbingband_index, j)*1.0 / uncertain_sf(absorbingband_index)) * 0.01
			
				if (set_bands(uncertainty_nonabsorbing_1621) <= 4) unc_reflectance(1) = max (VNIR_error, unc_reflectance(1))
				if (set_bands(uncertainty_nonabsorbing_1621) >= 5) unc_reflectance(1) = max (SWIR_error, unc_reflectance(1))
				if (set_bands(absorbingband_index) >= 5) unc_reflectance(2) = max(SWIR_error, unc_reflectance(2))

				call getuncertainties(optical_thickness_1621_final(i,j), &
									effective_radius_1621_final(i,j),         &
									liquid_water_path_1621(i,j), &
									phase,                     & 
									R1R2wavelengthIdx,         &
									unc_reflectance, &
									albedo_holder,             &
									transmittance_twoway(uncertainty_nonabsorbing_1621), &
									transmittance_twoway(absorbingband_index), &
									meandelta_trans(uncertainty_nonabsorbing_1621), &
									meandelta_trans(absorbingband_index), &
									transmittance_stddev(uncertainty_nonabsorbing_1621), &
									transmittance_stddev(absorbingband_index), &
									emission_pw, emission_Tc, sigma_R37_pw,&
									unc_tau_1621_real ,    &
									unc_re1621_real,   &
									unc_lwp1621_real, i, j)

			

				optical_thickness_1621_error(i, j) = nint(unc_tau_1621_real / unc_scale_factor)
				effective_radius_1621_error(i, j) = nint(unc_re1621_real / unc_scale_factor)
				liquid_water_path_1621_error(i,j) = nint(unc_lwp1621_real / unc_scale_factor)

			else!}{
				optical_thickness_1621_error(i, j) = fillvalue_int2
				effective_radius_1621_error(i, j) = fillvalue_int2
				liquid_water_path_1621_error(i,j) = fillvalue_int2
			endif!}

			if ( unc_tau_1621_real .lt. epsilon(unc_tau_1621_real) .or.  &
				unc_re1621_real .lt. epsilon(unc_re1621_real)  .or. &
				unc_lwp1621_real .lt. epsilon(unc_lwp1621_real)) then!{ 

				optical_thickness_1621_error(i, j) = fillvalue_int2
				effective_radius_1621_error(i, j) = fillvalue_int2
				liquid_water_path_1621_error(i,j) = fillvalue_int2
			endif!}


			if (DO_CSR) then 
	
! let us remember that band measurements are being overscanned
					if (iterationX == 1) then 
						if (i==1) then 
							istart = i
							iend = i+1
						endif
					else
						if (i==1) then 
							istart = i
							iend = i+2
						endif
					endif
	
					if (i > 1 .and. i < xdimension) then 
						istart = i-1
						iend = i+1
					endif
					
					if (iterationX < number_of_iterationsX) then 
						 if (i == xdimension) then 
							istart = i-1
							iend = i+1
						 endif
					else 
						if (i == xdimension) then 
							istart = i-1
							iend = i
						endif
					endif
					
					
					if (j == 1) then 
						jstart = j
						jend = j+1
					endif
					if (j >=2 .and. j <= (ydimension-1)) then 
						jstart = j-1
						jend = j+1
					endif
					if (j == ydimension) then 
						jstart = j-1
						jend = j
					endif

					! Check if 1km visible reflectance threshold or VIS/NIR ratio tests are applied and cloudy. Clear sky restoral
					! Part V (CSR=3) test will only be applied over ocean if vis1km_test = .true. (i.e., either one of visible
					! reflectance or VIS/NIR ratio tests are applied and cloudy).
					! KGM 3-4-13
					vis1km_test = .false.
					if ((cloudmask(i,j)%applied_visiblereflectance==1 .and. cloudmask(i,j)%test_visiblereflectance==1) .or. &
						(cloudmask(i,j)%applied_visnirratio==1 .and. cloudmask(i,j)%test_visnirratio==1) ) vis1km_test = .true.

					call cloudiness_test (cloudmask(i,j),           &
                                  cloudsummary(i,j),        &
                                  temp_meas, &
                                  band_measurements(istart:iend,:,jstart:jend), &
                                  sunglint_dust_test,       &
                                  lowvariability_confidence_test, &
                                  CSR_flag_array(i,j), latitude(i,j), &
                                  cloud_height_method(i,j), vis1km_test)	! KGM 3-4-13 GW 3.28.13
 
					if (CSR_flag_array(i,j) == 2 .and. cloudsummary(i,j)%ocean_surface) then 
					
						call compute_aod(i, j, scattering_angle, corr_meas, cur_wind_speed, aod550)

						! if the aerosol optical depth is too much then it's probably a cloud
						! and we can keep the retrieval, however we will mark it as a potentially
						! problematic cloud. 
						if (aod550 > 0.95) then ! aod550 is a ln(aod+0.01) quantity
							CSR_flag_array(i,j) = 0
						endif
					endif

			endif
			
		else
			retrievalstatus = 1
			call assign_retrieval_error(i,j)
			retrieval_failcount = retrieval_failcount + 1
		endif





	enddo
enddo

	call cleanup_retrieval

!	print*, "NUM_INTERP:", count_interpolations

!	print*, "num_sza: ", cnt_sza
!	print*, "num_vza: ", cnt_vza
!	print*, "num_raz: ", cnt_raz
!	print*, "num_scat: ", cnt_scat
!	print*, "num_cm_switch: ", cnt_cm_switch
!	print*, "num_wspeed: ", cnt_wspeed
!	print*, "** wind speed only: ", cnt_wspeed_only

! Now that we've done the clear sky restoral, we remove the edges of the clouds (ED)

	if (DO_CSR) then !{

		call remove_edge_scenes(cloudsummary, &
                           CSR_flag_array, &
                           xdimension, ydimension,  &
                           status)


!   reset cloudsummary variables for pixels "cleared" by CSR or ED, this step necessary
!     so that "cleared" pixels in cloud optical properties SDS get set to clear, and  
!     so that pertainent QA can be properly identified when set   GTA 6/7/05

!		endif


		where(CSR_flag_array == 2) 
			cloudsummary%cloudobserved = .false.
			cloudsummary%watercloud = .false.
			cloudsummary%icecloud = .false.
			cloudsummary%unknowncloud = .false.
			
			optical_thickness_final = fillvalue_real
			optical_thickness_16_final = fillvalue_real
			optical_thickness_37_final = fillvalue_real
			optical_thickness_1621_final = fillvalue_real
			effective_radius_16_final = fillvalue_real
			effective_radius_21_final = fillvalue_real
			effective_radius_37_final = fillvalue_real
			effective_radius_1621_final = fillvalue_real
			liquid_water_path = fillvalue_real
			liquid_water_path_16 = fillvalue_real
			liquid_water_path_37 = fillvalue_real
			liquid_water_path_1621 = fillvalue_real
			optical_thickness_error = fillvalue_int2
			effective_radius_21_error = fillvalue_int2
			effective_radius_16_error = fillvalue_int2
			effective_radius_37_error = fillvalue_int2
			liquid_water_path_error = fillvalue_int2
			liquid_water_path_16_error = fillvalue_int2
			liquid_water_path_37_error = fillvalue_int2
			cloud_layer_flag = 0
			ml_test_flag = 0
			optical_thickness_1621_error = fillvalue_int2
			optical_thickness_16_error = fillvalue_int2
			optical_thickness_37_error = fillvalue_int2
			effective_radius_1621_error = fillvalue_int2
			liquid_water_path_1621_error = fillvalue_int2
			
			failure_metric%tau = fillvalue_int2  
			failure_metric%re = fillvalue_int2  
			failure_metric%RSS = fillvalue_int2  

			failure_metric_16%tau = fillvalue_int2  
			failure_metric_16%re = fillvalue_int2  
			failure_metric_16%RSS = fillvalue_int2  

			failure_metric_1621%tau = fillvalue_int2  
			failure_metric_1621%re = fillvalue_int2  
			failure_metric_1621%RSS = fillvalue_int2  

			failure_metric_37%tau = fillvalue_int2  
			failure_metric_37%re = fillvalue_int2  
			failure_metric_37%RSS = fillvalue_int2  
			
			irw_temperature = fillvalue_real
			precip_water_094 = fillvalue_real
			
       end where
                 
	endif !}


! Now split off the PCL retrievals
! if it's an edge pixel or a 250m variable pixels then 
! split off the retrieval to the PCL storage and get rid of the value in 
! main retrieval arrays
    do j=1, ydimension
    	do i=1, xdimension
    	
    		if (CSR_flag_array(i,j) == 1 .or. CSR_flag_array(i,j) == 3) then 
    		

				if (optical_thickness_final(i,j) > 0.) &
					optical_thickness_final_PCL(i,j) = nint(optical_thickness_final(i,j) / retr_scale_factor)
				optical_thickness_final(i,j) = fillvalue_real

				if (optical_thickness_1621_final(i,j) > 0.) &
					optical_thickness_1621_final_PCL(i,j) = nint(optical_thickness_1621_final(i,j) / retr_scale_factor)
				optical_thickness_1621_final(i,j) = fillvalue_real

				if (optical_thickness_16_final(i,j) > 0.) &
					optical_thickness_16_final_PCL(i,j) = nint(optical_thickness_16_final(i,j) / retr_scale_factor)
				optical_thickness_16_final(i,j) = fillvalue_real

				if (optical_thickness_37_final(i,j) > 0.) &
					optical_thickness_37_final_PCL(i,j) = nint(optical_thickness_37_final(i,j) / retr_scale_factor)
				optical_thickness_37_final(i,j) = fillvalue_real

				if (effective_radius_16_final(i,j) > 0.) &
					effective_radius_16_final_PCL(i,j) = nint(effective_radius_16_final(i,j) / retr_scale_factor)
				effective_radius_16_final(i,j) = fillvalue_real

				if (effective_radius_37_final(i,j) > 0.) &
					effective_radius_37_final_PCL(i,j) = nint(effective_radius_37_final(i,j) / retr_scale_factor)
				effective_radius_37_final(i,j) = fillvalue_real

				if (effective_radius_21_final(i,j) > 0.) &
					effective_radius_21_final_PCL(i,j) = nint(effective_radius_21_final(i,j) / retr_scale_factor)
				effective_radius_21_final(i,j) = fillvalue_real

				if (effective_radius_1621_final(i,j) > 0.) &
					effective_radius_1621_final_PCL(i,j) = nint(effective_radius_1621_final(i,j) / retr_scale_factor)
				effective_radius_1621_final(i,j) = fillvalue_real

				if (liquid_water_path_16(i,j) > 0.) &
					liquid_water_path_16_PCL(i,j) = nint(liquid_water_path_16(i,j))
				liquid_water_path_16(i,j) = fillvalue_real

				if (liquid_water_path_37(i,j) > 0.) &
					liquid_water_path_37_PCL(i,j) = nint(liquid_water_path_37(i,j))
				liquid_water_path_37(i,j) = fillvalue_real

				if (liquid_water_path(i,j) > 0.) &
					liquid_water_path_PCL(i,j) = nint(liquid_water_path(i,j))
				liquid_water_path(i,j) = fillvalue_real

				if (liquid_water_path_1621(i,j) > 0.) &
					liquid_water_path_1621_PCL(i,j) = nint(liquid_water_path_1621(i,j))
				liquid_water_path_1621(i,j) = fillvalue_real
    	
    	
    		endif
    	
    	end do
    end do



! now we need to compute the inventory metadata that relates to the cloudiness percentage
! water cloud percentage and ice cloud percentage. We need to aggregate the little suckers

   do i=1, xdimension
      do j = 1, ydimension

! first of all we need to make sure that we have a cloud
         if (solar_zenith_angle(i,j) <= solar_zenith_threshold .and. &
              cloudsummary(i,j)%cloudobserved) then 
            IM_cloudy_count = IM_cloudy_count + 1
            
! now that we're sure, we can count the ice and water cloud pixels
            if (cloudsummary(i,j)%watercloud) then 
               IM_water_cloud_count = IM_water_cloud_count + 1
            endif
            if (cloudsummary(i,j)%icecloud) then 
               IM_ice_cloud_count = IM_ice_cloud_count + 1
            endif

            if (cloudsummary(i,j)%unknowncloud) then 
               IM_undet_count = IM_undet_count + 1
            endif

          endif

      end do
   end do


!	optical_thickness_final = abovecloud_watervapor
!	effective_radius_16_final = cloud_top_pressure

#ifndef RETRIEVE
	effective_radius_21_final(:,:) = aod550_store(:,:)
#endif
	
! 	print*, abovecloud_watervapor(14,884), cloud_top_pressure(14,884), surface_temperature(14,884)
 
!	print*, optical_thickness_final(19, 1992)
!	print*, 		effective_radius_16_final (19, 1992)
!	print*, 		effective_radius_21_final (19, 1992)
!	print*, 		effective_radius_37_final(19, 1992)
	


 end subroutine scienceinterface

subroutine compute_aod(x, y, scat_ang, corr_meas, ws, aod550)

	use core_arrays
	use mod06_run_settings
	use science_parameters, only: d2r
	use nonscience_parameters
	use modis_numerical_module, only: linearinterpolation
	
	integer, intent(in) :: x, y
	real, intent(in) :: scat_ang, ws
	real, intent(inout) :: aod550
	real, dimension(:), intent(in) :: corr_meas

	external ffnet_terra, ffnet_aqua

	real*8 :: input(15), output(1)
	real*8  :: ga, sza, vza, saz, vaz, raz, sca
	real :: check_val, temp_16
	
	sza = solar_zenith_angle(x,y)
	vza = sensor_zenith_angle(x,y)
	saz = solar_azimuth_angle(x,y)
	vaz = sensor_azimuth_angle(x,y)
	raz = relative_azimuth_angle(x,y)
	
	ga = cos(sza*d2r) * cos(vza*d2r) + sin(sza*d2r) * sin(vza*d2r) *  cos(raz*d2r)
	
	sca = cos(scat_ang*d2r)*1.0d0
	saz = cos(saz*d2r)*1.0d0
	sza = cos(sza*d2r)*1.0d0
	vaz = cos(vaz*d2r)*1.0d0
	vza = cos(vza*d2r)*1.0d0
				
	input = (/	corr_meas(band_0047)*1.0d0, &
			  	corr_meas(band_0055)*1.0d0, &
			  	corr_meas(band_0065)*1.0d0, &
			  	corr_meas(band_0086)*1.0d0, &
			  	corr_meas(band_0124)*1.0d0, &
			  	corr_meas(band_0163)*1.0d0, &
			  	corr_meas(band_0213)*1.0d0, &
				sca, ga, saz, sza, vaz, vza, 0.0d0, &
				ws*1.0d0 /)
	
	if ( corr_meas(band_0163) < 0.) then 
	
		temp_16 = linearinterpolation ( (/1.24, 2.13/), (/corr_meas(band_0124), corr_meas(band_0213)/), 1.63)
		input(6) = temp_16*1.0d0
			
	endif

	if (platform_name(1:4) == "Aqua") call ffnet_aqua(input, output)
	if (platform_name(1:5) == "Terra") call ffnet_terra(input, output)

	aod550 = output(1) ! ln(aod+0.01)

end subroutine compute_aod


subroutine assign_retrieval_error(xpoint, ypoint)
   use GeneralAuxType
   use nonscience_parameters
   use core_arrays
   integer,           intent(in)       :: xpoint,ypoint

   optical_thickness_final(xpoint, ypoint)     = fillvalue_real
   optical_thickness_16_final(xpoint, ypoint)     = fillvalue_real
   optical_thickness_37_final(xpoint, ypoint)     = fillvalue_real
   optical_thickness_1621_final(xpoint, ypoint) = fillvalue_real
   effective_radius_16_final(xpoint, ypoint)   = fillvalue_real
   effective_radius_21_final(xpoint, ypoint)   = fillvalue_real
   effective_radius_1621_final(xpoint, ypoint) = fillvalue_real
   effective_radius_37_final(xpoint, ypoint)   = fillvalue_real
   liquid_water_path(xpoint, ypoint) = fillvalue_real
   liquid_water_path_16(xpoint, ypoint) = fillvalue_real
   liquid_water_path_37(xpoint, ypoint) = fillvalue_real
   liquid_water_path_1621(xpoint, ypoint) = fillvalue_real
   optical_thickness_error(xpoint, ypoint) = fillvalue_int2
   optical_thickness_16_error(xpoint, ypoint) = fillvalue_int2
   optical_thickness_37_error(xpoint, ypoint) = fillvalue_int2
   effective_radius_21_error(xpoint, ypoint) = fillvalue_int2
   effective_radius_16_error(xpoint, ypoint) = fillvalue_int2
   effective_radius_37_error(xpoint, ypoint) = fillvalue_int2
   liquid_water_path_error(xpoint, ypoint) = fillvalue_int2
   liquid_water_path_16_error(xpoint, ypoint) = fillvalue_int2
   liquid_water_path_37_error(xpoint, ypoint) = fillvalue_int2
   cloud_layer_flag(xpoint, ypoint) = 0
   ml_test_flag(xpoint, ypoint) = 0
   CSR_flag_array(xpoint, ypoint) = 0 
   optical_thickness_1621_error(xpoint, ypoint) = fillvalue_int2
   effective_radius_1621_error(xpoint, ypoint) = fillvalue_int2
   liquid_water_path_1621_error(xpoint, ypoint) = fillvalue_int2

end subroutine assign_retrieval_error

subroutine init_science_arrays
	use GeneralAuxType
  	use nonscience_parameters
  	use core_arrays

  	optical_thickness_final = fillvalue_real
  	optical_thickness_16_final = fillvalue_real
  	optical_thickness_37_final = fillvalue_real
  	optical_thickness_1621_final = fillvalue_real
  	effective_radius_16_final = fillvalue_real
  	effective_radius_21_final = fillvalue_real
  	effective_radius_37_final = fillvalue_real
  	effective_radius_1621_final = fillvalue_real
  	liquid_water_path = fillvalue_real
  	liquid_water_path_16 = fillvalue_real
  	liquid_water_path_37 = fillvalue_real
  	liquid_water_path_1621 = fillvalue_real

  	optical_thickness_final_PCL = fillvalue_int2
  	optical_thickness_16_final_PCL = fillvalue_int2
  	optical_thickness_37_final_PCL = fillvalue_int2
  	optical_thickness_1621_final_PCL = fillvalue_int2
  	effective_radius_16_final_PCL = fillvalue_int2
  	effective_radius_21_final_PCL = fillvalue_int2
  	effective_radius_37_final_PCL = fillvalue_int2
  	effective_radius_1621_final_PCL = fillvalue_int2
  	liquid_water_path_PCL = fillvalue_int2
  	liquid_water_path_16_PCL = fillvalue_int2
  	liquid_water_path_37_PCL = fillvalue_int2
  	liquid_water_path_1621_PCL = fillvalue_int2

  	optical_thickness_error = fillvalue_int2
  	optical_thickness_16_error = fillvalue_int2
  	optical_thickness_37_error = fillvalue_int2
  	effective_radius_21_error = fillvalue_int2
  	effective_radius_16_error = fillvalue_int2
  	effective_radius_37_error = fillvalue_int2
  	liquid_water_path_error = fillvalue_int2
  	liquid_water_path_16_error = fillvalue_int2
  	liquid_water_path_37_error = fillvalue_int2
 	cloud_layer_flag = 0
  	ml_test_flag = 0
  	CSR_flag_array = 0 
  	precip_water_094 = fillvalue_real
  	optical_thickness_1621_error = fillvalue_int2
  	effective_radius_1621_error = fillvalue_int2
  	liquid_water_path_1621_error = fillvalue_int2
  	irw_temperature = fillvalue_real
  
	failure_metric(:,:)%tau = fillvalue_int2  
	failure_metric(:,:)%re = fillvalue_int2  
	failure_metric(:,:)%RSS = fillvalue_int2  

	failure_metric_16(:,:)%tau = fillvalue_int2  
	failure_metric_16(:,:)%re = fillvalue_int2  
	failure_metric_16(:,:)%RSS = fillvalue_int2  

	failure_metric_1621(:,:)%tau = fillvalue_int2  
	failure_metric_1621(:,:)%re = fillvalue_int2  
	failure_metric_1621(:,:)%RSS = fillvalue_int2  

	failure_metric_37(:,:)%tau = fillvalue_int2  
	failure_metric_37(:,:)%re = fillvalue_int2  
	failure_metric_37(:,:)%RSS = fillvalue_int2  

  	atm_corr_refl = fillvalue_real
  	

end subroutine init_science_arrays

end module modis_science_module

