module corescience_module

 implicit none
 private

 public :: corescience

 contains


subroutine corescience(xpoint,ypoint,process,measurements, Tc_liquid, Tc_ice, debug, na_band_used, nearest_liq, nearest_ice, &
									RSS_liq, RSS_ice, alt_ray_liq, alt_ray_ice, status)
   
   use GeneralAuxType
   use core_arrays
   use libraryarrays
   use libraryinterpolates
   use modis_cloudstructure, only: processflag
   use modis_numerical_module, only: linearinterpolation
   use get_retrieval_uncertainty
   use science_parameters
   use nonscience_parameters
   use mod06_run_settings
   use retrieval_solution_logic
   use retrieval_prep_logic
   
   use specific_other

   implicit none
 
   logical,           intent(in)       :: debug
   type(processflag), intent(in)       :: process
   real, dimension(:), intent(in) :: measurements
   real, intent(inout) :: Tc_liquid, Tc_ice, alt_ray_liq, alt_ray_ice
   integer,           intent(in)       :: xpoint,ypoint
   integer,           intent(out)      :: status
   integer, intent(out) :: na_band_used

	integer, dimension(:), intent(inout) :: nearest_liq, nearest_ice
	real, dimension(:), intent(inout) :: RSS_liq, RSS_ice

   type(cloudphase) :: local_process
   integer                             :: nonabsorbingband_index, absorbingband_index, maxradii
   real                                :: thermal_trans_1way, thermal_trans_2way
   real, allocatable                   :: optical_thickness_vector(:), residual(:)
   real        ::  retrievalopticalthickness,    &
                   retrievalopticalthickness16,  &
                   retrievalopticalthickness37,  &
                   retrievalopticalthickness1621,&
                   retrievalradius21,            &
                   retrievalradius16,            &
                   retrievalradius37,            &
                   retrievalradius1621 

	integer :: idx21, idx37, channel_37, idx11, i, idx_alb37, idx_alb16, quality_in
	real :: curTc, newTc, ray_temp
	logical :: PRN, use_nearest
	integer :: iii

   status = 0

   thermal_trans_1way   = thermal_correction_oneway(1)   ! only the 3.7um is used
   thermal_trans_2way   = thermal_correction_twoway(1)   ! only the 3.7um is used
 
	call get_band_idx(idx21, idx37, channel_37, idx11, idx_alb37, idx_alb16)
 
   nonabsorbingband_index = band_0065
   processing_information(xpoint, ypoint)%band_used_for_optical_thickness =  0

   if (process%cloudobserved .and. (process%ocean_surface .or. process%coastal_surface)) then!{
      if( process%snowice_surface) then!{
          ! 1.2 um
          nonabsorbingband_index = band_0124
      else!}{
          ! .86um
          nonabsorbingband_index = band_0086
      endif!}
	  if (nonabsorbingband_index > 3) then 
		processing_information(xpoint, ypoint)%band_used_for_optical_thickness = nonabsorbingband_index-1	  
	  else 
		processing_information(xpoint, ypoint)%band_used_for_optical_thickness = nonabsorbingband_index	  
	  endif
   endif!}
 
   if( process%cloudobserved .and. (process%land_surface .or. process%desert_surface) ) then!{      
      if( process%snowice_surface) then!{
          nonabsorbingband_index = band_0124
      else!}{
         ! .68um 
          nonabsorbingband_index = band_0065
      endif!} 
	  if (nonabsorbingband_index > 3) then 
		processing_information(xpoint, ypoint)%band_used_for_optical_thickness = nonabsorbingband_index-1	  
	  else 
		processing_information(xpoint, ypoint)%band_used_for_optical_thickness = nonabsorbingband_index	  
	  endif
   endif!}

   if( process%cloudobserved .and.process%snowice_surface) then!{
      ! 1.2 um
	   nonabsorbingband_index = band_0124
	  if (nonabsorbingband_index > 3) then 
		processing_information(xpoint, ypoint)%band_used_for_optical_thickness = nonabsorbingband_index-1	  
	  else 
		processing_information(xpoint, ypoint)%band_used_for_optical_thickness = nonabsorbingband_index	  
	  endif
   endif!}
 
	
	na_band_used = nonabsorbingband_index
 
 	local_process%watercloud = 0
	local_process%icecloud = 0

 
! initialize all retrievals to fillvalue   
	retrievalopticalthickness = fillvalue_real
	retrievalopticalthickness1621 = fillvalue_real
	retrievalopticalthickness16 = fillvalue_real
	retrievalopticalthickness37  = fillvalue_real
	retrievalradius16 = fillvalue_real
	retrievalradius21 = fillvalue_real
	retrievalradius1621 = fillvalue_real
	retrievalradius37 = fillvalue_real


	nearest_liq = 0
	nearest_ice = 0
	
	local_process%watercloud = 1


!  set surface QA

!  if the measured surface reflectance in the non absorbing bands is greater than
!  the surface albedo , and we have a cloudy pixel then!{ we can try a retrieval

   if (process%cloudobserved ) then!{

!    check for .86um saturation
!    if the .86um saturates (wisconsin reader sets band_meas to neg. if saturation) or there can be bad noisy neg. reflectance on low end,
!    try the other .65um band (and change the QA flag)
     if (nonabsorbingband_index == band_0086 .and. measurements(nonabsorbingband_index) < 0.) then!{
       if (measurements(band_0065) > 0.) nonabsorbingband_index = band_0065
	   	na_band_used = nonabsorbingband_index
        processing_information(xpoint, ypoint)%band_used_for_optical_thickness = nonabsorbingband_index
     endif!}


	nonabsorbingband_index = na_band_used

! this is for rotten snow albedoes
	if (albedo_real4( nonabsorbingband_index) < 0. .or. &
		albedo_real4( band_0163) < 0. .or. &
		albedo_real4( band_0213) < 0.) then 

		retrievalopticalthickness = fillvalue_real
		retrievalopticalthickness1621 = fillvalue_real
		retrievalopticalthickness16 = fillvalue_real
		retrievalopticalthickness37  = fillvalue_real
		retrievalradius16 = fillvalue_real
		retrievalradius21 = fillvalue_real
		retrievalradius1621 = fillvalue_real
		retrievalradius37 = fillvalue_real
		
		return 
	
	endif


!     if the shortwave saturates, try the other shortwave band 


        allocate(optical_thickness_vector(number_waterradii))
        allocate(residual(number_waterradii))
		
		allocate(reflibA(number_taus+1, number_waterradii), reflibB(number_taus+1, number_waterradii))
		allocate(rad37lib(number_taus+1, number_waterradii))
		
		rad37lib    = -999.0
		
        call vis_nonabsorbing_science(measurements(nonabsorbingband_index), &
                              nonabsorbingband_index,                 &
                              albedo_real4(nonabsorbingband_index), &
                              library_taus,                           &
                              water_radii,                            &
                              spherical_albedo_water,                 &
                              int_fluxdownwater_sensor,               &
                              int_fluxdownwater_solar,                & 
                              int_reflectance_water,                  &
                              sensor_zenith_angle(xpoint,ypoint),     &
                              solar_zenith_angle(xpoint,ypoint),      &
                              relative_azimuth_angle(xpoint,ypoint),  &
                              cloud_top_pressure(xpoint,ypoint),      &
                              local_process,                                &
                              optical_thickness_vector)

        absorbingband_index = band_0163 ! 1.6um band                      
        if(measurements(absorbingband_index) > 0. ) then!{ 
           
           call vis_absorbing_science(optical_thickness_vector,          &
                measurements(absorbingband_index), &
                absorbingband_index,                    &
                albedo_real4(idx_alb16), &
                library_taus,                           &
                water_radii,                            &
                spherical_albedo_water,                 &
                int_fluxdownwater_sensor,               &
                int_fluxdownwater_solar,                &
                int_reflectance_water,                  &
                residual,maxradii, debug)
           
           if (maxradii > 2) then!{
              call solveretrieval(residual(1:maxradii),                   &
                   optical_thickness_vector(1:maxradii),   &
                   water_radii(1:maxradii),                &
                   retrievalradius16,                      &
                   retrievalopticalthickness16,            &
                   debug, use_nearest, quality_in)
			if (use_nearest) then 
				nearest_liq(1) = 1
				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											water_radii, &
											retrievalopticalthickness16, &
											retrievalradius16, RSS_liq(re16), &
											.true., ray_temp, quality_in )		
			endif
		
		
           else!}{
              retrievalradius16 = fillvalue_real
              retrievalopticalthickness16 = fillvalue_real
           endif!}

        else!}{
           retrievalradius16 = fillvalue_real
           retrievalopticalthickness16 = fillvalue_real
        endif!}


#ifdef SEV_PR06OD
          retrievalradius21 = fillvalue_real
		  retrievalopticalthickness = retrievalopticalthickness16		
#else
        absorbingband_index = band_0213 ! 2.1um
        call vis_absorbing_science(optical_thickness_vector,          &
                              measurements(absorbingband_index), &
                              idx21,                    &
                              albedo_real4(absorbingband_index), &
                              library_taus,                           &
                              water_radii,                            &
                              spherical_albedo_water,                 &
                              int_fluxdownwater_sensor,               &
                              int_fluxdownwater_solar,                &
                              int_reflectance_water,                  &
                              residual, maxradii, debug)


        if (maxradii > 2) then!{  
          call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                             water_radii(1:maxradii), &
                             retrievalradius21, &
                             retrievalopticalthickness, &
                   debug, use_nearest, quality_in)
				   				   
			if (use_nearest) then 
				nearest_liq(2) = 1
				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											water_radii, &
											retrievalopticalthickness, &
											retrievalradius21, RSS_liq(re21), &
											.true., alt_ray_liq, quality_in)
				

			endif
		
			


        else!}{
          retrievalradius21 = fillvalue_real
          retrievalopticalthickness = fillvalue_real
        endif!}
#endif

        absorbingband_index = band_0370 ! 3.7um

! sep, 4 May: absorbingband_index - 1 being passed in for all libarary indices because 
!             band_0370 = 7 and library index corresponding to this band is 6. In the 
!             libraries, band_0935 is index 7.
!
! wind, 7 Dec: absorbingband_index - 1 has to be passed in for albedo as well. 
!

			call nir_absorbing_science(platform_name,&
                              optical_thickness_vector,               &
                              measurements(absorbingband_index), &
                              idx37,                &
                              albedo_real4( idx_alb37), &
                              xpoint, ypoint, &
                              cloud_top_temperature(xpoint, ypoint), &
                              thermal_trans_1way,                     &
                              thermal_trans_2way,                     &
                              library_taus,                           &
                              water_radii,                            &
                              spherical_albedo_water,                 &
                              int_fluxdownwater_sensor,               &
                              int_fluxdownwater_solar,                &
                              int_fluxupwater_sensor,                 &
                              int_reflectance_water,                  &
							  int_cloud_emissivity_water, &
							  int_surface_emissivity_water, &
                              residual, maxradii, channel_37, emission_uncertainty_pw_liq, &
                              emission_uncertainty_Tc_liq, sigma_R37_PW_liq,debug)

        if ( maxradii > 2 ) then!{
          call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                             water_radii(1:maxradii), &
                             retrievalradius37, &
                             retrievalopticalthickness37, &
                   debug, use_nearest, quality_in)
				   				   
			if (use_nearest) then 
			
				nearest_liq(3) = 1
				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											water_radii, &
											retrievalopticalthickness37, &
											retrievalradius37, RSS_liq(re37), &
											.true., ray_temp, quality_in,     &
											CH37_IDX = idx37, CTopT = cloud_top_temperature(xpoint, ypoint), &
											CH37_NUM =channel_37 , platFormName= platform_name)

			endif

        else!}{
          retrievalradius37 = fillvalue_real
          retrievalopticalthickness37 = fillvalue_real
        endif!}


! now we iterate the retrieval
#ifndef AMS_INST

! if the cloud too thick, don't bother, the iteration will not accomplish absolutely anything
		if ( retrievalopticalthickness37 > 0. &
				.and. retrievalradius37 > 0. .and. &
				irw_temperature(xpoint, ypoint) > 0. .and. nearest_liq(3) == 0 ) then 

		! if the cloud is too thick, then don't bother doing anything, but still set the temperature
			if (retrievalopticalthickness37 >= 10.) then 
				Tc_liquid = cloud_top_temperature(xpoint, ypoint)

			else

			curTc = irw_temperature(xpoint, ypoint)
			

			do i=1, 5

		
				if (retrievalopticalthickness37 < 0. .or. retrievalradius37 < 0. ) then 
					retrievalradius37 = fillvalue_real
					retrievalopticalthickness37 = fillvalue_real
					Tc_liquid = curTc
					exit
				endif
				
				
				
				call calculate_new_Tc(platform_name, &
									irw_temperature(xpoint, ypoint), &
									surface_temperature(xpoint, ypoint), &
									1.- surface_emissivity_land(xpoint, ypoint, 2), &
									idx11, &
									retrievalopticalthickness37, &
									retrievalradius37, &
									library_taus,                           &
									water_radii,                            &
									spherical_albedo_water,   &
									int_fluxdownwater_sensor,  &
									int_fluxupwater_sensor,   &
									int_cloud_emissivity_water, &
									int_surface_emissivity_water, &
									newTc, .false.)
									
				if (newTc < 0.) then 
					Tc_liquid	= curTc		 		
					exit
				endif

				call nir_absorbing_science(platform_name,&
                              optical_thickness_vector,               &
                              measurements(absorbingband_index), &
                              idx37,                &
                              albedo_real4( idx_alb37), &
							  xpoint, ypoint, & 
                              newTc,  &
                              thermal_trans_1way,                     &
                              thermal_trans_2way,                     &
                              library_taus,                           &
                              water_radii,                            &
                              spherical_albedo_water,                 &
                              int_fluxdownwater_sensor,               &
                              int_fluxdownwater_solar,                &
                              int_fluxupwater_sensor,                 &
                              int_reflectance_water,                  &
							  int_cloud_emissivity_water, &
							  int_surface_emissivity_water, &
                              residual, maxradii, channel_37, emission_uncertainty_pw_liq, &
                              emission_uncertainty_Tc_liq, sigma_R37_PW_liq,debug)
			
			
				if ( maxradii > 2 ) then!{
					call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                             water_radii(1:maxradii), &
                             retrievalradius37, &
                             retrievalopticalthickness37, &
                   debug, use_nearest, quality_in)
			if (use_nearest) then 
				nearest_liq(3) = 1			
				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											water_radii, &
											retrievalopticalthickness37, &
											retrievalradius37, RSS_liq(re37),&
											.true., ray_temp, quality_in,     &
											CH37_IDX = idx37, CTopT = cloud_top_temperature(xpoint, ypoint), &
											CH37_NUM =channel_37 , platFormName= platform_name)

			endif
				else!}{
					retrievalradius37 = fillvalue_real
					retrievalopticalthickness37 = fillvalue_real
					exit
				endif!}



				if ( abs(curTc - newTc) < 0.01 .or. retrievalradius37 < 0.) then 
						Tc_liquid = newTc
						exit
				endif

				curTc = newTc


			end do

			endif

		endif
#endif		

        retrievalradius1621 = fillvalue_real
        retrievalopticalthickness1621 = fillvalue_real
		
#ifndef SEV_PR06OD		
        if( (process%ocean_surface .or. process%snowice_surface) .and. &
            measurements(band_0163) > 0. ) then!{
           nonabsorbingband_index = band_0163
           absorbingband_index = band_0213
 
		   
		call vis_nonabsorbing_science(measurements(nonabsorbingband_index), &
                nonabsorbingband_index,                 &
                albedo_real4(nonabsorbingband_index), &
                library_taus,                           &
                water_radii,                            &
                spherical_albedo_water,                 &
                int_fluxdownwater_sensor,               &
                int_fluxdownwater_solar,                &
                int_reflectance_water,                  &
                sensor_zenith_angle(xpoint,ypoint),     &
                solar_zenith_angle(xpoint,ypoint),      &
                relative_azimuth_angle(xpoint,ypoint),  &
                cloud_top_pressure(xpoint,ypoint),      &
                local_process,                                &
                optical_thickness_vector)
                
           call vis_absorbing_science(optical_thickness_vector,          &
                measurements(absorbingband_index), &
                idx21,                    &
                albedo_real4(absorbingband_index), &
                library_taus,                           &
                water_radii,                            &
                spherical_albedo_water,                 &
                int_fluxdownwater_sensor,               &
                int_fluxdownwater_solar,                &
                int_reflectance_water,                  &
                residual,maxradii, debug)
           if (maxradii > 2) then!{
              call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                   water_radii(1:maxradii), &
                   retrievalradius1621, &
                   retrievalopticalthickness1621, &
                   debug, use_nearest, quality_in)
			if (use_nearest) then 
				nearest_liq(4) = 1

				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											water_radii, &
											retrievalopticalthickness1621, &
											retrievalradius1621, RSS_liq(re1621), &
											.true., ray_temp, quality_in)
			
			endif
           endif!}

        endif!}
#endif



!       End setting of water path for water phase clouds
		deallocate(rad37lib)
        deallocate(optical_thickness_vector, residual, reflibA, reflibB)

	optical_thickness_liquid         = retrievalopticalthickness
	optical_thickness_1621_liquid    = retrievalopticalthickness1621
	optical_thickness_16_liquid    = retrievalopticalthickness16
	optical_thickness_37_liquid    = retrievalopticalthickness37
	effective_radius_16_liquid       = retrievalradius16
	effective_radius_21_liquid       = retrievalradius21
	effective_radius_1621_liquid     = retrievalradius1621
	effective_radius_37_liquid       = retrievalradius37
	


!	if (nearest_liq(1) == 1) effective_radius_16_liquid(xpoint, ypoint) = fillvalue_real
!	if (nearest_liq(2) == 1) then 
!		effective_radius_21_liquid(xpoint, ypoint) = fillvalue_real
!		optical_thickness_liquid(xpoint, ypoint) = fillvalue_real
!	endif
!	if (nearest_liq(3) == 1) effective_radius_37_liquid(xpoint, ypoint) = fillvalue_real
!	if (nearest_liq(4) == 1) then 
!		effective_radius_1621_liquid(xpoint, ypoint) = fillvalue_real
!		optical_thickness_1621_liquid(xpoint, ypoint) = fillvalue_real
!	endif

! now reset everything so the ice cloud can be processed
		local_process%watercloud = 0
		local_process%icecloud = 1

		nonabsorbingband_index = na_band_used
        absorbingband_index = band_0163
		
        allocate(optical_thickness_vector(number_iceradii))
        allocate(residual(number_iceradii))
		allocate(reflibA(number_taus+1, number_iceradii), reflibB(number_taus+1, number_iceradii))
		allocate(rad37lib(number_taus+1, number_iceradii))
		
		rad37lib    = -999.0

! initialize all retrievals to fillvalue   
		retrievalopticalthickness = fillvalue_real
		retrievalopticalthickness1621 = fillvalue_real
		retrievalopticalthickness16 = fillvalue_real
		retrievalopticalthickness37  = fillvalue_real
		retrievalradius16 = fillvalue_real
		retrievalradius21 = fillvalue_real
		retrievalradius1621 = fillvalue_real
		retrievalradius37 = fillvalue_real


        call vis_nonabsorbing_science(measurements(nonabsorbingband_index), &
                              nonabsorbingband_index,                &
                              albedo_real4(nonabsorbingband_index), &
                              library_taus,                          &
                              ice_radii,                             &
                              spherical_albedo_ice,                  &
                              int_fluxdownice_sensor,                &
                              int_fluxdownice_solar,                 &
                              int_reflectance_ice,                   &
                              sensor_zenith_angle(xpoint,ypoint),    &
                              solar_zenith_angle(xpoint,ypoint),     &
                              relative_azimuth_angle(xpoint,ypoint), &
                              cloud_top_pressure(xpoint,ypoint),     &
                              local_process,                               &
                              optical_thickness_vector)



        if (measurements(absorbingband_index) > 0.) then!{ 
           call vis_absorbing_science(optical_thickness_vector,         &
                measurements(absorbingband_index), &
                absorbingband_index,                   &
                albedo_real4(idx_alb16), &
                library_taus,                          &
                ice_radii,                             &
                spherical_albedo_ice,                  &
                int_fluxdownice_sensor,                &
                int_fluxdownice_solar,                 &
                int_reflectance_ice,                   &
                residual,maxradii, debug)
           if ( maxradii > 2) then!{
              call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                   ice_radii(1:maxradii), &
                   retrievalradius16, &
                   retrievalopticalthickness16, &
                   debug, use_nearest, quality_in)
			if (use_nearest) then 
				nearest_ice(1) = 1
			
				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											ice_radii, &
											retrievalopticalthickness16, &
											retrievalradius16, RSS_ice(re16), &
											.false., ray_temp, quality_in)
			
			endif
           else!}{
              retrievalradius16 = fillvalue_real
              retrievalopticalthickness16 = fillvalue_real
           endif!}
           
        else!}{
           retrievalradius16 = fillvalue_real
           retrievalopticalthickness16 = fillvalue_real
        endif!}


#ifdef SEV_PR06OD
          retrievalradius21 = fillvalue_real
          retrievalopticalthickness = retrievalopticalthickness16

#else
        absorbingband_index = band_0213
        call vis_absorbing_science(optical_thickness_vector,         &
                              measurements(absorbingband_index), &
                              idx21,                   &
                              albedo_real4(absorbingband_index), &
                              library_taus,                          &
                              ice_radii,                             &
                              spherical_albedo_ice,                  &
                              int_fluxdownice_sensor,                &
                              int_fluxdownice_solar,                 &
                              int_reflectance_ice,                   &
                              residual,maxradii, debug)



        if ( maxradii > 2) then!{
        
        
          call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                             ice_radii(1:maxradii), &
                             retrievalradius21, &
                             retrievalopticalthickness, &
                   debug, use_nearest, quality_in)


			if (use_nearest) then 
				nearest_ice(2) = 1
				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											ice_radii, &
											retrievalopticalthickness, &
											retrievalradius21, RSS_ice(re21), &
											.false., alt_ray_ice, quality_in)

			
			endif
        else!}{
          retrievalradius21 = fillvalue_real
          retrievalopticalthickness = fillvalue_real
        endif!}
#endif

        absorbingband_index = band_0370

! sep, 4 May: absorbingband_index - 1 being passed in for all libarary indices because 
!             band_0370 = 7 and library index corresponding to this band is 6. In the 
!             libraries, band_0935 is index 7.
!
! wind, 7 Dec: absorbingband_index - 1 has to be passed in for albedo as well. 
!

        call nir_absorbing_science(platform_name, &
                               optical_thickness_vector,              &
                              measurements(absorbingband_index), &
                              idx37,               &
                              albedo_real4( idx_alb37), &
								xpoint, ypoint, &
                              cloud_top_temperature(xpoint, ypoint), &
                              thermal_trans_1way,                    &
                              thermal_trans_2way,                    &
                              library_taus,                          &
                              ice_radii,                             &
                              spherical_albedo_ice,                  &
                              int_fluxdownice_sensor,                &
                              int_fluxdownice_solar,                 &
                              int_fluxupice_sensor,                  &
                              int_reflectance_ice,                   &
							  int_cloud_emissivity_ice, &
							  int_surface_emissivity_ice, &
                              residual, maxradii, channel_37, emission_uncertainty_pw_ice, &
                              emission_uncertainty_Tc_ice, sigma_R37_PW_ice, debug)
		

        if ( maxradii > 2 ) then!{

          call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                             ice_radii(1:maxradii), &
                             retrievalradius37, &
                             retrievalopticalthickness37, &
                   debug, use_nearest, quality_in)
			if (use_nearest) then 
				nearest_ice(3) = 1

				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											ice_radii, &
											retrievalopticalthickness37, &
											retrievalradius37, RSS_ice(re37), &
											.false., ray_temp, quality_in,     &
											CH37_IDX = idx37, CTopT = cloud_top_temperature(xpoint, ypoint), &
											CH37_NUM =channel_37 , platFormName= platform_name)

			endif


        else!}{
          retrievalradius37 = fillvalue_real
          retrievalopticalthickness37 = fillvalue_real
        endif!}
	


											   					   
! now we iterate the retrieval
#ifndef AMS_INST
		if ( retrievalopticalthickness37 > 0. &
				.and. retrievalradius37 > 0. .and. &
				irw_temperature(xpoint, ypoint) > 0. .and. nearest_ice(3) == 0 ) then 

		! if the cloud is too thick, then don't bother doing anything, but still set the temperature
			if (retrievalopticalthickness37 > 10.) then 
				Tc_ice = cloud_top_temperature(xpoint, ypoint) 
			else


			curTc = irw_temperature(xpoint, ypoint)


		
			do i=1, 5
				
				if (retrievalopticalthickness37 < 0. .or. retrievalradius37 < 0. ) then 
					retrievalradius37 = fillvalue_real
					retrievalopticalthickness37 = fillvalue_real
					Tc_ice = curTc
					exit
				endif
				
								
				call calculate_new_Tc(platform_name, &
									irw_temperature(xpoint, ypoint), &
									surface_temperature(xpoint, ypoint), &
									1.- surface_emissivity_land(xpoint, ypoint, 2), &
									idx11, &
									retrievalopticalthickness37, &
									retrievalradius37, &
									library_taus,                           &
									ice_radii,                            &
									spherical_albedo_ice,   &
									int_fluxdownice_sensor,  &
									int_fluxupice_sensor,   &
									int_cloud_emissivity_ice, &
									int_surface_emissivity_ice, &
									newTc, .false.)
				if (newTc < 0.) then 
					Tc_ice	= curTc				
					exit
				endif
																														
        call nir_absorbing_science(platform_name, &
                               optical_thickness_vector,              &
                              measurements(absorbingband_index), &
                              idx37,               &
                              albedo_real4( idx_alb37), &
								xpoint, ypoint, &
							  newTc, &
                              thermal_trans_1way,                    &
                              thermal_trans_2way,                    &
                              library_taus,                          &
                              ice_radii,                             &
                              spherical_albedo_ice,                  &
                              int_fluxdownice_sensor,                &
                              int_fluxdownice_solar,                 &
                              int_fluxupice_sensor,                  &
                              int_reflectance_ice,                   &
							  int_cloud_emissivity_ice, &
							  int_surface_emissivity_ice, &
                              residual, maxradii, channel_37, emission_uncertainty_pw_ice, &
                              emission_uncertainty_Tc_ice, sigma_R37_PW_ice,debug)
		

        if ( maxradii > 2 ) then!{
          call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                             ice_radii(1:maxradii), &
                             retrievalradius37, &
                             retrievalopticalthickness37, &
                   debug, use_nearest, quality_in)
			if (use_nearest) then 
				nearest_ice(3) = 1
			
				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &											
											ice_radii, &
											retrievalopticalthickness37, &
											retrievalradius37, RSS_ice(re37), &
											.false., ray_temp, quality_in,     &
											CH37_IDX = idx37, CTopT = cloud_top_temperature(xpoint, ypoint), &
											CH37_NUM =channel_37 , platFormName= platform_name)

			endif



        else!}{
          retrievalradius37 = fillvalue_real
          retrievalopticalthickness37 = fillvalue_real
        endif!}


				if ( abs(curTc - newTc) < 0.01 .or. retrievalradius37 < 0.) then 
						Tc_ice = newTc
						exit
				endif


				curTc = newTc


			end do

			endif
								  
		endif
#endif		

		
        retrievalradius1621 = fillvalue_real
        retrievalopticalthickness1621 = fillvalue_real
		
#ifndef SEV_PR06OD		
        if( (process%ocean_surface .or. process%snowice_surface) .and. &
            measurements(band_0163) > 0. ) then!{
        nonabsorbingband_index = band_0163
        absorbingband_index = band_0213
		
				
        call vis_nonabsorbing_science(measurements(nonabsorbingband_index), &
                              nonabsorbingband_index,                &
                              albedo_real4(nonabsorbingband_index), &
                              library_taus,                          &
                              ice_radii,                             &
                              spherical_albedo_ice,                  &
                              int_fluxdownice_sensor,                &
                              int_fluxdownice_solar,                 &
                              int_reflectance_ice,                   &
                              sensor_zenith_angle(xpoint,ypoint),    &
                              solar_zenith_angle(xpoint,ypoint),     &
                              relative_azimuth_angle(xpoint,ypoint), &
                              cloud_top_pressure(xpoint,ypoint),     &
                              local_process,                               &
                              optical_thickness_vector)

        call vis_absorbing_science(optical_thickness_vector,         &
                              measurements(absorbingband_index), &
                              idx21,                   &
                              albedo_real4(absorbingband_index), &
                              library_taus,                          &
                              ice_radii,                             &
                              spherical_albedo_ice,                  &
                              int_fluxdownice_sensor,                &
                              int_fluxdownice_solar,                 &
                              int_reflectance_ice,                   &
                              residual,maxradii,debug)

        if (maxradii > 2 ) then!{
          call solveretrieval(residual(1:maxradii), optical_thickness_vector(1:maxradii), &
                             ice_radii(1:maxradii), &
                             retrievalradius1621, &
                             retrievalopticalthickness1621, &
                   debug, use_nearest, quality_in)
			if (use_nearest) then 
			
				nearest_ice(4) = 1
				call solveretrieval_nearest(xpoint,ypoint,measurements(nonabsorbingband_index), &
											measurements(absorbingband_index), &
											(/nonabsorbingband_index, absorbingband_index/), &
											ice_radii, &
											retrievalopticalthickness1621, &
											retrievalradius1621, RSS_ice(re1621), &
											.false., ray_temp, quality_in)
			
			endif
        endif!}
        endif!}
#endif
		deallocate(rad37lib)		
        deallocate(optical_thickness_vector, residual, reflibA, reflibB)

   else!}{
      ! there is no cloud
      retrievalopticalthickness = fillvalue_real
      retrievalopticalthickness1621 = fillvalue_real
      retrievalopticalthickness16 = fillvalue_real
      retrievalopticalthickness37  = fillvalue_real
      retrievalradius16 = fillvalue_real
      retrievalradius21 = fillvalue_real
      retrievalradius1621 = fillvalue_real
      retrievalradius37 = fillvalue_real
      cloud_layer_flag(xpoint, ypoint) = 0
      status = 1
   endif!}



! so the statistics are computed properly. The retrieval is only successful for statistical purposes if the main
! re retrieval is successful. 
   if (retrievalradius21 /= fillvalue_real) then 
      status = 1
   endif

   !assign the retrievals to the relevant arrays
	optical_thickness_ice         = retrievalopticalthickness
	optical_thickness_1621_ice    = retrievalopticalthickness1621
	optical_thickness_16_ice    = retrievalopticalthickness16
	optical_thickness_37_ice    = retrievalopticalthickness37
	effective_radius_16_ice       = retrievalradius16
	effective_radius_21_ice       = retrievalradius21
	effective_radius_1621_ice     = retrievalradius1621
	effective_radius_37_ice       = retrievalradius37


! reset the extra retrievals for now. 
!	if (nearest_ice(1) == 1) effective_radius_16_ice(xpoint, ypoint) = fillvalue_real
!	if (nearest_ice(2) == 1) then 
!		effective_radius_21_ice(xpoint, ypoint) = fillvalue_real
!		optical_thickness_ice(xpoint, ypoint) = fillvalue_real
!	endif
!	if (nearest_ice(3) == 1) effective_radius_37_ice(xpoint, ypoint) = fillvalue_real
!	if (nearest_ice(4) == 1) then 
!		effective_radius_1621_ice(xpoint, ypoint) = fillvalue_real
!		optical_thickness_1621_ice(xpoint, ypoint) = fillvalue_real
!	endif




end subroutine corescience




end module corescience_module
