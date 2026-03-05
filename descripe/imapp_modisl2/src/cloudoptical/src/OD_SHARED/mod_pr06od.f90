 subroutine mod_pr06od( statistics, errorlevel)

!-----------------------------------------------------------------------
!f90 mod_pr06od
!
!Description:
!
!  retrieve cloud optical and microphysical properties from MODIS radiation
!  measurements
!
!input parameters:
!
!output parameters:
!
!revision history:
!
! v.1 July 2001 Initial work mag gray@climate.gsfc.nasa.gov
!
!team-unique header:
!
! Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA 
!
!references and credits:
!
! Mark Gray
! gray@climate.gsfc.nasa.gov
! EmergentIT 
! Code 913, NASA GSFC 
! Greenbelt, MD 20771
!
! Nakajima, T. and M.D. King, 1990: determination of the
!    optical thickness and effective particle radius of
!    clouds from reflected solar radiation measurements
!    part i: theory. j. atmos. sci., 47, 1878-1893.
!
!design note:
!
!end
!----------------------------------------------------------------------


   use GeneralAuxType
   use core_arrays
   use nonscience_parameters
   use global_model_grids,only:grids_are_read  
   use ancillary_module
#if VIIRS_INST
   use shared_io_module_hdf5
   use modis_cloudhandler_iff
   use mod06_run_settings, only: IFF_yes
#else
   use shared_io_module
#endif
   use specific_other
   use modis_io_module
   use modis_frontend_module
   use modis_cloudhandler
   use modis_science_module
   use libraryarrays
#ifndef AMS_INST   
   use retrieval_irw
   use FASCODE_routines
#endif
   use MOD06AlbedoEcoModule
   use names, only: MY_UNIT_LUN, MY_TEXT_FILE
   use science_parameters, only: no_valid_data
   
   implicit none

!	include "hdf.f90"
!	include "dffunc.f90"

   integer, parameter :: MODFILLEN = 10
   real, intent(out)         :: statistics(7)
   integer, intent(out)      :: errorlevel

   integer, dimension (2)    :: start, edge, stride
   integer, dimension (2)    :: localstart, localedge, localstride
   integer, dimension(2) :: meas_start, meas_edge
 
   integer                   :: status, arguments, tilesize
   integer                   :: debug_start_iteration = 0, debug_stop_iteration = 0
   character(MAX_FILE_NAME_LEN)            ::  water_library, libnames_water(3), libnames_ice(3)

   character(MAX_FILE_NAME_LEN)            ::  ice_library, libnames_water_sdev(3), libnames_ice_sdev(3)

   character(MAX_FILE_NAME_LEN)            ::  phase_library
   character(MAX_FILE_NAME_LEN)            :: l1b_name(10), cloudmask_name,&
                                transmittance_library, &
                                geolocation_name, mod06_name , &
                                gdas_name, gdas_name2, ozone_name, ncepice_name, nise_name

   character(MAX_FILE_NAME_LEN)            :: surfacealbedo_lib_659,      &
                                surfacealbedo_lib_858,      &
                                surfacealbedo_lib_124,      &
                                surfacealbedo_lib_164,      &
                                surfacealbedo_lib_21,       &
                                ecosystem_data_name,        &
                                snowicealbedo_data_name, &
								emissivity_name

   real                      :: threshold_solar_zenith,      &
                                threshold_sensor_zenith,     &
                                threshold_relative_azimuth 

   integer                   :: l1b_filedata(MODFILLEN),        &
                                cloudmask_filedata(1),  &
                                geolocation_filedata(1), &
                                mod06_filedata(1), &
                                angle_filedata(1)

   logical                   :: debug

   integer, external         :: iargc
	integer :: start_iteration

	integer :: file_id, var_id, err_code, mystart(2), mystride(2), myedge(2)

!	integer*1, dimension(:,:), allocatable :: eco, snow, sfc
!	real, dimension(:,:), allocatable :: albedo

	integer :: start_time, end_time, crate, cmax, total_start, total_end
	integer :: st1, et1, cr1, cm1, anc_filedata(1)

	integer :: start_iterationX, start_iterationY


	call system_clock(total_start, crate, cmax)
	print*, "system speed: ", crate, " cycles/sec"


	IM_cloudy_count = 0
	IM_water_cloud_count = 0
	IM_ice_cloud_count = 0
	IM_undet_count = 0
	
	Statistics_1km%retrieval_fraction = 0.
	Statistics_1km%land_fraction = 0.
	Statistics_1km%water_fraction = 0.
	Statistics_1km%snow_fraction = 0.
	Statistics_1km%cloud_fraction = 0.
	Statistics_1km%water_cloud_fraction = 0.
	Statistics_1km%ice_cloud_fraction = 0.
	Statistics_1km%mean_liquid_tau = 0.
	Statistics_1km%mean_ice_tau = 0.
	Statistics_1km%mean_liquid_re21 = 0.
	Statistics_1km%mean_ice_re21 = 0.
	Statistics_1km%ctp_liquid = 0.
	Statistics_1km%ctp_ice = 0.
	Statistics_1km%ctt_ice = 0.
	Statistics_1km%ctt_liquid = 0.
	Statistics_1km%ctp_undetermined = 0.
	Statistics_1km%ctt_undetermined = 0.
	



   status = success
   errorlevel = 0
   statistics(:) = .0

!  get runtime setup parameters, array sizes, lat/long-pixel limits
   if (status == success) then
      call initialize_run (   start, edge, stride,       &
                              tilesize,                  &
                              l1b_name,                   &
                              cloudmask_name,            &
                              mod06_name,                &
                              geolocation_name,          &
                              gdas_name,                 &
                              gdas_name2,                 &
							  ozone_name, &
                              ncepice_name,              &
                              nise_name,                 &
							  water_library,              &
							  libnames_water, &
							  libnames_water_sdev, &
                              ice_library,                &
							  libnames_ice, &
							  libnames_ice_sdev, &
							  phase_library,			&
                              transmittance_library,     &
                              surfacealbedo_lib_659,      &
                              surfacealbedo_lib_858,      &
                              surfacealbedo_lib_124,      &
                              surfacealbedo_lib_164,      &
                              surfacealbedo_lib_21,       &
                              ecosystem_data_name,        &
                              snowicealbedo_data_name,    &
							  emissivity_name, &
                              threshold_solar_zenith,    &
                              threshold_sensor_zenith,   &
                              threshold_relative_azimuth,&
                              status )
     if (status /= success) then
       call local_message_handler('Problem reported in initialize_run. Check previous message/s', status,'mod_pr06od')
     endif
   endif


!  check that all files exist and are openable
!  if not there is a complete algorithm failure
   if( status  == success) then
      call check_datasources(l1b_name,                   &
                             cloudmask_name,            &
                             geolocation_name,          &
                             gdas_name,                 &
                             gdas_name2,                 &
							 ozone_name, &
                             ncepice_name,              &
                             nise_name,                 &
                             mod06_name,           &
							  water_library,              &
							  libnames_water, &
							  libnames_water_sdev, &
                              ice_library,                &
							  libnames_ice, &
							  libnames_ice_sdev, &
							 phase_library,			& 
                             transmittance_library,     &
                             surfacealbedo_lib_659,     &
                             surfacealbedo_lib_858,     &
                             surfacealbedo_lib_124,     &
                             surfacealbedo_lib_164,     &
                             surfacealbedo_lib_21,      &
                             ecosystem_data_name,       &
                             snowicealbedo_data_name,   & 
							 emissivity_name, &
                             status)
     if (status /= success) then
       call local_message_handler('Problem reported in check_datasources. Check previous message/s', status,'mod_pr06od')
     endif
   endif
   

!  open core hdf files for reading and chec 

   if (status == success) then
     call openclose_files (l1b_name,            &
                           cloudmask_name,     &
                           geolocation_name,   &
                           mod06_name,         &
                           'open',             &
                           l1b_filedata,       &
                           cloudmask_filedata, &
                           geolocation_filedata,&  
                           mod06_filedata,     &
                           status)
                           
     if (status /= success) then
       call local_message_handler('Problem reported in openclose_files . Check previous message/s', status,'mod_pr06od')
     endif
   endif

!  check the size of the array/s to be processed 
   if (status == success) then
     call check_datasize(l1b_filedata, start, stride, edge, status)
     if (status /= success) then
       call local_message_handler('Problem reported in check_datasize . Check previous message/s', status,'mod_pr06od')
     endif
   endif
   
!  if the total number of lines to read is greater than the 
!  tilesize (number of lines to process/ iteration) then
!  we need to iterate.  

!     read libraries in preparation for later interpolation
      if (status == success ) then
         call readlibraries_base(transmittance_library,	water_library,              &
                              ice_library, & 
							  libnames_ice(1),     &
							  phase_library,debug,status)
         if (status /= success) then
           call local_message_handler('Problem reported in readlibraries',status,'mod_pr06od')
         endif
      endif

#ifndef AMS_INST
	call init_FASCODE
	call init_irw
#endif

! for inventory metadata
   total_number_of_pixels = edge(2)*edge(1)

	mystart = 0
	mystride = 1
	myedge(1) = edge(1)
	myedge(2) = edge(2)
	
!	allocate(eco(edge(1), edge(2)), snow(edge(1), edge(2)), sfc(edge(1), edge(2)), albedo(edge(1), edge(2)))
	
	
#ifdef SEV_PR06OD
! because of SEVIRI's angle space, it's not possible to retrieve a single data stripe
! we have to use square blocks in order to keep the memory from going out of control 
   number_of_iterationsX = ceiling(real(edge(1)) / real(tilesize)) 
   number_of_iterationsY = ceiling(real(edge(2)) / real(tilesize))
#else
   number_of_iterationsX = ceiling(real(edge(1)) / real(tilesize))
   number_of_iterationsY = 1
#endif


   grids_are_read = .false.
   NISE_is_read = .false.
   snow_stats_are_read = .false.

   ! Default behavior for setting debug <start/stop> _iteration is to pick only the (10 scan line)
   ! iteration that bounds debugPxlYcoord_1to20L0inclusiv.  Specify additional iterations here. 
   if(debugPxlYcoord_1to20L0inclusiv .ne. 0) then
      debug_start_iteration = floor((real(debugPxlXcoord_1to1354inclusiv-1)/real(tilesize))) + 1
      debug_stop_iteration = debug_start_iteration
   endif
   debugPRN = .false.


#ifdef VIIRS_INST
	if (IFF_yes) then 
		anc_filedata(1) = geolocation_filedata(1)
	else
		anc_filedata(1) = cloudmask_filedata(1)
	endif
#else
	anc_filedata(1) = mod06_filedata(1)
#endif

	print*, "NUM_ITER: X, Y, total: ", number_of_iterationsX, number_of_iterationsY, &
											number_of_iterationsX * number_of_iterationsY


	start_iterationX = 1
	start_iterationY = 1

!	open(unit=MY_UNIT_LUN, file = MY_TEXT_FILE)
   do iterationX = start_iterationX, number_of_iterationsX

	do iterationY = start_iterationY, number_of_iterationsY


      print*, "ITERATION X Y: ", iterationX, iterationY

	call system_clock(start_time, crate, cmax)
		

      localstart  = start
      localedge   = edge
      localstride = stride

	  meas_start = start
	  meas_edge = edge


! we read one line more than we need on each side so that CSR can work correctly. 
      if (edge(1) > tilesize .and. iterationX == 1) then
         localedge(1) = tilesize
		 meas_edge(1) = tilesize+1
      endif
      if (edge(1) > tilesize .and. iterationX > 1 ) then 
         localstart(1)  = start(1) + (iterationX-1)*tilesize*stride(1)
		 meas_start(1) = start(1) + (iterationX-1)*tilesize*stride(1)-1
         localedge(1)   = tilesize
		 meas_edge(1) = tilesize+2
      endif


	
      if ( edge(1) - iterationX*tilesize <= 0 ) then 
		localedge(1) = edge(1) - (iterationX-1)*tilesize
		meas_edge(1) = edge(1) - (iterationX-1)*tilesize+1
	  endif

#ifdef SEV_PR06OD 
      if (edge(2) > tilesize .and. iterationY == 1) then
         localedge(2) = tilesize
      endif
      if (edge(2) > tilesize .and. iterationY > 1 ) then 
         localstart(2)  = start(2) + (iterationY-1)*tilesize*stride(2)
		 localedge(2)   = tilesize
      endif


      if ( edge(2) - iterationY*tilesize < 0 ) then 
		localedge(2) = edge(2) - (iterationY-1)*tilesize
	  endif
	  
#endif

#if MAS_OD || SEV_PR06OD || MODIS_HKM || AMS_OD
		meas_start = localstart
		meas_edge = localedge
#endif

	print*, "Chunk coords X, X+wid, Y, Y+ht:", localstart(1), localstart(1) + localedge(1), &
			localstart(2), localstart(2) + localedge(2), meas_edge(1), meas_edge(2)
	
			
!     Setup and allocate all data arrays
      if (status == success ) then
         call allocate_arrays (localedge, &
							   meas_edge, &
							   start_iterationX, &
							   start_iterationY, &
                               status )
         if (status /= success) then
           call local_message_handler('Failure detected before allocate_modisarrays, check previous messages', &
                                      status, &
                                      'mod_pr06od')
         endif
      endif
    
      call init_qualitydata
	  
	  print*, "allocate: "
	  
      ! get data cube to be processed and ancillary data arrays
      if (status == success ) then
         call get_modis_data_cube (l1b_filedata, &
                                   geolocation_filedata, &
                                   localstart,  &
                                   localedge,   &
                                   localstride, &
								   meas_start, &
								   meas_edge, &
                                   iterationX,   &
                                   debug,       &
                                   status )
         if (status /= success) then
           call local_message_handler('Failure detected before get_modis_data_cube check previous messages',&
                                      status, &
                                      'mod_pr06od')
         endif
      endif
		print*, "read data"


! if it's an empty iteration, then don't even bother wasting time doing anything else
	if (no_valid_data == 1) then 
		  
         call output_retrieval(mod06_filedata, &
                               mod06_name,     &
							   iterationX, iterationY, &
							   number_of_iterationsX, number_of_iterationsY, &
                               localstart, localedge, localstride, &
                               status)  


		 cycle
	
	endif

!     get cloud decision info
      if (status == success ) then
      
      
#ifdef VIIRS_INST
		if (IFF_yes) then 
		
	        call modis_cloudprocessing_iff(cloudmask_filedata, mod06_filedata, &
                                    iterationX,   &
                                    localstart,&
                                    localstride,&
                                    localedge,&
                                    debug,&
                                    status)		
		else
		
	        call modis_cloudprocessing(cloudmask_filedata, mod06_filedata, &
                                    iterationX,   &
                                    localstart,&
                                    localstride,&
                                    localedge,&
                                    debug,&
                                    status)		
		endif
#else      
         call modis_cloudprocessing(cloudmask_filedata, mod06_filedata, &
                                    iterationX,   &
                                    localstart,&
                                    localstride,&
                                    localedge,&
                                    debug,&
                                    status)
#endif



         if (status /= success) then
           call local_message_handler('Failure detected before cloudprocessing. Check previous messages/s', &
                                      status, &
                                      'mod_pr06od')
         endif
      endif

!     get ancillary data, right now that means albedos
      if (status == success ) then 
         call set_ancillary(gdas_name, gdas_name2, ozone_name,  &
                            ncepice_name,              &
                            nise_name,                 &
                            mod06_name,           &
							anc_filedata, &
                            surfacealbedo_lib_659,     &
                            surfacealbedo_lib_858,     &
                            surfacealbedo_lib_124,     &
                            surfacealbedo_lib_164,     &
                            surfacealbedo_lib_21,      &
                            ecosystem_data_name,       &
                            snowicealbedo_data_name,   &
							emissivity_name, &
                            localstart, localstride, localedge, &
                            debug,                     &
                            status )
         if (status /= success) then
           call local_message_handler('Problem reported in set_ancillary. Check previous message/s', &
                                      status, &
                                      'mod_pr06od')
         endif
      endif
	print*, "ancillary set"

      if (status == 0 ) then
         call readlibraries_extra(water_library,              &
							libnames_water, &
							libnames_water_sdev, &
							ice_library,                &
 							libnames_ice, &
							libnames_ice_sdev, &
								debug,status)
								
         if (status /= success) then
           call local_message_handler('Problem reported in readlibraries',status,'mod_pr06od')
         endif
      endif

!     process modis data for cloud optical thickness and cloud top effective
!     particle radius
	  
	print*, "libraries read"
      if (status == success .and. ( (debug_start_iteration .eq. 0) .or. &
          ( debug_start_iteration .le. iterationX .and. iterationX .le. debug_stop_iteration ) ) ) then

         call scienceinterface( threshold_solar_zenith,      &
                                threshold_sensor_zenith,     &
                                threshold_relative_azimuth,  &
                                debug,                       &
                                status)

         if (status /= success) then
           call local_message_handler('Problem reported in scienceinterface',status,'mod_pr06od')
         endif
      endif
	print*, "done retrieval"
	
!   eco(localstart(1)+1:localstart(1) + tilesize , &
!             localstart(2)+1:localstart(2) + localedge(2) ) = eco_out(:,:)

!     send output retrieval and qa data to a hdf file
      if (status == success ) then
         call output_retrieval(mod06_filedata, &
                               mod06_name,     &
							   iterationX, iterationY, &
							   number_of_iterationsX, number_of_iterationsY, &
                               localstart, localedge, localstride, &
                               status)  
         if (status /= success) then
           call local_message_handler('Problem reported in output_retrieval',status,'mod_pr06od')
         endif
      endif

	  call system_clock(end_time, crate, cmax)
	  
	  print*, "iteration_time: ", end_time - start_time
#ifdef MODIS_INST
	  call aggregate_statistics_1km
#endif
	  
	 enddo 
   enddo

      if (status == success) then
         call deallocate_cleanup(status )
         if (status /= success) then
           call local_message_handler('Problem reported in deallocate_cleanup', status,'mod_pr06od')
         endif
      endif
   

   call set_processing_attributes(mod06_filedata) 
   
   

   if (status == success) then
     call openclose_files (l1b_name,            &
                           cloudmask_name,     &
                           geolocation_name,   &    
                           mod06_name,         &
                           'close',            &
                           l1b_filedata,       &
                           cloudmask_filedata, &
                           geolocation_filedata,&
                           mod06_filedata,     &
                           status)
     if (status /= success) then
       call local_message_handler('Problem reported in openclose_files', status,'mod_pr06od')
     endif
   endif

! here we do the final assignment of the Inventory metadata statistics

#ifdef  MODIS_INST

	if (total_number_of_pixels > 0) then 
	   statistics(2) = real(IM_cloudy_count) / real(total_number_of_pixels) * 100.
	   Statistics_1km%land_fraction = Statistics_1km%land_fraction / real(total_number_of_pixels) * 100.
	   Statistics_1km%water_fraction = Statistics_1km%water_fraction / real(total_number_of_pixels) * 100.
	   Statistics_1km%snow_fraction = Statistics_1km%snow_fraction / real(total_number_of_pixels) * 100.
	else
	   statistics(2) = -9999.
	   Statistics_1km%land_fraction = -9999.
	   Statistics_1km%water_fraction = -9999.
	   Statistics_1km%snow_fraction = -9999.
	endif

	if (IM_cloudy_count > 0) then 
	   statistics(1) = real(IM_successful_retrieval_count) / real(IM_cloudy_count) * 100.
	   statistics(3) = real(IM_water_cloud_count) / real(IM_cloudy_count) * 100.
	   statistics(4) = real(IM_ice_cloud_count) / real(IM_cloudy_count) * 100.
	else
	   statistics(1) = -9999.
	   statistics(3:4) = -9999.
	endif

   
    Statistics_1km%retrieval_fraction = statistics(1)
	Statistics_1km%cloud_fraction = statistics(2)
	Statistics_1km%water_cloud_fraction = statistics(3)	
	Statistics_1km%ice_cloud_fraction = statistics(4)

	if (IM_water_cloud_count > 0) then 
		Statistics_1km%mean_liquid_tau = Statistics_1km%mean_liquid_tau / real(IM_water_cloud_count) 
		Statistics_1km%mean_liquid_re21 = Statistics_1km%mean_liquid_re21 / real(IM_water_cloud_count) 
		Statistics_1km%ctp_liquid = Statistics_1km%ctp_liquid / real(IM_water_cloud_count) 
		Statistics_1km%ctt_liquid = Statistics_1km%ctt_liquid / real(IM_water_cloud_count) 
	else
		Statistics_1km%mean_liquid_tau = -9999.
		Statistics_1km%mean_liquid_re21 = -9999.
		Statistics_1km%ctp_liquid = -9999.
		Statistics_1km%ctt_liquid = -9999.
	endif

	if (IM_ice_cloud_count > 0) then 
		Statistics_1km%mean_ice_tau = Statistics_1km%mean_ice_tau / real(IM_ice_cloud_count) 
		Statistics_1km%mean_ice_re21 = Statistics_1km%mean_ice_re21 / real(IM_ice_cloud_count) 
		Statistics_1km%ctp_ice = Statistics_1km%ctp_ice / real(IM_ice_cloud_count) 
		Statistics_1km%ctt_ice = Statistics_1km%ctt_ice / real(IM_ice_cloud_count) 
	else
		Statistics_1km%mean_ice_tau = -9999.
		Statistics_1km%mean_ice_re21 = -9999.
		Statistics_1km%ctp_ice = -9999.
		Statistics_1km%ctt_ice = -9999.
	endif
   
	if (IM_undet_count > 0) then 
		Statistics_1km%ctp_undetermined = Statistics_1km%ctp_undetermined / real(IM_undet_count) 
		Statistics_1km%ctt_undetermined = Statistics_1km%ctt_undetermined / real(IM_undet_count) 
    else
    	Statistics_1km%ctp_undetermined = -9999.
    	Statistics_1km%ctt_undetermined = -9999.
    endif
    

	call write_statistics(mod06_name)

   
#endif
   
!    close (MY_UNIT_LUN)
   
   	call system_clock(total_end, crate, cmax)
	print*, "Total Execution Time is: ", total_end - total_start, " cycles"
    print*, "total execution time is: ", (total_end - total_start) / ( crate*1.0 ) / 60.0, " minutes"
   

 end subroutine mod_pr06od
