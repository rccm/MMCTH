module modis_frontend_module

 implicit none
 
 include "hdf.f90"
 include "dffunc.f90"
 
 private

 integer                                    :: index_solarzenith, index_sensorzenith,&
                                                index_relativeazimuth,index_solarzenith_flux

 real, dimension(:), allocatable    :: solarzenith_all, sensorzenith_all, &
                                            relativeazimuth_all, solarzenith_flux_all


 public :: initialize_run, check_datasources, readlibraries_base, readlibraries_extra, &
          allocate_arrays, init_qualitydata, deallocate_cleanup



 contains

 subroutine initialize_run( start, edge, stride,        &
                           tilesize,                   &
                           level1b_name,               &
                           cloudmask_name,             &
                           mod06_name,            &
                           geolocation_name,           &
                           gdas_name,                  &
                           gdas_name2,                  &
						   ozone_name, &
                           ncepice_name,               &
                           nise_name,                  &
                           water_library,              &
							libnames_water, &
							libnames_water_sdev, &
                           ice_library,                &
							libnames_ice, &
							libnames_ice_sdev, &
						   phase_library,				&
                           transmittance_library,      &
                           surfacealbedo_lib_659,      &
                           surfacealbedo_lib_858,      &
                           surfacealbedo_lib_124,      &
                           surfacealbedo_lib_164,      &
                           surfacealbedo_lib_21,       &
                           ecosystem_data_name,        &
                           snowicealbedo_data_name,    &
						   emissivity_name, &
                           threshold_solar_zenith,     &
                           threshold_sensor_zenith,    &
                           threshold_relative_azimuth, &
                           status )
   
   use core_arrays
   use mod06_run_settings
   use nonscience_parameters
   use names
   use specific_ancillary
   use specific_other

   implicit none


   integer,          intent (out)         :: status, tilesize 

	character(*), intent(out) :: libnames_water(3), libnames_ice(3), &
								 libnames_water_sdev(3), libnames_ice_sdev(3)
   character(*),     intent (out)         :: level1b_name(:),  cloudmask_name,&
                                             water_library, ice_library, &
											 phase_library, &
                                             transmittance_library, &
                                             geolocation_name, mod06_name, &
                                             gdas_name, gdas_name2, ozone_name, ncepice_name,  &
                                             nise_name, surfacealbedo_lib_659,      &
                                             surfacealbedo_lib_858,      &
                                             surfacealbedo_lib_124,      &
                                             surfacealbedo_lib_164,      &
                                             surfacealbedo_lib_21,       &
                                             ecosystem_data_name,        &         
                                             snowicealbedo_data_name, emissivity_name


   real,                     intent(out)  :: threshold_solar_zenith,      &
                                             threshold_sensor_zenith,     &
                                             threshold_relative_azimuth

   integer,  dimension (2),  intent (out) :: start,edge,stride 

   integer                                :: number_of_bands, checkvariable,lun, readstatus,&
                                             i,index, file_version, localstatus 

   character*(20) , parameter :: functionname = 'initialize_run'
   integer fileid, version, gpn_status


!  function definitions
   integer get_platform_name
   external get_platform_name

	status = 0

 
	call set_l1b_names(level1b_name)

	file_version = 1
	gpn_status = 0

#ifdef VIIRS_INST
	localstatus = get_platform_name(Acloudmask_name, platform_name)
#else
	localstatus = get_platform_name(level1b_name(1), platform_name)
#endif

	call get_channels
 
 
	cloudmask_name = Acloudmask_name
   
	water_library = Awater_library

	libnames_water = Alibnames_water
	libnames_ice = Alibnames_ice
	libnames_water_sdev = Alibnames_water_sdev
	libnames_ice_sdev = Alibnames_ice_sdev
	

	ice_library = Aice_library
	
	phase_library = Aphase_library

	transmittance_library = Atransmittance_library
		
	mod06_name = Amod06_name


	geolocation_name = Ageolocation_name

	gdas_name = Agdas_name
	gdas_name2 = Agdas_name2

	ozone_name = Aozone_name
	ncepice_name = Ancepice_name

	nise_name = Anise_name

	surfacealbedo_lib_659 = Asurfacealbedo_lib_659
	surfacealbedo_lib_858 = Asurfacealbedo_lib_858
	surfacealbedo_lib_124 = Asurfacealbedo_lib_124
	surfacealbedo_lib_164 = Asurfacealbedo_lib_164
	surfacealbedo_lib_21 = Asurfacealbedo_lib_21


	ecosystem_data_name = Aecosystem_data_name

	snowicealbedo_data_name = Asnowicealbedo_data_name

	emissivity_name = Aemissivity_name

	call get_data_dims(mod06_name, start, stride, edge)


   threshold_solar_zenith = set_threshold_solar_zenith
   threshold_sensor_zenith = set_threshold_sensor_zenith
   threshold_relative_azimuth = set_threshold_relative_azimuth

   number_of_bands = size(set_bands)
   allocate(bands(number_of_bands))
   bands(:) = set_bands(:)
  
   tilesize = set_tilesize


 end subroutine initialize_run

 subroutine check_datasources(l1bname,                   &
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
						  phase_library, &
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
    use nonscience_parameters
    implicit none

	character(*), intent(in) :: libnames_water(3), libnames_ice(3), &
								 libnames_water_sdev(3), libnames_ice_sdev(3)

	character(len=MAX_FILE_NAME_LEN), intent(in) :: l1bname(:)
    character(*), intent(in)  :: cloudmask_name, geolocation_name, gdas_name, gdas_name2, ncepice_name,  &
                                  nise_name, mod06_name, water_library, ice_library, phase_library, &
                                 transmittance_library,  surfacealbedo_lib_659, ozone_name, &
                                 surfacealbedo_lib_858, surfacealbedo_lib_124, surfacealbedo_lib_164, &
                                 surfacealbedo_lib_21, ecosystem_data_name, snowicealbedo_data_name, emissivity_name
    
    integer,dimension(35)      :: allstatus
    integer, intent(out)       :: status
	integer :: i

    status = 0
    allstatus = 0


    allstatus(1) = checkfile(l1bname(1))
    if (allstatus(1) /= success) then
      status = failure
      call local_message_handler ('Error noted: level 1B file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(2) = checkfile(cloudmask_name)
    if (allstatus(2) /= success) then
      status = failure
      call local_message_handler ('Error noted: cloud mask file. Check PCF and input files.', &
                                   status, &
                                   'check_datasources')
    endif
    allstatus(3) = checkfile(geolocation_name)
    if (allstatus(3) /= success) then
      status = failure
      call local_message_handler ('Error noted: geolocation file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(4) = checkfile(gdas_name)    
    if (allstatus(4) /= success) then
      status = failure
      call local_message_handler ('Error noted: GDAS file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(33) = checkfile(gdas_name2)    
    if (allstatus(33) /= success) then
      status = failure
      call local_message_handler ('Error noted: GDAS file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif


    allstatus(5) = checkfile(ncepice_name)    
    if (allstatus(5) /= success) then
      status = failure
      call local_message_handler ('Error noted: NCEP Ice file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif

#ifdef USE_TOAST	
    allstatus(6) = checkfile(ozone_name)    
    if (allstatus(6) /= success) then
      status = failure
      call local_message_handler ('Error noted: NCEP Ice file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
#else
	allstatus(6) = success
#endif
	
    allstatus(7) = checkfile(nise_name)    
    if (allstatus(7) /= success) then
      status = failure
      call local_message_handler ('Error noted: NISE  file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif

	allstatus(8) = success

    allstatus(9) = checkfile(mod06_name)    
    if (allstatus(9) /= success) then
      status = failure
      call local_message_handler ('Error noted: MOD06_L2 file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(10) = checkfile(water_library)    
    if (allstatus(10) /= success) then
      status = failure
      call local_message_handler ('Error noted: water library file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(11) = checkfile(ice_library)    
    if (allstatus(11) /= success) then
      status = failure
      call local_message_handler ('Error noted: ice library file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(12) = checkfile(transmittance_library)    
    if (allstatus(12) /= success) then
      status = failure
      call local_message_handler ('Error noted: transmittance file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(13) = checkfile(surfacealbedo_lib_659)     
    if (allstatus(13) /= success) then
      status = failure
      call local_message_handler ('Error noted: albedo (.659um) file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(14) = checkfile(surfacealbedo_lib_858)    
    if (allstatus(14) /= success) then
      status = failure
      call local_message_handler ('Error noted: albedo (.858um) file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(15) = checkfile(surfacealbedo_lib_124)    
    if (allstatus(15) /= success) then
      status = failure
      call local_message_handler ('Error noted: albedo (1.24um) file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(16) = checkfile(surfacealbedo_lib_164)    
    if (allstatus(16) /= success) then
      status = failure
      call local_message_handler ('Error noted: albedo (1.64um) file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(17) = checkfile(surfacealbedo_lib_21)    
    if (allstatus(17) /= success) then
      status = failure
      call local_message_handler ('Error noted: albedo (2.1um) file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(18) = checkfile(ecosystem_data_name)    
    if (allstatus(18) /= success) then
      status = failure
      call local_message_handler ('Error noted: IGBP Ecosystem file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(19) = checkfile(snowicealbedo_data_name)    
    if (allstatus(19) /= success) then
      status = failure
      call local_message_handler ('Error noted: Static SnowIce Albedo file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
    allstatus(20) = checkfile(phase_library)    
    if (allstatus(20) /= success) then
      status = failure
      call local_message_handler ('Error noted: phase library file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
	
	do i=1, 3
		allstatus(20+i) = checkfile(libnames_water(i))    
		if (allstatus(20+i) /= success) then
			status = failure
			call local_message_handler ('Error noted: water library file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
		endif
		allstatus(23+i) = checkfile(libnames_water_sdev(i))    
		if (allstatus(23+i) /= success) then
			status = failure
			call local_message_handler ('Error noted: water library file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
		endif
		allstatus(26+i) = checkfile(libnames_ice(i))    
		if (allstatus(26+i) /= success) then
			status = failure
			call local_message_handler ('Error noted: water library file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
		endif
		allstatus(29+i) = checkfile(libnames_ice_sdev(i))    
		if (allstatus(29+i) /= success) then
			status = failure
			call local_message_handler ('Error noted: water library file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
		endif

	end do

    allstatus(34) = checkfile(emissivity_name)    
    if (allstatus(34) /= success) then
      status = failure
      call local_message_handler ('Error noted: phase library file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
	
	
	if (l1bname(2) /= "none" ) then 
		allstatus(35) = checkfile(l1bname(2))
		if (allstatus(1) /= success) then
			status = failure
			call local_message_handler ('Error noted: level 1B file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
		endif
	else
		allstatus(35) = 0
	endif
	
    status = maxval(allstatus)
   
    if (status /= success) then
      status = failure 
      call local_message_handler ('Error noted for a required input file. Check PCF and input files.' ,&
                                   status, &
                                   'check_datasources')
    endif
 end subroutine check_datasources
 integer function checkfile(name)

  implicit none
  
  character(*), intent(in)  :: name
  logical                   :: exist
  character(len=11)         :: readability

  checkfile = 0

  inquire( file = name,      &
           exist= exist,   &
           read = readability)
  if (.not. exist) then
         checkfile = 1
  endif
  if (readability == 'NO') then
         checkfile = 2
  endif
 end function checkfile



! this subroutine returns the needed array bounds for a resized array
 subroutine find_bounds(array, min_val, max_val, min_bnd, max_bnd) 
 
	implicit none
	
	real, dimension(:), intent(in) :: array
	real, intent(in) :: min_val, max_val
	integer, intent(inout) :: min_bnd, max_bnd
	
	integer :: i, N
	
	N = size(array)
	
	do i=1, N
		if (array(i) >= min_val) then 
			min_bnd = i-1
			exit
		endif
	end do
 
	if (min_bnd < 1) min_bnd = 1
	if (min_bnd > N) min_bnd = N
 
	do i=1, N
		if (array(i) >= max_val) then 
			max_bnd = i
			exit
		endif
	end do

	if (max_bnd >= N) max_bnd = N
	if (max_bnd <= 1) max_bnd = 1

 end subroutine find_bounds
 
!-----------------------------------------------------------------------
!f90 readlibraries
!
!Description:
!
! Read bidirectional reflectance and flux library values for all
! angles, as defined in a library description file.  Library arrays
! are allocated as required to read all library data.
!
!input parameters:
!
!output parameters:
!
!revision history:
!
! v2.1  July 2002 arrays are now allocated to only include the angles
!               needed.
!
! v2.0  June 2002 Changed the read routines to work with the HDF
!         version of the libraries. -- Gala. wind@climate.gsfc.nasa.gov
!
! v1.0  November 2001 Initial work mag gray@climate.gsfc.nasa.gov
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
! Gala Wind
! wind@climate.gsfc.nasa.gov
! L-3 Comm Analytics
! all same, all same
!
!design note:
!
!end
!----------------------------------------------------------------------

 subroutine readlibraries_base(trans_correct_lib, hdf_water_lib, hdf_ice_lib, &
							hdf_ice_lib_ws3, &
							hdf_phase_lib, debug,status)
   
   use GeneralAuxType
   use libraryarrays
   use libraryinterpolates
   use science_parameters
   use nonscience_parameters
   use hdf_mod
	use interpolate_libraries   
   implicit none  

   character(*), intent(in)                   :: hdf_water_lib, hdf_ice_lib,hdf_phase_lib, trans_correct_lib
   character(*), intent(in)                   :: hdf_ice_lib_ws3

   logical, intent(in)                        :: debug
   integer,intent (out)                       :: status

   integer                                    :: file_id, var_id, var_index,i, j
   integer									:: start(6), stride(6), edge(6)
   character(MAX_SDS_NAME_LEN)                              :: dummy_name
   integer                                    :: dummy_type, dummy_num_attr, dummy_rank
   integer                                    :: dim_sizes(6),dim_sizes4d(4)




   status = 0



!  read transmittance correction library

   start = 0
   stride = 1

   file_id = sfstart(trans_correct_lib, DFACC_READ)
   var_id = sfselect(file_id, sfn2index(file_id, "Transmittance"))
   status = sfginfo(var_id, dummy_name, dummy_rank, &
                      dim_sizes, dummy_type, dummy_num_attr)

	edge(1:4) = dim_sizes(1:4)

	allocate(transmit_correction_table(edge(1), edge(2), edge(3), edge(4)))
	allocate(transmit_stddev_table(edge(1), edge(2), edge(3), edge(4)))

   status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), transmit_correction_table)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "Standard_Deviation"))
   status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), transmit_stddev_table)
   status = sfendacc(var_id)

   status = sfend(file_id)



! we need to find angle range for the granule, so we can only read a limited amount from the library
! just enough to get the granule processed. 

! from the find_granule_angle_range routine we know just how wide the angle space we need. 
! it's fixed for sensor zenith and relative azimuth, but not for solar angle. However for MAS<TER>
! the relative azimuth space isn't fixed either and the sensor zenith limits are different. 


	number_wind_speed = 3

! start with the ice library first
	file_id = sfstart(hdf_ice_lib, DFACC_READ)

! get the number of wavelengths, we never use them, so we don't need to actually read them. 
	var_id = sfselect(file_id, sfn2index(file_id, "Wavelengths"))
	status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)
	number_wavelengths = dim_sizes(1)
	status = sfendacc(var_id)

	start = 0
	stride = 1

! get the radii
    var_id = sfselect(file_id, sfn2index(file_id, "ParticleRadius"))
    status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)
    number_iceradii = dim_sizes(1)
    allocate (ice_radii(number_iceradii))
    edge(1) = number_iceradii
    status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), ice_radii)
    status = sfendacc(var_id)

  ! get the taus
     var_id = sfselect(file_id, sfn2index(file_id, "OpticalThickness"))
     status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)
     number_taus = dim_sizes(1)
     allocate (library_taus(number_taus))
	 edge(1) = number_taus
	 status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), library_taus)
     status = sfendacc(var_id)

! get the full 1D angle arrays, they are small, we'll hang on to them and do resizing after we read them
! in the read_libraries_extra() subroutine

! now get the relative azimuth
    var_id = sfselect(file_id, sfn2index(file_id, "ReflectanceRelativeAzimuth"))
    status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)

	index_relativeazimuth = dim_sizes(1)
	allocate (relativeazimuth_all(index_relativeazimuth))
    edge(1) = index_relativeazimuth
    status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), relativeazimuth_all)
    status = sfendacc(var_id)

! now get the solar zenith 
    var_id = sfselect(file_id, sfn2index(file_id, "ReflectanceSolarZenith"))
    status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)

	index_solarzenith = dim_sizes(1)
	allocate (solarzenith_all(index_solarzenith))
    edge(1) = index_solarzenith
    status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), solarzenith_all)
    status = sfendacc(var_id)
	
! now get the sensor zenith 
    var_id = sfselect(file_id, sfn2index(file_id, "ReflectanceSensorZenith"))
    status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)

	index_sensorzenith = dim_sizes(1)
	allocate (sensorzenith_all(index_sensorzenith))
    edge(1) = index_sensorzenith
    status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), sensorzenith_all)
    status = sfendacc(var_id)

! now get the flux angle 
    var_id = sfselect(file_id, sfn2index(file_id, "FluxSolarZenith"))
    status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)

	index_solarzenith_flux = dim_sizes(1)
	allocate (solarzenith_flux_all(index_solarzenith))
    edge(1) = index_solarzenith_flux
    status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), solarzenith_flux_all)
    status = sfendacc(var_id)

! Read the 2D arrays 


   start = 0
   stride = 1
   edge(1) = number_wavelengths
   edge(2) = number_iceradii

   allocate(extinction_ice(number_wavelengths,number_iceradii))
   allocate(singlescattering_ice(number_wavelengths,number_iceradii))
   allocate(asymmetry_ice(number_wavelengths,number_iceradii))
   allocate(truncation_factor_ice(number_wavelengths,number_iceradii))
   allocate(phase_fun_norm_constant_ice(number_wavelengths,number_iceradii))

   var_id = sfselect(file_id, sfn2index(file_id, "ExtinctionCoefficient"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), extinction_ice)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "SingleScatterAlbedo"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), singlescattering_ice)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "Phase: AsymmetryParameter"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), asymmetry_ice)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "TruncationFactor"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), truncation_factor_ice)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "PhaseFuncNormConstant"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), phase_fun_norm_constant_ice)
   status = sfendacc(var_id)


! Read the 3D arrays

	
   start = 0
   stride = 1
   edge(1) = number_taus
   edge(2) = number_wavelengths
   edge(3) = number_iceradii
   
   allocate(spherical_albedo_ice(number_taus,number_wavelengths,number_iceradii))
   var_id = sfselect(file_id, sfn2index(file_id, "SphericalAlbedo"))
   status = sfrdata(var_id, start(1:3),stride(1:3),edge(1:3), spherical_albedo_ice)
   status = sfendacc(var_id)


! Read the 4D arrays
! Now we finally read the 6D reflectance array

 
   status = sfend(file_id)


!  Now that we're done with ice, get started with water...

!  Read all the water stuff there is to read.
   file_id = sfstart(hdf_water_lib, DFACC_READ)


! get the radii
    var_id = sfselect(file_id, sfn2index(file_id, "ParticleRadius"))
    status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)
    number_waterradii = dim_sizes(1)
    allocate (water_radii(number_waterradii))
    edge(1) = number_waterradii
    status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), water_radii)
    status = sfendacc(var_id)


! Read the 2D arrays 


   start = 0
   stride = 1
   edge(1) = number_wavelengths
   edge(2) = number_waterradii

   allocate(extinction_water(number_wavelengths,number_waterradii))
   allocate(singlescattering_water(number_wavelengths,number_waterradii))
   allocate(asymmetry_water(number_wavelengths, number_waterradii))
   allocate(truncation_factor_water(number_wavelengths,number_waterradii))
   allocate(phase_fun_norm_constant_water(number_wavelengths,number_waterradii))

   var_id = sfselect(file_id, sfn2index(file_id, "ExtinctionCoefficient"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), extinction_water)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "SingleScatterAlbedo"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), singlescattering_water)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "Phase: AsymmetryParameter"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), asymmetry_water)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "TruncationFactor"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), truncation_factor_water)
   status = sfendacc(var_id)

   var_id = sfselect(file_id, sfn2index(file_id, "PhaseFuncNormConstant"))
   status = sfrdata(var_id, start(1:2),stride(1:2),edge(1:2), phase_fun_norm_constant_water)
   status = sfendacc(var_id)


! Read the 3D arrays

	
   start = 0
   stride = 1
   edge(1) = number_taus
   edge(2) = number_wavelengths
   edge(3) = number_waterradii
   
   allocate(spherical_albedo_water(number_taus,number_wavelengths,number_waterradii))
   var_id = sfselect(file_id, sfn2index(file_id, "SphericalAlbedo"))
   status = sfrdata(var_id, start(1:3),stride(1:3),edge(1:3), spherical_albedo_water)
   status = sfendacc(var_id)


 
   status = sfend(file_id)



   allocate(int_reflectance_water(number_taus+1, number_wavelengths,number_waterradii))
   allocate(int_reflectance_ice(number_taus+1, number_wavelengths,number_iceradii))

   allocate(int_reflectance_water_sdev(number_taus+1, number_wavelengths,number_waterradii))
   allocate(int_reflectance_ice_sdev(number_taus+1, number_wavelengths,number_iceradii))

   allocate(int_refl_water_sdev_wspeed(number_taus+1, number_wavelengths,number_waterradii, number_wind_speed))
   allocate(int_refl_ice_sdev_wspeed(number_taus+1, number_wavelengths,number_iceradii, number_wind_speed))

   allocate(int_reflectance_water_wspeed(number_taus+1, number_wavelengths,number_waterradii, number_wind_speed))
   allocate(int_reflectance_ice_wspeed(number_taus+1, number_wavelengths,number_iceradii, number_wind_speed))


   allocate(int_fluxupwater_sensor(number_taus, number_wavelengths,number_waterradii))
   allocate(int_fluxdownwater_solar(number_taus, number_wavelengths,number_waterradii))
   allocate(int_fluxdownwater_sensor(number_taus, number_wavelengths,number_waterradii))
   allocate(int_fluxupice_sensor(number_taus, number_wavelengths,number_iceradii))
   allocate(int_fluxdownice_solar(number_taus, number_wavelengths,number_iceradii))
   allocate(int_fluxdownice_sensor(number_taus, number_wavelengths,number_iceradii))

   allocate(int_cloud_emissivity_water(number_taus+1, 2, number_waterradii))
   allocate(int_cloud_emissivity_ice(number_taus+1, 2, number_iceradii))
   allocate(int_surface_emissivity_water(number_taus+1, 2, number_waterradii))
   allocate(int_surface_emissivity_ice(number_taus+1, 2, number_iceradii))

   allocate(int_cloud_emissivity_water_sdev(number_taus+1, 2, number_waterradii))
   allocate(int_cloud_emissivity_ice_sdev(number_taus+1, 2, number_iceradii))
   allocate(int_surface_emissivity_water_sdev(number_taus+1, 2, number_waterradii))
   allocate(int_surface_emissivity_ice_sdev(number_taus+1, 2, number_iceradii))

   allocate(int_cloud_emis_water_wspeed(number_taus+1, 2, number_waterradii, number_wind_speed))
   allocate(int_cloud_emis_ice_wspeed(number_taus+1, 2, number_iceradii, number_wind_speed))
   allocate(int_surface_emis_water_wspeed(number_taus+1, 2, number_waterradii, number_wind_speed))
   allocate(int_surface_emis_ice_wspeed(number_taus+1, 2, number_iceradii, number_wind_speed))

   allocate(int_cloud_emis_water_sdev_wspeed(number_taus+1, 2, number_waterradii, number_wind_speed))
   allocate(int_cloud_emis_ice_sdev_wspeed(number_taus+1, 2, number_iceradii, number_wind_speed))
   allocate(int_surface_emis_water_sdev_wspeed(number_taus+1, 2, number_waterradii, number_wind_speed))
   allocate(int_surface_emis_ice_sdev_wspeed(number_taus+1, 2, number_iceradii, number_wind_speed))


	allocate(rayleigh_liq(number_waterradii))
	allocate(rayleigh_ice(number_iceradii))
	


! Now for the phase function information

     file_id = sfstart(hdf_phase_lib, DFACC_READ)

	start = 0
	stride = 1
   
     ! get the ice phase func info
     var_id = sfselect(file_id, sfn2index(file_id, "ScatAnglesIce"))
     status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)
     number_phase_angles_ice = dim_sizes(1)
     allocate (phase_angles_ice(number_phase_angles_ice))
     edge(1) = number_phase_angles_ice
     status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), phase_angles_ice)
     status = sfendacc(var_id)
   
     var_id = sfselect(file_id, sfn2index(file_id, "ScatAnglesWater"))
     status = sfginfo(var_id, dummy_name, dummy_rank, dim_sizes, dummy_type, dummy_num_attr)
     number_phase_angles_water = dim_sizes(1)
     allocate (phase_angles_water(number_phase_angles_water))
     edge(1) = number_phase_angles_water
     status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), phase_angles_water)
     status = sfendacc(var_id)
		
	 
	allocate(phase_funcs_water(number_phase_angles_water, number_wavelengths, number_waterradii))
	allocate(phase_funcs_ice(number_phase_angles_ice, number_wavelengths, number_iceradii))

	start = 0
	stride = 1
	edge(1) = number_phase_angles_water
	edge(2) = number_wavelengths
	edge(3) = number_waterradii
	

	var_id = sfselect(file_id, sfn2index(file_id, "WaterPhaseFuncVals"))
	status = sfrdata(var_id, start(1:3), stride(1:3), edge(1:3), phase_funcs_water)
	status = sfendacc(var_id)
	
	
	edge(1) = number_phase_angles_ice
	edge(2) = number_wavelengths
	edge(3) = number_iceradii
	
	var_id = sfselect(file_id, sfn2index(file_id, "IcePhaseFuncVals"))
	status = sfrdata(var_id, start(1:3), stride(1:3), edge(1:3), phase_funcs_ice)
	status = sfendacc(var_id)
	
	status = sfend(file_id)

! now for the aerosol and rayleigh property information

	allocate(rayleigh_tau(number_wavelengths))
	allocate(aerosol_tau(number_wavelengths))
	allocate(aerosol_asym(number_wavelengths))
	allocate(aerosol_ssa(number_wavelengths))


    file_id = sfstart(hdf_ice_lib_ws3, DFACC_READ)
	
	 start = 0
	 stride = 1
	 edge(1) = number_wavelengths
	
	var_id = sfselect(file_id, sfn2index(file_id, "RayleighOpticalThickness"))
	status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), rayleigh_tau)
	status = sfendacc(var_id)
	
	var_id = sfselect(file_id, sfn2index(file_id, "AerosolOpticalThickness"))
	status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), aerosol_tau)
	status = sfendacc(var_id)
	
	var_id = sfselect(file_id, sfn2index(file_id, "AerosolAsymParameter"))
	status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), aerosol_asym)
	status = sfendacc(var_id)
	
	var_id = sfselect(file_id, sfn2index(file_id, "AerosolSSAlbedo"))
	status = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), aerosol_ssa)
	status = sfendacc(var_id)

	 status = sfend(file_id)

	call lib_init


 end subroutine readlibraries_base
 
 
 subroutine readlibraries_extra(hdf_water_lib, libnames_water, libnames_water_sdev, &
								hdf_ice_lib, libnames_ice, libnames_ice_sdev, debug,status)
   
   use GeneralAuxType
   use libraryarrays
   use libraryinterpolates
   use science_parameters
   use nonscience_parameters
   use mod06_run_settings
   use hdf_mod
   implicit none  

   character(*), intent(in)                   :: hdf_water_lib
   character(*), intent(in)					:: hdf_ice_lib
	character(*), intent(in) :: libnames_water(3), libnames_ice(3), &
								libnames_water_sdev(3), libnames_ice_sdev(3)

   logical, intent(in)                        :: debug
   integer,intent (out)                       :: status

   integer                                    :: file_id, var_id, var_index,i, j
   integer									:: start(6), stride(6), edge(6)
   character(MAX_SDS_NAME_LEN)                              :: dummy_name
   integer                                    :: dummy_type, dummy_num_attr, dummy_rank
   integer                                    :: dim_sizes(6),dim_sizes4d(4)

	real, dimension(:,:,:,:,:,:), allocatable :: temp_refl
	real, dimension(:,:,:,:), allocatable :: temp_emis

	integer :: max_bnd, min_bnd
	integer :: sensor_min, sensor_max, solar_min, solar_max, sza_min, vza_min, ra_min
	real, dimension(:,:,:,:), allocatable :: temp_flux, temp_flux2

	integer :: start_time, end_time, crate, cmax

	real :: deg_to_rad, cos_vza_min, cos_vza_max, cos_sza_min, cos_sza_max
	integer :: i1, i2, i3, i4, i5, i6, nws

	status = 0

	deg_to_rad = d2r
	

! we need to find angle range for the granule, so we can only read a limited amount from the library
! just enough to get the granule processed. 

! from the find_granule_angle_range routine we know just how wide the angle space we need. 
! it's fixed for sensor zenith and relative azimuth, but not for solar angle. However for MAS<TER>
! the relative azimuth space isn't fixed either and the sensor zenith limits are different. 

! start with the ice library first
	file_id = sfstart(hdf_ice_lib, DFACC_READ)

	start = 0
	stride = 1


! now get the relative azimuth and resize it as needed
	if (allocated(library_relative_azimuth)) deallocate(library_relative_azimuth)	
! resize the array	
	call find_bounds(relativeazimuth_all, min_rel_azimuth, max_rel_azimuth, min_bnd, max_bnd)
	number_relazimuth = max_bnd - min_bnd + 1

	ra_min = min_bnd
	allocate(library_relative_azimuth(number_relazimuth))
	library_relative_azimuth(1:number_relazimuth) = relativeazimuth_all(min_bnd : max_bnd)


! now get the solar zenith and resize it as needed
	if (allocated(library_solar_zenith)) deallocate(library_solar_zenith)
! resize the array	

	cos_sza_max = cos(max_solar_zenith*deg_to_rad)
	cos_sza_min = cos(min_solar_zenith*deg_to_rad)

	call find_bounds(solarzenith_all, cos_sza_max, cos_sza_min, min_bnd, max_bnd)
	number_solarzenith = max_bnd - min_bnd + 1
		
	sza_min = min_bnd
	allocate(library_solar_zenith(number_solarzenith))
	library_solar_zenith(1:number_solarzenith) = solarzenith_all(min_bnd : max_bnd)
	
! now get the sensor zenith and resize it as needed
	if (allocated(library_sensor_zenith)) deallocate(library_sensor_zenith)
! resize the array	

	cos_vza_max = cos(max_sensor_zenith*deg_to_rad)
	cos_vza_min = cos(min_sensor_zenith*deg_to_rad)

	call find_bounds(sensorzenith_all, cos_vza_max, cos_vza_min, min_bnd, max_bnd)

        number_sensorzenith = max_bnd - min_bnd + 1
	vza_min = min_bnd
	allocate(library_sensor_zenith(number_sensorzenith))
	library_sensor_zenith(1:number_sensorzenith) = sensorzenith_all(min_bnd : max_bnd)

! now get the flux angle and resize it as needed
	if (allocated(library_fluxsensorzenith)) deallocate(library_fluxsensorzenith)
	if (allocated(library_fluxsolarzenith)) deallocate(library_fluxsolarzenith)
! resize the array once for solar and once for sensor
	call find_bounds(solarzenith_flux_all, cos_vza_max, cos_vza_min, min_bnd, max_bnd)
	number_fluxsensorzenith = max_bnd - min_bnd + 1
	sensor_min = min_bnd
	sensor_max = max_bnd
	allocate(library_fluxsensorzenith(number_fluxsensorzenith))
	library_fluxsensorzenith(1:number_fluxsensorzenith) = solarzenith_flux_all(min_bnd : max_bnd)

	call find_bounds(solarzenith_flux_all, cos_sza_max, cos_sza_min, min_bnd, max_bnd)
	number_fluxsolarzenith = max_bnd - min_bnd + 1
	solar_min = min_bnd
	solar_max = max_bnd
	allocate(library_fluxsolarzenith(number_fluxsolarzenith))
	library_fluxsolarzenith(1:number_fluxsolarzenith) = solarzenith_flux_all(min_bnd : max_bnd)


! Read the 3D arrays


! Read the 4D arrays

	start = 0
	stride = 1	

	allocate(temp_flux2 (number_iceradii, number_wavelengths, number_taus, index_solarzenith_flux))

   if (allocated(flux_up_ice_solar)) deallocate(flux_up_ice_solar)
   if (allocated(flux_up_ice_sensor)) deallocate(flux_up_ice_sensor)
   if (allocated(flux_down_ice_solar)) deallocate(flux_down_ice_solar)
   if (allocated(flux_down_ice_sensor)) deallocate(flux_down_ice_sensor)


   allocate(flux_up_ice_solar (number_fluxsolarzenith, number_taus,number_wavelengths,number_iceradii))
   allocate(flux_down_ice_solar(number_fluxsolarzenith, number_taus,number_wavelengths,number_iceradii)) 
   allocate(flux_up_ice_sensor(number_fluxsensorzenith, number_taus,number_wavelengths,number_iceradii))
   allocate(flux_down_ice_sensor(number_fluxsensorzenith, number_taus,number_wavelengths,number_iceradii)) 

!	edge(1) = index_solarzenith_flux
!	edge(2) = number_taus
!	edge(3) = number_wavelengths
!	edge(4) = number_iceradii

	edge(1) = number_iceradii
	edge(2) = number_wavelengths
	edge(3) = number_taus
	edge(4) = index_solarzenith_flux

	
	
   var_id = sfselect(file_id, sfn2index(file_id, "ReflectedFlux"))
   status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_flux2)
   status = sfendacc(var_id)

   allocate(temp_flux (index_solarzenith_flux, number_taus,number_wavelengths,number_iceradii))
	do i1=1, index_solarzenith_flux
		do i2=1, number_taus
			do i3=1, number_wavelengths
				do i4=1, number_iceradii
				
					temp_flux(i1, i2, i3, i4) = temp_flux2(i4, i3, i2, i1)
				
				end do
			end do
		end do
	end do

	
	flux_up_ice_solar(:,:,:,:) = temp_flux(solar_min:solar_max, :,:,:)
	flux_up_ice_sensor(:,:,:,:) = temp_flux(sensor_min:sensor_max, :,:,:)

   var_id = sfselect(file_id, sfn2index(file_id, "TransmittedFlux"))
   status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_flux2)
   status = sfendacc(var_id)

	do i1=1, index_solarzenith_flux
		do i2=1, number_taus
			do i3=1, number_wavelengths
				do i4=1, number_iceradii
				
					temp_flux(i1, i2, i3, i4) = temp_flux2(i4, i3, i2, i1)
				
				end do
			end do
		end do
	end do


	deallocate(temp_flux2)

	flux_down_ice_solar(:,:,:,:) = temp_flux(solar_min:solar_max, :,:,:)
	flux_down_ice_sensor(:,:,:,:) = temp_flux(sensor_min:sensor_max, :,:,:)
	
	deallocate(temp_flux)

! Now we finally read the 6D reflectance array
! but we don't read the 11um reflectance, as there sure isn't any. 
! we do however need the 11um flux to do emissivity over land surfaces.
	if (allocated(library_reflectance_ice)) deallocate(library_reflectance_ice)


	edge(6) = number_sensorzenith
	edge(5) = number_solarzenith
	edge(4) = number_relazimuth
	edge(3) = number_taus
	edge(2) = number_wavelengths
	edge(1) = number_iceradii
	
   allocate(library_reflectance_ice( number_sensorzenith, number_solarzenith, number_relazimuth, &
								number_taus+1, number_wavelengths,number_iceradii, number_wind_speed+1  ))

	allocate(temp_refl( edge(1), edge(2), edge(3), edge(4), edge(5), edge(6) ))


	start = 0
	start(6) = vza_min-1
	start(5) = sza_min-1
	start(4) = ra_min-1


   var_id = sfselect(file_id, sfn2index(file_id, "MultiScatBDReflectance"))
   status = sfrdata(var_id, start,stride,edge, temp_refl) 
   status = sfendacc(var_id)

	do i1=1, number_sensorzenith
		do i2=1, number_solarzenith
			do i3=1, number_relazimuth
				do i4=1, number_taus
					do i5=1, number_wavelengths
						do i6=1, number_iceradii
						
							library_reflectance_ice(i1, i2, i3, i4, i5, i6, 1) = temp_refl(i6, i5, i4, i3, i2, i1)
						
						end do
					end do
				end do
			end do
		end do
	end do


	temp_refl = 0.


	if (allocated(library_reflectance_ice_sdev)) deallocate(library_reflectance_ice_sdev)
		allocate(library_reflectance_ice_sdev( number_sensorzenith, number_solarzenith, number_relazimuth, &
									number_taus+1, number_wavelengths,number_iceradii, number_wind_speed+1  ))
							 
	var_id = sfselect(file_id, sfn2index(file_id, "StdDevMultiScatBDReflectance"))
	status = sfrdata(var_id, start,stride,edge, temp_refl) 
	status = sfendacc(var_id)
								 
	nws = number_wind_speed + 1
	
	do i1=1, number_sensorzenith
		do i2=1, number_solarzenith
			do i3=1, number_relazimuth
				do i4=1, number_taus
					do i5=1, number_wavelengths
						do i6=1, number_iceradii
						
							library_reflectance_ice_sdev(i1, i2, i3, i4+1, i5, i6, nws) = temp_refl(i6, i5, i4, i3, i2, i1)
						
						end do
					end do
				end do
			end do
		end do
	end do

   status = sfend(file_id)

	deallocate(temp_refl)
								 		

	if (allocated(cloud_emissivity_ice)) deallocate(cloud_emissivity_ice)
	if (allocated(surface_emissivity_ice)) deallocate(surface_emissivity_ice)
	
	allocate(cloud_emissivity_ice( number_sensorzenith, number_taus+1, 2, &
			 number_iceradii, number_wind_speed))
	allocate(surface_emissivity_ice(number_sensorzenith, number_taus+1, 2, &
			 number_iceradii, number_wind_speed))


	if (allocated(cloud_emissivity_ice_sdev)) deallocate(cloud_emissivity_ice_sdev)
	if (allocated(surface_emissivity_ice_sdev)) deallocate(surface_emissivity_ice_sdev)
	

	allocate(cloud_emissivity_ice_sdev(number_sensorzenith, number_taus+1, 2, &
			 number_iceradii, number_wind_speed))
	allocate(surface_emissivity_ice_sdev(number_sensorzenith, number_taus+1, 2, &
			 number_iceradii, number_wind_speed))



! wind speed 3
	allocate(temp_refl( number_iceradii, &
						number_wavelengths, &
						number_taus+1, &
						number_relazimuth, &
						number_solarzenith, &
						number_sensorzenith))
	
!	number_sensorzenith, number_solarzenith, number_relazimuth, number_taus+1, &
!                                    number_wavelengths,number_iceradii  ))
									
	allocate(temp_emis( number_iceradii, &
							2,			&
						number_taus + 1, &
						number_sensorzenith))
	
!number_sensorzenith, number_taus+1, 2, number_iceradii))

if (DO_COX_MUNK) then ! DO_COX_MUNK

	do i=1, 3


		file_id = sfstart(libnames_ice(i), DFACC_READ)
	
		edge(6) = number_sensorzenith
		edge(5) = number_solarzenith
		edge(4) = number_relazimuth
		edge(3) = number_taus+1
		edge(2) = number_wavelengths
		edge(1) = number_iceradii

		start = 0
		start(6) = vza_min-1
		start(5) = sza_min-1
		start(4) = ra_min-1

		var_id = sfselect(file_id, sfn2index(file_id, "MultiScatBDReflectance"))
		status = sfrdata(var_id, start,stride,edge, temp_refl) 
		status = sfendacc(var_id)
	
!		library_reflectance_ice(:,:,:, :,:,:, i+1) = temp_refl(:,:,:, :,:,:)

	do i1=1, number_sensorzenith
		do i2=1, number_solarzenith
			do i3=1, number_relazimuth
				do i4=1, number_taus+1
					do i5=1, number_wavelengths
						do i6=1, number_iceradii
						
							library_reflectance_ice(i1, i2, i3, i4, i5, i6, i+1) = temp_refl(i6, i5, i4, i3, i2, i1)
						
						end do
					end do
				end do
			end do
		end do
	end do

	
		edge(4) = number_sensorzenith
		edge(3) = number_taus+1
		edge(2) = 2
		edge(1) = number_iceradii
	
		start = 0
		start(4) = vza_min-1
	
		var_id = sfselect(file_id, sfn2index(file_id, "CloudEmissivity"))
		status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_emis)
		status = sfendacc(var_id)
	
!		cloud_emissivity_ice(:,:,:,:,i) = temp_emis(:,:,:,:)

	do i1=1, number_sensorzenith
				do i4=1, number_taus+1
					do i5=1, 2
						do i6=1, number_iceradii
						
							cloud_emissivity_ice(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
						end do
					end do
				end do
	end do

	
		var_id = sfselect(file_id, sfn2index(file_id, "SurfaceEmissivity"))
		status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_emis)
		status = sfendacc(var_id)

!		surface_emissivity_ice(:,:,:,:,i) = temp_emis(:,:,:,:)

	do i1=1, number_sensorzenith
				do i4=1, number_taus+1
					do i5=1, 2
						do i6=1, number_iceradii
						
							surface_emissivity_ice(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
						end do
					end do
				end do
	end do


		
		status = sfend(file_id)

		file_id = sfstart(libnames_ice_sdev(i), DFACC_READ)
	
		edge(6) = number_sensorzenith
		edge(5) = number_solarzenith
		edge(4) = number_relazimuth
		edge(3) = number_taus+1
		edge(2) = number_wavelengths
		edge(1) = number_iceradii

		start = 0
		start(6) = vza_min-1
		start(5) = sza_min-1
		start(4) = ra_min-1

		var_id = sfselect(file_id, sfn2index(file_id, "StdDevMultiScatBDReflectance"))
		status = sfrdata(var_id, start,stride,edge, temp_refl) 
		status = sfendacc(var_id)
	
!		library_reflectance_ice_sdev(:,:,:, :,:,:,i) = temp_refl(:,:,:, :,:,:)


	do i1=1, number_sensorzenith
		do i2=1, number_solarzenith
			do i3=1, number_relazimuth
				do i4=1, number_taus+1
					do i5=1, number_wavelengths
						do i6=1, number_iceradii
						
							library_reflectance_ice_sdev(i1, i2, i3, i4, i5, i6, i) = temp_refl(i6, i5, i4, i3, i2, i1)
						
						end do
					end do
				end do
			end do
		end do
	end do

	
		edge(4) = number_sensorzenith
		edge(3) = number_taus+1
		edge(2) = 2
		edge(1) = number_iceradii
	
		start = 0
		start(4) = vza_min-1
	
		var_id = sfselect(file_id, sfn2index(file_id, "CloudEmissivity"))
		status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_emis)
		status = sfendacc(var_id)
	
!		cloud_emissivity_ice_sdev(:,:,:,:,i) = temp_emis(:,:,:,:)

	do i1=1, number_sensorzenith
				do i4=1, number_taus+1
					do i5=1, 2
						do i6=1, number_iceradii
						
							cloud_emissivity_ice_sdev(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
						end do
					end do
				end do
	end do


	
		var_id = sfselect(file_id, sfn2index(file_id, "SurfaceEmissivity"))
		status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_emis)
		status = sfendacc(var_id)

!		surface_emissivity_ice_sdev(:,:,:,:, i) = temp_emis(:,:,:,:)

	do i1=1, number_sensorzenith
				do i4=1, number_taus+1
					do i5=1, 2
						do i6=1, number_iceradii
						
							surface_emissivity_ice_sdev(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
						end do
					end do
				end do
	end do

		
		status = sfend(file_id)

		


	end do

endif ! DO_COX_MUNK

	deallocate(temp_refl, temp_emis)

!  Now that we're done with ice, get started with water...

!  Read all the water stuff there is to read.
   file_id = sfstart(hdf_water_lib, DFACC_READ)


! Read the 3D arrays



! Read the 4D arrays

	start = 0
	stride = 1	

	allocate(temp_flux2 (number_waterradii, number_wavelengths, number_taus, index_solarzenith_flux))

   if (allocated(flux_up_water_solar)) deallocate(flux_up_water_solar)
   if (allocated(flux_up_water_sensor)) deallocate(flux_up_water_sensor)
   if (allocated(flux_down_water_solar)) deallocate(flux_down_water_solar)
   if (allocated(flux_down_water_sensor)) deallocate(flux_down_water_sensor)


   allocate(flux_up_water_solar (number_fluxsolarzenith, number_taus,number_wavelengths,number_waterradii))
   allocate(flux_down_water_solar(number_fluxsolarzenith, number_taus,number_wavelengths,number_waterradii)) 
   allocate(flux_up_water_sensor(number_fluxsensorzenith, number_taus,number_wavelengths,number_waterradii))
   allocate(flux_down_water_sensor(number_fluxsensorzenith, number_taus,number_wavelengths,number_waterradii)) 

!	edge(1) = index_solarzenith_flux
!	edge(2) = number_taus
!	edge(3) = number_wavelengths
!	edge(4) = number_waterradii

	edge(1) = number_waterradii
	edge(2) = number_wavelengths
	edge(3) = number_taus
	edge(4) = index_solarzenith_flux

	
   allocate(temp_flux (index_solarzenith_flux, number_taus,number_wavelengths,number_waterradii))

	
   var_id = sfselect(file_id, sfn2index(file_id, "ReflectedFlux"))
   status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_flux2)
   status = sfendacc(var_id)
	
	
	do i1=1, index_solarzenith_flux
		do i2=1, number_taus
			do i3=1, number_wavelengths
				do i4=1, number_waterradii
				
					temp_flux(i1, i2, i3, i4) = temp_flux2(i4, i3, i2, i1)
				
				end do
			end do
		end do
	end do
	
	
	
	
	flux_up_water_solar(:,:,:,:) = temp_flux(solar_min:solar_max, :,:,:)
	flux_up_water_sensor(:,:,:,:) = temp_flux(sensor_min:sensor_max, :,:,:)

   var_id = sfselect(file_id, sfn2index(file_id, "TransmittedFlux"))
   status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_flux2)
   status = sfendacc(var_id)


	do i1=1, index_solarzenith_flux
		do i2=1, number_taus
			do i3=1, number_wavelengths
				do i4=1, number_waterradii
				
					temp_flux(i1, i2, i3, i4) = temp_flux2(i4, i3, i2, i1)
				
				end do
			end do
		end do
	end do

	deallocate(temp_flux2)


	flux_down_water_solar(:,:,:,:) = temp_flux(solar_min:solar_max, :,:,:)
	flux_down_water_sensor(:,:,:,:) = temp_flux(sensor_min:sensor_max, :,:,:)
	

	deallocate(temp_flux)

! Now we finally read the 6D reflectance array


	if (allocated(library_reflectance_water)) deallocate(library_reflectance_water)
	
	
	edge(6) = number_sensorzenith
	edge(5) = number_solarzenith
	edge(4) = number_relazimuth
	edge(3) = number_taus
	edge(2) = number_wavelengths
	edge(1) = number_waterradii
	
   allocate(library_reflectance_water( number_sensorzenith, number_solarzenith, number_relazimuth, &
								number_taus+1, number_wavelengths,number_waterradii, number_wind_speed+1  ))	

	allocate(temp_refl( edge(1), edge(2), edge(3), edge(4), edge(5), edge(6) ))


	start = 0
	start(6) = vza_min-1
	start(5) = sza_min-1
	start(4) = ra_min-1

   var_id = sfselect(file_id, sfn2index(file_id, "MultiScatBDReflectance"))
   status = sfrdata(var_id, start,stride,edge, temp_refl) 
   status = sfendacc(var_id)


!	library_reflectance_water(:,:,1:number_taus,:,:,:, 1 ) = temp_refl(:,:,:, :,:,:)

	do i1=1, number_sensorzenith
		do i2=1, number_solarzenith
			do i3=1, number_relazimuth
				do i4=1, number_taus
					do i5=1, number_wavelengths
						do i6=1, number_waterradii
						
							library_reflectance_water(i1, i2, i3, i4, i5, i6, 1) = temp_refl(i6, i5, i4, i3, i2, i1)
						
						end do
					end do
				end do
			end do
		end do
	end do


	 
	temp_refl = 0.

	if (allocated(library_reflectance_water_sdev)) deallocate(library_reflectance_water_sdev)

		allocate(library_reflectance_water_sdev(number_sensorzenith, number_solarzenith, number_relazimuth, number_taus+1, &
                                 number_wavelengths,number_waterradii, number_wind_speed+1  ))
					 
	var_id = sfselect(file_id, sfn2index(file_id, "StdDevMultiScatBDReflectance"))
	status = sfrdata(var_id, start,stride,edge, temp_refl) 
	status = sfendacc(var_id)
								 
	nws = number_wind_speed + 1
	
	do i1=1, number_sensorzenith
		do i2=1, number_solarzenith
			do i3=1, number_relazimuth
				do i4=1, number_taus
					do i5=1, number_wavelengths
						do i6=1, number_waterradii
						
							library_reflectance_water_sdev(i1, i2, i3, i4+1, i5, i6, nws) = temp_refl(i6, i5, i4, i3, i2, i1)
						
						end do
					end do
				end do
			end do
		end do
	end do

   status = sfend(file_id)

	deallocate(temp_refl)
								 
								 
								 

	if (allocated(cloud_emissivity_water)) deallocate(cloud_emissivity_water)
	if (allocated(surface_emissivity_water)) deallocate(surface_emissivity_water)
	
	
	
	allocate(cloud_emissivity_water(number_sensorzenith, number_taus+1, 2, &
			 number_waterradii, number_wind_speed))
	allocate(surface_emissivity_water(number_sensorzenith, number_taus+1, 2, &
			 number_waterradii, number_wind_speed))

	if (allocated(cloud_emissivity_water_sdev)) deallocate(cloud_emissivity_water_sdev)
	if (allocated(surface_emissivity_water_sdev)) deallocate(surface_emissivity_water_sdev)
	allocate(cloud_emissivity_water_sdev(number_sensorzenith, number_taus+1, 2, &
			 number_waterradii, number_wind_speed))
	allocate(surface_emissivity_water_sdev(number_sensorzenith, number_taus+1, 2, &
			 number_waterradii, number_wind_speed))


! wind speed 3


!	allocate(temp_refl(number_sensorzenith, number_solarzenith, number_relazimuth, number_taus+1, &
!                                   number_wavelengths,number_waterradii  ))

!	allocate(temp_emis(number_sensorzenith, number_taus+1, 2, number_waterradii))


	allocate(temp_refl( number_waterradii, &
						number_wavelengths, &
						number_taus+1, &
						number_relazimuth, &
						number_solarzenith, &
						number_sensorzenith))
										
	allocate(temp_emis( number_waterradii, &
							2,			&
						number_taus + 1, &
						number_sensorzenith))
	


if (DO_COX_MUNK) then ! DO_COX_MUNK
	do i=1, 3
		
!		call system_clock(start_time, crate, cmax)

		file_id = sfstart(libnames_water(i), DFACC_READ)
	
		edge(6) = number_sensorzenith
		edge(5) = number_solarzenith
		edge(4) = number_relazimuth
		edge(3) = number_taus+1
		edge(2) = number_wavelengths
		edge(1) = number_waterradii

		start(6) = vza_min-1
		start(5) = sza_min-1
		start(4) = ra_min-1

		var_id = sfselect(file_id, sfn2index(file_id, "MultiScatBDReflectance"))
		status = sfrdata(var_id, start,stride,edge, temp_refl) 
		status = sfendacc(var_id)
	
!		library_reflectance_water(:,:,:,:,:,:, i+1)  = temp_refl(:,:,:, :,:,:)

	do i1=1, number_sensorzenith
		do i2=1, number_solarzenith
			do i3=1, number_relazimuth
				do i4=1, number_taus+1
					do i5=1, number_wavelengths
						do i6=1, number_waterradii
						
							library_reflectance_water(i1, i2, i3, i4, i5, i6, i+1) = temp_refl(i6, i5, i4, i3, i2, i1)
						
						end do
					end do
				end do
			end do
		end do
	end do

	
		edge(4) = number_sensorzenith
		edge(3) = number_taus+1
		edge(2) = 2
		edge(1) = number_waterradii
	
		start = 0
		start(4) = vza_min-1
	
		var_id = sfselect(file_id, sfn2index(file_id, "CloudEmissivity"))
		status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_emis)
		status = sfendacc(var_id)
	
!		cloud_emissivity_water(:,:,:,:, i) = temp_emis(:,:,:,:)

	do i1=1, number_sensorzenith
				do i4=1, number_taus+1
					do i5=1, 2
						do i6=1, number_waterradii
						
							cloud_emissivity_water(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
						end do
					end do
				end do
	end do

	
		var_id = sfselect(file_id, sfn2index(file_id, "SurfaceEmissivity"))
		status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_emis)
		status = sfendacc(var_id)

!		surface_emissivity_water(:,:,:,:, i) = temp_emis(:,:,:,:)

	do i1=1, number_sensorzenith
				do i4=1, number_taus+1
					do i5=1, 2
						do i6=1, number_waterradii
						
							surface_emissivity_water(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
						end do
					end do
				end do
	end do
	
		status = sfend(file_id)

		file_id = sfstart(libnames_water_sdev(i), DFACC_READ)
	
		edge(6) = number_sensorzenith
		edge(5) = number_solarzenith
		edge(4) = number_relazimuth
		edge(3) = number_taus+1
		edge(2) = number_wavelengths
		edge(1) = number_waterradii

		start(6) = vza_min-1
		start(5) = sza_min-1
		start(4) = ra_min-1

		var_id = sfselect(file_id, sfn2index(file_id, "StdDevMultiScatBDReflectance"))
		status = sfrdata(var_id, start,stride,edge, temp_refl) 
		status = sfendacc(var_id)
	
!		library_reflectance_water_sdev(:,:,:,:,:,:,i)  = temp_refl(:,:,:, :,:,:)

	do i1=1, number_sensorzenith
		do i2=1, number_solarzenith
			do i3=1, number_relazimuth
				do i4=1, number_taus+1
					do i5=1, number_wavelengths
						do i6=1, number_waterradii
						
							library_reflectance_water_sdev(i1, i2, i3, i4, i5, i6, i) = temp_refl(i6, i5, i4, i3, i2, i1)
						
						end do
					end do
				end do
			end do
		end do
	end do

	
		edge(4) = number_sensorzenith
		edge(3) = number_taus+1
		edge(2) = 2
		edge(1) = number_waterradii
	
		start = 0
		start(4) = vza_min-1
	
		var_id = sfselect(file_id, sfn2index(file_id, "CloudEmissivity"))
		status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_emis)
		status = sfendacc(var_id)
	
!		cloud_emissivity_water_sdev(:,:,:,:,i) = temp_emis(:,:,:,:)

	do i1=1, number_sensorzenith
				do i4=1, number_taus+1
					do i5=1, 2
						do i6=1, number_waterradii
						
							cloud_emissivity_water_sdev(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
						end do
					end do
				end do
	end do
	
		var_id = sfselect(file_id, sfn2index(file_id, "SurfaceEmissivity"))
		status = sfrdata(var_id, start(1:4),stride(1:4),edge(1:4), temp_emis)
		status = sfendacc(var_id)

!		surface_emissivity_water_sdev(:,:,:,:,i) = temp_emis(:,:,:,:)

	do i1=1, number_sensorzenith
				do i4=1, number_taus+1
					do i5=1, 2
						do i6=1, number_waterradii
						
							surface_emissivity_water_sdev(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
						end do
					end do
				end do
	end do
	
		status = sfend(file_id)

!		call system_clock(end_time, crate, cmax)
!		print*, "time elapsed: ", end_time - start_time

	end do
endif ! DO_COX_MUNK


	deallocate(temp_refl, temp_emis)


!	print*, "lamb: ", library_reflectance_water( 1, 1, 1, :, 4, 6, 1)
!	print*, "ws3: ",library_reflectance_water( 1, 1, 1, :, 4, 6, 2)
!	print*, "ws7: ",library_reflectance_water( 1, 1, 1, :, 4, 6, 3)
!	print*, "ws15: ",library_reflectance_water( 1, 1, 1, :, 4, 6, 4)

!	print*, "lamb: ", library_reflectance_ice( 1, 1, 1, :, 4, 6, 1)
!	print*, "ws3: ",library_reflectance_ice( 1, 1, 1, :, 4, 6, 2)
!	print*, "ws7: ",library_reflectance_ice( 1, 1, 1, :, 4, 6, 3)
!	print*, "ws15: ",library_reflectance_ice( 1, 1, 1, :, 4, 6, 4)

end subroutine readlibraries_extra
 
 
 

 subroutine allocate_arrays ( edge, meas_edge, st_iterX, st_iterY, status ) 

   use GeneralAuxType
   use core_arrays
   use nonscience_parameters
   use libraryarrays
   use hdf_mod
   use mod06_run_settings
   use science_parameters
   use specific_other
   
   implicit none
   
   integer, dimension(:), intent(in)     ::  edge, meas_edge
   integer,               intent (out)   :: status
   integer, intent(in) :: st_iterX, st_iterY

   integer                               :: checkvariable
   integer                               :: xdimension, meas_xdimension, ydimension, number_of_bands

   integer                               :: model_layers
   integer :: i, j
   
   logical       :: allocation_status, array_size_change
   
   status = 0
   model_layers = model_levels
   number_of_bands = size(bands)


   meas_xdimension = meas_edge(1)
   xdimension = edge(1)
   ydimension =  edge(2)
 
!  Use optical_thickness array as a size/allocation test to save time.
   allocation_status = allocated(optical_thickness_16_final) 

   array_size_change = .false.
   
   if (allocation_status) then
      if (size(optical_thickness_16_final,1) /= xdimension .or. &
          size(optical_thickness_16_final,2) /= ydimension .or. &
		  size(band_measurements, 1) /= meas_xdimension) then
         array_size_change  = .true.
      endif
   endif

   if (iterationX == st_iterX .and. iterationY == st_iterY) then 
   
    allocate(model_info(grid_xsize, grid_ysize))
   
	do i=1, grid_xsize
		do j=1, grid_ysize
		
			allocate(model_info(i,j)%mixr_profile(model_layers))
			allocate(model_info(i,j)%temp_profile(model_layers))
			allocate(model_info(i,j)%height_profile(model_layers))
			allocate(model_info(i,j)%pressure_profile(model_layers))
		
		end do
	end do
   
   endif


!  Core retrieval arrays   
   if(array_size_change) then

		deallocate(cloud_height_method)
		deallocate(irw_temperature)
		deallocate(optical_thickness_final, stat = checkvariable)

#ifndef SEVIRI_INST
		deallocate(optical_thickness_1621_final, stat = checkvariable)
		deallocate(effective_radius_21_final, stat = checkvariable)
		deallocate(effective_radius_1621_final, stat = checkvariable)
		deallocate(liquid_water_path_1621, stat = checkvariable)
		deallocate(liquid_water_path, stat = checkvariable)

		deallocate(optical_thickness_final_PCL, stat = checkvariable)
		deallocate(optical_thickness_1621_final_PCL, stat = checkvariable)
		deallocate(effective_radius_21_final_PCL, stat = checkvariable)
		deallocate(effective_radius_1621_final_PCL, stat = checkvariable)
		deallocate(liquid_water_path_1621_PCL, stat = checkvariable)
		deallocate(liquid_water_path_PCL, stat = checkvariable)

		deallocate(optical_thickness_error, stat = checkvariable)
		deallocate(effective_radius_21_error, stat = checkvariable)

		deallocate(liquid_water_path_error, stat = checkvariable)
		deallocate(optical_thickness_1621_error, stat = checkvariable)
		deallocate(effective_radius_1621_error, stat = checkvariable)
		deallocate(liquid_water_path_1621_error, stat = checkvariable)
		
		deallocate(cloud_mask_SPI, stat = checkvariable)
	
		deallocate(failure_metric, stat = checkvariable)
		deallocate(failure_metric_1621, stat = checkvariable)

		deallocate(atm_corr_refl, stat = checkvariable)

#endif

		deallocate(optical_thickness_16_final, stat = checkvariable)
		deallocate(optical_thickness_37_final, stat = checkvariable)
		deallocate(effective_radius_16_final, stat = checkvariable)
		deallocate(effective_radius_37_final, stat = checkvariable)
		deallocate(liquid_water_path_16, stat = checkvariable)
		deallocate(liquid_water_path_37, stat = checkvariable)

		deallocate(optical_thickness_16_final_PCL, stat = checkvariable)
		deallocate(optical_thickness_37_final_PCL, stat = checkvariable)
		deallocate(effective_radius_16_final_PCL, stat = checkvariable)
		deallocate(effective_radius_37_final_PCL, stat = checkvariable)
		deallocate(liquid_water_path_16_PCL, stat = checkvariable)
		deallocate(liquid_water_path_37_PCL, stat = checkvariable)

		deallocate(optical_thickness_16_error, stat = checkvariable)
		deallocate(optical_thickness_37_error, stat = checkvariable)
		deallocate(effective_radius_16_error, stat = checkvariable)
		deallocate(effective_radius_37_error, stat = checkvariable)

		deallocate(liquid_water_path_16_error, stat = checkvariable)
		deallocate(liquid_water_path_37_error, stat = checkvariable)
		deallocate(cloud_layer_flag, stat = checkvariable)
		deallocate(ml_test_flag, stat = checkvariable)
		deallocate(CSR_flag_array, stat = checkvariable)
		deallocate(precip_water_094, stat = checkvariable)

!******************
!		deallocate(clear_sky_btemp)
!		deallocate(clear_sky_rad)

		deallocate(latitude, stat = checkvariable)
		deallocate(longitude, stat = checkvariable)
		deallocate(sensor_zenith_angle, stat = checkvariable)
		deallocate(solar_zenith_angle, stat = checkvariable)
		deallocate(relative_azimuth_angle, stat = checkvariable)
		deallocate(band_measurements, stat = checkvariable)
		deallocate(band_uncertainty, stat = checkvariable)
	    deallocate(sensor_azimuth_angle, stat = checkvariable)
   		deallocate(solar_azimuth_angle, stat = checkvariable)
	
		deallocate(surface_albedo, stat = checkvariable)
		deallocate(surface_temperature, stat = checkvariable)
		deallocate(surface_emissivity_land, stat = checkvariable)
		deallocate(cloud_top_temperature, stat = checkvariable)
		deallocate(cloud_top_pressure, stat = checkvariable)
	    deallocate(cloud_phase_infrared, stat = checkvariable)
		deallocate(abovecloud_watervapor, stat = checkvariable)
		deallocate(column_ozone, stat = checkvariable)
		deallocate(cloudmask, stat = checkvariable)
		
		deallocate(cloud_mask_SPI, stat = checkvariable)

		deallocate(failure_metric_16, stat = checkvariable)
		deallocate(failure_metric_37, stat = checkvariable)	
		
		deallocate(cloudsummary, stat = checkvariable)


		call deallocate_extra
		deallocate(processing_information, stat = checkvariable)

	endif


#ifndef SEVIRI_INST
	if (.not. allocated(optical_thickness_final)) then 

		allocate (optical_thickness_final(xdimension,ydimension), stat = checkvariable)	 
		allocate (optical_thickness_1621_final(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_21_final(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_1621_final(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_1621(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path(xdimension,ydimension), stat = checkvariable)

		allocate (optical_thickness_final_PCL(xdimension,ydimension), stat = checkvariable)	 
		allocate (optical_thickness_1621_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_21_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_1621_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_1621_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_PCL(xdimension,ydimension), stat = checkvariable)

		allocate (optical_thickness_error(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_21_error(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_error(xdimension,ydimension), stat = checkvariable)
		allocate (optical_thickness_1621_error(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_1621_error(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_1621_error(xdimension,ydimension), stat = checkvariable)

		allocate(failure_metric(xdimension, ydimension), stat = checkvariable)
		allocate(failure_metric_1621(xdimension, ydimension), stat = checkvariable)

 		allocate(atm_corr_refl(set_albedo_bands, xdimension, ydimension), stat=checkvariable)
		allocate (cloud_mask_SPI(2,xdimension,ydimension), stat = checkvariable)

#else

	if (.not. allocated(optical_thickness_16_final)) then 

#endif

	
		allocate (optical_thickness_16_final(xdimension,ydimension), stat = checkvariable)	 
		allocate (optical_thickness_37_final(xdimension,ydimension), stat = checkvariable)	 
		allocate (effective_radius_16_final(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_37_final(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_16(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_37(xdimension,ydimension), stat = checkvariable)

		allocate (optical_thickness_16_final_PCL(xdimension,ydimension), stat = checkvariable)	 
		allocate (optical_thickness_37_final_PCL(xdimension,ydimension), stat = checkvariable)	 
		allocate (effective_radius_16_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_37_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_16_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_37_PCL(xdimension,ydimension), stat = checkvariable)
				
		allocate (optical_thickness_16_error(xdimension,ydimension), stat = checkvariable)
		allocate (optical_thickness_37_error(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_16_error(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_37_error(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_16_error(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_37_error(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_layer_flag(xdimension,ydimension), stat = checkvariable)
		allocate (ml_test_flag(xdimension,ydimension), stat = checkvariable)
		allocate (CSR_flag_array(xdimension,ydimension), stat = checkvariable)
		allocate(precip_water_094(xdimension,ydimension), stat = checkvariable)
 
		allocate(failure_metric_16(xdimension, ydimension), stat = checkvariable)
		allocate(failure_metric_37(xdimension, ydimension), stat = checkvariable)
 
 		 
!  Core input data arrays

		allocate (latitude(xdimension,ydimension), stat = checkvariable)
		allocate (longitude(xdimension,ydimension), stat = checkvariable)
		allocate (sensor_zenith_angle(xdimension,ydimension), stat = checkvariable)
		allocate (solar_zenith_angle(xdimension,ydimension), stat = checkvariable)
		allocate (sensor_azimuth_angle(xdimension,ydimension), stat = checkvariable)
		allocate (solar_azimuth_angle(xdimension,ydimension), stat = checkvariable)
		allocate (relative_azimuth_angle(xdimension,ydimension), stat = checkvariable)

		allocate (band_measurements(meas_xdimension,number_of_bands, ydimension), stat = checkvariable)
		allocate (band_uncertainty(meas_xdimension,set_albedo_bands, ydimension), stat = checkvariable)


!  Ancillary data array

		allocate (surface_albedo(xdimension,ydimension,set_albedo_bands), stat = checkvariable)
		allocate (surface_temperature(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_height_method(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_phase_infrared(xdimension,ydimension), stat = checkvariable)
		allocate (irw_temperature(xdimension,ydimension), stat = checkvariable)
		allocate (surface_emissivity_land(xdimension,ydimension, 2), stat = checkvariable)
		allocate (cloud_top_temperature(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_top_pressure(xdimension,ydimension), stat = checkvariable)
		allocate (abovecloud_watervapor(xdimension,ydimension), stat = checkvariable)
		allocate (column_ozone(xdimension,ydimension), stat = checkvariable)
		allocate (cloudmask(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_mask_SPI(2,xdimension,ydimension), stat = checkvariable)

!******************
!		allocate (clear_sky_btemp(xdimension, ydimension, 2))
!		allocate (clear_sky_rad(xdimension, ydimension, 2))
		

!  QA and Processing arrays

		allocate (cloudsummary(xdimension,ydimension), stat = checkvariable)

!  atmospheric correction

!		allocate (transmittance_stddev(xdimension,ydimension,7), stat = checkvariable)
!		allocate (transmittance_twoway(xdimension,ydimension,7), stat = checkvariable)
!		allocate (meandelta_trans(xdimension,ydimension,7), stat = checkvariable)
!		allocate (thermal_correction_twoway(xdimension,ydimension,2), stat = checkvariable)
!		allocate (thermal_correction_oneway(xdimension,ydimension,2), stat = checkvariable)
!		allocate (thermal_correction_bands(2), stat = checkvariable)
  
		call allocate_extra(xdimension, ydimension)

		allocate (processing_information(xdimension,ydimension), stat = checkvariable)
		

    endif
  
end subroutine allocate_arrays


integer function findpoint( vector, value)
   use GeneralAuxType

   implicit none

   real(single) ,intent(in) :: vector(:)
   real(single) ,intent(in) :: value
   integer                  ::temp(1)

   real(single) ,allocatable   :: localvector(:)

   allocate(localvector(size(vector)))

   localvector = vector
   call realsingle_s_where_equal(localvector,value)

   temp = maxloc(localvector)
   findpoint = temp(1)
   deallocate(localvector)

end function findpoint


subroutine init_qualitydata
   use GeneralAuxType
   use core_arrays, only: processing_information

!   product quality and retrieval processing QA flags
    processing_information(:,:)%optical_thickness_GC = 0_integer_onebyte
    processing_information(:,:)%optical_thickness_outofbounds = 0_integer_onebyte
    processing_information(:,:)%effective_radius_GC = 0_integer_onebyte
    processing_information(:,:)%water_path_GC = 0_integer_onebyte
    processing_information(:,:)%rayleigh_correction = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_1621 = 0_integer_onebyte

    processing_information(:,:)%path_and_outcome_16 = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_16_PCL = 0_integer_onebyte

    processing_information(:,:)%path_and_outcome_37 = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_37_PCL = 0_integer_onebyte

    processing_information(:,:)%path_and_outcome_PCL = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_1621_PCL = 0_integer_onebyte

    processing_information(:,:)%band_used_for_optical_thickness = 0_integer_onebyte
    processing_information(:,:)%optical_thickness_1621_GC = 0_integer_onebyte
    processing_information(:,:)%effective_radius_1621_GC = 0_integer_onebyte
    processing_information(:,:)%water_path_1621_GC= 0_integer_onebyte
    processing_information(:,:)%multi_layer_cloud = 0_integer_onebyte
    processing_information(:,:)%CSR_flag = 0_integer_onebyte
    processing_information(:,:)%ml_test_mark = 0_integer_onebyte

#ifdef SEVIRI_INST
    processing_information(:,:)%Tc_override = 0_integer_onebyte
#endif

end subroutine init_qualitydata

subroutine deallocate_cleanup(status)

	use core_arrays
	use libraryarrays
	use libraryinterpolates
	use specific_other
	use interpolate_libraries
	use science_parameters

	integer, intent(inout) :: status
	integer :: checkvariable

	integer :: i, j

	status = 0

   
	do i=1, grid_xsize
		do j=1, grid_ysize
		
			deallocate(model_info(i,j)%mixr_profile)
			deallocate(model_info(i,j)%temp_profile)
			deallocate(model_info(i,j)%height_profile)
			deallocate(model_info(i,j)%pressure_profile)
		
		end do
	end do

    deallocate(model_info)

#ifndef SEVIRI_INST
		deallocate(optical_thickness_final, stat = checkvariable)
		deallocate(optical_thickness_1621_final, stat = checkvariable)
		deallocate(effective_radius_21_final, stat = checkvariable)
		deallocate(effective_radius_1621_final, stat = checkvariable)
		deallocate(liquid_water_path_1621, stat = checkvariable)
		deallocate(liquid_water_path, stat = checkvariable)

		deallocate(optical_thickness_final_PCL, stat = checkvariable)
		deallocate(optical_thickness_1621_final_PCL, stat = checkvariable)
		deallocate(effective_radius_21_final_PCL, stat = checkvariable)
		deallocate(effective_radius_1621_final_PCL, stat = checkvariable)
		deallocate(liquid_water_path_1621_PCL, stat = checkvariable)
		deallocate(liquid_water_path_PCL, stat = checkvariable)

		deallocate(optical_thickness_error, stat = checkvariable)
		deallocate(effective_radius_21_error, stat = checkvariable)

		deallocate(liquid_water_path_error, stat = checkvariable)
		deallocate(optical_thickness_1621_error, stat = checkvariable)
		deallocate(effective_radius_1621_error, stat = checkvariable)
		deallocate(liquid_water_path_1621_error, stat = checkvariable)
		
		deallocate(cloud_mask_SPI, stat = checkvariable)
	
		deallocate(failure_metric, stat = checkvariable)
		deallocate(failure_metric_1621, stat = checkvariable)

		deallocate(atm_corr_refl, stat = checkvariable)

#endif


	
!  Core retrieval arrays   

		deallocate(optical_thickness_16_final, stat = checkvariable)
		deallocate(optical_thickness_37_final, stat = checkvariable)
		deallocate(effective_radius_16_final, stat = checkvariable)
		deallocate(effective_radius_37_final, stat = checkvariable)

     deallocate(liquid_water_path_16, stat = checkvariable)
     deallocate(liquid_water_path_37, stat = checkvariable)


		deallocate(optical_thickness_16_final_PCL, stat = checkvariable)
		deallocate(optical_thickness_37_final_PCL, stat = checkvariable)
		deallocate(effective_radius_16_final_PCL, stat = checkvariable)
		deallocate(effective_radius_37_final_PCL, stat = checkvariable)
		deallocate(liquid_water_path_16_PCL, stat = checkvariable)
		deallocate(liquid_water_path_37_PCL, stat = checkvariable)



     deallocate(optical_thickness_16_error, stat = checkvariable)
     deallocate(optical_thickness_37_error, stat = checkvariable)
     deallocate(effective_radius_16_error, stat = checkvariable)
     deallocate(effective_radius_37_error, stat = checkvariable)
	 deallocate(liquid_water_path_16_error, stat = checkvariable)
	 deallocate(liquid_water_path_37_error, stat = checkvariable)

     deallocate(cloud_layer_flag, stat = checkvariable)
     deallocate(ml_test_flag, stat = checkvariable)
     deallocate(CSR_flag_array, stat = checkvariable)
     deallocate(precip_water_094, stat = checkvariable)
     
	deallocate(failure_metric_16, stat = checkvariable)
	deallocate(failure_metric_37, stat = checkvariable)
     	
!  Core input data arrays
     deallocate(latitude, stat = checkvariable)
     deallocate(longitude, stat = checkvariable)
     deallocate(sensor_zenith_angle, stat = checkvariable)
     deallocate(solar_zenith_angle, stat = checkvariable)
     deallocate(relative_azimuth_angle, stat = checkvariable)
     deallocate(band_measurements, stat = checkvariable)
     deallocate(band_uncertainty, stat = checkvariable)
	    deallocate(sensor_azimuth_angle, stat = checkvariable)
   		deallocate(solar_azimuth_angle, stat = checkvariable)
	 
     deallocate(surface_albedo, stat = checkvariable)
     deallocate(surface_temperature, stat = checkvariable)
     deallocate(surface_emissivity_land, stat = checkvariable)
     deallocate(cloud_top_temperature, stat = checkvariable)
     deallocate(cloud_top_pressure, stat = checkvariable)
	 deallocate(cloud_phase_infrared, stat = checkvariable)
     deallocate(abovecloud_watervapor, stat = checkvariable)
     deallocate(column_ozone, stat = checkvariable)
     deallocate(cloudmask, stat = checkvariable)
	 deallocate(cloud_mask_SPI, stat = checkvariable)
!     deallocate(model_height_profile, stat = checkvariable)
!     deallocate(model_temp_profile, stat = checkvariable)
!     deallocate(model_water_profile, stat = checkvariable)


!***********
!	deallocate(clear_sky_btemp)
!	deallocate(clear_sky_rad)

!  QA and Processing arrays
     deallocate(cloudsummary, stat = checkvariable)
     deallocate(processing_information, stat = checkvariable)

!  atmospheric correction

!     deallocate(transmittance_stddev, stat = checkvariable)
!     deallocate(transmittance_twoway, stat = checkvariable)
!     deallocate(meandelta_trans, stat = checkvariable)
!     deallocate(thermal_correction_twoway, stat = checkvariable)
!     deallocate(thermal_correction_oneway, stat = checkvariable)
!     deallocate(thermal_correction_bands, stat = checkvariable)

		deallocate(cloud_height_method)
		deallocate(irw_temperature)

! library arrays

     deallocate(transmit_correction_table)
     deallocate(transmit_stddev_table)
 
	 deallocate(ice_radii)
	 deallocate(library_taus)
	 deallocate(water_radii)

!  finally. Let's allocate the arrays we'll need later 
     deallocate(extinction_water, asymmetry_water, singlescattering_water, &
					truncation_factor_water, phase_fun_norm_constant_water)
     deallocate(spherical_albedo_water)
   
!  ice arrays, library, not interpolated


 	 deallocate(extinction_ice, asymmetry_ice, singlescattering_ice, &
					truncation_factor_ice, phase_fun_norm_constant_ice)
     deallocate(spherical_albedo_ice)
   

	 deallocate(rayleigh_tau, aerosol_tau, aerosol_asym, aerosol_ssa)
   
  
!  arrays for library interpolation
     deallocate(int_reflectance_water, int_reflectance_water_sdev)
	 deallocate( int_fluxupwater_sensor, int_fluxdownwater_solar,int_fluxdownwater_sensor)
     deallocate(int_reflectance_ice , int_reflectance_ice_sdev )
	 deallocate ( int_fluxupice_sensor, int_fluxdownice_solar,int_fluxdownice_sensor)
	
     deallocate(int_cloud_emis_water_wspeed, int_surface_emis_water_wspeed)
	 deallocate(int_cloud_emis_ice_wspeed, int_surface_emis_ice_wspeed)

   deallocate(int_reflectance_water_wspeed)
   deallocate(int_reflectance_ice_wspeed)
   deallocate(int_refl_water_sdev_wspeed)
   deallocate(int_refl_ice_sdev_wspeed)

	 deallocate(int_cloud_emissivity_ice, int_surface_emissivity_ice)
	 deallocate(int_cloud_emissivity_ice_sdev, int_surface_emissivity_ice_sdev)



     deallocate(int_cloud_emis_water_sdev_wspeed, int_surface_emis_water_sdev_wspeed)
	 deallocate(int_cloud_emis_ice_sdev_wspeed, int_surface_emis_ice_sdev_wspeed)

	deallocate(rayleigh_liq, rayleigh_ice)
	
! phase function information
	 deallocate(phase_angles_water)
	 deallocate(phase_angles_ice)
     deallocate(phase_funcs_water)
     deallocate(phase_funcs_ice)
	 
	 call deallocate_extra
	 call lib_clean

	if (no_valid_data == 0) then 

  	 deallocate(library_solar_zenith)
	 deallocate(library_sensor_zenith)
	 deallocate(library_relative_azimuth)
	 deallocate(library_fluxsolarzenith)
     deallocate(library_fluxsensorzenith)

     deallocate( solarzenith_all, sensorzenith_all, &
                    relativeazimuth_all, solarzenith_flux_all)
	
	 deallocate(cloud_emissivity_ice, surface_emissivity_ice)
	 deallocate(cloud_emissivity_ice_sdev, surface_emissivity_ice_sdev)

     deallocate(flux_up_ice_solar, flux_down_ice_solar)
     deallocate(flux_up_ice_sensor, flux_down_ice_sensor)
     deallocate(library_reflectance_ice)  
     deallocate(library_reflectance_ice_sdev)  

     deallocate(flux_up_water_solar, flux_down_water_solar)
     deallocate(flux_up_water_sensor, flux_down_water_sensor)
 
	 deallocate(library_reflectance_water)
	 deallocate(library_reflectance_water_sdev)

	 deallocate(cloud_emissivity_water, surface_emissivity_water)
	 deallocate(cloud_emissivity_water_sdev, surface_emissivity_water_sdev)

	endif

end subroutine deallocate_cleanup


end module modis_frontend_module
