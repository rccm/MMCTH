 module names
 
  use mod06_run_settings
 
  implicit none
  
  integer :: MYYEAR, MYDAY, MYTIME, MYMONTH, MYDATE, MYMSG
  character(len = 2) :: MYMISSION

  character(len = 1000) :: Alevel1b_name(set_number_of_bands),   &
                           Acloudmask_name,             &
                           Amod06_name,            &
                           Ageolocation_name,           &
                           Agdas_name,                  &
                           Agdas_name2,                  &
						   Aozone_name, &
                           Ancepice_name,               &
                           Anise_name, ASHIS_name, Aangle_name              
	character(len = 1000) :: Awater_library,              &
                           Aice_library, ADEM_name, AHKM_name
						   
	character(len =1000) :: Alibnames_ice(3), Alibnames_ice_sdev(3)					   
	character(len =1000) :: Alibnames_water(3), Alibnames_water_sdev(3)					   
	character(len = 1000) :: Aphase_library,				&
                           Atransmittance_library,      &
                           Asurfacealbedo_lib_659,      &
                           Asurfacealbedo_lib_858,      &
                           Asurfacealbedo_lib_124,      &
                           Asurfacealbedo_lib_164,      &
                           Asurfacealbedo_lib_21,       &
                           Aecosystem_data_name,        &
                           Asnowicealbedo_data_name
	character(len=500) :: ACT_lib_path, Aemissivity_name, ACSRBias_name, Atmp_dirname

	character(len=500) :: MY_TEXT_FILE
	integer, parameter :: MY_UNIT_LUN = 20
														    
 end module names
 
 
