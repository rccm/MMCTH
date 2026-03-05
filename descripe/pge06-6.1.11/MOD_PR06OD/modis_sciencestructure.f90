module modis_sciencestructure
  use GeneralAuxType

  type qualityanalysis

!   product quality and retrieval processing QA flags   
    integer(integer_onebyte) :: optical_thickness_GC
    integer(integer_onebyte) :: optical_thickness_outofbounds
    integer(integer_onebyte) :: effective_radius_GC
    integer(integer_onebyte) :: water_path_GC
    integer(integer_onebyte) :: rayleigh_correction
! since they have to go together as it is, let them sit together in one byte instead of two
    integer(integer_onebyte) :: path_and_outcome
    integer(integer_onebyte) :: path_and_outcome_1621

    integer*1 :: path_and_outcome_16
    integer*1 :: path_and_outcome_16_PCL

    integer*1 :: path_and_outcome_37
    integer*1 :: path_and_outcome_37_PCL
    
    integer*1 :: path_and_outcome_PCL
    integer*1 :: path_and_outcome_1621_PCL
    
    integer(integer_onebyte) :: band_used_for_optical_thickness
    integer(integer_onebyte) :: optical_thickness_1621_GC
    integer(integer_onebyte) :: effective_radius_1621_GC
    integer(integer_onebyte) :: water_path_1621_GC
    integer(integer_onebyte) :: multi_layer_cloud
    integer(integer_onebyte) :: CSR_flag
	integer(integer_onebyte) :: ml_test_mark
	
#ifdef SEVIRI_INST
	integer :: Tc_override
#endif	
	
  end type qualityanalysis
  
  type stat_type
  
	real :: retrieval_fraction 
	real :: land_fraction 
	real :: water_fraction 
	real :: snow_fraction 
	real :: cloud_fraction 
	real :: water_cloud_fraction 
  	real :: ice_cloud_fraction 
	real*8 :: mean_liquid_tau 
	real*8 :: mean_ice_tau
	real*8 :: mean_liquid_re21 
	real*8 :: mean_ice_re21
	real*8 :: ctp_liquid 
	real*8 :: ctp_ice
	real*8 :: ctp_undetermined
	real*8 :: ctt_liquid 
	real*8 :: ctt_ice
	real*8 :: ctt_undetermined
  
  end type stat_type
  
  type failed_type
  
  	integer*2 :: tau
  	integer*2 :: re
  	integer*2 :: RSS
  
  end type failed_type
  
  
  
end module modis_sciencestructure
