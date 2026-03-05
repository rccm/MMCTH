 module specific_other
 
 implicit none
 
 contains

  subroutine set_l1b_names(level1b_name)
  
   use names

   character(*), intent(inout) :: level1b_name(:)

   
! if HKM name is garbage, so be it, no harm, no foul, we replace it with "none"
	level1b_name(1) = Alevel1b_name(1)
	level1b_name(2) = AHKM_name
	if (index(level1b_name(2), "02HKM") == 0) level1b_name(2) = "none"
      
  end subroutine set_l1b_names

  subroutine set_IRW_channel(IRW_channel)
	use mod06_run_settings
	
	integer, intent(inout) :: IRW_channel
	
	IRW_channel = 1
  
  end subroutine set_IRW_channel

 subroutine set_esfc(os_flag_in, x, y, esfc, os_flag )
 
	use libraryinterpolates
	use science_parameters
	use core_arrays
 
	logical, intent(in) :: os_flag_in
	real, dimension(:), intent(inout) :: esfc
	logical, intent(inout) :: os_flag
	integer, intent(in):: x, y
	
	if (os_flag_in) then 
		os_flag = .true.
	else
		os_flag = .false.
	endif

	if (os_flag .and. COX_MUNK) then 
		esfc(1) = int_surface_emissivity_water(1,2,1)
		esfc(2) = int_surface_emissivity_water(1,1,1)
	else
		esfc(1) = surface_emissivity_land(x,y,2)
		esfc(2) = surface_emissivity_land(x,y,1)
	endif
 
 end subroutine set_esfc



! this subroutine is intentionally left blank
  subroutine set_cox_munk_albedo(albedo, lib_albedo)
 
	use mod06_run_settings
 
	real, dimension(:), intent(in) :: albedo
	real, dimension(:), intent(in) :: lib_albedo
	
 
 end subroutine set_cox_munk_albedo

  subroutine get_band_idx(idx21, idx37, channel_37, idx_11, idx_alb37, idx_alb16) 
	
	use mod06_run_settings
 
	integer, intent(inout) :: idx21, idx37, channel_37, idx_11, idx_alb37, idx_alb16
 
	idx21 = band_0213
	idx37 = band_0370-1
	idx_11 = band_0370 !(it's 7 in MODIS)
	channel_37 = set_bands(band_0370)
	idx_alb37 = band_0370 - 1
	idx_alb16 = band_0163
 
 end subroutine get_band_idx

 ! this subroutine is intentionally left blank
 subroutine get_channels
	
	use mod06_run_settings
 
 end subroutine get_channels

 ! this subroutine is intentionally left blank
 subroutine allocate_extra(xdim, ydim)
	
	integer, intent(in) :: xdim, ydim
 
 end subroutine allocate_extra

 ! this subroutine is intentionally left blank
 subroutine deallocate_extra
	
	use core_arrays
 
 end subroutine deallocate_extra

 subroutine get_data_dims(filename, start, stride, edge)

	use mod06_run_settings
 
   integer, dimension(:), intent(inout) :: start, stride, edge
   character(len=*), intent(in) :: filename

   start = set_start
   stride = set_stride
   edge = set_edge


 end subroutine get_data_dims


 ! this subroutine is intentionally left blank
 subroutine set_process_time(file_id)

	integer, intent(in) :: file_id
 
 end subroutine set_process_time

 subroutine set_PH_desert(surface, R040, thres)
 
	logical, intent(in) :: surface
	real, intent(in) :: R040
	real, intent(inout) :: thres
	
	if (surface .and. R040 < 0.5) thres = 9999.
	

 end subroutine set_PH_desert

 logical function set_ice_ratio(ice_ratio)
 
	real, intent(in) :: ice_ratio
	
	if (ice_ratio < 1.3) then 
		set_ice_ratio = .true.
	else
		set_ice_ratio = .false.
	endif
 
 end function set_ice_ratio
 
 
! this subroutine is intentionally left blank
 subroutine set_albedo
 
 end subroutine set_albedo

! this subroutine is intentionally left blank
 subroutine set_processing_extra(file_id)
 
	integer, intent(in) :: file_id
	
 end subroutine set_processing_extra
 

 subroutine check_datasize(l1b_filedata, start, stride, edge, status)

!   use GeneralAuxType
   use nonscience_parameters

   implicit none

	include "hdf.f90"
	include "dffunc.f90"

   integer, dimension(:), intent(in)                   :: l1b_filedata
   integer, dimension (2), intent(inout) :: start, edge, stride
   integer, intent(out)                  :: status

   integer      :: Scans_Per_Granule
   character*40 :: att_N, dtype
   integer      :: RTN

   integer :: attr_id, file_id

   status = success

!  get number of scans
   att_N = 'Number of Scans'
   dtype = 'INTEGER*4'
  
	file_id = l1b_filedata(1)
	attr_id = sffattr(file_id, att_N)
	
	RTN = sfrattr(file_id, attr_id, Scans_Per_Granule)
	print*, Scans_Per_Granule
  
  
   if (rtn /= 0) then
     call MODIS_SMF_SETDYNAMICMSG(1, &
                                  'MAPI function GMFIN for Number of Scans  failed', &
                                  'check_datasize')
     status = failure
   endif

   edge(2) = Scans_Per_Granule * 10

 end subroutine check_datasize
 
subroutine aggregate_statistics_1km

	use modis_sciencestructure
	use core_arrays

	implicit none
	
	integer :: i, j, wid, ht
	
	
	wid = size(optical_thickness_final, 1)
	ht = size(optical_thickness_final, 2)
	
	do j=1, ht
		do i=1, wid
		
			if (.not. cloudsummary(i,j)%ocean_surface) &
									Statistics_1km%land_fraction = Statistics_1km%land_fraction + 1
			if (cloudsummary(i,j)%snowice_surface) Statistics_1km%snow_fraction = Statistics_1km%snow_fraction + 1
			if (cloudsummary(i,j)%ocean_surface) Statistics_1km%water_fraction = Statistics_1km%water_fraction + 1
		
			
			if (cloudsummary(i,j)%watercloud) then
				if (optical_thickness_final(i,j) > 0.) &
					Statistics_1km%mean_liquid_tau = Statistics_1km%mean_liquid_tau + optical_thickness_final(i,j)
				if (effective_radius_21_final(i,j) > 0.) &
					Statistics_1km%mean_liquid_re21 = Statistics_1km%mean_liquid_re21 + effective_radius_21_final(i,j)
				if (cloud_top_pressure(i,j) > 0.) & 
					Statistics_1km%ctp_liquid = Statistics_1km%ctp_liquid + cloud_top_pressure(i,j)
				if (cloud_top_temperature(i,j) > 0.) &
					Statistics_1km%ctt_liquid = Statistics_1km%ctt_liquid + cloud_top_temperature(i,j)
			endif
			
			if (cloudsummary(i,j)%icecloud) then 
				if (optical_thickness_final(i,j) > 0.) &
					Statistics_1km%mean_ice_tau = Statistics_1km%mean_ice_tau + optical_thickness_final(i,j)
				if (effective_radius_21_final(i,j) > 0.) &
					Statistics_1km%mean_ice_re21 = Statistics_1km%mean_ice_re21 + effective_radius_21_final(i,j)
				if (cloud_top_pressure(i,j) > 0.) & 
					Statistics_1km%ctp_ice = Statistics_1km%ctp_ice + cloud_top_pressure(i,j)
				if (cloud_top_temperature(i,j) > 0.) &
					Statistics_1km%ctt_ice = Statistics_1km%ctt_ice + cloud_top_temperature(i,j)
			endif
			
			if (cloudsummary(i,j)%unknowncloud) then 
				if (cloud_top_pressure(i,j) > 0.) & 
					Statistics_1km%ctp_undetermined = Statistics_1km%ctp_undetermined + cloud_top_pressure(i,j)
				if (cloud_top_temperature(i,j) > 0.) &
					Statistics_1km%ctt_undetermined = Statistics_1km%ctt_undetermined + cloud_top_temperature(i,j)			
			endif
		
		end do
	end do
	

end subroutine aggregate_statistics_1km

 
 
 end module specific_other