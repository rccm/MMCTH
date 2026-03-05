 module modis_io_module

 use shared_io_module

 implicit none

 
! include "hdf.f90"
! include "dffunc.f90"
 
 
 private

 public :: output_retrieval, get_modis_data_cube, write_statistics, writeqaarray_toolkit, write_int2_array

 integer*2, dimension(:), allocatable :: outputbuffer_twobyte
 
 integer :: out_xsize, out_ysize


 contains

  subroutine write_statistics(filename)

	use core_arrays

	character(*), intent(in) :: filename

	real :: data_to_write(17)
	integer :: file_id, err_code, vref, vid
	real, parameter :: fill_val = -999.
	character(len=100) :: ln, un	
	character(len=85) :: desc(19)
	
	data_to_write(1) = Statistics_1km%retrieval_fraction
	data_to_write(2) = Statistics_1km%land_fraction 
	data_to_write(3) = Statistics_1km%water_fraction 
	data_to_write(4) = Statistics_1km%snow_fraction 
	data_to_write(5) = Statistics_1km%cloud_fraction 
	data_to_write(6) = Statistics_1km%water_cloud_fraction 
  	data_to_write(7) = Statistics_1km%ice_cloud_fraction 
	data_to_write(8) = Statistics_1km%mean_liquid_tau 
	data_to_write(9) = Statistics_1km%mean_ice_tau
	data_to_write(10) = Statistics_1km%mean_liquid_re21 
	data_to_write(11) = Statistics_1km%mean_ice_re21
	data_to_write(12) = Statistics_1km%ctp_liquid 
	data_to_write(13) = Statistics_1km%ctp_ice
	data_to_write(14) = Statistics_1km%ctp_undetermined
	data_to_write(15) = Statistics_1km%ctt_liquid 
	data_to_write(16) = Statistics_1km%ctt_ice
	data_to_write(17) = Statistics_1km%ctt_undetermined
	
	file_id = hopen(filename, DFACC_WRITE, 0)
	
	err_code = vfstart(file_id)
	
	vref = vsffnd(file_id, "Statistics_1km")	
	vid = vsfatch(file_id, vref, "w")

	err_code = vsfwrt(vid, data_to_write, 17, 0)

	err_code = vsfdtch(vid)
	err_code = vfend(file_id)
	err_code = hclose(file_id)

! duplicate same information as an SDS
! we can't use MOD_PR06CR because CR automatically turns any 1D SDS into VData
! VData makes it difficult to see SDS attributes in common display tools
! without the attributes the data is useless
	file_id = sfstart(filename, DFACC_WRITE)
	vid = sfselect(file_id, sfn2index(file_id, "Statistics_1km_sds"))
	if (vid == -1) vid = sfcreate(file_id, "Statistics_1km_sds", DFNT_FLOAT, 1, (/ 17 /))
	err_code = sfwdata(vid, (/ 0 /), (/ 1 /), (/ 17 /), data_to_write)

	ln = "Granule Statistics for parameters at 1x1 resolution"
	un = "see description attribute"
	err_code = sfsattr(vid, "long_name", DFNT_CHAR, len(trim(ln)), ln )
	err_code = sfsattr(vid, "units", DFNT_CHAR, len(trim(un)), un )
	err_code = sfsattr(vid, "_FillValue", DFNT_FLOAT, 1, fill_val)

	desc(1:19)(1:85) = ' '

	desc(1) = " "
	desc(2) = "Statistics_1km:"
    desc(3) = "  1. Successful Retrieval Rate (%)"
    desc(4) = "  2. Land Cover Fraction (%)"
    desc(5) = "  3. Water Cover Fraction (%)"
   	desc(6) = "  4. Snow Cover Fraction (%)"
   	desc(7) = "  5. Cloud Cover Fraction (%)"
   	desc(8) = "  6. Water Cloud Detected (%)"
  	desc(9) = "  7. Ice Cloud Detected (%)"
   	desc(10) = "  8. Mean of Water Cloud Optical Thickness"
   	desc(11) = "  9. Mean of Ice Cloud Optical Thickness "
   	desc(12) = "  10. Mean of Water Cloud Effective Particle Radius (microns)"
   	desc(13) = "  11. Mean of Ice Cloud Effective Diameter (microns)"
   	desc(14) = "  12. Mean Liquid Water Cloud Top Pressure (mb)"
   	desc(15) = "  13. Mean Ice Cloud Top Pressure (mb)"
   	desc(16) = "  14. Mean Undetermined Cloud Top Pressure (mb)"
   	desc(17) = "  15. Mean Liquid Water Cloud Top Temperature (K)"
   	desc(18) = "  16. Mean Ice Cloud Top Temperature (K) "
   	desc(19) = "  17. Mean Undetermined Cloud Top Temperature (K)"
	
	desc(1:19)(85:85) = char(10)
	
	err_code = sfsattr(vid, "description", DFNT_CHAR, 19*85, desc)	

	err_code = sfendacc(vid)
	err_code = sfend(file_id)



  end subroutine write_statistics

 subroutine write_float_array(filedata, sdsname, start, stride, edge, data_to_write, status) 
 
	use nonscience_parameters
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(inout) :: data_to_write
 
 
	integer :: file_id, var_id, err_code, attr_id
	
	real :: fill_val
	
	fill_val = fillvalue_real
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
		
	attr_id = sffattr(var_id, "_FillValue")
	if (attr_id == -1) then
		attr_id = sfsattr(var_id, "_FillValue", DFNT_FLOAT, 1, fill_val)
	endif
	
	err_code = sfwdata(var_id, start, stride, edge, data_to_write)
	
	err_code = sfendacc(var_id)
	
	status = 0

 end subroutine write_float_array	 


 subroutine write_int2_array(filedata, sdsname, start, stride, edge, data_to_write, status) 
 
	use nonscience_parameters
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	integer*2, dimension(:,:,:), intent(inout) :: data_to_write
 
 
	integer :: file_id, var_id, err_code, attr_id
	
	real :: fill_val
	
	fill_val = fillvalue_real
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
			
	err_code = sfwdata(var_id, start, stride, edge, data_to_write)
	status = err_code
	
	err_code = sfendacc(var_id)
	

 end subroutine write_int2_array	 



 subroutine output_retrieval(mapi_filedata, &
                            filename,      &
                            currentscanX, currentscanY, nscansX, nscansY,  &
                            start, edge, stride, &
                            status)
   use GeneralAuxType
   use core_arrays
   use nonscience_parameters
   use mod06_run_settings
   use libraryarrays
   implicit none

   integer,      intent(in)          :: mapi_filedata(:), currentscanX, currentscanY, nscansX, nscansY
   integer,      intent(in)          :: start(:), edge(:), stride(:)
   character(*), intent (in)         :: filename

   integer,      intent (out)        :: status

   integer                              :: buffer_xsize, buffer_ysize, i,j, count
   integer(integer_twobyte)             :: fillint_twobyte
   integer(integer_twobyte),allocatable :: outputbuffer(:)
   real(double)                         :: scale, add_offset 
   integer :: localstart(3), localedge(3), localstride(3), xsize, ysize
   integer(integer_onebyte),allocatable :: quality_assurance_1km(:,:,:)
   real, allocatable :: retrieval_diff(:,:,:)
   
   integer*1, allocatable, dimension(:,:) :: cloud_phase_COP
   
   status = success
   xsize = size(optical_thickness_final,1)
   ysize = size(optical_thickness_final,2)

	out_xsize = xsize
	out_ysize = ysize

   allocate(quality_assurance_1km(9, xsize,ysize ))

   call set_quality_data(xsize, ysize)

   call convert_binary_qa(quality_assurance_1km,  &
                          status)

   call writeqaarray_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           quality_assurance_1km,   &
                           'Quality_Assurance_1km', &
                           status )
   deallocate(quality_assurance_1km)


   allocate(outputbuffer_twobyte(xsize*ysize))


   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           optical_thickness_final,   &
                           'Cloud_Optical_Thickness', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           optical_thickness_final_PCL,   &
                           'Cloud_Optical_Thickness_PCL', &
                           status )

   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           optical_thickness_1621_final,   &
                           'Cloud_Optical_Thickness_1621', &
                           status ) 
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           optical_thickness_1621_final_PCL,   &
                           'Cloud_Optical_Thickness_1621_PCL', &
                           status )

   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           optical_thickness_16_final,   &
                           'Cloud_Optical_Thickness_16', &
                           status ) 
                           
                           
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           optical_thickness_16_final_PCL,   &
                           'Cloud_Optical_Thickness_16_PCL', &
                           status )

   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           optical_thickness_37_final,   &
                           'Cloud_Optical_Thickness_37', &
                           status ) 
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           optical_thickness_37_final_PCL,   &
                           'Cloud_Optical_Thickness_37_PCL', &
                           status )


   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           effective_radius_21_final,   &
                           'Cloud_Effective_Radius', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           effective_radius_21_final_PCL,   &
                           'Cloud_Effective_Radius_PCL', &
                           status )


   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           effective_radius_16_final,   &
                           'Cloud_Effective_Radius_16', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           effective_radius_16_final_PCL,   &
                           'Cloud_Effective_Radius_16_PCL', &
                           status )

   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           effective_radius_37_final,   &
                           'Cloud_Effective_Radius_37', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           effective_radius_37_final_PCL,   &
                           'Cloud_Effective_Radius_37_PCL', &
                           status )


   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           effective_radius_1621_final,   &
                           'Cloud_Effective_Radius_1621', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           effective_radius_1621_final_PCL,   &
                           'Cloud_Effective_Radius_1621_PCL', &
                           status )

   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           liquid_water_path,   &
                           'Cloud_Water_Path', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           liquid_water_path_PCL,   &
                           'Cloud_Water_Path_PCL', &
                           status )
   
   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           liquid_water_path_16,   &
                           'Cloud_Water_Path_16', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           liquid_water_path_16_PCL,   &
                           'Cloud_Water_Path_16_PCL', &
                           status )

   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           liquid_water_path_37,   &
                           'Cloud_Water_Path_37', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           liquid_water_path_37_PCL,   &
                           'Cloud_Water_Path_37_PCL', &
                           status )
   
   
   call writefloatarray_toolkit(mapi_filedata,  &
                           filename,            &
                           currentscanX, currentscanY,         &
                           liquid_water_path_1621,   &
                          'Cloud_Water_Path_1621', &
                           status )
   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           liquid_water_path_1621_PCL,   &
                           'Cloud_Water_Path_1621_PCL', &
                           status )


   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           optical_thickness_error,   &
                           'Cloud_Optical_Thickness_Uncertainty', &
                           status )

   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           effective_radius_21_error,   &
                           'Cloud_Effective_Radius_Uncertainty', &
                           status )


   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           effective_radius_16_error,   &
                           'Cloud_Effective_Radius_Uncertainty_16', &
                           status )


   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           effective_radius_37_error,   &
                           'Cloud_Effective_Radius_Uncertainty_37', &
                           status )



   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           liquid_water_path_error,   &
                           'Cloud_Water_Path_Uncertainty', &
                           status )

   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           liquid_water_path_16_error,   &
                           'Cloud_Water_Path_Uncertainty_16', &
                           status )


   call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           liquid_water_path_37_error,   &
                           'Cloud_Water_Path_Uncertainty_37', &
                           status )


    call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                            optical_thickness_1621_error,   &
                            'Cloud_Optical_Thickness_Uncertainty_1621', &
                            status )
 
    call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                            optical_thickness_16_error,   &
                            'Cloud_Optical_Thickness_Uncertainty_16', &
                            status )
 
    call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                            optical_thickness_37_error,   &
                            'Cloud_Optical_Thickness_Uncertainty_37', &
                            status )
 

    call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                            effective_radius_1621_error,   &
                            'Cloud_Effective_Radius_Uncertainty_1621', &
                            status )
 
    call writeint2array_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                            liquid_water_path_1621_error,   &
                            'Cloud_Water_Path_Uncertainty_1621', &
                            status )

    call writefloatarray_toolkit(mapi_filedata,  &
                            filename,            &
                           currentscanX, currentscanY,         &
                            precip_water_094,   &
                            'Above_Cloud_Water_Vapor_094', &
                            status )


    call writefloatarray_toolkit(mapi_filedata,  &
                            filename,            &
                           currentscanX, currentscanY,         &
                            irw_temperature,   &
                            'IRW_Low_Cloud_Temperature_From_COP', &
                            status )

	deallocate(outputbuffer_twobyte)
	

	call write_failed_array(mapi_filedata, &
							currentscanX, currentscanY, &
							failure_metric, &
							'Retrieval_Failure_Metric', &
							status)

	call write_failed_array(mapi_filedata, &
							currentscanX, currentscanY, &
							failure_metric_16, &
							'Retrieval_Failure_Metric_16', &
							status)

	call write_failed_array(mapi_filedata, &
							currentscanX, currentscanY, &
							failure_metric_37, &
							'Retrieval_Failure_Metric_37', &
							status)

	call write_failed_array(mapi_filedata, &
							currentscanX, currentscanY, &
							failure_metric_1621, &
							'Retrieval_Failure_Metric_1621', &
							status)
 
   call writebytearray_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           cloud_layer_flag,   &
                           'Cloud_Multi_Layer_Flag', &
                           status )


	allocate(cloud_phase_COP(xsize,ysize))
	cloud_phase_COP = processing_information%path_and_outcome
	where(cloud_phase_COP > 4) 
		cloud_phase_COP = cloud_phase_COP - 8
	end where

   call writebytearray_toolkit(mapi_filedata,  &
                           currentscanX, currentscanY,         &
                           cloud_phase_COP,   &
                           'Cloud_Phase_Optical_Properties', &
                           status )

	deallocate(cloud_phase_COP)

!    call write3Dfloatarray(mapi_filedata,  &
!                            filename,            &
!                           currentscanX, currentscanY,         &
!                            cloud_mask_SPI,   &
!                            'Cloud_Mask_SPI', &
!                            status )

    call write3Dfloatarray(mapi_filedata,  &
                            filename,            &
                           currentscanX, currentscanY,         &
                            atm_corr_refl,   &
                            'Atm_Corr_Refl', &
                            status )

	if (currentscanX == nscansX .and. currentscanY == nscansY) then 
	
		localstart = 0
		localstride = 1
		localedge(1) = number_wavelengths
		localedge(2) = number_iceradii
	
		call write_float_array(mapi_filedata, "Extinction_Efficiency_Ice", &
			localstart(1:2), localstride(1:2), localedge(1:2), extinction_ice, status)	
		call write_float_array(mapi_filedata, "Single_Scatter_Albedo_Ice", &
			localstart(1:2), localstride(1:2), localedge(1:2), singlescattering_ice, status)
		call write_float_array(mapi_filedata, "Asymmetry_Parameter_Ice", &
			localstart(1:2), localstride(1:2), localedge(1:2), asymmetry_ice, status)
	
		localedge(2) = number_waterradii

		call write_float_array(mapi_filedata, "Extinction_Efficiency_Liq", &
			localstart(1:2), localstride(1:2), localedge(1:2), extinction_water, status)	
		call write_float_array(mapi_filedata, "Single_Scatter_Albedo_Liq", &
			localstart(1:2), localstride(1:2), localedge(1:2), singlescattering_water, status)
		call write_float_array(mapi_filedata, "Asymmetry_Parameter_Liq", &
			localstart(1:2), localstride(1:2), localedge(1:2), asymmetry_water, status)
	
	
	endif

 end subroutine output_retrieval
 
 
 subroutine read_L1B (level1b_filedata,  &
                      band, Cal_type_is_refl, xdimension, ydimension, start_val,&
                      level1b_buffer, uncertainty_buffer, uncertain_spec, uncertain_scale)

	use nonscience_parameters
					  				  
	implicit none
	
	include "hdf.f90"
	include "dffunc.f90"
	
	integer, dimension(:), intent(in) :: level1b_filedata
	logical, intent(in) ::  Cal_type_is_refl
	integer, intent(in) ::  band, xdimension, ydimension, start_val
	integer*1, dimension(:,:), intent(inout) :: uncertainty_buffer
	real, dimension(:,:), intent(inout) :: level1b_buffer
	real, intent(inout) :: uncertain_spec, uncertain_scale
	
	integer :: file_id, var_id, err_code, start(3), stride(3), edge(3), attr_id
	real*4 :: scale_factors(20), offsets(20)
	real*4 :: spec_unc(20), unc_scale(20)
	character(len=200) :: SDS_name, SDS_unc_name
	integer :: band_index
	integer*2 :: valid_range(2)
	integer*2, dimension(:,:), allocatable :: temp
	
	
        if (band >= 1 .and. band <=2) then 
			SDS_name = "EV_250_Aggr1km_RefSB"
			SDS_unc_name = "EV_250_Aggr1km_RefSB_Uncert_Indexes"  
		endif
        if (band >= 3 .and. band <=7) then
			SDS_name = "EV_500_Aggr1km_RefSB" 
			SDS_unc_name = "EV_500_Aggr1km_RefSB_Uncert_Indexes"          
		endif 
		if (band >=8 .and. band<=19 .or. band == 26) then
			SDS_name = "EV_1KM_RefSB" 
			SDS_unc_name = "EV_1KM_RefSB_Uncert_Indexes"
		endif
		if (band >= 20 .and. band <= 36 .and. band /=26) then 
			SDS_name = "EV_1KM_Emissive"
			SDS_unc_name = "EV_1KM_Emissive_Uncert_Indexes"
		endif
! first we need to find where to start
        if (band >=1 .and. band <= 2) band_index = band-1  
        if (band >=3 .and. band <= 7) band_index = band-3 
        if (band >=8 .and. band <=12) band_index = band-8  
        if (band >=13 .and. band <=15) band_index = 2*band-21 
        if (band >=15 .and. band <=19) band_index = band-6 
        if (band >=20 .and. band <=25) band_index = band-20  
        if (band >=27 .and. band <=36) band_index = band-21 
        if (band == 26) band_index = 14
		
!	start(1) = 0
!	start(2) = (scan_number-1)*10
	start(2) = 0
	start(1) = start_val !(scan_number-1)*99
	start(3) = band_index
	
	stride = 1
	
	edge(1) = xdimension
	edge(2) = ydimension
	edge(3) = 1
	
	allocate(temp(xdimension, ydimension))
	
	file_id = level1b_filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, SDS_name))
	
	err_code = sfrdata(var_id, start, stride, edge, temp)
	
	attr_id = sffattr(var_id, "valid_range")
	err_code = sfrattr(var_id, attr_id, valid_range)
	
	
	if (Cal_type_is_refl) then 
		attr_id = sffattr(var_id, "reflectance_scales")
		err_code = sfrattr(var_id, attr_id, scale_factors)
		attr_id = sffattr(var_id, "reflectance_offsets")
		err_code = sfrattr(var_id, attr_id, offsets)
            
	else
		attr_id = sffattr(var_id, "radiance_scales")
		err_code = sfrattr(var_id, attr_id, scale_factors)
		attr_id = sffattr(var_id, "radiance_offsets")
		err_code = sfrattr(var_id, attr_id, offsets)
        
	endif

	
	where (temp >= valid_range(1) .and. temp <= valid_range(2))
		level1b_buffer = (temp*1.0 - offsets(band_index+1)) *scale_factors(band_index+1)
	elsewhere 
		level1b_buffer = fillvalue_real
	end where
	
	
	deallocate(temp)
	
	err_code = sfendacc(var_id)
	
	var_id = sfselect(file_id, sfn2index(file_id, SDS_unc_name))
	
	attr_id = sffattr(var_id, "specified_uncertainty")
	err_code = sfrattr(var_id, attr_id, spec_unc)
	attr_id = sffattr(var_id, "scaling_factor")
	err_code = sfrattr(var_id, attr_id, unc_scale)
	
	err_code = sfrdata(var_id, start, stride, edge, uncertainty_buffer)
	
	err_code = sfendacc(var_id)
					  
					  
	uncertain_spec = spec_unc(band_index+1)
	uncertain_scale = unc_scale(band_index+1)				  
					  
 end subroutine read_L1B
 
 
 
 subroutine get_modis_data_cube ( level1b_filedata,       &
                                 geolocation_filedata, &
                                 start, edge, stride, meas_start, meas_edge, scan_number, debug, status)

!  Core sds of interest in the MODIS level 1B file read to fill the following:
!   latitude (MOD03)
!   longitude (MOD03)
!   band_measurements
!   solar_zenith_angle (MOD03)
!   sensor_zenith_angle (MOD03)
!   solar_azimuth_angle (MOD03)
!   sensor_azimuth_angle (MOD03)


   use GeneralAuxType
   use hdf_mod
   use core_arrays
   use science_parameters
   use nonscience_parameters
   use mod06_run_settings

   implicit none


   integer, dimension (2), intent (in)       :: start, edge, stride, meas_start, meas_edge
   integer, dimension(:), intent(in)         :: level1b_filedata, geolocation_filedata
   integer, intent(in)                       :: scan_number
  
   logical,                intent (in)       :: debug
   integer,                intent (out)      :: status

   integer                                   :: numberofbands ,checkvariable,  i, j, k

   integer                                   :: xdimension, meas_xdimension, ydimension
   real,dimension(:,:),allocatable           :: level1b_buffer, sza_temp
   integer(integer_onebyte),dimension(:,:), allocatable :: uncertainty_buffer
   logical                                   :: errorflag, useoffset
   real :: unc_spec, unc_scale
	integer :: unc_idx

   logical ::   Cal_type_is_refl


   status = success
   numberofbands = size(bands)

!  get level 1b data
 
   meas_xdimension = meas_edge(1)
   xdimension =   edge(1)
   ydimension =   edge(2)
 
   allocate(level1b_buffer(meas_xdimension, ydimension))
   allocate(uncertainty_buffer(meas_xdimension, ydimension))
   allocate(sza_temp(meas_xdimension, ydimension))
   
   

   solar_constant_37 = 10.9295 ! average of Terra and Aqua values from table in Platnick and Fontenla (2008)
!	solar_constant_37 = 8.82301 ! from formula in Platnick and Fontenla (2008)
	
   do i = 1, numberofbands
     if(bands(i) < 20 .or. bands(i) == 26) then
       Cal_type_is_refl = .true.
     else
       Cal_type_is_refl = .false.
     endif

     call read_L1B(level1b_filedata,     &
                   bands(i), Cal_type_is_refl,       &
                   meas_xdimension,ydimension,&
				   meas_start(1), &
                   level1b_buffer,       &
                   uncertainty_buffer, unc_spec, unc_scale )
     band_measurements(:,i,:) = level1b_buffer(:,:)

	unc_idx = 0

	if (bands(i) == 1) unc_idx = band_0065
	if (bands(i) == 2) unc_idx = band_0086
	if (bands(i) == 5) unc_idx = band_0124
	if (bands(i) == 6) unc_idx = band_0163
	if (bands(i) == 7) unc_idx = band_0213
	if (bands(i) == channel_37um) unc_idx = band_0370 - 1
	
	if (unc_idx /= 0) then 
	
		band_uncertainty(:,unc_idx,:) = uncertainty_buffer(:,:)
		spec_uncertain(unc_idx) = unc_spec
		uncertain_sf(unc_idx) = unc_scale

	endif

   enddo

   deallocate(level1b_buffer, uncertainty_buffer)

	no_valid_data = 0

!  get full resolution Latitude and Longitude arrays

	call read_float_array(geolocation_filedata, "Latitude", start, stride, edge, latitude, status)

	call read_float_array(geolocation_filedata, "Longitude", start, stride, edge, longitude, status) 


	useoffset = .false.
	call read_int_array(geolocation_filedata, "SensorZenith", start, stride, edge, useoffset, sensor_zenith_angle, status)
	
	call read_int_array(geolocation_filedata, "SensorAzimuth", start, stride, edge, useoffset, sensor_azimuth_angle, status) 
							  
							  
	call read_int_array(geolocation_filedata, "SolarZenith", meas_start, stride, meas_edge, useoffset, sza_temp, status)

	call read_int_array(geolocation_filedata, "SolarAzimuth", start, stride, edge, useoffset, solar_azimuth_angle, status) 



   do i = 1, numberofbands
     if(bands(i) < 20 .or. bands(i) == 26) then
        where (band_measurements(:,i,:) > 0.) 
           band_measurements(:,i,:) = band_measurements(:,i,:) / cos(d2r*sza_temp)
        end where
     endif
   enddo
   
		
	if (scan_number == 1) then 
		solar_zenith_angle(:,:) = sza_temp(1:edge(1), :)
	else
		solar_zenith_angle(:,:) = sza_temp(2:edge(1)+1, :)
	endif

	deallocate(sza_temp)


!  calculate the relative azimuth

   relative_azimuth_angle =  solar_azimuth_angle + 180. - sensor_azimuth_angle

	do j = 1,  ydimension
		do i = 1, xdimension
	 
       if (relative_azimuth_angle(i,j) <= 0.) relative_azimuth_angle(i,j) = -relative_azimuth_angle(i,j) 
       if (relative_azimuth_angle(i,j) > 180.) relative_azimuth_angle(i,j) = 360. - relative_azimuth_angle(i,j)  
     enddo
   enddo

   relative_azimuth_angle = abs(relative_azimuth_angle)
   where(relative_azimuth_angle > 180.) relative_azimuth_angle =360. - relative_azimuth_angle
   
   
   
   max_rel_azimuth = 0.
   min_rel_azimuth = 999999.
   
   
   max_sensor_zenith = 0.
   min_sensor_zenith = 99999.
   
   max_solar_zenith = 0.
   min_solar_zenith = 99999.
   
   do j=1, ydimension
		do i=1, xdimension
		
		if (relative_azimuth_angle(i,j) < min_rel_azimuth .and. relative_azimuth_angle(i,j) >= 0.) &
			min_rel_azimuth = relative_azimuth_angle(i,j)
		if (relative_azimuth_angle(i,j) > max_rel_azimuth .and. relative_azimuth_angle(i,j) >= 0.) &
			max_rel_azimuth = relative_azimuth_angle(i,j)
						
		if (solar_zenith_angle(i,j) < min_solar_zenith .and. solar_zenith_angle(i,j) >= 0. ) &
			min_solar_zenith = solar_zenith_angle(i,j)
		if (solar_zenith_angle(i,j) > max_solar_zenith .and. solar_zenith_angle(i,j) >= 0.) &
			max_solar_zenith = solar_zenith_angle(i,j)

					

	end do
   end do


! the sensor zenith is constant along a column
	do j=1, ydimension
		
		do i=1, xdimension
			if (sensor_zenith_angle(i,j) < 0.) exit ! bad data, get a different line
			if (sensor_zenith_angle(i,j) < min_sensor_zenith .and. sensor_zenith_angle(i,j) >= 0.) &
				min_sensor_zenith = sensor_zenith_angle(i,j)
			if (sensor_zenith_angle(i,j) > max_sensor_zenith .and. sensor_zenith_angle(i,j) >= 0.) &
				max_sensor_zenith = sensor_zenith_angle(i,j)
		end do
		
		if (i >= xdimension) then
			exit ! we are done, a good line of data
		endif

	end do
   	
	
 end subroutine get_modis_data_cube

 
 
 
 
 subroutine writeqaarray_toolkit(mapi_filedata, &
                                   currentscanX, currentscanY,  &
                                   array,         &
                                   sdsname,       &
                                   status)

  use GeneralAuxType
  use nonscience_parameters
  use mod06_run_settings

  implicit none

  include "hdf.f90"
  include "dffunc.f90"


  integer(integer_onebyte), intent(in)           :: array(:,:,:)
  integer, intent(in)        :: mapi_filedata(:), currentscanX, currentscanY
  character(*), intent (in)  :: sdsname

  integer,      intent (out)        :: status 

  real(double)     :: scale_factor, add_offset
  integer          :: buffer_xsize, buffer_ysize, i,j, count, k, buffer_zsize
  integer          :: rank, dimsize(3), nms
  integer          :: localstart(3), localedge(3), localstride(3)
  character(len=MAX_SDS_NAME_LEN)     :: groupname, data_type, attribute
!  integer(integer_onebyte),allocatable :: outputbuffer_onebyte(:)

  integer :: file_id, var_id, err_code

  buffer_zsize = size(array,1)
  buffer_xsize = size(array,2)
  buffer_ysize = size(array,3)

  localstride = 1

  rank = 3

  localstart(1) = 0
  localstart(2) = (currentscanX-1)*set_tilesize
  localstart(3) = (currentscanY-1)*set_tilesize

  localedge(1) = buffer_zsize
  localedge(2) = buffer_xsize
  localedge(3) = buffer_ysize



  file_id = mapi_filedata(1)

  var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)) )

  err_code = sfwdata(var_id, localstart, localstride, localedge, array)
  
  err_code = sfendacc(var_id)

  status = success

 end subroutine writeqaarray_toolkit


 subroutine writeint2array_toolkit(mapi_filedata, &
									currentscanX, currentscanY, &
                                   array,         &
                                   sdsname,       &
                                   status)
  use nonscience_parameters
  use mod06_run_settings

  implicit none

	include "hdf.f90"
	include "dffunc.f90"

  integer*2, intent(in)           :: array(:,:)
  integer, intent(in)        :: mapi_filedata(:), currentscanX, currentscanY
  character(*), intent (in)  :: sdsname

  integer,      intent (out)        :: status 

  integer          :: localstart(3), localedge(3), localstride(3)

  integer :: file_id, var_id, err_code
  
  localstride = 1
  
  file_id = mapi_filedata(1)
  var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)) )

  localstart(2) = (currentscanY-1)*set_tilesize
  localstart(1) = (currentscanX-1)*set_tilesize
  localstart(3) = 0

  localedge(1) = out_xsize
  localedge(2) = out_ysize
  localedge(3) = 1
	
  status = success
  
  err_code = sfwdata(var_id, localstart, localstride, localedge, array) 
  if (err_code /= 0) status = failure 
  
  err_code = sfendacc(var_id)
  

 
 end subroutine writeint2array_toolkit


 subroutine write_failed_array(mapi_filedata, &
									currentscanX, currentscanY, &
                                   failed_data,         &
                                   sdsname,       &
                                   status)
  use nonscience_parameters
  use mod06_run_settings
  use modis_sciencestructure

  implicit none

	include "hdf.f90"
	include "dffunc.f90"

  type(failed_type), intent(in)           :: failed_data(:,:)
  integer, intent(in)        :: mapi_filedata(:), currentscanX, currentscanY
  character(*), intent (in)  :: sdsname

  integer,      intent (out)        :: status 

  integer          :: localstart(3), localedge(3), localstride(3)

  integer :: file_id, var_id, err_code
  integer*2, dimension(:,:,:), allocatable :: temp_buffer
  
  allocate(temp_buffer(3, out_xsize, out_ysize))
  
  temp_buffer(1,:,:) = failed_data(:,:)%tau
  temp_buffer(2,:,:) = failed_data(:,:)%re
  temp_buffer(3,:,:) = failed_data(:,:)%RSS
  
  localstride = 1
  
  file_id = mapi_filedata(1)
  var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)) )

  localstart(1) = 0
  localstart(2) = (currentscanX-1)*set_tilesize
  localstart(3) = (currentscanY-1)*set_tilesize

  localedge(1) = 3
  localedge(2) = out_xsize
  localedge(3) = out_ysize
	
  status = success
  
  err_code = sfwdata(var_id, localstart, localstride, localedge, temp_buffer) 
  if (err_code /= 0) status = failure 
  
  err_code = sfendacc(var_id)
  
  deallocate(temp_buffer)
 
 end subroutine write_failed_array






 subroutine writebytearray_toolkit(mapi_filedata, &
									currentscanX, currentscanY, &
                                   array,         &
                                   sdsname,       &
                                   status)
  use nonscience_parameters
  use mod06_run_settings

  implicit none

	include "hdf.f90"
	include "dffunc.f90"

  integer*1, intent(in)           :: array(:,:)
  integer, intent(in)        :: mapi_filedata(:), currentscanX, currentscanY
  character(*), intent (in)  :: sdsname

  integer,      intent (out)        :: status 

  integer          :: localstart(3), localedge(3), localstride(3)

  integer :: file_id, var_id, err_code
  
  localstride = 1
  
  file_id = mapi_filedata(1)
  var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)) )

  localstart(2) = (currentscanY-1)*set_tilesize
  localstart(1) = (currentscanX-1)*set_tilesize
  localstart(3) = 0

  localedge(1) = out_xsize
  localedge(2) = out_ysize
  localedge(3) = 1
	
  status = success
  
  err_code = sfwdata(var_id, localstart, localstride, localedge, array) 
  if (err_code /= 0) status = failure 
  
  err_code = sfendacc(var_id)
  

 
 end subroutine writebytearray_toolkit






 subroutine writefloatarray_toolkit(mapi_filedata, &
                                   filename,      &
									currentscanX, currentscanY, &
                                   array,         &
                                   sdsname,       &
                                   status)

  use GeneralAuxType
  use nonscience_parameters
  use mod06_run_settings

  implicit none

	include "hdf.f90"
	include "dffunc.f90"

  real, intent(in)           :: array(:,:)
  integer, intent(in)        :: mapi_filedata(:), currentscanX, currentscanY
  character(*), intent (in)  :: filename, sdsname

  integer,      intent (out)        :: status 

  real(double)     :: scale_factor, add_offset
  integer          :: buffer_xsize, buffer_ysize, i,j, count
  integer          :: rank, dimsize(2), nms
  integer          :: localstart(3), localedge(3), localstride(3)
  character(len=MAX_SDS_NAME_LEN)     :: groupname, data_type, attribute

  integer :: file_id, var_id, attr_id, err_code
  character(len=400) :: dummy_name
  integer :: dtype, dummy, dims(10)
  real, allocatable :: outputbuffer_real4byte(:)
  
  localstride = 1
  

  file_id = mapi_filedata(1)
  var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)) )


  buffer_xsize = out_xsize
  buffer_ysize = out_ysize
  rank = 2

  
	attr_id = sffattr(var_id, "scale_factor")
	err_code = sfrattr(var_id, attr_id, scale_factor)

	attr_id = sffattr(var_id, "add_offset")
	err_code = sfrattr(var_id, attr_id, add_offset)


   err_code = sfginfo(var_id, dummy_name, dummy, dims, dtype, dummy)

  localstart(2) = (currentscanY-1)*set_tilesize
  localstart(1) = (currentscanX-1)*set_tilesize
  localstart(3) = 0
   
  localedge(1) = buffer_xsize
  localedge(2) = buffer_ysize
  localedge(3) = 1

  status = failure

   if (dtype == DFNT_INT16 .or. dtype == DFNT_UINT16) data_type = 'INTEGER*2'
   if (dtype == DFNT_FLOAT32) data_type = 'REAL*4'


  if (data_type(1:9) == 'INTEGER*2') then
  status = success
  count = 0
  do j = 1, buffer_ysize
     do i = 1, buffer_xsize
       count = count + 1
       if (array(i,j) < 0. ) then
          outputbuffer_twobyte(count) = -9999
       else
          outputbuffer_twobyte(count) = nint(array(i,j)/scale_factor+add_offset)
          ! guard against potentially zeroing-out a very small quantity. Basically indicate that it's something other than 0. 
		  if (outputbuffer_twobyte(count) == 0 .and. add_offset == 0.) outputbuffer_twobyte(count) = 1
       endif
     enddo
  enddo

  localstart(2) = (currentscanY-1)*set_tilesize
  localstart(1) = (currentscanX-1)*set_tilesize
  localstart(3) = 0

  localedge(1) = buffer_xsize
  localedge(2) = buffer_ysize
  localedge(3) = 1


   err_code = sfwdata(var_id, localstart, localstride, localedge, outputbuffer_twobyte)
	
  endif

  if (data_type(1:6) == 'REAL*4') then
    allocate(outputbuffer_real4byte(buffer_xsize*buffer_ysize))
    outputbuffer_real4byte(:) = 0.0
    status = success
    count = 0
    do j = 1, buffer_ysize
      do i = 1, buffer_xsize
        count = count + 1
        if (array(i,j) < 0. ) then
            outputbuffer_real4byte(count) = -999.
        else
            outputbuffer_real4byte(count) = array(i,j)/scale_factor+add_offset
        endif
      enddo
    enddo

    err_code = sfwdata(var_id, localstart, localstride, localedge, outputbuffer_real4byte)
	
  endif
 
 
 
 err_code = sfendacc(var_id)
 
 end subroutine writefloatarray_toolkit
 
 
 
 subroutine write3Dfloatarray(mapi_filedata, &
                                   filename,      &
									currentscanX, currentscanY, &
                                   array,         &
                                   sdsname,       &
                                   status)

  use GeneralAuxType
  use nonscience_parameters
  use mod06_run_settings

  implicit none

	include "hdf.f90"
	include "dffunc.f90"

  real, intent(in)           :: array(:,:,:)
  integer, intent(in)        :: mapi_filedata(:), currentscanX, currentscanY
  character(*), intent (in)  :: filename, sdsname

  integer,      intent (out)        :: status 

  real(double)     :: scale_factor, add_offset
  integer          :: buffer_xsize, buffer_ysize, buffer_zsize, i,j, k, count
  integer          :: rank, dimsize(2), nms
  integer          :: localstart(3), localedge(3), localstride(3)
  character(len=MAX_SDS_NAME_LEN)     :: groupname, data_type, attribute
  integer(integer_twobyte),allocatable :: outputbuffer_2b(:)
  integer(integer_onebyte),allocatable :: outputbuffer_onebyte(:)
  real, allocatable :: outputbuffer_real4byte(:)

  integer :: file_id, var_id, attr_id, err_code
  character(len=400) :: dummy_name
  integer :: dtype, dummy, dims(10)
  
	character(len = 256) :: attr_string, description
	integer*2 :: fillValue_integer2, range_integer2(2)
	integer*1 :: fillValue_integer1, range_integer1(2)
  
 
	buffer_zsize = size(array,1)
	buffer_xsize = size(array,2)
	buffer_ysize = size(array,3)
	rank = 3

  
  
  localstride = 1
  

  localstart(1) = 0
  localstart(2) = (currentscanX-1)*set_tilesize
  localstart(3) = (currentscanY-1)*set_tilesize
   
  localedge(1) = buffer_zsize
  localedge(2) = buffer_xsize
  localedge(3) = buffer_ysize


  file_id = mapi_filedata(1)
  var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)) )

	attr_id = sffattr(var_id, "scale_factor")
	err_code = sfrattr(var_id, attr_id, scale_factor)

	attr_id = sffattr(var_id, "add_offset")
	err_code = sfrattr(var_id, attr_id, add_offset)


   err_code = sfginfo(var_id, dummy_name, dummy, dims, dtype, dummy)


  status = failure

  if (dtype == DFNT_INT16 .or. dtype == DFNT_UINT16) data_type = 'INTEGER*2'
  if (dtype == DFNT_FLOAT32) data_type = 'REAL*4'


  if (data_type(1:9) == 'INTEGER*2') then
  allocate(outputbuffer_2b(buffer_xsize*buffer_ysize*buffer_zsize))
  status = success
  count = 0
  do j = 1, buffer_ysize
     do i = 1, buffer_xsize
		do k=1, buffer_zsize
       count = count + 1
       if (array(k,i,j) < 0. ) then
          outputbuffer_2b(count) = -9999
       else
          outputbuffer_2b(count) = nint(array(k,i,j)/scale_factor+add_offset)
       endif
     enddo
  enddo
  end do

   err_code = sfwdata(var_id, localstart, localstride, localedge, outputbuffer_2b)
	
  deallocate(outputbuffer_2b)
  endif 

  if (data_type(1:6) == 'REAL*4') then
    allocate(outputbuffer_real4byte(buffer_zsize*buffer_xsize*buffer_ysize))
    outputbuffer_real4byte(:) = 0.5
    status = success
    count = 0
    do j = 1, buffer_ysize
       do i = 1, buffer_xsize
         do k = 1, buffer_zsize
           count = count + 1
           if (array(k,i,j) == -999. ) then
             outputbuffer_real4byte(count) = -999.
           else
             outputbuffer_real4byte(count) = array(k,i,j)/scale_factor+add_offset
           endif
         enddo
      enddo
    enddo

    err_code = sfwdata(var_id, localstart, localstride, localedge, outputbuffer_real4byte)


    deallocate(outputbuffer_real4byte)
  endif

  
  err_code = sfendacc(var_id)
 
 end subroutine write3Dfloatarray
  
 end module modis_io_module
