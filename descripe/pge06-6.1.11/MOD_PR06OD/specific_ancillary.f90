module specific_ancillary

	implicit none

contains

! This subroutine is basically a stub as far as MODIS is concerned
 subroutine get_specific_ancillary(	mod06input_filedata,   &
                                mod06input_name,       &
								pressure_name, &
								temperature_name, &
								ctm_name, &
								sfc_temp_name, &
								cpi_name, &
                                mod06_start,       &
                                mod06_stride,      &
                                mod06_edge, EXECUTE_STANDARD, REPLACE_ALBEDO)
						 
   integer,     intent(inout),dimension(:)  :: mod06_start, mod06_stride, mod06_edge, mod06input_filedata
   character(*), intent(in) :: mod06input_name
	logical, intent(inout) :: EXECUTE_STANDARD, REPLACE_ALBEDO
 
	character(*), intent(inout) :: pressure_name, temperature_name, ctm_name, sfc_temp_name, cpi_name
 
#ifdef CT_1KM

	pressure_name = "cloud_top_pressure_1km"
	temperature_name = "cloud_top_temperature_1km"
	ctm_name = "cloud_top_method_1km"
	sfc_temp_name = "surface_temperature_1km"
	cpi_name = "Cloud_Phase_Infrared_1km"

    mod06_start  = mod06_start
    mod06_edge   = mod06_edge

#else

	pressure_name = "Cloud_Top_Pressure"
	temperature_name = "Cloud_Top_Temperature"
	ctm_name = "Cloud_Height_Method"
	sfc_temp_name = "Surface_Temperature"
	cpi_name = "Cloud_Phase_Infrared"

    mod06_start  = floor(real(mod06_start)/5.)
    mod06_edge   = floor(real(mod06_edge)/5.)


#endif

    mod06_stride = 1

	EXECUTE_STANDARD = .true.
	REPLACE_ALBEDO = .false.

 
 end subroutine get_specific_ancillary


 subroutine get_cloud_top_properties(mapi_filedata, &
                                    filename, &
									pressure_sds, &
									temperature_sds, &
									ctm_sds, &
									sfc_temp_sds, &
									cpi_sds, &
                                    start,    &
                                    stride,   &
                                    edge,     &
									start_1km, &
									edge_1km, &
                                    cloud_top_pressure, &
                                    cloud_top_temperature, &
									cloud_height_method, &
									surface_temperature, &
									cloud_phase_infrared, &
									EXECUTE_STANDARD, &
                                    status)

  use hdf_mod
  use GeneralAuxType
   use nonscience_parameters
   use mod06_run_settings
  use shared_io_module,only: read_int_array, read_byte_array
  implicit none

  character(*), intent(in)                  :: filename, pressure_sds, temperature_sds, ctm_sds, sfc_temp_sds, cpi_sds
  integer,      intent(in)                  :: mapi_filedata(:)
  integer,      intent(inout), dimension(:) :: start, stride, edge, start_1km, edge_1km
  real(single), dimension(:,:), intent(out) :: cloud_top_pressure, cloud_top_temperature, surface_temperature
  integer*1, dimension(:,:), intent(out ) :: cloud_height_method, cloud_phase_infrared
  integer,      intent(out)                 :: status
  logical, intent(in) :: EXECUTE_STANDARD

  integer(integer_twobyte)                  :: fill_value
  integer                                   :: i, j, local_i, local_j, dimsizes(2)
  real(single), dimension(:,:), allocatable :: ctp_temp_array, ctt_temp_array, sfc_temp_array
  integer*1, dimension(:,:), allocatable :: ctm_temp_array, cpi_temp_array
  logical :: offset
  integer :: localstride(2)
  

  cloud_top_pressure = 0.
  cloud_top_temperature = 0.
  cloud_height_method = 0
  cloud_phase_infrared = 0.

	offset = .true.

#ifdef CT_1KM

  call read_int_array(mapi_filedata, pressure_sds, start, stride, edge, offset, cloud_top_pressure, status)
  call read_int_array(mapi_filedata, temperature_sds, start, stride, edge, offset, cloud_top_temperature, status)
  call read_int_array(mapi_filedata, sfc_temp_sds, start, stride, edge, offset, surface_temperature, status)   
  call read_byte_array(mapi_filedata, ctm_sds, start, stride, edge, cloud_height_method, status)
  call read_byte_array(mapi_filedata, cpi_sds, start, stride, edge, cloud_phase_infrared, status)


#else

  allocate(ctp_temp_array(edge(1), edge(2)))
  allocate(ctt_temp_array(edge(1), edge(2)))
  allocate(ctm_temp_array(edge(1), edge(2)))
  allocate(sfc_temp_array(edge(1), edge(2)))
  allocate(cpi_temp_array(edge(1), edge(2)))

  call read_int_array(mapi_filedata, pressure_sds, start, stride, edge, offset, ctp_temp_array, status)
  call read_int_array(mapi_filedata, temperature_sds, start, stride, edge, offset, ctt_temp_array, status)
  call read_int_array(mapi_filedata, sfc_temp_sds, start, stride, edge, offset, sfc_temp_array, status)
  call read_byte_array(mapi_filedata, ctm_sds, start, stride, edge, ctm_temp_array, status)
  call read_byte_array(mapi_filedata, cpi_sds, start, stride, edge, cpi_temp_array, status)

	dimsizes(1) = size(cloud_top_pressure, 1)
	dimsizes(2) = size(cloud_top_pressure, 2)

  if (status /= success) then
    call local_message_handler('Problem reported, see earlier message/s',status,'get_cloud_top_properties') 
  endif 


  do i = 1, dimsizes(1)
    local_i = (i - 1) / 5 + 1
    if (local_i > edge(1) ) local_i = edge(1)
    do j = 1, dimsizes(2)
       local_j = (j - 1) / 5 + 1
	   if (local_j > edge(2)) local_j = edge(2)
       cloud_top_pressure(i,j) =  ctp_temp_array(local_i, local_j)
       cloud_top_temperature(i,j) = ctt_temp_array(local_i, local_j)
       cloud_height_method(i,j) = ctm_temp_array(local_i, local_j)
       surface_temperature(i,j) = sfc_temp_array(local_i, local_j)
       cloud_phase_infrared(i,j) = cpi_temp_array(local_i, local_j)
	   	   
    enddo
  enddo
  deallocate( ctt_temp_array, ctp_temp_array, ctm_temp_array, sfc_temp_array, cpi_temp_array )

#endif


 end subroutine get_cloud_top_properties


end module specific_ancillary
