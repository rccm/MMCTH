module ancillary_module

 implicit none
 private

 public :: set_ancillary, get_above_cloud_properties, given_P_get_T

contains

subroutine set_ancillary( &
                         gdas_name,                 &
						 gdas_name2, &
						 ozone_name, &
                         ncepice_name,              &
                         nise_name,                 &
                         mod06input_name,           &
                         mod06input_filedata,       &
                         surfacealbedo_lib_659,      &
                         surfacealbedo_lib_858,      &
                         surfacealbedo_lib_124,      &
                         surfacealbedo_lib_164,      &
                         surfacealbedo_lib_21,       &
                         ecosystem_data_name,        &
                         snowicealbedo_data_name,    &
						 emissivity_name, &
                         mod06_start, mod06_stride, mod06_edge, &
                         debug,                     &
                         status)

   use nonscience_parameters
   use names
   use GeneralAuxType
   use core_arrays
   use science_parameters

#ifndef SEV_PR06CT

   use mod06_run_settings
   use specific_other
   use specific_ancillary

#endif

   implicit none
   
   include "hdf.f90"
   include "dffunc.f90"
   
   logical,     intent(in)  :: debug
   integer,     intent(inout),dimension(:)  :: mod06_start, mod06_stride, mod06_edge,mod06input_filedata
   character(*),intent(in)  :: gdas_name, gdas_name2, ozone_name, ncepice_name,  &
                               nise_name,   &
                               mod06input_name,  &
                               surfacealbedo_lib_659,      &
                               surfacealbedo_lib_858,      &
                               surfacealbedo_lib_124,      &
                               surfacealbedo_lib_164,      &
                               surfacealbedo_lib_21,       &
                               ecosystem_data_name,        &
                               snowicealbedo_data_name, emissivity_name

   integer,     intent(out) :: status

   real,allocatable,dimension(:,:)   ::   icec
   integer                                   :: xsize, ysize, levels, i,j
   integer, allocatable, dimension(:)         :: mod06_edge_5km, mod06_stride_5km, mod06_start_5km

	character(len=200) :: pressure_sdsname, temperature_sdsname, ctm_sdsname, sfc_temp_sdsname, cpi_sdsname
	logical :: EXECUTE_STANDARD, REPLACE_ALBEDO
	integer :: dim_size


	integer :: ancillary_filedata(1)

   status = 0
   levels = model_levels

	dim_size = size(mod06_edge)
   allocate(mod06_edge_5km(dim_size))
   allocate(mod06_stride_5km(dim_size))
   allocate(mod06_start_5km(dim_size))

   xsize = size(latitude,1)
   ysize = size(latitude,2)
   allocate(icec(xsize, ysize))
   status = 0

 
	cloud_height_method = 0
	cloud_phase_infrared = 0
	cloud_top_pressure = fillvalue_real
	cloud_top_temperature = fillvalue_real
	irw_temperature = fillvalue_real
 
   call read_ancillary_grids(  &
                  latitude,      &
                  longitude,     &
                  gdas_name,     &
				  gdas_name2, &
				  ozone_name, &
                  ncepice_name,  &
                  column_ozone,         &
                  icec)

     
	 EXECUTE_STANDARD = .true.
	 REPLACE_ALBEDO = .false.
	 
	 mod06_start_5km = mod06_start
	 mod06_stride_5km = mod06_stride
	 mod06_edge_5km = mod06_edge
	 
#ifndef SEV_PR06CT
		
	 
	 call get_specific_ancillary(mod06input_filedata,   &
                                mod06input_name,       &
								pressure_sdsname, &
								temperature_sdsname, &
								ctm_sdsname, &
								sfc_temp_sdsname, &
								cpi_sdsname, &
                                mod06_start_5km,       &
                                mod06_stride_5km,      &
                                mod06_edge_5km, EXECUTE_STANDARD, REPLACE_ALBEDO)
	 
	call get_cloud_top_properties(mod06input_filedata,   &
                                mod06input_name,       &
								pressure_sdsname, &
								temperature_sdsname, &
								ctm_sdsname, &
								sfc_temp_sdsname, &
								cpi_sdsname, &
                                mod06_start_5km,       &
                                mod06_stride_5km,      &
                                mod06_edge_5km,        &
								mod06_start, &
								mod06_edge,            &
                                cloud_top_pressure,    &
                                cloud_top_temperature, &
								cloud_height_method, &
								surface_temperature, &
								cloud_phase_infrared, &
								EXECUTE_STANDARD, &
                                status)
 
		
#ifndef SEV_PR06OD
  call get_surface_albedo(latitude, &
                          longitude,&
                          nise_name, &
                          surfacealbedo_lib_659,      &
                          surfacealbedo_lib_858,      &
                          surfacealbedo_lib_124,      &
                          surfacealbedo_lib_164,      &
                          surfacealbedo_lib_21,       &
                          ecosystem_data_name,        &
                          snowicealbedo_data_name,    &
						  emissivity_name, &
                          debug, &
                          icec, &
                          surface_albedo, &
						  surface_emissivity_land, &
                          status)
#else


  call get_seviri_albedo(latitude, &
                          longitude,&
						  mod06_start, mod06_stride, mod06_edge, &
                          nise_name, &
                          surfacealbedo_lib_659,      &
                          ecosystem_data_name,        &
                          snowicealbedo_data_name,    &
						  emissivity_name, &
                          debug, &
                          icec, &
                          surface_albedo, &
						  surface_emissivity_land, &
                          status)


#endif

	if (REPLACE_ALBEDO) then 
		call set_albedo
	endif

#endif
   

  deallocate( icec)
  deallocate( mod06_start_5km, mod06_stride_5km, mod06_edge_5km)
  
end subroutine set_ancillary

 subroutine get_bilinear_idx(lat, lon, i, j, idx_x, idx_y)
 
	real, intent(in) :: lat, lon
	integer, intent(in ) :: i,j
	integer, intent(inout) :: idx_x(2), idx_y(2)
	
	 real :: x,y, x0, dx, y0, dy
	
      x = min( max( lon,  -179.99 ), 179.99 )
      if( x > -999. .and. x < 0.0 ) x = lon+ 360.0

      y = min( max( lat, -89.99 ), 89.99 )
      y0 = 90.0
      dy = -1.0

		if ( x < (i*1.) ) then !+0.5) then 
			idx_x(1) = i-1
			idx_x(2) = i
		else
			idx_x(1) = i
			idx_x(2) = i+1
		endif

		if (idx_x(1) < 0 .or. idx_x(1) > 359) idx_x(1) = 0
		if (idx_x(2) < 0 .or. idx_x(2) > 359) idx_x(2) = 0
		
		if ( (y-y0)/dy < (j*1.) ) then !+0.5) then 
			idx_y(1) = j-1
			idx_y(2) = j
		else
			idx_y(1) = j
			idx_y(2) = j+1
		endif

		if (idx_y(1) < 0 ) idx_y(1) = 0
		if (idx_y(1) > 180) idx_y(1) = 180

		if (idx_y(2) < 0 ) idx_y(2) = 0
		if (idx_y(2) > 180) idx_y(2) = 180


		idx_x = idx_x + 1
		idx_y = idx_y + 1
 
 end subroutine get_bilinear_idx


real function get_W(RH, T, P)

! After Rogers and Yau. 

	real*8, intent(in) :: RH, T, P

	real, parameter :: c = 4187.
	real, parameter :: cpv = 1870.
	real, parameter :: L0 = 2.501e6
	real, parameter :: T0 = 273.
	
	real :: L
	real :: esat
	real :: ws
	real :: W
	
	if (T < 235.) then 
		L = 2.83e6   ! latent heat of ice, doesn't change very much
	else
		L = L0 - (c - cpv) * (T - T0) ! latent heat of water
	endif
	
	esat = 6.11657*exp((L/461.51)*(1/273. - 1/T) ) ! saturation vapor pressure in mb

	ws = esat * 0.622 / (P - esat) ! saturation mixing ratio in kg/kg

	W = RH / 100. * ws ! mixing ratio in kg/kg
	
	get_W = W * 1000. ! mixing ratio in g/kg
	


end function get_W

real function get_W_from_Td( Td, T, P)

	real*8, intent(in) :: Td, T, P

	real, parameter :: c = 4187.
	real, parameter :: cpv = 1870.
	real, parameter :: L0 = 2.501e6
	real, parameter :: T0 = 273.
	
	real :: L
	real :: esat, esatW, RH
	real :: ws
	real :: W
	
	if (T < 235.) then 
		L = 2.83e6   ! latent heat of ice, doesn't change very much
	else
		L = L0 - (c - cpv) * (T - T0) ! latent heat of water
	endif
	
	esat = 6.11657*exp((L/461.51)*(1/273. - 1/T) ) ! saturation vapor pressure in mb
	esatW = 6.11657*exp((L/461.51)*(1/273. - 1/Td) ) ! saturation vapor pressure in mb

	RH = esat / esatW * 100.

	ws = esat * 0.622 / (P - esat) ! saturation mixing ratio in kg/kg

	W = RH / 100. * ws ! mixing ratio in kg/kg
	
	get_W_from_Td = W * 1000. ! mixing ratio in g/kg

end function get_W_from_Td

  subroutine given_P_get_T(P, model_point, T)

	use global_model_grids
	
	real, intent(in)  ::  P
	real, intent(inout) :: T
	type(ancillary_type) :: model_point
    
	integer :: k, lev_start
	real :: factor
   
   	lev_start = model_point%trop_level-2
   
  	do k=lev_start, model_point%surface_level
  	  if ( P < model_point%pressure_profile(k)) &
  	    exit
  	end do
    
  	if (k <= lev_start) then 
  	  T = model_point%temp_profile(lev_start)
  	elseif (k >=  model_point%surface_level) then 
  	  T = model_point%temp_profile(model_point%surface_level)
  	else 
		T = model_point%temp_profile(k)
	endif

  end subroutine given_P_get_T



integer function find_trop(temp_prof, p_prof, sfc_lev)

	real, dimension(:), intent(in) :: temp_prof, p_prof
	integer*1, intent(in) :: sfc_lev

	real :: xmin
	integer :: ilev, imin
	real, parameter :: PTOP = 100.
	
	xmin = 999999.
	imin = 1
	
	do ilev = 1, sfc_lev-5
		if (temp_prof(ilev) < xmin .and. p_prof(ilev) > PTOP) then
			xmin = temp_prof(ilev)
			imin = ilev
		endif
	end do


! don't allow trop height > 400 mb (level 71)
	find_trop = -99
	do ilev = imin, 71
	
		if (temp_prof(ilev-1) >= temp_prof(ilev) .and. &
			temp_prof(ilev+1) > temp_prof(ilev)) then 
				find_trop = ilev
				exit
		endif
	end do

	if (find_trop == -99) find_trop = imin

end function find_trop





subroutine read_ancillary_grids(  lat, lon, ncepgdas_name,ncepgdas_name2, ozone_name, ice_name,  &
                                 ozone,  icec )
  use GeneralAuxType
  use nonscience_parameters
  use global_model_grids
  use core_arrays, only : model_info, cloudmask
  use science_parameters
  use modis_numerical_module, only: linearinterpolation, bilinear_interpolation
  use names, only: MYTIME, MYMONTH, Atmp_dirname, MYYEAR, MYDAY
  implicit none

	include "hdf.f90"
	include "dffunc.f90"

!  include 'Atmos_AncData.f90.inc'

  !-----------------------------------------------------------------------
  ! !f90
  !
  ! !description:
  !      retrieve ancillary data items for a given set of 
  !      latitude and longitude pairs.
  !
  ! !input parameters:
  !      lat       latitude (degrees, -90s to +90.0n)
  !      lon       longitude (degrees, -180w to +180e, greenwich=0)
  !
  ! !output parameters:
  !      pres      array of pressure levels (hpa)
  !      temp      array of atmospheric temperatures (k) at pres(0:15)
  !      mixr      array of water vapor mixing ratios (g/kg) at pres(0:15)
  !      land      land mask (0=water, 1=land)
  !      sfctmp    surface temperature (k)
  !      prmsl     pressure (hpa) at mean sea level
  !      pwat      precipitable water (g/cm**2)
  !      ugrd      surface wind u component (m/s)
  !      vgrd      surface wind v component (m/s)
  !      Ozone     total column Ozone (dobson units)
  !      icec      ice concentration (fraction)
  !
  ! !notes
  !      Modified from the MOD_PR06OD Collection 4 ReadNCEP.f algorithm written
  !      by Liam Gumley and Jason Li.

  real(single), intent(in)     :: lat(:,:),lon(:,:)
  character(*), intent(in)     :: ncepgdas_name, ncepgdas_name2, ozone_name, ice_name

  real(single), intent(inout)     :: icec(:,:), ozone(:,:)

  ! ... parameter definitions
  integer     npoints_x,        npoints_y  
  parameter ( npoints_x = 360,  npoints_y = 180 )

  real        missing
  parameter ( missing = -999.0 )

  ! ... declaration of local variables and  arrays
  character*8    esdt_name_ncep, esdt_name_ozone, esdt_name
  character*7    esdt_name_ice
  character*160  errmsg
  character*13   ncepgdastemp_name, ozonetemp_name
  character*12   icetemp_name
  character*400  temp_output_name, temp_input_name


  integer header( 0:7 ), i, ios, j, k, level, lun, pcfnum, reclen, status, ii, jj

  integer data_xsize, data_ysize, grid_i, grid_j, kstart

  real(single)  ::    satmix, x, x0, xlon, dx, y, y0, dy

        real :: yy2(2), yy(2), xx(2)
		integer :: idx_x(2), idx_y(2)

#ifdef USE_GDAS
	real :: met_grid( 0:359, 0:180, 0:53 ), met_grid2( 0:359, 0:180, 0:53 )
	integer, parameter :: num_gdas_vars = 54
#else
	real :: met_grid( 0:359, 0:180, 0:44 ), met_grid2( 0:359, 0:180, 0:44 )
	integer, parameter :: num_gdas_vars = 45
#endif

	integer met_year, met_month, met_day, met_hour, &
          ice_year, ice_month, ice_day, ice_hour, &
          ozn_year, ozn_month, ozn_day, ozn_hour

	integer met_date(4), met_date2(4), ice_date(4),  ozn_date(4)

	logical met_success, ice_success, ozn_success


	real :: Ts2m, W2m, Pmsl, Ts

  ! ... declaration of external functions
  integer modis_grib_driver 
  external modis_grib_driver 

	integer make_profile_101
	external make_profile_101
	
	integer profile_to_101
	external profile_to_101

	integer height_profile
	external height_profile

	real*8 :: model_lat
	real*8 :: heights(model_levels)

#ifdef USE_GDAS
	integer, parameter :: model_coarse = 27
	integer, parameter :: wind_start = 47
	integer, parameter :: rh_start = 6
#else
	integer, parameter :: model_coarse = 21
	integer, parameter :: wind_start = 40	
	integer, parameter :: rh_start = 1
#endif

	real :: a_factor
	integer*1 :: sfc_level
	real :: frac_time

    real*8 :: p(model_coarse), p_source(model_coarse)

	real*8 :: mixr(model_coarse), temp(model_coarse)
	
	real*8 :: pressures(model_levels), temp_hires(model_levels), mixr_hires(model_levels)
	real*8 :: pressures_source(model_levels)

	real :: Ts_forint(2,2), lat_set(2), lon_set(2)

	integer :: file_id, var_id, err_code, start(2), stride(2), edge(2)

	character*30 :: name_tag
	character*4 :: ttag
	character*3 :: dtag

#ifdef USE_GDAS
	p_source = (/ 10., 20., 30., 50., 70., 100., 150., 200., 250., 300., 350., 400., 450., 500., &
			550., 600., 650., 700., 750., 800., 850., 900., 925., 950., 975., 1000., 1100. /)
#else
	p_source = (/ 10., 20., 30., 50., 70., 100., 150., 200., 250., 300.,       400.,       500., &
			      600.,       700.,       800., 850., 900., 925., 950.,       1000., 1100. /)

#endif

  !-----------------------------------------------------------------------
  !      begin executable code 
  !-----------------------------------------------------------------------
  ! ... read input data files if this is the first call

	lun = 222


  data_xsize = size(lat,1)
  data_ysize = size(lat,2)

  if(.not. grids_are_read) then
	grids_are_read = .true. 
	! ..... set data ingest success/fail flags
	met_success  = .false.
	ice_success  = .false.
	ozn_success  = .false.

	do i = 1, 4
		met_date( i )  = int( missing )
		ice_date( i )  = int( missing )
		ozn_date( i )  = int( missing )
	end do


  !-----------------------------------------------------------------------
  !      get ncep meteorological data
  !-----------------------------------------------------------------------
  ! ...   unpack grib met file and write to binary file
	errmsg    = ' '
	ESDT_name = 'GDAS_0ZF'
  
	temp_input_name = trim(ncepgdas_name) // char(0)
#ifdef SEVIRI_INST

	if (MYDAY < 10) then 
		write(dtag, '(a2i1)') "00", MYDAY
	else if (MYDAY >= 10 .and. MYDAY < 100) then
		write(dtag, '(a1i2)') "0", MYDAY
	else
		write(dtag, '(i3)') MYDAY
	endif
	if (MYTIME < 10) then 
		write(ttag, '(a3i1)') "000", MYTIME
	else if (MYTIME >= 10 .and. MYTIME < 100) then
		write(ttag, '(a2i2)') "00", MYTIME
	else if (MYTIME >= 100 .and. MYTIME < 1000) then
		write(ttag, '(a1i3)') "0", MYTIME
	else 
		write(ttag, '(i4)') MYTIME
		
	endif
	write(name_tag, '(i4)') MYYEAR
	name_tag = trim(name_tag) // trim(dtag) // trim(ttag)
	
	temp_output_name = trim(Atmp_dirname) //  "/" // trim(name_tag) // "_1.bin" // char(0)
#else
	temp_output_name = trim(ncepgdas_name) // ".bin" // char(0)
#endif  


	status    = modis_grib_driver( temp_input_name,  &
                                 ESDT_name, errmsg, &
                                 met_year, met_month, &
                                 met_day, met_hour )



	if ( status /= success ) then
		status = failure
		call local_message_handler ('error reported from grib driver, ncep read, Check PCF file', &
                                status, &
                                'read_ancillary_grids')

    ! ....  open unpacked met file
	else

		reclen = 360*181*num_gdas_vars*4

		open (unit = lun, file = temp_output_name, form = 'unformatted', &
	       access = 'direct', status = 'old', &
             recl = reclen, iostat = ios )

		if ( status /= 0 ) then
			status = 2
			call local_message_handler ('error reported from openr, ncep read, Check PCF file', &
                                status, &
                                'read_ancillary_grids')
		else 
		!  read the unpacked met file
			if ( status == 0 ) then
				ios= 0
				read( lun, rec = 1, iostat = ios ) met_grid


				if ( ios /= 0 ) then
					level = 2
				else
					met_success   = .true.
					met_date( 1 ) = met_year
					met_date( 2 ) = met_month
					met_date( 3 ) = met_day
					met_date( 4 ) = met_hour
					
				endif

				close(lun)
			endif
			
		endif
		
	endif

	errmsg    = ' '
	ESDT_name = 'GDAS_0ZF'
  
	temp_input_name = trim(ncepgdas_name2) // char(0)
#ifdef SEVIRI_INST
	temp_output_name = trim(Atmp_dirname) //  "/" // trim(name_tag) // "_2.bin" // char(0)
#else
	temp_output_name = trim(ncepgdas_name2) // ".bin" // char(0)
#endif  
  
	status    = modis_grib_driver( temp_input_name, &
                                 ESDT_name, errmsg, &
                                 met_year, met_month, &
                                 met_day, met_hour )
 
	if ( status /= success ) then
		status = failure
		call local_message_handler ('error reported from grib driver, ncep read, Check PCF file', &
                                status, &
                                'read_ancillary_grids')

    ! ....  open unpacked met file
	else

		reclen = 360*181*num_gdas_vars*4
		open (unit = lun, file = temp_output_name, form = 'unformatted', &
	       access = 'direct', status = 'old', &
             recl = reclen, iostat = ios )

		if ( status /= 0 ) then
			status = 2
			call local_message_handler ('error reported from openr, ncep read, Check PCF file', &
                                status, &
                                'read_ancillary_grids')
		else 
		!  read the unpacked met file
			if ( status == 0 ) then
				ios= 0
				read( lun, rec = 1, iostat = ios ) met_grid2

				if ( ios /= 0 ) then
					level = 2
				else
					met_success   = .true.
					met_date2( 1 ) = met_year
					met_date2( 2 ) = met_month
					met_date2( 3 ) = met_day
					met_date2( 4 ) = met_hour 
					if (met_date2(4) == 0) met_date2(4) = 24
				endif

				close(lun)
			endif
			
		endif
		
	endif

! capture month of the granule
	MYMONTH = met_date(2)
	
#ifdef USE_TOAST
!-----------------------------------------------------------------------
!     Get TOAST ozone concentration
!-----------------------------------------------------------------------
! ...   Unpack grib TOAST file and write to binary file
	errmsg    = ' '
	ESDT_name = 'OZ_DAILY'
  
	temp_input_name = trim(ozone_name) // char(0)
#ifdef SEVIRI_INST
	temp_output_name = trim(Atmp_dirname) //  "/" // trim(name_tag) // "_o3.bin" // char(0)
#else
	temp_output_name = trim(ozone_name) // ".bin" // char(0)
#endif  
  
	status    = modis_grib_driver( temp_input_name,  &
                                 ESDT_name, errmsg, &
                                 ozn_year, ozn_month, &
                                 ozn_day, ozn_hour )
	if ( status /= success ) then
		status = failure
		call local_message_handler ('error reported from grib driver, TOAST read, Check PCF file', &
                                status, &
                                'read_ancillary_grids')
	else

!   Open unpacked TOAST file
		reclen = 360*181*4
	
		open (unit = lun, file = temp_output_name, form = 'unformatted', &
             access = 'direct', status = 'old', &
             recl = reclen, iostat = ios )
	
	

		if ( status /= success ) then
			status = failure
			call local_message_handler ('error reported opening temp ozone file, Check PCF file', &
                                  status, &
                                  'read_ancillary_grids')
		else
!     Read the unpacked ozone file
			if ( status == 0 ) then
				ios  = 0
				read( lun, rec = 1, iostat = ios )ozn_grid
				if ( ios /= success ) then
					status = success
					call local_message_handler ('error reported reading temp sea ice file, Check PCF file', &
                                status, &
                                'read_ancillary_grids')
				else
					ozn_success   = .true.
					ozn_date( 1 ) = ozn_year
					ozn_date( 2 ) = ozn_month
					ozn_date( 3 ) = ozn_day
					ozn_date( 4 ) = ozn_hour
				endif
				close(lun)
			endif
			
		endif
		
	endif	

#endif	

#ifndef SEV_PR06CT
!-----------------------------------------------------------------------
!     Get SSM/I sea ice concentration
!-----------------------------------------------------------------------
! ...   Unpack grib sea ice file and write to binary file
	errmsg    = ' '
	ESDT_name = 'SEA_ICE'
  
	temp_input_name = trim(ice_name) // char(0)
#ifdef SEVIRI_INST
	temp_output_name = trim(Atmp_dirname) //  "/" // trim(name_tag) // "_ice.bin" // char(0)
#else
	temp_output_name = trim(ice_name) // ".bin" // char(0)
#endif  
  
	status    = modis_grib_driver( temp_input_name,  &
                                 ESDT_name, errmsg, &
                                 ice_year, ice_month, &
                                 ice_day, ice_hour )
	if ( status /= success ) then
		status = failure
		call local_message_handler ('error reported from grib driver, NCEP sea ice read, Check PCF file', &
                                status, &
                                'read_ancillary_grids')
	else

!   Open unpacked ice file
		reclen = 720*360*4
	
		open (unit = lun, file = temp_output_name, form = 'unformatted', &
             access = 'direct', status = 'old', &
             recl = reclen, iostat = ios )
	
	

		if ( status /= success ) then
			status = failure
			call local_message_handler ('error reported opening temp sea ice, Check PCF file', &
                                  status, &
                                  'read_ancillary_grids')
		else
!     Read the unpacked ice file
			if ( status == 0 ) then
				ios  = 0
				read( lun, rec = 1, iostat = ios )ice_grid
				if ( ios /= success ) then
					status = success
					call local_message_handler ('error reported reading temp sea ice file, Check PCF file', &
                                status, &
                                'read_ancillary_grids')
				else
					ice_success   = .true.
					ice_date( 1 ) = ice_year
					ice_date( 2 ) = ice_month
					ice_date( 3 ) = ice_day
					ice_date( 4 ) = ice_hour
				endif
				close(lun)
			endif
			
		endif
		
	endif	
#else
! no use for this in SEV_PR06CT
	ice_grid = 0.0
#endif

        xx(1) = met_date(4) * 1.0
        xx(2) = met_date2(4) * 1.0

		frac_time = MYTIME / 100 + mod(MYTIME, 100) / 60.
				
	! calculate the pressure levels for the 101-level profile, courtesy UW-Madison
    status = make_profile_101(pressures_source)

	do i=1, grid_xsize
		ii = i-1

		do j=1, grid_ysize
	
			jj = j-1

			p = p_source
			pressures = pressures_source
    
			model_lat = 89.5 - jj*1.0
			
			yy(1) = met_grid(ii,jj,wind_start)
			yy(2) = met_grid2(ii,jj, wind_start) 

			yy2(1) = met_grid(ii,jj,wind_start+1)
			yy2(2) = met_grid2(ii,jj,wind_start+1)       

			model_info(i,j)%wind_speed = sqrt( linearinterpolation( xx, yy, frac_time) **2 + &
									linearinterpolation( xx, yy2, frac_time) **2 )
									
									
			kstart = 0
			do k = 1, model_coarse-1
				temp(k) = linearinterpolation( xx, &
									(/ met_grid(ii,jj,kstart), met_grid2(ii,jj, kstart) /), frac_time) !met_grid( i, j, k )
				kstart = kstart + 1
			end do

			kstart = model_coarse-1
			do k = rh_start, model_coarse-1
				mixr(k) = linearinterpolation( xx, &
						(/ met_grid(ii,jj, kstart), met_grid2(ii,jj, kstart) /), frac_time)  !met_grid( i, j, k + 16 )
				if (mixr(k) < 0.) mixr(k) = 0.
				kstart = kstart + 1 
			end do

!     convert relative humidity profile (%) to mixing ratio (g/kg)

			do k = rh_start, model_coarse-1
				mixr(k) = get_W (mixr(k), temp(k), p(k))
			end do



!     extrapolate mixing ratio profile from 100 hPa to 10 hPa
#ifdef USE_GDAS
			do k = 1, 5
				mixr(k) = max( mixr(6), 0.003d0 ) * ( p( k ) / 100.0 )**3 ! was 21 
				mixr(k) = max( mixr(k), 0.003d0 )
			end do
#endif


            model_info(i,j)%Ps = linearinterpolation( xx, (/ met_grid(ii,jj,wind_start+2), &
													met_grid2(ii,jj, wind_start+2) /), frac_time) * 0.01 !met_grid( i, j, 76 )
			
! get the surface parameters and convert the RH at 2m to mixing ratio
! this Ts2m is really Ts, the surface temperature. Wisconsin switched to using TSFC instead of T2M
			Ts2m = linearinterpolation( xx, (/ met_grid(ii,jj,wind_start+3), &
											met_grid2(ii,jj, wind_start+3) /), frac_time)
			
! if we have ECMWF this W2m is not an RH, but a dew point temperature.
			W2m = linearinterpolation( xx, (/ met_grid(ii,jj,wind_start+4),&
											met_grid2(ii,jj, wind_start+4) /), frac_time)

#ifdef USE_GDAS
			W2m = get_W(W2m*1.0d0, Ts2m*1.0d0, model_info(i,j)%Ps*1.0d0)
			Pmsl = linearinterpolation( xx, (/ met_grid(ii,jj,52), met_grid2(ii,jj, 52) /), frac_time)* 0.01
#else
			W2m = get_W_from_Td(W2m*1.0d0, Ts2m*1.0d0, model_info(i,j)%Ps*1.0d0)
			Pmsl = model_info(i,j)%Ps
#endif

			if (model_info(i,j)%Ps > 0. .and. p(model_coarse-1) <= model_info(i,j)%Ps) then 
				model_info(i,j)%surface_level = model_coarse
			else
				model_info(i,j)%surface_level = 0
			endif

			model_info(i,j)%Ts = Ts2m
#ifdef USE_GDAS
			model_info(i,j)%col_o3 = linearinterpolation( xx, (/ met_grid(ii,jj,53), met_grid2(ii,jj, 53) /), frac_time)
#else
			model_info(i,j)%col_o3 = 0.0
#endif
			
			
			kstart = model_coarse / 2
			do k=kstart, model_coarse
			
				if (model_info(i,j)%Ps > 0. .and. p(k) > model_info(i,j)%Ps) then 
					if (model_info(i,j)%surface_level == 0) then
						
						model_info(i,j)%surface_level = k
						temp(k) = Ts2m
						if ( (model_info(i,j)%Ps - p(k-1)) < 5. .or. &
							 (p(k) - model_info(i,j)%Ps) < 5. ) then 
							p(k) = (p(k) + p(k-1)) / 2.
						else
							p(k) = model_info(i,j)%Ps
						endif
						mixr(k) = W2m
					
					else
					
						temp(k) = Ts2m
						mixr(k) = W2m
						p(k) = p_source(k-1)
					
					endif
				
				endif

			end do

! add the surface level into the coarse profile
		   if (model_info(i,j)%surface_level /= model_coarse) then 
				p(model_coarse) = p_source(model_coarse)
			else
				p(model_coarse) = model_info(i,j)%Ps
			endif
			temp(model_coarse) = Ts2m
			mixr(model_coarse) = W2m



! now we determine the lowest valid level of the new high-res profile
			kstart = model_levels / 2
			do k = kstart, model_levels
				if (pressures(k) >= model_info(i,j)%Ps) then
					model_info(i,j)%surface_level = k
					exit
				endif
			end do


	

			
! interpolate the profile to 101 levels, courtesy UW-Madison

			status = profile_to_101(p, temp, mixr, model_coarse, model_lat, pressures, temp_hires, mixr_hires)



			sfc_level = model_info(i,j)%surface_level
			a_factor = (model_info(i,j)%Ps - pressures(sfc_level-1)) / &
						( pressures(sfc_level) - pressures(sfc_level-1))
			temp_hires(sfc_level) = temp_hires(sfc_level-1) + a_factor * (temp_hires(sfc_level) - temp_hires(sfc_level-1))
			mixr_hires(sfc_level) = mixr_hires(sfc_level-1) + a_factor * (mixr_hires(sfc_level) - mixr_hires(sfc_level-1))
			
			pressures(sfc_level) = model_info(i,j)%Ps

			model_info(i,j)%temp_profile = temp_hires
			model_info(i,j)%mixr_profile = mixr_hires
							
! calculate the height profile, courtesy UW-Madison			

			status = height_profile(pressures, temp_hires, mixr_hires, &
								heights, model_levels, Pmsl*1.0d0)


			model_info(i,j)%height_profile = heights	
			model_info(i,j)%pressure_profile = pressures				


	
			model_info(i,j)%trop_level = find_trop (model_info(i,j)%temp_profile, &
											model_info(i,j)%pressure_profile, sfc_level)

			
!			if (i==1 .and. j==91) then 
!				do k=1, model_levels
!					print*, k, model_info(i,j)%pressure_profile(k), model_info(i,j)%mixr_profile(k), model_info(i,j)%temp_profile(k)
!				end do
!			end if
	

			
		end do
	end do



 endif 
	


  
! get lat/long data from grids

  do grid_i = 1, data_xsize
    do grid_j = 1, data_ysize

! don't waste time processing and setting ancillary if there is no cloud. Why bother?
      if (lon(grid_i,grid_j) <= -999. .or. lat(grid_i,grid_j) <= -999.0 &
#ifdef SEVIRI_INST   
			.or. cloudmask(grid_i,grid_j)%probablyclear_95 .or. cloudmask(grid_i,grid_j)%probablyclear_99 ) then
#else
		) then 
#endif		
		icec(grid_i, grid_j)   = -999.
		ozone(grid_i, grid_j)  = -999.
		
		cycle
		
      endif
	
	  call get_model_idx(lat(grid_i,grid_j), lon(grid_i,grid_j), i, j)

	  i = i-1
	  j = j-1

	  call get_bilinear_idx(lat(grid_i,grid_j), lon(grid_i,grid_j), i, j, idx_x, idx_y)
	

                if (idx_x(1)-1 < 180) then
                        lon_set(1) = (idx_x(1) - 1) * 1.0 !+ 0.5
                else
                        lon_set(1) = (idx_x(1) - 1) * 1.0 - 360. !+ 0.5
                endif

                if (idx_x(2)-1 < 180) then
                        lon_set(2) = (idx_x(2) - 1) * 1.0 !+ 0.5
                else
                        lon_set(2) = (idx_x(2) - 1) * 1.0 - 360. !+ 0.5
                endif


                if (lon(grid_i, grid_j) > 179 ) then !.5) then
                        lon_set(1) = 179 !.5
                        lon_set(2) = 180 !.5
                endif

                if (lon(grid_i, grid_j) < -179 ) then !.5) then
                        lon_set(1) = -180 !.5
                        lon_set(2) = -179 !.5
                endif


                lat_set = 90-(idx_y-1)*1.0 ! - 0.5
                if (lat(grid_i, grid_j) > 89 ) then !.5) then
                        lat_set(1) = 90.0
                        lat_set(2) = 89 !.5
                endif

                if (lat(grid_i, grid_j) < -89 ) then !.5) then
                        lat_set(1) = -89. !.5
                        lat_set(2) = -90.
                endif

#ifdef USE_TOAST
		Ts_forint(1,1) = ozn_grid(idx_x(1)-1, idx_y(1)-1)
		Ts_forint(1,2) = ozn_grid(idx_x(1)-1, idx_y(2)-1)
		Ts_forint(2,1) = ozn_grid(idx_x(2)-1, idx_y(1)-1)
		Ts_forint(2,2) = ozn_grid(idx_x(2)-1, idx_y(2)-1)


#else
		Ts_forint(1,1) = model_info(idx_x(1), idx_y(1))%col_o3
		Ts_forint(1,2) = model_info(idx_x(1), idx_y(2))%col_o3
		Ts_forint(2,1) = model_info(idx_x(2), idx_y(1))%col_o3
		Ts_forint(2,2) = model_info(idx_x(2), idx_y(2))%col_o3

	
#endif		

      ozone(grid_i, grid_j) = bilinear_interpolation(  lon_set, &
												lat_set, lon(grid_i, grid_j), lat(grid_i, grid_j), Ts_forint, 1)

	  
!     compute cell coordinates in ice grid and save ice data
      x = min( max( lon(grid_i,grid_j), -179.99 ), 179.99 )
      if( x .lt. 0.0 ) x = lon(grid_i,grid_j) + 360.0
      x0 = 0.25
      dx = 0.5
      i = int( ( x - x0 + 0.5*dx ) / dx )
      if( i .eq. 720 ) i = 0

      y = min( max( lat(grid_i,grid_j), -89.99 ), 89.99 )
      y0 = 89.75
      dy = -0.5
      j = int( ( y - y0 + 0.5*dy ) / dy )

      icec(grid_i, grid_j) = ice_grid( i, j )

    enddo
  enddo
  
  
 end subroutine read_ancillary_grids

#ifndef SEVIRI_INST
 subroutine get_surface_albedo(latitude, &
                              longitude,&
                              nise_filename, &
                              surfacealbedo_lib_659,      &
                              surfacealbedo_lib_858,      &
                              surfacealbedo_lib_124,      &
                              surfacealbedo_lib_164,      &
                              surfacealbedo_lib_21,       &
                              IGBPfilename,        &
                              snowicealbedo_data_name,    &
							  emissivity_name, &
                              debug, &
                              icec, &
                              surface_albedo, surface_emissivity, &
                              status)
!-----------------------------------------------------------------------
!f90 modisretrieval
!
!Description:
!
! get surface albedos for all required bands
!
!input parameters:
!
!output parameters:
!
!revision history:
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
!
!design note:
!
!end
!----------------------------------------------------------------------

! collection 4 albedo determination
  use GeneralAuxType
  use nonscience_parameters
  use modis_albedo, only : getAlbedoEco
  use core_arrays, only: cloudsummary, cloudmask
  use names, only : MYDAY
  use mod06_run_settings
  implicit none

!  include 'PGS_PC.f'
!  include 'PGS_TD.f'
!  include 'PGS_SMF.f'
  
  real, dimension(:,:), intent(in)    :: latitude, longitude
  character(*),intent(in)             :: IGBPfilename,  &
                                         nise_filename, &
                                         surfacealbedo_lib_659,      &
                                         surfacealbedo_lib_858,      &
                                         surfacealbedo_lib_124,      &
                                         surfacealbedo_lib_164,      &
                                         surfacealbedo_lib_21,       &
                                         snowicealbedo_data_name, emissivity_name
  logical, intent(in)                 :: debug 
  real(single), dimension(:,:), intent(in) ::  icec
  real, dimension(:,:,:), intent(out) ::  surface_emissivity
  integer*2, dimension(:,:,:), intent(out) :: surface_albedo
  integer, intent(out)                :: status 

  integer   (kind = integer_fourbyte),  parameter       :: MaxFileNameLength  = 2025
  integer*1,       &
              allocatable, dimension(:,:)               :: ecosystem, sfc_info
  integer   (integer_fourbyte)                          :: julianday,i,j,NumAlbWavelengths

  character (len = MaxFileNameLength)                   :: AlbedoWavelengthFN(set_albedo_bands)
  character (len = MaxFileNameLength)                   :: SnowAlbedoFN
  character (len = 10)							        :: WavelengthText(set_albedo_bands)
  real      											:: Wavelength(set_albedo_bands)
  logical							                    :: ProcessLandAlbedo(set_albedo_bands)

  !Timestamp variables:
  integer (kind = integer_fourbyte), parameter          :: LUN_TimeStamp = 10258
  character (len = 40)                                  :: dateTime_a, dateTime_b
  integer (kind = integer_fourbyte)                     :: PGS_PC_GetConfigData, &
                                                           PGS_TD_asciitime_atob
  character (len = 3)                                   :: DayOfYearString

	integer :: lat_wid, lat_ht, num_iter
	integer :: lat_start, lat_end
	
	lat_wid = size(latitude, 1)
	lat_ht = size(latitude, 2)




  allocate(ecosystem(lat_wid, lat_ht), sfc_info(lat_wid, lat_ht))


  NumAlbWavelengths = set_albedo_bands

  !Define the wavelengths:
  WavelengthText(1) = "0.659"
  WavelengthText(2) = "0.858"
  WavelengthText(3) = "1.24"
  WavelengthText(4) = "1.64"
  WavelengthText(5) = "2.13"
  WavelengthText(6) = "3.7"

  !Read in the HDF filenames containing the specified wavelengths from PCF file:
  ! Note to mark, these are read in from PCF file.
  AlbedoWavelengthFN(1) = surfacealbedo_lib_659 
  AlbedoWavelengthFN(2) = surfacealbedo_lib_858
  AlbedoWavelengthFN(3) = surfacealbedo_lib_124
  AlbedoWavelengthFN(4) = surfacealbedo_lib_164
  AlbedoWavelengthFN(5) = surfacealbedo_lib_21
  AlbedoWavelengthFN(6) = "notUsed"
  SnowAlbedoFN = snowicealbedo_data_name


   ecosystem = -99

  !Determine the Julian Day by reading in the time stamp and converting:
  !Get the Time Stamp from the pcf file:
!  Status = PGS_PC_GetConfigData(LUN_TimeStamp,dateTime_a)
!  if (status /= PGS_S_SUCCESS) then
!    call local_message_handler('Problem Reading in the PCF Time Stamp',status,'get_surface_albedo') 
!  endif 
  !Convert to day of year string:
!  Status = PGS_TD_asciitime_atob(dateTime_a,dateTime_b)
!  if (status /= PGS_S_SUCCESS) then
!    call local_message_handler('Problem converting time stamp to DOY',status,'get_surface_albedo') 
!  endif 
  !Extract the Day of Year:
!  DayOfYearString = dateTime_b(6:8)
  !Convert to integer:
!  read(DayOfYearString, *) JulianDay
  !If not successful, set JulianDay to 1 so that it does not crash:
	JulianDay = MYDAY
	
  if (JulianDay < 1 .or. JulianDay > 366) JulianDay = 1

  ProcessLandAlbedo (1) = .true.
  ProcessLandAlbedo (2) = .true.
  ProcessLandAlbedo (3) = .true.
  ProcessLandAlbedo (4) = .true.
  ProcessLandAlbedo (5) = .true.
  ProcessLandAlbedo (6) = .true.


  num_iter = lat_ht / 100 + 1

  do i=1, num_iter
  
		lat_start = (i-1)*100 + 1
		if (lat_start > lat_ht) lat_start = lat_ht

		if (i==num_iter) then 
			lat_end = lat_ht
		else
			lat_end = lat_start + 99
		endif
 
		call getAlbedoEco (latitude(:, lat_start:lat_end),          &
						longitude(:, lat_start:lat_end),         &
						julianday,         &
						AlbedoWavelengthFN,&
						IGBPfilename,      &
						SnowAlbedoFN,      &
						nise_filename,     &
						emissivity_name, &
						Debug,             &
						WavelengthText,    &
						ProcessLandAlbedo, &
						icec(:, lat_start:lat_end),              &
						Wavelength,        &
						surface_albedo(:, lat_start:lat_end, :),    &
						surface_emissivity(:, lat_start:lat_end, :), &
						Ecosystem(:, lat_start:lat_end),         &
						Status, sfc_info(:, lat_start:lat_end)             )

	end do
	
	
	where (surface_albedo > 1000) surface_albedo = -9999

  do i = 1, lat_wid
   do j = 1, lat_ht
	

	 cloudmask(i,j)%ocean_no_snow = 0
	 cloudmask(i,j)%ocean_snow = 0
	 cloudmask(i,j)%land_no_snow = 0
	 cloudmask(i,j)%land_snow = 0
		
	 if (sfc_info(i,j) == 0) then 
	 	cloudmask(i,j)%ocean_no_snow = 1
	 else if (sfc_info(i,j) == 1) then 
	 	cloudmask(i,j)%ocean_snow = 1
	 else if (sfc_info(i,j) == 2) then 
	 	cloudmask(i,j)%land_no_snow = 1
	 else if (sfc_info(i,j) == 3) then 
	 	cloudmask(i,j)%land_snow = 1
	 endif

     if (Ecosystem(i,j) == 0 )then

       cloudsummary(i,j)%land_surface  = .false.
       cloudsummary(i,j)%ocean_surface = .true.
       cloudsummary(i,j)%coastal_surface = .false.
       cloudsummary(i,j)%desert_surface = .false.
       cloudsummary(i,j)%snowice_surface = .false.
     elseif (Ecosystem(i,j) == 15) then

       cloudsummary(i,j)%land_surface  = .false.

       cloudsummary(i,j)%ocean_surface = .false.
       cloudsummary(i,j)%coastal_surface = .false.
       cloudsummary(i,j)%desert_surface = .false.
       cloudsummary(i,j)%snowice_surface = .true.
     else 
       cloudsummary(i,j)%land_surface  = .true.
       cloudsummary(i,j)%ocean_surface = .false.
       cloudsummary(i,j)%coastal_surface = .false.
       cloudsummary(i,j)%desert_surface = .false.
       cloudsummary(i,j)%snowice_surface = .false.
     endif

   enddo
  enddo
  deallocate(ecosystem, sfc_info)
  if (status /= success) then
    call local_message_handler('Problem reported in get_surface_albedo, see earlier message/s',status,'get_surface_albedo') 
  endif 
 end subroutine get_surface_albedo
#else


 subroutine get_seviri_albedo(latitude, &
                              longitude,&
							  start, stride, edge, &
                              nise_filename, &
                              surfacealbedo_lib_659,      &
                              IGBPfilename,        &
                              snowicealbedo_data_name,    &
							  emissivity_name, &
                              debug, &
                              icec, &
                              surface_albedo, surface_emissivity, &
                              status)
!-----------------------------------------------------------------------
!f90 modisretrieval
!
!Description:
!
! get surface albedos for all required bands
!
!input parameters:
!
!output parameters:
!
!revision history:
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
!
!design note:
!
!end
!----------------------------------------------------------------------

! collection 4 albedo determination
  use GeneralAuxType
  use nonscience_parameters
  use seviri_albedo, only : getAlbedoEco
  use core_arrays, only: cloudsummary, cloudmask
  use names, only : MYDAY
  use mod06_run_settings
  implicit none

!  include 'PGS_PC.f'
!  include 'PGS_TD.f'
!  include 'PGS_SMF.f'
  
  real, dimension(:,:), intent(in)    :: latitude, longitude
  character(*),intent(in)             :: IGBPfilename,  &
                                         nise_filename, &
                                         surfacealbedo_lib_659,      &
                                         snowicealbedo_data_name, emissivity_name
  logical, intent(in)                 :: debug 
  real(single), dimension(:,:), intent(in) ::  icec
  integer, intent(in), dimension(:) :: start, stride, edge
  real, dimension(:,:,:), intent(out) ::  surface_emissivity
  real, dimension(:,:,:), intent(out) :: surface_albedo
  integer, intent(out)                :: status 

  integer   (kind = integer_fourbyte),  parameter       :: MaxFileNameLength  = 2025
  integer*1,       &
              allocatable, dimension(:,:)               :: ecosystem
  integer   (integer_fourbyte)                          :: julianday,i,j,NumAlbWavelengths
  character (len = MaxFileNameLength)                   :: AlbedoWavelengthFN(set_albedo_bands)
  character (len = MaxFileNameLength)                   :: SnowAlbedoFN
  character (len = 10)							        :: WavelengthText(set_albedo_bands)
  real      											:: Wavelength(set_albedo_bands)
  logical							                    :: ProcessLandAlbedo(set_albedo_bands)

  !Timestamp variables:
  integer (kind = integer_fourbyte), parameter          :: LUN_TimeStamp = 10258
  character (len = 40)                                  :: dateTime_a, dateTime_b
  integer (kind = integer_fourbyte)                     :: PGS_PC_GetConfigData, &
                                                           PGS_TD_asciitime_atob
  character (len = 3)                                   :: DayOfYearString

	integer :: lat_wid, lat_ht, num_iter
	integer :: lat_start, lat_end
	
	lat_wid = size(latitude, 1)
	lat_ht = size(latitude, 2)




  allocate(ecosystem(lat_wid, lat_ht))


  NumAlbWavelengths = set_albedo_bands

  !Define the wavelengths:
  WavelengthText(1) = "0.659"
  WavelengthText(2) = "0.858"
  WavelengthText(3) = "1.24"
  WavelengthText(4) = "1.64"
  WavelengthText(5) = "2.13"
  WavelengthText(6) = "3.7"

  !Read in the HDF filenames containing the specified wavelengths from PCF file:
  ! Note to mark, these are read in from PCF file.
  ! for SEVIRI all albedo channels are contained in one file
  AlbedoWavelengthFN(1) = surfacealbedo_lib_659 
  AlbedoWavelengthFN(2) = surfacealbedo_lib_659
  AlbedoWavelengthFN(3) = "notUsed"
  AlbedoWavelengthFN(4) = surfacealbedo_lib_659
  AlbedoWavelengthFN(5) = "notUsed"
  AlbedoWavelengthFN(6) = "notUsed"
  SnowAlbedoFN = snowicealbedo_data_name


   ecosystem = -99

  !Determine the Julian Day by reading in the time stamp and converting:
  !Get the Time Stamp from the pcf file:
!  Status = PGS_PC_GetConfigData(LUN_TimeStamp,dateTime_a)
!  if (status /= PGS_S_SUCCESS) then
!    call local_message_handler('Problem Reading in the PCF Time Stamp',status,'get_surface_albedo') 
!  endif 
  !Convert to day of year string:
!  Status = PGS_TD_asciitime_atob(dateTime_a,dateTime_b)
!  if (status /= PGS_S_SUCCESS) then
!    call local_message_handler('Problem converting time stamp to DOY',status,'get_surface_albedo') 
!  endif 
  !Extract the Day of Year:
!  DayOfYearString = dateTime_b(6:8)
  !Convert to integer:
!  read(DayOfYearString, *) JulianDay
  !If not successful, set JulianDay to 1 so that it does not crash:
	JulianDay = MYDAY
	
  if (JulianDay < 1 .or. JulianDay > 366) JulianDay = 1

  ProcessLandAlbedo (1) = .true.
  ProcessLandAlbedo (2) = .true.
  ProcessLandAlbedo (3) = .false.
  ProcessLandAlbedo (4) = .true.
  ProcessLandAlbedo (5) = .false.
  ProcessLandAlbedo (6) = .false.

 
		call getAlbedoEco (latitude,          &
						longitude,         &
						julianday,         &
						start, stride, edge, &
						AlbedoWavelengthFN,&
						IGBPfilename,      &
						SnowAlbedoFN,      &
						nise_filename,     &
						emissivity_name, &
						Debug,             &
						WavelengthText,    &
						ProcessLandAlbedo, &
						icec,              &
						Wavelength,        &
						surface_albedo,    &
						surface_emissivity, &
						Ecosystem,         &
						Status             )


  do i = 1, lat_wid
   do j = 1, lat_ht
   
       cloudsummary(i,j)%land_surface  = .false.
       cloudsummary(i,j)%ocean_surface = .false.
       cloudsummary(i,j)%coastal_surface = .false.
       cloudsummary(i,j)%desert_surface = .false.
       cloudsummary(i,j)%snowice_surface = .false.
	
     if (Ecosystem(i,j) == 0 )then

       cloudsummary(i,j)%land_surface  = .false.
       cloudsummary(i,j)%ocean_surface = .true.
       cloudsummary(i,j)%coastal_surface = .false.
       cloudsummary(i,j)%desert_surface = .false.
       cloudsummary(i,j)%snowice_surface = .false.
     elseif (Ecosystem(i,j) == 15) then

       cloudsummary(i,j)%land_surface  = .false.

       cloudsummary(i,j)%ocean_surface = .false.
       cloudsummary(i,j)%coastal_surface = .false.
       cloudsummary(i,j)%desert_surface = .false.
       cloudsummary(i,j)%snowice_surface = .true.
     else if (Ecosystem(i,j) > 0 .and. Ecosystem(i,j) /= 15) then  
       cloudsummary(i,j)%land_surface  = .true.
       cloudsummary(i,j)%ocean_surface = .false.
       cloudsummary(i,j)%coastal_surface = .false.
       cloudsummary(i,j)%desert_surface = .false.
       cloudsummary(i,j)%snowice_surface = .false.
     endif

   enddo
  enddo
  deallocate(ecosystem)
  if (status /= success) then
    call local_message_handler('Problem reported in get_surface_albedo, see earlier message/s',status,'get_surface_albedo') 
  endif 
 end subroutine get_seviri_albedo





#endif


 subroutine get_above_cloud_properties(  pprof, wprof, sfc_lev,  &
                                   cloud_top_pressure,   &
                                   abovecloud_watervapor, &
                                   status)

   use nonscience_parameters

	implicit none

  real, intent(in), dimension(:) :: pprof, wprof
  integer*1, intent(in) :: sfc_lev
  real, intent(in)   :: cloud_top_pressure
  real, intent(out)   :: abovecloud_watervapor
  integer, intent(out)                       :: status

  integer    :: level_index, k, num_levels
  real       :: pw_layer, total_pw, slope, dummy
  real       ::  watervapor_lower, watervapor_upper

  abovecloud_watervapor = fillvalue_real

  if (cloud_top_pressure /= fillvalue_real) then    

	num_levels = sfc_lev + 1
	k = 2
	total_pw = 0.
	do while ( pprof(k) < cloud_top_pressure .and. k < num_levels)
		pw_layer = ( pprof(k) - pprof(k-1) ) *  &
             ( wprof(k-1) + wprof(k) ) * 0.5/ 980.616

		total_pw = total_pw + pw_layer
		k = k + 1
	enddo


	slope = ( wprof(k) -  wprof(k-1) ) / &
          ( pprof(k) - pprof(k-1))  

	watervapor_upper = wprof(k-1)
	watervapor_lower = slope * (cloud_top_pressure - pprof(k-1)) + &
	                         watervapor_upper
	pw_layer = (cloud_top_pressure - pprof(k-1)) * &
	           (watervapor_upper + watervapor_lower) * 0.5/980.616
	abovecloud_watervapor = total_pw + pw_layer

  endif
	 
 end subroutine get_above_cloud_properties


 end module ancillary_module
