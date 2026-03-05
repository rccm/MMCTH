subroutine allocate_arrays ( edge, status ) 
!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
!    allocate_arrays allocates data arrays (defined in core_arrays)
!    to hold MODIS L1B data.
!
! !INPUT PARAMETERS:
!
!    Type       Name             Description
!    ====       ====             ===========
!    INTEGER*4  edge             array, holds dimensions of L1B data
!                                 which are used as array sizes for
!                                 the allocation
!
! !OUTPUT PARAMETERS:
!
!    Type       Name             Description
!    ====       ====             ===========
!    INTEGER*4  status           return status
!
! !REVISION HISTORY:
!
!    Initial Version by Jennifer Wei, July 8, 2004
!    Jeremy Warner added prolog and several arrays, December 1, 2006
!                  and message handling, also removed
!                  unused calling arguments
!
! !TEAM-UNIQUE HEADER:
!
!    This software is developed by the Deep Blue Science Team
!    for the National Aeronautics and Space Administration,
!    Goddard Space Flight Center, under contract NAS5-02041.
!
! !REFERENCES AND CREDITS
!
! !DESIGN NOTES:
!
!
!   Externals:
!
!     MODIS_W_GENERIC            (MODIS_39500.f)
!
!   Functions:
!
!     MODIS_SMF_SETDYNAMICMSG
!
! !END
!-----------------------------------------------------------------------

   use GeneralAuxType
   use core_arrays
   use hdf
   
   implicit none
   
   include 'PGS_MODIS_39500.f'
   INCLUDE 'PGS_SMF.f'

   integer, dimension(:), intent(in)     :: edge
   integer,               intent (out)   :: status

   integer                               :: checkvariable, filelun, number, number_flux_solzen
   integer                               :: xdimension, ydimension, number_of_bands, nscan
   integer                               :: i,j,k
    
   logical         :: allocation_status, array_size_change
   character(256)  :: msg
   
   status = 0
   number_of_bands = size(bands)
   xdimension =  edge(1)
   ydimension =  edge(2)
	 nscan = ydimension/10
	 
!  Use latitude array as a size/allocation test to save time.
   allocation_status = allocated(latitude) 

   array_size_change = .false.
   
   if (allocation_status) then
      if (size(latitude,1) /= xdimension .or. &
          size(latitude,2) /= ydimension) then
         array_size_change  = .true.
      endif
   endif

  if(array_size_change .and. allocated(latitude)) then
     deallocate(latitude, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (latitude) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(latitude)) then
   allocate (latitude(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (latitude) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif

! Longitude 
  if(array_size_change .and. allocated(longitude)) then
     deallocate(longitude, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (longitude) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(longitude)) then
   allocate (longitude(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (longitude) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif

! landsea 
  if(array_size_change .and. allocated(landsea)) then
     deallocate(landsea, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (landsea) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(landsea)) then
   allocate (landsea(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (landsea) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif
	
! Sensor Zenith Angle   
  if(array_size_change .and. allocated(sensor_zenith_angle)) then
     deallocate(sensor_zenith_angle, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (sensor_zenith_angle) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(sensor_zenith_angle)) then
   allocate (sensor_zenith_angle(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (sensor_zenith_angle) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif
   
! Solar Zenith Angle
  if(array_size_change .and. allocated(solar_zenith_angle)) then
     deallocate(solar_zenith_angle, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (solar_zenith_angle) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(solar_zenith_angle)) then
   allocate (solar_zenith_angle(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (solar_zenith_angle) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif

! Solar Azimuth Angle
  if(array_size_change .and. allocated(solar_azimuth_angle)) then
     deallocate(solar_azimuth_angle, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (solar_azimuth_angle) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(solar_azimuth_angle)) then
   allocate (solar_azimuth_angle(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (solar_azimuth_angle) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif
  
! Sensor Azimuth Angle
    if(array_size_change .and. allocated(sensor_azimuth_angle)) then
     deallocate(sensor_azimuth_angle, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (sensor_azimuth_angle) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
! Sensor Azimuth Angle
  if (.not. allocated(sensor_azimuth_angle)) then 
   allocate (sensor_azimuth_angle(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (sensor_azimuth_angle) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif

   if(array_size_change .and. allocated(relative_azimuth_angle)) then
     deallocate(relative_azimuth_angle, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (relative_azimuth_angle) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
! Relative Azimuth Angle
  if (.not. allocated(relative_azimuth_angle)) then
   allocate (relative_azimuth_angle(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (relative_azimuth_angle) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif

! Added for polarization correction - 05 dec 2008 by CES
! note the dimension change for tt_inst2ecr
! T_inst2ECR
   array_size_change = .false.
   if(array_size_change .and. allocated(tt_inst2ecr)) then
     deallocate(tt_inst2ecr, stat = checkvariable)
     array_size_change = .true.
     if ( checkvariable /= 0 ) then
       msg = "Required array (tt_inst2ecr) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
   if (.not. allocated(tt_inst2ecr)) then
     allocate (tt_inst2ecr(1,3,nscan), stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (tt_inst2ecr) cannot be allocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
! Mirror_side
   array_size_change = .false.
   if(array_size_change .and. allocated(mirror_side)) then
     deallocate(mirror_side, stat = checkvariable)
     array_size_change = .true.
     if ( checkvariable /= 0 ) then
       msg = "Required array (mirror_side) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
!   write(*,*) 'mirror_side allocate? ', allocated(mirror_side)   
   if (.not. allocated(mirror_side)) then
     allocate (mirror_side(nscan), stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (mirror_side) cannot be allocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
! T_ECR
   array_size_change = .false.
   if(array_size_change .and. allocated(tt_ecr)) then
     deallocate(tt_ecr, stat = checkvariable)
     array_size_change = .true.
     if ( checkvariable /= 0 ) then
       msg = "Required array (tt_ecr) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
!   write(*,*) 'tt_ecr allocate? ', allocated(tt_ecr)   
   if (.not. allocated(tt_ecr)) then
     allocate (tt_ecr(3,nscan), stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (tt_ecr) cannot be allocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif

! EV center time
   array_size_change = .false.
   if(array_size_change .and. allocated(ev_cntr_time)) then
      deallocate(ev_cntr_time, stat = checkvariable)
      array_size_change = .true.
      if ( checkvariable /= 0 ) then
        msg = "Required array (ev_cntr_time) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
      endif
   endif
!   write(*,*) 'ev_cntr_time allocate?  ', allocated(ev_cntr_time)
   if (.not. allocated(ev_cntr_time)) then
     allocate (ev_cntr_time(nscan), stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (ev_cntr_time) cannot be allocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif


! band measurement from 36 channels....
! note the dimension order here
   if(array_size_change .and. allocated(band_measurements)) then
     deallocate(band_measurements, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (band_measurements) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(band_measurements)) then
   allocate (band_measurements(xdimension,number_of_bands, ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (band_measurements) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif

! Cloud_flag
  if(array_size_change .and. allocated(cloud_flag)) then
     deallocate(cloud_flag, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (cloud_flag) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(cloud_flag)) then
   allocate (cloud_flag(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (cloud_flag) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif

! Cloud_flag
  if(array_size_change .and. allocated(cirrus_flag)) then
     deallocate(cirrus_flag, stat = checkvariable)
     if ( checkvariable /= 0 ) then
       msg = "Required array (cirrus_flag) cannot be deallocated"
       call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
       status=1
     endif
   endif
   
  if (.not. allocated(cirrus_flag)) then
   allocate (cirrus_flag(xdimension,ydimension), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      msg = "Required array (cirrus_flag) cannot be allocated"
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'allocate_arrays')
      status=1
   endif
   endif

  
end subroutine allocate_arrays

