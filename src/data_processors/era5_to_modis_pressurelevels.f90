!=============================
! Linear interpolation function
!=============================
real(4) function linear_interp(target, x, y, n)
implicit none
integer, intent(in) :: n
real(4), intent(in) :: target, x(n), y(n)
integer :: i

do i = 1, n - 1
  if (x(i) <= target .and. target <= x(i + 1)) then
    linear_interp = y(i) + (y(i + 1) - y(i)) * &
                    (target - x(i)) / (x(i + 1) - x(i))
    return
  end if
end do

! If target is outside the range, return missing value
linear_interp = -999.0
end function linear_interp

!===============================================
! Main subroutine to interpolate to MODIS pressure levels
!===============================================
subroutine interpolate_to_pressure_levels(era5_levels, modis_levels, T_era5, W_era5, Z_era5, &
  surface_pressures, W_surface, T_skin, T_sst, T_modis, W_modis, Z_modis)
implicit none
! Declare arrays with assumed shape
real(4), intent(in) :: era5_levels(:)
real(4), intent(in) :: modis_levels(:)
real(4), intent(in) :: T_era5(:,:,:)
real(4), intent(in) :: W_era5(:,:,:)
real(4), intent(in) :: Z_era5(:,:,:)
real(4), intent(in) :: surface_pressures(:,:)
real(4), intent(in) :: W_surface(:,:)
real(4), intent(in) :: T_skin(:,:)
real(4), intent(in) :: T_sst(:,:)
 ! Output arrays
real(4), intent(inout) :: T_modis(:,:,:)
real(4), intent(inout) :: W_modis(:,:,:)
real(4), intent(inout) :: Z_modis(:,:,:)

integer :: num_era5_levels, num_modis_levels, lat_size, lon_size
integer :: i, j, k,insert_pos,ilev
real :: lowest_pressure  
real(4), allocatable :: log_interp_levels(:)
real(4), allocatable :: log_modis_levels(:)
real(4), allocatable :: T_values(:)
real(4), allocatable :: W_values(:)
real(4), allocatable :: Z_values(:)


real(4) :: surface_T

! Interface for linear interpolation function
interface
  real(4) function linear_interp(target, x, y, n)
    implicit none
    integer, intent(in) :: n
    real(4), intent(in) :: target, x(n), y(n)
  end function linear_interp
end interface

! Get dimensions from array sizes
num_era5_levels = size(era5_levels)
num_modis_levels = size(modis_levels)
lat_size = size(surface_pressures, 1)
lon_size = size(surface_pressures, 2)

! Allocate temporary arrays
allocate(log_interp_levels(num_era5_levels + 1))
allocate(log_modis_levels(num_modis_levels))
allocate(T_values(num_era5_levels + 1))
allocate(W_values(num_era5_levels + 1))
allocate(Z_values(num_era5_levels + 1))

! Compute the logarithm of MODIS pressure levels
do i = 1, num_modis_levels
  log_modis_levels(i) = log(modis_levels(i))
end do

! Loop over each spatial location (latitude and longitude)
do j = 1, lat_size
  do k = 1, lon_size
    ! Select the appropriate surface temperature
    surface_T = T_sst(j, k)
    if (surface_T < 0.0 .or. surface_T /= surface_T) then
      surface_T = T_skin(j, k)
    end if
    ! Combine ERA5 levels and surface pressure into one array
    insert_pos = num_era5_levels + 1
    do ilev = 1, num_era5_levels
      if (surface_pressures(j, k) >= era5_levels(i)) then
        insert_pos = ilev ! Position where surface_pressure should be inserted
        exit
      end if
    end do
    
    ! Build log_interp_levels with surface_pressure inserted
    log_interp_levels(1:insert_pos-1) = log(era5_levels(1:insert_pos-1))
    log_interp_levels(insert_pos) = log(surface_pressures(j, k))
    if (insert_pos <= num_era5_levels) then  ! Only assign upper part if not at end
      log_interp_levels(insert_pos+1:num_era5_levels+1) = log(era5_levels(insert_pos:num_era5_levels))
    end if
    

    T_values(1:insert_pos-1) = T_era5(1:insert_pos-1, j, k)
    T_values(insert_pos) = surface_T
    if (insert_pos <= num_era5_levels) then  ! Only assign upper part if not at end
      T_values(insert_pos+1:num_era5_levels+1) = T_era5(insert_pos:num_era5_levels, j, k)
    end if

    W_values(1:insert_pos-1) = W_era5(1:insert_pos-1, j, k)
    W_values(insert_pos) = W_surface(j, k)
    if (insert_pos <= num_era5_levels) then
      W_values(insert_pos+1:num_era5_levels+1) = W_era5(insert_pos:num_era5_levels, j, k)
    end if

    Z_values(1:insert_pos-1) = Z_era5(1:insert_pos-1, j, k)
    Z_values(insert_pos) = 0.0  ! Surface height
    if (insert_pos <= num_era5_levels) then
      Z_values(insert_pos+1:num_era5_levels+1) = Z_era5(insert_pos:num_era5_levels, j, k)
    end if
    ! Interpolate for each MODIS level
    !lowest_pressure = max(surface_pressures(j,k), era5_levels(num_era5_levels))
    do i = 1, num_modis_levels
      if (modis_levels(i) >surface_pressures(j,k) .or. modis_levels(i) < era5_levels(1) ) then
        ! Pressure level is below the surface; set to missing
        T_modis(i, j, k) = -999.0
        W_modis(i, j, k) = -999.0
        Z_modis(i, j, k) = -999.0
      else
        ! Interpolate between levels
        T_modis(i, j, k) = linear_interp(log(modis_levels(i)), log_interp_levels, T_values, num_era5_levels + 1)
        W_modis(i, j, k) = linear_interp(log(modis_levels(i)), log_interp_levels, W_values, num_era5_levels + 1)
        Z_modis(i, j, k) = linear_interp(log(modis_levels(i)), log_interp_levels, Z_values, num_era5_levels + 1)
      end if
    end do
  end do
end do

! Deallocate temporary arrays
deallocate(log_interp_levels)
deallocate(log_modis_levels)
deallocate(T_values)
deallocate(W_values)
deallocate(Z_values)
end subroutine interpolate_to_pressure_levels
