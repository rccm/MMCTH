!=============================
! Linear interpolation function
!=============================
real(4) function linear_interp(target, x, y, n)
implicit none

integer, intent(in) :: n
real(4), intent(in) :: target
real(4), intent(in) :: x(n), y(n)

integer :: i
real(4), parameter :: missing_value = -999.0

! This function assumes x is monotonically increasing.
! If target is outside the x range, return missing.

if (target < x(1) .or. target > x(n)) then
  linear_interp = missing_value
  return
end if

if (n == 1) then
  linear_interp = y(1)
  return
end if

if (target == x(n)) then
  linear_interp = y(n)
  return
end if

do i = 1, n - 1
  if (x(i) <= target .and. target <= x(i + 1)) then
    if (x(i + 1) == x(i)) then
      linear_interp = y(i)
      return
    end if
    linear_interp = y(i) + (y(i + 1) - y(i)) * &
                    (target - x(i)) / (x(i + 1) - x(i))
    return
  end if
end do

! Safety fallback. Should not be reached if x is monotonic.
linear_interp = missing_value

end function linear_interp


!===============================================
! Main subroutine to interpolate to MODIS pressure levels
!===============================================
subroutine interpolate_to_pressure_levels(era5_levels, modis_levels, T_era5, W_era5, Z_era5, &
  surface_pressures, W_surface, T_skin, T_sst, T_modis, W_modis, Z_modis)

implicit none

! Input arrays
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

! Local variables
integer :: num_era5_levels
integer :: num_modis_levels
integer :: lat_size
integer :: lon_size
integer :: i, j, k
integer :: ilev
integer :: n_profile

real(4) :: surface_T
real(4) :: surface_p
real(4), parameter :: missing_value = -999.0

real(4), allocatable :: log_interp_levels(:)
real(4), allocatable :: T_values(:)
real(4), allocatable :: W_values(:)
real(4), allocatable :: Z_values(:)

! Interface for linear interpolation function
interface
  real(4) function linear_interp(target, x, y, n)
    implicit none
    integer, intent(in) :: n
    real(4), intent(in) :: target
    real(4), intent(in) :: x(n), y(n)
  end function linear_interp
end interface

! Get dimensions
num_era5_levels = size(era5_levels)
num_modis_levels = size(modis_levels)
lat_size = size(surface_pressures, 1)
lon_size = size(surface_pressures, 2)

! Temporary arrays include ERA5 levels plus inserted surface level
allocate(log_interp_levels(num_era5_levels + 1))
allocate(T_values(num_era5_levels + 1))
allocate(W_values(num_era5_levels + 1))
allocate(Z_values(num_era5_levels + 1))

! Loop over each spatial location
do j = 1, lat_size
  do k = 1, lon_size

    surface_p = surface_pressures(j, k)

    if (surface_p /= surface_p .or. surface_p <= 0.0) then
      T_modis(:, j, k) = missing_value
      W_modis(:, j, k) = missing_value
      Z_modis(:, j, k) = missing_value
      cycle
    end if

    ! Select the appropriate surface temperature.
    ! Use SST when valid; otherwise use skin temperature.
    surface_T = T_sst(j, k)
    if (surface_T < 0.0 .or. surface_T /= surface_T) then
      surface_T = T_skin(j, k)
    end if

    ! -------------------------------------------------------
    ! Insert surface pressure into the ERA5 pressure grid.
    !
    ! Only levels above the local surface are retained. ERA5 pressure-level
    ! fields below terrain are extrapolated and should not affect retrievals.
    ! -------------------------------------------------------

    n_profile = 0
    do ilev = 1, num_era5_levels
      if (era5_levels(ilev) < surface_p) then
        n_profile = n_profile + 1
        log_interp_levels(n_profile) = log(era5_levels(ilev))
        T_values(n_profile) = T_era5(ilev, j, k)
        W_values(n_profile) = W_era5(ilev, j, k)
        Z_values(n_profile) = Z_era5(ilev, j, k)
      end if
    end do

    n_profile = n_profile + 1
    log_interp_levels(n_profile) = log(surface_p)
    T_values(n_profile) = surface_T
    W_values(n_profile) = W_surface(j, k)
    Z_values(n_profile) = 0.0

    ! -------------------------------------------------------
    ! Interpolate to every MODIS pressure level.
    !
    ! Levels below the local surface are set missing. Keeping them populated
    ! lets below-ground ERA5 fields leak into MISR CTH->CTP conversion and the
    ! CO2-slicing retrieval.
    ! -------------------------------------------------------

    do i = 1, num_modis_levels
      if (modis_levels(i) > surface_p .or. modis_levels(i) < era5_levels(1)) then
        T_modis(i, j, k) = missing_value
        W_modis(i, j, k) = missing_value
        Z_modis(i, j, k) = missing_value
      else
        T_modis(i, j, k) = linear_interp(log(modis_levels(i)), &
                                         log_interp_levels, T_values, &
                                         n_profile)

        W_modis(i, j, k) = linear_interp(log(modis_levels(i)), &
                                         log_interp_levels, W_values, &
                                         n_profile)

        Z_modis(i, j, k) = linear_interp(log(modis_levels(i)), &
                                         log_interp_levels, Z_values, &
                                         n_profile)
      end if
    end do

  end do
end do

! Deallocate temporary arrays
deallocate(log_interp_levels)
deallocate(T_values)
deallocate(W_values)
deallocate(Z_values)

end subroutine interpolate_to_pressure_levels
