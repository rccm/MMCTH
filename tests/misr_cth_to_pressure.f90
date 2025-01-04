! File: misr_cth_to_pressure.f90
subroutine height_to_log_pressure(geo_height, misr_cth, rows, cols, pressure_at_cth)
    implicit none
    integer, intent(in) :: rows, cols
    real(4), intent(in) :: geo_height(101, rows, cols)
    real(4), intent(in) :: misr_cth(rows, cols)
    real(4), intent(out) :: pressure_at_cth(rows, cols)

    integer :: i, j, k
    real(8) :: interp_log_p, h_below, h_above, p_below, p_above
    integer, parameter :: num_levels = 101
    real(8), parameter :: int_press_o(101) = (/4.9999999e-03, 1.6100001e-02, 3.8400002e-02, 7.6899998e-02, 1.3699999e-01, &
                  2.2440000e-01, 3.4540001e-01, 5.0639999e-01, 7.1399999e-01, 9.7530001e-01, &
                  1.2972000e+00, 1.6872000e+00, 2.1526000e+00, 2.7009001e+00, 3.3397999e+00, &
                  4.0770001e+00, 4.9204001e+00, 5.8776002e+00, 6.9566998e+00, 8.1654997e+00, &
                  9.5118999e+00, 1.1003800e+01, 1.2649200e+01, 1.4455900e+01, 1.6431801e+01, &
                  1.8584700e+01, 2.0922400e+01, 2.3452600e+01, 2.6182899e+01, 2.9121000e+01, &
                  3.2274399e+01, 3.5650501e+01, 3.9256599e+01, 4.3100101e+01, 4.7188202e+01, &
                  5.1527802e+01, 5.6125999e+01, 6.0989498e+01, 6.6125298e+01, 7.1539803e+01, &
                  7.7239601e+01, 8.3231003e+01, 8.9520401e+01, 9.6113800e+01, 1.0301720e+02, &
                  1.1023660e+02, 1.1777750e+02, 1.2564560e+02, 1.3384621e+02, 1.4238480e+02, &
                  1.5126640e+02, 1.6049590e+02, 1.7007840e+02, 1.8001830e+02, 1.9032030e+02, &
                  2.0098869e+02, 2.1202769e+02, 2.2344150e+02, 2.3523380e+02, 2.4740849e+02, &
                  2.5996909e+02, 2.7291910e+02, 2.8626169e+02, 3.0000000e+02, 3.1413690e+02, &
                  3.2867529e+02, 3.4361761e+02, 3.5896649e+02, 3.7472409e+02, 3.9089261e+02, &
                  4.0747379e+02, 4.2446979e+02, 4.4188190e+02, 4.5971179e+02, 4.7796069e+02, &
                  4.9662979e+02, 5.1571997e+02, 5.3523218e+02, 5.5516687e+02, 5.7552478e+02, &
                  5.9630621e+02, 6.1751123e+02, 6.3913977e+02, 6.6119202e+02, 6.8366730e+02, &
                  7.0656543e+02, 7.2988568e+02, 7.5362750e+02, 7.7778967e+02, 8.0237140e+02, &
                  8.2737128e+02, 8.5278802e+02, 8.7862012e+02, 9.0486591e+02, 9.3152362e+02, &
                  9.5859113e+02, 9.8606659e+02, 1.0139476e+03, 1.0422319e+03, 1.0709170e+03, &
                  1.1000000e+03 /)

    real(8) :: log_press(101)
    
    ! Compute log_press inside the subroutine
    do k = 1, num_levels
        log_press(k) = log(int_press_o(k))
    end do

    ! Initialize pressure_at_cth with NaN or -9999 as placeholder
    pressure_at_cth = -9999.0

    do i = 1, rows
        do j = 1, cols
            ! Skip invalid misr_cth values
            if (misr_cth(i, j) < 0) cycle
            
            ! Find the interpolation levels for each pixel
            do k = 1, num_levels-1
                if (misr_cth(i, j) <= geo_height(k, i, j) .and. misr_cth(i, j) >= geo_height(k+1, i, j)) then                
                    h_below = geo_height(k, i, j)
                    h_above = geo_height(k+1, i, j)
                    p_below = log_press(k)
                    p_above = log_press(k+1)
                    
                    ! Perform linear interpolation
                    interp_log_p = p_below + (misr_cth(i, j) - h_below) * (p_above - p_below) / (h_above - h_below)
                    pressure_at_cth(i, j) = exp(interp_log_p)
                end if
            end do
            
        end do
    end do
end subroutine height_to_log_pressure
