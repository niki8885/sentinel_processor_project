module stats_mod
    use iso_c_binding
    implicit none
    real(c_double), parameter :: NODATA = -9999.0d0

contains

    pure integer function flat3(t, r, c, rows, cols)
        integer, intent(in) :: t, r, c, rows, cols
        flat3 = (t - 1) * rows * cols + (r - 1) * cols + c
    end function flat3


    subroutine time_window_stats(arr, dates_days, n_times, rows, cols, &
            window_days, mean_out, std_out, slope_out) &
            bind(C, name = "time_window_stats")

        integer(c_int), intent(in), value :: n_times, rows, cols, window_days
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(in) :: dates_days(n_times)
        real(c_double), intent(out) :: mean_out(rows * cols)
        real(c_double), intent(out) :: std_out(rows * cols)
        real(c_double), intent(out) :: slope_out(rows * cols)

        integer :: r, c, t, pix, n_valid
        real(c_double) :: half_w, d_lo, d_hi, d_centre, val
        real(c_double) :: sum_v, sum_v2, mean_v, var_v
        real(c_double) :: sum_x, sum_y, sum_xx, sum_xy, denom, xc

        half_w = real(window_days, c_double) * 0.5d0

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                d_lo = huge(1.0d0)
                d_hi = -huge(1.0d0)
                do t = 1, n_times
                    val = arr(flat3(t, r, c, rows, cols))
                    if (abs(val - NODATA) > 1.0d-4) then
                        if (dates_days(t) < d_lo) d_lo = dates_days(t)
                        if (dates_days(t) > d_hi) d_hi = dates_days(t)
                    end if
                end do

                if (d_lo > d_hi) then
                    mean_out(pix) = NODATA
                    std_out(pix) = NODATA
                    slope_out(pix) = NODATA
                    cycle
                end if

                d_centre = 0.5d0 * (d_lo + d_hi)

                n_valid = 0
                sum_v = 0.0d0;  sum_v2 = 0.0d0
                sum_x = 0.0d0;  sum_y = 0.0d0
                sum_xx = 0.0d0;  sum_xy = 0.0d0

                do t = 1, n_times
                    if (abs(dates_days(t) - d_centre) > half_w) cycle
                    val = arr(flat3(t, r, c, rows, cols))
                    if (abs(val - NODATA) <= 1.0d-4) cycle

                    xc = dates_days(t) - d_centre
                    n_valid = n_valid + 1
                    sum_v = sum_v + val
                    sum_v2 = sum_v2 + val * val
                    sum_x = sum_x + xc
                    sum_y = sum_y + val
                    sum_xx = sum_xx + xc * xc
                    sum_xy = sum_xy + xc * val
                end do

                if (n_valid == 0) then
                    mean_out(pix) = NODATA
                    std_out(pix) = NODATA
                    slope_out(pix) = NODATA
                    cycle
                end if

                ! mean
                mean_v = sum_v / real(n_valid, c_double)
                mean_out(pix) = mean_v

                if (n_valid > 1) then
                    var_v = (sum_v2 - real(n_valid, c_double) * mean_v * mean_v) &
                            / real(n_valid - 1, c_double)
                    if (var_v < 0.0d0) var_v = 0.0d0
                    std_out(pix) = sqrt(var_v)
                else
                    std_out(pix) = 0.0d0
                end if

                denom = real(n_valid, c_double) * sum_xx - sum_x * sum_x
                if (abs(denom) > 1.0d-12) then
                    slope_out(pix) = (real(n_valid, c_double) * sum_xy &
                            - sum_x * sum_y) / denom
                else
                    slope_out(pix) = 0.0d0
                end if

            end do
        end do
    end subroutine time_window_stats


    subroutine anomaly_zscore(arr, mean_in, std_in, n_times, rows, cols, &
            zscore_out) &
            bind(C, name = "anomaly_zscore")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(in) :: mean_in(rows * cols)
        real(c_double), intent(in) :: std_in(rows * cols)
        real(c_double), intent(out) :: zscore_out(n_times * rows * cols)

        integer :: r, c, t, pix, idx
        real(c_double) :: val, mu, sigma

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c
                mu = mean_in(pix)
                sigma = std_in(pix)

                do t = 1, n_times
                    idx = flat3(t, r, c, rows, cols)
                    val = arr(idx)

                    if (abs(val - NODATA) <= 1.0d-4 .or. &
                            abs(mu - NODATA) <= 1.0d-4 .or. &
                            abs(sigma - NODATA) <= 1.0d-4 .or. &
                            sigma < 1.0d-12) then
                        zscore_out(idx) = NODATA
                    else
                        zscore_out(idx) = (val - mu) / sigma
                    end if
                end do
            end do
        end do
    end subroutine anomaly_zscore

    subroutine trend_theil_sen(arr, dates_days, n_times, rows, cols, &
            slope_out, intercept_out) &
            bind(C, name = "trend_theil_sen")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(in) :: dates_days(n_times)
        real(c_double), intent(out) :: slope_out(rows * cols)
        real(c_double), intent(out) :: intercept_out(rows * cols)

        integer :: r, c, t, pix, nv, k, i, j, npairs
        real(c_double), allocatable :: xv(:), yv(:), pairs(:)
        real(c_double) :: slope, intercept, med_x, med_y

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do

                if (nv == 0) then
                    slope_out(pix) = NODATA
                    intercept_out(pix) = NODATA
                    cycle
                end if

                if (nv == 1) then
                    slope_out(pix) = 0.0d0
                    do t = 1, n_times
                        if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) then
                            intercept_out(pix) = arr(flat3(t, r, c, rows, cols))
                            exit
                        end if
                    end do
                    cycle
                end if

                allocate(xv(nv), yv(nv))
                k = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) then
                        k = k + 1
                        xv(k) = dates_days(t)
                        yv(k) = arr(flat3(t, r, c, rows, cols))
                    end if
                end do

                ! pairwise slopes
                npairs = nv * (nv - 1) / 2
                allocate(pairs(npairs))
                k = 0
                do i = 1, nv - 1
                    do j = i + 1, nv
                        if (abs(xv(j) - xv(i)) > 1.0d-12) then
                            k = k + 1
                            pairs(k) = (yv(j) - yv(i)) / (xv(j) - xv(i))
                        end if
                    end do
                end do

                if (k == 0) then
                    slope = 0.0d0
                else
                    call insertion_sort(pairs, k)
                    slope = median_sorted(pairs, k)
                end if

                med_x = median_of(xv, nv)
                med_y = median_of(yv, nv)
                intercept = med_y - slope * med_x

                slope_out(pix) = slope
                intercept_out(pix) = intercept

                deallocate(xv, yv, pairs)
            end do
        end do
    end subroutine trend_theil_sen


    subroutine insertion_sort(a, n)
        integer, intent(in) :: n
        real(c_double), intent(inout) :: a(n)
        integer :: i, j
        real(c_double) :: key
        do i = 2, n
            key = a(i); j = i - 1
            do while (j >= 1 .and. a(j) > key)
                a(j + 1) = a(j); j = j - 1
            end do
            a(j + 1) = key
        end do
    end subroutine insertion_sort

    pure real(c_double) function median_sorted(a, n)
        integer, intent(in) :: n
        real(c_double), intent(in) :: a(n)
        integer :: half
        half = n / 2
        if (mod(n, 2) == 1) then
            median_sorted = a(half + 1)
        else
            median_sorted = 0.5d0 * (a(half) + a(half + 1))
        end if
    end function median_sorted

    real(c_double) function median_of(a, n)
        integer, intent(in) :: n
        real(c_double), intent(in) :: a(n)
        real(c_double), allocatable :: tmp(:)
        allocate(tmp(n))
        tmp = a
        call insertion_sort(tmp, n)
        median_of = median_sorted(tmp, n)
        deallocate(tmp)
    end function median_of

end module stats_mod