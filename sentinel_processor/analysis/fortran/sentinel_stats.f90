module sentinel_stats
    use iso_c_binding
    implicit none
    real(c_double), parameter :: NODATA = -9999.0d0
    real(c_double), parameter :: PI = 3.14159265358979323846d0

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


    subroutine valid_obs_count(arr, n_times, rows, cols, count_out, fraction_out) &
            bind(C, name = "valid_obs_count")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(out) :: count_out(rows * cols)
        real(c_double), intent(out) :: fraction_out(rows * cols)

        integer :: r, c, t, pix, nv

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c
                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do
                count_out(pix) = real(nv, c_double)
                fraction_out(pix) = real(nv, c_double) / real(n_times, c_double)
            end do
        end do
    end subroutine valid_obs_count


    subroutine temporal_gap_stats(arr, dates_days, n_times, rows, cols, &
            max_gap_out, mean_gap_out) &
            bind(C, name = "temporal_gap_stats")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(in) :: dates_days(n_times)
        real(c_double), intent(out) :: max_gap_out(rows * cols)
        real(c_double), intent(out) :: mean_gap_out(rows * cols)

        integer :: r, c, t, pix, nv, k
        real(c_double), allocatable :: xv(:)
        real(c_double) :: gap, max_gap, sum_gap

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do

                if (nv < 2) then
                    max_gap_out(pix) = NODATA
                    mean_gap_out(pix) = NODATA
                    cycle
                end if

                allocate(xv(nv))
                k = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) then
                        k = k + 1
                        xv(k) = dates_days(t)
                    end if
                end do

                max_gap = 0.0d0
                sum_gap = 0.0d0
                do k = 2, nv
                    gap = xv(k) - xv(k - 1)
                    if (gap > max_gap) max_gap = gap
                    sum_gap = sum_gap + gap
                end do

                max_gap_out(pix) = max_gap
                mean_gap_out(pix) = sum_gap / real(nv - 1, c_double)
                deallocate(xv)
            end do
        end do
    end subroutine temporal_gap_stats


    subroutine pixel_quantiles(arr, n_times, rows, cols, &
            p10, p25, p50, p75, p90) &
            bind(C, name = "pixel_quantiles")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(out) :: p10(rows * cols)
        real(c_double), intent(out) :: p25(rows * cols)
        real(c_double), intent(out) :: p50(rows * cols)
        real(c_double), intent(out) :: p75(rows * cols)
        real(c_double), intent(out) :: p90(rows * cols)

        integer :: r, c, t, pix, nv, k
        real(c_double), allocatable :: yv(:)

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do

                if (nv == 0) then
                    p10(pix) = NODATA; p25(pix) = NODATA; p50(pix) = NODATA
                    p75(pix) = NODATA; p90(pix) = NODATA
                    cycle
                end if

                allocate(yv(nv))
                k = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) then
                        k = k + 1
                        yv(k) = arr(flat3(t, r, c, rows, cols))
                    end if
                end do
                call insertion_sort(yv, nv)

                if (nv == 1) then
                    p10(pix) = yv(1); p25(pix) = yv(1); p50(pix) = yv(1)
                    p75(pix) = yv(1); p90(pix) = yv(1)
                else
                    p10(pix) = interp_quantile(yv, nv, 0.10d0)
                    p25(pix) = interp_quantile(yv, nv, 0.25d0)
                    p50(pix) = interp_quantile(yv, nv, 0.50d0)
                    p75(pix) = interp_quantile(yv, nv, 0.75d0)
                    p90(pix) = interp_quantile(yv, nv, 0.90d0)
                end if
                deallocate(yv)
            end do
        end do
    end subroutine pixel_quantiles


    subroutine pixel_iqr(arr, n_times, rows, cols, iqr_out, outlier_out) &
            bind(C, name = "pixel_iqr")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(out) :: iqr_out(rows * cols)
        real(c_double), intent(out) :: outlier_out(n_times * rows * cols)

        integer :: r, c, t, pix, nv, k, idx
        real(c_double), allocatable :: yv(:)
        real(c_double) :: q25, q75, iqr, lo, hi, val

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do

                if (nv < 2) then
                    iqr_out(pix) = NODATA
                    do t = 1, n_times
                        outlier_out(flat3(t, r, c, rows, cols)) = NODATA
                    end do
                    cycle
                end if

                allocate(yv(nv))
                k = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) then
                        k = k + 1
                        yv(k) = arr(flat3(t, r, c, rows, cols))
                    end if
                end do
                call insertion_sort(yv, nv)

                q25 = interp_quantile(yv, nv, 0.25d0)
                q75 = interp_quantile(yv, nv, 0.75d0)
                iqr = q75 - q25
                iqr_out(pix) = iqr
                lo = q25 - 1.5d0 * iqr
                hi = q75 + 1.5d0 * iqr

                do t = 1, n_times
                    idx = flat3(t, r, c, rows, cols)
                    val = arr(idx)
                    if (abs(val - NODATA) <= 1.0d-4) then
                        outlier_out(idx) = NODATA
                    else if (val < lo .or. val > hi) then
                        outlier_out(idx) = 1.0d0
                    else
                        outlier_out(idx) = 0.0d0
                    end if
                end do
                deallocate(yv)
            end do
        end do
    end subroutine pixel_iqr


    subroutine mann_kendall(arr, n_times, rows, cols, &
            s_out, vars_out, z_out, trend_out) &
            bind(C, name = "mann_kendall")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(out) :: s_out(rows * cols)
        real(c_double), intent(out) :: vars_out(rows * cols)
        real(c_double), intent(out) :: z_out(rows * cols)
        real(c_double), intent(out) :: trend_out(rows * cols)

        integer :: r, c, t, pix, nv, i, j, k
        real(c_double), allocatable :: yv(:)
        real(c_double) :: S, varS, Z, nf, tie_corr, prev, run

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do

                if (nv < 4) then
                    s_out(pix) = NODATA; vars_out(pix) = NODATA
                    z_out(pix) = NODATA; trend_out(pix) = NODATA
                    cycle
                end if

                allocate(yv(nv))
                k = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) then
                        k = k + 1
                        yv(k) = arr(flat3(t, r, c, rows, cols))
                    end if
                end do

                ! S statistic
                S = 0.0d0
                do i = 1, nv - 1
                    do j = i + 1, nv
                        if     (yv(j) > yv(i)) then; S = S + 1.0d0
                        else if(yv(j) < yv(i)) then; S = S - 1.0d0
                        end if
                    end do
                end do

                nf = real(nv, c_double)
                varS = nf * (nf - 1.0d0) * (2.0d0 * nf + 5.0d0) / 18.0d0

                call insertion_sort(yv, nv)
                i = 1
                tie_corr = 0.0d0
                do while (i <= nv)
                    prev = yv(i); run = 1.0d0; j = i + 1
                    do while (j <= nv .and. abs(yv(j) - prev) < 1.0d-14)
                        run = run + 1.0d0; j = j + 1
                    end do
                    if (run > 1.0d0) &
                            tie_corr = tie_corr + run * (run - 1.0d0) * (2.0d0 * run + 5.0d0) / 18.0d0
                    i = j
                end do
                varS = varS - tie_corr
                if (varS < 0.0d0) varS = 0.0d0

                if (varS < 1.0d-12) then
                    Z = 0.0d0
                else if (S > 0.0d0) then
                    Z = (S - 1.0d0) / sqrt(varS)
                else if (S < 0.0d0) then
                    Z = (S + 1.0d0) / sqrt(varS)
                else
                    Z = 0.0d0
                end if

                s_out(pix) = S
                vars_out(pix) = varS
                z_out(pix) = Z

                if     (Z >  1.96d0) then; trend_out(pix) = 1.0d0
                else if(Z < -1.96d0) then; trend_out(pix) = -1.0d0
                else;                      trend_out(pix) = 0.0d0
                end if

                deallocate(yv)
            end do
        end do
    end subroutine mann_kendall


    subroutine bfast_breakpoint(arr, dates_days, n_times, rows, cols, &
            break_day_out, magnitude_out, rss_ratio_out) &
            bind(C, name = "bfast_breakpoint")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(in) :: dates_days(n_times)
        real(c_double), intent(out) :: break_day_out(rows * cols)
        real(c_double), intent(out) :: magnitude_out(rows * cols)
        real(c_double), intent(out) :: rss_ratio_out(rows * cols)

        integer :: r, c, t, pix, nv, k, bp, best_bp
        real(c_double), allocatable :: xv(:), yv(:)
        real(c_double) :: rss_one, rss_best, rss_cur
        real(c_double) :: y1, y2

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do

                if (nv < 6) then
                    break_day_out(pix) = NODATA
                    magnitude_out(pix) = NODATA
                    rss_ratio_out(pix) = NODATA
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

                rss_one = ols_rss(xv, yv, 1, nv)

                rss_best = huge(1.0d0)
                best_bp = 3
                do bp = 3, nv - 3
                    rss_cur = ols_rss(xv, yv, 1, bp) + ols_rss(xv, yv, bp + 1, nv)
                    if (rss_cur < rss_best) then
                        rss_best = rss_cur
                        best_bp = bp
                    end if
                end do

                y1 = ols_predict(xv, yv, 1, best_bp, xv(best_bp))
                y2 = ols_predict(xv, yv, best_bp + 1, nv, xv(best_bp + 1))

                break_day_out(pix) = xv(best_bp)
                magnitude_out(pix) = y2 - y1
                if (rss_one > 1.0d-12) then
                    rss_ratio_out(pix) = rss_best / rss_one
                else
                    rss_ratio_out(pix) = 1.0d0
                end if

                deallocate(xv, yv)
            end do
        end do
    end subroutine bfast_breakpoint


    subroutine phenology_doy(arr, dates_days, n_times, rows, cols, &
            rising_pct, falling_pct, &
            sos_out, pos_out, eos_out) &
            bind(C, name = "phenology_doy")

        integer(c_int), intent(in), value :: n_times, rows, cols
        integer(c_int), intent(in), value :: rising_pct, falling_pct
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(in) :: dates_days(n_times)
        real(c_double), intent(out) :: sos_out(rows * cols)
        real(c_double), intent(out) :: pos_out(rows * cols)
        real(c_double), intent(out) :: eos_out(rows * cols)

        integer :: r, c, t, pix, nv, k, peak_k
        real(c_double), allocatable :: xv(:), yv(:)
        real(c_double) :: vmin, vmax, amp, r_thr, f_thr, peak_v
        real(c_double) :: sos, eos
        logical :: found_sos, found_eos

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do

                if (nv < 4) then
                    sos_out(pix) = NODATA; pos_out(pix) = NODATA
                    eos_out(pix) = NODATA; cycle
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

                vmin = minval(yv(1:nv)); vmax = maxval(yv(1:nv))
                amp = vmax - vmin

                if (amp < 1.0d-6) then
                    sos_out(pix) = NODATA; pos_out(pix) = NODATA
                    eos_out(pix) = NODATA; deallocate(xv, yv); cycle
                end if

                r_thr = vmin + real(rising_pct, c_double) / 100.0d0 * amp
                f_thr = vmin + real(falling_pct, c_double) / 100.0d0 * amp

                peak_v = yv(1); peak_k = 1
                do k = 2, nv
                    if (yv(k) > peak_v) then; peak_v = yv(k); peak_k = k;
                    end if
                end do
                pos_out(pix) = xv(peak_k)

                found_sos = .false.
                sos = xv(1)
                do k = 2, peak_k
                    if (yv(k - 1) < r_thr .and. yv(k) >= r_thr) then
                        sos = xv(k - 1) + (r_thr - yv(k - 1)) / &
                                max(yv(k) - yv(k - 1), 1.0d-12) * (xv(k) - xv(k - 1))
                        found_sos = .true.
                        exit
                    end if
                end do
                sos_out(pix) = merge(sos, NODATA, found_sos)

                found_eos = .false.
                eos = xv(nv)
                do k = nv, peak_k + 1, -1
                    if (yv(k) < f_thr .and. yv(k - 1) >= f_thr) then
                        eos = xv(k - 1) + (f_thr - yv(k - 1)) / &
                                (yv(k) - yv(k - 1)) * (xv(k) - xv(k - 1))
                        found_eos = .true.
                        exit
                    end if
                end do
                eos_out(pix) = merge(eos, NODATA, found_eos)

                deallocate(xv, yv)
            end do
        end do
    end subroutine phenology_doy


    subroutine pearson_map(arr_a, arr_b, n_times, rows, cols, r_out) &
            bind(C, name = "pearson_map")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: arr_a(n_times * rows * cols)
        real(c_double), intent(in) :: arr_b(n_times * rows * cols)
        real(c_double), intent(out) :: r_out(rows * cols)

        integer :: r, c, t, pix, n, idx
        real(c_double) :: va, vb, sum_a, sum_b
        real(c_double) :: mean_a, mean_b, cov, var_a, var_b, denom

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c
                n = 0
                sum_a = 0.0d0; sum_b = 0.0d0

                do t = 1, n_times
                    idx = flat3(t, r, c, rows, cols)
                    va = arr_a(idx); vb = arr_b(idx)
                    if (abs(va - NODATA) > 1.0d-4 .and. abs(vb - NODATA) > 1.0d-4) then
                        n = n + 1
                        sum_a = sum_a + va
                        sum_b = sum_b + vb
                    end if
                end do

                if (n < 3) then; r_out(pix) = NODATA; cycle;
                end if

                mean_a = sum_a / real(n, c_double)
                mean_b = sum_b / real(n, c_double)

                cov = 0.0d0; var_a = 0.0d0; var_b = 0.0d0
                do t = 1, n_times
                    idx = flat3(t, r, c, rows, cols)
                    va = arr_a(idx); vb = arr_b(idx)
                    if (abs(va - NODATA) > 1.0d-4 .and. abs(vb - NODATA) > 1.0d-4) then
                        cov = cov + (va - mean_a) * (vb - mean_b)
                        var_a = var_a + (va - mean_a) * (va - mean_a)
                        var_b = var_b + (vb - mean_b) * (vb - mean_b)
                    end if
                end do

                denom = sqrt(var_a * var_b)
                if (denom > 1.0d-12) then
                    r_out(pix) = cov / denom
                else
                    r_out(pix) = 0.0d0
                end if
            end do
        end do
    end subroutine pearson_map


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

    pure real(c_double) function interp_quantile(a, n, p)
        integer, intent(in) :: n
        real(c_double), intent(in) :: a(n)
        real(c_double), intent(in) :: p
        real(c_double) :: h
        integer :: lo
        h = p * real(n - 1, c_double)
        lo = int(h) + 1
        if (lo >= n) then
            interp_quantile = a(n)
        else
            interp_quantile = a(lo) + (h - real(lo - 1, c_double)) * (a(lo + 1) - a(lo))
        end if
    end function interp_quantile


    real(c_double) function ols_rss(x, y, lo, hi)
        real(c_double), intent(in) :: x(:), y(:)
        integer, intent(in) :: lo, hi
        integer :: i, nn
        real(c_double) :: sx, sy, sxx, sxy, denom, b, a_int, res
        nn = hi - lo + 1
        if (nn < 2) then; ols_rss = 0.0d0; return;
        end if
        sx = 0.0d0; sy = 0.0d0; sxx = 0.0d0; sxy = 0.0d0
        do i = lo, hi
            sx = sx + x(i); sy = sy + y(i)
            sxx = sxx + x(i) * x(i); sxy = sxy + x(i) * y(i)
        end do
        denom = real(nn, c_double) * sxx - sx * sx
        if (abs(denom) < 1.0d-12) then
            b = 0.0d0; a_int = sy / real(nn, c_double)
        else
            b = (real(nn, c_double) * sxy - sx * sy) / denom
            a_int = (sy - b * sx) / real(nn, c_double)
        end if
        ols_rss = 0.0d0
        do i = lo, hi
            res = y(i) - (a_int + b * x(i))
            ols_rss = ols_rss + res * res
        end do
    end function ols_rss

    real(c_double) function ols_predict(x, y, lo, hi, xq)
        real(c_double), intent(in) :: x(:), y(:), xq
        integer, intent(in) :: lo, hi
        integer :: i, nn
        real(c_double) :: sx, sy, sxx, sxy, denom, b, a_int
        nn = hi - lo + 1
        if (nn < 2) then; ols_predict = y(lo); return;
        end if
        sx = 0.0d0; sy = 0.0d0; sxx = 0.0d0; sxy = 0.0d0
        do i = lo, hi
            sx = sx + x(i); sy = sy + y(i)
            sxx = sxx + x(i) * x(i); sxy = sxy + x(i) * y(i)
        end do
        denom = real(nn, c_double) * sxx - sx * sx
        if (abs(denom) < 1.0d-12) then
            b = 0.0d0; a_int = sy / real(nn, c_double)
        else
            b = (real(nn, c_double) * sxy - sx * sy) / denom
            a_int = (sy - b * sx) / real(nn, c_double)
        end if
        ols_predict = a_int + b * xq
    end function ols_predict


    subroutine pixel_regression(y_arr, x_arr, n_times, rows, cols, &
            slope_out, intercept_out, r2_out) &
            bind(C, name = "pixel_regression")

        integer(c_int), intent(in), value :: n_times, rows, cols
        real(c_double), intent(in) :: y_arr(n_times * rows * cols)
        real(c_double), intent(in) :: x_arr(n_times)
        real(c_double), intent(out) :: slope_out(rows * cols)
        real(c_double), intent(out) :: intercept_out(rows * cols)
        real(c_double), intent(out) :: r2_out(rows * cols)

        integer :: r, c, t, pix, nv
        real(c_double) :: xi, yi, sx, sy, sxx, sxy, syy
        real(c_double) :: nf, denom, slope, intercept, ss_res, ss_tot, r2, mean_y

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                sx = 0.0d0;  sy = 0.0d0
                sxx = 0.0d0;  sxy = 0.0d0; syy = 0.0d0

                do t = 1, n_times
                    yi = y_arr(flat3(t, r, c, rows, cols))
                    if (abs(yi - NODATA) <= 1.0d-4) cycle
                    xi = x_arr(t)
                    nv = nv + 1
                    sx = sx + xi
                    sy = sy + yi
                    sxx = sxx + xi * xi
                    sxy = sxy + xi * yi
                    syy = syy + yi * yi
                end do

                if (nv < 2) then
                    slope_out(pix) = NODATA
                    intercept_out(pix) = NODATA
                    r2_out(pix) = NODATA
                    cycle
                end if

                nf = real(nv, c_double)
                denom = nf * sxx - sx * sx

                if (abs(denom) < 1.0d-12 * nf * nf) then
                    slope_out(pix) = 0.0d0
                    intercept_out(pix) = sy / nf
                    r2_out(pix) = 0.0d0
                    cycle
                end if

                slope = (nf * sxy - sx * sy) / denom
                intercept = (sy - slope * sx) / nf

                mean_y = sy / nf
                ss_tot = syy - nf * mean_y * mean_y

                ss_res = ss_tot - slope * (sxy - nf * (sx / nf) * mean_y)

                if (ss_tot < 1.0d-12) then
                    r2 = 0.0d0
                else
                    r2 = 1.0d0 - ss_res / ss_tot
                    if (r2 < 0.0d0) r2 = 0.0d0
                    if (r2 > 1.0d0) r2 = 1.0d0
                end if

                slope_out(pix) = slope
                intercept_out(pix) = intercept
                r2_out(pix) = r2

            end do
        end do
    end subroutine pixel_regression


    subroutine phenology_doy_v2(arr, dates_days, n_times, rows, cols, &
            rising_pct, falling_pct, &
            sos_out, eos_out, peak_doy_out, peak_val_out) &
            bind(C, name = "phenology_doy_v2")

        integer(c_int), intent(in), value :: n_times, rows, cols
        integer(c_int), intent(in), value :: rising_pct, falling_pct
        real(c_double), intent(in) :: arr(n_times * rows * cols)
        real(c_double), intent(in) :: dates_days(n_times)
        real(c_double), intent(out) :: sos_out(rows * cols)
        real(c_double), intent(out) :: eos_out(rows * cols)
        real(c_double), intent(out) :: peak_doy_out(rows * cols)
        real(c_double), intent(out) :: peak_val_out(rows * cols)

        integer :: r, c, t, pix, nv, k, peak_k
        real(c_double), allocatable :: xv(:), yv(:)
        real(c_double) :: vmin, vmax, amp, r_thr, f_thr, peak_v
        real(c_double) :: sos, eos
        logical :: found_sos, found_eos

        do r = 1, rows
            do c = 1, cols
                pix = (r - 1) * cols + c

                nv = 0
                do t = 1, n_times
                    if (abs(arr(flat3(t, r, c, rows, cols)) - NODATA) > 1.0d-4) &
                            nv = nv + 1
                end do

                if (nv < 5) then
                    sos_out(pix) = NODATA
                    eos_out(pix) = NODATA
                    peak_doy_out(pix) = NODATA
                    peak_val_out(pix) = NODATA
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

                vmin = minval(yv(1:nv));  vmax = maxval(yv(1:nv))
                amp = vmax - vmin

                if (amp < 1.0d-6) then
                    sos_out(pix) = NODATA
                    eos_out(pix) = NODATA
                    peak_doy_out(pix) = NODATA
                    peak_val_out(pix) = NODATA
                    deallocate(xv, yv);  cycle
                end if

                r_thr = vmin + real(rising_pct, c_double) / 100.0d0 * amp
                f_thr = vmin + real(falling_pct, c_double) / 100.0d0 * amp

                peak_v = yv(1);  peak_k = 1
                do k = 2, nv
                    if (yv(k) > peak_v) then
                        peak_v = yv(k);  peak_k = k
                    end if
                end do
                peak_doy_out(pix) = xv(peak_k)
                peak_val_out(pix) = peak_v

                found_sos = .false.
                sos = xv(1)
                do k = 2, peak_k
                    if (yv(k - 1) < r_thr .and. yv(k) >= r_thr) then
                        sos = xv(k - 1) + (r_thr - yv(k - 1)) / &
                                max(yv(k) - yv(k - 1), 1.0d-12) * (xv(k) - xv(k - 1))
                        found_sos = .true.
                        exit
                    end if
                end do
                sos_out(pix) = merge(sos, NODATA, found_sos)

                found_eos = .false.
                eos = xv(nv)
                do k = nv, peak_k + 1, -1
                    if (yv(k) < f_thr .and. yv(k - 1) >= f_thr) then
                        eos = xv(k - 1) + (f_thr - yv(k - 1)) / &
                                (yv(k) - yv(k - 1)) * (xv(k) - xv(k - 1))
                        found_eos = .true.
                        exit
                    end if
                end do
                eos_out(pix) = merge(eos, NODATA, found_eos)

                deallocate(xv, yv)
            end do
        end do
    end subroutine phenology_doy_v2


    subroutine savgol_smooth_stack(arr, n_times, rows, cols, window) &
            bind(C, name = "savgol_smooth_stack")

        integer(c_int), intent(in), value :: n_times, rows, cols, window
        real(c_double), intent(inout) :: arr(n_times * rows * cols)

        integer :: r, c, t, hw, half, j, jj
        real(c_double) :: s0, s1, s2, s3, s4, sy, sxy, sx2y, denom, a0, xj
        real(c_double), allocatable :: tmp(:), sg(:)

        hw = max(window / 2, 1)
        if (2 * hw + 1 < 3) hw = 1

        allocate(tmp(n_times), sg(n_times))

        do r = 1, rows
            do c = 1, cols
                do t = 1, n_times
                    tmp(t) = arr(flat3(t, r, c, rows, cols))
                end do

                if (all(abs(tmp - NODATA) <= 1.0d-4)) then
                    cycle
                end if

                do t = 1, n_times
                    half = min(hw, t - 1, n_times - t)
                    s0 = 0.0d0; s1 = 0.0d0; s2 = 0.0d0
                    s3 = 0.0d0; s4 = 0.0d0
                    sy = 0.0d0; sxy = 0.0d0; sx2y = 0.0d0
                    do j = -half, half
                        jj = t + j;  xj = real(j, c_double)
                        s0 = s0 + 1.0d0
                        s1 = s1 + xj
                        s2 = s2 + xj * xj
                        s3 = s3 + xj**3
                        s4 = s4 + xj**4
                        sy = sy + tmp(jj)
                        sxy = sxy + xj * tmp(jj)
                        sx2y = sx2y + xj * xj * tmp(jj)
                    end do
                    denom = s0 * (s2 * s4 - s3 * s3) &
                            - s1 * (s1 * s4 - s3 * s2) &
                            + s2 * (s1 * s3 - s2 * s2)
                    if (abs(denom) < 1.0d-12) then
                        sg(t) = tmp(t)
                    else
                        a0 = (sy * (s2 * s4 - s3 * s3) &
                                - sxy * (s1 * s4 - s3 * s2) &
                                + sx2y * (s1 * s3 - s2 * s2)) / denom
                        sg(t) = a0
                    end if
                end do

                do t = 1, n_times
                    if (abs(tmp(t) - NODATA) > 1.0d-4) then
                        arr(flat3(t, r, c, rows, cols)) = sg(t)
                    end if
                end do

            end do
        end do

        deallocate(tmp, sg)
    end subroutine savgol_smooth_stack

end module sentinel_stats