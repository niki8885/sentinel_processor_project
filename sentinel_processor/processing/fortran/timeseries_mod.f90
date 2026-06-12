! method  algorithm
! ------  -------------------------------------------------------------------
!   0     Linear interpolation   (piecewise, nearest-neighbour bounds)
!   1     Savitzky-Golay         (quadratic, window odd >= 3)
!   2     PCHIP cubic Hermite    (Fritsch-Carlson, no overshoot)
!   3     Holt double exponential smoothing (auto alpha/beta via grid search)
!   4     Gaussian-weighted moving average  (sigma = window/4, renormalised)

module timeseries_mod
    use iso_c_binding
    implicit none
    real(c_double), parameter :: NODATA = -9999.0d0

contains

    subroutine interpolate_gaps(arr, mask, n_times, rows, cols, method, window) &
            bind(C, name = "interpolate_gaps")

        integer(c_int), intent(in), value :: n_times, rows, cols, method, window
        real(c_double), intent(inout) :: arr(n_times * rows * cols)
        integer(c_int), intent(in) :: mask(n_times * rows * cols)

        integer :: r, c, base, npix
        npix = rows * cols

        do c = 1, cols
            do r = 1, rows
                base = (r - 1) * cols + (c - 1)
                select case (method)
                case (0); call fill_linear(arr, mask, n_times, npix, base)
                case (1); call fill_savgol(arr, mask, n_times, npix, base, window)
                case (2); call fill_pchip (arr, mask, n_times, npix, base)
                case (3); call fill_ets   (arr, mask, n_times, npix, base)
                case (4); call fill_gauss (arr, mask, n_times, npix, base, window)
                case default
                    call fill_linear(arr, mask, n_times, npix, base)
                end select
            end do
        end do
    end subroutine interpolate_gaps

    pure integer function ix(base, t, npix)
        integer, intent(in) :: base, t, npix
        ix = base + (t - 1) * npix + 1
    end function ix

    integer function valid_count(mask, n_times, npix, base)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base
        integer :: t
        valid_count = 0
        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) valid_count = valid_count + 1
        end do
    end function valid_count

    ! fill_edges_arr

    subroutine fill_edges_arr(arr, mask, n_times, npix, base)
        real(c_double), intent(inout) :: arr(*)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base
        integer :: t; real(c_double) :: v

        v = arr(ix(base, 1, npix))
        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) then; v = arr(ix(base, t, npix)); exit;
            end if
        end do
        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) exit
            arr(ix(base, t, npix)) = v
        end do

        v = arr(ix(base, n_times, npix))
        do t = n_times, 1, -1
            if (mask(ix(base, t, npix)) /= 0) then; v = arr(ix(base, t, npix)); exit;
            end if
        end do
        do t = n_times, 1, -1
            if (mask(ix(base, t, npix)) /= 0) exit
            arr(ix(base, t, npix)) = v
        end do
    end subroutine fill_edges_arr


    ! fill_edges_buf

    subroutine fill_edges_buf(buf, mask, n_times, npix, base)
        real(c_double), intent(inout) :: buf(n_times)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base
        integer :: t; real(c_double) :: v

        v = buf(1)
        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) then; v = buf(t); exit;
            end if
        end do
        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) exit
            buf(t) = v
        end do

        v = buf(n_times)
        do t = n_times, 1, -1
            if (mask(ix(base, t, npix)) /= 0) then; v = buf(t); exit;
            end if
        end do
        do t = n_times, 1, -1
            if (mask(ix(base, t, npix)) /= 0) exit
            buf(t) = v
        end do
    end subroutine fill_edges_buf


    ! writeback_gaps

    subroutine writeback_gaps(arr, mask, n_times, npix, base, out)
        real(c_double), intent(inout) :: arr(*)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base
        real(c_double), intent(in) :: out(n_times)
        integer :: t
        do t = 1, n_times
            if (mask(ix(base, t, npix)) == 0) arr(ix(base, t, npix)) = out(t)
        end do
    end subroutine writeback_gaps

    ! fill_linear_buf

    subroutine fill_linear_buf(buf, mask, n_times, npix, base)
        real(c_double), intent(inout) :: buf(n_times)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base
        integer :: t, lo, hi
        real(c_double) :: v_lo, v_hi, alpha

        call fill_edges_buf(buf, mask, n_times, npix, base)

        t = 1
        do while (t <= n_times)
            if (mask(ix(base, t, npix)) == 0) then
                lo = t - 1; hi = t
                do while (hi <= n_times)
                    if (mask(ix(base, hi, npix)) /= 0) exit; hi = hi + 1
                end do
                if (lo >= 1 .and. hi <= n_times) then
                    v_lo = buf(lo); v_hi = buf(hi)
                    do t = lo + 1, hi - 1
                        alpha = real(t - lo, c_double) / real(hi - lo, c_double)
                        buf(t) = v_lo + alpha * (v_hi - v_lo)
                    end do
                end if
                t = hi + 1
            else
                t = t + 1
            end if
        end do
    end subroutine fill_linear_buf

    ! METHOD 0 -- LINEAR

    subroutine fill_linear(arr, mask, n_times, npix, base)
        real(c_double), intent(inout) :: arr(*)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base
        integer :: t, lo, hi
        real(c_double) :: v_lo, v_hi, alpha

        if (valid_count(mask, n_times, npix, base) == 0) then
            do t = 1, n_times; arr(ix(base, t, npix)) = NODATA;
            end do; return
        end if

        call fill_edges_arr(arr, mask, n_times, npix, base)

        t = 1
        do while (t <= n_times)
            if (mask(ix(base, t, npix)) == 0) then
                lo = t - 1; hi = t
                do while (hi <= n_times)
                    if (mask(ix(base, hi, npix)) /= 0) exit; hi = hi + 1
                end do
                if (lo >= 1 .and. hi <= n_times) then
                    v_lo = arr(ix(base, lo, npix)); v_hi = arr(ix(base, hi, npix))
                    do t = lo + 1, hi - 1
                        alpha = real(t - lo, c_double) / real(hi - lo, c_double)
                        arr(ix(base, t, npix)) = v_lo + alpha * (v_hi - v_lo)
                    end do
                end if
                t = hi + 1
            else
                t = t + 1
            end if
        end do
    end subroutine fill_linear


    ! METHOD 1 -- SAVITZKY-GOLAY

    subroutine fill_savgol(arr, mask, n_times, npix, base, window)
        real(c_double), intent(inout) :: arr(*)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base, window
        integer :: t, hw, half, j, jj
        real(c_double), allocatable :: tmp(:), sg(:)
        real(c_double) :: s0, s1, s2, s3, s4, sy, sxy, sx2y, denom, a0, xj

        hw = max(window / 2, 1)
        if (2 * hw + 1 < 3) hw = 1

        if (valid_count(mask, n_times, npix, base) == 0) then
            do t = 1, n_times; arr(ix(base, t, npix)) = NODATA;
            end do; return
        end if

        allocate(tmp(n_times), sg(n_times))
        do t = 1, n_times; tmp(t) = arr(ix(base, t, npix));
        end do
        call fill_linear_buf(tmp, mask, n_times, npix, base)

        do t = 1, n_times
            half = min(hw, t - 1, n_times - t)
            s0 = 0.0d0;s1 = 0.0d0;s2 = 0.0d0;s3 = 0.0d0;s4 = 0.0d0
            sy = 0.0d0;sxy = 0.0d0;sx2y = 0.0d0
            do j = -half, half
                jj = t + j; xj = real(j, c_double)
                s0 = s0 + 1.0d0; s1 = s1 + xj; s2 = s2 + xj * xj
                s3 = s3 + xj**3; s4 = s4 + xj**4
                sy = sy + tmp(jj); sxy = sxy + xj * tmp(jj); sx2y = sx2y + xj * xj * tmp(jj)
            end do
            denom = s0 * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s3 * s2) + s2 * (s1 * s3 - s2 * s2)
            if (abs(denom) < 1.0d-12) then
                sg(t) = tmp(t)
            else
                a0 = (sy * (s2 * s4 - s3 * s3) - sxy * (s1 * s4 - s3 * s2) + sx2y * (s1 * s3 - s2 * s2)) / denom
                sg(t) = a0
            end if
        end do

        call writeback_gaps(arr, mask, n_times, npix, base, sg)
        deallocate(tmp, sg)
    end subroutine fill_savgol

    ! METHOD 2 -- PCHIP CUBIC HERMITE

    subroutine fill_pchip(arr, mask, n_times, npix, base)
        real(c_double), intent(inout) :: arr(*)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base
        integer :: t, k, nv, lo, hi
        real(c_double), allocatable :: xv(:), yv(:), dv(:), out(:)
        real(c_double) :: h, tt, tt2, tt3, h00, h10, h01, h11
        real(c_double) :: delta_k, alpha_fc, beta_fc, tau, w1, w2

        nv = valid_count(mask, n_times, npix, base)
        if (nv == 0) then
            do t = 1, n_times; arr(ix(base, t, npix)) = NODATA;
            end do; return
        end if
        if (nv == 1) then
            call fill_edges_arr(arr, mask, n_times, npix, base); return
        end if

        allocate(xv(nv), yv(nv), dv(nv), out(n_times))
        k = 0
        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) then
                k = k + 1; xv(k) = real(t, c_double); yv(k) = arr(ix(base, t, npix))
            end if
        end do

        do k = 1, nv - 1
            h = max(xv(k + 1) - xv(k), 1.0d-12)
            dv(k) = (yv(k + 1) - yv(k)) / h
        end do
        dv(nv) = dv(nv - 1)

        do k = 2, nv - 1
            w1 = 2.0d0 * (xv(k + 1) - xv(k)) + (xv(k) - xv(k - 1))
            w2 = (xv(k + 1) - xv(k)) + 2.0d0 * (xv(k) - xv(k - 1))
            dv(k) = (w1 * dv(k - 1) + w2 * dv(k)) / (w1 + w2)
        end do

        do k = 1, nv - 1
            delta_k = (yv(k + 1) - yv(k)) / max(abs(xv(k + 1) - xv(k)), 1.0d-12)
            if (abs(delta_k) < 1.0d-14) then
                dv(k) = 0.0d0; dv(k + 1) = 0.0d0
            else
                alpha_fc = dv(k) / delta_k
                beta_fc = dv(k + 1) / delta_k
                tau = alpha_fc * alpha_fc + beta_fc * beta_fc
                if (tau > 9.0d0) then
                    dv(k) = 3.0d0 * delta_k * alpha_fc / sqrt(tau)
                    dv(k + 1) = 3.0d0 * delta_k * beta_fc / sqrt(tau)
                end if
            end if
        end do

        do t = 1, n_times; out(t) = arr(ix(base, t, npix));
        end do

        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) cycle
            if (real(t, c_double) <= xv(1))  then; out(t) = yv(1);  cycle;
            end if
            if (real(t, c_double) >= xv(nv)) then; out(t) = yv(nv); cycle;
            end if

            lo = 1; hi = nv
            do while (hi - lo > 1)
                k = (lo + hi) / 2
                if (real(t, c_double) < xv(k)) then; hi = k;
                else; lo = k;
                end if
            end do

            h = xv(hi) - xv(lo)
            if (abs(h) < 1.0d-12) then; out(t) = yv(lo); cycle;
            end if
            tt = (real(t, c_double) - xv(lo)) / h
            tt2 = tt * tt; tt3 = tt2 * tt
            h00 = 2.0d0 * tt3 - 3.0d0 * tt2 + 1.0d0
            h10 = tt3 - 2.0d0 * tt2 + tt
            h01 = -2.0d0 * tt3 + 3.0d0 * tt2
            h11 = tt3 - tt2
            out(t) = h00 * yv(lo) + h10 * h * dv(lo) + h01 * yv(hi) + h11 * h * dv(hi)
        end do

        call fill_edges_buf(out, mask, n_times, npix, base)
        call writeback_gaps(arr, mask, n_times, npix, base, out)
        deallocate(xv, yv, dv, out)
    end subroutine fill_pchip

    ! METHOD 3 -- HOLT DOUBLE EXPONENTIAL SMOOTHING

    subroutine fill_ets(arr, mask, n_times, npix, base)
        real(c_double), intent(inout) :: arr(*)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base
        integer :: t, nv, k, ia, ib
        real(c_double), allocatable :: xv(:), yv(:), out(:)
        integer, allocatable :: tidx(:)
        real(c_double) :: alpha, beta, best_sse, sse, L, B, L_new
        real(c_double) :: a_try, b_try

        nv = valid_count(mask, n_times, npix, base)
        if (nv == 0) then
            do t = 1, n_times; arr(ix(base, t, npix)) = NODATA;
            end do; return
        end if
        if (nv <= 2) then
            call fill_linear(arr, mask, n_times, npix, base); return
        end if

        allocate(xv(nv), yv(nv), tidx(nv), out(n_times))
        k = 0
        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) then
                k = k + 1; tidx(k) = t
                xv(k) = real(t, c_double); yv(k) = arr(ix(base, t, npix))
            end if
        end do

        best_sse = huge(1.0d0); alpha = 0.3d0; beta = 0.1d0
        do ia = 1, 5
            do ib = 1, 5
                a_try = 0.1d0 + (ia - 1) * 0.2d0
                b_try = 0.05d0 + (ib - 1) * 0.1d0
                L = yv(1); B = (yv(2) - yv(1)) / max(xv(2) - xv(1), 1.0d0)
                sse = 0.0d0
                do k = 2, nv
                    L_new = a_try * yv(k) + (1.0d0 - a_try) * (L + B)
                    B = b_try * (L_new - L) + (1.0d0 - b_try) * B
                    L = L_new
                    sse = sse + (yv(k) - L)**2
                end do
                if (sse < best_sse) then; best_sse = sse; alpha = a_try; beta = b_try;
                end if
            end do
        end do

        do t = 1, n_times; out(t) = arr(ix(base, t, npix));
        end do
        L = yv(1); B = (yv(2) - yv(1)) / max(xv(2) - xv(1), 1.0d0)
        k = 1
        do t = 1, n_times
            if (mask(ix(base, t, npix)) /= 0) then
                if (k <= nv .and. k > 1) then
                    L_new = alpha * yv(k) + (1.0d0 - alpha) * (L + B)
                    B = beta * (L_new - L) + (1.0d0 - beta) * B
                    L = L_new
                end if
                if (k <= nv) k = k + 1
            else
                out(t) = L + B
            end if
        end do

        call fill_edges_buf(out, mask, n_times, npix, base)
        call writeback_gaps(arr, mask, n_times, npix, base, out)
        deallocate(xv, yv, tidx, out)
    end subroutine fill_ets

    ! METHOD 4 -- GAUSSIAN-WEIGHTED MOVING AVERAGE

    subroutine fill_gauss(arr, mask, n_times, npix, base, window)
        real(c_double), intent(inout) :: arr(*)
        integer(c_int), intent(in) :: mask(*)
        integer, intent(in) :: n_times, npix, base, window
        integer :: t, j, jj, hw, half
        real(c_double), allocatable :: tmp(:), gout(:)
        real(c_double) :: sigma, wsum, wval

        hw = max(window / 2, 1)
        if (2 * hw + 1 < 3) hw = 1
        sigma = real(2 * hw + 1, c_double) / 4.0d0

        if (valid_count(mask, n_times, npix, base) == 0) then
            do t = 1, n_times; arr(ix(base, t, npix)) = NODATA;
            end do; return
        end if

        allocate(tmp(n_times), gout(n_times))
        do t = 1, n_times; tmp(t) = arr(ix(base, t, npix));
        end do
        call fill_linear_buf(tmp, mask, n_times, npix, base)

        do t = 1, n_times
            half = min(hw, t - 1, n_times - t)
            wsum = 0.0d0; gout(t) = 0.0d0
            do j = -half, half
                jj = t + j
                wval = exp(-0.5d0 * (real(j, c_double) / sigma)**2)
                gout(t) = gout(t) + wval * tmp(jj); wsum = wsum + wval
            end do
            if (wsum > 1.0d-12) gout(t) = gout(t) / wsum
        end do

        call writeback_gaps(arr, mask, n_times, npix, base, gout)
        deallocate(tmp, gout)
    end subroutine fill_gauss

end module timeseries_mod