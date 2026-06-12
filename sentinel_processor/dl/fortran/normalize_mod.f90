module normalize_mod
    use iso_c_binding
    implicit none

    real(c_double), parameter :: NODATA = -9999.0d0
    real(c_double), parameter :: NODATA_TOL = 1.0d-3

contains

    pure logical function is_nodata(v)
        real(c_double), intent(in) :: v
        is_nodata = (abs(v - NODATA) < NODATA_TOL)
    end function is_nodata

    ! minmax_band

    subroutine minmax_band(arr, n, out_min, out_max, out) &
            bind(C, name = "minmax_band")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: arr(n)
        real(c_double), intent(out) :: out_min, out_max
        real(c_double), intent(out) :: out(n)

        real(c_double) :: vmin, vmax, range
        integer :: i
        logical :: found_valid

        vmin = huge(1.0d0)
        vmax = -huge(1.0d0)
        found_valid = .false.

        do i = 1, n
            if (.not. is_nodata(arr(i))) then
                if (arr(i) < vmin) vmin = arr(i)
                if (arr(i) > vmax) vmax = arr(i)
                found_valid = .true.
            end if
        end do

        if (.not. found_valid) then
            out_min = 0.0d0; out_max = 0.0d0
            do i = 1, n; out(i) = NODATA;
            end do
            return
        end if

        out_min = vmin
        out_max = vmax
        range = vmax - vmin

        if (range < 1.0d-12) then
            do i = 1, n
                if (is_nodata(arr(i))) then; out(i) = NODATA
                else;                        out(i) = 0.5d0
                end if
            end do
        else
            do i = 1, n
                if (is_nodata(arr(i))) then
                    out(i) = NODATA
                else
                    out(i) = (arr(i) - vmin) / range
                end if
            end do
        end if
    end subroutine minmax_band

    ! zscore_band

    subroutine zscore_band(arr, n, mean, std, out) &
            bind(C, name = "zscore_band")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in), value :: mean, std
        real(c_double), intent(in) :: arr(n)
        real(c_double), intent(out) :: out(n)

        integer :: i

        if (std < 1.0d-12) then
            do i = 1, n
                if (is_nodata(arr(i))) then; out(i) = NODATA
                else;                        out(i) = 0.0d0
                end if
            end do
            return
        end if

        do i = 1, n
            if (is_nodata(arr(i))) then
                out(i) = NODATA
            else
                out(i) = (arr(i) - mean) / std
            end if
        end do
    end subroutine zscore_band

    ! band_stats

    subroutine band_stats(arr, n, out_mean, out_std, out_min, out_max) &
            bind(C, name = "band_stats")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: arr(n)
        real(c_double), intent(out) :: out_mean, out_std, out_min, out_max

        integer :: i, cnt
        real(c_double) :: delta, mean_acc, m2_acc, vmin, vmax

        cnt = 0
        mean_acc = 0.0d0
        m2_acc = 0.0d0
        vmin = huge(1.0d0)
        vmax = -huge(1.0d0)

        do i = 1, n
            if (is_nodata(arr(i))) cycle
            cnt = cnt + 1
            delta = arr(i) - mean_acc
            mean_acc = mean_acc + delta / real(cnt, c_double)
            m2_acc = m2_acc + delta * (arr(i) - mean_acc)
            if (arr(i) < vmin) vmin = arr(i)
            if (arr(i) > vmax) vmax = arr(i)
        end do

        if (cnt == 0) then
            out_mean = 0.0d0; out_std = 0.0d0
            out_min = 0.0d0; out_max = 0.0d0
            return
        end if

        out_mean = mean_acc
        out_std = sqrt(max(m2_acc / real(cnt, c_double), 0.0d0))
        out_min = vmin
        out_max = vmax
    end subroutine band_stats

end module normalize_mod