module wavelet_mod
    use iso_c_binding
    implicit none

    integer, parameter :: HAAR_LEN = 2
    real(c_double), parameter :: HAAR_LO(2) = &
            (/  0.7071067811865476d0, 0.7071067811865476d0 /)
    real(c_double), parameter :: HAAR_HI(2) = &
            (/  0.7071067811865476d0, -0.7071067811865476d0 /)

    integer, parameter :: DB4_LEN = 4
    real(c_double), parameter :: DB4_LO(4) = (/ &
            0.4829629131445341d0, 0.8365163037378079d0, &
                    0.2241438680420134d0, -0.1294095225512604d0 /)
    real(c_double), parameter :: DB4_HI(4) = (/ &
            -0.1294095225512604d0, -0.2241438680420134d0, &
                    0.8365163037378079d0, -0.4829629131445341d0 /)

    integer, parameter :: DB6_LEN = 6
    real(c_double), parameter :: DB6_LO(6) = (/ &
            0.3326705529500825d0, 0.8068915093110924d0, &
                    0.4598775021184915d0, -0.1350110200102546d0, &
                    -0.0854412738820267d0, 0.0352262918857095d0 /)
    real(c_double), parameter :: DB6_HI(6) = (/ &
            0.0352262918857095d0, 0.0854412738820267d0, &
                    -0.1350110200102546d0, -0.4598775021184915d0, &
                    0.8068915093110924d0, -0.3326705529500825d0 /)

    integer, parameter :: COIF1_LEN = 6
    real(c_double), parameter :: COIF1_LO(6) = (/ &
            -0.0156557281357920d0, -0.0727326195125265d0, &
                    0.3848648468648578d0, 0.8525720202116004d0, &
                    0.3378976624574818d0, -0.0727326195125265d0 /)
    real(c_double), parameter :: COIF1_HI(6) = (/ &
            0.0727326195125265d0, 0.3378976624574818d0, &
                    -0.8525720202116004d0, 0.3848648468648578d0, &
                    0.0727326195125265d0, -0.0156557281357920d0 /)

    integer, parameter :: SYM4_LEN = 8
    real(c_double), parameter :: SYM4_LO(8) = (/ &
            -0.0757657147893551d0, -0.0296355276459541d0, &
                    0.4976186676324578d0, 0.8037387518052163d0, &
                    0.2978577956052774d0, -0.0992195435769354d0, &
                    -0.0126039672622612d0, 0.0322231006040427d0 /)
    real(c_double), parameter :: SYM4_HI(8) = (/ &
            0.0322231006040427d0, 0.0126039672622612d0, &
                    -0.0992195435769354d0, -0.2978577956052774d0, &
                    0.8037387518052163d0, -0.4976186676324578d0, &
                    -0.0296355276459541d0, 0.0757657147893551d0 /)

    integer, parameter :: SYM6_LEN = 12
    real(c_double), parameter :: SYM6_LO(12) = (/ &
            0.0154041093270274d0, 0.0034907120842175d0, &
                    -0.1179901111481906d0, -0.0483117425856330d0, &
                    0.4910559419267466d0, 0.7876411410301940d0, &
                    0.3379294217276218d0, -0.0726375227864625d0, &
                    -0.0210602925123006d0, 0.0447249017706658d0, &
                    0.0017677118642428d0, -0.0078007083250341d0 /)
    real(c_double), parameter :: SYM6_HI(12) = (/ &
            0.0078007083250341d0, 0.0017677118642428d0, &
                    -0.0447249017706658d0, -0.0210602925123006d0, &
                    0.0726375227864625d0, 0.3379294217276218d0, &
                    -0.7876411410301940d0, 0.4910559419267466d0, &
                    0.0483117425856330d0, -0.1179901111481906d0, &
                    -0.0034907120842175d0, 0.0154041093270274d0 /)

    real(c_double), parameter :: ZERO = 0.0d0
    real(c_double), parameter :: HALF = 0.5d0

contains

    pure integer function pmod(i, n)
        integer, intent(in) :: i, n
        integer :: tmp
        tmp = mod(i - 1, n)
        if (tmp < 0) tmp = tmp + n
        pmod = tmp + 1
    end function pmod


    subroutine select_filters(wavelet, lo, hi, flen)
        integer, intent(in) :: wavelet
        real(c_double), allocatable, intent(out) :: lo(:), hi(:)
        integer, intent(out) :: flen
        select case (wavelet)
        case (1); flen = DB4_LEN;   allocate(lo(flen), hi(flen)); lo = DB4_LO;   hi = DB4_HI
        case (2); flen = DB6_LEN;   allocate(lo(flen), hi(flen)); lo = DB6_LO;   hi = DB6_HI
        case (3); flen = COIF1_LEN; allocate(lo(flen), hi(flen)); lo = COIF1_LO; hi = COIF1_HI
        case (4); flen = SYM4_LEN;  allocate(lo(flen), hi(flen)); lo = SYM4_LO;  hi = SYM4_HI
        case (5); flen = SYM6_LEN;  allocate(lo(flen), hi(flen)); lo = SYM6_LO;  hi = SYM6_HI
        case default  ! 0 = Haar
            flen = HAAR_LEN;  allocate(lo(flen), hi(flen)); lo = HAAR_LO;  hi = HAAR_HI
        end select
    end subroutine select_filters


    subroutine dwt1d(x, n, lo, hi, flen, lo_out, hi_out, nout)
        integer, intent(in) :: n, flen, nout
        real(c_double), intent(in) :: x(n), lo(flen), hi(flen)
        real(c_double), intent(out) :: lo_out(nout), hi_out(nout)
        integer :: k, j, idx
        real(c_double) :: s_lo, s_hi
        do k = 1, nout
            s_lo = ZERO;  s_hi = ZERO
            do j = 1, flen
                idx = pmod(2 * k - 1 + j - 1, n)
                s_lo = s_lo + lo(j) * x(idx)
                s_hi = s_hi + hi(j) * x(idx)
            end do
            lo_out(k) = s_lo
            hi_out(k) = s_hi
        end do
    end subroutine dwt1d


    subroutine idwt1d(lo_in, hi_in, nout, lo, hi, flen, y, n)
        integer, intent(in) :: nout, flen, n
        real(c_double), intent(in) :: lo_in(nout), hi_in(nout)
        real(c_double), intent(in) :: lo(flen), hi(flen)
        real(c_double), intent(out) :: y(n)
        integer :: k, j, idx
        y(1:n) = ZERO
        do k = 1, nout
            do j = 1, flen
                idx = pmod(2 * k - 1 + j - 1, n)
                y(idx) = y(idx) + lo(j) * lo_in(k) + hi(j) * hi_in(k)
            end do
        end do
    end subroutine idwt1d

    subroutine dwt2d_one_level(buf, rT, cT, lo, hi, flen, out)
        integer, intent(in) :: rT, cT, flen
        real(c_double), intent(in) :: buf(rT * cT), lo(flen), hi(flen)
        real(c_double), intent(out) :: out(rT * cT)
        integer :: rh, ch, r, c, cc, rr
        real(c_double), allocatable :: row_lo(:), row_hi(:), tmp(:)
        real(c_double) :: rowvec(cT)

        rh = rT / 2;  ch = cT / 2
        allocate(row_lo(ch), row_hi(ch), tmp(rT * cT))

        do r = 1, rT
            do cc = 1, cT
                rowvec(cc) = buf((cc - 1) * rT + r)
            end do
            call dwt1d(rowvec, cT, lo, hi, flen, row_lo, row_hi, ch)
            do cc = 1, ch
                tmp((cc - 1) * rT + r) = row_lo(cc)
                tmp((ch + cc - 1) * rT + r) = row_hi(cc)
            end do
        end do

        do c = 1, cT
            block
                real(c_double) :: colvec(rT), col_lo(rh), col_hi(rh)
                do rr = 1, rT
                    colvec(rr) = tmp((c - 1) * rT + rr)
                end do
                call dwt1d(colvec, rT, lo, hi, flen, col_lo, col_hi, rh)
                do rr = 1, rh
                    out((c - 1) * rT + rr) = col_lo(rr)
                    out((c - 1) * rT + rh + rr) = col_hi(rr)
                end do
            end block
        end do

        deallocate(row_lo, row_hi, tmp)
    end subroutine dwt2d_one_level

    subroutine idwt2d_one_level(buf, rT, cT, lo, hi, flen, out)
        integer, intent(in) :: rT, cT, flen
        real(c_double), intent(in) :: buf(rT * cT), lo(flen), hi(flen)
        real(c_double), intent(out) :: out(rT * cT)
        integer :: rh, ch, r, c, cc, rr
        real(c_double), allocatable :: tmp(:)

        rh = rT / 2;  ch = cT / 2
        allocate(tmp(rT * cT))
        tmp(1:rT * cT) = ZERO

        do c = 1, cT
            block
                real(c_double) :: col_lo(rh), col_hi(rh), colvec(rT)
                do rr = 1, rh
                    col_lo(rr) = buf((c - 1) * rT + rr)
                    col_hi(rr) = buf((c - 1) * rT + rh + rr)
                end do
                call idwt1d(col_lo, col_hi, rh, lo, hi, flen, colvec, rT)
                do rr = 1, rT
                    tmp((c - 1) * rT + rr) = colvec(rr)
                end do
            end block
        end do

        out(1:rT * cT) = ZERO
        do r = 1, rT
            block
                real(c_double) :: row_lo(ch), row_hi(ch), rowvec(cT)
                do cc = 1, ch
                    row_lo(cc) = tmp((cc - 1) * rT + r)
                    row_hi(cc) = tmp((ch + cc - 1) * rT + r)
                end do
                call idwt1d(row_lo, row_hi, ch, lo, hi, flen, rowvec, cT)
                do cc = 1, cT
                    out((cc - 1) * rT + r) = rowvec(cc)
                end do
            end block
        end do

        deallocate(tmp)
    end subroutine idwt2d_one_level


    subroutine extract_tile(full_arr, full_rows, rT, cT, tile)
        integer, intent(in) :: full_rows, rT, cT
        real(c_double), intent(in) :: full_arr(*)
        real(c_double), intent(out) :: tile(rT * cT)
        integer :: r, c
        do c = 1, cT
            do r = 1, rT
                tile((c - 1) * rT + r) = full_arr((c - 1) * full_rows + r)
            end do
        end do
    end subroutine extract_tile

    subroutine insert_tile(tile, full_rows, rT, cT, full_arr)
        integer, intent(in) :: full_rows, rT, cT
        real(c_double), intent(in) :: tile(rT * cT)
        real(c_double), intent(inout) :: full_arr(*)
        integer :: r, c
        do c = 1, cT
            do r = 1, rT
                full_arr((c - 1) * full_rows + r) = tile((c - 1) * rT + r)
            end do
        end do
    end subroutine insert_tile


    subroutine dwt2d(arr, rows, cols, levels, wavelet, coeffs_out) &
            bind(C, name = "dwt2d")

        integer(c_int), intent(in), value :: rows, cols, levels, wavelet
        real(c_double), intent(in) :: arr(rows * cols)
        real(c_double), intent(out) :: coeffs_out(rows * cols)

        integer :: lv, rT, cT, flen
        real(c_double), allocatable :: lo(:), hi(:), tile_in(:), tile_out(:)

        call select_filters(wavelet, lo, hi, flen)
        coeffs_out(1:rows * cols) = arr(1:rows * cols)
        rT = rows;  cT = cols
        do lv = 1, levels
            allocate(tile_in(rT * cT), tile_out(rT * cT))
            call extract_tile(coeffs_out, rows, rT, cT, tile_in)
            call dwt2d_one_level(tile_in, rT, cT, lo, hi, flen, tile_out)
            call insert_tile(tile_out, rows, rT, cT, coeffs_out)
            deallocate(tile_in, tile_out)
            rT = rT / 2;  cT = cT / 2
        end do
        deallocate(lo, hi)
    end subroutine dwt2d

    subroutine idwt2d(coeffs, rows, cols, levels, wavelet, arr_out) &
            bind(C, name = "idwt2d")

        integer(c_int), intent(in), value :: rows, cols, levels, wavelet
        real(c_double), intent(in) :: coeffs(rows * cols)
        real(c_double), intent(out) :: arr_out(rows * cols)

        integer :: lv, rT, cT, flen
        real(c_double), allocatable :: lo(:), hi(:), tile_in(:), tile_out(:)

        call select_filters(wavelet, lo, hi, flen)
        arr_out(1:rows * cols) = coeffs(1:rows * cols)
        rT = rows;  cT = cols
        do lv = 1, levels - 1
            rT = rT / 2;  cT = cT / 2
        end do
        do lv = levels, 1, -1
            allocate(tile_in(rT * cT), tile_out(rT * cT))
            call extract_tile(arr_out, rows, rT, cT, tile_in)
            call idwt2d_one_level(tile_in, rT, cT, lo, hi, flen, tile_out)
            call insert_tile(tile_out, rows, rT, cT, arr_out)
            deallocate(tile_in, tile_out)
            rT = rT * 2;  cT = cT * 2
        end do
        deallocate(lo, hi)
    end subroutine idwt2d

    subroutine estimate_sigma(coeffs, n, sigma_out) &
            bind(C, name = "estimate_sigma")

        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: coeffs(n)
        real(c_double), intent(out) :: sigma_out

        real(c_double), allocatable :: tmp(:)
        real(c_double) :: med, t
        integer :: i, j

        allocate(tmp(n))
        tmp(1:n) = abs(coeffs(1:n))

        do i = 2, n
            t = tmp(i);  j = i - 1
            do while (j >= 1 .and. tmp(j) > t)
                tmp(j + 1) = tmp(j);  j = j - 1
            end do
            tmp(j + 1) = t
        end do

        if (mod(n, 2) == 0) then
            med = HALF * (tmp(n / 2) + tmp(n / 2 + 1))
        else
            med = tmp((n + 1) / 2)
        end if

        sigma_out = med / 0.6744897501960817d0
        deallocate(tmp)
    end subroutine estimate_sigma


    subroutine bayes_threshold(coeffs, n, sigma_n, threshold_out) &
            bind(C, name = "bayes_threshold")

        integer(c_int), intent(in), value :: n
        real(c_double), intent(in), value :: sigma_n
        real(c_double), intent(in) :: coeffs(n)
        real(c_double), intent(out) :: threshold_out

        real(c_double) :: var_y, var_s, sigma_s, max_abs
        integer :: i

        var_y = ZERO;  max_abs = ZERO
        do i = 1, n
            var_y = var_y + coeffs(i) * coeffs(i)
            if (abs(coeffs(i)) > max_abs) max_abs = abs(coeffs(i))
        end do
        var_y = var_y / real(n, c_double)

        var_s = max(var_y - sigma_n * sigma_n, ZERO)
        sigma_s = sqrt(var_s)

        if (sigma_s < 1.0d-12) then
            threshold_out = max_abs
        else
            threshold_out = sigma_n * sigma_n / sigma_s
        end if
    end subroutine bayes_threshold


    subroutine band_energy(coeffs, n, energy_out) &
            bind(C, name = "band_energy")

        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: coeffs(n)
        real(c_double), intent(out) :: energy_out

        integer :: i
        energy_out = ZERO
        do i = 1, n
            energy_out = energy_out + coeffs(i) * coeffs(i)
        end do
    end subroutine band_energy


    subroutine band_stats(coeffs, n, mean_out, var_out, l1_out, linf_out) &
            bind(C, name = "band_stats")

        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: coeffs(n)
        real(c_double), intent(out) :: mean_out, var_out, l1_out, linf_out

        integer :: i
        real(c_double) :: rn, diff

        rn = real(n, c_double)
        mean_out = ZERO;  l1_out = ZERO;  linf_out = ZERO

        do i = 1, n
            mean_out = mean_out + coeffs(i)
            l1_out = l1_out + abs(coeffs(i))
            if (abs(coeffs(i)) > linf_out) linf_out = abs(coeffs(i))
        end do
        mean_out = mean_out / rn
        l1_out = l1_out / rn

        var_out = ZERO
        do i = 1, n
            diff = coeffs(i) - mean_out
            var_out = var_out + diff * diff
        end do
        var_out = var_out / rn
    end subroutine band_stats


    subroutine dwt2d_batch(stack, rows, cols, n_bands, levels, wavelet, &
            coeffs_out) &
            bind(C, name = "dwt2d_batch")

        integer(c_int), intent(in), value :: rows, cols, n_bands, levels, wavelet
        real(c_double), intent(in) :: stack(rows * cols * n_bands)
        real(c_double), intent(out) :: coeffs_out(rows * cols * n_bands)

        integer :: b, npix, offset
        real(c_double), allocatable :: band_in(:), band_out(:)

        npix = rows * cols
        allocate(band_in(npix), band_out(npix))

        do b = 0, n_bands - 1
            offset = b * npix
            band_in(1:npix) = stack(offset + 1:offset + npix)
            call dwt2d(band_in, rows, cols, levels, wavelet, band_out)
            coeffs_out(offset + 1:offset + npix) = band_out(1:npix)
        end do

        deallocate(band_in, band_out)
    end subroutine dwt2d_batch


    subroutine idwt2d_batch(coeffs, rows, cols, n_bands, levels, wavelet, &
            arr_out) &
            bind(C, name = "idwt2d_batch")

        integer(c_int), intent(in), value :: rows, cols, n_bands, levels, wavelet
        real(c_double), intent(in) :: coeffs(rows * cols * n_bands)
        real(c_double), intent(out) :: arr_out(rows * cols * n_bands)

        integer :: b, npix, offset
        real(c_double), allocatable :: band_in(:), band_out(:)

        npix = rows * cols
        allocate(band_in(npix), band_out(npix))

        do b = 0, n_bands - 1
            offset = b * npix
            band_in(1:npix) = coeffs(offset + 1:offset + npix)
            call idwt2d(band_in, rows, cols, levels, wavelet, band_out)
            arr_out(offset + 1:offset + npix) = band_out(1:npix)
        end do

        deallocate(band_in, band_out)
    end subroutine idwt2d_batch


    subroutine dwt3d(arr, rows, cols, n_times, levels, wavelet, coeffs_out) &
            bind(C, name = "dwt3d")

        integer(c_int), intent(in), value :: rows, cols, n_times, levels, wavelet
        real(c_double), intent(in) :: arr(rows * cols * n_times)
        real(c_double), intent(out) :: coeffs_out(rows * cols * n_times)

        integer :: t, p, npix, lv, nt_work, flen
        real(c_double), allocatable :: lo(:), hi(:)
        real(c_double), allocatable :: band_in(:), band_out(:)
        real(c_double), allocatable :: tvec(:), tlo(:), thi(:)

        npix = rows * cols
        allocate(band_in(npix), band_out(npix))

        do t = 1, n_times
            band_in(1:npix) = arr((t - 1) * npix + 1:t * npix)
            call dwt2d(band_in, rows, cols, levels, wavelet, band_out)
            coeffs_out((t - 1) * npix + 1:t * npix) = band_out(1:npix)
        end do
        deallocate(band_in, band_out)

        call select_filters(wavelet, lo, hi, flen)
        nt_work = n_times

        do lv = 1, levels
            allocate(tlo(nt_work / 2), thi(nt_work / 2), tvec(nt_work))
            do p = 1, npix
                do t = 1, nt_work
                    tvec(t) = coeffs_out((t - 1) * npix + p)
                end do
                call dwt1d(tvec, nt_work, lo, hi, flen, tlo, thi, nt_work / 2)
                do t = 1, nt_work / 2
                    coeffs_out((t - 1) * npix + p) = tlo(t)
                    coeffs_out((t + nt_work / 2 - 1) * npix + p) = thi(t)
                end do
            end do
            deallocate(tlo, thi, tvec)
            nt_work = nt_work / 2
        end do
        deallocate(lo, hi)
    end subroutine dwt3d


    subroutine idwt3d(coeffs, rows, cols, n_times, levels, wavelet, arr_out) &
            bind(C, name = "idwt3d")

        integer(c_int), intent(in), value :: rows, cols, n_times, levels, wavelet
        real(c_double), intent(in) :: coeffs(rows * cols * n_times)
        real(c_double), intent(out) :: arr_out(rows * cols * n_times)

        integer :: t, p, npix, lv, nt_work, flen
        real(c_double), allocatable :: lo(:), hi(:)
        real(c_double), allocatable :: band_in(:), band_out(:)
        real(c_double), allocatable :: tvec(:), tlo(:), thi(:)

        npix = rows * cols
        arr_out(1:npix * n_times) = coeffs(1:npix * n_times)

        call select_filters(wavelet, lo, hi, flen)
        nt_work = n_times
        do lv = 1, levels - 1
            nt_work = nt_work / 2
        end do

        do lv = levels, 1, -1
            allocate(tlo(nt_work / 2), thi(nt_work / 2), tvec(nt_work))
            do p = 1, npix
                do t = 1, nt_work / 2
                    tlo(t) = arr_out((t - 1) * npix + p)
                    thi(t) = arr_out((t + nt_work / 2 - 1) * npix + p)
                end do
                call idwt1d(tlo, thi, nt_work / 2, lo, hi, flen, tvec, nt_work)
                do t = 1, nt_work
                    arr_out((t - 1) * npix + p) = tvec(t)
                end do
            end do
            deallocate(tlo, thi, tvec)
            nt_work = nt_work * 2
        end do
        deallocate(lo, hi)

        allocate(band_in(npix), band_out(npix))
        do t = 1, n_times
            band_in(1:npix) = arr_out((t - 1) * npix + 1:t * npix)
            call idwt2d(band_in, rows, cols, levels, wavelet, band_out)
            arr_out((t - 1) * npix + 1:t * npix) = band_out(1:npix)
        end do
        deallocate(band_in, band_out)
    end subroutine idwt3d

end module wavelet_mod