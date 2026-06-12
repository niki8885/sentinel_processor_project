module band_covariance_mod
    use iso_c_binding
    implicit none
    real(c_double), parameter :: NODATA = -9999.0d0

contains

    subroutine band_covariance(arr, rows, cols, n_bands, cov_out) &
            bind(C, name = "band_covariance")

        integer(c_int), intent(in), value :: rows, cols, n_bands
        real(c_double), intent(in) :: arr(n_bands * rows * cols)
        real(c_double), intent(out) :: cov_out(n_bands * n_bands)

        integer :: bi, bj, r, c, npix, idx_i, idx_j

        real(c_double), allocatable :: ksum_b(:), kcomp_b(:), mean_b(:)
        integer, allocatable :: cnt_b(:)

        real(c_double), allocatable :: cp_sum(:, :), cp_comp(:, :)
        integer, allocatable :: cp_cnt(:, :)

        real(c_double) :: val_i, val_j, di, dj, y, t

        npix = rows * cols

        allocate(ksum_b(n_bands), kcomp_b(n_bands), &
                mean_b(n_bands), cnt_b(n_bands))
        allocate(cp_sum(n_bands, n_bands), cp_comp(n_bands, n_bands), &
                cp_cnt(n_bands, n_bands))

        ksum_b = 0.0d0
        kcomp_b = 0.0d0
        cnt_b = 0
        cp_sum = 0.0d0
        cp_comp = 0.0d0
        cp_cnt = 0
        cov_out = 0.0d0

        do bi = 1, n_bands
            do r = 1, rows
                do c = 1, cols
                    idx_i = (bi - 1) * npix + (r - 1) * cols + c
                    val_i = arr(idx_i)
                    if (abs(val_i - NODATA) <= 1.0d-4) cycle

                    cnt_b(bi) = cnt_b(bi) + 1
                    y = val_i - kcomp_b(bi)
                    t = ksum_b(bi) + y
                    kcomp_b(bi) = (t - ksum_b(bi)) - y
                    ksum_b(bi) = t
                end do
            end do

            if (cnt_b(bi) > 0) then
                mean_b(bi) = ksum_b(bi) / real(cnt_b(bi), c_double)
            else
                mean_b(bi) = NODATA
            end if
        end do

        do bj = 1, n_bands
            do bi = 1, bj

                if (abs(mean_b(bi) - NODATA) <= 1.0d-4 .or. &
                        abs(mean_b(bj) - NODATA) <= 1.0d-4) cycle

                do r = 1, rows
                    do c = 1, cols
                        idx_i = (bi - 1) * npix + (r - 1) * cols + c
                        idx_j = (bj - 1) * npix + (r - 1) * cols + c
                        val_i = arr(idx_i)
                        val_j = arr(idx_j)

                        if (abs(val_i - NODATA) <= 1.0d-4) cycle
                        if (abs(val_j - NODATA) <= 1.0d-4) cycle

                        di = val_i - mean_b(bi)
                        dj = val_j - mean_b(bj)

                        cp_cnt(bi, bj) = cp_cnt(bi, bj) + 1

                        y = (di * dj) - cp_comp(bi, bj)
                        t = cp_sum(bi, bj) + y
                        cp_comp(bi, bj) = (t - cp_sum(bi, bj)) - y
                        cp_sum(bi, bj) = t
                    end do
                end do

            end do
        end do

        do bj = 1, n_bands
            do bi = 1, bj
                if (cp_cnt(bi, bj) > 1) then
                    cp_sum(bi, bj) = cp_sum(bi, bj) &
                            / real(cp_cnt(bi, bj) - 1, c_double)
                else
                    cp_sum(bi, bj) = NODATA
                end if

                cov_out((bi - 1) * n_bands + bj) = cp_sum(bi, bj)
                if (bi /= bj) then
                    cov_out((bj - 1) * n_bands + bi) = cp_sum(bi, bj)
                end if
            end do
        end do

        deallocate(ksum_b, kcomp_b, mean_b, cnt_b)
        deallocate(cp_sum, cp_comp, cp_cnt)

    end subroutine band_covariance

end module band_covariance_mod