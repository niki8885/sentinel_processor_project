module texture_mod
    use iso_c_binding, only : c_double, c_int
    implicit none
    private

    integer, parameter :: N_LEVELS = 64
    real(c_double), parameter :: NODATA = -9999.0d0

    public :: compute_glcm

contains

    ! quantise_band

    subroutine quantise_band(arr, rows, cols, qlev)
        integer(c_int), intent(in) :: rows, cols
        real(c_double), intent(in) :: arr(rows, cols)
        integer, intent(out) :: qlev(rows, cols)

        real(c_double) :: bmin, bmax, brange, inv_range
        integer :: r, c, lev

        bmin = huge(bmin)
        bmax = -huge(bmax)
        do c = 1, cols
            do r = 1, rows
                if (arr(r, c) > NODATA + 1.0d0) then
                    if (arr(r, c) < bmin) bmin = arr(r, c)
                    if (arr(r, c) > bmax) bmax = arr(r, c)
                end if
            end do
        end do

        brange = bmax - bmin
        if (brange < 1.0d-12) brange = 1.0d-12
        inv_range = real(N_LEVELS - 1, c_double) / brange

        do c = 1, cols
            do r = 1, rows
                if (arr(r, c) <= NODATA + 1.0d0) then
                    qlev(r, c) = 0
                else
                    lev = int((arr(r, c) - bmin) * inv_range) + 1
                    if (lev < 1)        lev = 1
                    if (lev > N_LEVELS) lev = N_LEVELS
                    qlev(r, c) = lev
                end if
            end do
        end do
    end subroutine quantise_band


    subroutine angle_offsets(angle_deg, distance, dr, dc)
        integer(c_int), intent(in) :: angle_deg, distance
        integer, intent(out) :: dr, dc

        select case (angle_deg)
        case (0)
            dr = 0;  dc = distance
        case (45)
            dr = -distance; dc = distance
        case (90)
            dr = -distance; dc = 0
        case (135)
            dr = -distance; dc = -distance
        case default
            dr = 0;  dc = distance
        end select
    end subroutine angle_offsets


    subroutine glcm_features_single_angle(&
            qlev, rows, cols, &
            r0, c0, half_w, &
            dr, dc, &
            energy, contrast, homogeneity)

        integer(c_int), intent(in) :: rows, cols
        integer, intent(in) :: qlev(rows, cols)
        integer, intent(in) :: r0, c0, half_w, dr, dc
        real(c_double), intent(out) :: energy, contrast, homogeneity

        real(c_double) :: glcm(N_LEVELS, N_LEVELS)
        real(c_double) :: norm, p, diff2
        integer :: r, c, r1, c1, i, j, n_pairs, di, dj

        glcm = 0.0d0
        n_pairs = 0

        do c = c0 - half_w, c0 + half_w
            do r = r0 - half_w, r0 + half_w
                r1 = r + dr
                c1 = c + dc
                if (r  < 1 .or. r  > rows) cycle
                if (c  < 1 .or. c  > cols) cycle
                if (r1 < 1 .or. r1 > rows) cycle
                if (c1 < 1 .or. c1 > cols) cycle
                if (qlev(r, c)  == 0) cycle
                if (qlev(r1, c1) == 0) cycle

                i = qlev(r, c)
                j = qlev(r1, c1)
                glcm(i, j) = glcm(i, j) + 1.0d0
                glcm(j, i) = glcm(j, i) + 1.0d0   ! symmetrise
                n_pairs = n_pairs + 1
            end do
        end do

        norm = sum(glcm)
        if (norm < 1.0d-12) then
            energy = 0.0d0
            contrast = 0.0d0
            homogeneity = 0.0d0
            return
        end if
        glcm = glcm / norm

        energy = 0.0d0
        contrast = 0.0d0
        homogeneity = 0.0d0

        do j = 1, N_LEVELS
            do i = 1, N_LEVELS
                p = glcm(i, j)
                di = i - j
                diff2 = real(di * di, c_double)
                energy = energy + p * p
                contrast = contrast + diff2 * p
                homogeneity = homogeneity + p / (1.0d0 + diff2)
            end do
        end do
    end subroutine glcm_features_single_angle


    ! compute_glcm  (C-callable entry point)
    !
    ! arr            – (rows, cols) band, Fortran column-major order
    ! rows, cols     – raster dimensions
    ! window         – neighbourhood side length (e.g. 7 → 7×7 patch)
    ! distance       – pixel offset for co-occurrence (usually 1)
    ! angle_deg      – 0 | 45 | 90 | 135 ; -1 = isotropic average
    ! energy_out     – (rows, cols) output, filled with NODATA at borders
    ! contrast_out   – (rows, cols) output
    ! homogeneity_out– (rows, cols) output


    subroutine compute_glcm(&
            arr, rows, cols, &
            window, distance, angle_deg, &
            energy_out, contrast_out, homogeneity_out) &
            bind(C, name = "compute_glcm")

        integer(c_int), value, intent(in) :: rows, cols
        integer(c_int), value, intent(in) :: window, distance, angle_deg
        real(c_double), intent(in) :: arr(rows, cols)
        real(c_double), intent(out) :: energy_out(rows, cols)
        real(c_double), intent(out) :: contrast_out(rows, cols)
        real(c_double), intent(out) :: homogeneity_out(rows, cols)

        integer :: qlev(rows, cols)
        integer :: half_w, border, r, c
        integer :: dr, dc
        real(c_double) :: en, co, ho
        real(c_double) :: en_sum, co_sum, ho_sum

        integer, parameter :: N_ANGLES = 4
        integer :: angles(N_ANGLES)
        data angles / 0, 45, 90, 135 /
        integer :: a_idx

        half_w = window / 2

        border = half_w + distance

        ! Step 1: quantise to N_LEVELS grey levels

        call quantise_band(arr, rows, cols, qlev)

        ! Step 2: per-pixel GLCM

        do c = 1, cols
            do r = 1, rows

                ! Fill border pixels with NODATA
                if (r <= border .or. r > rows - border .or. &
                        c <= border .or. c > cols - border) then
                    energy_out(r, c) = NODATA
                    contrast_out(r, c) = NODATA
                    homogeneity_out(r, c) = NODATA
                    cycle
                end if

                if (angle_deg == -1) then
                    en_sum = 0.0d0
                    co_sum = 0.0d0
                    ho_sum = 0.0d0
                    do a_idx = 1, N_ANGLES
                        call angle_offsets(&
                                int(angles(a_idx), c_int), &
                                int(distance, c_int), dr, dc)
                        call glcm_features_single_angle(&
                                qlev, rows, cols, &
                                r, c, half_w, dr, dc, &
                                en, co, ho)
                        en_sum = en_sum + en
                        co_sum = co_sum + co
                        ho_sum = ho_sum + ho
                    end do
                    energy_out(r, c) = en_sum / real(N_ANGLES, c_double)
                    contrast_out(r, c) = co_sum / real(N_ANGLES, c_double)
                    homogeneity_out(r, c) = ho_sum / real(N_ANGLES, c_double)

                else
                    call angle_offsets(angle_deg, distance, dr, dc)
                    call glcm_features_single_angle(&
                            qlev, rows, cols, &
                            r, c, half_w, dr, dc, &
                            en, co, ho)
                    energy_out(r, c) = en
                    contrast_out(r, c) = co
                    homogeneity_out(r, c) = ho
                end if

            end do
        end do

    end subroutine compute_glcm

end module texture_mod