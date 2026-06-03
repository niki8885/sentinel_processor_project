module filters_mod
  use iso_c_binding
  implicit none

  real(c_double), parameter :: PI = 3.14159265358979323846d0

contains

  ! GENERIC 2-D CONVOLUTION

  subroutine convolve2d(src, rows, cols, kernel, krows, kcols, dst) &
      bind(C, name="convolve2d")

    integer(c_int), intent(in), value :: rows, cols, krows, kcols
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(in)        :: kernel(krows * kcols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: r, c, kr, kc, sr, sc
    integer        :: hr, hc
    real(c_double) :: acc

    hr = krows / 2
    hc = kcols / 2

    do c = 1, cols
      do r = 1, rows
        acc = 0.0d0
        do kc = 1, kcols
          sc = c + (kc - 1 - hc)
          if (sc < 1 .or. sc > cols) cycle
          do kr = 1, krows
            sr = r + (kr - 1 - hr)
            if (sr < 1 .or. sr > rows) cycle
            acc = acc + kernel((kc-1)*krows + kr) * src((sc-1)*rows + sr)
          end do
        end do
        dst((c-1)*rows + r) = acc
      end do
    end do

  end subroutine convolve2d


  ! SEPARABLE GAUSSIAN BLUR

  subroutine gaussian_blur(src, rows, cols, sigma, kradius, dst) &
      bind(C, name="gaussian_blur")

    integer(c_int), intent(in), value :: rows, cols, kradius
    real(c_double), intent(in), value :: sigma
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: ksize
    integer        :: r, c, k, sc, sr
    real(c_double) :: acc, wsum, w
    real(c_double), allocatable :: kernel(:), tmp(:)

    ksize = 2 * kradius + 1
    allocate(kernel(ksize), tmp(rows * cols))

    wsum = 0.0d0
    do k = 1, ksize
      w = exp(-0.5d0 * real(k - 1 - kradius, c_double)**2 / (sigma * sigma))
      kernel(k) = w
      wsum = wsum + w
    end do
    kernel = kernel / wsum   ! normalise

    do c = 1, cols
      do r = 1, rows
        acc = 0.0d0
        do k = 1, ksize
          sc = c + (k - 1 - kradius)
          sc = max(1, min(sc, cols))
          acc = acc + kernel(k) * src((sc-1)*rows + r)
        end do
        tmp((c-1)*rows + r) = acc
      end do
    end do

    do c = 1, cols
      do r = 1, rows
        acc = 0.0d0
        do k = 1, ksize
          sr = r + (k - 1 - kradius)
          sr = max(1, min(sr, rows))
          acc = acc + kernel(k) * tmp((c-1)*rows + sr)
        end do
        dst((c-1)*rows + r) = acc
      end do
    end do

    deallocate(kernel, tmp)

  end subroutine gaussian_blur


  ! SOBEL EDGE MAGNITUDE
  subroutine sobel_magnitude(src, rows, cols, use_l2, dst) &
      bind(C, name="sobel_magnitude")

    integer(c_int), intent(in), value :: rows, cols, use_l2
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: r, c, ri, ci
    real(c_double) :: gx, gy
    real(c_double) :: px(-1:1, -1:1)

    do c = 1, cols
      do r = 1, rows
        do ci = -1, 1
          do ri = -1, 1
            px(ri, ci) = src( &
              (max(1,min(c+ci,cols))-1)*rows + max(1,min(r+ri,rows)) )
          end do
        end do

        gx = -px(-1,-1) + px(-1,1) &
             - 2.0d0*px(0,-1) + 2.0d0*px(0,1) &
             - px(1,-1) + px(1,1)

        gy = -px(-1,-1) - 2.0d0*px(-1,0) - px(-1,1) &
             + px(1,-1) + 2.0d0*px(1,0) + px(1,1)

        if (use_l2 == 1) then
          dst((c-1)*rows + r) = sqrt(gx*gx + gy*gy)
        else
          dst((c-1)*rows + r) = abs(gx) + abs(gy)
        end if
      end do
    end do

  end subroutine sobel_magnitude

  ! SOBEL EDGE DIRECTION

  subroutine sobel_direction(src, rows, cols, dst) &
      bind(C, name="sobel_direction")

    integer(c_int), intent(in), value :: rows, cols
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: r, c, ri, ci
    real(c_double) :: gx, gy
    real(c_double) :: px(-1:1, -1:1)

    do c = 1, cols
      do r = 1, rows
        do ci = -1, 1
          do ri = -1, 1
            px(ri, ci) = src( &
              (max(1,min(c+ci,cols))-1)*rows + max(1,min(r+ri,rows)) )
          end do
        end do

        gx = -px(-1,-1) + px(-1,1) &
             - 2.0d0*px(0,-1) + 2.0d0*px(0,1) &
             - px(1,-1) + px(1,1)

        gy = -px(-1,-1) - 2.0d0*px(-1,0) - px(-1,1) &
             + px(1,-1) + 2.0d0*px(1,0) + px(1,1)

        dst((c-1)*rows + r) = atan2(gy, gx)
      end do
    end do

  end subroutine sobel_direction


  ! LAPLACIAN FILTER

  subroutine laplacian(src, rows, cols, connectivity, dst) &
      bind(C, name="laplacian")

    integer(c_int), intent(in), value :: rows, cols, connectivity
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: r, c
    real(c_double) :: centre, top, bot, lft, rgt
    real(c_double) :: tl, tr, bl, br
    real(c_double) :: diag_w, centre_w

    if (connectivity == 8) then
      diag_w   = 1.0d0
      centre_w = -8.0d0
    else
      diag_w   = 0.0d0
      centre_w = -4.0d0
    end if

    do c = 1, cols
      do r = 1, rows
        centre = src((c-1)*rows + r)
        top    = src((c-1)*rows + max(r-1, 1))
        bot    = src((c-1)*rows + min(r+1, rows))
        lft    = src((max(c-1,1)-1)*rows + r)
        rgt    = src((min(c+1,cols)-1)*rows + r)

        tl     = src((max(c-1,1)-1)*rows + max(r-1,1))
        tr     = src((max(c-1,1)-1)*rows + min(r+1,rows))
        bl     = src((min(c+1,cols)-1)*rows + max(r-1,1))
        br     = src((min(c+1,cols)-1)*rows + min(r+1,rows))

        dst((c-1)*rows + r) = centre_w * centre &
            + top + bot + lft + rgt  &
            + diag_w * (tl + tr + bl + br)
      end do
    end do

  end subroutine laplacian

  ! UNSHARP MASK

  subroutine unsharp_mask(src, rows, cols, sigma, kradius, amount, threshold, dst) &
      bind(C, name="unsharp_mask")

    integer(c_int), intent(in), value :: rows, cols, kradius
    real(c_double), intent(in), value :: sigma, amount, threshold
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: i, n
    real(c_double), allocatable :: blurred(:)
    real(c_double) :: residual

    n = rows * cols
    allocate(blurred(n))

    call gaussian_blur(src, rows, cols, sigma, kradius, blurred)

    do i = 1, n
      residual = src(i) - blurred(i)
      if (abs(residual) >= threshold) then
        dst(i) = src(i) + amount * residual
      else
        dst(i) = src(i)
      end if
    end do

    deallocate(blurred)

  end subroutine unsharp_mask


  ! MEDIAN FILTER

  subroutine median_filter(src, rows, cols, radius, dst) &
      bind(C, name="median_filter")

    integer(c_int), intent(in), value :: rows, cols, radius
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: r, c, ri, ci, nr, nc, wsize, mid, i, j
    real(c_double), allocatable :: win(:)
    real(c_double) :: tmp

    wsize = (2*radius+1) * (2*radius+1)
    allocate(win(wsize))

    do c = 1, cols
      do r = 1, rows
        k_cnt: block
          integer :: cnt
          cnt = 0
          do ci = -radius, radius
            nc = max(1, min(c + ci, cols))
            do ri = -radius, radius
              nr = max(1, min(r + ri, rows))
              cnt = cnt + 1
              win(cnt) = src((nc-1)*rows + nr)
            end do
          end do

          do i = 2, cnt
            tmp = win(i)
            j   = i - 1
            do while (j >= 1 .and. win(j) > tmp)
              win(j+1) = win(j)
              j = j - 1
            end do
            win(j+1) = tmp
          end do

          mid = cnt / 2 + 1
          dst((c-1)*rows + r) = win(mid)
        end block k_cnt
      end do
    end do

    deallocate(win)

  end subroutine median_filter


  ! BILATERAL FILTER

  subroutine bilateral_filter(src, rows, cols, sigma_s, sigma_r, kradius, dst) &
      bind(C, name="bilateral_filter")

    integer(c_int), intent(in), value :: rows, cols, kradius
    real(c_double), intent(in), value :: sigma_s, sigma_r
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: r, c, ri, ci, nr, nc
    real(c_double) :: acc, wsum, ws, wr, w, centre, d_spatial, d_range
    real(c_double) :: inv_s2, inv_r2

    inv_s2 = 1.0d0 / (2.0d0 * sigma_s * sigma_s)
    inv_r2 = 1.0d0 / (2.0d0 * sigma_r * sigma_r)

    do c = 1, cols
      do r = 1, rows
        centre = src((c-1)*rows + r)
        acc    = 0.0d0
        wsum   = 0.0d0

        do ci = -kradius, kradius
          nc = max(1, min(c + ci, cols))
          d_spatial = real(ci, c_double)
          ws = d_spatial * d_spatial
          do ri = -kradius, kradius
            nr = max(1, min(r + ri, rows))
            d_spatial = real(ri, c_double)
            d_range   = src((nc-1)*rows + nr) - centre
            w = exp(-(ws + d_spatial*d_spatial) * inv_s2 &
                    - d_range * d_range * inv_r2)
            acc  = acc  + w * src((nc-1)*rows + nr)
            wsum = wsum + w
          end do
        end do

        dst((c-1)*rows + r) = acc / wsum
      end do
    end do

  end subroutine bilateral_filter


  ! MORPHOLOGICAL EROSION  (rectangular structuring element, radius r)

  subroutine morpho_erode(src, rows, cols, radius, dst) &
      bind(C, name="morpho_erode")

    integer(c_int), intent(in), value :: rows, cols, radius
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: r, c, ri, ci, nr, nc
    real(c_double) :: mn, v

    do c = 1, cols
      do r = 1, rows
        mn = huge(mn)
        do ci = -radius, radius
          nc = max(1, min(c + ci, cols))
          do ri = -radius, radius
            nr = max(1, min(r + ri, rows))
            v  = src((nc-1)*rows + nr)
            if (v < mn) mn = v
          end do
        end do
        dst((c-1)*rows + r) = mn
      end do
    end do

  end subroutine morpho_erode


  ! MORPHOLOGICAL DILATION  (rectangular structuring element, radius r)

  subroutine morpho_dilate(src, rows, cols, radius, dst) &
      bind(C, name="morpho_dilate")

    integer(c_int), intent(in), value :: rows, cols, radius
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: r, c, ri, ci, nr, nc
    real(c_double) :: mx, v

    do c = 1, cols
      do r = 1, rows
        mx = -huge(mx)
        do ci = -radius, radius
          nc = max(1, min(c + ci, cols))
          do ri = -radius, radius
            nr = max(1, min(r + ri, rows))
            v  = src((nc-1)*rows + nr)
            if (v > mx) mx = v
          end do
        end do
        dst((c-1)*rows + r) = mx
      end do
    end do

  end subroutine morpho_dilate

  ! MORPHOLOGICAL WHITE TOP-HAT

  subroutine top_hat_white(src, rows, cols, radius, dst) &
      bind(C, name="top_hat_white")

    integer(c_int), intent(in), value :: rows, cols, radius
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: n, i
    real(c_double), allocatable :: eroded(:), opened(:)

    n = rows * cols
    allocate(eroded(n), opened(n))

    call morpho_erode  (src,    rows, cols, radius, eroded)
    call morpho_dilate (eroded, rows, cols, radius, opened)

    do i = 1, n
      dst(i) = src(i) - opened(i)
    end do

    deallocate(eroded, opened)

  end subroutine top_hat_white


  ! MORPHOLOGICAL BLACK TOP-HAT

  subroutine top_hat_black(src, rows, cols, radius, dst) &
      bind(C, name="top_hat_black")

    integer(c_int), intent(in), value :: rows, cols, radius
    real(c_double), intent(in)        :: src(rows * cols)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: n, i
    real(c_double), allocatable :: dilated(:), closed(:)

    n = rows * cols
    allocate(dilated(n), closed(n))

    call morpho_dilate (src,     rows, cols, radius, dilated)
    call morpho_erode  (dilated, rows, cols, radius, closed)

    do i = 1, n
      dst(i) = closed(i) - src(i)
    end do

    deallocate(dilated, closed)

  end subroutine top_hat_black

end module filters_mod