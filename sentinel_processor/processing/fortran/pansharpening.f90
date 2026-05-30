module pansharpening_mod
  use iso_c_binding
  implicit none

contains

  subroutine bilinear_upsample(src, ms_rows, ms_cols, dst, rows, cols)
    integer,          intent(in)  :: ms_rows, ms_cols, rows, cols
    real(c_double),   intent(in)  :: src(ms_rows * ms_cols)
    real(c_double),   intent(out) :: dst(rows * cols)

    integer        :: r, c, r0, c0, r1, c1
    real(c_double) :: fr, fc, v00, v01, v10, v11
    real(c_double) :: scale_r, scale_c

    scale_r = real(ms_rows - 1, c_double) / real(max(rows - 1, 1), c_double)
    scale_c = real(ms_cols - 1, c_double) / real(max(cols - 1, 1), c_double)

    do c = 1, cols
      do r = 1, rows
        fr = (r - 1) * scale_r
        fc = (c - 1) * scale_c
        r0 = int(fr) + 1
        c0 = int(fc) + 1
        r1 = min(r0 + 1, ms_rows)
        c1 = min(c0 + 1, ms_cols)
        fr = fr - int(fr)
        fc = fc - int(fc)
        v00 = src((c0 - 1) * ms_rows + r0)
        v10 = src((c0 - 1) * ms_rows + r1)
        v01 = src((c1 - 1) * ms_rows + r0)
        v11 = src((c1 - 1) * ms_rows + r1)
        dst((c - 1) * rows + r) = &
          (1.0d0 - fr) * (1.0d0 - fc) * v00 + &
          fr            * (1.0d0 - fc) * v10 + &
          (1.0d0 - fr) * fc            * v01 + &
          fr            * fc            * v11
      end do
    end do
  end subroutine bilinear_upsample


  pure function arr_mean(a, n) result(m)
    integer,        intent(in) :: n
    real(c_double), intent(in) :: a(n)
    real(c_double)             :: m
    integer :: i
    m = 0.0d0
    do i = 1, n; m = m + a(i); end do
    m = m / real(n, c_double)
  end function arr_mean


  pure function arr_std(a, n, mean) result(s)
    integer,        intent(in) :: n
    real(c_double), intent(in) :: a(n), mean
    real(c_double)             :: s
    integer        :: i
    real(c_double) :: acc
    acc = 0.0d0
    do i = 1, n; acc = acc + (a(i) - mean)**2; end do
    s = sqrt(acc / real(n, c_double) + 1.0d-12)
  end function arr_std


  pure function arr_dot(a, b, n) result(d)
    integer,        intent(in) :: n
    real(c_double), intent(in) :: a(n), b(n)
    real(c_double)             :: d
    integer :: i
    d = 0.0d0
    do i = 1, n; d = d + a(i) * b(i); end do
  end function arr_dot


  subroutine gram_schmidt_sharpen( &
      pan_data,  rows,    cols,    &
      ms_data,   ms_rows, ms_cols, &
      n_bands,                     &
      out_data)                    &
    bind(C, name="gram_schmidt_sharpen")

    integer(c_int),  intent(in), value :: rows, cols, ms_rows, ms_cols, n_bands
    real(c_double),  intent(in)        :: pan_data(rows * cols)
    real(c_double),  intent(in)        :: ms_data(ms_rows * ms_cols * n_bands)
    real(c_double),  intent(out)       :: out_data(rows * cols * n_bands)

    integer        :: npx, b, i
    real(c_double), allocatable :: up(:,:)
    real(c_double), allocatable :: gs(:,:)
    real(c_double), allocatable :: proj(:)
    real(c_double), allocatable :: pan_m(:)
    real(c_double) :: mu_syn, sig_syn, mu_pan, sig_pan
    real(c_double) :: alpha, beta_coeff

    npx = rows * cols
    allocate(up(npx, n_bands), gs(npx, n_bands), proj(npx), pan_m(npx))

    do b = 1, n_bands
      call bilinear_upsample( &
        ms_data((b-1)*ms_rows*ms_cols + 1 : b*ms_rows*ms_cols), &
        ms_rows, ms_cols, up(:, b), rows, cols)
    end do

    gs(:, 1) = 0.0d0
    do b = 1, n_bands
      do i = 1, npx
        gs(i, 1) = gs(i, 1) + up(i, b)
      end do
    end do
    gs(:, 1) = gs(:, 1) / real(n_bands, c_double)

    do b = 2, n_bands
      gs(:, b) = up(:, b)
      do i = 1, b - 1
        alpha = arr_dot(gs(:, b), gs(:, i), npx) / &
                (arr_dot(gs(:, i), gs(:, i), npx) + 1.0d-12)
        gs(:, b) = gs(:, b) - alpha * gs(:, i)
      end do
    end do

    mu_syn  = arr_mean(gs(:, 1), npx)
    sig_syn = arr_std (gs(:, 1), npx, mu_syn)
    mu_pan  = arr_mean(pan_data, npx)
    sig_pan = arr_std (pan_data, npx, mu_pan)
    do i = 1, npx
      pan_m(i) = mu_syn + (pan_data(i) - mu_pan) * (sig_syn / (sig_pan + 1.0d-12))
    end do

    do b = 1, n_bands
      beta_coeff = arr_dot(up(:, b), gs(:, 1), npx) / &
                   (arr_dot(gs(:, 1), gs(:, 1), npx) + 1.0d-12)
      do i = 1, npx
        out_data((b-1)*npx + i) = up(i, b) + beta_coeff * (pan_m(i) - gs(i, 1))
      end do
    end do

    deallocate(up, gs, proj, pan_m)
  end subroutine gram_schmidt_sharpen


  subroutine ihs_sharpen( &
      pan_data,  rows,    cols,    &
      ms_data,   ms_rows, ms_cols, &
      out_data)                    &
    bind(C, name="ihs_sharpen")

    integer(c_int),  intent(in), value :: rows, cols, ms_rows, ms_cols
    real(c_double),  intent(in)        :: pan_data(rows * cols)
    real(c_double),  intent(in)        :: ms_data(ms_rows * ms_cols * 3)
    real(c_double),  intent(out)       :: out_data(rows * cols * 3)

    integer        :: npx, i
    real(c_double), allocatable :: r(:), g(:), b(:)
    real(c_double), allocatable :: intens(:), v1(:), v2(:)
    real(c_double), allocatable :: pan_m(:)
    real(c_double) :: mu_i, sig_i, mu_p, sig_p
    real(c_double) :: sv6, sv2, delta

    real(c_double), parameter :: INV3  = 1.0d0 / 3.0d0
    real(c_double), parameter :: SQ6   = 2.449489742783178d0
    real(c_double), parameter :: SQ2   = 1.4142135623730951d0

    npx = rows * cols
    allocate(r(npx), g(npx), b(npx), intens(npx), v1(npx), v2(npx), pan_m(npx))

    call bilinear_upsample(ms_data(1              : ms_rows*ms_cols  ), ms_rows, ms_cols, r, rows, cols)
    call bilinear_upsample(ms_data(ms_rows*ms_cols+1 : 2*ms_rows*ms_cols), ms_rows, ms_cols, g, rows, cols)
    call bilinear_upsample(ms_data(2*ms_rows*ms_cols+1 : 3*ms_rows*ms_cols), ms_rows, ms_cols, b, rows, cols)

    sv6 = 1.0d0 / SQ6
    sv2 = 1.0d0 / SQ2
    do i = 1, npx
      intens(i) = INV3 * (r(i) + g(i) + b(i))
      v1(i)     = sv6  * (-r(i) - g(i) + 2.0d0 * b(i))
      v2(i)     = sv2  * ( r(i) - g(i))
    end do

    mu_i  = arr_mean(intens, npx)
    sig_i = arr_std (intens, npx, mu_i)
    mu_p  = arr_mean(pan_data, npx)
    sig_p = arr_std (pan_data, npx, mu_p)
    do i = 1, npx
      pan_m(i) = mu_i + (pan_data(i) - mu_p) * (sig_i / (sig_p + 1.0d-12))
    end do

    do i = 1, npx
      delta = pan_m(i) - intens(i)
      out_data(          i) = r(i) + delta
      out_data(  npx  +  i) = g(i) + delta
      out_data(2*npx  +  i) = b(i) + delta
    end do

    deallocate(r, g, b, intens, v1, v2, pan_m)
  end subroutine ihs_sharpen


  subroutine wavelet_sharpen( &
      pan_data,  rows,    cols,    &
      ms_data,   ms_rows, ms_cols, &
      n_bands,                     &
      out_data)                    &
    bind(C, name="wavelet_sharpen")

    integer(c_int),  intent(in), value :: rows, cols, ms_rows, ms_cols, n_bands
    real(c_double),  intent(in)        :: pan_data(rows * cols)
    real(c_double),  intent(in)        :: ms_data(ms_rows * ms_cols * n_bands)
    real(c_double),  intent(out)       :: out_data(rows * cols * n_bands)

    integer :: npx, b, hr, hc
    real(c_double), allocatable :: up(:)
    real(c_double), allocatable :: pan2d(:,:)
    real(c_double), allocatable :: ms2d(:,:)
    real(c_double), allocatable :: pan_ll(:,:), pan_lh(:,:), pan_hl(:,:), pan_hh(:,:)
    real(c_double), allocatable :: ms_ll(:,:),  ms_lh(:,:),  ms_hl(:,:),  ms_hh(:,:)
    real(c_double), allocatable :: rec(:,:)

    npx = rows * cols
    hr  = rows / 2
    hc  = cols / 2

    allocate(up(npx))
    allocate(pan2d(rows, cols), ms2d(rows, cols))
    allocate(pan_ll(hr, hc), pan_lh(hr, hc), pan_hl(hr, hc), pan_hh(hr, hc))
    allocate(ms_ll(hr, hc),  ms_lh(hr, hc),  ms_hl(hr, hc),  ms_hh(hr, hc))
    allocate(rec(rows, cols))

    pan2d = reshape(pan_data, [rows, cols])
    call haar_fwd_2d(pan2d, rows, cols, hr, hc, pan_ll, pan_lh, pan_hl, pan_hh)

    do b = 1, n_bands
      call bilinear_upsample( &
        ms_data((b-1)*ms_rows*ms_cols + 1 : b*ms_rows*ms_cols), &
        ms_rows, ms_cols, up, rows, cols)
      ms2d = reshape(up, [rows, cols])
      call haar_fwd_2d(ms2d, rows, cols, hr, hc, ms_ll, ms_lh, ms_hl, ms_hh)

      ms_ll = pan_ll

      call haar_inv_2d(ms_ll, ms_lh, ms_hl, ms_hh, hr, hc, rows, cols, rec)
      out_data((b-1)*npx + 1 : b*npx) = reshape(rec, [npx])
    end do

    deallocate(up, pan2d, ms2d, pan_ll, pan_lh, pan_hl, pan_hh, &
               ms_ll, ms_lh, ms_hl, ms_hh, rec)
  end subroutine wavelet_sharpen


  subroutine haar_fwd_2d(img, rows, cols, hr, hc, ll, lh, hl, hh)
    integer,        intent(in)  :: rows, cols, hr, hc
    real(c_double), intent(in)  :: img(rows, cols)
    real(c_double), intent(out) :: ll(hr, hc), lh(hr, hc), hl(hr, hc), hh(hr, hc)

    real(c_double), allocatable :: tmp(:, :)
    integer :: r, c, ri, ci
    real(c_double) :: a, b_val

    allocate(tmp(hr, cols))

    do c = 1, cols
      do r = 1, hr
        ri = 2 * r - 1
        a     = img(ri,   c)
        b_val = img(ri+1, c)
        tmp(r, c) = (a + b_val) * 0.5d0   ! low
      end do
    end do

    do r = 1, hr
      do c = 1, hc
        ci = 2 * c - 1
        ll(r, c) = (tmp(r, ci) + tmp(r, ci+1)) * 0.5d0
        lh(r, c) = (tmp(r, ci) - tmp(r, ci+1)) * 0.5d0
      end do
    end do

    do c = 1, cols
      do r = 1, hr
        ri = 2 * r - 1
        tmp(r, c) = (img(ri, c) - img(ri+1, c)) * 0.5d0  ! high
      end do
    end do
    do r = 1, hr
      do c = 1, hc
        ci = 2 * c - 1
        hl(r, c) = (tmp(r, ci) + tmp(r, ci+1)) * 0.5d0
        hh(r, c) = (tmp(r, ci) - tmp(r, ci+1)) * 0.5d0
      end do
    end do

    deallocate(tmp)
  end subroutine haar_fwd_2d


  subroutine haar_inv_2d(ll, lh, hl, hh, hr, hc, rows, cols, out)
    integer,        intent(in)  :: hr, hc, rows, cols
    real(c_double), intent(in)  :: ll(hr, hc), lh(hr, hc), hl(hr, hc), hh(hr, hc)
    real(c_double), intent(out) :: out(rows, cols)

    real(c_double), allocatable :: row_l(:,:), row_h(:,:)
    integer :: r, c

    allocate(row_l(hr, cols), row_h(hr, cols))

    do r = 1, hr
      do c = 1, hc
        row_l(r, 2*c-1) = ll(r, c) + lh(r, c)
        row_l(r, 2*c  ) = ll(r, c) - lh(r, c)
        row_h(r, 2*c-1) = hl(r, c) + hh(r, c)
        row_h(r, 2*c  ) = hl(r, c) - hh(r, c)
      end do
    end do

    do c = 1, cols
      do r = 1, hr
        out(2*r-1, c) = row_l(r, c) + row_h(r, c)
        out(2*r,   c) = row_l(r, c) - row_h(r, c)
      end do
    end do

    deallocate(row_l, row_h)
  end subroutine haar_inv_2d

end module pansharpening_mod