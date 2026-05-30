module raster_ops_mod
  use iso_c_binding
  implicit none

contains

  subroutine rgb_to_luminance(src, rows, cols, n_bands, dst) &
    bind(C, name="rgb_to_luminance")

    integer(c_int), intent(in), value :: rows, cols, n_bands
    real(c_double), intent(in)        :: src(rows * cols * n_bands)
    real(c_double), intent(out)       :: dst(rows * cols)

    integer        :: npx, i, b
    real(c_double) :: w(n_bands)

    npx = rows * cols

    if (n_bands == 1) then
      dst = src(1:npx)
      return
    end if

    if (n_bands == 3) then
      w(1) = 0.299d0   ! R
      w(2) = 0.587d0   ! G
      w(3) = 0.114d0   ! B
    else
      do b = 1, n_bands
        w(b) = 1.0d0 / real(n_bands, c_double)
      end do
    end if

    dst = 0.0d0
    do b = 1, n_bands
      do i = 1, npx
        dst(i) = dst(i) + w(b) * src((b-1)*npx + i)
      end do
    end do
  end subroutine rgb_to_luminance


  subroutine align_bands(src, src_rows, src_cols, dst, dst_rows, dst_cols) &
    bind(C, name="align_bands")

    integer(c_int), intent(in), value :: src_rows, src_cols, dst_rows, dst_cols
    real(c_double), intent(in)        :: src(src_rows * src_cols)
    real(c_double), intent(out)       :: dst(dst_rows * dst_cols)

    integer        :: r, c, sr, sc
    real(c_double) :: scale_r, scale_c

    scale_r = real(src_rows, c_double) / real(dst_rows, c_double)
    scale_c = real(src_cols, c_double) / real(dst_cols, c_double)

    do c = 1, dst_cols
      sc = min(int((c - 1) * scale_c) + 1, src_cols)
      do r = 1, dst_rows
        sr = min(int((r - 1) * scale_r) + 1, src_rows)
        dst((c-1)*dst_rows + r) = src((sc-1)*src_rows + sr)
      end do
    end do
  end subroutine align_bands


  subroutine clip_box_indices( &
      origin_x, origin_y, pixel_w, pixel_h, &
      rows, cols,                             &
      min_lon, max_lon, min_lat, max_lat,     &
      row_min, row_max, col_min, col_max)     &
    bind(C, name="clip_box_indices")

    real(c_double), intent(in), value :: origin_x, origin_y, pixel_w, pixel_h
    integer(c_int), intent(in), value :: rows, cols
    real(c_double), intent(in), value :: min_lon, max_lon, min_lat, max_lat
    integer(c_int), intent(out)       :: row_min, row_max, col_min, col_max

    integer        :: r0, r1, c0, c1
    real(c_double) :: inv_h, inv_w

    inv_w = 1.0d0 / pixel_w
    inv_h = 1.0d0 / pixel_h

    c0 = int((min_lon - origin_x) * inv_w)
    c1 = int((max_lon - origin_x) * inv_w) + 1

    if (pixel_h < 0.0d0) then
      r0 = int((max_lat - origin_y) * inv_h)
      r1 = int((min_lat - origin_y) * inv_h) + 1
    else
      r0 = int((min_lat - origin_y) * inv_h)
      r1 = int((max_lat - origin_y) * inv_h) + 1
    end if

    col_min = max(c0 + 1, 1)
    col_max = min(c1,     cols)
    row_min = max(r0 + 1, 1)
    row_max = min(r1,     rows)

  end subroutine clip_box_indices


  subroutine band_stats(arr, n, mean_out, std_out, min_out, max_out) &
    bind(C, name="band_stats")

    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)        :: arr(n)
    real(c_double), intent(out)       :: mean_out, std_out, min_out, max_out

    integer        :: i
    real(c_double) :: acc, acc2

    if (n <= 0) then
      mean_out = 0.0d0; std_out = 0.0d0
      min_out  = 0.0d0; max_out = 0.0d0
      return
    end if

    acc  = arr(1)
    acc2 = 0.0d0
    min_out = arr(1)
    max_out = arr(1)

    do i = 2, n
      acc  = acc + arr(i)
      if (arr(i) < min_out) min_out = arr(i)
      if (arr(i) > max_out) max_out = arr(i)
    end do
    mean_out = acc / real(n, c_double)

    do i = 1, n
      acc2 = acc2 + (arr(i) - mean_out)**2
    end do
    std_out = sqrt(acc2 / real(n, c_double) + 1.0d-12)

  end subroutine band_stats


  subroutine normalize_band(arr, n, src_min, src_max, out_min, out_max) &
    bind(C, name="normalize_band")

    integer(c_int), intent(in), value :: n
    real(c_double), intent(in), value :: src_min, src_max, out_min, out_max
    real(c_double), intent(inout)     :: arr(n)

    integer        :: i
    real(c_double) :: scale, range

    range = src_max - src_min
    if (abs(range) < 1.0d-12) then
      arr = out_min
      return
    end if
    scale = (out_max - out_min) / range
    do i = 1, n
      arr(i) = out_min + (arr(i) - src_min) * scale
    end do

  end subroutine normalize_band


  subroutine reproject_nearest( &
      src,    src_rows, src_cols,               &
      src_ox, src_oy,   src_pw, src_ph,         &
      dst,    dst_rows, dst_cols,               &
      dst_ox, dst_oy,   dst_pw, dst_ph)         &
    bind(C, name="reproject_nearest")

    integer(c_int), intent(in), value :: src_rows, src_cols, dst_rows, dst_cols
    real(c_double), intent(in), value :: src_ox, src_oy, src_pw, src_ph
    real(c_double), intent(in), value :: dst_ox, dst_oy, dst_pw, dst_ph
    real(c_double), intent(in)        :: src(src_rows * src_cols)
    real(c_double), intent(out)       :: dst(dst_rows * dst_cols)

    integer        :: dr, dc, sr, sc
    real(c_double) :: lon, lat

    do dc = 1, dst_cols
      lon = dst_ox + (dc - 1) * dst_pw
      sc  = nint((lon - src_ox) / src_pw) + 1
      sc  = max(1, min(sc, src_cols))
      do dr = 1, dst_rows
        lat = dst_oy + (dr - 1) * dst_ph
        sr  = nint((lat - src_oy) / src_ph) + 1
        sr  = max(1, min(sr, src_rows))
        dst((dc-1)*dst_rows + dr) = src((sc-1)*src_rows + sr)
      end do
    end do

  end subroutine reproject_nearest

end module raster_ops_mod