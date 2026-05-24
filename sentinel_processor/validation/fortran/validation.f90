module validation_mod
  use iso_c_binding
  implicit none

contains

  subroutine validate_scl(scl_values, n, max_cloud_thr, &
      confidence_score, cloud_ratio, snow_ratio, water_excluded, &
      issues_buf, max_issue_len, n_issues) &
    bind(C, name="validate_scl")

    integer(c_int), intent(in),  value :: n, max_issue_len
    integer(c_int), intent(in) :: scl_values(n)
    real(c_double), intent(in),  value :: max_cloud_thr
    real(c_double), intent(out) :: confidence_score, cloud_ratio, snow_ratio
    integer(c_int), intent(out) :: water_excluded, n_issues
    character(c_char), intent(out) :: issues_buf(max_issue_len * 4)

    integer :: i, filtered_total, bad_pixels, snow_pixels
    real(c_double) :: total_d

    filtered_total = 0; bad_pixels = 0; snow_pixels = 0

    do i = 1, n
      if (scl_values(i) == 0 .or. scl_values(i) == 6) cycle
      filtered_total = filtered_total + 1
      if (scl_values(i) == 3 .or. scl_values(i) == 8 .or. &
          scl_values(i) == 9 .or. scl_values(i) == 10) bad_pixels = bad_pixels + 1
      if (scl_values(i) == 11) snow_pixels = snow_pixels + 1
    end do

    water_excluded = 1
    total_d = real(filtered_total, c_double)

    if (filtered_total == 0) then
      cloud_ratio = 0.0d0; snow_ratio = 0.0d0
    else
      cloud_ratio = real(bad_pixels,  c_double) / total_d
      snow_ratio  = real(snow_pixels, c_double) / total_d
    end if

    if (filtered_total == 0) then
      confidence_score = 0.0d0
    else if (snow_ratio > 0.5d0) then
      confidence_score = 0.0d0
    else if (cloud_ratio > max_cloud_thr) then
      confidence_score = 0.0d0
    else if (cloud_ratio < 0.1d0) then
      confidence_score = 1.0d0
    else if (cloud_ratio < 0.3d0) then
      confidence_score = 0.75d0
    else if (cloud_ratio < 0.4d0) then
      confidence_score = 0.5d0
    else
      confidence_score = 0.0d0
    end if

    n_issues = 0
    call zero_buf(issues_buf, max_issue_len * 4)
    if (filtered_total == 0) &
      call write_issue(issues_buf, n_issues, max_issue_len, "No valid pixels after filtering")
    if (cloud_ratio > 0.3d0) &
      call write_issue(issues_buf, n_issues, max_issue_len, "High cloud cover")
    if (snow_ratio > 0.5d0) &
      call write_issue(issues_buf, n_issues, max_issue_len, "Excessive snow/ice")

  end subroutine validate_scl


  subroutine check_radiometry(pixels, n, result) &
    bind(C, name="check_radiometry")

    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: pixels(n)
    integer(c_int), intent(out) :: result

    real(c_double), parameter :: LIMIT = 15000.0d0
    integer :: saturated, i

    saturated = 0
    do i = 1, n
      if (pixels(i) > LIMIT) saturated = saturated + 1
    end do

    result = merge(1, 0, real(saturated, c_double) / real(n, c_double) < 0.01d0)

  end subroutine check_radiometry

  subroutine zero_buf(buf, total_len)
    character(c_char), intent(out) :: buf(total_len)
    integer, intent(in) :: total_len
    integer :: i
    do i = 1, total_len; buf(i) = c_null_char; end do
  end subroutine zero_buf

  subroutine write_issue(buf, n_issues, max_len, msg)
    character(c_char), intent(inout) :: buf(max_len * 4)
    integer, intent(inout) :: n_issues
    integer, intent(in) :: max_len
    character(*), intent(in) :: msg
    integer :: offset, j, msg_len
    if (n_issues >= 4) return
    offset = n_issues * max_len + 1
    msg_len = min(len_trim(msg), max_len - 1)
    do j = 1, msg_len; buf(offset + j - 1) = msg(j:j); end do
    buf(offset + msg_len) = c_null_char
    n_issues = n_issues + 1
  end subroutine write_issue

  function ieee_nan() result(nan_val)
    real(c_double) :: nan_val
    nan_val = 0.0d0; nan_val = nan_val / nan_val
  end function ieee_nan

end module validation_mod
