module indices_mod
  use iso_c_binding
  implicit none

contains

  ! -----------------------------------------------------------------------
  ! Normalised-difference helper:  (a - b) / (a + b), NaN-safe
  ! -----------------------------------------------------------------------
  pure function norm_diff(a, b) result(r)
    real(c_double), intent(in) :: a, b
    real(c_double) :: r
    real(c_double) :: denom
    denom = a + b
    if (abs(denom) < 1.0d-9) then
      r = 0.0d0
    else
      r = (a - b) / denom
    end if
  end function norm_diff


  ! -----------------------------------------------------------------------
  ! NDVI  – Normalised Difference Vegetation Index
  !         (NIR - Red) / (NIR + Red)   →  B08, B04
  ! -----------------------------------------------------------------------
  subroutine compute_ndvi(nir, red, n, out) &
      bind(C, name="compute_ndvi")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: nir(n), red(n)
    real(c_double), intent(out) :: out(n)
    integer :: i
    do i = 1, n
      out(i) = norm_diff(nir(i), red(i))
    end do
  end subroutine compute_ndvi


  ! -----------------------------------------------------------------------
  ! EVI   – Enhanced Vegetation Index
  !         2.5*(NIR-Red) / (NIR + 6*Red - 7.5*Blue + 1)  →  B08, B04, B02
  ! -----------------------------------------------------------------------
  subroutine compute_evi(nir, red, blue, n, out) &
      bind(C, name="compute_evi")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: nir(n), red(n), blue(n)
    real(c_double), intent(out) :: out(n)
    real(c_double) :: denom
    integer :: i
    do i = 1, n
      denom = nir(i) + 6.0d0*red(i) - 7.5d0*blue(i) + 1.0d0
      if (abs(denom) < 1.0d-9) then
        out(i) = 0.0d0
      else
        out(i) = 2.5d0 * (nir(i) - red(i)) / denom
      end if
    end do
  end subroutine compute_evi


  ! -----------------------------------------------------------------------
  ! SAVI  – Soil Adjusted Vegetation Index
  !         (NIR-Red)*(1+L) / (NIR+Red+L)   →  B08, B04
  ! -----------------------------------------------------------------------
  subroutine compute_savi(nir, red, n, out) &
      bind(C, name="compute_savi")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: nir(n), red(n)
    real(c_double), intent(out) :: out(n)
    real(c_double), parameter :: L = 0.5d0
    real(c_double) :: denom
    integer :: i
    do i = 1, n
      denom = nir(i) + red(i) + L
      if (abs(denom) < 1.0d-9) then
        out(i) = 0.0d0
      else
        out(i) = (nir(i) - red(i)) * (1.0d0 + L) / denom
      end if
    end do
  end subroutine compute_savi


  ! -----------------------------------------------------------------------
  ! NDWI  – Normalised Difference Water Index
  !         (Green - NIR) / (Green + NIR)   →  B03, B08
  ! -----------------------------------------------------------------------
  subroutine compute_ndwi(green, nir, n, out) &
      bind(C, name="compute_ndwi")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: green(n), nir(n)
    real(c_double), intent(out) :: out(n)
    integer :: i
    do i = 1, n
      out(i) = norm_diff(green(i), nir(i))
    end do
  end subroutine compute_ndwi


  ! -----------------------------------------------------------------------
  ! MNDWI – Modified NDWI
  !         (Green - SWIR1) / (Green + SWIR1)  →  B03, B11
  ! -----------------------------------------------------------------------
  subroutine compute_mndwi(green, swir1, n, out) &
      bind(C, name="compute_mndwi")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: green(n), swir1(n)
    real(c_double), intent(out) :: out(n)
    integer :: i
    do i = 1, n
      out(i) = norm_diff(green(i), swir1(i))
    end do
  end subroutine compute_mndwi


  ! -----------------------------------------------------------------------
  ! NDBI  – Normalised Difference Built-up Index
  !         (SWIR1 - NIR) / (SWIR1 + NIR)   →  B11, B08
  ! -----------------------------------------------------------------------
  subroutine compute_ndbi(swir1, nir, n, out) &
      bind(C, name="compute_ndbi")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: swir1(n), nir(n)
    real(c_double), intent(out) :: out(n)
    integer :: i
    do i = 1, n
      out(i) = norm_diff(swir1(i), nir(i))
    end do
  end subroutine compute_ndbi


  ! -----------------------------------------------------------------------
  ! NBR   – Normalised Burn Ratio
  !         (NIR - SWIR2) / (NIR + SWIR2)   →  B08, B12
  ! -----------------------------------------------------------------------
  subroutine compute_nbr(nir, swir2, n, out) &
      bind(C, name="compute_nbr")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: nir(n), swir2(n)
    real(c_double), intent(out) :: out(n)
    integer :: i
    do i = 1, n
      out(i) = norm_diff(nir(i), swir2(i))
    end do
  end subroutine compute_nbr


  ! -----------------------------------------------------------------------
  ! NDSI  – Normalised Difference Snow Index
  !         (Green - SWIR1) / (Green + SWIR1)  →  B03, B11
  !   (same formula as MNDWI; semantic alias kept for clarity)
  ! -----------------------------------------------------------------------
  subroutine compute_ndsi(green, swir1, n, out) &
      bind(C, name="compute_ndsi")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: green(n), swir1(n)
    real(c_double), intent(out) :: out(n)
    integer :: i
    do i = 1, n
      out(i) = norm_diff(green(i), swir1(i))
    end do
  end subroutine compute_ndsi


  ! -----------------------------------------------------------------------
  ! CIG   – Chlorophyll Index Green
  !         (NIR / Green) - 1              →  B08, B03
  ! -----------------------------------------------------------------------
  subroutine compute_cig(nir, green, n, out) &
      bind(C, name="compute_cig")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: nir(n), green(n)
    real(c_double), intent(out) :: out(n)
    integer :: i
    do i = 1, n
      if (abs(green(i)) < 1.0d-9) then
        out(i) = 0.0d0
      else
        out(i) = (nir(i) / green(i)) - 1.0d0
      end if
    end do
  end subroutine compute_cig


  ! -----------------------------------------------------------------------
  ! ARVI  – Atmospherically Resistant Vegetation Index
  !         (NIR - (2*Red - Blue)) / (NIR + (2*Red - Blue))  →  B08, B04, B02
  ! -----------------------------------------------------------------------
  subroutine compute_arvi(nir, red, blue, n, out) &
      bind(C, name="compute_arvi")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in)  :: nir(n), red(n), blue(n)
    real(c_double), intent(out) :: out(n)
    real(c_double) :: rb
    integer :: i
    do i = 1, n
      rb = 2.0d0*red(i) - blue(i)
      out(i) = norm_diff(nir(i), rb)
    end do
  end subroutine compute_arvi

end module indices_mod