import numpy as np
import pytest

from sentinel_processor.validation._fortran_bridge import (
    call_check_dimensions,
    call_check_radiometry,
    call_validate_scl,
    validate_file,
)

MIN_SIDE = 32
MAX_SIDE = 10_980
MAX_ASPECT = 4.0
RADIO_LIMIT = 15_000.0
RADIO_MAX_FRAC = 0.01


class TestCheckDimensions:

    def test_valid_square(self):
        r = call_check_dimensions(512, 512)
        assert r["passed"] is True
        assert r["issues"] == []

    def test_valid_rectangle(self):
        r = call_check_dimensions(1000, 500)
        assert r["passed"] is True

    def test_exact_minimum(self):
        """MIN_SIDE × MIN_SIDE must pass."""
        r = call_check_dimensions(MIN_SIDE, MIN_SIDE)
        assert r["passed"] is True

    def test_exact_maximum(self):
        """MAX_SIDE × MAX_SIDE must pass."""
        r = call_check_dimensions(MAX_SIDE, MAX_SIDE)
        assert r["passed"] is True

    def test_too_small_rows(self):
        r = call_check_dimensions(MIN_SIDE - 1, 512)
        assert r["passed"] is False
        assert any("small" in issue.lower() for issue in r["issues"])

    def test_too_small_cols(self):
        r = call_check_dimensions(512, MIN_SIDE - 1)
        assert r["passed"] is False

    def test_too_small_both(self):
        r = call_check_dimensions(1, 1)
        assert r["passed"] is False

    def test_too_large_rows(self):
        r = call_check_dimensions(MAX_SIDE + 1, 512)
        assert r["passed"] is False
        assert any("large" in issue.lower() for issue in r["issues"])

    def test_too_large_cols(self):
        r = call_check_dimensions(512, MAX_SIDE + 1)
        assert r["passed"] is False

    def test_aspect_ratio_exactly_at_limit(self):
        """rows/cols = 4.0 — exactly at limit, should still pass."""
        r = call_check_dimensions(512, 128)  # 512/128 = 4.0
        assert r["passed"] is True

    def test_aspect_ratio_exceeds_limit(self):
        """rows/cols > 4.0 — degenerate shape."""
        r = call_check_dimensions(513, 128)
        assert r["passed"] is False
        assert any("aspect" in issue.lower() or "degenerate" in issue.lower()
                   for issue in r["issues"])

    def test_issues_is_list(self):
        r = call_check_dimensions(64, 64)
        assert isinstance(r["issues"], list)

    def test_multiple_failures_reported(self):
        """Both too-small AND bad aspect: at least 2 issues (or 1 combined)."""
        r = call_check_dimensions(1, 1)
        assert r["passed"] is False
        assert len(r["issues"]) >= 1


# call_validate_scl

class TestValidateScl:
    """
    SCL class reference (Sentinel-2):
      0  = No Data       → filtered out (excluded)
      1  = Saturated/Def → counted as bad? No — only 3,8,9,10 are bad
      2  = Dark Area
      3  = Cloud Shadow  → bad
      4  = Vegetation
      5  = Bare Soil
      6  = Water         → filtered out (excluded)
      7  = Unclassified
      8  = Cloud med     → bad
      9  = Cloud high    → bad
      10 = Thin Cirrus   → bad
      11 = Snow/Ice      → snow counter
    """

    def test_returns_all_keys(self):
        r = call_validate_scl([4, 4, 4, 4])
        assert set(r) == {
            "confidence_score", "cloud_ratio", "snow_ratio",
            "water_excluded", "issues",
        }

    def test_water_excluded_always_true(self):
        """Fortran always sets water_excluded=1."""
        r = call_validate_scl([4, 5, 2])
        assert r["water_excluded"] is True

    def test_clear_scene_confidence_1(self):
        """All vegetation — cloud_ratio < 0.1 -> confidence = 1.0."""
        scl = [4] * 1000
        r = call_validate_scl(scl)
        assert r["confidence_score"] == pytest.approx(1.0)
        assert r["cloud_ratio"] == pytest.approx(0.0)
        assert r["issues"] == []

    def test_confidence_075_when_cloud_between_10_and_30_pct(self):
        """cloud_ratio = 0.20 -> confidence = 0.75."""
        scl = [8] * 20 + [4] * 80  # 20% cloud
        r = call_validate_scl(scl, max_cloud_threshold=0.30)
        assert r["cloud_ratio"] == pytest.approx(0.20)
        assert r["confidence_score"] == pytest.approx(0.75)

    def test_confidence_05_when_cloud_between_30_and_40_pct(self):
        """cloud_ratio = 0.35, threshold=0.40 -> not > threshold → falls to
        the 0.3 < cloud < 0.4 branch -> confidence = 0.5."""
        scl = [9] * 35 + [4] * 65
        r = call_validate_scl(scl, max_cloud_threshold=0.40)
        assert r["cloud_ratio"] == pytest.approx(0.35)
        assert r["confidence_score"] == pytest.approx(0.5)

    def test_confidence_0_when_cloud_exceeds_threshold(self):
        """cloud_ratio > max_cloud_threshold -> confidence = 0."""
        scl = [10] * 40 + [4] * 60
        r = call_validate_scl(scl, max_cloud_threshold=0.30)
        assert r["confidence_score"] == pytest.approx(0.0)

    def test_confidence_0_when_cloud_above_40_pct(self):
        """cloud_ratio = 0.45 -> confidence = 0 regardless of threshold."""
        scl = [8] * 45 + [4] * 55
        r = call_validate_scl(scl, max_cloud_threshold=0.50)
        assert r["confidence_score"] == pytest.approx(0.0)

    def test_custom_threshold_respected(self):
        """With threshold=0.50, 40% cloud should still get confidence=0.5."""
        scl = [8] * 35 + [4] * 65
        r = call_validate_scl(scl, max_cloud_threshold=0.50)
        assert r["confidence_score"] == pytest.approx(0.5)

    def test_snow_ratio_computed(self):
        scl = [11] * 30 + [4] * 70
        r = call_validate_scl(scl)
        assert r["snow_ratio"] == pytest.approx(0.30)

    def test_excessive_snow_zeroes_confidence(self):
        """snow_ratio > 0.5 -> confidence = 0."""
        scl = [11] * 60 + [4] * 40
        r = call_validate_scl(scl)
        assert r["confidence_score"] == pytest.approx(0.0)
        assert any("snow" in i.lower() for i in r["issues"])

    def test_snow_below_threshold_does_not_zero(self):
        scl = [11] * 40 + [4] * 60
        r = call_validate_scl(scl)
        assert r["confidence_score"] > 0.0

    def test_nodata_pixels_excluded_from_ratio(self):
        """SCL=0 pixels are excluded; 10% of remaining are cloud."""
        scl = [0] * 500 + [8] * 50 + [4] * 450  # 50/500 valid = 10% cloud
        r = call_validate_scl(scl)
        assert r["cloud_ratio"] == pytest.approx(0.10, abs=1e-6)

    def test_water_pixels_excluded_from_ratio(self):
        scl = [6] * 500 + [8] * 50 + [4] * 450
        r = call_validate_scl(scl)
        assert r["cloud_ratio"] == pytest.approx(0.10, abs=1e-6)

    def test_all_nodata_returns_zero_confidence(self):
        scl = [0] * 100
        r = call_validate_scl(scl)
        assert r["confidence_score"] == pytest.approx(0.0)
        assert any("no valid" in i.lower() for i in r["issues"])

    def test_all_water_returns_zero_confidence(self):
        scl = [6] * 100
        r = call_validate_scl(scl)
        assert r["confidence_score"] == pytest.approx(0.0)

    def test_accepts_numpy_array(self):
        scl = np.array([4, 4, 8, 4], dtype=np.int32)
        r = call_validate_scl(scl)
        assert isinstance(r["cloud_ratio"], float)

    def test_all_bad_classes_counted(self):
        """Classes 3, 8, 9, 10 all count as bad pixels."""
        scl = [3, 8, 9, 10, 4, 4, 4, 4, 4, 4]  # 40% bad
        r = call_validate_scl(scl)
        assert r["cloud_ratio"] == pytest.approx(0.40, abs=1e-6)

    def test_high_cloud_issue_message(self):
        """cloud_ratio > 0.3 -> 'High cloud cover' in issues."""
        scl = [8] * 40 + [4] * 60
        r = call_validate_scl(scl, max_cloud_threshold=0.50)
        assert any("cloud" in i.lower() for i in r["issues"])


# call_check_radiometry

class TestCheckRadiometry:

    def test_all_normal_pixels_pass(self):
        pixels = np.full(1000, 5000.0)
        assert call_check_radiometry(pixels) is True

    def test_exactly_at_limit_passes(self):
        """pixel == 15000 is NOT > 15000 → not saturated."""
        pixels = np.full(100, 15_000.0)
        assert call_check_radiometry(pixels) is True

    def test_one_pixel_above_limit_in_large_array_passes(self):
        """1/1000 = 0.1% saturated < 1% threshold → passes."""
        pixels = np.full(1000, 5000.0)
        pixels[0] = 15_001.0
        assert call_check_radiometry(pixels) is True

    def test_exactly_1pct_saturated_passes(self):
        """Exactly 10/1000 = 1.0% but condition is < 0.01 so fails."""
        pixels = np.full(1000, 5000.0)
        pixels[:10] = 15_001.0  # 1.0% — NOT < 0.01 → fails
        assert call_check_radiometry(pixels) is False

    def test_just_below_1pct_passes(self):
        """9/1000 = 0.9% < 1% → passes."""
        pixels = np.full(1000, 5000.0)
        pixels[:9] = 15_001.0
        assert call_check_radiometry(pixels) is True

    def test_high_saturation_fails(self):
        """50% saturated → fails."""
        pixels = np.full(100, 5000.0)
        pixels[:50] = 20_000.0
        assert call_check_radiometry(pixels) is False

    def test_all_saturated_fails(self):
        pixels = np.full(100, 20_000.0)
        assert call_check_radiometry(pixels) is False

    def test_accepts_list_input(self):
        pixels = [1000.0] * 50
        assert call_check_radiometry(pixels) is True

    def test_returns_bool(self):
        assert isinstance(call_check_radiometry([5000.0] * 10), bool)


# validate_file


class TestValidateFile:

    def test_missing_file_returns_error_dict(self):
        r = validate_file("/nonexistent/path/scl.tif")
        assert r["passed"] is False
        assert "error" in r
        assert r["file"] == "/nonexistent/path/scl.tif"

    def test_error_dict_has_file_key(self):
        r = validate_file("does_not_exist.nc")
        assert "file" in r

    def test_min_confidence_respected(self, tmp_path):
        pytest.importorskip("rioxarray")
        pytest.importorskip("rasterio")
        import rasterio
        from rasterio.transform import from_bounds

        tif = tmp_path / "scl_test.tif"
        # 64×64 pixels, all vegetation (SCL=4) -> confidence=1.0
        data = np.full((1, 64, 64), 4, dtype=np.uint8)
        transform = from_bounds(0, 0, 1, 1, 64, 64)
        with rasterio.open(
                str(tif), "w",
                driver="GTiff", height=64, width=64,
                count=1, dtype="uint8",
                crs="EPSG:4326", transform=transform,
        ) as dst:
            dst.write(data)

        # min_confidence=0.99 -> should pass (confidence=1.0)
        r = validate_file(str(tif), min_confidence=0.99)
        assert r["passed"] is True
        assert r["confidence_score"] == pytest.approx(1.0)

        # min_confidence higher than possible -> fails
        r2 = validate_file(str(tif), min_confidence=1.01)
        assert r2["passed"] is False