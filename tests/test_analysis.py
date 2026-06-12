from __future__ import annotations
from datetime import datetime, timedelta, timezone
from unittest.mock import patch
import numpy as np
import pytest

try:
    from sentinel_processor.analysis._sentinel_stats_bridge import (
        NODATA,
        anomaly_zscore,
        bfast_breakpoint,
        mann_kendall,
        pearson_map,
        phenology_doy,
        phenology_metrics,
        pixel_iqr,
        pixel_quantiles,
        pixel_regression,
        save_phenology,
        temporal_gap_stats,
        time_window_stats,
        trend_theil_sen,
        valid_obs_count,
    )

    _BRIDGE_AVAILABLE = True
except (ImportError, FileNotFoundError, OSError):
    _BRIDGE_AVAILABLE = False

pytestmark = pytest.mark.fortran
skip_no_lib = pytest.mark.skipif(
    not _BRIDGE_AVAILABLE,
    reason="libsentinel_stats not compiled — run python build_sentinel_stats.py",
)


def _dates(n: int, start: datetime | None = None, step_days: float = 16.0) -> list[datetime]:
    t0 = start or datetime(2023, 1, 1, tzinfo=timezone.utc)
    return [t0 + timedelta(days=step_days * i) for i in range(n)]


def _stack(n: int, rows: int, cols: int, seed: int = 0) -> np.ndarray:
    """Random float64 stack (n, rows, cols) with values in [0, 1]."""
    return np.random.default_rng(seed).random((n, rows, cols))


def _nodata_stack(n: int, rows: int, cols: int) -> np.ndarray:
    """Stack where every pixel is NODATA."""
    return np.full((n, rows, cols), NODATA, dtype=np.float64)


def _linear_pixel(n: int, slope: float, intercept: float,
                  dates: list[datetime]) -> np.ndarray:
    """1×1 stack with a perfect linear trend."""
    d = np.array([(dt - dates[0]).total_seconds() / 86400.0 for dt in dates])
    vals = intercept + slope * d
    return vals.reshape(n, 1, 1)


# NODATA constant

class TestNodata:

    def test_nodata_is_minus_9999(self):
        assert NODATA == -9999.0

    def test_nodata_type_float(self):
        assert isinstance(NODATA, float)


# valid_obs_count


@skip_no_lib
class TestValidObsCount:

    def test_all_valid_count_equals_n(self):
        arr = _stack(10, 8, 8)
        res = valid_obs_count(arr)
        assert res["count"].shape == (8, 8)
        np.testing.assert_array_equal(res["count"], 10.0)

    def test_all_valid_fraction_is_one(self):
        arr = _stack(10, 8, 8)
        res = valid_obs_count(arr)
        np.testing.assert_allclose(res["fraction"], 1.0)

    def test_all_nodata_count_is_zero(self):
        arr = _nodata_stack(10, 4, 4)
        res = valid_obs_count(arr)
        np.testing.assert_array_equal(res["count"], 0.0)

    def test_all_nodata_fraction_is_zero(self):
        arr = _nodata_stack(10, 4, 4)
        res = valid_obs_count(arr)
        np.testing.assert_allclose(res["fraction"], 0.0)

    def test_partial_nodata_counts_correctly(self):
        arr = _stack(10, 4, 4)
        arr[3:, 0, 0] = NODATA  # 7 NODATA → 3 valid
        res = valid_obs_count(arr)
        assert res["count"][0, 0] == pytest.approx(3.0)
        assert res["fraction"][0, 0] == pytest.approx(3.0 / 10.0)

    def test_output_shapes(self):
        arr = _stack(5, 16, 24)
        res = valid_obs_count(arr)
        assert res["count"].shape == (16, 24)
        assert res["fraction"].shape == (16, 24)

    def test_fraction_in_0_1_range(self):
        arr = _stack(8, 6, 6)
        arr[::2, :, :] = NODATA
        res = valid_obs_count(arr)
        assert res["fraction"].min() >= 0.0
        assert res["fraction"].max() <= 1.0

    def test_wrong_ndim_raises(self):
        with pytest.raises(ValueError, match="3-D"):
            valid_obs_count(np.ones((4, 4)))


# temporal_gap_stats


@skip_no_lib
class TestTemporalGapStats:

    def test_uniform_16day_gaps(self):
        n = 8
        dates = _dates(n, step_days=16.0)
        arr = _stack(n, 4, 4)
        res = temporal_gap_stats(arr, dates)
        np.testing.assert_allclose(res["max_gap"], 16.0, atol=0.01)
        np.testing.assert_allclose(res["mean_gap"], 16.0, atol=0.01)

    def test_single_valid_obs_returns_nodata(self):
        n = 6
        dates = _dates(n)
        arr = np.full((n, 1, 1), NODATA, dtype=np.float64)
        arr[2, 0, 0] = 0.5  # only one valid
        res = temporal_gap_stats(arr, dates)
        assert res["max_gap"][0, 0] == NODATA
        assert res["mean_gap"][0, 0] == NODATA

    def test_all_nodata_returns_nodata(self):
        n = 5
        dates = _dates(n)
        arr = _nodata_stack(n, 3, 3)
        res = temporal_gap_stats(arr, dates)
        assert np.all(res["max_gap"] == NODATA)
        assert np.all(res["mean_gap"] == NODATA)

    def test_max_gap_ge_mean_gap(self):
        n = 10
        dates = _dates(n, step_days=5.0)
        arr = _stack(n, 8, 8)
        arr[3, :, :] = NODATA  # create a longer gap on one step
        res = temporal_gap_stats(arr, dates)
        valid = res["max_gap"] != NODATA
        assert np.all(res["max_gap"][valid] >= res["mean_gap"][valid] - 1e-6)

    def test_output_shapes(self):
        n = 6
        dates = _dates(n)
        arr = _stack(n, 12, 10)
        res = temporal_gap_stats(arr, dates)
        assert res["max_gap"].shape == (12, 10)
        assert res["mean_gap"].shape == (12, 10)

    def test_date_length_mismatch_raises(self):
        arr = _stack(5, 4, 4)
        dates = _dates(3)
        with pytest.raises(ValueError):
            temporal_gap_stats(arr, dates)


# pixel_quantiles


@skip_no_lib
class TestPixelQuantiles:

    def test_output_keys(self):
        arr = _stack(10, 4, 4)
        res = pixel_quantiles(arr)
        assert set(res.keys()) == {"p10", "p25", "p50", "p75", "p90"}

    def test_all_nodata_returns_nodata(self):
        arr = _nodata_stack(10, 3, 3)
        res = pixel_quantiles(arr)
        assert np.all(res["p50"] == NODATA)

    def test_constant_series_all_quantiles_equal(self):
        arr = np.full((10, 4, 4), 0.5, dtype=np.float64)
        res = pixel_quantiles(arr)
        for key in ("p10", "p25", "p50", "p75", "p90"):
            np.testing.assert_allclose(res[key], 0.5, atol=1e-6)

    def test_quantile_ordering(self):
        arr = _stack(20, 8, 8, seed=7)
        res = pixel_quantiles(arr)
        valid = res["p10"] != NODATA
        assert np.all(res["p10"][valid] <= res["p25"][valid] + 1e-9)
        assert np.all(res["p25"][valid] <= res["p50"][valid] + 1e-9)
        assert np.all(res["p50"][valid] <= res["p75"][valid] + 1e-9)
        assert np.all(res["p75"][valid] <= res["p90"][valid] + 1e-9)

    def test_p50_approximates_median(self):
        rng = np.random.default_rng(3)
        arr = rng.random((20, 1, 1))
        res = pixel_quantiles(arr)
        expected = float(np.median(arr[:, 0, 0]))
        assert res["p50"][0, 0] == pytest.approx(expected, abs=0.05)

    def test_output_shapes(self):
        arr = _stack(8, 14, 18)
        res = pixel_quantiles(arr)
        for v in res.values():
            assert v.shape == (14, 18)

    def test_wrong_ndim_raises(self):
        with pytest.raises(ValueError, match="3-D"):
            pixel_quantiles(np.ones((10, 10)))


# pixel_iqr


@skip_no_lib
class TestPixelIQR:

    def test_output_keys(self):
        arr = _stack(10, 4, 4)
        res = pixel_iqr(arr)
        assert set(res.keys()) == {"iqr", "outlier"}

    def test_iqr_shape(self):
        arr = _stack(10, 6, 8)
        res = pixel_iqr(arr)
        assert res["iqr"].shape == (6, 8)

    def test_outlier_shape(self):
        arr = _stack(10, 6, 8)
        res = pixel_iqr(arr)
        assert res["outlier"].shape == (10, 6, 8)

    def test_constant_series_iqr_zero(self):
        arr = np.full((10, 3, 3), 0.5, dtype=np.float64)
        res = pixel_iqr(arr)
        np.testing.assert_allclose(res["iqr"], 0.0, atol=1e-6)

    def test_all_nodata_iqr_is_nodata(self):
        arr = _nodata_stack(10, 4, 4)
        res = pixel_iqr(arr)
        assert np.all(res["iqr"] == NODATA)

    def test_extreme_outlier_flagged(self):
        arr = np.full((10, 1, 1), 0.5, dtype=np.float64)
        arr[0, 0, 0] = 1000.0  # extreme outlier
        res = pixel_iqr(arr)
        assert res["outlier"][0, 0, 0] == pytest.approx(1.0)

    def test_inlier_not_flagged(self):
        arr = _stack(20, 1, 1, seed=0)  # values in [0, 1], no wild outliers
        res = pixel_iqr(arr)
        # centre values should not be flagged as outliers
        assert res["outlier"][10, 0, 0] in (0.0, NODATA)

    def test_iqr_non_negative(self):
        arr = _stack(15, 8, 8)
        res = pixel_iqr(arr)
        valid = res["iqr"] != NODATA
        assert np.all(res["iqr"][valid] >= 0.0)


# time_window_stats


@skip_no_lib
class TestTimeWindowStats:

    def test_output_keys(self):
        n = 8
        dates = _dates(n)
        arr = _stack(n, 4, 4)
        res = time_window_stats(arr, dates, window_days=60)
        assert set(res.keys()) == {"mean", "std", "slope"}

    def test_output_shapes(self):
        n = 8
        dates = _dates(n)
        arr = _stack(n, 6, 10)
        res = time_window_stats(arr, dates, window_days=120)
        assert res["mean"].shape == (6, 10)
        assert res["std"].shape == (6, 10)
        assert res["slope"].shape == (6, 10)

    def test_all_nodata_returns_nodata(self):
        n = 6
        dates = _dates(n)
        arr = _nodata_stack(n, 3, 3)
        res = time_window_stats(arr, dates, window_days=60)
        assert np.all(res["mean"] == NODATA)

    def test_constant_series_std_zero(self):
        n = 8
        dates = _dates(n)
        arr = np.full((n, 4, 4), 0.5, dtype=np.float64)
        res = time_window_stats(arr, dates, window_days=200)
        np.testing.assert_allclose(res["std"], 0.0, atol=1e-6)

    def test_constant_series_slope_zero(self):
        n = 8
        dates = _dates(n)
        arr = np.full((n, 4, 4), 0.5, dtype=np.float64)
        res = time_window_stats(arr, dates, window_days=200)
        np.testing.assert_allclose(res["slope"], 0.0, atol=1e-8)

    def test_linear_pixel_slope_recovered(self):
        n = 10
        dates = _dates(n, step_days=16.0)
        slope = 0.003
        arr = _linear_pixel(n, slope=slope, intercept=0.2, dates=dates)
        res = time_window_stats(arr, dates, window_days=200)
        assert res["slope"][0, 0] == pytest.approx(slope, rel=1e-4)

    def test_date_length_mismatch_raises(self):
        arr = _stack(5, 4, 4)
        with pytest.raises(ValueError):
            time_window_stats(arr, _dates(3), window_days=30)

    def test_window_days_zero_raises(self):
        arr = _stack(5, 4, 4)
        with pytest.raises(ValueError, match="window_days"):
            time_window_stats(arr, _dates(5), window_days=0)


# anomaly_zscore


@skip_no_lib
class TestAnomalyZscore:

    def test_output_shape(self):
        n = 8
        dates = _dates(n)
        arr = _stack(n, 6, 8)
        z = anomaly_zscore(arr, dates, window_days=120)
        assert z.shape == (n, 6, 8)

    def test_all_nodata_returns_nodata(self):
        n = 6
        dates = _dates(n)
        arr = _nodata_stack(n, 3, 3)
        z = anomaly_zscore(arr, dates, window_days=60)
        assert np.all(z == NODATA)

    def test_constant_series_is_nodata_or_zero(self):
        """std = 0 -> z-score undefined -> NODATA."""
        n = 8
        dates = _dates(n)
        arr = np.full((n, 3, 3), 0.5, dtype=np.float64)
        z = anomaly_zscore(arr, dates, window_days=200)
        assert np.all((z == NODATA) | (z == 0.0))

    def test_mean_of_z_scores_near_zero(self):
        """For a series without outliers, mean |z| should be small."""
        rng = np.random.default_rng(0)
        n = 20
        dates = _dates(n)
        arr = rng.normal(0.5, 0.1, (n, 8, 8))
        z = anomaly_zscore(arr, dates, window_days=400)
        valid = z != NODATA
        assert np.abs(z[valid].mean()) < 1.0

    def test_precomputed_background_accepted(self):
        n = 8
        dates = _dates(n)
        arr = _stack(n, 4, 4)
        bg = time_window_stats(arr, dates, window_days=120)
        z1 = anomaly_zscore(arr, dates, background=bg)
        z2 = anomaly_zscore(arr, dates, window_days=120)
        np.testing.assert_allclose(z1, z2, atol=1e-10)

    def test_date_mismatch_raises(self):
        arr = _stack(5, 4, 4)
        with pytest.raises(ValueError):
            anomaly_zscore(arr, _dates(3), window_days=60)


# trend_theil_sen


@skip_no_lib
class TestTrendTheilSen:

    def test_output_keys(self):
        n = 8
        dates = _dates(n)
        arr = _stack(n, 4, 4)
        res = trend_theil_sen(arr, dates)
        assert set(res.keys()) == {"slope", "intercept"}

    def test_output_shapes(self):
        n = 8
        dates = _dates(n)
        arr = _stack(n, 6, 10)
        res = trend_theil_sen(arr, dates)
        assert res["slope"].shape == (6, 10)
        assert res["intercept"].shape == (6, 10)

    def test_all_nodata_returns_nodata(self):
        n = 6
        dates = _dates(n)
        arr = _nodata_stack(n, 3, 3)
        res = trend_theil_sen(arr, dates)
        assert np.all(res["slope"] == NODATA)

    def test_linear_slope_recovered(self):
        n = 10
        dates = _dates(n, step_days=16.0)
        slope = 0.002
        arr = _linear_pixel(n, slope=slope, intercept=0.3, dates=dates)
        res = trend_theil_sen(arr, dates)
        assert res["slope"][0, 0] == pytest.approx(slope, rel=1e-4)

    def test_flat_series_slope_zero(self):
        n = 10
        dates = _dates(n)
        arr = np.full((n, 3, 3), 0.6, dtype=np.float64)
        res = trend_theil_sen(arr, dates)
        np.testing.assert_allclose(res["slope"], 0.0, atol=1e-8)

    def test_robust_to_single_outlier(self):
        """Theil-Sen slope should be close to true slope even with an outlier."""
        n = 12
        dates = _dates(n, step_days=16.0)
        slope = 0.001
        arr = _linear_pixel(n, slope=slope, intercept=0.2, dates=dates)
        arr[5, 0, 0] = 10.0  # extreme outlier
        res = trend_theil_sen(arr, dates)
        assert res["slope"][0, 0] == pytest.approx(slope, abs=0.002)


# mann_kendall


@skip_no_lib
class TestMannKendall:

    def test_output_keys(self):
        n = 8
        dates = _dates(n)
        arr = _stack(n, 4, 4)
        res = mann_kendall(arr, dates)
        assert set(res.keys()) == {"S", "varS", "Z", "trend"}

    def test_output_shapes(self):
        n = 8
        dates = _dates(n)
        arr = _stack(n, 6, 10)
        res = mann_kendall(arr, dates)
        for v in res.values():
            assert v.shape == (6, 10)

    def test_fewer_than_4_obs_returns_nodata(self):
        n = 3
        dates = _dates(n)
        arr = _stack(n, 3, 3)
        res = mann_kendall(arr, dates)
        assert np.all(res["trend"] == NODATA)

    def test_all_nodata_returns_nodata(self):
        n = 8
        dates = _dates(n)
        arr = _nodata_stack(n, 3, 3)
        res = mann_kendall(arr, dates)
        assert np.all(res["trend"] == NODATA)

    def test_monotone_increasing_detected(self):
        n = 12
        dates = _dates(n, step_days=16.0)
        arr = _linear_pixel(n, slope=0.01, intercept=0.1, dates=dates)
        res = mann_kendall(arr, dates)
        assert res["trend"][0, 0] == pytest.approx(1.0)

    def test_monotone_decreasing_detected(self):
        n = 12
        dates = _dates(n, step_days=16.0)
        arr = _linear_pixel(n, slope=-0.01, intercept=0.9, dates=dates)
        res = mann_kendall(arr, dates)
        assert res["trend"][0, 0] == pytest.approx(-1.0)

    def test_trend_values_in_allowed_set(self):
        n = 10
        dates = _dates(n)
        arr = _stack(n, 8, 8, seed=5)
        res = mann_kendall(arr, dates)
        allowed = {-1.0, 0.0, 1.0, NODATA}
        unique = set(np.unique(res["trend"]))
        assert unique.issubset(allowed)

    def test_vars_non_negative(self):
        n = 10
        dates = _dates(n)
        arr = _stack(n, 4, 4)
        res = mann_kendall(arr, dates)
        valid = res["varS"] != NODATA
        assert np.all(res["varS"][valid] >= 0.0)


# bfast_breakpoint

@skip_no_lib
class TestBfastBreakpoint:

    def test_output_keys(self):
        n = 12
        dates = _dates(n)
        arr = _stack(n, 4, 4)
        res = bfast_breakpoint(arr, dates)
        assert set(res.keys()) == {"break_day", "magnitude", "rss_ratio"}

    def test_output_shapes(self):
        n = 12
        dates = _dates(n)
        arr = _stack(n, 6, 8)
        res = bfast_breakpoint(arr, dates)
        for v in res.values():
            assert v.shape == (6, 8)

    def test_fewer_than_6_obs_returns_nodata(self):
        n = 4
        dates = _dates(n)
        arr = _stack(n, 3, 3)
        res = bfast_breakpoint(arr, dates)
        assert np.all(res["break_day"] == NODATA)

    def test_all_nodata_returns_nodata(self):
        n = 12
        dates = _dates(n)
        arr = _nodata_stack(n, 3, 3)
        res = bfast_breakpoint(arr, dates)
        assert np.all(res["break_day"] == NODATA)

    def test_rss_ratio_in_0_1_range(self):
        n = 12
        dates = _dates(n, step_days=16.0)
        arr = _stack(n, 4, 4)
        res = bfast_breakpoint(arr, dates)
        valid = res["rss_ratio"] != NODATA
        assert np.all(res["rss_ratio"][valid] >= 0.0)
        assert np.all(res["rss_ratio"][valid] <= 1.0 + 1e-6)

    def test_sharp_break_has_low_rss_ratio(self):
        """A step-change should give rss_ratio well below 1."""
        n = 12
        dates = _dates(n, step_days=16.0)
        arr = np.zeros((n, 1, 1), dtype=np.float64)
        arr[:6, 0, 0] = 0.2
        arr[6:, 0, 0] = 0.8
        res = bfast_breakpoint(arr, dates)
        assert res["rss_ratio"][0, 0] < 0.5

    def test_magnitude_sign_reflects_direction(self):
        """Jump up -> positive magnitude."""
        n = 12
        dates = _dates(n, step_days=16.0)
        arr = np.zeros((n, 1, 1), dtype=np.float64)
        arr[:6, 0, 0] = 0.2
        arr[6:, 0, 0] = 0.8
        res = bfast_breakpoint(arr, dates)
        assert res["magnitude"][0, 0] > 0.0


# pixel_regression


@skip_no_lib
class TestPixelRegression:

    def test_output_keys(self):
        n = 10
        x = np.arange(n, dtype=np.float64)
        arr = _stack(n, 4, 4)
        res = pixel_regression(arr, x)
        assert set(res.keys()) == {"slope", "intercept", "r2"}

    def test_output_shapes(self):
        n = 10
        x = np.arange(n, dtype=np.float64)
        arr = _stack(n, 6, 8)
        res = pixel_regression(arr, x)
        for v in res.values():
            assert v.shape == (6, 8)

    def test_perfect_linear_slope_and_r2(self):
        n = 10
        x = np.arange(n, dtype=np.float64)
        arr = np.zeros((n, 1, 1), dtype=np.float64)
        arr[:, 0, 0] = 2.0 * x + 1.0
        res = pixel_regression(arr, x)
        assert res["slope"][0, 0] == pytest.approx(2.0, rel=1e-5)
        assert res["intercept"][0, 0] == pytest.approx(1.0, rel=1e-5)
        assert res["r2"][0, 0] == pytest.approx(1.0, rel=1e-5)

    def test_flat_series_slope_zero_r2_zero(self):
        n = 10
        x = np.arange(n, dtype=np.float64)
        arr = np.full((n, 3, 3), 5.0, dtype=np.float64)
        res = pixel_regression(arr, x)
        np.testing.assert_allclose(res["slope"], 0.0, atol=1e-8)
        np.testing.assert_allclose(res["r2"], 0.0, atol=1e-8)

    def test_all_nodata_returns_nodata(self):
        n = 10
        x = np.arange(n, dtype=np.float64)
        arr = _nodata_stack(n, 3, 3)
        res = pixel_regression(arr, x)
        assert np.all(res["slope"] == NODATA)
        assert np.all(res["r2"] == NODATA)

    def test_r2_in_0_1_range(self):
        rng = np.random.default_rng(1)
        n = 15
        x = np.arange(n, dtype=np.float64)
        arr = rng.random((n, 8, 8))
        res = pixel_regression(arr, x)
        valid = res["r2"] != NODATA
        assert np.all(res["r2"][valid] >= 0.0)
        assert np.all(res["r2"][valid] <= 1.0 + 1e-9)

    def test_x_length_mismatch_raises(self):
        arr = _stack(5, 4, 4)
        x = np.arange(3, dtype=np.float64)
        with pytest.raises(ValueError):
            pixel_regression(arr, x)

    def test_negative_slope_recovered(self):
        n = 10
        x = np.arange(n, dtype=np.float64)
        arr = np.zeros((n, 1, 1), dtype=np.float64)
        arr[:, 0, 0] = 1.0 - 0.05 * x
        res = pixel_regression(arr, x)
        assert res["slope"][0, 0] == pytest.approx(-0.05, rel=1e-4)

    def test_one_nodata_pixel_rest_valid(self):
        n = 10
        x = np.arange(n, dtype=np.float64)
        arr = _stack(n, 2, 2)
        arr[:, 1, 1] = NODATA
        res = pixel_regression(arr, x)
        assert res["slope"][0, 0] != NODATA  # valid pixel computed
        assert res["slope"][1, 1] == NODATA  # all-NODATA pixel → NODATA


# pearson_map


@skip_no_lib
class TestPearsonMap:

    def test_output_shape(self):
        n = 8
        dates = _dates(n)
        a = _stack(n, 6, 8)
        b = _stack(n, 6, 8, seed=99)
        r = pearson_map(a, b, dates)
        assert r.shape == (6, 8)

    def test_perfect_positive_correlation(self):
        n = 10
        dates = _dates(n)
        arr = _stack(n, 4, 4)
        r = pearson_map(arr, arr * 2.0 + 1.0, dates)
        np.testing.assert_allclose(r, 1.0, atol=1e-6)

    def test_perfect_negative_correlation(self):
        n = 10
        dates = _dates(n)
        arr = _stack(n, 4, 4)
        r = pearson_map(arr, -arr, dates)
        np.testing.assert_allclose(r, -1.0, atol=1e-6)

    def test_all_nodata_returns_nodata(self):
        n = 6
        dates = _dates(n)
        a = _nodata_stack(n, 3, 3)
        b = _nodata_stack(n, 3, 3)
        r = pearson_map(a, b, dates)
        assert np.all(r == NODATA)

    def test_r_in_minus1_1_range(self):
        rng = np.random.default_rng(2)
        n = 12
        dates = _dates(n)
        a = rng.random((n, 8, 8))
        b = rng.random((n, 8, 8))
        r = pearson_map(a, b, dates)
        valid = r != NODATA
        assert np.all(r[valid] >= -1.0 - 1e-9)
        assert np.all(r[valid] <= 1.0 + 1e-9)

    def test_shape_mismatch_raises(self):
        n = 6
        dates = _dates(n)
        a = _stack(n, 4, 4)
        b = _stack(n, 6, 6)
        with pytest.raises(ValueError):
            pearson_map(a, b, dates)

    def test_fewer_than_3_joint_obs_returns_nodata(self):
        n = 6
        dates = _dates(n)
        a = np.full((n, 1, 1), NODATA, dtype=np.float64)
        b = np.full((n, 1, 1), NODATA, dtype=np.float64)
        # only 2 valid joint steps
        a[0, 0, 0] = 0.3
        b[0, 0, 0] = 0.4
        a[1, 0, 0] = 0.5
        b[1, 0, 0] = 0.6
        r = pearson_map(a, b, dates)
        assert r[0, 0] == NODATA


# phenology_doy


@skip_no_lib
class TestPhenologyDoy:

    def _bell_stack(self, n: int = 14, rows: int = 1, cols: int = 1) -> tuple:
        dates = _dates(n, step_days=30.0)
        curve = [0.1, 0.15, 0.2, 0.4, 0.6, 0.75, 0.85,
                 0.9, 0.8, 0.65, 0.45, 0.3, 0.2, 0.1]
        arr = np.zeros((n, rows, cols), dtype=np.float64)
        for t, v in enumerate(curve):
            arr[t, :, :] = v
        return arr, dates

    def test_output_keys(self):
        arr, dates = self._bell_stack()
        res = phenology_doy(arr, dates)
        assert set(res.keys()) == {"sos", "pos", "eos"}

    def test_pos_at_peak(self):
        arr, dates = self._bell_stack()
        d = np.array([(dt - dates[0]).total_seconds() / 86_400.0 for dt in dates])
        res = phenology_doy(arr, dates)
        assert res["pos"][0, 0] == pytest.approx(d[7], abs=1.0)

    def test_sos_before_pos(self):
        arr, dates = self._bell_stack()
        res = phenology_doy(arr, dates)
        sos, pos = res["sos"][0, 0], res["pos"][0, 0]
        if sos != NODATA and pos != NODATA:
            assert sos < pos

    def test_eos_after_pos(self):
        arr, dates = self._bell_stack()
        res = phenology_doy(arr, dates)
        pos, eos = res["pos"][0, 0], res["eos"][0, 0]
        if pos != NODATA and eos != NODATA:
            assert eos > pos

    def test_flat_series_returns_nodata(self):
        n = 10
        dates = _dates(n)
        arr = np.full((n, 1, 1), 0.5, dtype=np.float64)
        res = phenology_doy(arr, dates)
        assert res["pos"][0, 0] == NODATA

    def test_output_shapes(self):
        arr, dates = self._bell_stack(rows=4, cols=6)
        res = phenology_doy(arr, dates)
        for v in res.values():
            assert v.shape == (4, 6)


# phenology_metrics

@skip_no_lib
class TestPhenologyMetrics:

    def _ndvi_stack(self, n: int = 20) -> tuple:
        dates = _dates(n, step_days=18.0)
        curve = [0.1, 0.12, 0.18, 0.3, 0.5, 0.68, 0.78, 0.85,
                 0.88, 0.87, 0.82, 0.72, 0.58, 0.42, 0.3,
                 0.22, 0.16, 0.13, 0.11, 0.1]
        arr = np.zeros((n, 3, 3), dtype=np.float64)
        for t, v in enumerate(curve):
            arr[t, :, :] = v
        return arr, dates

    def test_output_keys(self):
        arr, dates = self._ndvi_stack()
        res = phenology_metrics(arr, dates)
        assert set(res.keys()) == {"sos_doy", "eos_doy", "peak_doy", "peak_val"}

    def test_output_shapes(self):
        arr, dates = self._ndvi_stack()
        res = phenology_metrics(arr, dates)
        for v in res.values():
            assert v.shape == (3, 3)

    def test_peak_val_in_0_1_for_ndvi(self):
        arr, dates = self._ndvi_stack()
        res = phenology_metrics(arr, dates)
        valid = res["peak_val"] != NODATA
        assert np.all(res["peak_val"][valid] <= 1.0 + 1e-6)
        assert np.all(res["peak_val"][valid] >= -1.0 - 1e-6)

    def test_sos_before_eos(self):
        arr, dates = self._ndvi_stack()
        res = phenology_metrics(arr, dates)
        valid = (res["sos_doy"] != NODATA) & (res["eos_doy"] != NODATA)
        assert np.all(res["sos_doy"][valid] < res["eos_doy"][valid])

    def test_fewer_than_5_valid_returns_nodata(self):
        n = 3
        dates = _dates(n)
        arr = np.full((n, 2, 2), 0.5, dtype=np.float64)
        res = phenology_metrics(arr, dates)
        assert np.all(res["peak_doy"] == NODATA)

    def test_smooth_false_accepted(self):
        arr, dates = self._ndvi_stack()
        res = phenology_metrics(arr, dates, smooth=False)
        assert res["peak_doy"].shape == (3, 3)

    def test_xarray_input_accepted(self):
        pytest.importorskip("xarray")
        import xarray as xr
        arr, dates = self._ndvi_stack()
        da = xr.DataArray(arr, dims=["time", "y", "x"])
        res = phenology_metrics(da, dates)
        assert res["peak_doy"].shape == (3, 3)

    def test_rising_pct_0_sos_at_first_obs(self):
        """rising_pct=0 means threshold = min → SOS should be very early."""
        arr, dates = self._ndvi_stack()
        res0 = phenology_metrics(arr, dates, rising_pct=0, smooth=False)
        res20 = phenology_metrics(arr, dates, rising_pct=20, smooth=False)
        valid = (res0["sos_doy"] != NODATA) & (res20["sos_doy"] != NODATA)
        assert np.all(res0["sos_doy"][valid] <= res20["sos_doy"][valid] + 1e-6)


# save_phenology


@skip_no_lib
class TestSavePhenology:

    def _make_metrics(self, rows: int = 4, cols: int = 6) -> dict:
        rng = np.random.default_rng(0)
        return {
            "sos_doy": rng.random((rows, cols)) * 30,
            "eos_doy": rng.random((rows, cols)) * 30 + 150,
            "peak_doy": rng.random((rows, cols)) * 30 + 80,
            "peak_val": rng.random((rows, cols)),
        }

    def test_four_tifs_written(self, tmp_path):
        pytest.importorskip("rasterio")
        metrics = self._make_metrics()
        saved = save_phenology(metrics, tmp_path)
        assert len(saved) == 4
        for p in saved.values():
            assert p.exists()
            assert p.suffix == ".tif"

    def test_output_dir_created(self, tmp_path):
        pytest.importorskip("rasterio")
        metrics = self._make_metrics()
        out = tmp_path / "new_subdir" / "phenology"
        save_phenology(metrics, out)
        assert out.is_dir()

    def test_returned_keys_match_expected(self, tmp_path):
        pytest.importorskip("rasterio")
        metrics = self._make_metrics()
        saved = save_phenology(metrics, tmp_path)
        assert set(saved.keys()) == {"sos_doy", "eos_doy", "peak_doy", "peak_val"}

    def test_round_trip_values(self, tmp_path):
        """Values written to GeoTIFF match the input array."""
        rasterio = pytest.importorskip("rasterio")
        metrics = self._make_metrics()
        saved = save_phenology(metrics, tmp_path)
        with rasterio.open(saved["peak_val"]) as src:
            data = src.read(1)
        np.testing.assert_allclose(
            data.astype(np.float64),
            metrics["peak_val"].astype(np.float64),
            atol=1e-6,
        )

    def test_missing_rasterio_raises_import_error(self, tmp_path):
        metrics = self._make_metrics()
        with patch.dict("sys.modules", {"rasterio": None}):
            with pytest.raises((ImportError, TypeError)):
                save_phenology(metrics, tmp_path)
