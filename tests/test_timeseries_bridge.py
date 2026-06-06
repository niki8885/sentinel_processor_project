import numpy as np
import pytest

from sentinel_processor.processing._timeseries_bridge import (
    NODATA,
    _METHOD_MAP,
    _prepare_arr,
    _prepare_mask,
    interpolate_gaps,
)

pytestmark = pytest.mark.fortran

ALL_METHODS = list(_METHOD_MAP.keys())


# ---------------------------------------------------------------------------
# Helpers shared across tests
# ---------------------------------------------------------------------------

def _ramp(n_times=10, rows=4, cols=4, seed=0):
    """Linear ramp along time axis, all pixels identical."""
    rng = np.random.default_rng(seed)
    base = np.linspace(100.0, 900.0, n_times)
    arr = np.broadcast_to(base[:, None, None], (n_times, rows, cols)).copy()
    return arr.astype(np.float64)


def _full_mask(shape):
    return np.ones(shape, dtype=np.int32)


def _gap_mask(shape, gap_times):
    m = np.ones(shape, dtype=np.int32)
    for t in gap_times:
        m[t] = 0
    return m


# ---------------------------------------------------------------------------
# _prepare_arr
# ---------------------------------------------------------------------------

class TestPrepareArr:

    def test_returns_float64(self):
        arr = np.ones((5, 4, 4), dtype=np.float32)
        _, buf = _prepare_arr(arr)
        assert buf.dtype == np.float64

    def test_returns_c_contiguous(self):
        arr = np.asfortranarray(np.ones((5, 4, 4)))
        _, buf = _prepare_arr(arr)
        assert buf.flags["C_CONTIGUOUS"]

    def test_returns_copy(self):
        arr = np.ones((5, 4, 4), dtype=np.float64)
        _, buf = _prepare_arr(arr)
        buf[0, 0, 0] = 999.0
        assert arr[0, 0, 0] == 1.0

    def test_values_preserved(self):
        arr = np.arange(60, dtype=np.float32).reshape(3, 4, 5)
        _, buf = _prepare_arr(arr)
        np.testing.assert_allclose(buf, arr.astype(np.float64))


# ---------------------------------------------------------------------------
# _prepare_mask
# ---------------------------------------------------------------------------

class TestPrepareMask:

    def test_returns_int32(self):
        m = np.ones((5, 4, 4), dtype=bool)
        _, buf = _prepare_mask(m, (5, 4, 4))
        assert buf.dtype == np.int32

    def test_shape_mismatch_raises(self):
        m = np.ones((5, 4, 4), dtype=np.int32)
        with pytest.raises(ValueError, match="shape"):
            _prepare_mask(m, (6, 4, 4))

    def test_returns_c_contiguous(self):
        m = np.asfortranarray(np.ones((5, 4, 4), dtype=np.int32))
        _, buf = _prepare_mask(m, (5, 4, 4))
        assert buf.flags["C_CONTIGUOUS"]

    def test_bool_mask_converted(self):
        m = np.array([[[True, False]], [[False, True]]], dtype=bool)
        _, buf = _prepare_mask(m, m.shape)
        assert buf[0, 0, 0] == 1
        assert buf[0, 0, 1] == 0


# ---------------------------------------------------------------------------
# interpolate_gaps — validation (Python layer, no Fortran needed)
# ---------------------------------------------------------------------------

class TestInterpolateGapsValidation:

    def test_wrong_ndim_raises(self):
        arr  = np.ones((5, 4))
        mask = np.ones((5, 4), dtype=np.int32)
        with pytest.raises(ValueError, match="3-D"):
            interpolate_gaps(arr, mask)

    def test_unknown_method_raises(self):
        arr  = np.ones((5, 4, 4))
        mask = np.ones((5, 4, 4), dtype=np.int32)
        with pytest.raises(ValueError, match="Unknown method"):
            interpolate_gaps(arr, mask, method="spline")

    def test_mask_shape_mismatch_raises(self):
        arr  = np.ones((5, 4, 4))
        mask = np.ones((4, 4, 4), dtype=np.int32)
        with pytest.raises(ValueError, match="shape"):
            interpolate_gaps(arr, mask)

    def test_original_array_not_mutated(self):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[3, 4])
        orig = arr.copy()
        interpolate_gaps(arr, mask, method="linear")
        np.testing.assert_array_equal(arr, orig)

    def test_even_window_forced_odd(self):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[4])
        out_even = interpolate_gaps(arr, mask, method="savgol", window=4)
        out_odd  = interpolate_gaps(arr, mask, method="savgol", window=5)
        np.testing.assert_allclose(out_even, out_odd)

    def test_window_below_3_forced_to_3(self):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[4])
        out_1 = interpolate_gaps(arr, mask, method="savgol", window=1)
        out_3 = interpolate_gaps(arr, mask, method="savgol", window=3)
        np.testing.assert_allclose(out_1, out_3)

    def test_returns_float64(self):
        arr  = _ramp().astype(np.float32)
        mask = _full_mask(arr.shape)
        out  = interpolate_gaps(arr, mask, method="linear")
        assert out.dtype == np.float64

    def test_output_shape_preserved(self):
        arr  = _ramp(n_times=8, rows=6, cols=7)
        mask = _full_mask(arr.shape)
        out  = interpolate_gaps(arr, mask)
        assert out.shape == arr.shape


# ---------------------------------------------------------------------------
# All methods: shared contract
# ---------------------------------------------------------------------------

class TestAllMethodsContract:

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_no_gaps_unchanged(self, method):
        """When mask is all-valid, output equals input."""
        arr  = _ramp()
        mask = _full_mask(arr.shape)
        out  = interpolate_gaps(arr, mask, method=method)
        np.testing.assert_allclose(out, arr)

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_all_masked_returns_nodata(self, method):
        arr  = _ramp()
        mask = np.zeros(arr.shape, dtype=np.int32)
        out  = interpolate_gaps(arr, mask, method=method)
        np.testing.assert_array_equal(out, NODATA)

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_no_nan_in_output(self, method):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[2, 5, 7])
        out  = interpolate_gaps(arr, mask, method=method)
        assert not np.isnan(out).any()

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_no_inf_in_output(self, method):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[0, 9])
        out  = interpolate_gaps(arr, mask, method=method)
        assert not np.isinf(out).any()

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_valid_positions_unchanged(self, method):
        """Positions where mask==1 must keep their original values."""
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[3, 6])
        out  = interpolate_gaps(arr, mask, method=method)
        valid = mask.astype(bool)
        np.testing.assert_allclose(out[valid], arr[valid])

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_single_valid_obs_constant_fill(self, method):
        arr  = _ramp()
        mask = np.zeros(arr.shape, dtype=np.int32)
        mask[5] = 1          # only one valid time step
        out  = interpolate_gaps(arr, mask, method=method)
        expected = arr[5]    # all pixels should equal t=5 value
        for t in range(arr.shape[0]):
            if mask[t, 0, 0] == 0:
                np.testing.assert_allclose(out[t], expected, atol=1e-6)

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_output_shape_matches_input(self, method):
        arr  = _ramp(n_times=12, rows=5, cols=7)
        mask = _gap_mask(arr.shape, gap_times=[1, 4, 9])
        out  = interpolate_gaps(arr, mask, method=method)
        assert out.shape == arr.shape

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_leading_gap_filled_with_first_valid(self, method):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[0, 1, 2])
        out  = interpolate_gaps(arr, mask, method=method)
        first_valid = arr[3]
        np.testing.assert_allclose(out[0], first_valid, atol=1e-4)

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_trailing_gap_filled_with_last_valid(self, method):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[8, 9])
        out  = interpolate_gaps(arr, mask, method=method)
        last_valid = arr[7]
        np.testing.assert_allclose(out[9], last_valid, atol=1e-4)


# ---------------------------------------------------------------------------
# METHOD 0 — linear
# ---------------------------------------------------------------------------

class TestLinear:

    def test_midpoint_of_ramp(self):
        """Gap at t=5 in a linear ramp -> must equal the ramp value exactly."""
        arr  = _ramp(n_times=10)
        mask = _gap_mask(arr.shape, gap_times=[5])
        out  = interpolate_gaps(arr, mask, method="linear")
        np.testing.assert_allclose(out[5], arr[5], atol=1e-6)

    def test_interior_gap_interpolated_linearly(self):
        """t=2 gap between t=1=200 and t=3=400 -> t=2 must equal 300."""
        arr       = np.zeros((5, 2, 2))
        arr[1]    = 200.0
        arr[3]    = 400.0
        arr[4]    = 500.0
        mask      = np.ones((5, 2, 2), dtype=np.int32)
        mask[2]   = 0
        out = interpolate_gaps(arr, mask, method="linear")
        np.testing.assert_allclose(out[2], 300.0, atol=1e-6)

    def test_multi_step_gap_linear(self):
        """Three consecutive gaps between t=0=0 and t=4=400."""
        arr     = np.zeros((5, 1, 1))
        arr[0]  = 0.0
        arr[4]  = 400.0
        mask    = np.ones((5, 1, 1), dtype=np.int32)
        mask[1] = mask[2] = mask[3] = 0
        out = interpolate_gaps(arr, mask, method="linear")
        np.testing.assert_allclose(out[1, 0, 0], 100.0, atol=1e-6)
        np.testing.assert_allclose(out[2, 0, 0], 200.0, atol=1e-6)
        np.testing.assert_allclose(out[3, 0, 0], 300.0, atol=1e-6)

    def test_constant_series_gap_stays_constant(self):
        arr  = np.full((8, 3, 3), 42.0)
        mask = _gap_mask(arr.shape, gap_times=[2, 3, 5])
        out  = interpolate_gaps(arr, mask, method="linear")
        np.testing.assert_allclose(out, 42.0, atol=1e-6)


# ---------------------------------------------------------------------------
# METHOD 1 — savgol
# ---------------------------------------------------------------------------

class TestSavgol:

    def test_ramp_reconstructed(self):
        """S-G on a linear ramp must recover the missing value accurately."""
        arr  = _ramp(n_times=12)
        mask = _gap_mask(arr.shape, gap_times=[5])
        out  = interpolate_gaps(arr, mask, method="savgol", window=5)
        np.testing.assert_allclose(out[5], arr[5], atol=5.0)

    def test_output_in_plausible_range(self):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[3, 7])
        out  = interpolate_gaps(arr, mask, method="savgol", window=5)
        assert out.min() > arr.min() - 200
        assert out.max() < arr.max() + 200

    def test_window_7_closer_to_trend_than_window_3(self):
        """window=7 S-G should produce filled values closer to the underlying
        linear trend than window=3 on a heavily noisy series (averaged over
        many pixels to get a stable estimate)."""
        rng   = np.random.default_rng(42)
        rows, cols = 16, 16
        base  = np.linspace(0.0, 1000.0, 20)
        trend = np.broadcast_to(base[:, None, None], (20, rows, cols)).copy()
        noise = rng.normal(0, 80, (20, rows, cols))
        arr   = trend + noise
        mask  = _gap_mask(arr.shape, gap_times=[5, 10, 15])
        out3  = interpolate_gaps(arr, mask, method="savgol", window=3)
        out7  = interpolate_gaps(arr, mask, method="savgol", window=7)
        gap_idx = [5, 10, 15]
        mae3 = np.mean([np.abs(out3[t] - trend[t]) for t in gap_idx])
        mae7 = np.mean([np.abs(out7[t] - trend[t]) for t in gap_idx])
        assert mae7 < mae3


# ---------------------------------------------------------------------------
# METHOD 2 — pchip
# ---------------------------------------------------------------------------

class TestPchip:

    def test_no_overshoot_on_step(self):
        """PCHIP must not exceed the boundary values on a step-like series."""
        arr       = np.zeros((8, 2, 2))
        arr[:4]   = 100.0
        arr[4:]   = 800.0
        mask      = np.ones((8, 2, 2), dtype=np.int32)
        mask[3]   = 0       # gap right at the step
        out = interpolate_gaps(arr, mask, method="pchip")
        assert out[3].min() >= 99.0
        assert out[3].max() <= 801.0

    def test_ramp_exact(self):
        arr  = _ramp(n_times=10)
        mask = _gap_mask(arr.shape, gap_times=[4])
        out  = interpolate_gaps(arr, mask, method="pchip")
        np.testing.assert_allclose(out[4], arr[4], atol=1e-4)

    def test_monotone_between_knots(self):
        """Values between two monotone knots must not exceed them."""
        arr     = np.zeros((6, 1, 1))
        arr[0]  = 200.0
        arr[5]  = 800.0
        mask    = np.ones((6, 1, 1), dtype=np.int32)
        mask[1] = mask[2] = mask[3] = mask[4] = 0
        out = interpolate_gaps(arr, mask, method="pchip")
        for t in range(1, 5):
            assert out[t, 0, 0] >= 200.0 - 1e-3
            assert out[t, 0, 0] <= 800.0 + 1e-3

    def test_requires_2_knots_for_cubic(self):
        """With only 1 valid observation, falls back to constant fill."""
        arr  = _ramp()
        mask = np.zeros(arr.shape, dtype=np.int32)
        mask[4] = 1
        out = interpolate_gaps(arr, mask, method="pchip")
        assert not np.isnan(out).any()
        assert not np.isinf(out).any()


# ---------------------------------------------------------------------------
# METHOD 3 — ets
# ---------------------------------------------------------------------------

class TestEts:

    def test_falls_back_to_linear_for_two_valid(self):
        """With <= 2 valid obs ETS falls back to linear — must not crash."""
        arr  = _ramp()
        mask = np.zeros(arr.shape, dtype=np.int32)
        mask[2] = mask[7] = 1
        out = interpolate_gaps(arr, mask, method="ets")
        assert not np.isnan(out).any()

    def test_trending_series_gap_in_range(self):
        """Gap in a steadily rising series must be near the trend."""
        arr     = np.zeros((10, 1, 1))
        for t in range(10):
            arr[t, 0, 0] = t * 100.0
        mask    = _gap_mask(arr.shape, gap_times=[5])
        out     = interpolate_gaps(arr, mask, method="ets")
        assert 0.0 <= out[5, 0, 0] <= 1000.0

    def test_no_crash_on_flat_series(self):
        arr  = np.full((8, 3, 3), 500.0)
        mask = _gap_mask(arr.shape, gap_times=[3])
        out  = interpolate_gaps(arr, mask, method="ets")
        assert not np.isnan(out).any()


# ---------------------------------------------------------------------------
# METHOD 4 — gauss
# ---------------------------------------------------------------------------

class TestGauss:

    def test_ramp_approx_recovered(self):
        arr  = _ramp(n_times=12)
        mask = _gap_mask(arr.shape, gap_times=[5])
        out  = interpolate_gaps(arr, mask, method="gauss", window=5)
        np.testing.assert_allclose(out[5], arr[5], atol=30.0)

    def test_wider_window_smoother(self):
        rng  = np.random.default_rng(7)
        base = np.linspace(0, 1000, 20)
        arr  = base[:, None, None] + rng.normal(0, 80, (20, 4, 4))
        mask = _gap_mask(arr.shape, gap_times=[6, 12])
        out5  = interpolate_gaps(arr, mask, method="gauss", window=5)
        out11 = interpolate_gaps(arr, mask, method="gauss", window=11)
        std5  = np.std([out5[t]  for t in [6, 12]])
        std11 = np.std([out11[t] for t in [6, 12]])
        assert std11 <= std5 + 1.0

    def test_constant_series_stays_constant(self):
        arr  = np.full((10, 3, 3), 250.0)
        mask = _gap_mask(arr.shape, gap_times=[2, 7])
        out  = interpolate_gaps(arr, mask, method="gauss", window=7)
        np.testing.assert_allclose(out, 250.0, atol=1e-4)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:

    def test_single_time_step_all_valid(self):
        arr  = np.ones((1, 4, 4)) * 500.0
        mask = _full_mask(arr.shape)
        out  = interpolate_gaps(arr, mask, method="linear")
        np.testing.assert_allclose(out, 500.0)

    def test_single_time_step_all_masked(self):
        arr  = np.ones((1, 4, 4)) * 500.0
        mask = np.zeros((1, 4, 4), dtype=np.int32)
        out  = interpolate_gaps(arr, mask, method="linear")
        np.testing.assert_array_equal(out, NODATA)

    def test_large_array_no_crash(self):
        arr  = _ramp(n_times=20, rows=64, cols=64)
        mask = _gap_mask(arr.shape, gap_times=[5, 10, 15])
        out  = interpolate_gaps(arr, mask, method="pchip")
        assert out.shape == arr.shape
        assert not np.isnan(out).any()

    def test_all_gaps_except_first_and_last(self):
        arr     = np.zeros((6, 2, 2))
        arr[0]  = 100.0
        arr[5]  = 600.0
        mask    = np.zeros((6, 2, 2), dtype=np.int32)
        mask[0] = mask[5] = 1
        out = interpolate_gaps(arr, mask, method="linear")
        np.testing.assert_allclose(out[0], 100.0, atol=1e-6)
        np.testing.assert_allclose(out[5], 600.0, atol=1e-6)
        assert not np.isnan(out).any()

    def test_float32_input_accepted(self):
        arr  = _ramp().astype(np.float32)
        mask = _gap_mask(arr.shape, gap_times=[3])
        out  = interpolate_gaps(arr, mask, method="linear")
        assert out.dtype == np.float64

    def test_bool_mask_accepted(self):
        arr  = _ramp()
        mask = _gap_mask(arr.shape, gap_times=[4]).astype(bool)
        out  = interpolate_gaps(arr, mask, method="linear")
        assert not np.isnan(out).any()

    def test_per_pixel_independence(self):
        """Different pixels can have gaps at different positions."""
        arr        = np.zeros((6, 2, 1))
        arr[:, 0, 0] = np.linspace(0, 500, 6)
        arr[:, 1, 0] = np.linspace(500, 0, 6)
        mask       = np.ones((6, 2, 1), dtype=np.int32)
        mask[2, 0, 0] = 0   # gap only in pixel 0
        mask[4, 1, 0] = 0   # gap only in pixel 1
        out = interpolate_gaps(arr, mask, method="linear")
        assert not np.isnan(out).any()
        np.testing.assert_allclose(out[2, 0, 0], arr[2, 0, 0], atol=1e-4)
        np.testing.assert_allclose(out[4, 1, 0], arr[4, 1, 0], atol=1e-4)

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_two_time_steps_no_gap(self, method):
        arr  = np.array([[[100.0]], [[200.0]]])
        mask = _full_mask(arr.shape)
        out  = interpolate_gaps(arr, mask, method=method)
        np.testing.assert_allclose(out, arr)

    @pytest.mark.parametrize("method", ALL_METHODS)
    def test_nodata_constant_is_minus_9999(self, method):
        arr  = _ramp()
        mask = np.zeros(arr.shape, dtype=np.int32)
        out  = interpolate_gaps(arr, mask, method=method)
        assert (out == NODATA).all()