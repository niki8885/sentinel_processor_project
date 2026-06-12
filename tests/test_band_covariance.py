from __future__ import annotations
import numpy as np
import pytest

try:
    from sentinel_processor.analysis._band_covariance_bridge import (
        NODATA,
        band_covariance,
    )

    _BRIDGE_AVAILABLE = True
except (ImportError, FileNotFoundError, OSError):
    _BRIDGE_AVAILABLE = False

pytestmark = pytest.mark.fortran
skip_no_lib = pytest.mark.skipif(
    not _BRIDGE_AVAILABLE,
    reason="libband_covariance not compiled — run python build_band_covariance.py",
)


def _cube(n_bands: int, rows: int, cols: int, seed: int = 0) -> np.ndarray:
    """Random float64 cube (n_bands, rows, cols) with values in [0, 1]."""
    return np.random.default_rng(seed).random((n_bands, rows, cols))


def _nodata_cube(n_bands: int, rows: int, cols: int) -> np.ndarray:
    """Cube where every pixel is NODATA."""
    return np.full((n_bands, rows, cols), NODATA, dtype=np.float64)


def _linear_cube(n_bands: int, rows: int, cols: int) -> np.ndarray:
    """Cube with a linear ramp per band — each band is perfectly correlated
    with its index, so we know the expected covariance exactly."""
    rng = np.random.default_rng(42)
    base = rng.random((rows, cols))
    arr = np.stack([base * (b + 1) for b in range(n_bands)], axis=0).astype(np.float64)
    return arr


class TestNodata:

    def test_nodata_is_minus_9999(self):
        assert NODATA == -9999.0

    def test_nodata_type_float(self):
        assert isinstance(NODATA, float)


@skip_no_lib
class TestBandCovarianceShape:

    def test_square_output_shape(self):
        arr = _cube(4, 16, 16)
        out = band_covariance(arr)
        assert out.shape == (4, 4)

    def test_shape_single_band(self):
        arr = _cube(1, 8, 8)
        out = band_covariance(arr)
        assert out.shape == (1, 1)

    def test_shape_many_bands(self):
        arr = _cube(12, 20, 20)
        out = band_covariance(arr)
        assert out.shape == (12, 12)

    def test_non_square_spatial_extent(self):
        arr = _cube(6, 10, 40)
        out = band_covariance(arr)
        assert out.shape == (6, 6)

    def test_output_dtype_float64(self):
        arr = _cube(4, 8, 8)
        out = band_covariance(arr)
        assert out.dtype == np.float64

    def test_input_not_mutated(self):
        arr = _cube(4, 8, 8)
        arr_copy = arr.copy()
        band_covariance(arr)
        np.testing.assert_array_equal(arr, arr_copy)


@skip_no_lib
class TestBandCovarianceSymmetry:

    def test_symmetric_clean_data(self):
        arr = _cube(5, 30, 30)
        out = band_covariance(arr)
        np.testing.assert_allclose(out, out.T, atol=1e-12)

    def test_symmetric_with_nodata(self):
        arr = _cube(4, 20, 20)
        arr[1, :10, :10] = NODATA
        out = band_covariance(arr)
        valid = out != NODATA
        # only compare entries that are both defined
        both = valid & valid.T
        np.testing.assert_allclose(out[both], out.T[both], atol=1e-12)

    def test_diagonal_non_negative(self):
        """Variance is always >= 0."""
        arr = _cube(5, 20, 20)
        out = band_covariance(arr)
        diag = np.diag(out)
        valid = diag != NODATA
        assert np.all(diag[valid] >= 0.0)

    def test_positive_semi_definite(self):
        """All eigenvalues of a clean covariance matrix must be >= 0."""
        arr = _cube(5, 50, 50)
        out = band_covariance(arr)
        eigvals = np.linalg.eigvalsh(out)
        assert np.all(eigvals >= -1e-10)


@skip_no_lib
class TestBandCovarianceCorrectness:

    def test_matches_numpy_cov_clean(self):
        """Full agreement with numpy.cov on clean data."""
        n_bands, rows, cols = 4, 50, 60
        arr = _cube(n_bands, rows, cols, seed=7)
        got = band_covariance(arr)
        want = np.cov(arr.reshape(n_bands, -1))
        np.testing.assert_allclose(got, want, atol=1e-10)

    def test_matches_numpy_cov_large(self):
        n_bands, rows, cols = 6, 100, 100
        arr = _cube(n_bands, rows, cols, seed=42)
        got = band_covariance(arr)
        want = np.cov(arr.reshape(n_bands, -1))
        np.testing.assert_allclose(got, want, atol=1e-9)

    def test_single_band_variance_matches_numpy(self):
        arr = _cube(1, 40, 40, seed=3)
        got = band_covariance(arr)
        want = np.var(arr, ddof=1)
        assert got[0, 0] == pytest.approx(want, rel=1e-9)

    def test_two_band_off_diagonal_matches_numpy(self):
        arr = _cube(2, 30, 30, seed=5)
        got = band_covariance(arr)
        flat = arr.reshape(2, -1)
        want = np.cov(flat)
        assert got[0, 1] == pytest.approx(want[0, 1], rel=1e-9)
        assert got[1, 0] == pytest.approx(want[1, 0], rel=1e-9)

    def test_perfectly_correlated_bands(self):
        """Band 1 = 2 * Band 0 + 1  →  Pearson r should be exactly 1."""
        arr = _cube(2, 40, 40, seed=9)
        arr[1] = arr[0] * 2.0 + 1.0
        out = band_covariance(arr)
        # corr(0,1) = cov(0,1) / sqrt(var0 * var1)
        corr = out[0, 1] / np.sqrt(out[0, 0] * out[1, 1])
        assert corr == pytest.approx(1.0, abs=1e-9)

    def test_constant_band_variance_zero(self):
        """A constant band has zero variance."""
        arr = _cube(3, 20, 20)
        arr[1, :, :] = 0.5
        out = band_covariance(arr)
        assert out[1, 1] == pytest.approx(0.0, abs=1e-10)

    def test_float32_input_accepted(self):
        """float32 input should be cast internally and give a close result."""
        arr64 = _cube(3, 20, 20, seed=11)
        arr32 = arr64.astype(np.float32)
        got64 = band_covariance(arr64)
        got32 = band_covariance(arr32)
        np.testing.assert_allclose(got32, got64, atol=1e-5)


@skip_no_lib
class TestBandCovarianceNodata:

    def test_all_nodata_all_entries_nodata(self):
        arr = _nodata_cube(4, 10, 10)
        out = band_covariance(arr)
        assert np.all(out == NODATA)

    def test_single_dead_band_row_col_nodata(self):
        """A fully NODATA band produces NODATA in its row and column."""
        arr = _cube(4, 20, 20)
        arr[2, :, :] = NODATA
        out = band_covariance(arr)
        assert out[2, 2] == NODATA
        for j in range(4):
            if j != 2:
                assert out[2, j] == NODATA
                assert out[j, 2] == NODATA

    def test_partial_nodata_excludes_affected_pixels(self):
        """Masking some pixels in band 0 should change its variance."""
        rng = np.random.default_rng(0)
        arr = _cube(3, 40, 40)
        arr_nd = arr.copy()
        mask = rng.random((40, 40)) < 0.3
        arr_nd[0][mask] = NODATA
        out_full = band_covariance(arr)
        out_nd = band_covariance(arr_nd)
        assert out_nd[0, 0] != pytest.approx(out_full[0, 0], rel=1e-3)

    def test_partial_nodata_symmetry_preserved(self):
        arr = _cube(4, 20, 20)
        arr[0, :5, :5] = NODATA
        arr[3, 10:, 10:] = NODATA
        out = band_covariance(arr)
        for i in range(4):
            for j in range(4):
                if out[i, j] != NODATA and out[j, i] != NODATA:
                    assert out[i, j] == pytest.approx(out[j, i], abs=1e-12)

    def test_two_valid_obs_gives_defined_value(self):
        """Exactly 2 jointly valid pixels is the minimum for sample covariance."""
        arr = np.full((2, 1, 2), NODATA, dtype=np.float64)
        arr[0, 0, 0] = 0.2
        arr[1, 0, 0] = 0.8
        arr[0, 0, 1] = 0.4
        arr[1, 0, 1] = 0.6
        out = band_covariance(arr)
        assert out[0, 0] != NODATA
        assert out[1, 1] != NODATA
        assert out[0, 1] != NODATA

    def test_one_joint_obs_returns_nodata(self):
        """Fewer than 2 jointly valid pixels → NODATA."""
        arr = np.full((2, 1, 2), NODATA, dtype=np.float64)
        arr[0, 0, 0] = 0.3
        arr[1, 0, 0] = 0.7
        arr[0, 0, 1] = 0.5
        out = band_covariance(arr)
        assert out[0, 1] == NODATA
        assert out[1, 0] == NODATA

    def test_nodata_not_treated_as_valid_value(self):
        """Replacing all NODATA with a real value must change the result."""
        arr_nd = _nodata_cube(2, 10, 10)
        arr_nd[0, 0, 0] = 0.5
        arr_nd[1, 0, 0] = 0.5
        out_nd = band_covariance(arr_nd)

        arr_real = np.full((2, 10, 10), 0.5, dtype=np.float64)
        out_real = band_covariance(arr_real)
        assert out_nd[0, 0] == NODATA
        assert out_real[0, 0] == pytest.approx(0.0, abs=1e-10)


@skip_no_lib
class TestBandCovarianceKahan:

    def test_large_offset_does_not_break_variance(self):
        """Values clustered near a large constant — naive summation loses
        precision but Kahan should stay accurate."""
        rng = np.random.default_rng(0)
        small = rng.random((3, 80, 80)) * 1e-4  # tiny spread
        arr = small + 1e6  # large offset
        out = band_covariance(arr)
        want = np.cov(arr.reshape(3, -1))
        np.testing.assert_allclose(out, want, rtol=1e-6)

    def test_covariance_stable_for_high_pixel_count(self):
        """500×500 pixel cube — result must still match numpy.cov closely."""
        arr = _cube(4, 500, 500, seed=2)
        got = band_covariance(arr)
        want = np.cov(arr.reshape(4, -1))
        np.testing.assert_allclose(got, want, atol=1e-8)


@skip_no_lib
class TestBandCovarianceValidation:

    def test_2d_input_raises_value_error(self):
        with pytest.raises(ValueError, match="3-D"):
            band_covariance(np.ones((8, 8)))

    def test_1d_input_raises_value_error(self):
        with pytest.raises(ValueError, match="3-D"):
            band_covariance(np.ones(8))

    def test_4d_input_raises_value_error(self):
        with pytest.raises(ValueError, match="3-D"):
            band_covariance(np.ones((2, 4, 4, 4)))

    def test_error_message_includes_actual_shape(self):
        with pytest.raises(ValueError, match=r"\(8, 8\)"):
            band_covariance(np.ones((8, 8)))


@skip_no_lib
class TestBandCovariancePCA:

    def test_eigh_accepts_output(self):
        """np.linalg.eigh must not raise on the returned matrix."""
        arr = _cube(6, 40, 40)
        out = band_covariance(arr)
        eigvals, eigvecs = np.linalg.eigh(out)
        assert eigvecs.shape == (6, 6)

    def test_explained_variance_sums_to_one(self):
        arr = _cube(5, 50, 50)
        out = band_covariance(arr)
        eigvals = np.linalg.eigvalsh(out)
        eigvals = np.maximum(eigvals, 0.0)  # clamp numerical noise
        explained = eigvals / eigvals.sum()
        assert explained.sum() == pytest.approx(1.0, abs=1e-9)

    def test_first_pc_captures_most_variance(self):
        """PC1 should explain more variance than PC2."""
        arr = _cube(5, 50, 50)
        out = band_covariance(arr)
        eigvals = np.sort(np.linalg.eigvalsh(out))[::-1]
        assert eigvals[0] >= eigvals[1]

    def test_reconstruction_error_small(self):
        """Projecting to all PCs and back should recover original data."""
        n_bands, rows, cols = 4, 20, 20
        arr = _cube(n_bands, rows, cols)
        out = band_covariance(arr)
        _, eigvecs = np.linalg.eigh(out)
        flat = arr.reshape(n_bands, -1)
        projected = eigvecs.T @ flat
        reconstructed = eigvecs @ projected
        np.testing.assert_allclose(reconstructed, flat, atol=1e-10)
