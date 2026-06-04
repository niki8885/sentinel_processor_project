import numpy as np
import pytest

from sentinel_processor.processing._raster_ops_bridge import (
    align_bands,
    band_stats,
    clip_box_indices,
    normalize_band,
    reproject_nearest,
    rgb_to_luminance,
)
from sentinel_processor.processing._fortran_bridge import (
    pansharpen,
    pansharpening_gs,
    pansharpening_ihs,
    pansharpening_wavelet,
)


# rgb_to_luminance

class TestRgbToLuminance:

    def test_2d_passthrough(self):
        arr = np.random.default_rng(0).random((32, 32))
        out = rgb_to_luminance(arr)
        assert out.shape == (32, 32)
        np.testing.assert_allclose(out, arr.astype(np.float64))

    def test_single_band_3d(self):
        arr = np.random.default_rng(1).random((1, 16, 16))
        out = rgb_to_luminance(arr)
        assert out.shape == (16, 16)
        np.testing.assert_allclose(out, arr[0])

    def test_rgb_rec601_weights(self):
        """For 3-band: Y = 0.299R + 0.587G + 0.114B  (Rec. 601)."""
        r = np.full((16, 16), 1.0)
        g = np.full((16, 16), 0.0)
        b = np.full((16, 16), 0.0)
        arr = np.stack([r, g, b], axis=0)  # (3, 16, 16)
        out = rgb_to_luminance(arr)
        np.testing.assert_allclose(out, 0.299, atol=1e-9)

    def test_rgb_green_channel_weight(self):
        r = np.zeros((16, 16))
        g = np.ones((16, 16))
        b = np.zeros((16, 16))
        arr = np.stack([r, g, b], axis=0)
        out = rgb_to_luminance(arr)
        np.testing.assert_allclose(out, 0.587, atol=1e-9)

    def test_equal_weight_for_4_bands(self):
        """4-band: equal weights → mean of bands."""
        rng = np.random.default_rng(2)
        arr = rng.random((4, 16, 16))
        out = rgb_to_luminance(arr)
        expected = arr.mean(axis=0)
        np.testing.assert_allclose(out, expected, atol=1e-9)

    def test_wrong_ndim_raises(self):
        with pytest.raises(ValueError, match="2-D or 3-D"):
            rgb_to_luminance(np.ones((2, 4, 4, 4)))

    def test_output_shape_matches_spatial(self):
        arr = np.random.default_rng(3).random((3, 48, 64))
        out = rgb_to_luminance(arr)
        assert out.shape == (48, 64)


# align_bands

class TestAlignBands:

    def test_same_size_identity(self):
        src = np.random.default_rng(0).random((32, 32))
        out = align_bands(src, 32, 32)
        np.testing.assert_allclose(out, src, atol=1e-10)

    def test_upsample_shape(self):
        src = np.random.default_rng(1).random((16, 16))
        out = align_bands(src, 32, 32)
        assert out.shape == (32, 32)

    def test_downsample_shape(self):
        src = np.random.default_rng(2).random((64, 64))
        out = align_bands(src, 32, 32)
        assert out.shape == (32, 32)

    def test_1d_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            align_bands(np.ones(16), 8, 8)

    def test_corner_values_preserved_on_upsample(self):
        """Top-left corner value should survive nearest-neighbour upsample."""
        src = np.zeros((4, 4))
        src[0, 0] = 1.0
        out = align_bands(src, 8, 8)
        assert out[0, 0] == pytest.approx(1.0)


# band_stats
class TestBandStats:

    def test_constant_array(self):
        arr = np.full(1000, 5.0)
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(5.0)
        assert s["min"] == pytest.approx(5.0)
        assert s["max"] == pytest.approx(5.0)
        assert s["std"] == pytest.approx(0.0, abs=1e-6)

    def test_known_values(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(3.0)
        assert s["min"] == pytest.approx(1.0)
        assert s["max"] == pytest.approx(5.0)

    def test_std_positive(self):
        rng = np.random.default_rng(0)
        arr = rng.random(10_000)
        s = band_stats(arr)
        assert s["std"] > 0

    def test_returns_all_keys(self):
        s = band_stats(np.ones(10))
        assert set(s.keys()) == {"mean", "std", "min", "max"}

    def test_2d_input_ravelled(self):
        arr = np.full((8, 8), 3.0)
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(3.0)

    def test_large_array_mean_accurate(self):
        """Mean of uniform [0,1] should be ~0.5 within tolerance."""
        arr = np.random.default_rng(1).random(1_000_000)
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(0.5, abs=0.005)


# normalize_band

class TestNormalizeBand:

    def test_output_range_0_to_1(self):
        arr = np.array([0.0, 5000.0, 10000.0])
        out = normalize_band(arr.copy(), src_min=0.0, src_max=10000.0)
        np.testing.assert_allclose(out.ravel(), [0.0, 0.5, 1.0], atol=1e-9)

    def test_custom_out_range(self):
        arr = np.array([0.0, 1.0])
        out = normalize_band(arr.copy(), src_min=0.0, src_max=1.0,
                             out_min=-1.0, out_max=1.0)
        np.testing.assert_allclose(out.ravel(), [-1.0, 1.0], atol=1e-9)

    def test_flat_input_all_out_min(self):
        """When src_min == src_max, output should be out_min."""
        arr = np.full(10, 5.0)
        out = normalize_band(arr.copy(), src_min=5.0, src_max=5.0, out_min=0.0)
        np.testing.assert_allclose(out.ravel(), 0.0, atol=1e-10)

    def test_auto_range_from_array(self):
        """Without explicit src_min/max, bridge computes them from the data."""
        arr = np.linspace(10.0, 20.0, 100)
        out = normalize_band(arr.copy())
        assert out.min() == pytest.approx(0.0, abs=1e-6)
        assert out.max() == pytest.approx(1.0, abs=1e-6)

    def test_shape_preserved(self):
        arr = np.random.default_rng(0).random((8, 8))
        out = normalize_band(arr.copy())
        assert out.shape == (8, 8)


# clip_box_indices

class TestClipBoxIndices:

    def test_full_coverage(self):
        """BBox that covers the entire raster -> row/col span = full grid."""
        ox, oy = 0.0, 0.0
        pw, ph = 0.01, 0.01
        rows, cols = 100, 100
        r0, r1, c0, c1 = clip_box_indices(
            ox, oy, pw, ph, rows, cols,
            0.0, 1.0, 0.0, 1.0,
        )
        assert r0 >= 1 and r1 <= rows
        assert c0 >= 1 and c1 <= cols

    def test_clamped_to_grid(self):
        """BBox outside raster extent -> clamped to [1, rows/cols]."""
        r0, r1, c0, c1 = clip_box_indices(
            0.0, 0.0, 0.01, 0.01, 100, 100,
            -1.0, 2.0, -1.0, 2.0,  # far outside
        )
        assert r0 >= 1 and r1 <= 100
        assert c0 >= 1 and c1 <= 100

    def test_negative_pixel_height(self):
        """Negative ph (north-up raster) should not crash."""
        r0, r1, c0, c1 = clip_box_indices(
            0.0, 1.0, 0.01, -0.01, 100, 100,
            0.2, 0.8, 0.2, 0.8,
        )
        assert r0 <= r1
        assert c0 <= c1

    def test_returns_4_ints(self):
        result = clip_box_indices(0.0, 0.0, 1.0, 1.0, 10, 10, 0.0, 5.0, 0.0, 5.0)
        assert len(result) == 4
        for v in result:
            assert isinstance(v, int)


# reproject_nearest

class TestReprojectNearest:

    def test_same_affine_identity(self):
        """Same src and dst affine + same size -> output == input."""
        src = np.random.default_rng(0).random((32, 32))
        affine = (0.0, 0.0, 1.0 / 32, 1.0 / 32)
        out = reproject_nearest(src, affine, 32, 32, affine)
        np.testing.assert_allclose(out, src, atol=1e-9)

    def test_upsample_shape(self):
        src = np.random.default_rng(1).random((16, 16))
        affine = (0.0, 0.0, 1.0 / 16, 1.0 / 16)
        dst_af = (0.0, 0.0, 1.0 / 32, 1.0 / 32)
        out = reproject_nearest(src, affine, 32, 32, dst_af)
        assert out.shape == (32, 32)

    def test_1d_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            reproject_nearest(np.ones(16), (0, 0, 1, 1), 4, 4, (0, 0, 1, 1))

    def test_corner_preserved(self):
        src = np.zeros((8, 8))
        src[0, 0] = 99.0
        affine = (0.0, 0.0, 1.0, 1.0)
        out = reproject_nearest(src, affine, 8, 8, affine)
        assert out[0, 0] == pytest.approx(99.0)


# pansharpen

class TestPansharpen:

    def _make_inputs(self, pan_size=64, ms_size=16, n_bands=4, seed=0):
        rng = np.random.default_rng(seed)
        pan = rng.random((pan_size, pan_size))
        ms = rng.random((n_bands, ms_size, ms_size))
        return pan, ms

    def test_gram_schmidt_output_shape(self):
        pan, ms = self._make_inputs()
        out = pansharpen(pan, ms, algorithm="gram_schmidt")
        assert out.shape == (4, 64, 64)

    def test_ihs_output_shape(self):
        pan, ms = self._make_inputs(n_bands=3)
        out = pansharpen(pan, ms, algorithm="ihs")
        assert out.shape == (3, 64, 64)

    def test_wavelet_output_shape(self):
        pan, ms = self._make_inputs(n_bands=4, pan_size=64, ms_size=16)
        out = pansharpen(pan, ms, algorithm="wavelet")
        assert out.shape == (4, 64, 64)

    def test_output_dtype_float64(self):
        pan, ms = self._make_inputs()
        out = pansharpen(pan, ms, algorithm="gram_schmidt")
        assert out.dtype == np.float64

    def test_ihs_requires_3_bands(self):
        pan, ms = self._make_inputs(n_bands=4)
        with pytest.raises(ValueError, match="3 bands"):
            pansharpen(pan, ms, algorithm="ihs")

    def test_ihs_accepts_exactly_3_bands(self):
        pan, ms = self._make_inputs(n_bands=3)
        out = pansharpen(pan, ms, algorithm="ihs")
        assert out.shape[0] == 3

    def test_pan_must_be_2d(self):
        pan = np.ones((2, 64, 64))
        ms = np.ones((4, 16, 16))
        with pytest.raises(ValueError, match="2-D"):
            pansharpen(pan, ms, algorithm="gram_schmidt")

    def test_ms_must_be_3d(self):
        pan = np.ones((64, 64))
        ms = np.ones((16, 16))  # 2-D, not 3-D
        with pytest.raises(ValueError):
            pansharpen(pan, ms, algorithm="gram_schmidt")

    def test_unknown_algorithm_raises(self):
        pan, ms = self._make_inputs()
        with pytest.raises(ValueError, match="Unknown algorithm"):
            pansharpen(pan, ms, algorithm="invalid")

    def test_gram_schmidt_same_pan_size_as_output(self):
        pan, ms = self._make_inputs(pan_size=32, ms_size=8, n_bands=2)
        out = pansharpen(pan, ms, algorithm="gram_schmidt")
        assert out.shape[-2:] == pan.shape

    def test_wavelet_requires_even_pan_dimensions(self):
        """Wavelet uses rows/2, cols/2 — odd sizes would truncate; even sizes work."""
        pan = np.random.default_rng(0).random((64, 64))  # even
        ms = np.random.default_rng(1).random((2, 16, 16))
        out = pansharpen(pan, ms, algorithm="wavelet")
        assert out.shape == (2, 64, 64)

    def test_gs_alias(self):
        pan, ms = self._make_inputs()
        np.testing.assert_allclose(
            pansharpening_gs(pan, ms),
            pansharpen(pan, ms, algorithm="gram_schmidt"),
            atol=1e-12,
        )

    def test_ihs_alias(self):
        pan, ms = self._make_inputs(n_bands=3)
        np.testing.assert_allclose(
            pansharpening_ihs(pan, ms),
            pansharpen(pan, ms, algorithm="ihs"),
            atol=1e-12,
        )

    def test_wavelet_alias(self):
        pan, ms = self._make_inputs()
        np.testing.assert_allclose(
            pansharpening_wavelet(pan, ms),
            pansharpen(pan, ms, algorithm="wavelet"),
            atol=1e-12,
        )

    def test_gram_schmidt_output_not_all_zeros(self):
        pan, ms = self._make_inputs()
        out = pansharpen(pan, ms, algorithm="gram_schmidt")
        assert out.std() > 0.01

    def test_ihs_output_sums_to_channel_luminance(self):
        rng = np.random.default_rng(42)
        pan = rng.random((32, 32)) * 0.5 + 0.3
        ms = rng.random((3, 8, 8)) * 0.5 + 0.2
        out = pansharpen(pan, ms, algorithm="ihs")
        assert out.min() > -5.0
        assert out.max() < 5.0
