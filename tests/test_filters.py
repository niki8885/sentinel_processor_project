import math
import numpy as np
import pytest

from sentinel_processor.filters._filters_bridge import (
    _to_f_f64,
    _empty_f,
    _prepare_2d,
    apply_to_bands,
    bilateral_filter,
    convolve2d,
    gaussian_blur,
    laplacian,
    median_filter,
    morpho_close,
    morpho_dilate,
    morpho_erode,
    morpho_open,
    sobel_direction,
    sobel_magnitude,
    top_hat_black,
    top_hat_white,
    unsharp_mask,
)
from sentinel_processor.filters.compute import (
    AVAILABLE_FILTERS,
    apply_filter,
    list_filters,
)


def _ramp(rows=16, cols=16):
    """Simple gradient ramp — predictable values for edge/blur tests."""
    r = np.linspace(0.0, 1.0, rows * cols).reshape(rows, cols)
    return r


def _flat(val=0.5, rows=16, cols=16):
    return np.full((rows, cols), val, dtype=np.float64)


def _rand(rows=32, cols=32, seed=0):
    return np.random.default_rng(seed).random((rows, cols))


class TestToFF64:

    def test_c_order_converted_to_f_order(self):
        a = np.ones((4, 4), dtype=np.float64, order="C")
        _, out = _to_f_f64(a)
        assert out.flags["F_CONTIGUOUS"]

    def test_already_f_order_no_copy(self):
        a = np.asfortranarray(np.ones((4, 4), dtype=np.float64))
        _, out = _to_f_f64(a)
        assert out.ctypes.data == a.ctypes.data

    def test_float32_converted(self):
        a = np.ones((4, 4), dtype=np.float32)
        _, out = _to_f_f64(a)
        assert out.dtype == np.float64

    def test_returns_pointer_and_array(self):
        ptr, arr = _to_f_f64(np.ones((2, 2)))
        assert arr.ndim == 2


class TestPrepare2d:

    def test_2d_accepted(self):
        a, r, c = _prepare_2d(np.ones((8, 8)))
        assert r == 8 and c == 8

    def test_1d_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            _prepare_2d(np.ones(8))

    def test_3d_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            _prepare_2d(np.ones((2, 8, 8)))


# Gaussian blur

class TestGaussianBlur:

    def test_flat_image_unchanged(self):
        src = _flat(0.5)
        out = gaussian_blur(src, sigma=2.0)
        np.testing.assert_allclose(out, src, atol=1e-10)

    def test_output_shape(self):
        src = _rand(64, 64)
        out = gaussian_blur(src, sigma=1.5)
        assert out.shape == src.shape

    def test_smooths_noise(self):
        rng = np.random.default_rng(0)
        src = rng.random((64, 64))
        out = gaussian_blur(src, sigma=3.0)
        assert out.std() < src.std()

    def test_larger_sigma_more_smooth(self):
        src = _rand(64, 64, seed=1)
        s1 = gaussian_blur(src, sigma=1.0).std()
        s3 = gaussian_blur(src, sigma=3.0).std()
        assert s3 < s1

    def test_explicit_kradius(self):
        src = _rand()
        out = gaussian_blur(src, sigma=1.0, kradius=5)
        assert out.shape == src.shape

    def test_default_kradius_equals_ceil3sigma(self):
        src = _rand(32, 32, seed=2)
        sigma = 1.5
        kr = math.ceil(3 * sigma)
        auto = gaussian_blur(src, sigma=sigma)
        manual = gaussian_blur(src, sigma=sigma, kradius=kr)
        np.testing.assert_allclose(auto, manual, atol=1e-12)


# Sobel magnitude

class TestSobelMagnitude:

    def test_flat_image_zero_edges(self):
        src = _flat(0.5, 32, 32)
        out = sobel_magnitude(src, norm="l2")
        interior = out[1:-1, 1:-1]
        np.testing.assert_allclose(interior, 0.0, atol=1e-10)

    def test_output_non_negative(self):
        src = _rand(32, 32)
        assert sobel_magnitude(src, norm="l2").min() >= 0.0
        assert sobel_magnitude(src, norm="l1").min() >= 0.0

    def test_l1_ge_l2(self):
        src = _rand(32, 32, seed=3)
        l1 = sobel_magnitude(src, norm="l1")
        l2 = sobel_magnitude(src, norm="l2")
        assert np.all(l1 >= l2 - 1e-10)

    def test_output_shape(self):
        src = _rand(48, 64)
        assert sobel_magnitude(src).shape == (48, 64)

    def test_step_edge_detected(self):
        src = np.zeros((32, 32))
        src[:, 16:] = 1.0
        out = sobel_magnitude(src)
        assert out[:, 16].mean() > out[:, 4].mean()


# Sobel direction

class TestSobelDirection:

    def test_output_range(self):
        """All values must be in [-π, π]."""
        src = _rand(32, 32)
        out = sobel_direction(src)
        assert out.min() >= -math.pi - 1e-6
        assert out.max() <= math.pi + 1e-6

    def test_output_shape(self):
        src = _rand(32, 48)
        assert sobel_direction(src).shape == (32, 48)

    def test_horizontal_edge_direction(self):
        """Horizontal step edge -> gradient should point vertically (≈ ±π/2)."""
        src = np.zeros((32, 32))
        src[16:, :] = 1.0
        out = sobel_direction(src)
        row_angles = out[15, 2:-2]
        assert np.abs(row_angles).mean() > 0.5


# Laplacian

class TestLaplacian:

    def test_flat_image_zero(self):
        src = _flat(0.5, 32, 32)
        out = laplacian(src, connectivity=4)
        np.testing.assert_allclose(out[1:-1, 1:-1], 0.0, atol=1e-10)

    def test_connectivity_4_vs_8_differ_on_diagonal(self):
        """On a diagonal ramp, 4- and 8-conn give different results."""
        src = _ramp(32, 32)
        out4 = laplacian(src, connectivity=4)
        out8 = laplacian(src, connectivity=8)
        assert not np.allclose(out4, out8)

    def test_invalid_connectivity_raises(self):
        with pytest.raises(ValueError, match="connectivity"):
            laplacian(_rand(), connectivity=6)

    def test_output_shape(self):
        src = _rand(24, 36)
        assert laplacian(src, connectivity=4).shape == (24, 36)


# Unsharp mask

class TestUnsharpMask:

    def test_flat_image_unchanged(self):
        """On flat input, residual = 0 -> output == input."""
        src = _flat(0.5)
        out = unsharp_mask(src, sigma=1.0, amount=2.0)
        np.testing.assert_allclose(out, src, atol=1e-10)

    def test_sharpening_increases_local_contrast(self):
        """After unsharp mask, local std should be >= original."""
        src = _rand(32, 32, seed=4)
        out = unsharp_mask(src, sigma=1.0, amount=2.0)
        assert out.std() >= src.std() - 1e-10

    def test_threshold_suppresses_noise(self):
        """High threshold: near-flat areas should be unchanged."""
        src = _flat(0.5) + np.random.default_rng(5).random((16, 16)) * 0.001
        out = unsharp_mask(src, sigma=1.0, amount=10.0, threshold=0.1)
        np.testing.assert_allclose(out, src, atol=0.01)

    def test_output_shape(self):
        assert unsharp_mask(_rand()).shape == _rand().shape


# Median filter

class TestMedianFilter:

    def test_flat_image_unchanged(self):
        src = _flat(0.7)
        out = median_filter(src, radius=2)
        np.testing.assert_allclose(out, src, atol=1e-10)

    def test_removes_salt_pepper(self):
        """Single hot pixel should be removed by median r=1."""
        src = _flat(0.5, 16, 16)
        src[8, 8] = 1.0  # salt pixel
        out = median_filter(src, radius=1)
        assert out[8, 8] < 0.6  # should be suppressed

    def test_radius_1_vs_3(self):
        """Larger radius smooths more — interior of output should be flatter."""
        src = _rand(32, 32, seed=6)
        out1 = median_filter(src, radius=1)
        out3 = median_filter(src, radius=3)
        assert out3[4:-4, 4:-4].std() <= out1[4:-4, 4:-4].std() + 1e-10

    def test_output_shape(self):
        assert median_filter(_rand(24, 36), radius=2).shape == (24, 36)


# Bilateral filter

class TestBilateralFilter:

    def test_flat_image_unchanged(self):
        src = _flat(0.5, 16, 16)
        out = bilateral_filter(src, sigma_s=2.0, sigma_r=0.1)
        np.testing.assert_allclose(out, src, atol=1e-8)

    def test_output_shape(self):
        assert bilateral_filter(_rand()).shape == _rand().shape

    def test_edge_preserving(self):
        """Edge pixels should retain larger values than interior smoothed pixels."""
        src = np.zeros((32, 32))
        src[:, 16:] = 1.0  # hard step edge
        out = bilateral_filter(src, sigma_s=3.0, sigma_r=0.05)

        left_mean = out[:, 4].mean()
        right_mean = out[:, 28].mean()
        assert right_mean - left_mean > 0.5


# Morphological filters

class TestMorphoErode:

    def test_flat_unchanged(self):
        src = _flat(0.5)
        out = morpho_erode(src, radius=2)
        np.testing.assert_allclose(out, src, atol=1e-10)

    def test_erode_reduces_max(self):
        src = _rand(32, 32, seed=7)
        out = morpho_erode(src, radius=2)
        assert out.max() <= src.max() + 1e-10

    def test_erode_reduces_mean(self):
        src = _rand(32, 32, seed=8)
        out = morpho_erode(src, radius=1)
        assert out.mean() <= src.mean() + 1e-10

    def test_output_shape(self):
        assert morpho_erode(_rand()).shape == _rand().shape


class TestMorphoDilate:

    def test_flat_unchanged(self):
        src = _flat(0.5)
        out = morpho_dilate(src, radius=2)
        np.testing.assert_allclose(out, src, atol=1e-10)

    def test_dilate_increases_min(self):
        src = _rand(32, 32, seed=9)
        out = morpho_dilate(src, radius=2)
        assert out.min() >= src.min() - 1e-10

    def test_dilate_duality_with_erode(self):
        """Duality: dilate(src) = 1 - erode(1 - src) for binary image."""
        src = (_rand(32, 32, seed=10) > 0.5).astype(np.float64)
        dil = morpho_dilate(src, radius=2)
        er = morpho_erode(1.0 - src, radius=2)
        np.testing.assert_allclose(dil, 1.0 - er, atol=1e-10)


class TestMorphoOpenClose:

    def test_open_removes_small_bright_blob(self):
        src = _flat(0.0, 32, 32)
        src[15, 15] = 1.0  # single bright pixel
        out = morpho_open(src, radius=2)
        assert out[15, 15] < 0.1  # removed

    def test_close_fills_small_dark_hole(self):
        src = _flat(1.0, 32, 32)
        src[15, 15] = 0.0  # single dark hole
        out = morpho_close(src, radius=2)
        assert out[15, 15] > 0.9  # filled

    def test_open_is_erode_then_dilate(self):
        src = _rand(32, 32, seed=11)
        manual = morpho_dilate(morpho_erode(src, 1), 1)
        fn_out = morpho_open(src, 1)
        np.testing.assert_allclose(fn_out, manual, atol=1e-10)

    def test_close_is_dilate_then_erode(self):
        src = _rand(32, 32, seed=12)
        manual = morpho_erode(morpho_dilate(src, 1), 1)
        fn_out = morpho_close(src, 1)
        np.testing.assert_allclose(fn_out, manual, atol=1e-10)


# Top-hat filters

class TestTopHat:

    def test_white_tophat_nonnegative(self):
        """White top-hat is always >= 0."""
        src = _rand(32, 32, seed=13)
        out = top_hat_white(src, radius=3)
        assert out.min() >= -1e-10

    def test_black_tophat_nonnegative(self):
        src = _rand(32, 32, seed=14)
        out = top_hat_black(src, radius=3)
        assert out.min() >= -1e-10

    def test_white_tophat_flat_is_zero(self):
        src = _flat(0.5)
        out = top_hat_white(src, radius=3)
        np.testing.assert_allclose(out, 0.0, atol=1e-10)

    def test_black_tophat_flat_is_zero(self):
        src = _flat(0.5)
        out = top_hat_black(src, radius=3)
        np.testing.assert_allclose(out, 0.0, atol=1e-10)

    def test_white_tophat_detects_small_bright(self):
        """Bright pixel smaller than SE radius should appear in white top-hat."""
        src = _flat(0.0, 32, 32)
        src[16, 16] = 1.0
        out = top_hat_white(src, radius=5)
        assert out[16, 16] > 0.5

    def test_white_tophat_is_src_minus_opening(self):
        src = _rand(32, 32, seed=15)
        manual = src - morpho_dilate(morpho_erode(src, 3), 3)
        out = top_hat_white(src, radius=3)
        np.testing.assert_allclose(out, manual, atol=1e-10)

    def test_black_tophat_is_closing_minus_src(self):
        src = _rand(32, 32, seed=16)
        manual = morpho_erode(morpho_dilate(src, 3), 3) - src
        out = top_hat_black(src, radius=3)
        np.testing.assert_allclose(out, manual, atol=1e-10)


# Generic convolve2d

class TestConvolve2d:

    def test_identity_kernel(self):
        """3×3 kernel with 1 in center -> output equals input."""
        k = np.zeros((3, 3))
        k[1, 1] = 1.0
        src = _rand(16, 16)
        out = convolve2d(src, k)
        np.testing.assert_allclose(out[1:-1, 1:-1], src[1:-1, 1:-1], atol=1e-10)

    def test_mean_kernel_3x3(self):
        """3×3 box mean of flat image = same value."""
        k = np.ones((3, 3)) / 9.0
        src = _flat(0.6)
        out = convolve2d(src, k)
        np.testing.assert_allclose(out[1:-1, 1:-1], 0.6, atol=1e-10)

    def test_kernel_must_be_2d(self):
        with pytest.raises(ValueError, match="2-D"):
            convolve2d(_rand(), np.ones(9))

    def test_output_shape(self):
        src = _rand(24, 36)
        k = np.ones((3, 5)) / 15.0
        assert convolve2d(src, k).shape == (24, 36)

    def test_zero_kernel_gives_zeros(self):
        src = _rand(16, 16)
        k = np.zeros((3, 3))
        np.testing.assert_allclose(convolve2d(src, k), 0.0, atol=1e-10)

    def test_1x1_kernel_passthrough(self):
        src = _rand(16, 16)
        k = np.array([[1.0]])
        np.testing.assert_allclose(convolve2d(src, k), src, atol=1e-10)


# apply_to_bands  (multi-band wrapper)

class TestApplyToBands:

    def test_2d_passthrough(self):
        src = _rand()
        out = apply_to_bands(gaussian_blur, src, sigma=1.0)
        assert out.shape == src.shape

    def test_3d_applied_per_band(self):
        ms = np.stack([_rand(16, 16, seed=i) for i in range(4)], axis=0)
        assert ms.shape == (4, 16, 16)
        out = apply_to_bands(gaussian_blur, ms, sigma=1.0)
        assert out.shape == (4, 16, 16)

    def test_each_band_independent(self):
        """apply_to_bands result == stacking individual band calls."""
        ms = np.stack([_rand(16, 16, seed=i) for i in range(3)], axis=0)
        expected = np.stack([gaussian_blur(ms[b], sigma=1.5) for b in range(3)])
        result = apply_to_bands(gaussian_blur, ms, sigma=1.5)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_wrong_ndim_raises(self):
        with pytest.raises(ValueError):
            apply_to_bands(gaussian_blur, np.ones((2, 4, 4, 4)), sigma=1.0)


# apply_filter  (high-level registry)

class TestApplyFilter:

    def test_unknown_filter_raises(self):
        with pytest.raises(ValueError, match="Unknown filter"):
            apply_filter(_rand(), "nonexistent_filter")

    def test_all_registered_filters_callable(self):
        """Every filter in the registry must run without error on a small array."""
        src = _rand(16, 16)
        for name in AVAILABLE_FILTERS:
            if name == "convolve":
                result = apply_filter(src, name, kernel=np.ones((3, 3)) / 9)
            else:
                result = apply_filter(src, name)
            assert result.shape == src.shape, f"{name} changed shape"

    def test_3d_input_accepted(self):
        ms = np.stack([_rand(16, 16, seed=i) for i in range(3)], axis=0)
        out = apply_filter(ms, "gaussian", sigma=1.0)
        assert out.shape == ms.shape

    def test_1d_input_raises(self):
        with pytest.raises(ValueError):
            apply_filter(np.ones(100), "gaussian")

    def test_case_insensitive(self):
        src = _rand()
        out1 = apply_filter(src, "Gaussian", sigma=1.0)
        out2 = apply_filter(src, "gaussian", sigma=1.0)
        np.testing.assert_allclose(out1, out2)


class TestListFilters:

    def test_returns_dict(self):
        assert isinstance(list_filters(), dict)

    def test_all_filters_listed(self):
        assert set(list_filters().keys()) == set(AVAILABLE_FILTERS)

    def test_each_has_description(self):
        for name, entry in list_filters().items():
            assert "description" in entry
            assert len(entry["description"]) > 0

    def test_each_has_default_params(self):
        for name, entry in list_filters().items():
            assert "default_params" in entry
            assert isinstance(entry["default_params"], dict)
