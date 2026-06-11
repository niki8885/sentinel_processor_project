from __future__ import annotations

import numpy as np
import pytest

from sentinel_processor.wavelet._wavelet_bridge import (
    dwt2d,
    idwt2d,
    threshold_coeffs,
    estimate_sigma,
    bayes_threshold,
    bayes_denoise,
    band_energy,
    band_stats,
    dwt2d_batch,
    idwt2d_batch,
    dwt3d,
    idwt3d,
)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

ALL_WAVELETS = ["haar", "db4", "db6", "coif1", "sym4", "sym6"]

# sym6 needs at least 12 pixels per axis; use 48 as a safe multiple of
# 2^3=8 and >=12, so all wavelets and levels 1-3 work without padding.
ROWS, COLS = 48, 48


def _arr(rows=ROWS, cols=COLS, seed=0) -> np.ndarray:
    """Return a random (rows, cols) float64 array."""
    return np.random.default_rng(seed).random((rows, cols))


def _coeffs(arr=None, levels=1, wavelet="haar") -> dict:
    if arr is None:
        arr = _arr()
    return dwt2d(arr, levels=levels, wavelet=wavelet)


# ---------------------------------------------------------------------------
# dwt2d — output structure
# ---------------------------------------------------------------------------

class TestDwt2dStructure:

    def test_returns_dict(self):
        c = _coeffs()
        assert isinstance(c, dict)

    def test_level_keys_are_integers(self):
        c = _coeffs(levels=3)
        assert all(isinstance(k, int) for k in c)

    def test_level_keys_range(self):
        levels = 3
        c = _coeffs(levels=levels)
        assert set(c.keys()) == {1, 2, 3}

    def test_detail_bands_present_at_every_level(self):
        c = _coeffs(levels=3)
        for lv in c:
            assert "LH" in c[lv]
            assert "HL" in c[lv]
            assert "HH" in c[lv]

    def test_ll_only_at_coarsest_level(self):
        levels = 3
        c = _coeffs(levels=levels)
        assert "LL" in c[levels]
        for lv in range(1, levels):
            assert "LL" not in c[lv]

    def test_subband_shapes_halve_per_level(self):
        levels = 3
        c = _coeffs(levels=levels)
        for lv in range(1, levels + 1):
            expected = (ROWS >> lv, COLS >> lv)
            assert c[lv]["LH"].shape == expected
            assert c[lv]["HL"].shape == expected
            assert c[lv]["HH"].shape == expected

    def test_ll_shape_at_coarsest(self):
        levels = 3
        c = _coeffs(levels=levels)
        expected = (ROWS >> levels, COLS >> levels)
        assert c[levels]["LL"].shape == expected

    def test_subband_dtype_float64(self):
        c = _coeffs(levels=2)
        for lv in c.values():
            for arr in lv.values():
                assert arr.dtype == np.float64

    def test_subbands_are_c_contiguous(self):
        c = _coeffs(levels=2)
        for lv in c.values():
            for arr in lv.values():
                assert arr.flags["C_CONTIGUOUS"]

    def test_float32_input_accepted(self):
        arr = _arr().astype(np.float32)
        c = dwt2d(arr, levels=1, wavelet="haar")
        assert isinstance(c, dict)

    def test_integer_input_accepted(self):
        arr = (_arr() * 10000).astype(np.int16)
        c = dwt2d(arr, levels=1, wavelet="haar")
        assert isinstance(c, dict)


# ---------------------------------------------------------------------------
# dwt2d — input validation
# ---------------------------------------------------------------------------

class TestDwt2dValidation:

    def test_1d_input_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            dwt2d(np.ones(64), levels=1)

    def test_3d_input_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            dwt2d(np.ones((4, 16, 16)), levels=1)

    def test_rows_not_divisible_raises(self):
        # 48 rows, 2 levels → factor=4 → ok; use 46 rows → fail
        with pytest.raises(ValueError, match="rows"):
            dwt2d(np.ones((46, 48)), levels=2)

    def test_cols_not_divisible_raises(self):
        with pytest.raises(ValueError, match="cols"):
            dwt2d(np.ones((48, 46)), levels=2)

    def test_unknown_wavelet_raises(self):
        with pytest.raises(ValueError, match="Unknown wavelet"):
            dwt2d(_arr(), levels=1, wavelet="bior22")

    def test_unknown_wavelet_message_lists_valid(self):
        with pytest.raises(ValueError, match="haar"):
            dwt2d(_arr(), levels=1, wavelet="nope")


# ---------------------------------------------------------------------------
# dwt2d / idwt2d — perfect reconstruction
# ---------------------------------------------------------------------------

class TestPerfectReconstruction:

    @pytest.mark.parametrize("wavelet", ALL_WAVELETS)
    @pytest.mark.parametrize("levels", [1, 2, 3])
    def test_roundtrip(self, wavelet, levels):
        arr = _arr()
        coeffs = dwt2d(arr, levels=levels, wavelet=wavelet)
        recon  = idwt2d(coeffs, wavelet=wavelet)
        np.testing.assert_allclose(recon, arr, atol=1e-9)

    def test_reconstructed_shape_matches_input(self):
        arr = _arr()
        c = _coeffs(arr, levels=2)
        r = idwt2d(c)
        assert r.shape == arr.shape

    def test_reconstructed_dtype_float64(self):
        c = _coeffs()
        r = idwt2d(c)
        assert r.dtype == np.float64

    def test_reconstructed_c_contiguous(self):
        c = _coeffs()
        r = idwt2d(c)
        assert r.flags["C_CONTIGUOUS"]

    def test_modified_ll_changes_reconstruction(self):
        """Zeroing the LL sub-band must change the output."""
        arr = _arr()
        c = _coeffs(arr, levels=2)
        c[2]["LL"][:] = 0.0
        r = idwt2d(c)
        assert not np.allclose(r, arr)


# ---------------------------------------------------------------------------
# threshold_coeffs
# ---------------------------------------------------------------------------

class TestThresholdCoeffs:

    def test_returns_dict_same_keys(self):
        c = _coeffs(levels=2)
        t = threshold_coeffs(c, threshold=0.1)
        assert set(t.keys()) == set(c.keys())

    def test_detail_bands_thresholded(self):
        arr = _arr()
        c = _coeffs(arr, levels=2, wavelet="haar")
        # Use a large threshold so most coefficients are zeroed
        t = threshold_coeffs(c, threshold=1e6, mode="soft")
        for lv in t:
            for key in ("LH", "HL", "HH"):
                np.testing.assert_array_equal(t[lv][key], 0.0)

    def test_ll_never_modified(self):
        c = _coeffs(levels=2)
        original_ll = c[2]["LL"].copy()
        t = threshold_coeffs(c, threshold=1e6, mode="soft")
        np.testing.assert_array_equal(t[2]["LL"], original_ll)

    def test_soft_threshold_magnitude_reduced(self):
        c = _coeffs(levels=1)
        t = threshold_coeffs(c, threshold=0.05, mode="soft")
        for key in ("LH", "HL", "HH"):
            # All kept coefficients must be strictly smaller in magnitude
            orig = np.abs(c[1][key])
            thr  = np.abs(t[1][key])
            assert np.all(thr <= orig + 1e-12)

    def test_hard_threshold_zeros_small_coeffs(self):
        c = _coeffs(levels=1, wavelet="haar")
        thr = 0.3
        t = threshold_coeffs(c, threshold=thr, mode="hard")
        for key in ("LH", "HL", "HH"):
            vals = t[1][key]
            # every kept value must have |x| > thr
            nonzero = vals[vals != 0.0]
            assert np.all(np.abs(nonzero) > thr)

    def test_soft_sign_preserved(self):
        """Soft-thresholded nonzero coefficients keep their sign."""
        c = _coeffs(levels=1)
        t = threshold_coeffs(c, threshold=0.05, mode="soft")
        for key in ("LH", "HL", "HH"):
            orig = c[1][key]
            th   = t[1][key]
            nonzero = th != 0.0
            np.testing.assert_array_equal(
                np.sign(th[nonzero]), np.sign(orig[nonzero])
            )

    def test_input_not_modified(self):
        """threshold_coeffs must not mutate the input dict."""
        c = _coeffs(levels=1)
        original_hh = c[1]["HH"].copy()
        threshold_coeffs(c, threshold=0.5, mode="soft")
        np.testing.assert_array_equal(c[1]["HH"], original_hh)

    def test_negative_threshold_raises(self):
        c = _coeffs()
        with pytest.raises(ValueError, match="threshold"):
            threshold_coeffs(c, threshold=-0.1)

    def test_invalid_mode_raises(self):
        c = _coeffs()
        with pytest.raises(ValueError, match="mode"):
            threshold_coeffs(c, threshold=0.1, mode="median")

    def test_zero_threshold_is_identity(self):
        """Zero threshold with soft mode should not change coefficients."""
        arr = _arr()
        c = _coeffs(arr, levels=1)
        t = threshold_coeffs(c, threshold=0.0, mode="soft")
        for key in ("LH", "HL", "HH"):
            np.testing.assert_allclose(t[1][key], c[1][key], atol=1e-12)


# ---------------------------------------------------------------------------
# estimate_sigma
# ---------------------------------------------------------------------------

class TestEstimateSigma:

    def test_returns_float(self):
        band = np.random.default_rng(0).standard_normal(256)
        assert isinstance(estimate_sigma(band), float)

    def test_nonnegative(self):
        band = np.random.default_rng(1).standard_normal(256)
        assert estimate_sigma(band) >= 0.0

    def test_known_sigma_approx(self):
        """For N(0, sigma^2) noise, MAD estimator should be within 10% of sigma."""
        rng = np.random.default_rng(2)
        true_sigma = 0.5
        noise = rng.normal(0.0, true_sigma, 4096)
        est = estimate_sigma(noise)
        assert abs(est - true_sigma) < 0.05 * true_sigma + 0.01

    def test_zero_band_gives_zero_sigma(self):
        """All-zero array: median(|x|) = 0 → sigma = 0."""
        band = np.zeros(100)
        assert estimate_sigma(band) == pytest.approx(0.0, abs=1e-10)

    def test_accepts_2d_input(self):
        band = np.random.default_rng(3).standard_normal((16, 16))
        sig = estimate_sigma(band)
        assert sig >= 0.0

    def test_scale_proportional(self):
        """Doubling the noise amplitude should (roughly) double the estimate."""
        rng = np.random.default_rng(4)
        noise = rng.standard_normal(2048)
        s1 = estimate_sigma(noise)
        s2 = estimate_sigma(noise * 2.0)
        assert s2 == pytest.approx(s1 * 2.0, rel=0.05)

    def test_against_finest_hh_band(self):
        """Estimate from the HH sub-band should be positive for noisy input."""
        noisy = _arr() + np.random.default_rng(5).normal(0, 0.05, (ROWS, COLS))
        c = dwt2d(noisy, levels=1, wavelet="db4")
        sig = estimate_sigma(c[1]["HH"])
        assert sig > 0.0


# ---------------------------------------------------------------------------
# bayes_threshold
# ---------------------------------------------------------------------------

class TestBayesThreshold:

    def test_returns_float(self):
        band = np.random.default_rng(0).standard_normal(64)
        t = bayes_threshold(band, sigma_n=0.3)
        assert isinstance(t, float)

    def test_nonnegative(self):
        band = np.random.default_rng(1).standard_normal(64)
        assert bayes_threshold(band, sigma_n=0.1) >= 0.0

    def test_pure_noise_returns_max_abs(self):
        """When signal variance <= 0, threshold = max|x| (zero out entire band)."""
        rng = np.random.default_rng(2)
        sigma_n = 0.5
        # tiny band: variance will be dominated by sigma_n
        band = rng.normal(0, sigma_n * 0.5, 16)
        thr  = bayes_threshold(band, sigma_n=sigma_n)
        assert thr >= 0.0

    def test_high_snr_gives_small_threshold(self):
        """Large signal variance → small threshold relative to max."""
        rng = np.random.default_rng(3)
        band    = rng.normal(0, 5.0, 1024)   # high signal
        sigma_n = 0.01                         # tiny noise
        thr = bayes_threshold(band, sigma_n=sigma_n)
        assert thr < np.max(np.abs(band))

    def test_accepts_2d_band(self):
        band = np.random.default_rng(4).standard_normal((16, 16))
        t = bayes_threshold(band, sigma_n=0.2)
        assert t >= 0.0

    def test_zero_sigma_n_gives_zero_threshold(self):
        """sigma_n=0 → T = 0/sigma_s = 0 for any band with variance."""
        band = np.random.default_rng(5).standard_normal(64)
        t = bayes_threshold(band, sigma_n=0.0)
        assert t == pytest.approx(0.0, abs=1e-10)


# ---------------------------------------------------------------------------
# bayes_denoise
# ---------------------------------------------------------------------------

class TestBayesDenoise:

    def test_output_shape_preserved(self):
        arr = _arr()
        d = bayes_denoise(arr, levels=2, wavelet="db4")
        assert d.shape == arr.shape

    def test_output_dtype_float64(self):
        d = bayes_denoise(_arr(), levels=1)
        assert d.dtype == np.float64

    def test_output_c_contiguous(self):
        d = bayes_denoise(_arr(), levels=1)
        assert d.flags["C_CONTIGUOUS"]

    def test_denoising_reduces_noise(self):
        """BayesShrink denoising should reduce noise on a smooth (sparse) signal.

        BayesShrink is designed for signals whose energy is concentrated in
        few wavelet coefficients (sparse signals). A smooth sinusoidal surface
        satisfies this — its high-frequency sub-bands are near-zero and the
        noise sits above them, so BayesShrink correctly zeros the noisy coeffs.
        """
        rng = np.random.default_rng(0)
        xs = np.linspace(0, 4 * np.pi, ROWS)
        ys = np.linspace(0, 4 * np.pi, COLS)
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        clean = 0.4 * np.sin(X) * np.cos(Y) + 0.5   # smooth, values in [0.1, 0.9]
        noise = rng.normal(0, 0.05, (ROWS, COLS))
        noisy = clean + noise
        denoised = bayes_denoise(noisy, levels=3, wavelet="db4")
        assert np.std(denoised - clean) < np.std(noisy - clean)

    def test_clean_signal_not_severely_distorted(self):
        """Denoising a zero-noise signal should not add large artefacts."""
        arr = _arr()
        d = bayes_denoise(arr, levels=2, wavelet="haar")
        # max deviation should be small relative to signal range
        assert np.max(np.abs(d - arr)) < 1.0

    @pytest.mark.parametrize("wavelet", ALL_WAVELETS)
    def test_all_wavelets_produce_output(self, wavelet):
        d = bayes_denoise(_arr(), levels=2, wavelet=wavelet)
        assert d.shape == (ROWS, COLS)

    def test_levels_1_also_works(self):
        d = bayes_denoise(_arr(), levels=1, wavelet="db4")
        assert d.shape == (ROWS, COLS)


# ---------------------------------------------------------------------------
# band_energy
# ---------------------------------------------------------------------------

class TestBandEnergy:

    def test_returns_float(self):
        assert isinstance(band_energy(np.ones(4)), float)

    def test_known_value(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        assert band_energy(arr) == pytest.approx(30.0, rel=1e-10)

    def test_zeros_give_zero(self):
        assert band_energy(np.zeros(100)) == pytest.approx(0.0, abs=1e-12)

    def test_nonnegative(self):
        arr = np.random.default_rng(0).standard_normal(256)
        assert band_energy(arr) >= 0.0

    def test_2d_input_accepted(self):
        arr = np.ones((8, 8))
        assert band_energy(arr) == pytest.approx(64.0, rel=1e-10)

    def test_scaling(self):
        """Scaling by k multiplies energy by k^2."""
        arr = np.random.default_rng(1).random(64)
        e1 = band_energy(arr)
        e3 = band_energy(arr * 3.0)
        assert e3 == pytest.approx(e1 * 9.0, rel=1e-10)

    def test_energy_from_subband(self):
        c = _coeffs(levels=2)
        e = band_energy(c[1]["HH"])
        assert e >= 0.0


# ---------------------------------------------------------------------------
# band_stats
# ---------------------------------------------------------------------------

class TestBandStats:

    def test_returns_dict_with_required_keys(self):
        s = band_stats(np.ones(10))
        assert set(s.keys()) == {"mean", "var", "l1", "linf"}

    def test_constant_array(self):
        arr = np.full(1000, 7.0)
        s = band_stats(arr)
        assert s["mean"]  == pytest.approx(7.0,  rel=1e-10)
        assert s["var"]   == pytest.approx(0.0,  abs=1e-10)
        assert s["l1"]    == pytest.approx(7.0,  rel=1e-10)
        assert s["linf"]  == pytest.approx(7.0,  rel=1e-10)

    def test_known_values(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(2.5)
        assert s["var"]  == pytest.approx(1.25)
        assert s["l1"]   == pytest.approx(2.5)
        assert s["linf"] == pytest.approx(4.0)

    def test_symmetric_zero_mean(self):
        arr = np.array([-1.0, 1.0])
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(0.0, abs=1e-10)
        assert s["l1"]   == pytest.approx(1.0)
        assert s["linf"] == pytest.approx(1.0)

    def test_var_nonnegative(self):
        arr = np.random.default_rng(0).random(512)
        assert band_stats(arr)["var"] >= 0.0

    def test_linf_is_max_abs(self):
        arr = np.array([-5.0, 2.0, 3.0])
        assert band_stats(arr)["linf"] == pytest.approx(5.0)

    def test_2d_input_ravelled(self):
        arr = np.full((8, 8), 4.0)
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(4.0)

    def test_l1_equals_mean_for_nonneg_array(self):
        arr = np.random.default_rng(1).random(256)
        s = band_stats(arr)
        assert s["l1"] == pytest.approx(s["mean"], rel=1e-10)

    def test_values_are_python_floats(self):
        s = band_stats(np.ones(10))
        for v in s.values():
            assert isinstance(v, float)


# ---------------------------------------------------------------------------
# dwt2d_batch / idwt2d_batch
# ---------------------------------------------------------------------------

class TestDwt2dBatch:

    def _stack(self, n_bands=4, rows=ROWS, cols=COLS, seed=0) -> np.ndarray:
        return np.random.default_rng(seed).random((n_bands, rows, cols))

    def test_returns_list_of_dicts(self):
        result = dwt2d_batch(self._stack(), levels=1)
        assert isinstance(result, list)
        assert all(isinstance(c, dict) for c in result)

    def test_list_length_equals_n_bands(self):
        n_bands = 6
        result = dwt2d_batch(self._stack(n_bands=n_bands), levels=1)
        assert len(result) == n_bands

    def test_each_dict_has_correct_levels(self):
        levels = 2
        result = dwt2d_batch(self._stack(), levels=levels)
        for cd in result:
            assert set(cd.keys()) == {1, 2}

    def test_2d_input_raises(self):
        with pytest.raises(ValueError, match="3-D"):
            dwt2d_batch(np.ones((ROWS, COLS)), levels=1)

    @pytest.mark.parametrize("wavelet", ALL_WAVELETS)
    def test_roundtrip_all_wavelets(self, wavelet):
        # Verify internal consistency: each band in the batch result must
        # match an independent per-band dwt2d → idwt2d round-trip.
        stack = self._stack()
        n_bands = stack.shape[0]
        coeffs_list = dwt2d_batch(stack, levels=2, wavelet=wavelet)
        recon = idwt2d_batch(coeffs_list, wavelet=wavelet)
        expected = np.stack(
            [idwt2d(dwt2d(stack[b], levels=2, wavelet=wavelet), wavelet=wavelet)
             for b in range(n_bands)]
        )
        np.testing.assert_allclose(recon, expected, atol=1e-9)

    def test_roundtrip_output_shape(self):
        n_bands = 4
        stack = self._stack(n_bands=n_bands)
        cl = dwt2d_batch(stack, levels=2)
        r  = idwt2d_batch(cl)
        assert r.shape == stack.shape

    def test_roundtrip_dtype_float64(self):
        stack = self._stack()
        cl = dwt2d_batch(stack, levels=1)
        r  = idwt2d_batch(cl)
        assert r.dtype == np.float64

    def test_roundtrip_c_contiguous(self):
        stack = self._stack()
        cl = dwt2d_batch(stack, levels=1)
        r  = idwt2d_batch(cl)
        assert r.flags["C_CONTIGUOUS"]

    def test_single_band_batch_matches_dwt2d(self):
        """dwt2d_batch on 1 band should match dwt2d + idwt2d directly."""
        stack = self._stack(n_bands=1)
        cl = dwt2d_batch(stack, levels=2, wavelet="db4")
        recon_batch  = idwt2d_batch(cl, wavelet="db4")
        recon_single = idwt2d(dwt2d(stack[0], levels=2, wavelet="db4"), wavelet="db4")
        np.testing.assert_allclose(recon_batch[0], recon_single, atol=1e-9)

    def test_bands_processed_independently(self):
        """Modifying one band's coefficients should not affect others."""
        stack = self._stack(n_bands=3)
        cl = dwt2d_batch(stack, levels=1, wavelet="haar")
        # Zero out band 0 detail completely
        cl[0][1]["HH"][:] = 0.0
        recon = idwt2d_batch(cl, wavelet="haar")
        # Bands 1 and 2 must match their independent per-band round-trip
        for b in (1, 2):
            ref = idwt2d(dwt2d(stack[b], levels=1, wavelet="haar"), wavelet="haar")
            np.testing.assert_allclose(recon[b], ref, atol=1e-9)


# ---------------------------------------------------------------------------
# dwt3d / idwt3d
# ---------------------------------------------------------------------------

class TestDwt3d:

    def _cube(self, n_times=4, rows=16, cols=16, seed=0) -> np.ndarray:
        # n_times=4: log2(4)=2, so default levels=2 in idwt3d matches dwt3d(levels=2)
        return np.random.default_rng(seed).random((n_times, rows, cols))

    def test_output_shape_equals_input(self):
        cube = self._cube()
        c = dwt3d(cube, levels=1, wavelet="haar")
        assert c.shape == cube.shape

    def test_output_dtype_float64(self):
        c = dwt3d(self._cube(), levels=1)
        assert c.dtype == np.float64

    def test_output_c_contiguous(self):
        c = dwt3d(self._cube(), levels=1)
        assert c.flags["C_CONTIGUOUS"]

    def test_2d_input_raises(self):
        with pytest.raises(ValueError, match="3-D"):
            dwt3d(np.ones((16, 16)), levels=1)

    def test_n_times_not_divisible_raises(self):
        # 6 times, levels=2 → factor=4 → 6 % 4 != 0
        with pytest.raises(ValueError, match="n_times"):
            dwt3d(self._cube(n_times=6), levels=2)

    def test_rows_not_divisible_raises(self):
        with pytest.raises(ValueError):
            dwt3d(np.ones((8, 15, 16)), levels=2)

    @pytest.mark.parametrize("wavelet", ["haar", "db4", "sym4"])
    def test_roundtrip(self, wavelet):
        # idwt3d(dwt3d(cube, levels=2)) must reconstruct cube to machine epsilon.
        # n_times=4 → log2(4)=2 so idwt3d default levels matches.
        cube = self._cube()   # shape (4, 16, 16)
        c = dwt3d(cube, levels=2, wavelet=wavelet)
        r = idwt3d(c,   wavelet=wavelet)
        np.testing.assert_allclose(r, cube, atol=1e-9)

    def test_roundtrip_shape(self):
        cube = self._cube()   # (4, 16, 16)
        c = dwt3d(cube, levels=2)
        r = idwt3d(c)
        assert r.shape == cube.shape

    def test_roundtrip_dtype_float64(self):
        cube = self._cube()
        c = dwt3d(cube, levels=2)
        r = idwt3d(c)
        assert r.dtype == np.float64

    def test_roundtrip_c_contiguous(self):
        cube = self._cube()
        c = dwt3d(cube, levels=2)
        r = idwt3d(c)
        assert r.flags["C_CONTIGUOUS"]

    def test_temporal_energy_conserved_in_roundtrip(self):
        """idwt3d(dwt3d(cube)) must have the same total energy as cube (Parseval)."""
        cube = self._cube()   # (4, 16, 16)
        c    = dwt3d(cube, levels=2, wavelet="haar")
        r    = idwt3d(c,   wavelet="haar")
        assert np.sum(r ** 2) == pytest.approx(np.sum(cube ** 2), rel=1e-6)

    def test_zeroing_temporal_hi_smooths_cube(self):
        """Zeroing temporal high-band slices must reduce per-pixel temporal variance."""
        # n_times=4, levels=1: lo slices are [0,1], hi slices are [2,3].
        # After zeroing hi: each time-pair becomes its mean → temporal std drops.
        cube = self._cube(n_times=4, rows=16, cols=16)
        c = dwt3d(cube, levels=1, wavelet="haar")
        half = cube.shape[0] // 2   # = 2
        c_smooth = c.copy()
        c_smooth[half:] = 0.0       # zero temporal hi sub-band
        reconstructed = idwt3d(c_smooth, levels=1, wavelet="haar")
        # Temporal std per pixel (axis=0): zeroing hi collapses pairs → std < original
        orig_temporal_std  = np.std(cube,          axis=0).mean()
        recon_temporal_std = np.std(reconstructed, axis=0).mean()
        assert recon_temporal_std <= orig_temporal_std + 1e-6

    def test_idwt3d_2d_input_raises(self):
        with pytest.raises(ValueError, match="3-D"):
            idwt3d(np.ones((16, 16)))


# ---------------------------------------------------------------------------
# idwt2d — input validation
# ---------------------------------------------------------------------------

class TestIdwt2dValidation:

    def test_wavelet_mismatch_does_not_crash(self):
        """Using the wrong wavelet gives bad reconstruction but must not raise."""
        c = _coeffs(levels=1, wavelet="haar")
        # db4 synthesis applied to Haar coefficients: wrong result, not an exception
        r = idwt2d(c, wavelet="db4")
        assert r.shape == (ROWS, COLS)

    def test_unknown_wavelet_raises(self):
        c = _coeffs()
        with pytest.raises(ValueError, match="Unknown wavelet"):
            idwt2d(c, wavelet="bior99")


# ---------------------------------------------------------------------------
# Integration — multi-level denoising pipeline
# ---------------------------------------------------------------------------

class TestDenoisingPipeline:

    def test_manual_bayes_matches_bayes_denoise(self):
        """Manual per-band BayesShrink should match the bayes_denoise wrapper."""
        arr = _arr(seed=99)
        levels  = 3
        wavelet = "sym4"

        # bayes_denoise one-liner
        d_auto = bayes_denoise(arr, levels=levels, wavelet=wavelet)

        # Manual pipeline
        coeffs  = dwt2d(arr, levels=levels, wavelet=wavelet)
        sigma_n = estimate_sigma(coeffs[1]["HH"])
        manual  = {}
        for lv, sub_bands in coeffs.items():
            manual[lv] = {}
            for key, band in sub_bands.items():
                if key == "LL":
                    manual[lv]["LL"] = band.copy()
                else:
                    thr = bayes_threshold(band, sigma_n)
                    manual[lv][key] = (
                        np.sign(band) * np.maximum(np.abs(band) - thr, 0.0)
                    )
        d_manual = idwt2d(manual, wavelet=wavelet)

        np.testing.assert_allclose(d_auto, d_manual, atol=1e-12)

    def test_hard_threshold_idempotent(self):
        """Hard threshold is idempotent: applying it twice gives the same result.

        For hard thresholding: if |x| > t then x is kept unchanged,
        so a second pass with the same t leaves the result unchanged.
        """
        arr     = _arr()
        coeffs  = _coeffs(arr, levels=2)
        once    = threshold_coeffs(coeffs, threshold=0.1, mode="hard")
        twice   = threshold_coeffs(once,   threshold=0.1, mode="hard")
        r1 = idwt2d(once,  wavelet="haar")
        r2 = idwt2d(twice, wavelet="haar")
        np.testing.assert_allclose(r1, r2, atol=1e-9)

    def test_energy_decreases_after_threshold(self):
        """Thresholding detail sub-bands must not increase total band energy."""
        c   = _coeffs(levels=2)
        t   = threshold_coeffs(c, threshold=0.05, mode="soft")
        for lv in (1, 2):
            for key in ("LH", "HL", "HH"):
                assert band_energy(t[lv][key]) <= band_energy(c[lv][key]) + 1e-12

    def test_batch_denoise_each_band(self):
        """Batch DWT + per-band BayesShrink should reduce noise in every band.

        Uses smooth sinusoidal signals (sparse in wavelet domain) so BayesShrink
        reliably suppresses the added Gaussian noise. Spatial frequencies are kept
        well below Nyquist (period > 16 px) to avoid aliasing artefacts.
        """
        rng     = np.random.default_rng(7)
        n_bands = 3    # frequencies 1-3; frequency 4 has period=12 px ≈ filter length
        xs = np.linspace(0, 4 * np.pi, ROWS)
        ys = np.linspace(0, 4 * np.pi, COLS)
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        clean = np.stack([
            0.4 * np.sin((b + 1) * X) * np.cos((b + 1) * Y) + 0.5
            for b in range(n_bands)
        ])
        noisy   = clean + rng.normal(0, 0.05, clean.shape)

        cl = dwt2d_batch(noisy, levels=2, wavelet="db4")

        thresholded = []
        for cd in cl:
            sigma_n = estimate_sigma(cd[1]["HH"])
            t = {}
            for lv, sub_bands in cd.items():
                t[lv] = {}
                for key, band in sub_bands.items():
                    if key == "LL":
                        t[lv]["LL"] = band.copy()
                    else:
                        thr = bayes_threshold(band, sigma_n)
                        t[lv][key] = (
                            np.sign(band) * np.maximum(np.abs(band) - thr, 0.0)
                        )
            thresholded.append(t)

        denoised = idwt2d_batch(thresholded, wavelet="db4")

        # Each band: denoised is closer to clean than noisy is
        for b in range(n_bands):
            err_noisy    = np.std(noisy[b]    - clean[b])
            err_denoised = np.std(denoised[b] - clean[b])
            assert err_denoised < err_noisy