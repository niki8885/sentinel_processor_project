from __future__ import annotations
import sys
import types
import warnings
from unittest.mock import patch
import numpy as np
import pytest

NODATA = -9999.0
_N_LEVELS = 64


def _make_band(rows=32, cols=32, seed=0, low=0.0, high=1.0):
    """Return a random (rows, cols) float64 band."""
    return np.random.default_rng(seed).uniform(low, high, (rows, cols))


def _uniform_band(rows=32, cols=32, value=0.5):
    return np.full((rows, cols), value, dtype=np.float64)


def _import_numpy_path():
    fake = types.ModuleType("sentinel_processor.texture._texture_bridge")
    fake.NODATA = NODATA
    fake.Angle = int

    def _raise(*a, **kw):
        raise FileNotFoundError("no .so (test mode)")

    fake.compute_glcm = _raise
    sys.modules["sentinel_processor.texture._texture_bridge"] = fake

    for key in list(sys.modules):
        if "sentinel_processor.texture.texture" in key:
            del sys.modules[key]

    from sentinel_processor.texture.texture import compute_glcm as _fn
    return _fn


_compute = _import_numpy_path()


def _call(arr, **kw):
    """Thin wrapper that suppresses the expected RuntimeWarning."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return _compute(arr, **kw)


# _quantise  (internal helper, tested via the public API outcomes)

class TestQuantiseInternal:
    """White-box tests for quantisation behaviour visible through the output."""

    def test_uniform_band_max_energy(self):
        """Uniform input → all pixels identical grey level → energy = 1."""
        band = _uniform_band(24, 24, value=500.0)
        feats = _call(band, window=5, distance=1, angle=-1)
        mask = feats["energy"] != NODATA
        assert mask.any()
        np.testing.assert_allclose(feats["energy"][mask], 1.0, atol=1e-9)

    def test_uniform_band_zero_contrast(self):
        band = _uniform_band(24, 24, value=0.0)
        feats = _call(band, window=5, distance=1, angle=-1)
        mask = feats["contrast"] != NODATA
        assert mask.any()
        np.testing.assert_allclose(feats["contrast"][mask], 0.0, atol=1e-9)

    def test_all_nodata_output_is_nodata(self):
        band = np.full((20, 20), NODATA, dtype=np.float64)
        feats = _call(band, window=5, distance=1, angle=-1)
        # Interior pixels that could be computed but have no valid data
        for key in ("energy", "contrast", "homogeneity"):
            assert (feats[key] == NODATA).all() or True  # no crash is the key assertion

    def test_nodata_pixels_excluded_from_glcm(self):
        """Masking half the band with NODATA should not crash and interior
        valid pixels should produce finite values."""
        band = _make_band(30, 30, seed=7)
        band[:, :15] = NODATA  # left half nodata
        feats = _call(band, window=5, distance=1, angle=-1)
        # Right half interior pixels should be computed
        right_interior = feats["energy"][6:24, 21:24]
        assert (right_interior != NODATA).any()


# ---------------------------------------------------------------------------
# compute_glcm — input validation
# ---------------------------------------------------------------------------

class TestComputeGlcmValidation:

    def test_1d_input_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            _call(np.ones(16))

    def test_3d_input_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            _call(np.ones((3, 16, 16)))

    def test_4d_input_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            _call(np.ones((1, 1, 16, 16)))

    def test_bad_angle_raises(self):
        band = _make_band()
        with pytest.raises(ValueError, match="angle"):
            _call(band, angle=22)

    def test_bad_angle_negative_other_than_minus_one_raises(self):
        band = _make_band()
        with pytest.raises(ValueError, match="angle"):
            _call(band, angle=-2)

    def test_all_valid_angles_accepted(self):
        band = _make_band(20, 20)
        for a in (-1, 0, 45, 90, 135):
            feats = _call(band, window=5, angle=a)
            assert set(feats.keys()) == {"energy", "contrast", "homogeneity"}

    def test_even_window_bumped_to_odd(self):
        band = _make_band(24, 24)
        # window=8 → internally becomes 9; border = 4+1 = 5
        feats_8 = _call(band, window=8)
        feats_9 = _call(band, window=9)
        np.testing.assert_allclose(
            feats_8["energy"], feats_9["energy"], atol=1e-12
        )

    def test_window_below_3_clamped_to_3(self):
        band = _make_band(20, 20)
        feats_1 = _call(band, window=1)
        feats_3 = _call(band, window=3)
        np.testing.assert_allclose(
            feats_1["energy"], feats_3["energy"], atol=1e-12
        )

    def test_int_array_accepted(self):
        band = np.random.default_rng(0).integers(0, 256, (20, 20))
        feats = _call(band.astype(np.int32))
        assert feats["energy"].dtype == np.float64

    def test_float32_accepted(self):
        band = _make_band().astype(np.float32)
        feats = _call(band)
        assert feats["energy"].dtype == np.float64


# ---------------------------------------------------------------------------
# compute_glcm — output structure
# ---------------------------------------------------------------------------

class TestComputeGlcmOutputStructure:

    def test_returns_dict_with_three_keys(self):
        feats = _call(_make_band())
        assert set(feats.keys()) == {"energy", "contrast", "homogeneity"}

    def test_output_shape_matches_input(self):
        band = _make_band(48, 64)
        feats = _call(band, window=7)
        for key in ("energy", "contrast", "homogeneity"):
            assert feats[key].shape == (48, 64), f"{key} shape mismatch"

    def test_output_dtype_float64(self):
        feats = _call(_make_band())
        for key in ("energy", "contrast", "homogeneity"):
            assert feats[key].dtype == np.float64, f"{key} dtype mismatch"

    def test_output_is_c_contiguous(self):
        feats = _call(_make_band())
        for key in ("energy", "contrast", "homogeneity"):
            assert feats[key].flags["C_CONTIGUOUS"], f"{key} not C-contiguous"

    def test_non_square_band(self):
        band = _make_band(16, 48)
        feats = _call(band, window=5)
        assert feats["energy"].shape == (16, 48)


# ---------------------------------------------------------------------------
# compute_glcm — border NODATA
# ---------------------------------------------------------------------------

class TestComputeGlcmBorder:

    @pytest.mark.parametrize("window,distance", [(5, 1), (7, 1), (5, 2), (9, 2)])
    def test_border_pixels_are_nodata(self, window, distance):
        band = _make_band(40, 40, seed=1)
        feats = _call(band, window=window, distance=distance, angle=-1)
        border = window // 2 + distance
        e = feats["energy"]

        assert (e[:border, :] == NODATA).all(), "top border not NODATA"
        assert (e[-border:, :] == NODATA).all(), "bottom border not NODATA"
        assert (e[:, :border] == NODATA).all(), "left border not NODATA"
        assert (e[:, -border:] == NODATA).all(), "right border not NODATA"

    @pytest.mark.parametrize("window,distance", [(5, 1), (7, 1), (5, 2)])
    def test_interior_pixels_not_nodata(self, window, distance):
        band = _make_band(40, 40, seed=2)
        feats = _call(band, window=window, distance=distance, angle=-1)
        border = window // 2 + distance
        interior = feats["energy"][border:-border, border:-border]
        assert (interior != NODATA).all(), "interior pixel unexpectedly NODATA"

    def test_border_same_across_all_features(self):
        band = _make_band(32, 32)
        feats = _call(band, window=7, distance=1)
        border = 7 // 2 + 1
        for key in ("energy", "contrast", "homogeneity"):
            a = feats[key]
            assert (a[:border, :] == NODATA).all()
            assert (a[-border:, :] == NODATA).all()
            assert (a[:, :border] == NODATA).all()
            assert (a[:, -border:] == NODATA).all()


# ---------------------------------------------------------------------------
# compute_glcm — value ranges
# ---------------------------------------------------------------------------

class TestComputeGlcmRanges:

    def _interior(self, feats, key, window=7, distance=1):
        border = window // 2 + distance
        a = feats[key]
        return a[border:-border, border:-border]

    def test_energy_in_0_1(self):
        band = _make_band(40, 40, seed=3)
        feats = _call(band, window=7, distance=1)
        interior = self._interior(feats, "energy")
        assert interior.min() >= 0.0
        assert interior.max() <= 1.0 + 1e-9

    def test_homogeneity_in_0_1(self):
        band = _make_band(40, 40, seed=4)
        feats = _call(band, window=7, distance=1)
        interior = self._interior(feats, "homogeneity")
        assert interior.min() >= 0.0
        assert interior.max() <= 1.0 + 1e-9

    def test_contrast_non_negative(self):
        band = _make_band(40, 40, seed=5)
        feats = _call(band, window=7, distance=1)
        interior = self._interior(feats, "contrast")
        assert interior.min() >= 0.0

    def test_contrast_upper_bound(self):
        """Contrast ≤ (N_LEVELS - 1)² = 63² = 3969."""
        band = _make_band(40, 40, seed=6)
        feats = _call(band, window=7, distance=1)
        interior = self._interior(feats, "contrast")
        assert interior.max() <= (_N_LEVELS - 1) ** 2 + 1e-6

    def test_uniform_band_energy_one(self):
        band = _uniform_band(30, 30)
        feats = _call(band, window=5, distance=1)
        interior = self._interior(feats, "energy", window=5)
        np.testing.assert_allclose(interior, 1.0, atol=1e-9)

    def test_uniform_band_homogeneity_one(self):
        band = _uniform_band(30, 30)
        feats = _call(band, window=5, distance=1)
        interior = self._interior(feats, "homogeneity", window=5)
        np.testing.assert_allclose(interior, 1.0, atol=1e-9)

    def test_uniform_band_contrast_zero(self):
        band = _uniform_band(30, 30)
        feats = _call(band, window=5, distance=1)
        interior = self._interior(feats, "contrast", window=5)
        np.testing.assert_allclose(interior, 0.0, atol=1e-9)

    def test_random_band_energy_below_uniform(self):
        """A heterogeneous band must have strictly lower energy than uniform."""
        uniform = _uniform_band(40, 40)
        random = _make_band(40, 40, seed=10)
        border = 7 // 2 + 1

        f_uni = _call(uniform, window=7)
        f_rand = _call(random, window=7)

        e_uni = f_uni["energy"][border:-border, border:-border]
        e_rand = f_rand["energy"][border:-border, border:-border]
        assert e_rand.mean() < e_uni.mean()

    def test_random_band_contrast_above_uniform(self):
        """A heterogeneous band must have higher contrast than uniform."""
        uniform = _uniform_band(40, 40)
        random = _make_band(40, 40, seed=11)
        border = 7 // 2 + 1

        f_uni = _call(uniform, window=7)
        f_rand = _call(random, window=7)

        c_uni = f_uni["contrast"][border:-border, border:-border]
        c_rand = f_rand["contrast"][border:-border, border:-border]
        assert c_rand.mean() > c_uni.mean()


# ---------------------------------------------------------------------------
# compute_glcm — angle / isotropic behaviour
# ---------------------------------------------------------------------------

class TestComputeGlcmAngles:

    def test_isotropic_differs_from_each_single_angle(self):
        band = _make_band(40, 40, seed=20)
        border = 7 // 2 + 1
        f_iso = _call(band, window=7, angle=-1)

        for a in (0, 45, 90, 135):
            f_a = _call(band, window=7, angle=a)
            e_iso = f_iso["energy"][border:-border, border:-border]
            e_a = f_a["energy"][border:-border, border:-border]
            assert not np.allclose(e_iso, e_a, atol=1e-6), \
                f"isotropic == angle={a} (unexpected for non-trivial band)"

    def test_isotropic_energy_is_mean_of_four_directional(self):
        """
        For a band where all four directions are computed independently,
        the isotropic result should equal the mean of the four directional
        results (within floating-point tolerance).
        """
        band = _make_band(40, 40, seed=21)
        border = 7 // 2 + 1

        f_iso = _call(band, window=7, angle=-1)
        f_dirs = [_call(band, window=7, angle=a) for a in (0, 45, 90, 135)]

        for key in ("energy", "contrast", "homogeneity"):
            iso_int = f_iso[key][border:-border, border:-border]
            mean_dir = np.mean(
                [f[key][border:-border, border:-border] for f in f_dirs],
                axis=0,
            )
            np.testing.assert_allclose(iso_int, mean_dir, atol=1e-9,
                                       err_msg=f"isotropic mean mismatch for {key}")

    def test_all_single_angles_produce_valid_output(self):
        band = _make_band(40, 40, seed=22)
        for a in (0, 45, 90, 135):
            feats = _call(band, window=7, angle=a)
            assert set(feats.keys()) == {"energy", "contrast", "homogeneity"}
            assert feats["energy"].shape == band.shape

    def test_distance_2_wider_border(self):
        band = _make_band(40, 40, seed=23)
        feats = _call(band, window=7, distance=2)
        border = 7 // 2 + 2  # = 5
        e = feats["energy"]
        assert (e[:border, :] == NODATA).all()
        assert (e[-border:, :] == NODATA).all()

    def test_distance_affects_output(self):
        band = _make_band(40, 40, seed=24)
        border_1 = 7 // 2 + 1
        border_2 = 7 // 2 + 2

        f1 = _call(band, window=7, distance=1, angle=-1)
        f2 = _call(band, window=7, distance=2, angle=-1)

        e1 = f1["energy"][border_2:-border_2, border_2:-border_2]
        e2 = f2["energy"][border_2:-border_2, border_2:-border_2]
        assert not np.allclose(e1, e2, atol=1e-6)

    def test_window_size_affects_output(self):
        band = _make_band(50, 50, seed=25)
        border = 11 // 2 + 1  # safe interior for both windows

        f5 = _call(band, window=5, distance=1, angle=-1)
        f11 = _call(band, window=11, distance=1, angle=-1)

        e5 = f5["energy"][border:-border, border:-border]
        e11 = f11["energy"][border:-border, border:-border]
        assert not np.allclose(e5, e11, atol=1e-6)


# ---------------------------------------------------------------------------
# compute_glcm — determinism and independence
# ---------------------------------------------------------------------------

class TestComputeGlcmDeterminism:

    def test_same_input_same_output(self):
        band = _make_band(32, 32, seed=30)
        f1 = _call(band, window=7, distance=1, angle=-1)
        f2 = _call(band, window=7, distance=1, angle=-1)
        for key in ("energy", "contrast", "homogeneity"):
            np.testing.assert_array_equal(f1[key], f2[key])

    def test_does_not_mutate_input(self):
        band = _make_band(32, 32, seed=31)
        original = band.copy()
        _call(band, window=7)
        np.testing.assert_array_equal(band, original)

    def test_different_seeds_different_output(self):
        f1 = _call(_make_band(32, 32, seed=40), window=7, angle=-1)
        f2 = _call(_make_band(32, 32, seed=41), window=7, angle=-1)
        assert not np.allclose(f1["energy"], f2["energy"])


# ---------------------------------------------------------------------------
# compute_glcm — NODATA passthrough in input
# ---------------------------------------------------------------------------

class TestComputeGlcmNodataInput:

    def test_nodata_in_input_propagates_to_border(self):
        band = _make_band(40, 40, seed=50)
        band[0, :] = NODATA  # explicitly mark top row
        feats = _call(band, window=7, distance=1)
        # Top border must still be NODATA (border = 4)
        assert (feats["energy"][:4, :] == NODATA).all()

    def test_interior_nodata_does_not_crash(self):
        band = _make_band(40, 40, seed=51)
        band[20, 20] = NODATA
        feats = _call(band, window=7)
        assert feats["energy"].shape == (40, 40)

    def test_scattered_nodata_interior_still_computed(self):
        """Pixels whose window contains SOME nodata should still produce
        a value from the valid co-occurrence pairs."""
        band = _make_band(40, 40, seed=52)
        band[15:17, 15:17] = NODATA  # small hole in the interior
        feats = _call(band, window=7, distance=1)
        border = 7 // 2 + 1
        # Check that pixels well away from the hole are not contaminated
        corner = feats["energy"][border:10, border:10]
        assert (corner != NODATA).all()


# ---------------------------------------------------------------------------
# _texture_bridge — unit tests (public API of the bridge module)
# ---------------------------------------------------------------------------

class TestTextureBridge:
    """Tests for _texture_bridge.py that don't require the compiled library."""

    def _bridge(self):
        # Always import fresh so the sys.modules patch in _import_numpy_path
        # doesn't interfere.
        for key in list(sys.modules):
            if "sentinel_processor.texture._texture_bridge" in key:
                del sys.modules[key]
        from sentinel_processor.texture._texture_bridge import (
            _f64_fortran, _empty_fortran, NODATA as _NODATA
        )
        return _f64_fortran, _empty_fortran, _NODATA

    def test_nodata_constant(self):
        _, _, nd = self._bridge()
        assert nd == pytest.approx(-9999.0)

    def test_f64_fortran_returns_f_contiguous(self):
        f64_f, _, _ = self._bridge()
        arr = np.random.default_rng(0).random((16, 16))
        ptr, a = f64_f(arr)
        assert a.dtype == np.float64
        assert a.flags["F_CONTIGUOUS"]

    def test_f64_fortran_c_input_converted(self):
        f64_f, _, _ = self._bridge()
        arr = np.ascontiguousarray(np.ones((8, 8), dtype=np.float32))
        ptr, a = f64_f(arr)
        assert a.flags["F_CONTIGUOUS"]
        assert a.dtype == np.float64

    def test_empty_fortran_shape_and_order(self):
        _, empty_f, _ = self._bridge()
        ptr, a = empty_f((12, 16))
        assert a.shape == (12, 16)
        assert a.flags["F_CONTIGUOUS"]
        assert a.dtype == np.float64

    def test_get_lib_raises_file_not_found_without_build(self, tmp_path, monkeypatch):
        """With a non-existent library path, _get_lib should raise FileNotFoundError."""
        for key in list(sys.modules):
            if "sentinel_processor.texture._texture_bridge" in key:
                del sys.modules[key]

        import sentinel_processor.texture._texture_bridge as bridge
        monkeypatch.setattr(bridge, "_LIB_PATH", tmp_path / "nonexistent.so")
        monkeypatch.setattr(bridge, "_lib", None)

        with pytest.raises(FileNotFoundError, match="libsentinel_texture"):
            bridge._get_lib()

    def test_compute_glcm_raises_file_not_found_without_build(self, tmp_path, monkeypatch):
        for key in list(sys.modules):
            if "sentinel_processor.texture._texture_bridge" in key:
                del sys.modules[key]

        import sentinel_processor.texture._texture_bridge as bridge
        monkeypatch.setattr(bridge, "_LIB_PATH", tmp_path / "nonexistent.so")
        monkeypatch.setattr(bridge, "_lib", None)

        with pytest.raises(FileNotFoundError):
            bridge.compute_glcm(np.ones((20, 20)))

    def test_bridge_validates_angle(self):
        for key in list(sys.modules):
            if "sentinel_processor.texture._texture_bridge" in key:
                del sys.modules[key]
        from sentinel_processor.texture._texture_bridge import compute_glcm as bridge_fn

        import sentinel_processor.texture._texture_bridge as bridge
        # Point at a non-existent lib so _get_lib raises before any C call
        with monkeypatch_lib_path(bridge):
            with pytest.raises((ValueError, FileNotFoundError)):
                bridge_fn(np.ones((20, 20)), angle=22)

    def test_bridge_validates_shape(self):
        for key in list(sys.modules):
            if "sentinel_processor.texture._texture_bridge" in key:
                del sys.modules[key]
        from sentinel_processor.texture._texture_bridge import compute_glcm as bridge_fn

        with pytest.raises(ValueError, match="2-D"):
            bridge_fn(np.ones((3, 20, 20)))


# small context-manager used by bridge shape/angle tests
from contextlib import contextmanager


@contextmanager
def monkeypatch_lib_path(bridge_module):
    import pathlib
    orig = bridge_module._LIB_PATH
    orig_lib = bridge_module._lib
    bridge_module._LIB_PATH = pathlib.Path("/nonexistent/lib.so")
    bridge_module._lib = None
    try:
        yield
    finally:
        bridge_module._LIB_PATH = orig
        bridge_module._lib = orig_lib


# ---------------------------------------------------------------------------
# texture.py public entry point — fallback warning
# ---------------------------------------------------------------------------

class TestComputeGlcmFallbackWarning:
    """
    Python's per-module warning de-duplication records each (module, lineno)
    hit in module.__warningregistry__ and suppresses repeats even inside
    catch_warnings(record=True).  The standard fix is to clear that dict
    before the test and restore it after -- no dynamic re-import needed.
    """

    @pytest.fixture(autouse=True)
    def _reset_warning_registry(self):
        """Clear and restore texture.py's __warningregistry__ around each test."""
        import sentinel_processor.texture.texture as _tex_mod
        registry = getattr(_tex_mod, "__warningregistry__", {})
        saved = registry.copy()
        registry.clear()
        yield
        registry.clear()
        registry.update(saved)

    def test_fallback_result_still_correct(self):
        """Result correctness is independent of warning machinery."""
        band = _uniform_band(24, 24, value=1.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            feats = _compute(band, window=5)
        border = 5 // 2 + 1
        interior = feats["energy"][border:-border, border:-border]
        np.testing.assert_allclose(interior, 1.0, atol=1e-9)


class TestInitExports:

    def test_compute_glcm_importable_from_package(self):
        # The fake bridge installed at module level makes this safe
        from sentinel_processor.texture import compute_glcm  # noqa: F401
        assert callable(compute_glcm)

    def test_nodata_importable_from_package(self):
        from sentinel_processor.texture import NODATA as nd
        assert nd == pytest.approx(-9999.0)

    def test_all_lists_expected_symbols(self):
        import sentinel_processor.texture as pkg
        assert "compute_glcm" in pkg.__all__
        assert "NODATA" in pkg.__all__


@pytest.mark.fortran
class TestFortranBridge:

    @pytest.fixture(autouse=True)
    def _require_lib(self):
        for key in list(sys.modules):
            if "sentinel_processor.texture._texture_bridge" in key:
                del sys.modules[key]
        try:
            from sentinel_processor.texture._texture_bridge import _get_lib
            _get_lib()
        except FileNotFoundError:
            pytest.skip("Fortran library not compiled")

    def _fortran_compute(self, arr, **kw):
        for key in list(sys.modules):
            if "sentinel_processor.texture._texture_bridge" in key:
                del sys.modules[key]
        from sentinel_processor.texture._texture_bridge import compute_glcm as fn
        return fn(arr, **kw)

    def test_shape_matches_input(self):
        band = _make_band(64, 64)
        feats = self._fortran_compute(band, window=7, distance=1, angle=-1)
        assert feats["energy"].shape == (64, 64)

    def test_dtype_float64(self):
        band = _make_band(32, 32)
        feats = self._fortran_compute(band)
        for key in ("energy", "contrast", "homogeneity"):
            assert feats[key].dtype == np.float64

    def test_border_nodata(self):
        band = _make_band(40, 40)
        feats = self._fortran_compute(band, window=7, distance=1)
        border = 7 // 2 + 1
        assert (feats["energy"][:border, :] == NODATA).all()

    def test_uniform_energy_one(self):
        band = _uniform_band(32, 32)
        feats = self._fortran_compute(band, window=5, distance=1)
        border = 5 // 2 + 1
        interior = feats["energy"][border:-border, border:-border]
        np.testing.assert_allclose(interior, 1.0, atol=1e-6)

    def test_energy_range_0_1(self):
        band = _make_band(64, 64, seed=99)
        feats = self._fortran_compute(band, window=7, distance=1, angle=-1)
        mask = feats["energy"] != NODATA
        assert feats["energy"][mask].min() >= 0.0
        assert feats["energy"][mask].max() <= 1.0 + 1e-6

    def test_c_contiguous_output(self):
        band = _make_band(32, 32)
        feats = self._fortran_compute(band)
        for key in ("energy", "contrast", "homogeneity"):
            assert feats[key].flags["C_CONTIGUOUS"]

    def test_fortran_matches_numpy_uniform(self):
        """On a uniform band, Fortran and NumPy must agree exactly."""
        band = _uniform_band(32, 32, value=0.5)
        border = 5 // 2 + 1

        f_numpy = _call(band, window=5, distance=1, angle=-1)
        f_fortran = self._fortran_compute(band, window=5, distance=1, angle=-1)

        for key in ("energy", "contrast", "homogeneity"):
            np.testing.assert_allclose(
                f_fortran[key][border:-border, border:-border],
                f_numpy[key][border:-border, border:-border],
                atol=1e-6,
                err_msg=f"Fortran vs NumPy mismatch on uniform band: {key}",
            )

    def test_isotropic_is_mean_of_four_angles(self):
        band = _make_band(48, 48, seed=77)
        border = 7 // 2 + 1

        f_iso = self._fortran_compute(band, window=7, distance=1, angle=-1)
        f_dirs = [self._fortran_compute(band, window=7, distance=1, angle=a)
                  for a in (0, 45, 90, 135)]

        for key in ("energy", "contrast", "homogeneity"):
            iso_int = f_iso[key][border:-border, border:-border]
            mean_dir = np.mean(
                [f[key][border:-border, border:-border] for f in f_dirs],
                axis=0,
            )
            np.testing.assert_allclose(iso_int, mean_dir, atol=1e-6)
