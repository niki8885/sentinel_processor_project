from pathlib import Path

import numpy as np
import pytest

from sentinel_processor.indices._indices_bridge import (
    _ensure_f64c,
    _prepare,
    compute_arvi,
    compute_cig,
    compute_evi,
    compute_mndwi,
    compute_nbr,
    compute_ndbi,
    compute_ndsi,
    compute_ndvi,
    compute_ndwi,
    compute_savi,
)
from sentinel_processor.indices.compute import (
    AVAILABLE_INDICES,
    _BAND_ALIASES,
    _check_missing_bands,
    compute_indices,
    list_indices,
)


# _ensure_f64c / _prepare

class TestEnsureF64c:

    def test_list_input_converted(self):
        a = _ensure_f64c([1, 2, 3])
        assert a.dtype == np.float64

    def test_already_f64_no_copy(self):
        src = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        out = _ensure_f64c(src)
        assert out.ctypes.data == src.ctypes.data  # same buffer

    def test_float32_converted(self):
        src = np.array([1.0, 2.0], dtype=np.float32)
        out = _ensure_f64c(src)
        assert out.dtype == np.float64

    def test_2d_ravelled(self):
        src = np.ones((4, 4), dtype=np.float64)
        out = _ensure_f64c(src)
        assert out.ndim == 1
        assert out.size == 16

    def test_fortran_order_made_contiguous(self):
        src = np.asfortranarray(np.ones((4, 4), dtype=np.float64))
        out = _ensure_f64c(src)
        assert out.flags["C_CONTIGUOUS"]


class TestPrepare:

    def test_mismatched_sizes_raise(self):
        a = np.ones(10)
        b = np.ones(20)
        with pytest.raises(ValueError, match="same number"):
            _prepare(a, b)

    def test_matched_sizes_ok(self):
        a = np.ones(100)
        b = np.ones(100)
        arrs, n = _prepare(a, b)
        assert n == 100
        assert len(arrs) == 2

    def test_three_band_prepare(self):
        a = b = c = np.ones(50)
        arrs, n = _prepare(a, b, c)
        assert n == 50
        assert len(arrs) == 3


# Index formula correctness

class TestNDVI:

    def test_pure_vegetation(self):
        """NIR=1, RED=0 → NDVI = (1-0)/(1+0) = 1.0"""
        nir = np.array([1.0])
        red = np.array([0.0])
        result = compute_ndvi(nir, red)
        assert result[0] == pytest.approx(1.0)

    def test_pure_soil(self):
        """NIR=0.5, RED=0.5 → NDVI = 0"""
        result = compute_ndvi(np.array([0.5]), np.array([0.5]))
        assert result[0] == pytest.approx(0.0, abs=1e-9)

    def test_water_negative(self):
        """NIR=0.1, RED=0.6 → NDVI negative"""
        result = compute_ndvi(np.array([0.1]), np.array([0.6]))
        assert result[0] < 0

    def test_range_bulk(self):
        """All pixels in a random array must be in [-1, 1]."""
        rng = np.random.default_rng(0)
        nir = rng.random(10_000)
        red = rng.random(10_000)
        out = compute_ndvi(nir, red)
        assert out.min() >= -1.0 - 1e-6
        assert out.max() <= 1.0 + 1e-6

    def test_formula(self):
        """Spot-check formula: NDVI = (NIR - RED) / (NIR + RED)."""
        rng = np.random.default_rng(42)
        nir = rng.random(500) + 0.01
        red = rng.random(500) + 0.01
        expected = (nir - red) / (nir + red)
        result = compute_ndvi(nir, red)
        np.testing.assert_allclose(result, expected, atol=1e-9)

    def test_output_shape(self):
        nir = np.ones(256 * 256)
        red = np.ones(256 * 256) * 0.5
        assert compute_ndvi(nir, red).shape == (256 * 256,)


class TestEVI:

    def test_formula(self):
        """EVI = 2.5*(NIR-RED)/(NIR + 6*RED - 7.5*BLUE + 1)"""
        rng = np.random.default_rng(1)
        nir = rng.random(200) + 0.1
        red = rng.random(200) * 0.3
        blue = rng.random(200) * 0.1
        expected = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1.0)
        result = compute_evi(nir, red, blue)
        np.testing.assert_allclose(result, expected, atol=1e-9)


class TestSAVI:

    def test_formula(self):
        """SAVI = 1.5*(NIR-RED)/(NIR+RED+0.5)"""
        rng = np.random.default_rng(2)
        nir = rng.random(200) + 0.1
        red = rng.random(200) * 0.3
        expected = 1.5 * (nir - red) / (nir + red + 0.5)
        result = compute_savi(nir, red)
        np.testing.assert_allclose(result, expected, atol=1e-9)


class TestNDWI:

    def test_formula(self):
        """NDWI = (GREEN - NIR) / (GREEN + NIR)"""
        rng = np.random.default_rng(3)
        green = rng.random(200) + 0.1
        nir = rng.random(200) + 0.1
        expected = (green - nir) / (green + nir)
        result = compute_ndwi(green, nir)
        np.testing.assert_allclose(result, expected, atol=1e-9)


class TestMNDWI:

    def test_formula(self):
        """MNDWI = (GREEN - SWIR1) / (GREEN + SWIR1)"""
        rng = np.random.default_rng(4)
        green = rng.random(200) + 0.1
        swir1 = rng.random(200) + 0.1
        expected = (green - swir1) / (green + swir1)
        result = compute_mndwi(green, swir1)
        np.testing.assert_allclose(result, expected, atol=1e-9)


class TestNDBI:

    def test_formula(self):
        """NDBI = (SWIR1 - NIR) / (SWIR1 + NIR)"""
        rng = np.random.default_rng(5)
        swir1 = rng.random(200) + 0.1
        nir = rng.random(200) + 0.1
        expected = (swir1 - nir) / (swir1 + nir)
        result = compute_ndbi(swir1, nir)
        np.testing.assert_allclose(result, expected, atol=1e-9)


class TestNBR:

    def test_formula(self):
        """NBR = (NIR - SWIR2) / (NIR + SWIR2)"""
        rng = np.random.default_rng(6)
        nir = rng.random(200) + 0.1
        swir2 = rng.random(200) + 0.1
        expected = (nir - swir2) / (nir + swir2)
        result = compute_nbr(nir, swir2)
        np.testing.assert_allclose(result, expected, atol=1e-9)


class TestNDSI:

    def test_formula(self):
        """NDSI = (GREEN - SWIR1) / (GREEN + SWIR1)"""
        rng = np.random.default_rng(7)
        green = rng.random(200) + 0.1
        swir1 = rng.random(200) + 0.1
        expected = (green - swir1) / (green + swir1)
        result = compute_ndsi(green, swir1)
        np.testing.assert_allclose(result, expected, atol=1e-9)

    def test_ndsi_equals_mndwi(self):
        """NDSI and MNDWI share the same formula — results must be identical."""
        rng = np.random.default_rng(8)
        green = rng.random(300) + 0.1
        swir1 = rng.random(300) + 0.1
        np.testing.assert_allclose(
            compute_ndsi(green, swir1),
            compute_mndwi(green, swir1),
            atol=1e-12,
        )


class TestCIG:

    def test_formula(self):
        """CIG = (NIR / GREEN) - 1"""
        rng = np.random.default_rng(9)
        nir = rng.random(200) + 0.1
        green = rng.random(200) + 0.1
        expected = (nir / green) - 1.0
        result = compute_cig(nir, green)
        np.testing.assert_allclose(result, expected, atol=1e-9)


class TestARVI:

    def test_formula(self):
        """ARVI = (NIR - (2*RED - BLUE)) / (NIR + (2*RED - BLUE))"""
        rng = np.random.default_rng(10)
        nir = rng.random(200) + 0.2
        red = rng.random(200) * 0.3
        blue = rng.random(200) * 0.1
        rb = 2.0 * red - blue
        expected = (nir - rb) / (nir + rb)
        result = compute_arvi(nir, red, blue)
        np.testing.assert_allclose(result, expected, atol=1e-9)


class TestBandAliases:

    def test_human_names_resolve(self):
        assert _BAND_ALIASES["nir"] == "B08"
        assert _BAND_ALIASES["red"] == "B04"
        assert _BAND_ALIASES["green"] == "B03"
        assert _BAND_ALIASES["blue"] == "B02"
        assert _BAND_ALIASES["swir16"] == "B11"
        assert _BAND_ALIASES["swir22"] == "B12"

    def test_canonical_passthrough(self):
        for key in ["B02", "B03", "B04", "B08", "B11", "B12"]:
            assert _BAND_ALIASES[key] == key

    def test_lowercase_canonical_resolves(self):
        assert _BAND_ALIASES["b08"] == "B08"

    def test_alternative_spellings(self):
        assert _BAND_ALIASES["swir1"] == "B11"
        assert _BAND_ALIASES["swir2"] == "B12"
        assert _BAND_ALIASES["nir_broad"] == "B08"
        assert _BAND_ALIASES["nir_narrow"] == "B8A"


class TestCheckMissingBands:

    def test_all_present_no_missing(self):
        available = {"B08", "B04"}
        missing = _check_missing_bands(available, ["ndvi"])
        assert missing.get("ndvi") == [] or missing.get("ndvi") is None

    def test_missing_one_band(self):
        available = {"B08"}  # missing B04
        missing = _check_missing_bands(available, ["ndvi"])
        assert "B04" in missing["ndvi"]

    def test_unknown_index_reported(self):
        missing = _check_missing_bands(set(), ["xyz_unknown"])
        assert "xyz_unknown" in missing

    def test_multiple_indices_independent(self):
        available = {"B08", "B04", "B03"}
        missing = _check_missing_bands(available, ["ndvi", "ndwi", "evi"])
        assert missing.get("ndvi") == [] or "ndvi" not in missing
        assert missing.get("ndwi") == [] or "ndwi" not in missing
        assert "B02" in missing.get("evi", [])


class TestAvailableIndices:

    def test_all_expected_present(self):
        expected = {"ndvi", "evi", "savi", "ndwi", "mndwi",
                    "ndbi", "nbr", "ndsi", "cig", "arvi"}
        assert expected.issubset(set(AVAILABLE_INDICES))

    def test_sorted(self):
        assert AVAILABLE_INDICES == sorted(AVAILABLE_INDICES)


class TestListIndices:

    def test_returns_dict(self):
        r = list_indices()
        assert isinstance(r, dict)

    def test_each_entry_has_bands_required(self):
        for name, entry in list_indices().items():
            assert "bands_required" in entry, f"{name} missing bands_required"
            assert isinstance(entry["bands_required"], list)

    def test_each_entry_has_long_name(self):
        for name, entry in list_indices().items():
            assert "long_name" in entry
            assert len(entry["long_name"]) > 0


class TestComputeIndices:

    def test_missing_source_raises(self):
        with pytest.raises(FileNotFoundError):
            compute_indices("/nonexistent.nc", ["ndvi"])

    def test_empty_indices_raises(self, tmp_path):
        f = tmp_path / "dummy.tif"
        f.touch()
        with pytest.raises(ValueError, match="must not be empty"):
            compute_indices(str(f), [])

    def test_unknown_index_raises(self, tmp_path):
        f = tmp_path / "dummy.tif"
        f.touch()
        with pytest.raises(ValueError, match="Unknown index"):
            compute_indices(str(f), ["foobar"])

    def test_unsupported_format_raises(self, tmp_path):
        f = tmp_path / "dummy.tif"
        f.touch()
        with pytest.raises(ValueError, match="Unsupported output_format"):
            compute_indices(str(f), ["ndvi"], output_format="jpg")

    def test_unsupported_extension_raises(self, tmp_path):
        f = tmp_path / "dummy.hdf"
        f.write_bytes(b"")
        with pytest.raises(ValueError, match="Unsupported file format"):
            compute_indices(str(f), ["ndvi"])

    def test_tif_with_named_bands(self, tmp_path):
        pytest.importorskip("rioxarray")
        pytest.importorskip("rasterio")
        import rasterio
        from rasterio.transform import from_bounds

        rows, cols = 32, 32
        nir = np.full((rows, cols), 0.8, dtype=np.float32)
        red = np.full((rows, cols), 0.2, dtype=np.float32)
        transform = from_bounds(0, 0, 1, 1, cols, rows)

        tif = tmp_path / "scene.tif"
        with rasterio.open(
                str(tif), "w",
                driver="GTiff", height=rows, width=cols, count=2,
                dtype="float32", crs="EPSG:4326", transform=transform,
        ) as dst:
            dst.write(nir, 1)
            dst.write(red, 2)
            dst.update_tags(1, name="nir")
            dst.update_tags(2, name="red")
        result = compute_indices(str(tif), ["ndvi"], output_dir=str(tmp_path))
        assert isinstance(result, dict)

    def test_overwrite_false_skips_existing(self, tmp_path):
        pytest.importorskip("rioxarray")
        pytest.importorskip("rasterio")
        import rasterio
        from rasterio.transform import from_bounds

        rows, cols = 32, 32
        transform = from_bounds(0, 0, 1, 1, cols, rows)
        tif = tmp_path / "scene.tif"
        with rasterio.open(
                str(tif), "w",
                driver="GTiff", height=rows, width=cols, count=1,
                dtype="float32", crs="EPSG:4326", transform=transform,
        ) as dst:
            dst.write(np.ones((rows, cols), dtype=np.float32), 1)

        idx_dir = tmp_path / "indices"
        idx_dir.mkdir()
        existing = idx_dir / "indices_scene_ndvi.tif"
        existing.write_bytes(b"placeholder")

        result = compute_indices(
            str(tif), ["ndvi"],
            output_dir=str(idx_dir),
            overwrite=False,
        )
        assert isinstance(result, dict)


# compute_indices — NetCDF integration


def _write_scene_nc(tmp_path, var_names=("B04", "B08"), rows=24, cols=24):
    """Write a NetCDF dataset with named band variables; return its path."""
    pytest.importorskip("netCDF4")
    import xarray as xr

    rng = np.random.default_rng(11)
    xs = np.linspace(600000.0, 600000.0 + cols * 10.0, cols, endpoint=False)
    ys = np.linspace(5200000.0, 5200000.0 - rows * 10.0, rows, endpoint=False)
    ds = xr.Dataset(
        {
            name: xr.DataArray(
                rng.uniform(100.0, 4000.0, (rows, cols)),
                dims=["y", "x"],
                coords={"y": ys, "x": xs},
            )
            for name in var_names
        }
    )
    path = tmp_path / "scene_20260101T000000.nc"
    ds.to_netcdf(str(path))
    ds.close()
    return path


class TestComputeIndicesValidation:

    def test_missing_source_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            compute_indices(str(tmp_path / "nope.nc"), ["ndvi"])

    def test_empty_indices_raises(self, tmp_path):
        src = _write_scene_nc(tmp_path)
        with pytest.raises(ValueError, match="must not be empty"):
            compute_indices(str(src), [])

    def test_unknown_index_raises(self, tmp_path):
        src = _write_scene_nc(tmp_path)
        with pytest.raises(ValueError, match="Unknown index"):
            compute_indices(str(src), ["nope"])

    def test_bad_output_format_raises(self, tmp_path):
        src = _write_scene_nc(tmp_path)
        with pytest.raises(ValueError, match="output_format"):
            compute_indices(str(src), ["ndvi"], output_format="png")

    def test_unsupported_suffix_raises(self, tmp_path):
        bad = tmp_path / "scene.txt"
        bad.write_text("not a raster")
        with pytest.raises(ValueError, match="Unsupported file format"):
            compute_indices(str(bad), ["ndvi"])


class TestComputeIndicesNetCDF:

    def test_ndvi_written_as_nc(self, tmp_path):
        src = _write_scene_nc(tmp_path)
        out_dir = tmp_path / "out"
        result = compute_indices(
            str(src), ["ndvi"], output_dir=str(out_dir), output_format="nc",
        )
        assert "ndvi" in result
        assert Path(result["ndvi"]).exists()
        assert Path(result["ndvi"]).suffix == ".nc"

    def test_ndvi_written_as_tif(self, tmp_path):
        src = _write_scene_nc(tmp_path)
        out_dir = tmp_path / "out"
        result = compute_indices(
            str(src), ["ndvi"], output_dir=str(out_dir), output_format="tif",
        )
        assert "ndvi" in result
        assert Path(result["ndvi"]).suffix == ".tif"

    def test_ndvi_values_in_valid_range(self, tmp_path):
        import xarray as xr
        src = _write_scene_nc(tmp_path)
        result = compute_indices(
            str(src), ["ndvi"], output_dir=str(tmp_path / "out"),
            output_format="nc",
        )
        with xr.open_dataarray(result["ndvi"]) as da:
            vals = da.values
        valid = vals[np.abs(vals - (-9999.0)) > 1.0]
        assert np.all(valid >= -1.0 - 1e-6)
        assert np.all(valid <= 1.0 + 1e-6)

    def test_missing_band_index_skipped(self, tmp_path):
        # mndwi needs B11 which the scene lacks; ndvi is still computed
        src = _write_scene_nc(tmp_path)
        result = compute_indices(
            str(src), ["ndvi", "mndwi"], output_dir=str(tmp_path / "out"),
        )
        assert "ndvi" in result
        assert "mndwi" not in result

    def test_all_bands_missing_returns_empty(self, tmp_path):
        src = _write_scene_nc(tmp_path, var_names=("B04",))
        result = compute_indices(
            str(src), ["mndwi"], output_dir=str(tmp_path / "out"),
        )
        assert result == {}

    def test_scale_factor_applied(self, tmp_path):
        src = _write_scene_nc(tmp_path)
        result = compute_indices(
            str(src), ["ndvi"], output_dir=str(tmp_path / "out"),
            scale_factor=1e-4,
        )
        assert "ndvi" in result

    def test_default_output_dir_next_to_source(self, tmp_path):
        src = _write_scene_nc(tmp_path)
        result = compute_indices(str(src), ["ndvi"])
        assert Path(result["ndvi"]).parent == tmp_path / "indices"

    def test_multiple_indices_one_call(self, tmp_path):
        src = _write_scene_nc(tmp_path, var_names=("B03", "B04", "B08"))
        result = compute_indices(
            str(src), ["ndvi", "savi", "ndwi"], output_dir=str(tmp_path / "out"),
        )
        assert set(result) == {"ndvi", "savi", "ndwi"}

    def test_overwrite_false_keeps_existing(self, tmp_path):
        src = _write_scene_nc(tmp_path)
        out_dir = tmp_path / "out"
        first = compute_indices(str(src), ["ndvi"], output_dir=str(out_dir))
        mtime = Path(first["ndvi"]).stat().st_mtime_ns
        second = compute_indices(
            str(src), ["ndvi"], output_dir=str(out_dir), overwrite=False,
        )
        assert Path(second["ndvi"]).stat().st_mtime_ns == mtime


class TestComputeIndicesTifWithoutBandNames:

    def test_unnamed_bands_return_empty(self, tmp_path):
        # Plain GeoTIFF band indices (1..N) cannot be matched to B04/B08
        pytest.importorskip("rasterio")
        import rasterio
        from rasterio.transform import from_bounds

        rows, cols = 24, 24
        tif = tmp_path / "plain.tif"
        with rasterio.open(
                str(tif), "w",
                driver="GTiff", height=rows, width=cols, count=2,
                dtype="float32", crs="EPSG:4326",
                transform=from_bounds(0, 0, 1, 1, cols, rows),
        ) as dst:
            dst.write(np.ones((2, rows, cols), dtype=np.float32))

        result = compute_indices(
            str(tif), ["ndvi"], output_dir=str(tmp_path / "out"),
        )
        assert result == {}
