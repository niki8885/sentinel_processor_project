from __future__ import annotations
import math
import sys
import tempfile
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from sentinel_processor.dl._normalize_bridge import (
    minmax_band,
    zscore_band,
    band_stats,
)
from sentinel_processor.dl.presets import (
    get_preset,
    AVAILABLE_PRESETS,
    SSL4EO_S12_ALL,
)
from sentinel_processor.dl.normalize import (
    normalize_for_dl,
    save_stats,
    load_stats,
)
from sentinel_processor.dl.tiling import (
    extract_tiles,
    stitch_tiles,
    tile_grid_shape,
)

NODATA = -9999.0


class TestMinmaxBand:
    def test_basic_range(self):
        arr = np.array([0.0, 500.0, 1000.0])
        out, vmin, vmax = minmax_band(arr)
        assert vmin == pytest.approx(0.0)
        assert vmax == pytest.approx(1000.0)
        np.testing.assert_allclose(out, [0.0, 0.5, 1.0])

    def test_nodata_preserved(self):
        arr = np.array([NODATA, 0.0, 1000.0])
        out, _, _ = minmax_band(arr)
        assert out[0] == pytest.approx(NODATA)

    def test_constant_band(self):
        arr = np.full(10, 500.0)
        out, _, _ = minmax_band(arr)
        np.testing.assert_allclose(out, 0.5)

    def test_all_nodata(self):
        arr = np.full(5, NODATA)
        out, vmin, vmax = minmax_band(arr)
        assert vmin == pytest.approx(0.0)
        np.testing.assert_allclose(out, NODATA)

    def test_2d_input_flattened(self):
        arr = np.arange(1, 10, dtype=np.float64).reshape(3, 3)
        out, vmin, vmax = minmax_band(arr)
        assert out.shape == (9,)
        assert vmin == pytest.approx(1.0)
        assert vmax == pytest.approx(9.0)


class TestZscoreBand:
    def test_basic(self):
        arr = np.array([0.0, 1.0, 2.0])
        out = zscore_band(arr, mean=1.0, std=1.0)
        np.testing.assert_allclose(out, [-1.0, 0.0, 1.0])

    def test_zero_std_returns_zeros(self):
        arr = np.array([5.0, 5.0, 5.0])
        out = zscore_band(arr, mean=5.0, std=0.0)
        np.testing.assert_allclose(out, 0.0)

    def test_nodata_preserved(self):
        arr = np.array([NODATA, 2.0])
        out = zscore_band(arr, mean=1.0, std=1.0)
        assert out[0] == pytest.approx(NODATA)
        assert out[1] == pytest.approx(1.0)


class TestBandStats:
    def test_known_values(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(3.0)
        assert s["std"] == pytest.approx(math.sqrt(2.0))
        assert s["min"] == pytest.approx(1.0)
        assert s["max"] == pytest.approx(5.0)

    def test_nodata_excluded(self):
        arr = np.array([NODATA, 2.0, 4.0])
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(3.0)
        assert s["min"] == pytest.approx(2.0)
        assert s["max"] == pytest.approx(4.0)

    def test_all_nodata(self):
        arr = np.full(5, NODATA)
        s = band_stats(arr)
        assert s["mean"] == pytest.approx(0.0)
        assert s["std"] == pytest.approx(0.0)


class TestNormalizeMinmax:
    def test_3d_output_range(self):
        rng = np.random.default_rng(42)
        arr = rng.uniform(0, 10000, (4, 64, 64))
        out = normalize_for_dl(arr, method="minmax")
        assert out.shape == (4, 64, 64)
        for b in range(4):
            assert out[b].min() >= -1e-9
            assert out[b].max() <= 1.0 + 1e-9

    def test_preserves_shape_1d(self):
        arr = np.array([100.0, 200.0, 300.0, 400.0])
        out = normalize_for_dl(arr, method="minmax")
        assert out.shape == arr.shape

    def test_preserves_shape_2d(self):
        arr = np.arange(12, dtype=np.float64).reshape(3, 4)
        out = normalize_for_dl(arr, method="minmax")
        assert out.shape == arr.shape

    def test_custom_nodata_preserved(self):
        """Pixels equal to a non-default nodata must be masked and preserved.

        Regression: a custom nodata used to be stretched as ordinary data.
        """
        arr = np.array([[0.0, 5.0, 10.0, -1.0]])
        out = normalize_for_dl(arr, method="minmax", nodata=-1.0)
        assert out[0, 3] == pytest.approx(-1.0)
        np.testing.assert_allclose(out[0, :3], [0.0, 0.5, 1.0])

    def test_custom_nodata_zscore(self):
        arr = np.array([[100.0, 200.0, 300.0, -1.0]])
        stats = {"0": {"mean": 200.0, "std": 100.0}}
        out = normalize_for_dl(arr, method="zscore", stats=stats, nodata=-1.0)
        assert out[0, 3] == pytest.approx(-1.0)
        np.testing.assert_allclose(out[0, :3], [-1.0, 0.0, 1.0])


class TestNormalizeZscore:
    def test_scene_stats(self):
        rng = np.random.default_rng(0)
        arr = rng.normal(1000, 200, (3, 32, 32))
        out = normalize_for_dl(arr, method="zscore")
        for b in range(3):
            assert abs(out[b].mean()) < 0.1

    def test_provided_stats(self):
        arr = np.array([[100.0, 200.0, 300.0],
                        [400.0, 500.0, 600.0]], dtype=np.float64)
        stats = {
            "0": {"mean": 200.0, "std": 100.0},
            "1": {"mean": 500.0, "std": 100.0},
        }
        out = normalize_for_dl(arr, method="zscore", stats=stats)
        np.testing.assert_allclose(out[0], [-1.0, 0.0, 1.0])
        np.testing.assert_allclose(out[1], [-1.0, 0.0, 1.0])

    def test_missing_band_raises(self):
        arr = np.ones((2, 4, 4))
        stats = {"B04": {"mean": 100.0, "std": 10.0}}
        with pytest.raises(KeyError):
            normalize_for_dl(arr, method="zscore", bands=["B04", "B03"], stats=stats)


class TestNormalizePresets:
    @pytest.mark.parametrize("method,bands", [
        ("sentinel2_rgb", ["B04", "B03", "B02"]),
        ("sentinel2_rgbn", ["B04", "B03", "B02", "B08"]),
        ("sentinel2_all", ["B02", "B03", "B04", "B08", "B05", "B06", "B07", "B8A", "B11", "B12"]),
    ])
    def test_output_shape(self, method, bands):
        C = len(bands)
        arr = np.full((C, 16, 16), 1000.0)
        out = normalize_for_dl(arr, method=method, bands=bands)
        assert out.shape == (C, 16, 16)

    def test_missing_bands_arg_raises(self):
        arr = np.ones((3, 8, 8))
        with pytest.raises(ValueError, match="'bands' must be provided"):
            normalize_for_dl(arr, method="sentinel2_rgb")

    def test_unknown_method_raises(self):
        arr = np.ones((3, 8, 8))
        with pytest.raises(ValueError, match="Unknown normalisation method"):
            normalize_for_dl(arr, method="foobar")


class TestStatsIO:
    def test_roundtrip(self):
        stats = {
            "B04": {"mean": 1354.99, "std": 2173.26},
            "B03": {"mean": 1117.36, "std": 2065.40},
        }
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "mystats.json"
            save_stats(stats, p)
            assert p.exists()
            loaded = load_stats(p)
        assert loaded["B04"]["mean"] == pytest.approx(1354.99)
        assert loaded["B03"]["std"] == pytest.approx(2065.40)

    def test_extension_added_automatically(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "stats_no_ext"
            save_stats({"B04": {"mean": 1.0, "std": 2.0}}, p)
            assert (Path(tmp) / "stats_no_ext.json").exists()

    def test_load_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_stats("/tmp/nonexistent_xyz_12345.json")


class TestPresets:
    def test_all_presets_available(self):
        for name in AVAILABLE_PRESETS:
            preset = get_preset(name)
            for band, s in preset.items():
                assert "mean" in s and "std" in s

    def test_unknown_preset_raises(self):
        with pytest.raises(ValueError):
            get_preset("imagenet")

    def test_ssl4eo_has_10_bands(self):
        assert len(SSL4EO_S12_ALL) == 10


class TestExtractTiles:
    def test_basic_grid_size(self):
        arr = np.zeros((4, 512, 512))
        tiles, meta = extract_tiles(arr, tile_size=256, overlap=32)
        assert tiles.shape == (meta.grid_rows * meta.grid_cols, 4, 256, 256)

    def test_exact_fit_no_padding(self):
        arr = np.ones((2, 512, 512))
        tiles, meta = extract_tiles(arr, tile_size=256, overlap=0)
        assert meta.grid_rows == 2
        assert meta.grid_cols == 2
        assert meta.pad_h == 0
        assert meta.pad_w == 0

    def test_non_square_image(self):
        arr = np.zeros((3, 300, 400))
        tiles, meta = extract_tiles(arr, tile_size=128, overlap=16)
        assert tiles.ndim == 4
        assert tiles.shape[1] == 3
        assert tiles.shape[2] == tiles.shape[3] == 128

    def test_2d_input_squeezed(self):
        arr = np.ones((100, 100))
        tiles, meta = extract_tiles(arr, tile_size=64, overlap=8)
        assert tiles.ndim == 3

    def test_invalid_params(self):
        arr = np.zeros((1, 64, 64))
        with pytest.raises(ValueError):
            extract_tiles(arr, tile_size=0)
        with pytest.raises(ValueError):
            extract_tiles(arr, tile_size=64, overlap=-1)
        with pytest.raises(ValueError):
            extract_tiles(arr, tile_size=64, overlap=32)


class TestStitchTiles:
    def _roundtrip(self, C, H, W, tile_size, overlap):
        rng = np.random.default_rng(99)
        arr = rng.uniform(0, 1, (C, H, W))
        tiles, meta = extract_tiles(arr, tile_size=tile_size, overlap=overlap)
        recon = stitch_tiles(tiles, meta, blend=True)
        assert recon.shape == (C, H, W)
        np.testing.assert_allclose(recon, arr, atol=1e-10)

    def test_exact_fit(self):
        self._roundtrip(C=3, H=512, W=512, tile_size=256, overlap=0)

    def test_with_overlap(self):
        self._roundtrip(C=4, H=256, W=256, tile_size=128, overlap=16)

    def test_non_square(self):
        self._roundtrip(C=2, H=300, W=400, tile_size=128, overlap=16)

    def test_small_image(self):
        arr = np.random.rand(1, 50, 200)
        tiles, meta = extract_tiles(arr, tile_size=128, overlap=16)
        recon = stitch_tiles(tiles, meta)
        assert recon.shape == (1, 50, 200)

    def test_2d_roundtrip(self):
        arr = np.random.rand(64, 64)
        tiles, meta = extract_tiles(arr, tile_size=32, overlap=4)
        recon = stitch_tiles(tiles, meta)
        assert recon.shape == (64, 64)


class TestTileGridShape:
    def test_exact_fit(self):
        n, rows, cols = tile_grid_shape(512, 512, tile_size=256, overlap=0)
        assert rows == 2 and cols == 2 and n == 4

    def test_overlap_increases_grid(self):
        n0, r0, c0 = tile_grid_shape(512, 512, 256, overlap=0)
        n1, r1, c1 = tile_grid_shape(512, 512, 256, overlap=64)
        assert n1 >= n0


# compute_dataset_stats / _welford_combine / save & load stats


def _write_stats_scene(tmp_path, name, offset=0.0, rows=16, cols=16):
    """Write a NetCDF scene with B04/B08 variables; return its path."""
    pytest.importorskip("netCDF4")
    import xarray as xr

    rng = np.random.default_rng(int(offset) + 3)
    ds = xr.Dataset(
        {
            "B04": (("y", "x"), rng.uniform(100, 2000, (rows, cols)) + offset),
            "B08": (("y", "x"), rng.uniform(500, 4000, (rows, cols)) + offset),
        }
    )
    path = tmp_path / name
    ds.to_netcdf(str(path))
    ds.close()
    return path


class TestWelfordCombine:

    def test_two_batches_match_global_stats(self):
        from sentinel_processor.dl.normalize import _welford_combine

        rng = np.random.default_rng(0)
        a = rng.normal(100, 20, 1000)
        b = rng.normal(300, 50, 500)

        agg = {"n": 0, "mean": 0.0, "M2": 0.0,
               "min": float("inf"), "max": float("-inf")}
        for batch in (a, b):
            _welford_combine(
                agg, batch.size,
                float(batch.mean()), float(batch.std()),
                float(batch.min()), float(batch.max()),
            )

        full = np.concatenate([a, b])
        assert agg["n"] == full.size
        assert agg["mean"] == pytest.approx(full.mean())
        assert math.sqrt(agg["M2"] / agg["n"]) == pytest.approx(full.std())
        assert agg["min"] == pytest.approx(full.min())
        assert agg["max"] == pytest.approx(full.max())


class TestComputeDatasetStats:

    def test_stats_across_scenes(self, tmp_path):
        from sentinel_processor.dl.normalize import compute_dataset_stats

        p1 = _write_stats_scene(tmp_path, "s_20260101T000000.nc", offset=0.0)
        p2 = _write_stats_scene(tmp_path, "s_20260111T000000.nc", offset=500.0)
        stats = compute_dataset_stats([p1, p2], ["B04", "B08"])

        for b in ("B04", "B08"):
            assert stats[b]["n_pixels"] == 2 * 16 * 16
            assert stats[b]["std"] > 0
            assert stats[b]["min"] <= stats[b]["mean"] <= stats[b]["max"]

    def test_missing_band_yields_zero_pixels(self, tmp_path):
        from sentinel_processor.dl.normalize import compute_dataset_stats

        p1 = _write_stats_scene(tmp_path, "s_20260101T000000.nc")
        stats = compute_dataset_stats([p1], ["B11"])
        assert stats["B11"]["n_pixels"] == 0
        assert stats["B11"]["std"] == 1.0

    def test_unreadable_scene_skipped(self, tmp_path):
        from sentinel_processor.dl.normalize import compute_dataset_stats

        good = _write_stats_scene(tmp_path, "s_20260101T000000.nc")
        bad = tmp_path / "broken.nc"
        bad.write_bytes(b"not a netcdf file")
        stats = compute_dataset_stats([bad, good], ["B04"])
        assert stats["B04"]["n_pixels"] == 16 * 16

    def test_unsupported_format_skipped(self, tmp_path):
        from sentinel_processor.dl.normalize import compute_dataset_stats

        txt = tmp_path / "scene.txt"
        txt.write_text("nope")
        stats = compute_dataset_stats([txt], ["B04"])
        assert stats["B04"]["n_pixels"] == 0

    def test_custom_nodata_excluded(self, tmp_path):
        from sentinel_processor.dl.normalize import compute_dataset_stats
        pytest.importorskip("netCDF4")
        import xarray as xr

        arr = np.full((8, 8), 1000.0)
        arr[0, :4] = -1.0  # 4 nodata pixels under a custom sentinel
        ds = xr.Dataset({"B04": (("y", "x"), arr)})
        path = tmp_path / "s_20260101T000000.nc"
        ds.to_netcdf(str(path))
        ds.close()

        stats = compute_dataset_stats([path], ["B04"], nodata=-1.0)
        assert stats["B04"]["n_pixels"] == 60
        assert stats["B04"]["mean"] == pytest.approx(1000.0)


class TestStatsPersistence:

    def test_save_load_roundtrip(self, tmp_path):
        stats = {"B04": {"mean": 1.5, "std": 0.5}}
        path = tmp_path / "stats.json"
        save_stats(stats, path)
        assert load_stats(path) == stats

    def test_json_suffix_appended(self, tmp_path):
        save_stats({"a": 1}, tmp_path / "stats")
        assert (tmp_path / "stats.json").exists()

    def test_parent_dirs_created(self, tmp_path):
        save_stats({"a": 1}, tmp_path / "deep" / "dir" / "stats.json")
        assert (tmp_path / "deep" / "dir" / "stats.json").exists()

    def test_load_missing_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_stats(tmp_path / "absent.json")


class TestNormalizeErrors:

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match="Unknown normalisation method"):
            normalize_for_dl(np.ones((2, 4, 4)), method="banana")

    def test_4d_input_raises(self):
        with pytest.raises(ValueError, match="1-D, 2-D, or 3-D"):
            normalize_for_dl(np.ones((1, 2, 4, 4)), method="minmax")

    def test_preset_without_bands_raises(self):
        with pytest.raises(ValueError, match="'bands' must be provided"):
            normalize_for_dl(np.ones((3, 4, 4)), method="sentinel2_rgb")

    def test_preset_band_count_mismatch_raises(self):
        with pytest.raises(ValueError, match="does not match"):
            normalize_for_dl(
                np.ones((2, 4, 4)), method="sentinel2_rgb",
                bands=["B04", "B03", "B02"],
            )

    def test_preset_unknown_band_raises(self):
        with pytest.raises(KeyError, match="not in preset"):
            normalize_for_dl(
                np.ones((3, 4, 4)), method="sentinel2_rgb",
                bands=["B04", "B03", "B99"],
            )
