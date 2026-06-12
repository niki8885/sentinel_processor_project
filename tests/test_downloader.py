from __future__ import annotations
import datetime
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from sentinel_processor.input.downloader import (
    DownloadConfig,
    _apply_pansharpening,
    _extract_affine,
    _progress_done,
    _progress_start,
    _progress_update,
    _run_validation,
    _save_both,
    _write_report,
)


def _make_da(rows=64, cols=64, n_bands=3, band_names=None, crs="EPSG:32634"):
    """Small in-memory DataArray (band, y, x) with spatial metadata."""
    import xarray as xr
    import rioxarray  # noqa: F401

    if band_names is None:
        band_names = [f"B{i:02d}" for i in range(1, n_bands + 1)]

    data = np.random.default_rng(7).random((n_bands, rows, cols)).astype(np.float32)
    xs = np.linspace(600000.0, 600000.0 + cols * 10.0, cols, endpoint=False)
    ys = np.linspace(5200000.0, 5200000.0 - rows * 10.0, rows, endpoint=False)

    da = xr.DataArray(
        data,
        dims=["band", "y", "x"],
        coords={"band": band_names, "y": ys, "x": xs},
    )
    if crs:
        da = da.rio.write_crs(crs)
    return da


# DownloadConfig


class TestDownloadConfig:

    def test_resolved_end_defaults_to_now(self):
        cfg = DownloadConfig()
        now = datetime.datetime.now(datetime.UTC)
        assert abs((cfg.resolved_end() - now).total_seconds()) < 5

    def test_resolved_end_explicit(self):
        end = datetime.datetime(2026, 1, 1, tzinfo=datetime.UTC)
        cfg = DownloadConfig(end_date=end)
        assert cfg.resolved_end() == end

    def test_resolved_start_from_lookback(self):
        end = datetime.datetime(2026, 1, 31, tzinfo=datetime.UTC)
        cfg = DownloadConfig(end_date=end, lookback_days=30)
        assert cfg.resolved_start() == end - datetime.timedelta(days=30)

    def test_resolved_start_explicit(self):
        start = datetime.datetime(2025, 6, 1, tzinfo=datetime.UTC)
        cfg = DownloadConfig(start_date=start)
        assert cfg.resolved_start() == start

    def test_date_range_str_format(self):
        start = datetime.datetime(2026, 1, 1, tzinfo=datetime.UTC)
        end = datetime.datetime(2026, 2, 1, tzinfo=datetime.UTC)
        cfg = DownloadConfig(start_date=start, end_date=end)
        assert cfg.date_range_str() == "2026-01-01T00:00:00Z/2026-02-01T00:00:00Z"

    def test_band_keys_from_list(self):
        cfg = DownloadConfig(bands=["red", "nir"])
        assert cfg.band_keys() == ["red", "nir"]

    def test_tech_keys_none_is_empty(self):
        cfg = DownloadConfig(tech_bands=None)
        assert cfg.tech_keys() == []

    def test_tech_keys_from_list(self):
        cfg = DownloadConfig(tech_bands=["scl"])
        assert cfg.tech_keys() == ["scl"]

    def test_subdir_creates_directory(self, tmp_path):
        cfg = DownloadConfig(output_dir=str(tmp_path / "out"))
        sub = cfg.subdir("spectral")
        assert Path(sub).is_dir()
        assert Path(sub).name == "spectral"


# _extract_affine


class TestExtractAffine:

    def test_returns_origin_and_pixel_size(self):
        da = _make_da(rows=32, cols=32)
        affine = _extract_affine(da)
        assert affine is not None
        ox, oy, pw, ph = affine
        assert pw == pytest.approx(10.0)
        assert ph == pytest.approx(-10.0)
        assert ox == pytest.approx(600000.0)
        assert oy == pytest.approx(5200000.0)

    def test_too_few_coords_returns_none(self):
        import xarray as xr
        da = xr.DataArray(
            np.zeros((1, 1)),
            dims=["y", "x"],
            coords={"y": [0.0], "x": [0.0]},
        )
        assert _extract_affine(da) is None

    def test_degenerate_pixel_size_returns_none(self):
        import xarray as xr
        da = xr.DataArray(
            np.zeros((4, 4)),
            dims=["y", "x"],
            coords={"y": [0.0, 0.0, 0.0, 0.0], "x": [0.0, 0.0, 0.0, 0.0]},
        )
        assert _extract_affine(da) is None


# _save_both / _write_report


class TestSaveBoth:

    def test_writes_tif_and_nc(self, tmp_path):
        da = _make_da(rows=16, cols=16, n_bands=2)
        written = _save_both(da, str(tmp_path), "scene")
        assert (tmp_path / "scene.tif").exists()
        assert str(tmp_path / "scene.tif") in written
        # NetCDF is best-effort (backend may be missing)
        if len(written) == 2:
            assert (tmp_path / "scene.nc").exists()


class TestWriteReport:

    def test_writes_json_report(self, tmp_path):
        report = {"item_id": "abc", "passed": True}
        _write_report(report, str(tmp_path), "scene01")
        path = tmp_path / "scene01_report.json"
        assert path.exists()
        assert json.loads(path.read_text())["item_id"] == "abc"


# _run_validation (mocked STAC items)


def _has_validation_lib() -> bool:
    try:
        from sentinel_processor.validation._fortran_bridge import _get_lib
        _get_lib()
        return True
    except Exception:
        return False


skip_no_lib = pytest.mark.skipif(
    not _has_validation_lib(), reason="Fortran validation library not built"
)


class TestRunValidation:

    def test_missing_scl_asset_passes_with_warning(self):
        item = SimpleNamespace(id="item-1", assets={}, properties={})
        cfg = DownloadConfig()
        passes, report, scl = _run_validation(item, 19.0, 47.5, 0.05, None, cfg)
        assert passes is True
        assert scl is None
        assert "warning" in report

    @skip_no_lib
    def test_dimension_precheck_rejects_tiny_scene(self):
        scl_asset = SimpleNamespace(
            href="https://example.invalid/scl.tif",
            extra_fields={"proj:shape": [10, 10]},  # below MIN_SIDE=32
        )
        item = SimpleNamespace(id="item-2", assets={"scl": scl_asset}, properties={})
        cfg = DownloadConfig()
        passes, report, scl = _run_validation(item, 19.0, 47.5, 0.05, None, cfg)
        assert passes is False
        assert report["dimension_pass"] is False
        assert report["issues"]

    @skip_no_lib
    def test_unreadable_scl_passes_with_error(self):
        scl_asset = SimpleNamespace(
            href=str(Path("definitely") / "missing" / "scl.tif"),
            extra_fields={},
        )
        item = SimpleNamespace(id="item-3", assets={"scl": scl_asset}, properties={})
        cfg = DownloadConfig()
        passes, report, scl = _run_validation(item, 19.0, 47.5, 0.05, None, cfg)
        assert passes is True
        assert "error" in report


# _apply_pansharpening — degraded paths


class TestApplyPansharpening:

    def test_missing_pan_asset_returns_input(self):
        ds = _make_da(rows=16, cols=16, n_bands=3)
        item = SimpleNamespace(id="item-4", assets={})
        cfg = DownloadConfig()
        out = _apply_pansharpening(ds, item, 19.0, 47.5, 0.05, cfg)
        assert out is ds

    def test_unreadable_pan_returns_input(self):
        ds = _make_da(rows=16, cols=16, n_bands=3)
        pan_asset = SimpleNamespace(href=str(Path("missing") / "pan.tif"))
        item = SimpleNamespace(id="item-5", assets={"visual": pan_asset})
        cfg = DownloadConfig()
        out = _apply_pansharpening(ds, item, 19.0, 47.5, 0.05, cfg)
        assert out is ds


# progress helpers


class TestProgress:

    def test_progress_start_prints(self, capsys):
        _progress_start(3, enabled=True)
        assert "3 scene(s)" in capsys.readouterr().out

    def test_progress_start_disabled_silent(self, capsys):
        _progress_start(3, enabled=False)
        assert capsys.readouterr().out == ""

    def test_progress_update_prints_bar(self, capsys):
        _progress_update(1, 2, "scene_a", enabled=True)
        out = capsys.readouterr().out
        assert "50%" in out
        assert "scene_a" in out

    def test_progress_update_truncates_long_names(self, capsys):
        _progress_update(1, 1, "x" * 50, enabled=True)
        assert "..." in capsys.readouterr().out

    def test_progress_update_disabled_silent(self, capsys):
        _progress_update(1, 2, "scene", enabled=False)
        assert capsys.readouterr().out == ""

    def test_progress_done_prints_summary(self, capsys):
        _progress_done(2, {"a": ["f1", "f2"], "b": ["f3"]}, enabled=True)
        out = capsys.readouterr().out
        assert "2 scene(s)" in out
        assert "3 file(s)" in out

    def test_progress_done_disabled_silent(self, capsys):
        _progress_done(2, {}, enabled=False)
        assert capsys.readouterr().out == ""
