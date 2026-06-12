from __future__ import annotations
import json
import re
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pytest
from sentinel_processor.input.timeseries import (
    SceneInfo,
    StackResult,
    TimeSeriesConfig,
    _extract_affine,
    _resolve_bands,
    _scl_path_for_scene,
    _report_for_scene,
    _ts_from_path,
    stack_timeseries,
)


# Helpers shared across tests

def _make_da(rows=64, cols=64, n_bands=3, band_names=None, crs="EPSG:32634"):
    """Return a small in-memory DataArray (band, y, x) with spatial metadata."""
    import xarray as xr
    import rioxarray  # noqa: F401

    if band_names is None:
        band_names = [f"B{i:02d}" for i in range(1, n_bands + 1)]

    data = np.random.default_rng(42).random((n_bands, rows, cols)).astype(np.float32)
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


def _write_nc(da, path: Path) -> Path:
    """Save a DataArray to NetCDF and return the path."""
    da.to_netcdf(str(path))
    return path


def _make_scene_info(
        name="budapest_20260520T094746.nc",
        accepted=True,
        cloud=0.02,
        snow=0.0,
        conf=1.0,
        dim_pass=True,
        radio_pass=True,
        reject_reason="",
        report_source="sidecar_json",
):
    return SceneInfo(
        path=Path(name),
        timestamp=datetime(2026, 5, 20, 9, 47, 46, tzinfo=timezone.utc),
        accepted=accepted,
        cloud_fraction=cloud,
        snow_fraction=snow,
        confidence_score=conf,
        dimension_pass=dim_pass,
        radiometry_pass=radio_pass,
        reject_reason=reject_reason,
        report_source=report_source,
    )


# _ts_from_path

class TestTsFromPath:

    def test_standard_filename(self):
        p = Path("budapest_20260520T094746.nc")
        ts = _ts_from_path(p)
        assert ts == datetime(2026, 5, 20, 9, 47, 46, tzinfo=timezone.utc)

    def test_timestamp_in_middle_of_name(self):
        p = Path("loc_name_20251231T235959_extra.tif")
        ts = _ts_from_path(p)
        assert ts.year == 2025
        assert ts.month == 12
        assert ts.day == 31

    def test_no_timestamp_falls_back_to_mtime(self, tmp_path):
        f = tmp_path / "no_ts_file.nc"
        f.write_text("x")
        ts = _ts_from_path(f)
        assert isinstance(ts, datetime)
        assert ts.tzinfo is not None

    def test_returns_utc(self):
        p = Path("x_20260101T000000.nc")
        ts = _ts_from_path(p)
        assert ts.tzinfo == timezone.utc

    def test_sorts_correctly(self):
        names = [
            "budapest_20260528T094727.nc",
            "budapest_20260520T094746.nc",
            "budapest_20260523T095743.nc",
        ]
        paths = [Path(n) for n in names]
        paths.sort(key=_ts_from_path)
        stems = [p.stem for p in paths]
        assert stems[0].startswith("budapest_20260520")
        assert stems[-1].startswith("budapest_20260528")


# _extract_affine

class TestExtractAffine:

    def test_valid_da(self):
        da = _make_da(32, 32)
        result = _extract_affine(da)
        assert result is not None
        ox, oy, pw, ph = result
        assert pw > 0
        assert ph < 0

    def test_single_pixel_returns_none(self):
        import xarray as xr
        da = xr.DataArray(
            np.ones((1, 1, 1)),
            dims=["band", "y", "x"],
            coords={"band": ["B01"], "y": [0.0], "x": [0.0]},
        )
        assert _extract_affine(da) is None

    def test_missing_x_coord_returns_none(self):
        import xarray as xr
        da = xr.DataArray(np.ones((2, 4, 1)), dims=["band", "y", "x"])
        assert _extract_affine(da) is None


# _resolve_bands

class TestResolveBands:

    def test_canonical_names_pass_through(self):
        da = _make_da(band_names=["B02", "B03", "B04"])
        result = _resolve_bands(da)
        assert result == ["B02", "B03", "B04"]

    def test_human_readable_aliases_resolved(self):
        da = _make_da(n_bands=4, band_names=["blue", "green", "red", "nir"])
        result = _resolve_bands(da)
        assert result == ["B02", "B03", "B04", "B08"]

    def test_mixed_aliases(self):
        da = _make_da(band_names=["nir", "swir16", "swir22"])
        result = _resolve_bands(da)
        assert result == ["B08", "B11", "B12"]

    def test_no_band_coord_returns_empty(self):
        import xarray as xr
        da = xr.DataArray(np.ones((4, 4)), dims=["y", "x"])
        assert _resolve_bands(da) == []

    def test_unknown_alias_passes_through_as_is(self):
        da = _make_da(n_bands=2, band_names=["nir", "unknown_band"])
        result = _resolve_bands(da)
        assert "B08" in result
        assert "unknown_band" in result


# _scl_path_for_scene  and  _report_for_scene

class TestSclPathForScene:

    def test_finds_tif_in_scl_dir(self, tmp_path):
        scl = tmp_path / "scl_budapest_20260520T094746.tif"
        scl.write_bytes(b"")
        spectral = Path("data/spectral/budapest_20260520T094746.nc")
        result = _scl_path_for_scene(spectral, tmp_path)
        assert result == scl

    def test_finds_nc_in_scl_dir(self, tmp_path):
        scl = tmp_path / "scl_budapest_20260520T094746.nc"
        scl.write_bytes(b"")
        spectral = Path("data/spectral/budapest_20260520T094746.nc")
        result = _scl_path_for_scene(spectral, tmp_path)
        assert result == scl

    def test_returns_none_when_missing(self, tmp_path):
        spectral = Path("data/spectral/budapest_20260520T094746.nc")
        result = _scl_path_for_scene(spectral, tmp_path)
        assert result is None

    def test_no_scl_dir_tries_same_dir(self, tmp_path):
        scl = tmp_path / "scl_my_scene_20260520T094746.tif"
        scl.write_bytes(b"")
        spectral = tmp_path / "my_scene_20260520T094746.nc"
        result = _scl_path_for_scene(spectral, scl_dir=None)
        assert result == scl


class TestReportForScene:

    def test_finds_report_next_to_spectral(self, tmp_path):
        report = {"passed": True, "confidence_score": 1.0}
        rp = tmp_path / "budapest_20260520T094746_report.json"
        rp.write_text(json.dumps(report))
        spectral = tmp_path / "budapest_20260520T094746.nc"
        result = _report_for_scene(spectral, scl_dir=None)
        assert result is not None
        assert result["passed"] is True

    def test_finds_report_in_scl_dir(self, tmp_path):
        scl_dir = tmp_path / "technical"
        scl_dir.mkdir()
        report = {"confidence_score": 0.75}
        (scl_dir / "budapest_20260520T094746_report.json").write_text(json.dumps(report))
        spectral = tmp_path / "spectral" / "budapest_20260520T094746.nc"
        result = _report_for_scene(spectral, scl_dir=scl_dir)
        assert result is not None
        assert result["confidence_score"] == 0.75

    def test_returns_none_when_no_report(self, tmp_path):
        spectral = tmp_path / "budapest_20260520T094746.nc"
        assert _report_for_scene(spectral, scl_dir=None) is None

    def test_returns_none_on_corrupt_json(self, tmp_path):
        rp = tmp_path / "budapest_20260520T094746_report.json"
        rp.write_text("{not valid json")
        spectral = tmp_path / "budapest_20260520T094746.nc"
        assert _report_for_scene(spectral, scl_dir=None) is None


# SceneInfo

class TestSceneInfo:

    def test_as_dict_keys(self):
        s = _make_scene_info()
        d = s.as_dict()
        expected_keys = {
            "file", "timestamp", "accepted", "cloud_fraction", "snow_fraction",
            "confidence_score", "dimension_pass", "radiometry_pass",
            "reject_reason", "report_source",
        }
        assert set(d.keys()) == expected_keys

    def test_as_dict_timestamp_is_iso_string(self):
        s = _make_scene_info()
        d = s.as_dict()
        assert isinstance(d["timestamp"], str)
        datetime.fromisoformat(d["timestamp"])  # must not raise

    def test_as_dict_file_is_string(self):
        s = _make_scene_info()
        assert isinstance(s.as_dict()["file"], str)


# TimeSeriesConfig

class TestTimeSeriesConfig:

    def test_defaults(self):
        cfg = TimeSeriesConfig()
        assert cfg.align is True
        assert cfg.max_cloud_fraction == pytest.approx(0.30)
        assert cfg.min_confidence == pytest.approx(0.01)
        assert cfg.max_snow_fraction == pytest.approx(1.0)
        assert cfg.fill_rejected is False
        assert cfg.use_sidecar_report is True
        assert cfg.run_radiometry_check is True
        assert cfg.save_dir is None
        assert cfg.save_format == "nc"
        assert cfg.save_name is None

    def test_custom_values(self):
        cfg = TimeSeriesConfig(
            max_cloud_fraction=0.05,
            min_confidence=0.80,
            require_bands=["nir", "red"],
            fill_rejected=True,
            save_dir="/tmp/stacks",
            save_format="tif",
            save_name="test_stack",
        )
        assert cfg.max_cloud_fraction == pytest.approx(0.05)
        assert cfg.require_bands == ["nir", "red"]
        assert cfg.save_format == "tif"
        assert cfg.save_name == "test_stack"


# StackResult  (unit — no I/O)

class TestStackResult:

    def _make_result(self, n_accepted=4, n_rejected=1):
        import xarray as xr
        stack = xr.DataArray(
            np.zeros((n_accepted, 3, 8, 8), dtype=np.float32),
            dims=["time", "band", "y", "x"],
        )
        scenes = [_make_scene_info(accepted=True) for _ in range(n_accepted)]
        scenes += [_make_scene_info(accepted=False, reject_reason="cloud 45.0% > 30.0%")
                   for _ in range(n_rejected)]
        return StackResult(stack=stack, scenes=scenes)

    def test_n_accepted(self):
        r = self._make_result(4, 1)
        assert r.n_accepted == 4

    def test_n_rejected(self):
        r = self._make_result(4, 1)
        assert r.n_rejected == 1

    def test_summary_contains_checkmarks(self):
        r = self._make_result(2, 1)
        s = r.summary()
        assert "✓" in s
        assert "✗" in s

    def test_summary_contains_totals(self):
        r = self._make_result(3, 2)
        s = r.summary()
        assert "3" in s
        assert "5" in s

    def test_to_json(self, tmp_path):
        r = self._make_result(2, 1)
        p = tmp_path / "log.json"
        r.to_json(p)
        data = json.loads(p.read_text())
        assert "scenes" in data
        assert len(data["scenes"]) == 3

    def test_to_json_all_fields_present(self, tmp_path):
        r = self._make_result(1, 0)
        p = tmp_path / "log.json"
        r.to_json(p)
        scene = json.loads(p.read_text())["scenes"][0]
        for key in ("file", "timestamp", "accepted", "cloud_fraction",
                    "snow_fraction", "confidence_score", "dimension_pass",
                    "radiometry_pass", "reject_reason", "report_source"):
            assert key in scene, f"missing key: {key}"


# StackResult.save  (integration — writes real files)

@pytest.mark.integration
class TestStackResultSave:

    def _make_real_result(self, tmp_path, n_times=3):
        """Build a StackResult with actual xarray data that can be saved."""
        import xarray as xr
        import rioxarray  # noqa: F401

        times = np.array(
            ["2026-05-20T09:47:46", "2026-05-23T09:47:26", "2026-05-26T09:57:21"],
            dtype="datetime64[ns]",
        )[:n_times]

        xs = np.linspace(600000.0, 600110.0, 12, endpoint=False)
        ys = np.linspace(5200000.0, 5199880.0, 12, endpoint=False)

        data = np.random.default_rng(0).random(
            (n_times, 3, 12, 12)
        ).astype(np.float32)

        stack = xr.DataArray(
            data,
            dims=["time", "band", "y", "x"],
            coords={
                "time": times,
                "band": ["red", "green", "blue"],
                "y": ys,
                "x": xs,
            },
        )
        stack = stack.rio.write_crs("EPSG:32634")
        scenes = [_make_scene_info() for _ in range(n_times)]
        return StackResult(stack=stack, scenes=scenes)

    def test_save_nc_creates_file(self, tmp_path):
        r = self._make_real_result(tmp_path)
        paths = r.save(tmp_path, name="test_stack", fmt="nc")
        assert len(paths) == 1
        assert Path(paths[0]).exists()
        assert paths[0].endswith(".nc")

    def test_save_nc_filename(self, tmp_path):
        r = self._make_real_result(tmp_path)
        paths = r.save(tmp_path, name="my_location_stack", fmt="nc")
        assert "my_location_stack.nc" in paths[0]

    def test_save_nc_auto_name_derives_from_scenes(self, tmp_path):
        r = self._make_real_result(tmp_path)
        paths = r.save(tmp_path, fmt="nc")
        fname = Path(paths[0]).name
        assert fname.endswith(".nc")
        assert "_stack" in fname

    def test_save_nc_readable(self, tmp_path):
        pytest.importorskip("xarray")
        import xarray as xr
        r = self._make_real_result(tmp_path)
        paths = r.save(tmp_path, name="readable_test", fmt="nc")
        ds = xr.open_dataset(paths[0])
        assert "time" in ds.coords or len(ds.dims) > 0
        ds.close()

    def test_save_nc_time_axis_preserved(self, tmp_path):
        pytest.importorskip("xarray")
        import xarray as xr
        r = self._make_real_result(tmp_path, n_times=3)
        paths = r.save(tmp_path, name="time_test", fmt="nc")
        ds = xr.open_dataset(paths[0])
        # at least one dimension should have length 3
        assert any(v == 3 for v in ds.sizes.values())
        ds.close()

    def test_save_tif_creates_one_per_timestep(self, tmp_path):
        r = self._make_real_result(tmp_path, n_times=3)
        paths = r.save(tmp_path, name="scene", fmt="tif")
        assert len(paths) == 3
        for p in paths:
            assert Path(p).exists()
            assert p.endswith(".tif")

    def test_save_tif_filenames_contain_timestamp(self, tmp_path):
        r = self._make_real_result(tmp_path, n_times=2)
        paths = r.save(tmp_path, name="scene", fmt="tif")
        ts_re = re.compile(r"\d{8}T\d{6}")
        for p in paths:
            assert ts_re.search(Path(p).name), f"no timestamp in {Path(p).name}"

    def test_save_invalid_fmt_raises(self, tmp_path):
        r = self._make_real_result(tmp_path)
        with pytest.raises(ValueError, match="fmt"):
            r.save(tmp_path, fmt="xyz")

    def test_save_creates_output_dir(self, tmp_path):
        r = self._make_real_result(tmp_path)
        new_dir = tmp_path / "nested" / "output"
        r.save(new_dir, fmt="nc")
        assert new_dir.exists()

    def test_save_returns_absolute_paths(self, tmp_path):
        r = self._make_real_result(tmp_path)
        paths = r.save(tmp_path, fmt="nc")
        assert all(Path(p).is_absolute() for p in paths)


# stack_timeseries  (integration — writes and reads real files)

@pytest.mark.integration
class TestStackTimeseries:

    def _make_scene_files(self, tmp_path, n=3, band_names=None, rows=32, cols=32):
        import rasterio
        from rasterio.transform import from_origin
        import rioxarray  # noqa: F401

        if band_names is None:
            band_names = ["blue", "green", "red", "nir"]

        timestamps = [
            "20260520T094746",
            "20260523T094726",
            "20260526T095721",
            "20260528T094727",
            "20260530T095800",
        ][:n]

        tmp_path.mkdir(parents=True, exist_ok=True)
        transform = from_origin(600000.0, 5200000.0, 10.0, 10.0)
        paths = []

        for ts in timestamps:
            data = np.random.default_rng(int(ts[:8])).random(
                (len(band_names), rows, cols)
            ).astype(np.float32)
            fpath = tmp_path / f"budapest_{ts}.tif"
            with rasterio.open(
                    str(fpath), "w",
                    driver="GTiff",
                    height=rows, width=cols,
                    count=len(band_names),
                    dtype="float32",
                    crs="EPSG:32634",
                    transform=transform,
            ) as dst:
                dst.write(data)
                dst.update_tags(**{f"band_{i + 1}": name
                                   for i, name in enumerate(band_names)})
            paths.append(fpath)

        return paths

    def _write_reports(self, tmp_path, paths, cloud=0.02, conf=1.0):
        for p in paths:
            report = {
                "passed": True,
                "cloud_ratio": cloud,
                "snow_ratio": 0.0,
                "confidence_score": conf,
                "dimension_pass": True,
                "radiometry_pass": True,
                "issues": [],
            }
            rp = tmp_path / f"{p.stem}_report.json"
            rp.write_text(json.dumps(report))

    # basic functionality

    def test_returns_stack_result(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, scl_dir=None, progress=False)
        assert isinstance(result, StackResult)

    def test_stack_shape_time_first(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=3)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, progress=False)
        assert result.stack.dims[0] == "time"
        assert result.stack.shape[0] == 3

    def test_stack_dtype_float32(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, progress=False)
        assert result.stack.dtype == np.float32

    def test_time_coordinate_sorted(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=3)
        self._write_reports(tmp_path, paths)
        # pass in reverse order
        result = stack_timeseries(list(reversed(paths)), progress=False)
        times = result.stack.coords["time"].values
        assert list(times) == sorted(times)

    def test_band_coord_preserved(self, tmp_path):
        # tif files have integer band coords
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, progress=False)
        assert "band" in result.stack.dims

    def test_n_accepted_matches_all_clear_scenes(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=4)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, progress=False)
        assert result.n_accepted == 4
        assert result.n_rejected == 0

    # quality filtering

    def test_cloudy_scenes_rejected(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=3)
        # first two clean, third cloudy
        self._write_reports(tmp_path, paths[:2], cloud=0.02, conf=1.0)
        self._write_reports(tmp_path, paths[2:], cloud=0.60, conf=0.0)

        cfg = TimeSeriesConfig(max_cloud_fraction=0.10, min_confidence=0.01)
        result = stack_timeseries(paths, cfg=cfg, progress=False)
        assert result.n_accepted == 2
        assert result.n_rejected == 1

    def test_low_confidence_rejected(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=3)
        self._write_reports(tmp_path, paths[:2], conf=1.0)
        self._write_reports(tmp_path, paths[2:], cloud=0.25, conf=0.75)

        cfg = TimeSeriesConfig(min_confidence=0.80)
        result = stack_timeseries(paths, cfg=cfg, progress=False)
        assert result.n_rejected >= 1

    def test_reject_reason_populated(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths[:1], cloud=0.02, conf=1.0)
        self._write_reports(tmp_path, paths[1:], cloud=0.70, conf=0.0)

        cfg = TimeSeriesConfig(max_cloud_fraction=0.10)
        result = stack_timeseries(paths, cfg=cfg, progress=False)
        rejected = [s for s in result.scenes if not s.accepted]
        assert len(rejected) == 1
        assert "cloud" in rejected[0].reject_reason.lower()

    def test_all_rejected_raises_value_error(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths, cloud=0.80, conf=0.0)
        cfg = TimeSeriesConfig(max_cloud_fraction=0.05, min_confidence=0.50)
        with pytest.raises(ValueError, match="[Nn]o scenes"):
            stack_timeseries(paths, cfg=cfg, progress=False)

    def test_fill_rejected_keeps_time_axis_contiguous(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=3)
        self._write_reports(tmp_path, paths[:2], cloud=0.02, conf=1.0)
        self._write_reports(tmp_path, paths[2:], cloud=0.80, conf=0.0)

        cfg = TimeSeriesConfig(
            max_cloud_fraction=0.10,
            fill_rejected=True,
            nodata=float("nan"),
        )
        result = stack_timeseries(paths, cfg=cfg, progress=False)
        assert result.stack.shape[0] == 3  # all 3 time steps present
        assert result.n_rejected == 1

    def _make_nc_scenes(self, tmp_path, n=2, band_names=None):
        """Write named-band .nc files for require_bands tests.
        Separate helper because _make_scene_files uses .tif (no named coords).
        """
        import xarray as xr
        import rioxarray  # noqa: F401
        if band_names is None:
            band_names = ["red", "nir"]
        timestamps = ["20260520T094746", "20260523T094726"][:n]
        tmp_path.mkdir(parents=True, exist_ok=True)
        xs = np.linspace(600000.0, 600320.0, 32, endpoint=False)
        ys = np.linspace(5200000.0, 5199680.0, 32, endpoint=False)
        paths = []
        for ts in timestamps:
            data = np.random.default_rng(0).random(
                (len(band_names), 32, 32)).astype(np.float32)
            da = xr.DataArray(data, dims=["band", "y", "x"],
                              coords={"band": band_names, "y": ys, "x": xs})
            da = da.rio.write_crs("EPSG:32634")
            fpath = tmp_path / f"budapest_{ts}.nc"
            da.to_netcdf(str(fpath))
            paths.append(fpath)
        return paths

    def test_require_bands_rejects_missing(self, tmp_path):
        nc_dir = tmp_path / "nc_scenes"
        paths = self._make_nc_scenes(nc_dir, n=2, band_names=["red", "nir"])
        self._write_reports(nc_dir, paths)

        cfg = TimeSeriesConfig(require_bands=["red", "nir", "swir16"])
        with pytest.raises(ValueError, match="[Nn]o scenes"):
            stack_timeseries(paths, cfg=cfg, progress=False)

    def test_require_bands_alias_resolves(self, tmp_path):
        # files have "B08" (canonical) but we require "nir" (alias)
        nc_dir = tmp_path / "nc_scenes"
        paths = self._make_nc_scenes(nc_dir, n=2, band_names=["B04", "B08"])
        self._write_reports(nc_dir, paths)

        cfg = TimeSeriesConfig(require_bands=["red", "nir"])
        result = stack_timeseries(paths, cfg=cfg, progress=False)
        assert result.n_accepted == 2

    def test_use_sidecar_report_true_reads_json(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths, cloud=0.01, conf=1.0)
        cfg = TimeSeriesConfig(use_sidecar_report=True)
        result = stack_timeseries(paths, cfg=cfg, progress=False)
        assert all(s.report_source == "sidecar_json" for s in result.scenes)

    def test_use_sidecar_false_skips_json(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths)
        # no SCL files -> report_source will be "none"
        cfg = TimeSeriesConfig(use_sidecar_report=False)
        result = stack_timeseries(paths, cfg=cfg, progress=False)
        assert all(s.report_source == "none" for s in result.scenes)

    # empty source

    def test_empty_sources_raises(self):
        with pytest.raises(ValueError, match="empty"):
            stack_timeseries([], progress=False)

    # auto-save via config

    def test_save_dir_in_config_writes_file(self, tmp_path):
        scene_dir = tmp_path / "scenes"
        paths = self._make_scene_files(scene_dir, n=2)
        self._write_reports(scene_dir, paths)
        out_dir = tmp_path / "stacks"

        cfg = TimeSeriesConfig(save_dir=out_dir, save_format="tif", save_name="test")
        stack_timeseries(paths, cfg=cfg, progress=False)
        written = list(out_dir.glob("test_*.tif"))
        assert len(written) == 2

    def test_save_dir_none_does_not_write(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths)
        cfg = TimeSeriesConfig(save_dir=None)
        stack_timeseries(paths, cfg=cfg, progress=False)
        assert not (tmp_path / "stacks").exists()

    # result metadata

    def test_scenes_list_length_equals_sources(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=4)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, progress=False)
        assert len(result.scenes) == 4

    def test_scene_info_has_correct_timestamps(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, progress=False)
        for scene in result.scenes:
            assert scene.timestamp.year == 2026

    def test_summary_runs_without_error(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, progress=False)
        s = result.summary()
        assert isinstance(s, str)
        assert len(s) > 0

    def test_to_json_written_correctly(self, tmp_path):
        paths = self._make_scene_files(tmp_path, n=2)
        self._write_reports(tmp_path, paths)
        result = stack_timeseries(paths, progress=False)
        log_path = tmp_path / "log.json"
        result.to_json(log_path)
        data = json.loads(log_path.read_text())
        assert len(data["scenes"]) == 2


# _align_to_reference / _scl_metrics


class TestAlignToReference:

    def test_same_shape_returned_unchanged(self):
        from sentinel_processor.input.timeseries import _align_to_reference
        ref = _make_da(rows=32, cols=32)
        da = _make_da(rows=32, cols=32)
        assert _align_to_reference(da, ref) is da

    def test_same_crs_resampled_to_reference_grid(self):
        from sentinel_processor.input.timeseries import _align_to_reference
        ref = _make_da(rows=64, cols=64)
        da = _make_da(rows=32, cols=32)
        out = _align_to_reference(da, ref)
        assert out.shape[-2:] == (64, 64)
        np.testing.assert_array_equal(out.coords["x"].values,
                                      ref.coords["x"].values)


class TestSclMetrics:

    def _write_scl_tif(self, tmp_path, values):
        pytest.importorskip("rasterio")
        import rasterio
        from rasterio.transform import from_bounds

        rows, cols = values.shape
        path = tmp_path / "scl_scene_20260101T000000.tif"
        with rasterio.open(
                str(path), "w",
                driver="GTiff", height=rows, width=cols, count=1,
                dtype="int32", crs="EPSG:32634",
                transform=from_bounds(0, 0, cols, rows, cols, rows),
        ) as dst:
            dst.write(values.astype(np.int32), 1)
        return path

    def test_clear_scene_metrics(self, tmp_path):
        from sentinel_processor.input.timeseries import _scl_metrics
        scl = np.full((64, 64), 4)  # all vegetation
        path = self._write_scl_tif(tmp_path, scl)
        m = _scl_metrics(path, max_cloud_threshold=0.3)
        assert m["dimension_pass"] is True
        assert m["cloud_fraction"] == pytest.approx(0.0)
        assert m["confidence_score"] == pytest.approx(1.0)

    def test_cloudy_scene_metrics(self, tmp_path):
        from sentinel_processor.input.timeseries import _scl_metrics
        scl = np.full((64, 64), 9)  # all high-probability cloud
        path = self._write_scl_tif(tmp_path, scl)
        m = _scl_metrics(path, max_cloud_threshold=0.3)
        assert m["cloud_fraction"] == pytest.approx(1.0)
        assert m["confidence_score"] == pytest.approx(0.0)

    def test_tiny_scene_fails_dimension_check(self, tmp_path):
        from sentinel_processor.input.timeseries import _scl_metrics
        scl = np.full((8, 8), 4)  # below MIN_SIDE=32
        path = self._write_scl_tif(tmp_path, scl)
        m = _scl_metrics(path, max_cloud_threshold=0.3)
        assert m["dimension_pass"] is False
        assert m["confidence_score"] == 0.0

    def test_unreadable_scl_passes_with_warning(self, tmp_path):
        from sentinel_processor.input.timeseries import _scl_metrics
        m = _scl_metrics(tmp_path / "missing.tif", max_cloud_threshold=0.3)
        assert m["confidence_score"] == 1.0
        assert "warning" in m
