import numpy as np
import pytest

plotly = pytest.importorskip("plotly")
go = pytest.importorskip("plotly.graph_objects")

from sentinel_processor.visualisation.plot import (  # noqa: E402
    _extract_band,
    _load_da,
    _pstretch,
    DEFAULT_BAD_CLASSES,
    SCL_CLASSES,
    plot_band,
    plot_grid,
    plot_mask,
    plot_rgb,
    plot_timeseries,
)


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture()
def single_band_tif(tmp_path):
    """32×32 single-band GeoTIFF."""
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_bounds
    rng = np.random.default_rng(0)
    data = rng.random((1, 32, 32)).astype(np.float32)
    tf = from_bounds(0, 0, 1, 1, 32, 32)
    p = tmp_path / "band.tif"
    with rasterio.open(str(p), "w", driver="GTiff", height=32, width=32,
                       count=1, dtype="float32", crs="EPSG:4326",
                       transform=tf) as dst:
        dst.write(data)
    return p


@pytest.fixture()
def three_band_tif(tmp_path):
    """32×32 three-band GeoTIFF (R, G, B)."""
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_bounds
    rng = np.random.default_rng(1)
    data = rng.random((3, 32, 32)).astype(np.float32)
    tf = from_bounds(0, 0, 1, 1, 32, 32)
    p = tmp_path / "rgb.tif"
    with rasterio.open(str(p), "w", driver="GTiff", height=32, width=32,
                       count=3, dtype="float32", crs="EPSG:4326",
                       transform=tf) as dst:
        dst.write(data)
    return p


@pytest.fixture()
def scl_tif(tmp_path):
    """32×32 SCL GeoTIFF: mix of vegetation (4), cloud (8), nodata (0)."""
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_bounds
    data = np.full((1, 32, 32), 4, dtype=np.uint8)
    data[0, :8, :] = 8   # cloud rows
    data[0, :2, :] = 0   # nodata rows
    tf = from_bounds(0, 0, 1, 1, 32, 32)
    p = tmp_path / "scl.tif"
    with rasterio.open(str(p), "w", driver="GTiff", height=32, width=32,
                       count=1, dtype="uint8", crs="EPSG:4326",
                       transform=tf) as dst:
        dst.write(data)
    return p


@pytest.fixture()
def timeseries_stack(tmp_path):
    """Small (time=6, band=2, y=16, x=16) xarray DataArray with NDVI+EVI.

    Geographic extent: lon 19.0–19.1 / lat 47.4–47.5 (Budapest area).
    Time steps: monthly from 2024-04 to 2024-09 (all clear by default).
    """
    xr = pytest.importorskip("xarray")
    rng = np.random.default_rng(42)

    times = np.array(
        ["2024-04-01", "2024-05-01", "2024-06-01",
         "2024-07-01", "2024-08-01", "2024-09-01"],
        dtype="datetime64[ns]",
    )
    lons = np.linspace(19.0, 19.1, 16)
    lats = np.linspace(47.4, 47.5, 16)
    data = rng.uniform(0.1, 0.8, (6, 2, 16, 16)).astype(np.float32)

    da = xr.DataArray(
        data,
        dims=["time", "band", "y", "x"],
        coords={
            "time": times,
            "band": ["ndvi", "evi"],
            "y": lats,
            "x": lons,
        },
    )
    return da


@pytest.fixture()
def scl_dir_with_files(tmp_path, timeseries_stack):
    """Directory of 6 SCL files, one per time step in timeseries_stack.

    Steps 0 and 2 are cloudy (SCL=8 at centre pixel); rest are clear (SCL=4).
    """
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_bounds

    scl_dir = tmp_path / "technical"
    scl_dir.mkdir()

    timestamps = [
        "20240401T000000", "20240501T000000", "20240601T000000",
        "20240701T000000", "20240801T000000", "20240901T000000",
    ]
    cloudy_steps = {0, 2}   # steps that will carry SCL=8

    for i, ts in enumerate(timestamps):
        scl_val = np.uint8(8 if i in cloudy_steps else 4)
        data = np.full((1, 16, 16), scl_val, dtype=np.uint8)
        tf = from_bounds(19.0, 47.4, 19.1, 47.5, 16, 16)
        p = scl_dir / f"scl_budapest_{ts}.tif"
        with rasterio.open(str(p), "w", driver="GTiff", height=16, width=16,
                           count=1, dtype="uint8", crs="EPSG:4326",
                           transform=tf) as dst:
            dst.write(data)

    return scl_dir, cloudy_steps


# ── _pstretch ─────────────────────────────────────────────────────────────────

class TestPstretch:

    def test_output_range_0_to_1(self):
        arr = np.linspace(-100, 100, 1000)
        out = _pstretch(arr)
        assert out.min() >= 0.0
        assert out.max() <= 1.0

    def test_dtype_float32(self):
        arr = np.random.default_rng(0).random(100)
        out = _pstretch(arr)
        assert out.dtype == np.float32

    def test_all_same_returns_zeros(self):
        arr = np.full(100, 5.0)
        out = _pstretch(arr)
        np.testing.assert_allclose(out, 0.0, atol=1e-6)

    def test_custom_percentiles(self):
        arr = np.linspace(0, 100, 1000)
        out = _pstretch(arr, lo=0, hi=100)
        assert out.max() == pytest.approx(1.0, abs=1e-3)


# ── _load_da ──────────────────────────────────────────────────────────────────

class TestLoadDa:

    def test_loads_tif(self, single_band_tif):
        da = _load_da(single_band_tif)
        assert da.values is not None

    def test_missing_file_raises(self):
        with pytest.raises(Exception):
            _load_da("/nonexistent/path.tif")


# ── _extract_band ─────────────────────────────────────────────────────────────

class TestExtractBand:

    def test_none_squeezes(self, single_band_tif):
        da = _load_da(single_band_tif)
        arr = _extract_band(da, None)
        assert arr.ndim == 2

    def test_integer_index(self, three_band_tif):
        da = _load_da(three_band_tif)
        arr = _extract_band(da, 0)
        assert arr.ndim == 2

    def test_out_of_range_raises(self, single_band_tif):
        da = _load_da(single_band_tif)
        with pytest.raises(Exception):
            _extract_band(da, 99)


# ── plot_band ─────────────────────────────────────────────────────────────────

class TestPlotBand:

    def test_returns_figure(self, single_band_tif):
        fig = plot_band(str(single_band_tif))
        assert isinstance(fig, go.Figure)

    def test_has_one_trace(self, single_band_tif):
        fig = plot_band(str(single_band_tif))
        assert len(fig.data) == 1

    def test_trace_is_heatmap(self, single_band_tif):
        fig = plot_band(str(single_band_tif))
        assert isinstance(fig.data[0], go.Heatmap)

    def test_custom_colorscale(self, single_band_tif):
        fig = plot_band(str(single_band_tif), colorscale="Plasma")
        cs = fig.data[0].colorscale
        viridis_start = "#440154"
        if isinstance(cs, str):
            assert cs.lower() == "plasma"
        else:
            assert cs[0][1].lower() != viridis_start.lower()

    def test_custom_title(self, single_band_tif):
        fig = plot_band(str(single_band_tif), title="My Test Band")
        assert "My Test Band" in fig.layout.title.text

    def test_save_html(self, single_band_tif, tmp_path):
        out = tmp_path / "band.html"
        plot_band(str(single_band_tif), save_html=str(out))
        assert out.exists()
        assert out.stat().st_size > 0

    def test_save_creates_parent_dirs(self, single_band_tif, tmp_path):
        out = tmp_path / "nested" / "deep" / "band.html"
        plot_band(str(single_band_tif), save_html=str(out))
        assert out.exists()

    def test_band_int_selector(self, three_band_tif):
        fig = plot_band(str(three_band_tif), band=0)
        assert isinstance(fig, go.Figure)


# ── plot_rgb ──────────────────────────────────────────────────────────────────

class TestPlotRgb:

    def test_mode_a_returns_figure(self, three_band_tif):
        fig = plot_rgb(str(three_band_tif))
        assert isinstance(fig, go.Figure)

    def test_mode_a_trace_is_image(self, three_band_tif):
        fig = plot_rgb(str(three_band_tif))
        assert isinstance(fig.data[0], go.Image)

    def test_mode_b_explicit_bands(self, three_band_tif):
        fig = plot_rgb(str(three_band_tif), red_band=0, green_band=1, blue_band=2)
        assert isinstance(fig, go.Figure)

    def test_custom_title(self, three_band_tif):
        fig = plot_rgb(str(three_band_tif), title="RGB Test")
        assert "RGB Test" in fig.layout.title.text

    def test_save_html(self, three_band_tif, tmp_path):
        out = tmp_path / "rgb.html"
        plot_rgb(str(three_band_tif), save_html=str(out))
        assert out.exists()


# ── plot_grid ─────────────────────────────────────────────────────────────────

class TestPlotGrid:

    def test_returns_figure(self, single_band_tif):
        panels = [{"file": str(single_band_tif), "label": "Band 1"}]
        fig = plot_grid(panels)
        assert isinstance(fig, go.Figure)

    def test_trace_count_equals_panels(self, single_band_tif, three_band_tif):
        panels = [
            {"file": str(single_band_tif), "label": "A"},
            {"file": str(single_band_tif), "label": "B"},
            {"file": str(three_band_tif), "band": 0, "label": "C"},
        ]
        fig = plot_grid(panels, ncols=3)
        assert len(fig.data) == 3

    def test_empty_panels_raises(self):
        with pytest.raises(ValueError, match="empty"):
            plot_grid([])

    def test_ncols_capped_at_n(self, single_band_tif):
        panels = [{"file": str(single_band_tif)}]
        fig = plot_grid(panels, ncols=10)
        assert isinstance(fig, go.Figure)

    def test_per_panel_colorscale(self, single_band_tif):
        panels_hot = [{"file": str(single_band_tif), "colorscale": "Hot"}]
        panels_default = [{"file": str(single_band_tif)}]
        fig_hot = plot_grid(panels_hot)
        fig_default = plot_grid(panels_default)
        assert fig_hot.data[0].colorscale != fig_default.data[0].colorscale

    def test_same_file_cached(self, single_band_tif):
        panels = [
            {"file": str(single_band_tif), "label": "A"},
            {"file": str(single_band_tif), "label": "B"},
        ]
        fig = plot_grid(panels)
        assert len(fig.data) == 2

    def test_save_html(self, single_band_tif, tmp_path):
        out = tmp_path / "grid.html"
        plot_grid([{"file": str(single_band_tif)}], save_html=str(out))
        assert out.exists()


# ── plot_mask ─────────────────────────────────────────────────────────────────

class TestPlotMask:

    def test_returns_figure(self, scl_tif):
        fig = plot_mask(str(scl_tif))
        assert isinstance(fig, go.Figure)

    def test_three_traces(self, scl_tif):
        fig = plot_mask(str(scl_tif))
        assert len(fig.data) == 3

    def test_first_trace_is_heatmap(self, scl_tif):
        assert isinstance(plot_mask(str(scl_tif)).data[0], go.Heatmap)

    def test_second_trace_is_pie(self, scl_tif):
        assert isinstance(plot_mask(str(scl_tif)).data[1], go.Pie)

    def test_third_trace_is_bar(self, scl_tif):
        assert isinstance(plot_mask(str(scl_tif)).data[2], go.Bar)

    def test_return_mask_flag(self, scl_tif):
        result = plot_mask(str(scl_tif), return_mask=True)
        assert isinstance(result, tuple)
        fig, mask = result
        assert isinstance(fig, go.Figure)
        assert isinstance(mask, np.ndarray)

    def test_mask_shape_matches_scl(self, scl_tif):
        _, mask = plot_mask(str(scl_tif), return_mask=True)
        assert mask.shape == (32, 32)

    def test_mask_is_binary(self, scl_tif):
        _, mask = plot_mask(str(scl_tif), return_mask=True)
        assert set(np.unique(mask)).issubset({0, 1})

    def test_default_bad_classes_applied(self, scl_tif):
        """Rows 2-7 have SCL=8 (cloud) → must be flagged bad."""
        _, mask = plot_mask(str(scl_tif), return_mask=True)
        assert mask[2:8, :].all()

    def test_clear_rows_not_bad(self, scl_tif):
        _, mask = plot_mask(str(scl_tif), return_mask=True)
        assert mask[10:, :].sum() == 0

    def test_custom_bad_classes(self, scl_tif):
        _, mask = plot_mask(str(scl_tif), bad_classes=[0], return_mask=True)
        assert mask[0, 0] == 1
        assert mask[10, 0] == 0

    def test_save_html(self, scl_tif, tmp_path):
        out = tmp_path / "mask.html"
        plot_mask(str(scl_tif), save_html=str(out))
        assert out.exists()

    def test_custom_title(self, scl_tif):
        fig = plot_mask(str(scl_tif), title="Cloud Check")
        assert "Cloud Check" in fig.layout.title.text


# ── plot_timeseries ───────────────────────────────────────────────────────────

class TestPlotTimeseries:

    # ── return type & basic shape ─────────────────────────────────────────────

    def test_returns_figure(self, timeseries_stack):
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45)
        assert isinstance(fig, go.Figure)

    def test_has_traces(self, timeseries_stack):
        """At least one Scatter trace per band."""
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45)
        scatter = [t for t in fig.data if isinstance(t, go.Scatter)]
        assert len(scatter) >= 2   # at least one per band

    def test_line_traces_equal_band_count(self, timeseries_stack):
        """Each band contributes exactly one continuous line trace."""
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi", "evi"])
        lines = [t for t in fig.data if isinstance(t, go.Scatter) and t.mode == "lines"]
        assert len(lines) == 2

    # ── band selection ────────────────────────────────────────────────────────

    def test_single_band_selection(self, timeseries_stack):
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"])
        lines = [t for t in fig.data if isinstance(t, go.Scatter) and t.mode == "lines"]
        assert len(lines) == 1
        assert lines[0].name == "ndvi"

    def test_unknown_band_raises(self, timeseries_stack):
        with pytest.raises(ValueError, match="None of the requested bands"):
            plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["fake_band"])

    def test_default_bands_uses_all(self, timeseries_stack):
        """bands=None → all bands in stack plotted."""
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=None)
        lines = [t for t in fig.data if isinstance(t, go.Scatter) and t.mode == "lines"]
        assert len(lines) == 2

    # ── data values ───────────────────────────────────────────────────────────

    def test_y_values_finite_or_nan(self, timeseries_stack):
        """Line trace y-values should be finite floats or NaN (no inf)."""
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"])
        line = next(t for t in fig.data if isinstance(t, go.Scatter) and t.mode == "lines")
        y = np.array(line.y, dtype=np.float64)
        assert not np.any(np.isinf(y))

    def test_x_length_matches_time_steps(self, timeseries_stack):
        """Line trace x-axis length == number of time steps in the stack."""
        n_times = timeseries_stack.sizes["time"]
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"])
        line = next(t for t in fig.data if isinstance(t, go.Scatter) and t.mode == "lines")
        assert len(line.x) == n_times

    # ── extraction modes ──────────────────────────────────────────────────────

    def test_region_mode_returns_figure(self, timeseries_stack):
        fig = plot_timeseries(
            timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"], agg_bbox=0.01
        )
        assert isinstance(fig, go.Figure)

    def test_region_mode_same_time_length(self, timeseries_stack):
        n_times = timeseries_stack.sizes["time"]
        fig = plot_timeseries(
            timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"], agg_bbox=0.01
        )
        line = next(t for t in fig.data if isinstance(t, go.Scatter) and t.mode == "lines")
        assert len(line.x) == n_times

    # ── SCL cloud markers ─────────────────────────────────────────────────────

    def test_no_scl_no_cloud_markers(self, timeseries_stack):
        """Without scl_path there must be no red ✕ cloud trace."""
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"])
        cloud_traces = [
            t for t in fig.data
            if isinstance(t, go.Scatter)
            and getattr(t.marker, "symbol", None) == "x"
        ]
        assert len(cloud_traces) == 0

    def test_scl_dir_adds_cloud_markers(self, timeseries_stack, scl_dir_with_files):
        scl_dir, cloudy_steps = scl_dir_with_files
        fig = plot_timeseries(
            timeseries_stack,
            lon=19.05, lat=47.45,
            bands=["ndvi"],
            scl_path=scl_dir,
        )
        cloud_traces = [
            t for t in fig.data
            if isinstance(t, go.Scatter)
            and getattr(t.marker, "symbol", None) == "x"
        ]
        assert len(cloud_traces) >= 1

    def test_cloud_marker_count_matches_cloudy_steps(
        self, timeseries_stack, scl_dir_with_files
    ):
        """The cloud trace must contain exactly as many points as cloudy steps."""
        scl_dir, cloudy_steps = scl_dir_with_files
        fig = plot_timeseries(
            timeseries_stack,
            lon=19.05, lat=47.45,
            bands=["ndvi"],
            scl_path=scl_dir,
        )
        cloud_traces = [
            t for t in fig.data
            if isinstance(t, go.Scatter)
            and getattr(t.marker, "symbol", None) == "x"
        ]
        total_cloud_points = sum(len(t.x) for t in cloud_traces)
        assert total_cloud_points == len(cloudy_steps)

    def test_cloud_marker_colour_is_red(self, timeseries_stack, scl_dir_with_files):
        scl_dir, _ = scl_dir_with_files
        fig = plot_timeseries(
            timeseries_stack,
            lon=19.05, lat=47.45,
            bands=["ndvi"],
            scl_path=scl_dir,
        )
        cloud_traces = [
            t for t in fig.data
            if isinstance(t, go.Scatter)
            and getattr(t.marker, "symbol", None) == "x"
        ]
        assert len(cloud_traces) >= 1
        assert cloud_traces[0].marker.color == "#e63946"

    def test_scl_single_file_accepted(self, timeseries_stack, scl_tif):
        """A single SCL file path (not a dir) must not raise."""
        fig = plot_timeseries(
            timeseries_stack,
            lon=19.05, lat=47.45,
            bands=["ndvi"],
            scl_path=scl_tif,
        )
        assert isinstance(fig, go.Figure)

    # ── layout & style ────────────────────────────────────────────────────────

    def test_title_contains_coordinates(self, timeseries_stack):
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"])
        assert "19.05" in fig.layout.title.text or "47.45" in fig.layout.title.text

    def test_title_contains_band_name(self, timeseries_stack):
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"])
        assert "ndvi" in fig.layout.title.text

    def test_dark_background(self, timeseries_stack):
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"])
        assert fig.layout.paper_bgcolor == "#0f1117"

    def test_hovermode_unified(self, timeseries_stack):
        fig = plot_timeseries(timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"])
        assert fig.layout.hovermode == "x unified"

    # ── save_html ─────────────────────────────────────────────────────────────

    def test_save_html(self, timeseries_stack, tmp_path):
        out = tmp_path / "ts.html"
        plot_timeseries(
            timeseries_stack, lon=19.05, lat=47.45, bands=["ndvi"],
            save_html=str(out),
        )
        assert out.exists()
        assert out.stat().st_size > 0

    # ── input validation ──────────────────────────────────────────────────────

    def test_no_time_dim_raises(self):
        xr = pytest.importorskip("xarray")
        da = xr.DataArray(np.zeros((4, 16, 16)), dims=["band", "y", "x"])
        with pytest.raises(ValueError, match="time"):
            plot_timeseries(da, lon=19.05, lat=47.45)

    def test_2d_input_raises(self):
        xr = pytest.importorskip("xarray")
        da = xr.DataArray(np.zeros((16, 16)), dims=["y", "x"])
        with pytest.raises(ValueError):
            plot_timeseries(da, lon=19.05, lat=47.45)


# ── SCL metadata ──────────────────────────────────────────────────────────────

class TestSclClasses:

    def test_all_12_classes_defined(self):
        assert len(SCL_CLASSES) == 12

    def test_classes_0_to_11(self):
        assert set(SCL_CLASSES.keys()) == set(range(12))

    def test_each_has_label_and_color(self):
        for v, entry in SCL_CLASSES.items():
            assert "label" in entry, f"Class {v} missing label"
            assert "color" in entry, f"Class {v} missing color"

    def test_default_bad_classes_subset_of_scl(self):
        assert DEFAULT_BAD_CLASSES.issubset(set(SCL_CLASSES.keys()))