from __future__ import annotations
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    import plotly.graph_objects as go

logger = logging.getLogger(__name__)

SCL_CLASSES: dict[int, dict] = {
    0: {"label": "No Data", "color": "#000000"},
    1: {"label": "Saturated / Defective", "color": "#ff0000"},
    2: {"label": "Dark Area", "color": "#2f2f2f"},
    3: {"label": "Cloud Shadow", "color": "#643200"},
    4: {"label": "Vegetation", "color": "#00a000"},
    5: {"label": "Bare Soil", "color": "#ffe65a"},
    6: {"label": "Water", "color": "#0000ff"},
    7: {"label": "Unclassified", "color": "#808080"},
    8: {"label": "Cloud (med.)", "color": "#c0c0c0"},
    9: {"label": "Cloud (high)", "color": "#ffffff"},
    10: {"label": "Thin Cirrus", "color": "#64c8ff"},
    11: {"label": "Snow / Ice", "color": "#ff96ff"},
}

DEFAULT_BAD_CLASSES = {0, 1, 3, 8, 9, 10}

_CLEAR_COLOR = "#3cb371"
_BAD_COLOR = "#e63946"


def _load_da(path: str | Path) -> xr.DataArray:
    path = Path(path)
    if path.suffix.lower() in (".tif", ".tiff"):
        import rioxarray
        return rioxarray.open_rasterio(path)
    ds = xr.open_dataset(path)
    var = next((v for v in ds.data_vars if v != "spatial_ref"), list(ds.data_vars)[0])
    return ds[var]


def _extract_band(da: xr.DataArray, band: str | int | None) -> np.ndarray:
    if band is None:
        return np.asarray(da.squeeze(), dtype=np.float64)
    if isinstance(band, str):
        if "band" in da.coords:
            coord_vals = [str(v) for v in da.coords["band"].values]
            if band in coord_vals:
                return np.asarray(da.sel(band=band).squeeze(), dtype=np.float64)

        try:
            return np.asarray(da.isel(band=int(band)).squeeze(), dtype=np.float64)
        except Exception:
            raise ValueError(
                f"Band '{band}' not found. Available: {coord_vals if 'band' in da.coords else 'no band coord'}")
    return np.asarray(da.isel(band=band).squeeze(), dtype=np.float64)


def _pstretch(arr: np.ndarray, lo: float = 2, hi: float = 98) -> np.ndarray:
    valid = arr[np.isfinite(arr)]
    if valid.size == 0:
        return np.zeros_like(arr, dtype=np.float32)
    p_lo, p_hi = np.percentile(valid, [lo, hi])
    if p_hi <= p_lo:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - p_lo) / (p_hi - p_lo), 0.0, 1.0).astype(np.float32)


def _dark_layout(fig, title: str, height: int = 480) -> None:
    fig.update_layout(
        title=dict(text=title, font_size=14, font_color="#e0e0e0"),
        height=height,
        margin=dict(l=10, r=10, t=55, b=10),
        paper_bgcolor="#0f1117",
        plot_bgcolor="#0f1117",
        font=dict(color="#e0e0e0"),
    )
    fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
    fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False)


def _save(fig, path: str | Path | None) -> None:
    if path is None:
        return
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(p), include_plotlyjs="cdn")
    logger.info(f"[vis] → {p}")


def plot_band(
        source: str | Path,
        band: str | int | None = None,
        colorscale: str = "Viridis",
        title: str | None = None,
        save_html: str | Path | None = None,
) -> "go.Figure":
    import plotly.graph_objects as go

    da = _load_da(source)
    arr = _extract_band(da, band)
    arr_s = _pstretch(arr)

    band_label = str(band) if band is not None else Path(source).stem
    fig = go.Figure(go.Heatmap(
        z=arr_s,
        colorscale=colorscale,
        showscale=True,
        colorbar=dict(thickness=12, title=dict(text="scaled", side="right")),
        hovertemplate="x: %{x}<br>y: %{y}<br>val: %{z:.3f}<extra></extra>",
    ))
    _dark_layout(fig, title or f"{band_label} — {Path(source).name}")
    _save(fig, save_html)
    return fig


def plot_rgb(
        source: str | Path,
        red_band: str | int | None = None,
        green_band: str | int | None = None,
        blue_band: str | int | None = None,
        title: str | None = None,
        save_html: str | Path | None = None,
) -> "go.Figure":
    import plotly.graph_objects as go

    da = _load_da(source)

    if red_band is None and green_band is None and blue_band is None:
        arr = np.asarray(da, dtype=np.float64)
        if arr.ndim == 2:
            arr = np.stack([arr, arr, arr], axis=0)
        rgb_raw = np.moveaxis(arr[:3], 0, -1)
    else:
        r = _extract_band(da, red_band)
        g = _extract_band(da, green_band)
        b = _extract_band(da, blue_band)
        rgb_raw = np.stack([r, g, b], axis=-1)

    rgb_uint8 = (
            np.stack([_pstretch(rgb_raw[..., c]) for c in range(3)], axis=-1) * 255
    ).astype(np.uint8)

    fig = go.Figure(go.Image(z=rgb_uint8, hoverinfo="skip"))
    _dark_layout(fig, title or f"RGB — {Path(source).name}")
    _save(fig, save_html)
    return fig


def plot_grid(
        panels: Sequence[dict],
        ncols: int = 3,
        colorscale: str = "Viridis",
        title: str | None = None,
        save_html: str | Path | None = None,
) -> "go.Figure":
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    n = len(panels)
    if n == 0:
        raise ValueError("panels list is empty.")
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))

    subplot_titles = [
        p.get("label") or f"{Path(p['file']).stem} / {p.get('band', 'band')}"
        for p in panels
    ]
    subplot_titles += [""] * (nrows * ncols - n)

    fig = make_subplots(
        rows=nrows, cols=ncols,
        subplot_titles=subplot_titles,
        vertical_spacing=0.05,
        horizontal_spacing=0.03,
    )

    _cache: dict[str, xr.DataArray] = {}

    for idx, panel in enumerate(panels):
        row = idx // ncols + 1
        col = idx % ncols + 1
        fpath = str(panel["file"])
        band = panel.get("band")
        cs = panel.get("colorscale", colorscale)

        if fpath not in _cache:
            _cache[fpath] = _load_da(fpath)
        da = _cache[fpath]
        arr = _extract_band(da, band)
        arr_s = _pstretch(arr)

        show_scale = (idx == 0)
        fig.add_trace(
            go.Heatmap(
                z=arr_s,
                colorscale=cs,
                showscale=show_scale,
                colorbar=dict(thickness=10, len=0.3) if show_scale else None,
                hovertemplate="x: %{x}<br>y: %{y}<br>val: %{z:.3f}<extra></extra>",
            ),
            row=row, col=col,
        )

    _dark_layout(fig, title or "Scene Grid", height=280 * nrows + 80)
    _save(fig, save_html)
    return fig


def plot_mask(
        scl_source: str | Path,
        bad_classes: Sequence[int] | None = None,
        title: str | None = None,
        save_html: str | Path | None = None,
        return_mask: bool = False,
) -> "go.Figure | tuple[go.Figure, np.ndarray]":
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    bad_set = set(bad_classes) if bad_classes is not None else DEFAULT_BAD_CLASSES
    da = _load_da(scl_source)
    scl = np.asarray(da.squeeze(), dtype=np.int32)
    mask = np.isin(scl, list(bad_set)).astype(np.uint8)

    bad_pct = 100.0 * mask.sum() / mask.size
    good_pct = 100.0 - bad_pct

    unique_vals, counts = np.unique(scl, return_counts=True)
    class_labels = [SCL_CLASSES.get(int(v), {"label": f"SCL {v}"})["label"] for v in unique_vals]
    pct_vals = 100.0 * counts / counts.sum()
    is_bad = [int(v) in bad_set for v in unique_vals]

    fig = make_subplots(
        rows=1, cols=3,
        column_widths=[0.50, 0.20, 0.30],
        subplot_titles=["Binary Mask", "Clear / Bad", "SCL Classes"],
        specs=[[{"type": "heatmap"}, {"type": "pie"}, {"type": "bar"}]],
        horizontal_spacing=0.06,
    )

    fig.add_trace(
        go.Heatmap(
            z=mask,
            colorscale=[[0, _CLEAR_COLOR], [1, _BAD_COLOR]],
            zmin=0, zmax=1,
            showscale=False,
            hovertemplate="x: %{x}<br>y: %{y}<br>bad: %{z}<extra></extra>",
        ),
        row=1, col=1,
    )

    fig.add_trace(
        go.Pie(
            labels=["Clear", "Bad"],
            values=[good_pct, bad_pct],
            marker_colors=[_CLEAR_COLOR, _BAD_COLOR],
            textinfo="label+percent",
            hole=0.4,
        ),
        row=1, col=2,
    )

    bar_colors = [_BAD_COLOR if b else _CLEAR_COLOR for b in is_bad]
    fig.add_trace(
        go.Bar(
            x=pct_vals,
            y=class_labels,
            orientation="h",
            marker_color=bar_colors,
            text=[f"{p:.1f}%" for p in pct_vals],
            textposition="outside",
            hovertemplate="%{y}: %{x:.2f}%<extra></extra>",
        ),
        row=1, col=3,
    )

    stem = Path(scl_source).stem
    bad_list = sorted(bad_set)
    _dark_layout(
        fig,
        title or f"SCL Mask — {stem}  (bad classes: {bad_list})",
        height=440,
    )
    fig.update_layout(showlegend=False)

    _save(fig, save_html)
    if return_mask:
        return fig, mask
    return fig


_TS_CLOUDY_CLASSES = {0, 1, 3, 8, 9, 10}

_BAND_COLOURS = [
    "#4fc3f7",  # sky-blue   – first band
    "#81c784",  # green      – second band
    "#ffb74d",  # amber      – third band
    "#f06292",  # pink       – fourth band
    "#ce93d8",  # lavender   – fifth band
    "#80cbc4",  # teal       – sixth band
    "#fff176",  # yellow     – seventh band
]


def _reproject_lonlat_to_crs(
        lon: float,
        lat: float,
        crs,
) -> tuple[float, float]:
    """Convert (lon, lat) WGS-84 to the native CRS of *crs*.

    Falls back to the raw values when pyproj is unavailable or CRS is
    already geographic.
    """
    try:
        from pyproj import Transformer

        transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
        x, y = transformer.transform(lon, lat)
        return float(x), float(y)
    except Exception:
        return lon, lat


def _nearest_pixel_value(
        da: "xr.DataArray",
        lon: float,
        lat: float,
) -> float:
    """Return the value of the pixel nearest to (lon, lat).

    Handles projected CRS by reprojecting the query point first;
    falls back to direct coordinate selection otherwise.
    """

    crs = None
    try:
        crs = da.rio.crs
    except Exception:
        pass

    if crs is not None and not crs.is_geographic:
        qx, qy = _reproject_lonlat_to_crs(lon, lat, crs)
        coord_x, coord_y = "x", "y"
    else:
        qx, qy = lon, lat
        coord_x = "x" if "x" in da.coords else "lon"
        coord_y = "y" if "y" in da.coords else "lat"

    kwargs = {coord_x: qx, coord_y: qy}
    val = float(da.sel(method="nearest", **kwargs).values)
    return val


def _bbox_mean_value(
        da: "xr.DataArray",
        lon: float,
        lat: float,
        half_deg: float,
) -> float:
    """Return the spatial mean over a (lon±half_deg, lat±half_deg) bbox."""
    try:
        clipped = da.rio.clip_box(
            minx=lon - half_deg,
            miny=lat - half_deg,
            maxx=lon + half_deg,
            maxy=lat + half_deg,
            crs="EPSG:4326",
        )
        vals = np.asarray(clipped.values, dtype=np.float64)
        valid = vals[np.isfinite(vals) & (vals != NODATA)]
        return float(np.mean(valid)) if valid.size else float("nan")
    except Exception:
        return _nearest_pixel_value(da, lon, lat)


def _scl_is_cloudy(scl_path: "str | Path", lon: float, lat: float) -> bool:
    """Return True when the SCL pixel nearest (lon, lat) is cloud/shadow."""
    try:
        import rioxarray  # noqa: F401

        da = _load_da(scl_path).squeeze()
        val = int(round(_nearest_pixel_value(da, lon, lat)))
        return val in _TS_CLOUDY_CLASSES
    except Exception:
        return False


NODATA: float = -9999.0


def plot_timeseries(
        stack: "xr.DataArray | str | Path",
        lon: float,
        lat: float,
        bands: list[str] | None = None,
        scl_path: "str | Path | None" = None,
        agg_bbox: "float | None" = None,
        save_html: "str | Path | None" = None,
) -> "go.Figure":
    """Interactive pixel / region time-series chart.

    Parameters
    ----------
    stack : xr.DataArray or str/Path
        A ``(time, band, y, x)`` DataArray **or** a path to a NetCDF / GeoTIFF
        file that will be opened with :func:`_load_da`.
    lon, lat : float
        WGS-84 longitude / latitude of the point of interest.
    bands : list[str] | None
        Band names to plot (matched against the ``band`` coordinate).
        ``None`` → all bands in the stack.
    scl_path : str | Path | None
        Path to a *single* SCL scene **or** a directory that contains one
        SCL file per time step (naming convention ``scl_*<timestamp>*.tif``).
        When supplied, cloud / shadow observations are highlighted as red
        markers; clear observations get filled circles.
    agg_bbox : float | None
        When set, spatially average all pixels within a bounding box of
        ``±agg_bbox`` degrees around ``(lon, lat)`` instead of sampling the
        nearest single pixel.
    save_html : str | Path | None
        If given, write the figure as a self-contained HTML file (same
        behaviour as all other ``plot_*`` functions in this module).

    Returns
    -------
    plotly.graph_objects.Figure

    Notes
    -----
    *Single-pixel mode* (``agg_bbox=None``): extracts the nearest pixel to
    ``(lon, lat)`` at every time step.

    *Region mode* (``agg_bbox`` is a float): computes the spatial mean of
    all valid (finite, non-nodata) pixels inside the bounding box.

    When ``scl_path`` points to a **directory**, the function matches each
    time step to ``scl_<scene_stem>*.tif`` files inside that directory.
    When it points to a **single file**, that file is used for every time
    step (useful for quick single-scene checks).

    Cloud / shadow observations (SCL classes 0, 1, 3, 8, 9, 10) are shown
    as red ✕ markers on top of the line; clear observations get filled
    circles.  The SCL overlay is optional — omit ``scl_path`` to disable.
    """
    import plotly.graph_objects as go
    import xarray as xr

    if not isinstance(stack, xr.DataArray):
        stack = _load_da(stack)

    if stack.ndim == 2:
        raise ValueError(
            "plot_timeseries() requires a 3-D (time, …) or 4-D "
            "(time, band, y, x) DataArray."
        )

    if "time" not in stack.dims:
        raise ValueError("DataArray must have a 'time' dimension.")

    if stack.ndim == 3:
        stack = stack.expand_dims(dim={"band": ["value"]}, axis=1)

    if stack.dims != ("time", "band", "y", "x"):
        try:
            stack = stack.transpose("time", "band", "y", "x")
        except ValueError:
            pass

    if "band" in stack.coords:
        all_bands = [str(b) for b in stack.coords["band"].values]
    else:
        all_bands = [f"band_{i}" for i in range(stack.sizes.get("band", 1))]

    plot_bands = bands if bands is not None else all_bands

    times_raw = stack.coords["time"].values
    try:
        import pandas as pd
        times = pd.to_datetime(times_raw).to_pydatetime().tolist()
    except Exception:
        times = list(times_raw)

    n_times = len(times)

    cloudy_mask: list[bool] = [False] * n_times

    if scl_path is not None:
        scl_p = Path(scl_path)
        if scl_p.is_dir():
            scl_files = sorted(scl_p.glob("scl_*.tif")) + sorted(scl_p.glob("scl_*.nc"))
            for ti, t in enumerate(times):
                ts_str = (
                    t.strftime("%Y%m%dT%H%M%S")
                    if hasattr(t, "strftime") else str(t)[:19].replace("-", "").replace(":", "")
                )
                matched = [f for f in scl_files if ts_str in f.name]
                if matched:
                    cloudy_mask[ti] = _scl_is_cloudy(matched[0], lon, lat)
        else:
            single_cloudy = _scl_is_cloudy(scl_p, lon, lat)
            cloudy_mask = [single_cloudy] * n_times

    band_series: dict[str, list[float]] = {}
    extract_fn = (
        (lambda da: _bbox_mean_value(da, lon, lat, agg_bbox))
        if agg_bbox is not None
        else (lambda da: _nearest_pixel_value(da, lon, lat))
    )

    for band_name in plot_bands:
        if band_name in all_bands:
            bidx = all_bands.index(band_name)
        else:
            logger.warning(
                f"[plot_timeseries] Band '{band_name}' not found in stack "
                f"(available: {all_bands}); skipping."
            )
            continue

        vals: list[float] = []
        for ti in range(n_times):
            da_slice = stack.isel(time=ti, band=bidx)
            vals.append(extract_fn(da_slice))
        band_series[band_name] = vals

    if not band_series:
        raise ValueError(
            f"[plot_timeseries] None of the requested bands {plot_bands} "
            f"were found in the stack (available: {all_bands})."
        )

    fig = go.Figure()

    clear_idx = [i for i, c in enumerate(cloudy_mask) if not c]
    cloud_idx = [i for i, c in enumerate(cloudy_mask) if c]

    for b_idx, (band_name, vals) in enumerate(band_series.items()):
        colour = _BAND_COLOURS[b_idx % len(_BAND_COLOURS)]
        vals_arr = np.array(vals, dtype=np.float64)
        vals_arr[(vals_arr == NODATA) | ~np.isfinite(vals_arr)] = np.nan

        t_all = [times[i] for i in range(n_times)]

        fig.add_trace(go.Scatter(
            x=t_all,
            y=vals_arr.tolist(),
            mode="lines",
            name=band_name,
            line=dict(color=colour, width=2),
            legendgroup=band_name,
            showlegend=True,
            hovertemplate=(
                f"<b>{band_name}</b><br>"
                "Date: %{x|%Y-%m-%d}<br>"
                "Value: %{y:.4f}<extra></extra>"
            ),
        ))

        if clear_idx:
            t_clear = [times[i] for i in clear_idx]
            v_clear = [vals_arr[i] for i in clear_idx]
            fig.add_trace(go.Scatter(
                x=t_clear,
                y=v_clear,
                mode="markers",
                name=f"{band_name} (clear)",
                marker=dict(color=colour, size=7, symbol="circle"),
                legendgroup=band_name,
                showlegend=scl_path is not None,
                hovertemplate=(
                    f"<b>{band_name}</b> – clear<br>"
                    "Date: %{x|%Y-%m-%d}<br>"
                    "Value: %{y:.4f}<extra></extra>"
                ),
            ))

        # cloud/shadow markers (red X)
        if cloud_idx and scl_path is not None:
            t_cloud = [times[i] for i in cloud_idx]
            v_cloud = [vals_arr[i] for i in cloud_idx]
            first_cloud = b_idx == 0  # show legend entry once
            fig.add_trace(go.Scatter(
                x=t_cloud,
                y=v_cloud,
                mode="markers",
                name="cloud / shadow" if first_cloud else f"cloud ({band_name})",
                marker=dict(
                    color="#e63946",
                    size=10,
                    symbol="x",
                    line=dict(width=2, color="#e63946"),
                ),
                legendgroup="cloud",
                showlegend=first_cloud,
                hovertemplate=(
                    f"<b>{band_name}</b> – cloud/shadow<br>"
                    "Date: %{x|%Y-%m-%d}<br>"
                    "Value: %{y:.4f}<extra></extra>"
                ),
            ))

    mode_label = (
        f"region ±{agg_bbox}°" if agg_bbox is not None else "pixel"
    )
    band_label = ", ".join(band_series.keys())
    auto_title = (
        f"Time Series [{band_label}] — ({lat:.4f}, {lon:.4f})  [{mode_label}]"
    )

    fig.update_layout(
        title=dict(text=auto_title, font_size=14, font_color="#e0e0e0"),
        height=480,
        margin=dict(l=60, r=20, t=60, b=50),
        paper_bgcolor="#0f1117",
        plot_bgcolor="#1a1d27",
        font=dict(color="#e0e0e0"),
        legend=dict(
            bgcolor="rgba(15,17,23,0.8)",
            bordercolor="#333",
            borderwidth=1,
            font_size=11,
        ),
        hovermode="x unified",
        xaxis=dict(
            title="Date",
            showgrid=True,
            gridcolor="#2a2d3a",
            zeroline=False,
            tickfont=dict(size=11),
        ),
        yaxis=dict(
            title="Index value",
            showgrid=True,
            gridcolor="#2a2d3a",
            zeroline=False,
            tickfont=dict(size=11),
        ),
    )

    _save(fig, save_html)
    return fig
