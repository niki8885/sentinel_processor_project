from __future__ import annotations
import logging
from pathlib import Path
from typing import Sequence
import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)

SCL_CLASSES: dict[int, dict] = {
    0:  {"label": "No Data",               "color": "#000000"},
    1:  {"label": "Saturated / Defective", "color": "#ff0000"},
    2:  {"label": "Dark Area",             "color": "#2f2f2f"},
    3:  {"label": "Cloud Shadow",          "color": "#643200"},
    4:  {"label": "Vegetation",            "color": "#00a000"},
    5:  {"label": "Bare Soil",             "color": "#ffe65a"},
    6:  {"label": "Water",                 "color": "#0000ff"},
    7:  {"label": "Unclassified",          "color": "#808080"},
    8:  {"label": "Cloud (med.)",          "color": "#c0c0c0"},
    9:  {"label": "Cloud (high)",          "color": "#ffffff"},
    10: {"label": "Thin Cirrus",           "color": "#64c8ff"},
    11: {"label": "Snow / Ice",            "color": "#ff96ff"},
}

DEFAULT_BAD_CLASSES = {0, 1, 3, 8, 9, 10}

_CLEAR_COLOR = "#3cb371"
_BAD_COLOR   = "#e63946"

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
            raise ValueError(f"Band '{band}' not found. Available: {coord_vals if 'band' in da.coords else 'no band coord'}")
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

    da  = _load_da(source)
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
        band  = panel.get("band")
        cs    = panel.get("colorscale", colorscale)

        if fpath not in _cache:
            _cache[fpath] = _load_da(fpath)
        da  = _cache[fpath]
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
    da  = _load_da(scl_source)
    scl = np.asarray(da.squeeze(), dtype=np.int32)
    mask = np.isin(scl, list(bad_set)).astype(np.uint8)

    bad_pct  = 100.0 * mask.sum() / mask.size
    good_pct = 100.0 - bad_pct

    unique_vals, counts = np.unique(scl, return_counts=True)
    class_labels = [SCL_CLASSES.get(int(v), {"label": f"SCL {v}"})["label"] for v in unique_vals]
    class_colors = [SCL_CLASSES.get(int(v), {"color": "#888888"})["color"]  for v in unique_vals]
    pct_vals     = 100.0 * counts / counts.sum()
    is_bad       = [int(v) in bad_set for v in unique_vals]

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