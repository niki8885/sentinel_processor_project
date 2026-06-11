# Visualisation

Interactive Plotly figures for Sentinel-2 scenes, indices, and quality layers.
All functions return a `plotly.graph_objects.Figure` and optionally save a self-contained HTML file.

## Install

```bash
pip install sentinel-processor
```

`plotly` is included in the default dependencies — no extra install needed.

## Import

```python
from sentinel_processor.visualisation.plot import (
    plot_band,
    plot_rgb,
    plot_grid,
    plot_mask,
    plot_timeseries,
)
```

## Module layout

```
sentinel_processor/
└── visualisation/
    ├── __init__.py
    └── plot.py
```

## Output

All functions return a `plotly.graph_objects.Figure`.
Pass `save_html="path/to/file.html"` to write a self-contained interactive HTML file.
Call `.show()` to open in the browser, or use `.to_html()` / `.write_image()` for further export.

---

## API reference

### plot_band

```python
def plot_band(
    source: str | Path,
    band: str | int | None = None,
    colorscale: str = "Viridis",
    title: str | None = None,
    save_html: str | Path | None = None,
) -> go.Figure
```

Visualise a single band or index as an interactive heatmap.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source` | `str \| Path` | — | Path to `.nc` or `.tif` |
| `band` | `str \| int \| None` | `None` | Band name (`"nir"`, `"red"`, …), integer index, or `None` to squeeze (single-band files) |
| `colorscale` | `str` | `"Viridis"` | Any Plotly colorscale: `"RdYlGn"`, `"Greys"`, `"Plasma"`, … |
| `title` | `str \| None` | auto | Figure title |
| `save_html` | `str \| Path \| None` | `None` | Save path for HTML output |

---

### plot_rgb

```python
def plot_rgb(
    source: str | Path,
    red_band: str | int | None = None,
    green_band: str | int | None = None,
    blue_band: str | int | None = None,
    title: str | None = None,
    save_html: str | Path | None = None,
) -> go.Figure
```

RGB preview with automatic percentile stretch (2–98 %).

**Mode A — visual file** (pre-made 3-band overview from the downloader):
Pass only `source`. Bands are taken in order: index 0 = R, 1 = G, 2 = B.

**Mode B — raw spectral** (build RGB from any three bands):
Pass `source` + explicit band names. Supports true colour and any false-colour combination.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source` | `str \| Path` | — | Path to `.nc` or `.tif` |
| `red_band` | `str \| int \| None` | `None` | Red channel selector (mode B) |
| `green_band` | `str \| int \| None` | `None` | Green channel selector (mode B) |
| `blue_band` | `str \| int \| None` | `None` | Blue channel selector (mode B) |
| `title` | `str \| None` | auto | Figure title |
| `save_html` | `str \| Path \| None` | `None` | Save path for HTML output |

---

### plot_grid

```python
def plot_grid(
    panels: list[dict],
    ncols: int = 3,
    colorscale: str = "Viridis",
    title: str | None = None,
    save_html: str | Path | None = None,
) -> go.Figure
```

Grid of heatmaps from any mix of files and bands.

Each panel is a `dict`:

| Key | Type | Required | Description |
|---|---|---|---|
| `"file"` | `str \| Path` | ✓ | Path to `.nc` or `.tif` |
| `"band"` | `str \| int \| None` | — | Band name or index; `None` squeezes (default) |
| `"label"` | `str` | — | Subplot title; auto-generated if omitted |
| `"colorscale"` | `str` | — | Per-panel colorscale override |

Files are cached in memory — referencing the same file multiple times with different bands does not re-open it.

---

### plot_mask

```python
def plot_mask(
    scl_source: str | Path,
    bad_classes: Sequence[int] | None = None,
    title: str | None = None,
    save_html: str | Path | None = None,
    return_mask: bool = False,
) -> go.Figure | tuple[go.Figure, np.ndarray]
```

Binary cloud / artefact mask from an SCL layer. Produces three panels:

- **Binary heatmap** — green = clear, red = bad
- **Pie chart** — clear % vs bad %
- **SCL class breakdown** — horizontal bar chart; bad classes highlighted in red

| Parameter | Type | Default | Description |
|---|---|---|---|
| `scl_source` | `str \| Path` | — | Path to SCL `.nc` or `.tif` |
| `bad_classes` | `Sequence[int] \| None` | `{0,1,3,8,9,10}` | SCL values treated as bad |
| `title` | `str \| None` | auto | Figure title |
| `save_html` | `str \| Path \| None` | `None` | Save path for HTML output |
| `return_mask` | `bool` | `False` | If `True`, also returns `np.ndarray` (1 = bad, 0 = clear) |

---

### plot_timeseries

```python
def plot_timeseries(
    stack: xr.DataArray | str | Path,
    lon: float,
    lat: float,
    bands: list[str] | None = None,
    scl_path: str | Path | None = None,
    agg_bbox: float | None = None,
    save_html: str | Path | None = None,
) -> go.Figure
```

Interactive time-series chart for a specific pixel or small region.
The primary visual debugging tool for all time-series workflows.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `stack` | `xr.DataArray \| str \| Path` | — | `(time, band, y, x)` DataArray **or** path to `.nc` / `.tif` stack |
| `lon` | `float` | — | WGS-84 longitude of the query point |
| `lat` | `float` | — | WGS-84 latitude of the query point |
| `bands` | `list[str] \| None` | `None` | Band names to plot. `None` → all bands in the stack |
| `scl_path` | `str \| Path \| None` | `None` | SCL file **or** directory of per-scene SCL files; enables cloud / shadow markers |
| `agg_bbox` | `float \| None` | `None` | Half-width in degrees of a spatial averaging box. `None` = nearest pixel |
| `save_html` | `str \| Path \| None` | `None` | Save path for HTML output |

**Extraction modes**

- **Single-pixel** (`agg_bbox=None`): nearest-pixel selection.
  Projected CRS (e.g. UTM) is handled — the query point is reprojected automatically.
- **Region** (`agg_bbox=<float>`): clips a `±agg_bbox`° bbox and returns the spatial
  mean of all finite, non-nodata pixels. Useful for field-scale averages.

**Cloud / shadow overlay** (requires `scl_path`)

SCL classes `{0, 1, 3, 8, 9, 10}` are treated as contaminated:

- Contaminated observations → **red ✕ markers**
- Clear observations → **filled circles** in the band colour

`scl_path` can be:
- A **single SCL file** — applied to every time step
- A **directory** — matched per time step using `scl_<timestamp>*.tif` naming
  (the layout produced by `download_sentinel2`)

**Multiple bands** are plotted on the same figure with separate coloured traces
(sky-blue → green → amber → pink → lavender → teal → yellow, cycling).

---

## SCL class reference

| Value | Label | Default bad |
|---|---|---|
| 0 | No Data | ✓ |
| 1 | Saturated / Defective | ✓ |
| 2 | Dark Area | — |
| 3 | Cloud Shadow | ✓ |
| 4 | Vegetation | — |
| 5 | Bare Soil | — |
| 6 | Water | — |
| 7 | Unclassified | — |
| 8 | Cloud (medium prob.) | ✓ |
| 9 | Cloud (high prob.) | ✓ |
| 10 | Thin Cirrus | ✓ |
| 11 | Snow / Ice | — |

---

## Examples

### Single band or index

```python
from sentinel_processor.visualisation.plot import plot_band

# NIR band from raw spectral file
plot_band("data/spectral/budapest_20260526T095725.nc", band="nir").show()

# NDVI index (single-band .tif)
plot_band(
    "data/indices/indices_budapest_20260526T095725_ndvi.tif",
    colorscale="RdYlGn",
    save_html="data/vis/ndvi.html",
).show()
```

### RGB preview

```python
from sentinel_processor.visualisation.plot import plot_rgb

# from pre-downloaded visual file (mode A)
plot_rgb("data/visual/vis_budapest_20260526T095725.nc").show()

# true colour from raw spectral bands (mode B)
plot_rgb(
    "data/spectral/budapest_20260526T095725.nc",
    red_band="red", green_band="green", blue_band="blue",
).show()

# false colour: NIR / Red / Green
plot_rgb(
    "data/spectral/budapest_20260526T095725.nc",
    red_band="nir", green_band="red", blue_band="green",
    title="False Colour (NIR-R-G)",
    save_html="data/vis/false_colour.html",
).show()
```

### Grid

```python
from sentinel_processor.visualisation.plot import plot_grid

scene = "data/spectral/budapest_20260526T095725.nc"

plot_grid(
    [{"file": scene, "band": b, "label": b}
     for b in ["blue", "green", "red", "nir", "swir16", "swir22"]],
    ncols=3,
    save_html="data/vis/bands.html",
).show()
```

### Binary mask

```python
from sentinel_processor.visualisation.plot import plot_mask

plot_mask("data/technical/scl_budapest_20260526T095725.nc").show()

fig, mask_arr = plot_mask(
    "data/technical/scl_budapest_20260526T095725.nc",
    return_mask=True,
)
print(mask_arr.shape, mask_arr.sum(), "bad pixels")
```

### NDVI + EVI time series for a crop field with cloud markers

```python
from pathlib import Path
from sentinel_processor.visualisation.plot import plot_timeseries

# Maize field near Kecskemét, Hungary
LON, LAT = 19.688, 46.901

fig = plot_timeseries(
    stack="output/kecskémet_stack.nc",   # (time, band, y, x) NetCDF
    lon=LON,
    lat=LAT,
    bands=["ndvi", "evi"],               # both indices on one chart
    scl_path="output/technical/",        # per-scene SCL directory → cloud markers
    agg_bbox=0.001,                      # ~100 m spatial average for field scale
    save_html="output/crop_ts.html",
)
fig.show()
```

The chart renders two coloured lines (NDVI in sky-blue, EVI in green).
Cloud-contaminated acquisitions appear as **red ✕** marks; clear observations
get **filled circle** markers. Hover over any point for the exact date and value.

### Single-pixel inspection from an in-memory DataArray

```python
from sentinel_processor.visualisation.plot import plot_timeseries

# stack_result is a StackResult from stack_timeseries()
fig = plot_timeseries(
    stack=stack_result.stack,
    lon=19.05,
    lat=47.49,
    bands=["ndvi"],
)
fig.show()
```

### Region mean, no cloud overlay

```python
fig = plot_timeseries(
    stack="output/budapest_stack.nc",
    lon=19.040,
    lat=47.498,
    bands=["ndvi", "swir16"],
    agg_bbox=0.005,   # ~500 m spatial average
)
fig.show()
```

---

### Full pipeline

```python
import sentinel_processor as sp
from sentinel_processor.indices.compute import compute_indices
from sentinel_processor.visualisation.plot import (
    plot_band, plot_rgb, plot_grid, plot_mask, plot_timeseries,
)

# 1. download
results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=sp.DownloadConfig(bands=sp.SpectralBands.ALL, keep_items=6),
)

scene = "data/spectral/budapest_20260526T095725.nc"
scl   = "data/technical/scl_budapest_20260526T095725.nc"

# 2. compute indices
idx_paths = compute_indices(scene, ["ndvi", "evi", "ndwi"])

# 3. single-scene views
plot_rgb(scene, red_band="red", green_band="green", blue_band="blue",
         save_html="data/vis/rgb.html").show()
plot_band(scene, band="nir", save_html="data/vis/nir.html").show()
plot_mask(scl, save_html="data/vis/mask.html").show()

# 4. time series from a stack built by stack_timeseries()
from sentinel_processor.timeseries import stack_timeseries, TimeSeriesConfig
from pathlib import Path

spectral_files = sorted(Path("data/spectral").glob("budapest_*.nc"))
result = stack_timeseries(spectral_files, scl_dir="data/technical/")

plot_timeseries(
    stack=result.stack,
    lon=19.17, lat=47.56,
    bands=["ndvi", "evi"],
    scl_path="data/technical/",
    agg_bbox=0.002,
    save_html="data/vis/timeseries.html",
).show()
```

---

## Dark-theme colour reference

All functions share a consistent dark palette.

| Role | Hex |
|---|---|
| Paper / outer background | `#0f1117` |
| Plot area (time series) | `#1a1d27` |
| Grid lines | `#2a2d3a` |
| Cloud / shadow marker | `#e63946` |
| Clear vegetation mask | `#3cb371` |
| Bad-pixel mask | `#e63946` |

Band colours cycle through:
`#4fc3f7` (sky-blue) → `#81c784` (green) → `#ffb74d` (amber) →
`#f06292` (pink) → `#ce93d8` (lavender) → `#80cbc4` (teal) → `#fff176` (yellow).

---

## Running the tests

```bash
# from the repo root
pytest tests/test_visualisation.py -v

# run only the new time-series tests
pytest tests/test_visualisation.py -v -k "timeseries"

# skip if optional deps are absent (rasterio, plotly)
pytest tests/test_visualisation.py -v --ignore-glob="*integration*"
```

Expected output (all deps present):

```
tests/test_visualisation.py::TestPstretch::test_output_range_0_to_1         PASSED
tests/test_visualisation.py::TestPstretch::test_dtype_float32                PASSED
...
tests/test_visualisation.py::TestPlotTimeseries::test_returns_figure         PASSED
tests/test_visualisation.py::TestPlotTimeseries::test_line_traces_equal_band_count PASSED
tests/test_visualisation.py::TestPlotTimeseries::test_cloud_marker_count_matches_cloudy_steps PASSED
...
```