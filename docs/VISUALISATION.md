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
from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask
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

## API reference

### plot_band

```python
def plot_band(
    source: str | Path,
    band: str | int | None = None,
    colorscale: str = "Viridis",
    title: str | None = None,
    save_html: str | Path | None = None,
) -> go.Figure:
```

Visualise a single band or index as an interactive heatmap.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source` | `str \| Path` | — | Path to `.nc` or `.tif` |
| `band` | `str \| int \| None` | `None` | Band name (`"nir"`, `"red"`, …), integer index, or `None` to squeeze (single-band files) |
| `colorscale` | `str` | `"Viridis"` | Any Plotly colorscale: `"RdYlGn"`, `"Greys"`, `"Plasma"`, … |
| `title` | `str \| None` | auto | Figure title |
| `save_html` | `str \| Path \| None` | `None` | Save path for HTML output |

### plot_rgb

```python
def plot_rgb(
    source: str | Path,
    red_band: str | int | None = None,
    green_band: str | int | None = None,
    blue_band: str | int | None = None,
    title: str | None = None,
    save_html: str | Path | None = None,
) -> go.Figure:
```

RGB preview with automatic percentile stretch (2–98%).

Two modes — determined by whether band arguments are passed:

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

### plot_grid

```python
def plot_grid(
    panels: list[dict],
    ncols: int = 3,
    colorscale: str = "Viridis",
    title: str | None = None,
    save_html: str | Path | None = None,
) -> go.Figure:
```

Grid of heatmaps from any mix of files and bands — spectral, indices, visual, technical.

Each panel is a `dict` with these keys:

| Key | Type | Required | Description |
|---|---|---|---|
| `"file"` | `str \| Path` | ✓ | Path to `.nc` or `.tif` |
| `"band"` | `str \| int \| None` | — | Band name or index; `None` squeezes (default) |
| `"label"` | `str` | — | Subplot title; auto-generated if omitted |
| `"colorscale"` | `str` | — | Per-panel colorscale override |

| Parameter | Type | Default | Description |
|---|---|---|---|
| `panels` | `list[dict]` | — | Panel definitions (see above) |
| `ncols` | `int` | `3` | Number of columns in the grid |
| `colorscale` | `str` | `"Viridis"` | Default colorscale for all panels |
| `title` | `str \| None` | auto | Figure title |
| `save_html` | `str \| Path \| None` | `None` | Save path for HTML output |

Files are cached in memory — referencing the same file multiple times with different bands does not re-open it.

### plot_mask

```python
def plot_mask(
    scl_source: str | Path,
    bad_classes: Sequence[int] | None = None,
    title: str | None = None,
    save_html: str | Path | None = None,
    return_mask: bool = False,
) -> go.Figure | tuple[go.Figure, np.ndarray]:
```

Binary cloud / artefact mask from an SCL layer. Produces three panels:

- **Binary heatmap** — green = clear, red = bad
- **Pie chart** — clear % vs bad %
- **SCL class breakdown** — horizontal bar chart for classes present in the file; bad classes highlighted in red

| Parameter | Type | Default | Description |
|---|---|---|---|
| `scl_source` | `str \| Path` | — | Path to SCL `.nc` or `.tif` |
| `bad_classes` | `Sequence[int] \| None` | `{0,1,3,8,9,10}` | SCL values treated as bad |
| `title` | `str \| None` | auto | Figure title |
| `save_html` | `str \| Path \| None` | `None` | Save path for HTML output |
| `return_mask` | `bool` | `False` | If `True`, also returns `np.ndarray` (1 = bad, 0 = clear) |

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

## Examples

### Single band or index

```python
from sentinel_processor.visualisation.plot import plot_band

# NIR band from raw spectral file
plot_band("data/spectral/budapest_20260526T095725.nc", band="nir").show()

# NDVI index (single-band .tif — no band arg needed)
plot_band(
    "data/indices/indices_budapest_20260526T095725_ndvi.tif",
    colorscale="RdYlGn",
    save_html="data/vis/ndvi.html",
).show()

# SWIR1 with custom colorscale
plot_band("data/spectral/budapest_20260526T095725.nc", band="swir16", colorscale="Inferno").show()
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

# all 10 spectral bands
plot_grid(
    [
        {"file": scene, "band": b, "label": b}
        for b in ["blue", "green", "red", "nir",
                  "rededge1", "rededge2", "rededge3",
                  "nir08", "swir16", "swir22"]
    ],
    ncols=5,
    save_html="data/vis/all_bands.html",
).show()

# visual channels + indices
plot_grid(
    [
        {"file": "data/visual/vis_budapest_20260526T095725.nc", "band": 0, "label": "Red"},
        {"file": "data/visual/vis_budapest_20260526T095725.nc", "band": 1, "label": "Green"},
        {"file": "data/visual/vis_budapest_20260526T095725.nc", "band": 2, "label": "Blue"},
        {"file": "data/indices/indices_budapest_20260526T095725_ndvi.tif", "label": "NDVI", "colorscale": "RdYlGn"},
        {"file": "data/indices/indices_budapest_20260526T095725_ndwi.tif", "label": "NDWI", "colorscale": "Blues"},
        {"file": "data/indices/indices_budapest_20260526T095725_ndbi.tif", "label": "NDBI", "colorscale": "Reds"},
    ],
    ncols=3,
    save_html="data/vis/grid.html",
).show()
```

### Binary mask

```python
from sentinel_processor.visualisation.plot import plot_mask

# default bad classes (clouds, shadows, no-data)
plot_mask("data/technical/scl_budapest_20260526T095725.nc").show()

# only flag cloud and cirrus
plot_mask(
    "data/technical/scl_budapest_20260526T095725.nc",
    bad_classes=[8, 9, 10],
    save_html="data/vis/cloud_mask.html",
).show()

# also get the numpy mask array
fig, mask_arr = plot_mask(
    "data/technical/scl_budapest_20260526T095725.nc",
    return_mask=True,
)
print(mask_arr.shape, mask_arr.sum(), "bad pixels")
```

### Full pipeline

```python
import sentinel_processor as sp
from sentinel_processor.indices.compute import compute_indices
from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask

# 1. download
results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=sp.DownloadConfig(bands=sp.SpectralBands.ALL, keep_items=1),
)

scene = "data/spectral/budapest_20260526T095725.nc"
vis   = "data/visual/vis_budapest_20260526T095725.nc"
scl   = "data/technical/scl_budapest_20260526T095725.nc"

# 2. compute indices
idx_paths = compute_indices(scene, ["ndvi", "ndwi", "ndbi"])

# 3. visualise
plot_rgb(vis, save_html="data/vis/rgb.html").show()
plot_band(scene, band="nir", save_html="data/vis/nir.html").show()
plot_band(idx_paths["ndvi"], colorscale="RdYlGn", save_html="data/vis/ndvi.html").show()
plot_mask(scl, save_html="data/vis/mask.html").show()
plot_grid(
    [
        {"file": scene,               "band": "red", "label": "Red"},
        {"file": scene,               "band": "nir", "label": "NIR"},
        {"file": idx_paths["ndvi"],                  "label": "NDVI", "colorscale": "RdYlGn"},
        {"file": idx_paths["ndwi"],                  "label": "NDWI", "colorscale": "Blues"},
    ],
    ncols=4,
    save_html="data/vis/summary.html",
).show()
```