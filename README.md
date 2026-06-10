# sentinel-processor

Sentinel-2 L2A downloader and processing toolkit built on the [Element84 STAC API](https://earth-search.aws.element84.com/v1).

Downloads spectral bands, quality layers, and visual overviews for any coordinate. Validation, spectral index computation, convolution filters, pansharpening, time-series stacking, and band covariance are all backed by compiled Fortran kernels — Python handles I/O and orchestration, Fortran handles the pixels.

[![ko-fi](https://ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/Z8Z01TOFUW)
[![PyPI](https://img.shields.io/pypi/v/sentinel-processor)](https://pypi.org/project/sentinel-processor/)
[![Python](https://img.shields.io/badge/python-3.10%20|%203.11%20|%203.12-blue)](https://pypi.org/project/sentinel-processor/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE.md)

---

## Features

| Module | What it does |
|---|---|
| **Downloader** | STAC search → parallel band fetch → validation → save `.tif` / `.nc` |
| **Validation** | SCL cloud/snow analysis, radiometry check, dimension check (Fortran) |
| **Indices** | 10 spectral indices: NDVI, EVI, SAVI, NDWI, MNDWI, NDBI, NBR, NDSI, CIG, ARVI (Fortran) |
| **Filters** | 14 convolution and morphological filters: Gaussian, bilateral, Sobel, Laplacian, unsharp mask, median, erosion, dilation, top-hat, arbitrary kernel (Fortran) |
| **Pansharpening** | Gram-Schmidt, IHS, Wavelet — inject PAN detail into MS bands (Fortran) |
| **Time series** | Quality-filtered temporal stack builder with cloud/snow filtering, alignment, and save (Fortran validation + raster ops) |
| **Gap filling** | Fill cloud-masked holes in time stacks: linear, Savitzky-Golay, PCHIP, Holt ETS, Gaussian (Fortran) |
| **Texture** | GLCM texture features per pixel: energy, contrast, homogeneity — window, distance, and angle-configurable (Fortran) |
| **Analysis** | 14 per-pixel temporal statistics: coverage, gap stats, quantiles, IQR outlier mask, rolling mean/std/slope, z-score anomaly, Mann-Kendall trend test, Theil-Sen robust slope, BFAST structural break, OLS regression with R², phenology (SOS/EOS/peak), Pearson correlation (Fortran) |
| **Covariance** | Per-band covariance matrix — single-pass Kahan-compensated algorithm for PCA, feature reduction, and Mahalanobis anomaly detection (Fortran) |
| **Visualisation** | Interactive Plotly figures: band heatmap, RGB composite, grid, SCL mask |

---

## Installation

```bash
pip install sentinel-processor
```

For NetCDF output:

```bash
pip install "sentinel-processor[netcdf]"
```

All optional extras:

```bash
pip install "sentinel-processor[all]"
```

### Fortran libraries

The Fortran kernels must be compiled once before validation, indices, filters, pansharpening, time-series alignment, and covariance are available. Without them the downloader still works; set `validate=False` and skip Fortran-dependent calls.

**Linux / macOS**
```bash
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/validation/fortran/libsentinel_validation.so \
  sentinel_processor/validation/fortran/validation.f90

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/indices/fortran/libsentinel_indices.so \
  sentinel_processor/indices/fortran/indices_mod.f90

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/processing/fortran/libsentinel_raster_ops.so \
  sentinel_processor/processing/fortran/raster_ops.f90

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/processing/fortran/libsentinel_processing.so \
  sentinel_processor/processing/fortran/pansharpening.f90

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/processing/fortran/libsentinel_timeseries.so \
  sentinel_processor/processing/fortran/timeseries_mod.f90

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/filters/fortran/libsentinel_filters.so \
  sentinel_processor/filters/fortran/filters.f90

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/analysis/fortran/libsentinel_stats.so \
  sentinel_processor/analysis/fortran/sentinel_stats.f90

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/analysis/fortran/libband_covariance.so \
  sentinel_processor/analysis/fortran/band_covariance.f90

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/texture/fortran/libsentinel_texture.so \
  sentinel_processor/texture/fortran/texture_mod.f90
```

**Windows** (MSYS2 UCRT64 — do **not** use `-static-libgfortran` on GCC 16+)

```bat
gfortran -O2 -shared -o sentinel_processor\validation\fortran\libsentinel_validation.dll sentinel_processor\validation\fortran\validation.f90
gfortran -O2 -shared -o sentinel_processor\indices\fortran\libsentinel_indices.dll sentinel_processor\indices\fortran\indices_mod.f90
gfortran -O2 -shared -o sentinel_processor\processing\fortran\libsentinel_raster_ops.dll sentinel_processor\processing\fortran\raster_ops.f90
gfortran -O2 -shared -o sentinel_processor\processing\fortran\libsentinel_processing.dll sentinel_processor\processing\fortran\pansharpening.f90
gfortran -O2 -shared -o sentinel_processor\processing\fortran\libsentinel_timeseries.dll sentinel_processor\processing\fortran\timeseries_mod.f90
gfortran -O2 -shared -o sentinel_processor\filters\fortran\libsentinel_filters.dll sentinel_processor\filters\fortran\filters.f90
gfortran -O2 -shared -o sentinel_processor\analysis\fortran\libsentinel_stats.dll sentinel_processor\analysis\fortran\sentinel_stats.f90
gfortran -O2 -shared -o sentinel_processor\analysis\fortran\libband_covariance.dll sentinel_processor\analysis\fortran\band_covariance.f90
gfortran -O2 -shared -o sentinel_processor\texture\fortran\libsentinel_texture.dll sentinel_processor\texture\fortran\texture_mod.f90
```

After compiling on Windows, copy the MSYS2 runtime DLLs next to each `.dll`:

```bat
for %d in (validation indices processing filters analysis texture) do (
  copy C:\msys64\ucrt64\bin\libgfortran-5.dll    sentinel_processor\%d\fortran\
  copy C:\msys64\ucrt64\bin\libgcc_s_seh-1.dll   sentinel_processor\%d\fortran\
  copy C:\msys64\ucrt64\bin\libwinpthread-1.dll   sentinel_processor\%d\fortran\
)
```

---

## Quick start

```python
import sentinel_processor as sp
from sentinel_processor.indices.compute import compute_indices
from sentinel_processor.filters import apply_filter
from sentinel_processor.input.timeseries import stack_timeseries, TimeSeriesConfig
from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask
from pathlib import Path

# 1. Download
results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=sp.DownloadConfig(
        keep_items  = 10,
        save_report = True,   # enables fast sidecar path in stack_timeseries
    ),
)

scene = "data/spectral/budapest_20260526T095725.nc"
vis   = "data/visual/vis_budapest_20260526T095725.nc"
scl   = "data/technical/scl_budapest_20260526T095725.nc"

# 2. Spectral indices
idx = compute_indices(scene, ["ndvi", "ndwi", "ndbi"])

# 3. Filters
import rioxarray
nir = rioxarray.open_rasterio(scene).sel(band="nir").squeeze().values
nir_smooth = apply_filter(nir, "bilateral",      sigma_s=2.0, sigma_r=0.08)
nir_edges  = apply_filter(nir, "sobel_magnitude")
nir_sharp  = apply_filter(nir, "unsharp_mask",   sigma=1.5, amount=1.2)

# 4. Time series stack
result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg = TimeSeriesConfig(
        max_cloud_fraction = 0.10,
        min_confidence     = 0.75,
        save_dir           = "data/stacks",
    ),
)
print(result.summary())
da = result.stack   # xr.DataArray  (time, band, y, x)  float32

# 5. Gap filling
import numpy as np
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps
from sentinel_processor.analysis._sentinel_stats_bridge import (
    time_window_stats, anomaly_zscore, mann_kendall,
    phenology_metrics, save_phenology, NODATA,
)

arr  = da.values.astype(np.float64)
mask = np.isfinite(arr).astype(np.int32)
filled = np.empty_like(arr)
for b in range(arr.shape[1]):
    filled[:, b] = interpolate_gaps(arr[:, b], mask[:, b], method="pchip")

# 6. Temporal analysis
nir_idx = list(da.coords["band"].values).index("nir")
nir     = filled[:, nir_idx]
dates   = [s.timestamp for s in result.scenes]

mk  = mann_kendall(nir, dates)
bg  = time_window_stats(nir, dates, window_days=90)
z   = anomaly_zscore(nir, dates, background=bg)

# NDVI and phenology
ri      = list(da.coords["band"].values).index("red")
ndvi    = (nir - filled[:, ri]) / (nir + filled[:, ri] + 1e-9) / 10_000.0
metrics = phenology_metrics(ndvi, dates)
save_phenology(metrics, "data/analysis/phenology", crs="EPSG:32633")

# 7. Per-band covariance → PCA
from sentinel_processor.analysis._band_covariance_bridge import band_covariance

stack = filled.transpose(1, 0, 2, 3)   # (n_bands, time, rows, cols) → treat time as pixels
# or pass a single (n_bands, rows, cols) scene cube:
scene_cube = filled[0]                  # (n_bands, rows, cols)
cov = band_covariance(scene_cube)       # (n_bands, n_bands)
eigenvalues, eigenvectors = np.linalg.eigh(cov)

# 8. Visualise
plot_rgb(vis).show()
plot_band(idx["ndvi"], colorscale="RdYlGn").show()
plot_mask(scl).show()
plot_grid([
    {"file": scene,       "band": "nir", "label": "NIR"},
    {"file": idx["ndvi"],                "label": "NDVI", "colorscale": "RdYlGn"},
    {"file": idx["ndwi"],                "label": "NDWI", "colorscale": "Blues"},
], ncols=3).show()
```

---

## Modules

| Module | Description | Docs |
|---|---|---|
| `sentinel_processor` | Download, STAC search, validation | [DOWNLOADER.md](docs/DOWNLOADER.md) |
| `sentinel_processor.validation` | SCL + radiometry quality checks (Fortran) | [VALIDATION.md](docs/VALIDATION.md) |
| `sentinel_processor.indices` | Spectral index computation (Fortran) | [INDICES.md](docs/INDICES.md) |
| `sentinel_processor.filters` | Convolution and morphological filters (Fortran) | [FILTERS.md](docs/FILTERS.md) |
| `sentinel_processor.processing` | Pansharpening + raster ops (Fortran) | [PANSHARPENING.md](docs/PANSHARPENING.md) · [RASTER_OPS.md](docs/RASTER_OPS.md) |
| `sentinel_processor.input.timeseries` | Quality-filtered temporal stack builder | [TIMESERIES.md](docs/TIMESERIES.md) |
| `sentinel_processor.processing._timeseries_bridge` | Gap filling for time stacks (Fortran) | [TIMESERIES.md](docs/TIMESERIES.md#gap-filling) |
| `sentinel_processor.texture` | GLCM texture features: energy, contrast, homogeneity (Fortran) | [TEXTURE.md](docs/TEXTURE.md) |
| `sentinel_processor.analysis` | 14 per-pixel temporal statistics: trend, anomaly, phenology, regression (Fortran) | [ANALYSIS.md](docs/ANALYSIS.md) |
| `sentinel_processor.analysis._band_covariance_bridge` | Per-band covariance matrix for PCA and anomaly detection (Fortran) | [COVARIANCE.md](docs/COVARIANCE.md) |
| `sentinel_processor.visualisation` | Interactive Plotly figures | [VISUALISATION.md](docs/VISUALISATION.md) |

---

## Usage

### Download with defaults

```python
import sentinel_processor as sp

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="my_field")],
)
```

Downloads the 10 most recent cloud-free scenes within ~5 km of the point.
Saves to `data/spectral/`, `data/technical/`, `data/visual/`.

### Custom download config

```python
import datetime, sentinel_processor as sp

cfg = sp.DownloadConfig(
    bands               = sp.SpectralBands.VEGETATION,
    tech_bands          = sp.TechnicalLayers.SCL,
    visual              = False,
    output_dir          = "/mnt/sentinel",
    start_date          = datetime.datetime(2025, 3, 1, tzinfo=datetime.UTC),
    end_date            = datetime.datetime(2025, 6, 1, tzinfo=datetime.UTC),
    keep_items          = 3,
    max_cloud_threshold = 0.15,
    min_confidence      = 0.75,
    save_report         = True,
)
results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=cfg,
)
```

### Time series

Build a quality-filtered temporal stack from any number of downloaded scenes.

```python
from pathlib import Path
from sentinel_processor.input.timeseries import stack_timeseries, TimeSeriesConfig

result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg = TimeSeriesConfig(
        max_cloud_fraction = 0.10,   # same scale as DownloadConfig.max_cloud_threshold
        min_confidence     = 0.75,
        require_bands      = ["red", "nir"],   # reject scenes missing these bands
        save_dir           = "data/stacks",    # auto-save on completion
    ),
)
print(result.summary())    # per-scene quality table
da = result.stack          # xr.DataArray (time, band, y, x) float32
```

**Quality pipeline** — identical thresholds to the downloader:

```
sidecar *_report.json   ← fast path when save_report=True in DownloadConfig
    │  not found ↓
SCL file  →  check_dimensions → validate_scl → check_radiometry  (Fortran)
    │
    └── apply TimeSeriesConfig thresholds → accept or reject scene
```

**Save options:**

```python
# Single NetCDF-4 (default)
result.save("data/stacks", fmt="nc")        # → data/stacks/budapest_stack.nc

# One GeoTIFF per time step
result.save("data/stacks", fmt="tif")       # → data/stacks/budapest_stack_20260520T094746.tif …

# Explicit name
result.save("data/stacks", name="may_2026", fmt="nc")

# Auto-save via config
cfg = TimeSeriesConfig(save_dir="data/stacks", save_format="nc", save_name="may_2026")
```

**Keep rejected scenes as nodata to preserve a contiguous time axis:**

```python
cfg = TimeSeriesConfig(max_cloud_fraction=0.10, fill_rejected=True)
result = stack_timeseries(sources, scl_dir="data/technical", cfg=cfg)
# result.stack.shape[0] == total scenes including rejected (filled with NaN)
```

### Gap filling

Fill cloud-masked holes in the stack using Fortran-accelerated interpolation.

```python
import numpy as np
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps

arr  = result.stack.values.astype(np.float64)  # (time, band, y, x)
mask = np.isfinite(arr).astype(np.int32)        # 1 = valid, 0 = gap

filled = np.empty_like(arr)
for b in range(arr.shape[1]):
    filled[:, b] = interpolate_gaps(arr[:, b], mask[:, b], method="pchip")
```

| `method` | Best for |
|---|---|
| `"linear"` | Short gaps, fast baseline |
| `"savgol"` | Noisy series, preserves peaks (use `window=5..11`) |
| `"pchip"` | NDVI / EVI / LAI — monotone cubic, no overshoot |
| `"ets"` | Series with persistent seasonal trend |
| `"gauss"` | Smooth phenology curves (use `window=9..15`) |

### Temporal analysis

14 Fortran-accelerated per-pixel statistics on a `(n_times, rows, cols)` stack.

```python
import numpy as np
from sentinel_processor.analysis._sentinel_stats_bridge import (
    valid_obs_count, temporal_gap_stats,
    time_window_stats, anomaly_zscore,
    mann_kendall, trend_theil_sen,
    bfast_breakpoint, pixel_regression,
    phenology_metrics, save_phenology,
    pearson_map, NODATA,
)

dates = [s.timestamp for s in result.scenes]

# coverage and data quality
oc = valid_obs_count(nir)
gs = temporal_gap_stats(nir, dates)
print(f"Mean coverage: {oc['fraction'].mean():.1%}  Max gap: {gs['max_gap'].max():.0f} d")

# rolling statistics + per-scene anomaly z-scores
bg = time_window_stats(nir, dates, window_days=90)
z  = anomaly_zscore(nir, dates, background=bg)
drought = (z < -2.0) & (z != NODATA)   # (n_times, rows, cols)

# trend significance and magnitude
mk = mann_kendall(nir, dates)
ts = trend_theil_sen(nir, dates)
print(f"Increasing pixels: {(mk['trend'] == 1.0).sum():,}")

# structural break detection (fire, deforestation)
bp = bfast_breakpoint(nir, dates)
disturbed = (bp["rss_ratio"] < 0.7) & (bp["magnitude"] < -0.05)

# phenology from NDVI (normalised [0, 1])
metrics = phenology_metrics(ndvi, dates, smooth=True)
saved   = save_phenology(metrics, "data/analysis/phenology", crs="EPSG:32633")
```

See [ANALYSIS.md](docs/ANALYSIS.md) for full API reference, method selection guide, and all 14 functions.

### Band covariance

Per-band covariance matrix for PCA-based sharpening, feature reduction, and Mahalanobis anomaly detection.

```python
import numpy as np
from sentinel_processor.analysis._band_covariance_bridge import band_covariance

# cube: (n_bands, rows, cols) float64
cov = band_covariance(cube)                        # (n_bands, n_bands)

# PCA
eigenvalues, eigenvectors = np.linalg.eigh(cov)
explained = eigenvalues / eigenvalues.sum()
print(f"PC1 explains {explained[-1]:.1%} of variance")

# Mahalanobis anomaly detection
prec  = np.linalg.inv(cov)
flat  = cube.reshape(cube.shape[0], -1)            # (n_bands, pixels)
mu    = flat.mean(axis=1, keepdims=True)
delta = flat - mu
d2    = np.einsum("ij,jk,ki->i", delta.T, prec, delta.T.T)
anomaly_map = np.sqrt(d2).reshape(cube.shape[1:])  # (rows, cols)
```

See [COVARIANCE.md](docs/COVARIANCE.md) for full API reference and examples.

### Spectral indices

```python
from sentinel_processor.indices.compute import compute_indices, list_indices

for name, info in list_indices().items():
    print(f"{name:8s}  {info['bands_required']}")

results = compute_indices(
    source="data/spectral/scene.nc",
    indices=["ndvi", "evi", "ndwi", "ndbi", "nbr"],
    output_dir="data/indices",
)
# {"ndvi": "data/indices/indices_scene_ndvi.tif", ...}
```

Supported: `ndvi`, `evi`, `savi`, `ndwi`, `mndwi`, `ndbi`, `nbr`, `ndsi`, `cig`, `arvi`.

### Filters

```python
import numpy as np
from sentinel_processor.filters import apply_filter, apply_filter_da, list_filters

band   = np.random.rand(512, 512)
smooth = apply_filter(band, "gaussian",       sigma=1.5)
clean  = apply_filter(band, "median",         radius=2)
edges  = apply_filter(band, "sobel_magnitude")
sharp  = apply_filter(band, "unsharp_mask",   sigma=1.5, amount=1.2)
feats  = apply_filter(band, "top_hat_white",  radius=5)
custom = apply_filter(band, "convolve",       kernel=np.ones((3, 3)) / 9)

# xarray DataArray — preserves CRS and coordinates
result_da = apply_filter_da(da, "bilateral", sigma_s=3.0, sigma_r=0.08)

# Multi-band (n_bands, rows, cols) — applied per band
ms = np.random.rand(4, 256, 256)
ms_smooth = apply_filter(ms, "gaussian", sigma=1.0)
```

Available: `gaussian`, `bilateral`, `median`, `sobel_magnitude`, `sobel_direction`, `laplacian`, `unsharp_mask`, `erode`, `dilate`, `open`, `close`, `top_hat_white`, `top_hat_black`, `convolve`.

### Pansharpening

```python
import sentinel_processor as sp

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=sp.DownloadConfig(
        bands                = sp.SpectralBands.RGB_NIR,
        pansharpen_algorithm = "gram_schmidt",  # "ihs" (RGB only) | "wavelet"
    ),
)
```

### Validate an SCL file

```python
from sentinel_processor import validate_file

report = validate_file(
    "data/technical/scl_budapest_20260526T095725.tif",
    max_cloud_threshold=0.30,
    min_confidence=0.50,
)
print(report["passed"])            # True
print(report["cloud_ratio"])       # 0.04
print(report["confidence_score"])  # 1.0
```

### Visualise

```python
from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask

plot_band("data/spectral/scene.nc",      band="nir",   colorscale="Plasma").show()
plot_rgb("data/visual/vis_scene.nc").show()
plot_rgb("data/spectral/scene.nc",       "nir", "red", "green").show()  # false colour
plot_mask("data/technical/scl_scene.nc", bad_classes=[8, 9, 10]).show()
plot_grid([
    {"file": "data/spectral/scene.nc", "band": "red",  "label": "Red"},
    {"file": "data/spectral/scene.nc", "band": "nir",  "label": "NIR"},
    {"file": "data/indices/scene_ndvi.tif",             "label": "NDVI", "colorscale": "RdYlGn"},
], ncols=3).show()
```

---

## Output structure

```
<output_dir>/
├── spectral/
│   ├── <name>_<timestamp>.tif
│   ├── <name>_<timestamp>.nc            ← requires [netcdf]
│   ├── <name>_<timestamp>_report.json   ← if save_report=True
│   └── indices/
│       └── indices_<name>_<timestamp>_<index>.tif
├── technical/
│   ├── scl_<name>_<timestamp>.tif/.nc
│   ├── aot_<name>_<timestamp>.tif/.nc
│   └── wvp_<name>_<timestamp>.tif/.nc
├── visual/
│   └── vis_<name>_<timestamp>.tif/.nc
└── stacks/                              ← written by stack_timeseries
    └── <name>_stack.nc / <name>_stack_<timestamp>.tif
analysis_output/                         ← written by analysis functions
    ├── coverage.tif · max_gap.tif · slope_*.tif · r2_*.tif
    ├── mk_trend_*.tif · pearson_*.tif
    └── phenology/
        └── sos_doy.tif · eos_doy.tif · peak_doy.tif · peak_val.tif
```

---

## Band presets

| Preset | Bands |
|---|---|
| `SpectralBands.RGB` | blue, green, red |
| `SpectralBands.RGB_NIR` | blue, green, red, NIR |
| `SpectralBands.VEGETATION` | red, NIR, rededge 1-2-3 |
| `SpectralBands.AGRICULTURE` | blue, green, red, NIR, rededge1, SWIR1, SWIR2 |
| `SpectralBands.ALL_10M` | blue, green, red, NIR |
| `SpectralBands.ALL_20M` | rededge 1-3, NIR narrow, SWIR 1-2 |
| `SpectralBands.ALL` | all 10 m + 20 m bands |
| `TechnicalLayers.SCL` | Scene Classification Layer |
| `TechnicalLayers.ALL` | SCL + AOT + WVP |

Custom list: `DownloadConfig(bands=["red", "nir", "swir16"])`

---

## Project layout

```
sentinel_processor/
├── __init__.py
├── config.py
├── analysis/
│   ├── __init__.py
│   ├── _sentinel_stats_bridge.py
│   ├── _band_covariance_bridge.py
│   └── fortran/
│       ├── sentinel_stats.f90
│       └── band_covariance.f90
├── filters/
│   ├── __init__.py
│   ├── _filters_bridge.py
│   ├── compute.py
│   └── fortran/filters.f90
├── indices/
│   ├── __init__.py
│   ├── _indices_bridge.py
│   ├── compute.py
│   └── fortran/indices_mod.f90
├── input/
│   ├── __init__.py
│   ├── downloader.py
│   └── timeseries.py
├── texture/
│   ├── __init__.py
│   ├── texture.py
│   ├── _texture_bridge.py
│   └── fortran/texture_mod.f90
├── processing/
│   ├── __init__.py
│   ├── _fortran_bridge.py
│   ├── _raster_ops_bridge.py
│   ├── _timeseries_bridge.py
│   └── fortran/pansharpening.f90 · raster_ops.f90 · timeseries_mod.f90
├── utils/
│   ├── __init__.py
│   └── data_utils.py
├── validation/
│   ├── __init__.py
│   ├── _fortran_bridge.py
│   └── fortran/validation.f90
└── visualisation/
    ├── __init__.py
    └── plot.py

tests/
├── conftest.py
├── test_analysis.py
├── test_band_covariance.py
├── test_validation.py
├── test_indices.py
├── test_filters.py
├── test_processing.py
├── test_texture.py
├── test_timeseries.py
├── test_timeseries_bridge.py
└── test_visualisation.py

docs/
├── ANALYSIS.md
├── COVARIANCE.md
├── DOWNLOADER.md
├── FILTERS.md
├── INDICES.md
├── PANSHARPENING.md
├── RASTER_OPS.md
├── TEXTURE.md
├── TIMESERIES.md
├── VALIDATION.md
└── VISUALISATION.md
```

---

## Development

```bash
git clone https://github.com/niki8885/sentinel-processor
cd sentinel-processor
pip install -e ".[dev,netcdf]"
pytest tests/ -v
```

Run only fast unit tests (no I/O, no Fortran required):

```bash
pytest tests/ -m "not integration and not fortran" -v
```

See [CONTRIBUTING.md](.github/CONTRIBUTING.md) for full contribution guidelines.

---

## Requirements

- Python ≥ 3.10
- `pystac-client`, `rioxarray`, `xarray`, `numpy`, `rasterio`, `plotly`
- `netCDF4` or `h5netcdf` for `.nc` output (optional)
- `gfortran` ≥ 9 to build the Fortran kernels

---

## Maintainer

Nikita Manaenkov — [nick.maanenkov@gmail.com](mailto:nick.maanenkov@gmail.com) · [@niki8885](https://github.com/niki8885)