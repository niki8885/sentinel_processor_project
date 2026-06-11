# sentinel-processor

Sentinel-2 L2A downloader and processing toolkit built on the [Element84 STAC API](https://earth-search.aws.element84.com/v1).

Downloads spectral bands, quality layers, and visual overviews for any coordinate. Validation, spectral index computation, convolution filters, pansharpening, wavelet denoising, time-series stacking, and band covariance are all backed by compiled Fortran kernels — Python handles I/O and orchestration, Fortran handles the pixels.

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
| **Wavelet** | Multi-level 2-D/3-D DWT with 6 orthogonal wavelets, BayesShrink denoising, batch multi-spectral processing, sub-band energy and statistics (Fortran) |
| **Time series** | Quality-filtered temporal stack builder with cloud/snow filtering, alignment, and save (Fortran validation + raster ops) |
| **Gap filling** | Fill cloud-masked holes in time stacks: linear, Savitzky-Golay, PCHIP, Holt ETS, Gaussian (Fortran) |
| **Texture** | GLCM texture features per pixel: energy, contrast, homogeneity — window, distance, and angle-configurable (Fortran) |
| **Analysis** | 14 per-pixel temporal statistics: coverage, gap stats, quantiles, IQR outlier mask, rolling mean/std/slope, z-score anomaly, Mann-Kendall trend test, Theil-Sen robust slope, BFAST structural break, OLS regression with R², phenology (SOS/EOS/peak), Pearson correlation (Fortran) |
| **Covariance** | Per-band covariance matrix — single-pass Kahan-compensated algorithm for PCA, feature reduction, and Mahalanobis anomaly detection (Fortran) |
| **Visualisation** | Interactive Plotly figures: band heatmap, RGB composite, multi-panel grid, SCL mask, pixel/region time series with cloud markers |

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

The Fortran kernels must be compiled once before validation, indices, filters, pansharpening, wavelet, time-series alignment, and covariance are available. Without them the downloader still works; set `validate=False` and skip Fortran-dependent calls.

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

gfortran -O2 -shared -fPIC \
  -o sentinel_processor/wavelet/fortran/libsentinel_wavelet.so \
  sentinel_processor/wavelet/fortran/wavelet_mod.f90
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
gfortran -O2 -shared -o sentinel_processor\wavelet\fortran\libsentinel_wavelet.dll sentinel_processor\wavelet\fortran\wavelet_mod.f90
```

After compiling on Windows, copy the MSYS2 runtime DLLs next to each `.dll`:

```bat
for %d in (validation indices processing filters analysis texture wavelet) do (
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
from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask, plot_timeseries
from pathlib import Path

# 1. Download
results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=sp.DownloadConfig(
        keep_items  = 10,
        save_report = True,
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

# 4. Wavelet denoising
import numpy as np
from sentinel_processor.wavelet._wavelet_bridge import (
    dwt2d, idwt2d, threshold_coeffs,
    bayes_denoise,
    dwt2d_batch, idwt2d_batch,
    dwt3d, idwt3d,
    band_energy, band_stats,
)

# One-liner BayesShrink denoising
nir_f64    = nir.astype(np.float64)
nir_clean  = bayes_denoise(nir_f64, levels=3, wavelet="db4")

# Manual control: forward DWT → inspect sub-bands → threshold → inverse
coeffs = dwt2d(nir_f64, levels=3, wavelet="sym4")
for lv in sorted(coeffs):
    e = band_energy(coeffs[lv]["HH"])
    s = band_stats(coeffs[lv]["HH"])
    print(f"L{lv} HH  energy={e:.2f}  linf={s['linf']:.4f}")
thresholded = threshold_coeffs(coeffs, threshold=0.02, mode="soft")
nir_denoised = idwt2d(thresholded, wavelet="sym4")

# Batch denoising for a multi-spectral stack (n_bands, rows, cols)
ms = rioxarray.open_rasterio(scene).values.astype(np.float64)  # (n_bands, rows, cols)
coeffs_list  = dwt2d_batch(ms, levels=2, wavelet="db4")
recon_stack  = idwt2d_batch(coeffs_list, wavelet="db4")

# 3-D DWT on a temporal cube  (n_times, rows, cols)
result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg = TimeSeriesConfig(
        max_cloud_fraction = 0.10,
        min_confidence     = 0.75,
        save_dir           = "data/stacks",
    ),
)
da = result.stack   # xr.DataArray  (time, band, y, x)  float32
nir_idx  = list(da.coords["band"].values).index("nir")
nir_cube = da.values[:, nir_idx].astype(np.float64)  # (n_times, rows, cols)

# Apply separable 3-D DWT — spatial + temporal at once
# n_times must be a power of 2; pad if needed
coeffs_vol    = dwt3d(nir_cube, levels=2, wavelet="haar")
nir_cube_back = idwt3d(coeffs_vol, wavelet="haar")

# 5. Gap filling
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps
from sentinel_processor.analysis._sentinel_stats_bridge import (
    time_window_stats, anomaly_zscore, mann_kendall,
    phenology_metrics, save_phenology, NODATA,
)

arr    = da.values.astype(np.float64)
mask   = np.isfinite(arr).astype(np.int32)
filled = np.empty_like(arr)
for b in range(arr.shape[1]):
    filled[:, b] = interpolate_gaps(arr[:, b], mask[:, b], method="pchip")

# 6. Temporal analysis
nir    = filled[:, nir_idx]
dates  = [s.timestamp for s in result.scenes]

mk  = mann_kendall(nir, dates)
bg  = time_window_stats(nir, dates, window_days=90)
z   = anomaly_zscore(nir, dates, background=bg)

ri      = list(da.coords["band"].values).index("red")
ndvi    = (nir - filled[:, ri]) / (nir + filled[:, ri] + 1e-9) / 10_000.0
metrics = phenology_metrics(ndvi, dates)
save_phenology(metrics, "data/analysis/phenology", crs="EPSG:32633")

# 7. Per-band covariance → PCA
from sentinel_processor.analysis._band_covariance_bridge import band_covariance

scene_cube  = filled[0]            # (n_bands, rows, cols)
cov         = band_covariance(scene_cube)
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

## Modules

| Module | Description | Docs |
|---|---|---|
| `sentinel_processor` | Download, STAC search, validation | [DOWNLOADER.md](docs/DOWNLOADER.md) |
| `sentinel_processor.indices` | Spectral index computation | [INDICES.md](docs/INDICES.md) |
| `sentinel_processor.filters` | Spatial filters | [FILTERS.md](docs/FILTERS.md) |
| `sentinel_processor.wavelet` | Wavelet DWT, denoising, sub-band analysis | [WAVELET.md](docs/WAVELET.md) |
| `sentinel_processor.processing` | Pansharpening, raster ops, gap filling | [PANSHARPENING.md](docs/PANSHARPENING.md) |
| `sentinel_processor.input.timeseries` | Time-series stacking | [TIMESERIES.md](docs/TIMESERIES.md) |
| `sentinel_processor.analysis` | Temporal stats, phenology, trends | [ANALYSIS.md](docs/ANALYSIS.md) |
| `sentinel_processor.analysis` (covariance) | Band covariance, PCA support | [COVARIANCE.md](docs/COVARIANCE.md) |
| `sentinel_processor.texture` | GLCM texture features | [TEXTURE.md](docs/TEXTURE.md) |
| `sentinel_processor.visualisation` | Plotly visualisation | [VISUALISATION.md](docs/VISUALISATION.md) |

---

## Code examples

### Spectral indices

```python
from sentinel_processor.indices.compute import compute_indices, list_indices

print(list_indices())
# ['arvi', 'cig', 'evi', 'mndwi', 'nbr', 'ndbi', 'ndsi', 'ndvi', 'ndwi', 'savi']

paths = compute_indices("data/spectral/scene.nc", ["ndvi", "evi", "ndwi"])
# paths == {'ndvi': Path('...'), 'evi': Path('...'), 'ndwi': Path('...')}
```

### Filters

```python
import numpy as np
from sentinel_processor.filters import apply_filter, apply_filter_da
import xarray as xr

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

### Wavelet

```python
import numpy as np
from sentinel_processor.wavelet._wavelet_bridge import (
    dwt2d, idwt2d, threshold_coeffs,
    estimate_sigma, bayes_threshold, bayes_denoise,
    band_energy, band_stats,
    dwt2d_batch, idwt2d_batch,
    dwt3d, idwt3d,
)

band = np.random.rand(256, 256)

# --- BayesShrink denoising (one-liner) ---
denoised = bayes_denoise(band, levels=3, wavelet="db4")

# --- Manual pipeline ---
coeffs = dwt2d(band, levels=3, wavelet="sym4")
sigma_n = estimate_sigma(coeffs[1]["HH"])      # noise from finest HH sub-band
for lv in sorted(coeffs):
    for key in ("LH", "HL", "HH"):
        thr = bayes_threshold(coeffs[lv][key], sigma_n)
        coeffs[lv][key] = np.sign(coeffs[lv][key]) * np.maximum(
            np.abs(coeffs[lv][key]) - thr, 0.0
        )
denoised = idwt2d(coeffs, wavelet="sym4")

# --- Sub-band statistics ---
e = band_energy(coeffs[1]["HH"])
s = band_stats(coeffs[2]["LH"])  # {'mean', 'var', 'l1', 'linf'}

# --- Multi-spectral batch ---
ms = np.random.rand(6, 256, 256)                    # (n_bands, rows, cols)
cl = dwt2d_batch(ms, levels=2, wavelet="db4")       # list of dicts, one per band
ms_back = idwt2d_batch(cl, wavelet="db4")           # (n_bands, rows, cols)

# --- 3-D (spatial + temporal) ---
cube = np.random.rand(8, 64, 64)                    # (n_times, rows, cols)
c3d  = dwt3d(cube, levels=2, wavelet="haar")
back = idwt3d(c3d, wavelet="haar")
```

Supported wavelets: `haar`, `db4`, `db6`, `coif1`, `sym4`, `sym6`. All orthogonal — perfect reconstruction to float64 machine epsilon.

See [docs/WAVELET.md](docs/WAVELET.md) for the full API reference and wavelet selection guide.

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
from sentinel_processor.visualisation.plot import (
    plot_band, plot_rgb, plot_grid, plot_mask, plot_timeseries,
)

plot_band("data/spectral/scene.nc",      band="nir",   colorscale="Plasma").show()
plot_rgb("data/visual/vis_scene.nc").show()
plot_rgb("data/spectral/scene.nc",       "nir", "red", "green").show()  # false colour
plot_mask("data/technical/scl_scene.nc", bad_classes=[8, 9, 10]).show()
plot_grid([
    {"file": "data/spectral/scene.nc", "band": "red",  "label": "Red"},
    {"file": "data/spectral/scene.nc", "band": "nir",  "label": "NIR"},
    {"file": "data/indices/scene_ndvi.tif",             "label": "NDVI", "colorscale": "RdYlGn"},
], ncols=3).show()

# Time series — pixel or region, with optional cloud markers
plot_timeseries(
    stack="data/stacks/budapest_stack.nc",
    lon=19.17, lat=47.56,
    bands=["ndvi", "evi"],
    scl_path="data/technical/",   # directory → per-scene cloud markers
    agg_bbox=0.002,               # ~200 m spatial average
    save_html="data/vis/timeseries.html",
).show()
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
├── visualisation/
│   ├── __init__.py
│   └── plot.py
└── wavelet/                        ← NEW
    ├── __init__.py
    ├── _wavelet_bridge.py
    └── fortran/
        └── wavelet_mod.f90

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
├── test_visualisation.py
└── test_wavelet.py                 ← NEW

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
├── VISUALISATION.md
└── WAVELET.md                      ← NEW
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