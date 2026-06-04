# sentinel-processor

Sentinel-2 L2A downloader and processing toolkit built on the [Element84 STAC API](https://earth-search.aws.element84.com/v1).

Downloads spectral bands, quality layers, and visual overviews for any coordinate. Validation, spectral index computation, convolution filters, and pansharpening are all backed by compiled Fortran kernels — Python handles I/O and orchestration, Fortran handles the pixels.

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

The Fortran kernels must be compiled once before validation, indices, filters, and pansharpening are available. Without them, set `validate=False` and skip `compute_indices` / `apply_filter` — the downloader works normally regardless.

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
  -o sentinel_processor/filters/fortran/libsentinel_filters.so \
  sentinel_processor/filters/fortran/filters.f90
```

**Windows** (MSYS2 UCRT64 — do **not** use `-static-libgfortran` on GCC 16+)

```bat
gfortran -O2 -shared -o sentinel_processor\validation\fortran\libsentinel_validation.dll sentinel_processor\validation\fortran\validation.f90
gfortran -O2 -shared -o sentinel_processor\indices\fortran\libsentinel_indices.dll sentinel_processor\indices\fortran\indices_mod.f90
gfortran -O2 -shared -o sentinel_processor\processing\fortran\libsentinel_raster_ops.dll sentinel_processor\processing\fortran\raster_ops.f90
gfortran -O2 -shared -o sentinel_processor\processing\fortran\libsentinel_processing.dll sentinel_processor\processing\fortran\pansharpening.f90
gfortran -O2 -shared -o sentinel_processor\filters\fortran\libsentinel_filters.dll sentinel_processor\filters\fortran\filters.f90
```

After compiling on Windows, copy the MSYS2 runtime DLLs next to each `.dll`:

```bat
for %d in (validation indices processing filters) do (
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
from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask

# 1. Download
results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
)

scene = "data/spectral/budapest_20260526T095725.nc"
vis   = "data/visual/vis_budapest_20260526T095725.nc"
scl   = "data/technical/scl_budapest_20260526T095725.nc"

# 2. Spectral indices
idx = compute_indices(scene, ["ndvi", "ndwi", "ndbi"])

# 3. Filters
import rioxarray
nir = rioxarray.open_rasterio(scene).sel(band="nir").squeeze().values
nir_smooth  = apply_filter(nir, "bilateral",     sigma_s=2.0, sigma_r=0.08)
nir_edges   = apply_filter(nir, "sobel_magnitude")
nir_sharp   = apply_filter(nir, "unsharp_mask",  sigma=1.5, amount=1.2)

# 4. Visualise
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
)
results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=cfg,
)
```

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
│   ├── <name>_<timestamp>_report.json
│   └── indices/
│       └── indices_<name>_<timestamp>_<index>.tif
├── technical/
│   ├── scl_<name>_<timestamp>.tif/.nc
│   ├── aot_<name>_<timestamp>.tif/.nc
│   └── wvp_<name>_<timestamp>.tif/.nc
└── visual/
    └── vis_<name>_<timestamp>.tif/.nc
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
│   └── downloader.py
├── processing/
│   ├── _fortran_bridge.py
│   ├── _raster_ops_bridge.py
│   └── fortran/pansharpening.f90 · raster_ops.f90
├── utils/
│   └── data_utils.py
├── validation/
│   ├── _fortran_bridge.py
│   └── fortran/validation.f90
└── visualisation/
    └── plot.py

tests/
├── conftest.py
├── test_validation.py
├── test_indices.py
├── test_filters.py
├── test_processing.py
└── test_visualisation.py

docs/
├── DOWNLOADER.md
├── VALIDATION.md
├── INDICES.md
├── FILTERS.md
├── PANSHARPENING.md
├── RASTER_OPS.md
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
