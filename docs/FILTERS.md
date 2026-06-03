# Filters

## Overview

The filters module applies convolution and morphological operations to single-band
raster arrays. All kernels are implemented in Fortran and exposed via
`libsentinel_filters`:

```
input band  (rows × cols, float64)
    │
    ▼
category choice
    ├── gaussian_blur()       ──  separable Gaussian, two-pass O(n·r)   [smoothing]
    ├── bilateral_filter()    ──  edge-preserving spatial + range weight  [smoothing]
    ├── median_filter()       ──  box median, insertion sort              [smoothing]
    ├── sobel_magnitude()     ──  Gx/Gy edge magnitude (L1 or L2)        [edges]
    ├── sobel_direction()     ──  edge orientation in radians             [edges]
    ├── laplacian()           ──  4- or 8-connectivity discrete Laplacian [edges]
    ├── unsharp_mask()        ──  high-boost sharpening with noise gate   [sharpening]
    ├── morpho_erode()        ──  rectangular erosion                     [morphology]
    ├── morpho_dilate()       ──  rectangular dilation                    [morphology]
    ├── morpho_open()         ──  erode → dilate                         [morphology]
    ├── morpho_close()        ──  dilate → erode                         [morphology]
    ├── top_hat_white()       ──  src − opening(src)                     [morphology]
    ├── top_hat_black()       ──  closing(src) − src                     [morphology]
    └── convolve2d()          ──  arbitrary user kernel                  [custom]
    │
    ▼
filtered band  (rows × cols, float64)
```

`apply_filter` is the single entry point for all filters.
`apply_filter_da` wraps xarray DataArrays and preserves CRS and coordinates.

## Install

```bat
# Windows
mkdir "...\sentinel_processor\filters\fortran"
gfortran -O2 -shared -o "...\sentinel_processor\filters\fortran\libsentinel_filters.dll" "...\sentinel_processor\filters\fortran\filters.f90"
```

```bash
# Linux / macOS
gfortran -O2 -march=native -shared -fPIC \
  -o sentinel_processor/filters/fortran/libsentinel_filters.so \
  sentinel_processor/filters/fortran/filters.f90
```

If the library is missing, every call raises `FileNotFoundError` with the
build command printed in the message.

## Module layout

```
sentinel_processor/
└── filters/
    ├── __init__.py           ← public re-exports
    ├── compute.py            ← apply_filter, apply_filter_da, list_filters
    ├── _filters_bridge.py    ← ctypes bridge
    └── fortran/
        ├── filters.f90
        ├── libsentinel_filters.dll / .so
        └── (runtime DLLs on Windows — see Install)
```

## API reference

### apply_filter

```python
from sentinel_processor.filters import apply_filter

result = apply_filter(arr, filter_name, **kwargs)
```

| Parameter | Type | Description |
|---|---|---|
| `arr` | `np.ndarray` | 2-D `(rows, cols)` or 3-D `(n_bands, rows, cols)` float64 array |
| `filter_name` | `str` | One of the names in `AVAILABLE_FILTERS` |
| `**kwargs` | — | Filter-specific parameters (see each filter below) |

Returns `np.ndarray` with the same spatial shape as `arr`, dtype `float64`.
For 3-D input the filter is applied independently to each band.

### apply_filter_da

```python
from sentinel_processor.filters import apply_filter_da

result_da = apply_filter_da(da, filter_name, **kwargs)
```

Wraps `apply_filter` for xarray DataArrays. Preserves all coordinates, dims,
and CRS (via rioxarray when installed). Adds `filter_applied` and parameter
keys to `.attrs`.

### list_filters

```python
from sentinel_processor.filters import list_filters

print(list_filters())
# {'bilateral': {'description': '...', 'default_params': {...}}, ...}
```

Returns a dict of all registered filters with descriptions and default parameters.

### AVAILABLE_FILTERS

```python
from sentinel_processor.filters import AVAILABLE_FILTERS
# ['bilateral', 'close', 'convolve', 'dilate', 'erode', 'gaussian',
#  'laplacian', 'median', 'open', 'sobel_direction', 'sobel_magnitude',
#  'top_hat_black', 'top_hat_white', 'unsharp_mask']
```

## Filters

### gaussian (`gaussian`)

Separable two-pass Gaussian blur. O(n · kradius) — much faster than a 2-D
kernel for large sigma.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sigma` | `float` | `1.0` | Standard deviation in pixels |
| `kradius` | `int \| None` | `ceil(3·sigma)` | Kernel half-width; full kernel = 2·kradius+1 taps |

Pass 1 convolves horizontally, Pass 2 vertically — the result is identical to
a 2-D Gaussian with the same sigma.

### bilateral (`bilateral`)

Edge-preserving smoothing combining a spatial Gaussian weight (sigma_s) with
an intensity range weight (sigma_r). Pixels with large intensity difference
from the centre receive near-zero weight regardless of distance.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sigma_s` | `float` | `2.0` | Spatial sigma in pixels |
| `sigma_r` | `float` | `0.1` | Range sigma (intensity units) |
| `kradius` | `int \| None` | `ceil(2·sigma_s)` | Kernel half-width |

Typical `sigma_r` values: `0.05–0.15` for data in `[0, 1]`; `300–1000` for
raw Sentinel-2 uint16 reflectance `[0, 10000]`.

Cost is O(n · kradius²) — keep kradius ≤ 10 for interactive use.

### median (`median`)

Box median filter. Removes salt-and-pepper noise while preserving edges
better than Gaussian blur.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `radius` | `int` | `1` | Half-width; full window = (2·radius+1)² pixels |

Uses insertion sort on the neighbourhood window — optimal for small windows
(radius ≤ 5). For larger windows consider pre-processing with `gaussian` first.

### sobel_magnitude (`sobel_magnitude`)

Sobel edge magnitude from the Gx and Gy gradient operators.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `norm` | `"l1" \| "l2"` | `"l2"` | `"l2"` → √(Gx²+Gy²); `"l1"` → \|Gx\|+\|Gy\| (faster) |

Output range depends on input scale. Normalise the input to `[0, 1]` first
for consistent thresholding across scenes.

Kernels used:

```
Gx = [-1  0 +1]    Gy = [-1 -2 -1]
     [-2  0 +2]         [ 0  0  0]
     [-1  0 +1]         [+1 +2 +1]
```

### sobel_direction (`sobel_direction`)

Edge orientation angle in radians computed as `atan2(Gy, Gx)`.
Output is in `[-π, π]`.

Useful for line detection and texture orientation analysis.
No parameters.

### laplacian (`laplacian`)

Discrete Laplacian — second-order derivative. Highlights rapid intensity
changes; output is positive on bright edges, negative on dark edges.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `connectivity` | `4 \| 8` | `4` | 4 → cross kernel `(0,1,0 / 1,−4,1 / 0,1,0)`; 8 → full 3×3 kernel with centre −8 |

The 8-connectivity version responds to diagonal edges; use 4 for isotropic
response on axis-aligned features.

### unsharp_mask (`unsharp_mask`)

High-boost sharpening via residual injection:

```
dst = src + amount × (src − gaussian_blur(src, sigma))
```

A noise gate suppresses sharpening where the residual is below `threshold`,
preventing amplification of sensor noise.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sigma` | `float` | `1.0` | Blur sigma for the low-frequency estimate |
| `amount` | `float` | `1.0` | Sharpening strength; typical `0.5–2.0` |
| `kradius` | `int \| None` | `ceil(3·sigma)` | Gaussian kernel half-width |
| `threshold` | `float` | `0.0` | Minimum \|residual\| to apply sharpening |

### erode (`erode`)

Morphological erosion with a rectangular structuring element.
Replaces each pixel with the minimum of its neighbourhood — shrinks bright
regions, removes small bright blobs.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `radius` | `int` | `1` | Half-width of the structuring element |

### dilate (`dilate`)

Morphological dilation with a rectangular structuring element.
Replaces each pixel with the maximum of its neighbourhood — expands bright
regions, fills small dark holes.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `radius` | `int` | `1` | Half-width of the structuring element |

### open (`open`)

Morphological opening: `erode` then `dilate`. Removes small bright blobs
smaller than the structuring element without shifting larger structures.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `radius` | `int` | `1` | Half-width of the structuring element |

### close (`close`)

Morphological closing: `dilate` then `erode`. Fills small dark holes
smaller than the structuring element.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `radius` | `int` | `1` | Half-width of the structuring element |

### top_hat_white (`top_hat_white`)

White top-hat transform: `dst = src − opening(src)`.

Extracts small bright features (thinner than 2·radius pixels) against a
varying background — useful for road extraction, building detection, and
removing large-scale brightness gradients before thresholding.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `radius` | `int` | `3` | Structuring element half-width |

### top_hat_black (`top_hat_black`)

Black top-hat transform: `dst = closing(src) − src`.

Extracts small dark features (valleys, shadows, water channels) against a
varying background.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `radius` | `int` | `3` | Structuring element half-width |

### convolve (`convolve`)

Arbitrary 2-D convolution with a user-supplied kernel. Zero-padding at borders.

| Parameter | Type | Description |
|---|---|---|
| `kernel` | `np.ndarray` | 2-D array of shape `(krows, kcols)`. Both dimensions should be odd for symmetric padding. The kernel is applied as-is (no flipping). |

Both `krows` and `kcols` can be any positive integer.

## Array layout

All arrays passed to Fortran are converted to **column-major (Fortran) order**
before the ctypes call. Input arrives as a C-order `(rows, cols)` numpy array
and is transposed to Fortran order by `_to_f_f64`. The output buffer is
allocated in Fortran order and converted back to C-contiguous by `_out_2d`
before being returned.

For 3-D multi-band arrays `(n_bands, rows, cols)`, `apply_to_bands` iterates
over the band axis and calls the Fortran routine once per band. No band axis
is ever passed to Fortran directly — all Fortran subroutines operate on flat
2-D single-band arrays.

## Examples

### Single band

```python
import numpy as np
from sentinel_processor.filters import apply_filter

band = ...  # np.ndarray shape (rows, cols), e.g. from ds["nir"].values

# Gaussian blur
blurred = apply_filter(band, "gaussian", sigma=1.5)

# Sobel edge magnitude
edges = apply_filter(band, "sobel_magnitude")

# Median — salt-and-pepper removal
clean = apply_filter(band, "median", radius=2)

# Bilateral — smooth noise, keep field boundaries sharp
smooth = apply_filter(band, "bilateral", sigma_s=3.0, sigma_r=0.08)

# Unsharp mask — increase apparent resolution
sharp = apply_filter(band, "unsharp_mask", sigma=1.5, amount=1.2, threshold=0.01)

# White top-hat — extract small bright objects on varying background
feats = apply_filter(band, "top_hat_white", radius=5)
```

### Multi-band stack

```python
import numpy as np
from sentinel_processor.filters import apply_filter

ms = ...  # np.ndarray shape (n_bands, rows, cols)

# Applied independently to each band
ms_blurred = apply_filter(ms, "gaussian", sigma=1.5)
ms_edges   = apply_filter(ms, "sobel_magnitude")
```

### xarray DataArray — preserves CRS and coordinates

```python
import xarray as xr
from sentinel_processor.filters import apply_filter_da

da = ...  # xarray.DataArray from rioxarray.open_rasterio(...)

sharpened_da = apply_filter_da(da, "unsharp_mask", sigma=1.5, amount=1.2)
sharpened_da.rio.to_raster("sharpened_nir.tif")
```

### Arbitrary kernel

```python
import numpy as np
from sentinel_processor.filters import apply_filter

# 3×3 mean (box blur)
kernel = np.ones((3, 3)) / 9.0
result = apply_filter(band, "convolve", kernel=kernel)

# Horizontal Prewitt edge
kernel = np.array([[-1, -1, -1],
                   [ 0,  0,  0],
                   [ 1,  1,  1]], dtype=np.float64)
result = apply_filter(band, "convolve", kernel=kernel)

# 5×5 sharpening kernel
kernel = np.array([[ 0,  0, -1,  0,  0],
                   [ 0, -1, -2, -1,  0],
                   [-1, -2, 17, -2, -1],
                   [ 0, -1, -2, -1,  0],
                   [ 0,  0, -1,  0,  0]], dtype=np.float64) / 5.0
result = apply_filter(band, "convolve", kernel=kernel)
```

### Full pipeline with indices

```python
import sentinel_processor as sp
from sentinel_processor.indices.compute import compute_indices
from sentinel_processor.filters import apply_filter_da
import xarray as xr
import rioxarray

scene = "data/spectral/budapest_20260528T094727.nc"

# 1. Compute NDVI
idx_paths = compute_indices(scene, ["ndvi"])

# 2. Open NDVI raster
ndvi_da = rioxarray.open_rasterio(idx_paths["ndvi"]).squeeze()

# 3. Denoise with bilateral filter
ndvi_clean = apply_filter_da(ndvi_da, "bilateral", sigma_s=2.0, sigma_r=0.05)

# 4. Extract vegetation patches via white top-hat
veg_patches = apply_filter_da(ndvi_clean, "top_hat_white", radius=7)

# 5. Save
ndvi_clean.rio.to_raster("data/filtered/ndvi_bilateral.tif")
veg_patches.rio.to_raster("data/filtered/ndvi_tophat.tif")
```

### Quick visual test

```python
# test_filters.py — run from project root
# Opens budapest scene, applies all 16 filter variants,
# saves interactive HTML grid and opens it in the browser.
from pathlib import Path
import numpy as np, xarray as xr, plotly.graph_objects as go
from plotly.subplots import make_subplots
from sentinel_processor.filters import apply_filter

NC_FILE  = Path(r"data\spectral\budapest_20260528T094727.nc")
OUT_HTML = NC_FILE.parent.parent / "vis" / "filter_test.html"
OUT_HTML.parent.mkdir(parents=True, exist_ok=True)

ds = xr.open_dataset(NC_FILE)
# ... (see test_filters.py in project root for full loader)
```

See `test_filters.py` in the project root for the complete runnable script
that produced the visualisation below.

## Filter selection guide

| Goal | Recommended filter | Key parameters |
|---|---|---|
| Remove sensor noise before index computation | `bilateral` | `sigma_s=2`, `sigma_r=0.05–0.1` |
| Remove salt-and-pepper artefacts | `median` | `radius=1–2` |
| Smooth before morphological ops | `gaussian` | `sigma=1.0–2.0` |
| Detect field/urban boundaries | `sobel_magnitude` | `norm="l2"` |
| Detect thin linear features (roads, rivers) | `laplacian` | `connectivity=8` |
| Increase apparent spatial detail | `unsharp_mask` | `sigma=1.5, amount=1.0–2.0` |
| Extract small bright objects (buildings, boats) | `top_hat_white` | `radius=3–7` |
| Extract water channels, shadows | `top_hat_black` | `radius=3–7` |
| Remove large-scale illumination gradient | `top_hat_white` | `radius=15–30` |
| Custom edge / texture kernel | `convolve` | `kernel=<ndarray>` |