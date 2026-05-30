# Raster Ops

## Overview

`raster_ops.f90` contains fast Fortran routines for CPU-bound raster operations
that previously ran in Python or through the GDAL/rioxarray stack. They are called
from the downloader after network I/O completes and before files are saved.

All routines are compiled into `libsentinel_raster_ops` and wrapped by
`processing/_raster_ops_bridge.py`.

```
downloaded band (numpy array)
    │
    ▼
_raster_ops_bridge.py  ──  convert to Fortran-order float64, pass pointer
    │
    ▼
libsentinel_raster_ops  ──  Fortran routine executes
    │
    ▼
result returned as C-contiguous numpy array
```

Every public function has a **NumPy fallback** — if the `.dll`/`.so` is not found,
the same result is produced in pure Python/NumPy without raising an error.

## Install

```bat
# Windows (MSYS2 UCRT64 terminal or any shell with gfortran on PATH)
gfortran -O2 -shared -o "...\sentinel_processor\processing\fortran\libsentinel_raster_ops.dll" "...\sentinel_processor\processing\fortran\raster_ops.f90"
```

```bash
# Linux / macOS
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/processing/fortran/libsentinel_raster_ops.so \
  sentinel_processor/processing/fortran/raster_ops.f90
```

## Module layout

```
sentinel_processor/
└── processing/
    ├── _raster_ops_bridge.py
    └── fortran/
        ├── raster_ops.f90
        └── libsentinel_raster_ops.dll / .so
```

## Array layout convention

All arrays are passed as **Fortran-order (column-major) float64** pointers.
The bridge converts automatically:

```python
# 3-D input (n_bands, rows, cols)  →  transpose to (rows, cols, n_bands) + order='F'
src_f = np.asfortranarray(arr.transpose(1, 2, 0), dtype=np.float64)

# 2-D input (rows, cols)  →  order='F' in-place
src_f = np.asfortranarray(arr, dtype=np.float64)

# Output buffer allocated Fortran-order, transposed back after call
out_f = np.empty(shape, dtype=np.float64, order='F')
```

Indices inside Fortran are 1-based; the flat offset for pixel `(r, c)` is
`(c-1)*rows + r`.

## Fortran routines

### rgb_to_luminance

Collapse a multi-band array to a single luminance plane.

**Signature (Fortran)**
```fortran
subroutine rgb_to_luminance(src, rows, cols, n_bands, dst)
```

| Argument | Type | Description |
|---|---|---|
| `src` | `real(c_double)(rows*cols*n_bands)` | Input, band-major flat column-major |
| `rows`, `cols` | `integer(c_int)` | Spatial dimensions |
| `n_bands` | `integer(c_int)` | Number of input bands |
| `dst` | `real(c_double)(rows*cols)` | Output luminance plane |

**Weights**

| `n_bands` | Method |
|---|---|
| 1 | Direct copy |
| 3 | Rec. 601: `Y = 0.299·R + 0.587·G + 0.114·B` |
| other | Equal weights: `1 / n_bands` per band |

**Python call**
```python
from sentinel_processor.processing._raster_ops_bridge import rgb_to_luminance

lum = rgb_to_luminance(arr)   # arr: (n_bands, rows, cols)
# returns: (rows, cols), float64
```

---

### align_bands

Nearest-neighbour resample a single band to a target pixel grid.
Used when two bands cover the same geographic extent but have different
pixel counts (e.g. 20 m band aligned to 10 m reference).

**Signature (Fortran)**
```fortran
subroutine align_bands(src, src_rows, src_cols, dst, dst_rows, dst_cols)
```

| Argument | Type | Description |
|---|---|---|
| `src` | `real(c_double)(src_rows*src_cols)` | Input band, column-major |
| `src_rows`, `src_cols` | `integer(c_int)` | Source dimensions |
| `dst` | `real(c_double)(dst_rows*dst_cols)` | Output band, column-major |
| `dst_rows`, `dst_cols` | `integer(c_int)` | Target dimensions |

Scale factors: `scale_r = src_rows / dst_rows`, `scale_c = src_cols / dst_cols`.
Each destination pixel maps to `src[floor((r-1)*scale_r), floor((c-1)*scale_c)]`.

**Python call**
```python
from sentinel_processor.processing._raster_ops_bridge import align_bands

out = align_bands(src, dst_rows=1100, dst_cols=1100)
# src: (rows, cols) 2-D array
# returns: (dst_rows, dst_cols), float64
```

---

### clip_box_indices

Compute 1-based pixel index bounds for a lon/lat bounding box given
a raster's affine parameters. Pure arithmetic — no array data is read.

**Signature (Fortran)**
```fortran
subroutine clip_box_indices(origin_x, origin_y, pixel_w, pixel_h, &
    rows, cols, min_lon, max_lon, min_lat, max_lat, &
    row_min, row_max, col_min, col_max)
```

| Argument | Type | Description |
|---|---|---|
| `origin_x`, `origin_y` | `real(c_double)` | Top-left corner coordinates |
| `pixel_w` | `real(c_double)` | Pixel width (positive) |
| `pixel_h` | `real(c_double)` | Pixel height (negative for north-up) |
| `rows`, `cols` | `integer(c_int)` | Full raster dimensions |
| `min_lon` … `max_lat` | `real(c_double)` | Bounding box |
| `row_min`, `row_max`, `col_min`, `col_max` | `integer(c_int)` | Output indices, 1-based, clamped |

**Python call**
```python
from sentinel_processor.processing._raster_ops_bridge import clip_box_indices

row_min, row_max, col_min, col_max = clip_box_indices(
    origin_x=600000.0, origin_y=5500000.0,
    pixel_w=10.0,       pixel_h=-10.0,
    rows=10980,         cols=10980,
    min_lon=19.12,      max_lon=19.22,
    min_lat=47.51,      max_lat=47.61,
)
```

---

### band_stats

Compute mean, std, min, max in a **single pass** over a flat array.
Faster than calling `np.mean`, `np.std`, `np.min`, `np.max` separately
because it avoids four separate traversals of the data.

**Signature (Fortran)**
```fortran
subroutine band_stats(arr, n, mean_out, std_out, min_out, max_out)
```

| Argument | Type | Description |
|---|---|---|
| `arr` | `real(c_double)(n)` | Input array (flat) |
| `n` | `integer(c_int)` | Element count |
| `mean_out`, `std_out`, `min_out`, `max_out` | `real(c_double)` | Output scalars |

Standard deviation uses the population formula with a `1e-12` stability term.

**Python call**
```python
from sentinel_processor.processing._raster_ops_bridge import band_stats

stats = band_stats(arr)   # arr: any shape, flattened internally
# returns: {"mean": ..., "std": ..., "min": ..., "max": ...}
```

---

### normalize_band

Linear stretch to `[out_min, out_max]` **in-place**. No intermediate array
is allocated — operates directly on the caller's buffer.

**Signature (Fortran)**
```fortran
subroutine normalize_band(arr, n, src_min, src_max, out_min, out_max)
```

| Argument | Type | Description |
|---|---|---|
| `arr` | `real(c_double)(n)` | In/out array |
| `src_min`, `src_max` | `real(c_double)` | Input range |
| `out_min`, `out_max` | `real(c_double)` | Output range |

If `abs(src_max - src_min) < 1e-12` (flat band), the output is filled with
`out_min`.

**Python call**
```python
from sentinel_processor.processing._raster_ops_bridge import normalize_band

arr = normalize_band(arr)                          # auto src_min/max, out [0,1]
arr = normalize_band(arr, src_min=0, src_max=10000, out_min=0.0, out_max=1.0)
```

---

### reproject_nearest

Full nearest-neighbour reproject between two grids described by their affine
parameters. Replaces `rioxarray.reproject_match` for the same-CRS case —
no GDAL warp overhead, direct pixel-address arithmetic.

**Signature (Fortran)**
```fortran
subroutine reproject_nearest(src, src_rows, src_cols, src_ox, src_oy, src_pw, src_ph, &
    dst, dst_rows, dst_cols, dst_ox, dst_oy, dst_pw, dst_ph)
```

| Argument | Type | Description |
|---|---|---|
| `src` | `real(c_double)(src_rows*src_cols)` | Source raster, column-major |
| `src_rows`, `src_cols` | `integer(c_int)` | Source dimensions |
| `src_ox`, `src_oy` | `real(c_double)` | Source top-left corner |
| `src_pw`, `src_ph` | `real(c_double)` | Source pixel size (pw > 0, ph < 0) |
| `dst` | `real(c_double)(dst_rows*dst_cols)` | Output raster, column-major |
| `dst_*` | — | Same fields for destination grid |

For each destination pixel `(dr, dc)`:
```
lon = dst_ox + (dc-1) * dst_pw
lat = dst_oy + (dr-1) * dst_ph
sc  = round((lon - src_ox) / src_pw) + 1   (clamped)
sr  = round((lat - src_oy) / src_ph) + 1   (clamped)
dst[dr, dc] = src[sr, sc]
```

**Python call**
```python
from sentinel_processor.processing._raster_ops_bridge import reproject_nearest

out = reproject_nearest(
    src        = band_20m,                        # (rows_20m, cols_20m)
    src_affine = (ox, oy, 20.0, -20.0),
    dst_rows   = 1100,
    dst_cols   = 1100,
    dst_affine = (ox, oy, 10.0, -10.0),
)
# returns: (1100, 1100), float64
```

**When it is used automatically**

In the downloader, `reproject_nearest` replaces `reproject_match` when:
- Source and reference CRS are identical
- Both grids have uniform, extractable affine parameters
- Source band is 2-D after squeezing

Otherwise the code falls back to `rioxarray.reproject_match`.

## Where each routine is called in the downloader

| Location in `downloader.py` | Routine | Replaces |
|---|---|---|
| `_clip()` — same-CRS alignment to reference grid | `reproject_nearest` | `da.rio.reproject_match(reference_da)` |
| `_download_item()` — band stack alignment loop | `reproject_nearest` | `da.rio.reproject_match(reference_da)` / `da.interp(..., method="nearest")` |
| `_apply_pansharpening()` — PAN preparation | `rgb_to_luminance` | `np.tensordot(weights, pan_np, ...)` |

## Low-level example

```python
import numpy as np
from sentinel_processor.processing._raster_ops_bridge import (
    rgb_to_luminance,
    align_bands,
    band_stats,
    normalize_band,
    reproject_nearest,
)

# Collapse visual RGB → luminance PAN
visual = np.random.randint(0, 10000, (3, 1100, 1100)).astype(np.float64)
pan = rgb_to_luminance(visual)        # (1100, 1100)

# Align a 20 m band to a 10 m grid
band_20m = np.random.rand(550, 550)
band_10m = align_bands(band_20m, dst_rows=1100, dst_cols=1100)

# Stats in one pass
stats = band_stats(band_10m)
print(stats)  # {"mean": ..., "std": ..., "min": ..., "max": ...}

# Normalize to [0, 1]
norm = normalize_band(band_10m.copy(), out_min=0.0, out_max=1.0)

# Reproject with full affine info
out = reproject_nearest(
    src        = band_20m,
    src_affine = (399960.0, 5100060.0,  20.0, -20.0),
    dst_rows   = 1100,
    dst_cols   = 1100,
    dst_affine = (399960.0, 5100060.0,  10.0, -10.0),
)
```