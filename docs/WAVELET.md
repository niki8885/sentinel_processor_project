# Wavelet

## Overview

The wavelet module applies multi-level 2-D (and 3-D) Discrete Wavelet Transforms
to single-band rasters, multi-spectral stacks, and temporal image cubes.
All transforms are implemented in Fortran and exposed via `libsentinel_wavelet`:

```
input band  (rows × cols, float64)
    │
    ▼
dwt2d()            ──  forward N-level 2-D DWT → coefficient dict    [transform]
    │
    ▼
operation choice
    ├── threshold_coeffs()   ──  soft / hard wavelet thresholding      [denoising]
    ├── bayes_denoise()      ──  full BayesShrink pipeline             [denoising]
    ├── band_energy()        ──  Parseval energy of a sub-band         [analysis]
    └── band_stats()         ──  mean, variance, L1, L∞ of a sub-band [analysis]
    │
    ▼
idwt2d()           ──  inverse N-level 2-D DWT → reconstructed band  [transform]
```

For multi-spectral stacks and temporal cubes:

```
stack  (n_bands, rows, cols)         cube  (n_times, rows, cols)
    │                                    │
    ▼                                    ▼
dwt2d_batch()  ──  per-band DWT      dwt3d()  ──  separable 3-D DWT
    │                                    │         (2-D spatial + 1-D temporal)
    ▼                                    ▼
idwt2d_batch()                       idwt3d()
```

Noise estimation helpers (`estimate_sigma`, `bayes_threshold`) are also exported
for use outside the denoising pipeline.

## Install

```bat
# Windows
mkdir "...\sentinel_processor\wavelet\fortran"
gfortran -O2 -shared ^
  -o "...\sentinel_processor\wavelet\fortran\libsentinel_wavelet.dll" ^
     "...\sentinel_processor\wavelet\fortran\wavelet_mod.f90"
```

```bash
# Linux / macOS
gfortran -O2 -march=native -shared -fPIC \
  -o sentinel_processor/wavelet/fortran/libsentinel_wavelet.so \
     sentinel_processor/wavelet/fortran/wavelet_mod.f90
```

If the library is missing, every call raises `FileNotFoundError` with the
build command printed in the message.

## Module layout

```
sentinel_processor/
└── wavelet/
    ├── __init__.py           ← public re-exports
    ├── _wavelet_bridge.py    ← ctypes bridge (all public functions)
    └── fortran/
        ├── wavelet_mod.f90
        ├── libsentinel_wavelet.dll / .so
        └── (runtime DLLs on Windows — see Install)
```

## API reference

### dwt2d

```python
from sentinel_processor.wavelet import dwt2d

coeffs = dwt2d(arr, levels=3, wavelet="db4")
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `arr` | `np.ndarray` | — | 2-D `(rows, cols)` float64 array. `rows` and `cols` must each be divisible by `2^levels`. |
| `levels` | `int` | `1` | Number of decomposition levels (≥ 1). |
| `wavelet` | `str` | `"haar"` | Wavelet family; see [Supported wavelets](#supported-wavelets). |

Returns a nested dict `{ level: {'LL', 'LH', 'HL', 'HH'} }`.
`'LL'` (coarsest approximation) is present only at `level == levels`.
All sub-band arrays are C-order float64. Level `1` is the finest scale,
level `levels` is the coarsest.

### idwt2d

```python
from sentinel_processor.wavelet import idwt2d

reconstructed = idwt2d(coeffs, wavelet="db4")
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coeffs` | `dict` | — | Dict returned by `dwt2d`, optionally with thresholded sub-bands. |
| `wavelet` | `str` | `"haar"` | Must match the forward transform. |

Returns `(rows, cols)` float64 C-order array.
Perfect reconstruction error is at float64 machine epsilon (≤ 1e-9) for all
supported orthogonal wavelets.

### threshold_coeffs

```python
from sentinel_processor.wavelet import threshold_coeffs

thresholded = threshold_coeffs(coeffs, threshold=0.05, mode="soft")
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coeffs` | `dict` | — | Dict from `dwt2d`. |
| `threshold` | `float` | — | Non-negative threshold value. |
| `mode` | `"soft" \| "hard"` | `"soft"` | `"soft"` shrinks coefficients toward zero by `threshold`; `"hard"` zeroes all coefficients with `\|x\| ≤ threshold`. |

The `'LL'` approximation sub-band is never modified.
Returns a new dict with the same structure; input is not modified.

### estimate_sigma

```python
from sentinel_processor.wavelet import estimate_sigma

sigma_n = estimate_sigma(coeffs[1]["HH"])
```

| Parameter | Type | Description |
|---|---|---|
| `band` | `np.ndarray` | Any flat or 2-D sub-band array. Recommended: finest-scale `HH` sub-band. |

Returns `float`. Applies the Donoho-Johnstone MAD estimator:
`σ = median(|x|) / 0.6745`.

This estimator is robust to outliers and is the standard preprocessing step
before BayesShrink thresholding.

### bayes_threshold

```python
from sentinel_processor.wavelet import bayes_threshold

thr = bayes_threshold(coeffs[1]["LH"], sigma_n)
```

| Parameter | Type | Description |
|---|---|---|
| `band` | `np.ndarray` | 2-D or 1-D sub-band array. |
| `sigma_n` | `float` | Noise sigma from `estimate_sigma`. |

Returns `float`. Computes the per-sub-band BayesShrink optimal threshold:

```
σ_y² = E[x²]
σ_s  = sqrt(max(σ_y² − σ_n², 0))
T    = σ_n² / σ_s      (or max|x| when σ_s ≈ 0 → sub-band is pure noise)
```

Call separately per sub-band for adaptive per-scale thresholding, or use
the convenience wrapper `bayes_denoise` for the full pipeline.

### bayes_denoise

```python
from sentinel_processor.wavelet import bayes_denoise

denoised = bayes_denoise(arr, levels=3, wavelet="db4")
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `arr` | `np.ndarray` | — | 2-D `(rows, cols)` float64 array. |
| `levels` | `int` | `3` | Decomposition depth. |
| `wavelet` | `str` | `"db4"` | Wavelet family. |

Full pipeline in one call:

1. `dwt2d` → coefficient dict
2. `estimate_sigma` on the finest-scale `HH` sub-band
3. `bayes_threshold` per detail sub-band (LH, HL, HH at every level)
4. Soft-threshold all detail sub-bands
5. `idwt2d` → denoised image

Returns `(rows, cols)` float64 array.

### band_energy

```python
from sentinel_processor.wavelet import band_energy

e = band_energy(coeffs[1]["HH"])
```

| Parameter | Type | Description |
|---|---|---|
| `band` | `np.ndarray` | Any flat or 2-D sub-band array. |

Returns `float`. Computes the squared L2 norm (Parseval energy):
`Σ x²` over all elements.

Useful for scale-space energy analysis and compression ratio estimation.

### band_stats

```python
from sentinel_processor.wavelet import band_stats

stats = band_stats(coeffs[2]["LH"])
# {'mean': ..., 'var': ..., 'l1': ..., 'linf': ...}
```

| Parameter | Type | Description |
|---|---|---|
| `band` | `np.ndarray` | Any flat or 2-D sub-band array. |

Returns `dict` with keys:

| Key | Description |
|---|---|
| `mean` | Arithmetic mean |
| `var` | Population variance |
| `l1` | Mean absolute value (`Σ\|x\| / n`) |
| `linf` | Maximum absolute value |

### dwt2d_batch

```python
from sentinel_processor.wavelet import dwt2d_batch

coeffs_list = dwt2d_batch(stack, levels=2, wavelet="db4")
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `stack` | `np.ndarray` | — | `(n_bands, rows, cols)` float64 array. |
| `levels` | `int` | `1` | Decomposition depth. |
| `wavelet` | `str` | `"haar"` | Wavelet family. |

Applies `dwt2d` to each band independently.
Returns `list[dict]` of length `n_bands`, one coefficient dict per band (band-first order).

The Fortran call processes all bands in a single `dwt2d_batch` invocation —
no Python loop over bands.

### idwt2d_batch

```python
from sentinel_processor.wavelet import idwt2d_batch

reconstructed = idwt2d_batch(coeffs_list, wavelet="db4")
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coeffs_list` | `list[dict]` | — | List from `dwt2d_batch`. |
| `wavelet` | `str` | `"haar"` | Must match the forward transform. |

Returns `(n_bands, rows, cols)` float64 C-order array.

### dwt3d

```python
from sentinel_processor.wavelet import dwt3d

coeffs_vol = dwt3d(cube, levels=2, wavelet="haar")
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `arr` | `np.ndarray` | — | `(n_times, rows, cols)` float64 array. `rows`, `cols`, and `n_times` must each be divisible by `2^levels`. |
| `levels` | `int` | `1` | Decomposition depth (applied to both spatial and temporal axes). |
| `wavelet` | `str` | `"haar"` | Wavelet family. |

Separable 3-D DWT:

1. `dwt2d` applied to every time slice independently (spatial sub-bands)
2. 1-D DWT applied along the time axis for every spatial pixel (temporal sub-bands)

Returns `(n_times, rows, cols)` float64 coefficient array.
Temporal low sub-bands occupy slices `0 .. n_times//2 − 1`;
high sub-bands occupy slices `n_times//2 .. n_times − 1` (per level).

### idwt3d

```python
from sentinel_processor.wavelet import idwt3d

reconstructed = idwt3d(coeffs_vol, wavelet="haar")
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coeffs` | `np.ndarray` | — | `(n_times, rows, cols)` coefficient array from `dwt3d`. |
| `wavelet` | `str` | `"haar"` | Must match the forward transform. |

Returns `(n_times, rows, cols)` float64 C-order array.

## Supported wavelets

| Key | Name | Filter length | Notes |
|---|---|---|---|
| `"haar"` | Haar / Daubechies-1 | 2 | Fastest; piecewise-constant; strong blocking at high compression |
| `"db4"` | Daubechies-4 | 4 | Good default; 2 vanishing moments; compact support |
| `"db6"` | Daubechies-6 | 6 | Smoother than DB4; 3 vanishing moments |
| `"coif1"` | Coiflet-1 | 6 | Vanishing moments on both lo and hi filters; near-symmetric |
| `"sym4"` | Symlet-4 | 8 | Near-symmetric version of DB4; less phase distortion |
| `"sym6"` | Symlet-6 | 12 | Smoothest in the set; best frequency selectivity; needs larger arrays |

All wavelets are orthogonal. Synthesis = transpose of analysis — perfect
reconstruction to float64 machine epsilon (≤ 1e-9) is guaranteed.

**Array-size constraints by wavelet:**

| Wavelet | Min spatial dimension | Min `n_times` (3-D) |
|---|---|---|
| `haar`, `db4` | `2^levels` | `2^levels` |
| `db6`, `coif1` | `max(6, 2^levels)` | `max(6, 2^levels)` |
| `sym4` | `max(8, 2^levels)` | `max(8, 2^levels)` |
| `sym6` | `max(12, 2^levels)` | `max(12, 2^levels)` |

Practical safe minimum for 3 levels with `sym6`: **rows ≥ 12, cols ≥ 12**.

## Coefficient dict layout

`dwt2d` returns `{ level: sub_bands }` where:

```
level = 1      ← finest detail (full spatial resolution / 2)
level = 2
...
level = L      ← coarsest scale
```

Each `sub_bands` dict contains:

| Key | Spatial position | Description |
|---|---|---|
| `'LL'` | top-left quadrant | Coarse approximation (only at `level == L`) |
| `'LH'` | top-right quadrant | Horizontal detail (vertical edges) |
| `'HL'` | bottom-left quadrant | Vertical detail (horizontal edges) |
| `'HH'` | bottom-right quadrant | Diagonal detail |

Shape of each sub-band array at level `l`:
`(rows / 2^l,  cols / 2^l)`, C-order float64.

The flat in-memory layout mirrors the standard quad-tree tiling
used by `wavelet_mod.f90` — the Python conversion is zero-copy when
the array is already F-contiguous.

## Array layout

All arrays passed to Fortran are converted to **column-major (Fortran) order**
before the ctypes call. Input arrives as a C-order numpy array and is
transposed by `_to_f64_f`. Output buffers are allocated in Fortran order and
converted back to C-contiguous before being returned.

3-D stacks `(n_bands, rows, cols)` are transposed to `(rows·cols, n_bands)`
F-order so that Fortran reads spatial pixels as the fast-varying index.
Temporal cubes `(n_times, rows, cols)` are laid out as
`(rows·cols, n_times)` F-order; time slice `t` starts at element
`(t−1) × rows × cols`.

## Examples

### Denoise a single band

```python
import numpy as np
from sentinel_processor.wavelet import dwt2d, idwt2d, threshold_coeffs

band = ...  # np.ndarray shape (rows, cols)

# 3-level soft threshold with DB4
coeffs     = dwt2d(band, levels=3, wavelet="db4")
thresholded = threshold_coeffs(coeffs, threshold=0.02, mode="soft")
denoised   = idwt2d(thresholded, wavelet="db4")
```

### BayesShrink — one-liner

```python
from sentinel_processor.wavelet import bayes_denoise

denoised = bayes_denoise(band, levels=3, wavelet="db4")
```

### BayesShrink — manual per-sub-band control

```python
from sentinel_processor.wavelet import (
    dwt2d, idwt2d, estimate_sigma, bayes_threshold, threshold_coeffs
)

coeffs  = dwt2d(band, levels=3, wavelet="sym4")
sigma_n = estimate_sigma(coeffs[1]["HH"])

thresholded = {}
for lv, sub_bands in coeffs.items():
    thresholded[lv] = {}
    for key, sb in sub_bands.items():
        if key == "LL":
            thresholded[lv]["LL"] = sb.copy()
        else:
            thr = bayes_threshold(sb, sigma_n)
            thresholded[lv][key] = (
                np.sign(sb) * np.maximum(np.abs(sb) - thr, 0.0)
            )

denoised = idwt2d(thresholded, wavelet="sym4")
```

### Inspect sub-band energy across scales

```python
from sentinel_processor.wavelet import dwt2d, band_energy, band_stats

coeffs = dwt2d(band, levels=4, wavelet="db4")

for lv in sorted(coeffs):
    for key in ["LH", "HL", "HH"]:
        if key not in coeffs[lv]:
            continue
        e = band_energy(coeffs[lv][key])
        s = band_stats(coeffs[lv][key])
        print(f"L{lv} {key}  energy={e:.4f}  linf={s['linf']:.4f}")
```

### Denoise a multi-spectral stack

```python
from sentinel_processor.wavelet import dwt2d_batch, idwt2d_batch, band_energy

ms = ...  # np.ndarray shape (n_bands, rows, cols)

# Forward — one Fortran call for all bands
coeffs_list = dwt2d_batch(ms, levels=2, wavelet="db4")

# Threshold per band and per sub-band
import numpy as np
from sentinel_processor.wavelet import estimate_sigma, bayes_threshold

thresholded_list = []
for coeffs in coeffs_list:
    sigma_n = estimate_sigma(coeffs[1]["HH"])
    thr_coeffs = {}
    for lv, sub_bands in coeffs.items():
        thr_coeffs[lv] = {}
        for key, sb in sub_bands.items():
            if key == "LL":
                thr_coeffs[lv]["LL"] = sb.copy()
            else:
                thr = bayes_threshold(sb, sigma_n)
                thr_coeffs[lv][key] = (
                    np.sign(sb) * np.maximum(np.abs(sb) - thr, 0.0)
                )
    thresholded_list.append(thr_coeffs)

# Inverse — one Fortran call
ms_denoised = idwt2d_batch(thresholded_list, wavelet="db4")
# shape: (n_bands, rows, cols)
```

### 3-D DWT on a Sentinel time stack

```python
import numpy as np
from sentinel_processor.wavelet import dwt3d, idwt3d, threshold_coeffs

# cube: (n_times, rows, cols) — e.g. 8 monthly composites
cube = ...  # np.ndarray shape (8, 256, 256)

coeffs_vol = dwt3d(cube, levels=2, wavelet="haar")

# Zero out the temporal high-frequency sub-band (temporal smoothing)
n_times = cube.shape[0]
smoothed_vol = coeffs_vol.copy()
half = n_times // 2
smoothed_vol[half:] = 0.0   # zero temporal hi-band slices at coarsest level

reconstructed = idwt3d(smoothed_vol, wavelet="haar")
# result: temporally smoothed cube, same shape (8, 256, 256)
```

### Full pipeline with index computation

```python
import sentinel_processor as sp
from sentinel_processor.indices.compute import compute_indices
from sentinel_processor.wavelet import bayes_denoise
import xarray as xr, rioxarray, numpy as np

scene = "data/spectral/budapest_20260528T094727.nc"

# 1. Compute NDVI
idx_paths = compute_indices(scene, ["ndvi"])

# 2. Open NDVI raster
ndvi = rioxarray.open_rasterio(idx_paths["ndvi"]).squeeze().values.astype(np.float64)

# 3. Pad to power-of-2 multiples if necessary
from math import ceil
rows, cols = ndvi.shape
r2 = 2**3;  pr = ceil(rows / r2) * r2;  pc = ceil(cols / r2) * r2
padded = np.pad(ndvi, [(0, pr-rows), (0, pc-cols)], mode="reflect")

# 4. BayesShrink denoising
denoised = bayes_denoise(padded, levels=3, wavelet="sym4")[:rows, :cols]

# 5. Save
import rasterio, rasterio.transform as rt
with rasterio.open(idx_paths["ndvi"]) as src:
    meta = src.meta.copy()
with rasterio.open("data/filtered/ndvi_wavelet.tif", "w", **meta) as dst:
    dst.write(denoised.astype(np.float32), 1)
```

## Wavelet selection guide

| Goal | Recommended wavelet | Notes |
|---|---|---|
| Fast preview / prototype | `haar` | Cheapest; visible blocking at high thresholds |
| General denoising | `db4` | Best speed/quality balance; standard default |
| Smooth natural imagery | `sym4` or `sym6` | Less ringing than DB4; `sym6` needs arrays ≥ 12 px |
| Spectral indices, phenology | `coif1` | Vanishing moments on both filters; good for band ratios |
| Compression / coding | `db6` or `sym6` | Better frequency separation |
| Temporal cube analysis | `haar` | Exact Haar lifting; interpreting temporal sub-bands is simplest |

## Denoising mode guide

| Scenario | Recommended mode | Threshold |
|---|---|---|
| Sensor noise (Gaussian, known sigma) | `bayes_denoise` | Automatic (BayesShrink) |
| Salt-and-pepper / impulse noise | `threshold_coeffs` + `mode="hard"` | Tune manually |
| Compression artefact removal | `threshold_coeffs` + `mode="soft"` | 0.01–0.05 for normalised reflectance |
| Edge preservation critical | `sym4` or `sym6` + soft threshold | Lower threshold than for Haar |
| Aggressive smoothing | `haar` + `mode="hard"` | 0.05–0.2 |