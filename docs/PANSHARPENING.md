# Pansharpening

## Overview

Pansharpening increases the spatial resolution of multispectral bands by injecting
detail from a higher-resolution panchromatic (PAN) source. Three algorithms are
implemented in Fortran and exposed via `libsentinel_processing`:

```
MS bands (coarser resolution)
    │
    ▼
bilinear_upsample()     ──  upsample each MS band to PAN grid  [pansharpening]
    │
    ▼
rgb_to_luminance()      ──  collapse visual RGB → single PAN   [raster_ops]
    │
    ▼
algorithm choice
    ├── gram_schmidt_sharpen()   ──  any number of bands
    ├── ihs_sharpen()            ──  exactly 3 bands
    └── wavelet_sharpen()        ──  any number of bands
    │
    ▼
sharpened MS stack at PAN resolution  →  saved to spectral/
```

Pansharpening is **disabled by default** (`pansharpen_algorithm = None`).
Enable it by setting `pansharpen_algorithm` in `DownloadConfig`.

## Install

```bat
# Windows
gfortran -O2 -shared -o "...\sentinel_processor\processing\fortran\libsentinel_processing.dll" "...\sentinel_processor\processing\fortran\pansharpening.f90"
```

```bash
# Linux / macOS
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/processing/fortran/libsentinel_processing.so \
  sentinel_processor/processing/fortran/pansharpening.f90
```

`libsentinel_raster_ops` is also required for the PAN preparation step
(`rgb_to_luminance`). See [RASTER_OPS.md](RASTER_OPS.md) for its build command.
If either library is missing, pansharpening is skipped with a warning.

## Module layout

```
sentinel_processor/
└── processing/
    ├── _fortran_bridge.py        ← pansharpening bridge
    ├── _raster_ops_bridge.py     ← raster ops bridge
    └── fortran/
        ├── pansharpening.f90
        ├── raster_ops.f90
        ├── libsentinel_processing.dll / .so
        └── libsentinel_raster_ops.dll / .so
```

## API reference

### DownloadConfig — pansharpening fields

| Parameter | Type | Default | Description |
|---|---|---|---|
| `pansharpen_algorithm` | `"gram_schmidt" \| "ihs" \| "wavelet" \| None` | `None` | Algorithm to use; `None` disables pansharpening |
| `pansharpen_pan_key` | `str` | `"visual"` | STAC asset key used as the PAN source |

### pansharpen

```python
from sentinel_processor.processing._fortran_bridge import pansharpen

sharpened = pansharpen(
    pan,        # np.ndarray, shape (rows, cols)
    ms,         # np.ndarray, shape (n_bands, ms_rows, ms_cols)
    algorithm,  # "gram_schmidt" | "ihs" | "wavelet"
)
# returns np.ndarray, shape (n_bands, rows, cols), dtype float64
```

### Convenience wrappers

```python
from sentinel_processor.processing._fortran_bridge import (
    pansharpening_gs,
    pansharpening_ihs,
    pansharpening_wavelet,
)
```

## Algorithms

### Gram-Schmidt (`gram_schmidt`)

Standard remote-sensing spectral sharpening. Works with any number of bands.

Steps:
1. Upsample each MS band to PAN resolution (bilinear, in Fortran)
2. Compute a synthetic low-resolution PAN as the mean of all upsampled bands
3. Orthogonalise the band stack via Gram-Schmidt, starting from the synthetic PAN
4. Histogram-match the real PAN to the synthetic PAN (preserves radiometric scale)
5. Replace the first GS component with the matched PAN
6. Back-project to recover sharpened bands

Good default choice. Preserves relative spectral relationships well.

### IHS (`ihs`)

Intensity-Hue-Saturation sharpening. **Requires exactly 3 bands.**

Steps:
1. Upsample MS bands to PAN resolution
2. Forward IHS transform: `I = (b1+b2+b3)/3`, `v1`, `v2` are orthogonal chrominance axes
3. Histogram-match PAN to the I channel
4. Inject matched PAN as the new I
5. Inverse IHS transform

Fastest option for RGB triplets.

### Wavelet (`wavelet`)

Haar lifting-scheme wavelet injection. Works with any number of bands.

Steps (per band):
1. Upsample MS band to PAN resolution
2. One-level 2D forward Haar DWT on PAN → LL, LH, HL, HH subbands
3. One-level 2D forward Haar DWT on upsampled MS band
4. Replace the LL subband of MS with the LL subband of PAN
5. Inverse DWT → sharpened band

Retains original high-frequency texture from the MS bands while injecting
PAN spatial structure into the low-frequency subband.

## PAN source on Sentinel-2

Sentinel-2 does not have a dedicated panchromatic band. The `visual` STAC asset
(the default `pansharpen_pan_key`) is an RGB composite at 10 m. The bridge
collapses it to a single luminance channel via `rgb_to_luminance` (Fortran):

| Channels in asset | Collapse method |
|---|---|
| 1 | Used directly |
| 3 (RGB) | Rec. 601: `Y = 0.299·R + 0.587·G + 0.114·B` |
| > 3 | Equal-weight mean: `1 / n_bands` per band |

See [RASTER_OPS.md](RASTER_OPS.md#rgb_to_luminance) for implementation details.

## Array layout

All arrays passed to Fortran are converted to **column-major (Fortran) order**
before the ctypes call. Multispectral arrays arrive in Python as
`(n_bands, rows, cols)` and are transposed to `(rows, cols, n_bands)` + `order='F'`
so that Fortran reads the spatial dimensions as the fast-varying indices.
The output buffer is allocated in the same order and transposed back to
`(n_bands, rows, cols)` C-contiguous before being returned.

## Examples

### Enable during download

```python
import sentinel_processor as sp

# Gram-Schmidt — any number of bands
cfg = sp.DownloadConfig(
    bands                = sp.SpectralBands.RGB_NIR,
    pansharpen_algorithm = "gram_schmidt",
)

# IHS — exactly 3 bands required
cfg = sp.DownloadConfig(
    bands                = sp.SpectralBands.RGB,
    pansharpen_algorithm = "ihs",
)

# Wavelet — any number of bands
cfg = sp.DownloadConfig(
    bands                = sp.SpectralBands.ALL,
    pansharpen_algorithm = "wavelet",
)

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=cfg,
)
```

### Low-level call

```python
import numpy as np
from sentinel_processor.processing._fortran_bridge import pansharpen

pan = np.random.rand(1000, 1000)       # high-res PAN
ms  = np.random.rand(10, 250, 250)     # 10 MS bands at coarser resolution

sharpened = pansharpen(pan, ms, algorithm="gram_schmidt")
print(sharpened.shape)   # (10, 1000, 1000)
```

## Output

Pansharpened bands are saved in place of the original spectral stack — same
filenames, same output directory (`spectral/`), but at PAN resolution.
The `_report.json` sidecar (when `save_report=True`) gains extra fields:

```json
{
  "item_id": "S2B_34TCT_20260501_1_L2A",
  "passed": true,
  "pansharpen_algorithm": "gram_schmidt",
  "pan_key": "visual",
  ...
}
```