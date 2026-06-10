# Texture

Compute GLCM-based texture features (energy, contrast, homogeneity) from
single Sentinel-2 bands. Texture features are widely used in land-cover
classification and urban mapping; the Fortran backend makes full-tile
computation practical where a naive Python/NumPy loop would be prohibitively
slow.

## Install

```bash
pip install sentinel-processor
```

## Import

```python
from sentinel_processor.texture import compute_glcm
```

Or explicitly:

```python
from sentinel_processor.texture.texture import compute_glcm
```

## Module layout

```
sentinel_processor/
└── texture/
    ├── __init__.py
    ├── texture.py               ← public API + NumPy fallback
    ├── _texture_bridge.py       ← ctypes bridge
    └── fortran/
        ├── texture_mod.f90
        └── libsentinel_texture.dll / .so
```

## Build

```bash
# Linux / macOS
gfortran -O2 -shared -fPIC \
    -o sentinel_processor/texture/fortran/libsentinel_texture.so \
    sentinel_processor/texture/fortran/texture_mod.f90

# Windows (MSYS2 / MinGW-w64)
gfortran -O2 -shared \
    -o sentinel_processor\texture\fortran\libsentinel_texture.dll \
    sentinel_processor\texture\fortran\texture_mod.f90
```

If the library is not compiled, `compute_glcm` falls back to a pure-NumPy
implementation and emits a `RuntimeWarning`. Output is identical; only
performance differs.

## API reference

### compute_glcm

```python
def compute_glcm(
    arr:      np.ndarray,   # (rows, cols)  float64-compatible
    window:   int = 7,      # neighbourhood side length (odd, >= 3)
    distance: int = 1,      # pixel offset for co-occurrence pair
    angle:    int = -1,     # 0 | 45 | 90 | 135 | -1 (isotropic)
) -> dict[str, np.ndarray]:
```

Computes three Haralick texture features from the Grey-Level Co-occurrence
Matrix (GLCM) of a 2-D raster band. Returns a `dict` with keys `"energy"`,
`"contrast"`, and `"homogeneity"`, each an `(rows, cols)` float64 array.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `arr` | `np.ndarray` | — | 2-D array `(rows, cols)`, any numeric dtype. Internally cast to `float64`. May contain `NODATA=-9999` values; those pixels are excluded from co-occurrence counts and produce `NODATA` in all output maps. |
| `window` | `int` | `7` | Side length of the local patch used to build each per-pixel GLCM (e.g. `7` → 7×7 neighbourhood). Must be odd and ≥ 3; even values are incremented by 1 automatically. Larger windows produce smoother feature maps but blur spatial edges. Recommended range: 5–15. |
| `distance` | `int` | `1` | Pixel offset for the co-occurrence pair. Use `1` for most classification workflows. Increase to `2`–`3` to capture coarser texture scales (e.g. urban block structure). |
| `angle` | `int` | `-1` | Direction of the co-occurrence offset. One of `{0, 45, 90, 135}` for a single direction, or `-1` for the isotropic average of all four. See angle table below. |

**Angle table**

| `angle` | Direction | Offset `(Δrow, Δcol)` |
|---|---|---|
| `0` | Horizontal → | `(0, +d)` |
| `45` | Diagonal ↗ | `(−d, +d)` |
| `90` | Vertical ↑ | `(−d, 0)` |
| `135` | Anti-diagonal ↖ | `(−d, −d)` |
| `−1` | Isotropic | Mean of all four |

**Returns**

`dict[str, np.ndarray]` with keys:

| Key | Feature | Range | Interpretation |
|---|---|---|---|
| `"energy"` | Angular Second Moment | [0, 1] | High = uniform / repetitive texture; low = heterogeneous surface |
| `"contrast"` | Local intensity variation | [0, (N−1)²] | High = strong grey-level differences between neighbours (edges, rough surfaces) |
| `"homogeneity"` | Inverse Difference Moment | [0, 1] | High = similar grey levels between neighbours (smooth / homogeneous regions) |

Border pixels within `window // 2 + distance` of any edge are set to
`NODATA = -9999.0`.

**Raises**

- `ValueError` — `arr` is not 2-D, or `angle` is not in `{−1, 0, 45, 90, 135}`.
- `FileNotFoundError` — shared library absent and no fallback is available (should not occur; the NumPy fallback always covers this case).

## Fortran internals

**Quantisation** — Grey levels are mapped to 64 bins once per tile before
the per-pixel loop. The GLCM matrix is always 64×64 regardless of the
original radiometric depth, keeping the inner loop cache-friendly.

**Symmetry** — Each co-occurrence pair `(i, j)` increments both `glcm(i,j)`
and `glcm(j,i)`, producing the symmetric normalised matrix used in the
standard Haralick definition.

**Isotropic mode** — The four directional GLCMs are averaged inside Fortran
before being written to the output arrays. No extra output storage is needed.

**Border filling** — Pixels whose neighbourhood or co-occurrence offset would
fall outside the raster are written as `NODATA` without entering the GLCM
loop, avoiding all out-of-bounds reads.

## Examples

### Minimal

```python
import numpy as np
from sentinel_processor.texture import compute_glcm

band = np.load("B04_10m.npy")                   # (rows, cols) float64
features = compute_glcm(band)

energy      = features["energy"]
contrast    = features["contrast"]
homogeneity = features["homogeneity"]
```

### Single direction

```python
# Horizontal co-occurrence only
features_h = compute_glcm(band, window=7, distance=1, angle=0)

# Vertical co-occurrence only
features_v = compute_glcm(band, window=7, distance=1, angle=90)
```

### Larger window and distance for coarse urban texture

```python
# 11×11 neighbourhood, 2-pixel offset — captures building block scale
features = compute_glcm(band, window=11, distance=2, angle=-1)
```

### Handle NODATA input

```python
import numpy as np
from sentinel_processor.texture import compute_glcm, NODATA

band = np.load("B04_10m.npy").astype(np.float64)
band[cloud_mask] = NODATA           # mark cloud pixels before passing in

features = compute_glcm(band, window=7, distance=1, angle=-1)

# Border and cloud pixels are NODATA in the output too
valid = features["energy"] != NODATA
print(f"Valid pixels: {valid.sum()} / {band.size}")
```

### All three features as a stacked array

```python
feats = compute_glcm(band, window=7)
texture_stack = np.stack(
    [feats["energy"], feats["contrast"], feats["homogeneity"]],
    axis=0,
)   # shape: (3, rows, cols)
```

### Multi-band stack

```python
import xarray as xr

band_names = ["energy", "contrast", "homogeneity"]
arrays = {}

for b_name in ["red", "nir", "swir16"]:
    band_data = scene.sel(band=b_name).values.astype(np.float64)
    feats     = compute_glcm(band_data, window=7, distance=1, angle=-1)
    for feat_name in band_names:
        arrays[f"{b_name}_{feat_name}"] = feats[feat_name]
```

### Integrate with raster_ops luminance

```python
from sentinel_processor.processing._raster_ops_bridge import rgb_to_luminance
from sentinel_processor.texture import compute_glcm

# Collapse RGB visual → single luminance band, then compute texture
rgb   = scene.sel(band=["red", "green", "blue"]).values  # (3, rows, cols)
lum   = rgb_to_luminance(rgb)                            # (rows, cols)
feats = compute_glcm(lum, window=9, distance=1, angle=-1)
```

### Integrate with timeseries stack

```python
import numpy as np
from sentinel_processor import stack_timeseries, TimeSeriesConfig
from sentinel_processor.texture import compute_glcm

result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg     = TimeSeriesConfig(max_cloud_fraction=0.15),
)

# Compute texture for every accepted time step on the NIR band
stack  = result.stack                           # (time, band, y, x)
nir_ts = stack.sel(band="nir").values           # (time, y, x)

energy_ts = np.stack(
    [compute_glcm(nir_ts[t].astype(np.float64))["energy"]
     for t in range(nir_ts.shape[0])],
    axis=0,
)   # (time, y, x)
```

### Window selection guide

```
Fine texture, 10 m resolution (urban, cropfield edges)  →  window=5 or 7
Medium texture, 10–20 m resolution                      →  window=7 or 9
Coarse texture (urban blocks, forest patches)           →  window=11–15,  distance=2–3
Smoothest possible feature map                          →  window=15,     angle=-1
```

## Notes

- **Input dtype** — Any numeric array is accepted; it is cast to `float64` internally. Passing `float32` arrays does not allocate a second copy if the bridge handles the cast, but the returned arrays are always `float64`.
- **Memory** — Each call allocates three `(rows, cols)` output arrays plus one Fortran-order copy of the input. For a 10 980×10 980 full Sentinel-2 tile (≈ 2.4 GB at `float64`) this is around 10 GB peak. Process sub-tiles or a clipped AOI if memory is constrained.
- **CRS / geotransform** — `compute_glcm` operates on raw pixel values only. CRS and affine transform are not consumed or written; wrap the outputs in `xarray.DataArray` or `rioxarray` manually if spatial reference is needed.
- **Fortran libraries** — if `libsentinel_texture.so/.dll` is not compiled, the module falls back to a pure-NumPy implementation automatically. The fallback is correct but significantly slower on full tiles. Build instructions: see the Build section above.