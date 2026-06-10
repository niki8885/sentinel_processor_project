# Band Covariance

Fortran-accelerated per-band covariance matrix for `(n_bands, rows, cols)` raster cubes.
Mathematical foundation for PCA-based sharpening, spectral feature reduction, and
multi-band anomaly detection.

## Install

```bash
pip install sentinel-processor
```

Build the Fortran library once (required before `band_covariance` can be called):

```bash
# Linux / macOS
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/analysis/fortran/libband_covariance.so \
  sentinel_processor/analysis/fortran/band_covariance.f90

# Windows (MSYS2 / MinGW ucrt64)
gfortran -O2 -shared \
  -static-libgfortran -static-libgcc \
  -o sentinel_processor\analysis\fortran\libband_covariance.dll \
  sentinel_processor\analysis\fortran\band_covariance.f90
```

Or use the build script:

```bash
python build_band_covariance.py          # Release
python build_band_covariance.py --debug  # -g -O0 -fcheck=all
python build_band_covariance.py --check  # check dependencies only, no compile
```

## Import

```python
from sentinel_processor.analysis._band_covariance_bridge import (
    band_covariance,
    NODATA,
)
```

## Module layout

```
sentinel_processor/
└── analysis/
    ├── __init__.py
    ├── _band_covariance_bridge.py
    └── fortran/
        ├── band_covariance.f90
        ├── libband_covariance.dll   ← Windows
        └── libband_covariance.so    ← Linux / macOS
```

## Array layout

All functions share the same conventions:

| Convention | Detail |
|---|---|
| Input cube | `(n_bands, rows, cols)` `float64`, C-contiguous |
| Output matrix | `(n_bands, n_bands)` `float64` |
| Missing value | `NODATA = -9999.0` — used in both input and output |

Arrays of any numeric dtype are accepted and cast to `float64` internally.
The input array is never mutated — all functions return new buffers.

## NODATA handling

`NODATA` values (`-9999.0`) are excluded from all statistics.
For an off-diagonal entry `[i, j]` a pixel contributes only when **both** band `i`
and band `j` carry a non-NODATA value at that spatial position.
Output entries are `NODATA` where the result is undefined (fewer than 2 jointly
valid observations). Check for `NODATA` explicitly — do not rely on `np.nan`:

```python
valid = result != NODATA
mean_variance = result[np.eye(n_bands, dtype=bool) & valid].mean()
```

---

## API reference

### band_covariance

```python
def band_covariance(arr: np.ndarray) -> np.ndarray:
```

Per-band sample covariance matrix over all spatial pixels.

Computes the `(n_bands × n_bands)` covariance matrix of band reflectances across
the full spatial extent of the image. A single-pass Fortran kernel avoids loading
the full array multiple times, making this suitable for large tiles.

Kahan compensated summation is used in both the mean and cross-product accumulation
steps to minimise floating-point cancellation errors.

**Parameters**

| Parameter | Type | Description |
|---|---|---|
| `arr` | `np.ndarray` `(n_bands, rows, cols)` | Multi-band raster cube |

**Returns** `np.ndarray`, shape `(n_bands, n_bands)`, dtype `float64` —
symmetric sample covariance matrix (divided by N − 1).
Entry `[i, j]` is `NODATA` when fewer than 2 pixels are jointly valid in bands `i`
and `j`.

**Algorithm**

| Pass | Operation |
|---|---|
| Pass 1 | Per-band pixel sums with Kahan compensation → band means |
| Pass 2 | Centred cross-products `(val_i − mean_i)(val_j − mean_j)` with Kahan compensation, upper triangle only |
| Finalise | Divide by `joint_count − 1` → sample covariance; mirror to lower triangle |

**Raises** `ValueError` when `arr` is not exactly 3-D.
**Raises** `FileNotFoundError` if the compiled Fortran shared library is not present.

---

## Fortran subroutine

The C-exported subroutine lives in `band_covariance.f90` in module `band_covariance_mod`.

| Subroutine | Python function | Signature |
|---|---|---|
| `band_covariance` | `band_covariance` | `(arr, rows, cols, n_bands, cov_out)` |

The internal flat C/row-major array layout:

```
arr[b, r, c]  →  flat index  =  (b-1) * rows * cols  +  (r-1) * cols  +  c   (1-based in Fortran)
```

Output `cov_out` is row-major:

```
cov_out[i, j]  →  flat index  =  (i-1) * n_bands  +  j
```

---

## Examples

### PCA via eigendecomposition

```python
import numpy as np
from sentinel_processor.analysis._band_covariance_bridge import band_covariance

# stack: (n_bands, rows, cols) float64
cov = band_covariance(stack)                        # (n_bands, n_bands)
eigenvalues, eigenvectors = np.linalg.eigh(cov)

# eigenvectors[:, -1] is the first principal component direction
pc1_direction = eigenvectors[:, -1]

# project the stack onto the first two PCs
flat = stack.reshape(stack.shape[0], -1)            # (n_bands, pixels)
scores = eigenvectors[:, -2:].T @ flat              # (2, pixels)
pc_image = scores.reshape(2, *stack.shape[1:])      # (2, rows, cols)
```

### Mahalanobis distance for anomaly detection

```python
import numpy as np
from sentinel_processor.analysis._band_covariance_bridge import band_covariance, NODATA

cov  = band_covariance(stack)                       # (n_bands, n_bands)

# replace NODATA sentinel with NaN before inversion
cov_clean = np.where(cov == NODATA, np.nan, cov)
prec = np.linalg.inv(cov_clean)                     # precision matrix

# Mahalanobis distance per pixel
flat  = stack.reshape(stack.shape[0], -1)           # (n_bands, pixels)
mu    = np.nanmean(flat, axis=1, keepdims=True)
delta = flat - mu                                   # (n_bands, pixels)
d2    = np.einsum("ij,jk,ki->i", delta.T, prec, delta.T.T)  # (pixels,)
mahal = np.sqrt(d2).reshape(stack.shape[1:])        # (rows, cols)

anomaly = mahal > 3.0                               # ~3-sigma threshold
print(f"Anomalous pixels: {anomaly.sum():,}")
```

### Spectral correlation matrix

```python
import numpy as np
from sentinel_processor.analysis._band_covariance_bridge import band_covariance, NODATA

cov = band_covariance(stack)                        # (n_bands, n_bands)

# convert covariance to correlation
valid = cov != NODATA
diag  = np.diag(cov)                               # per-band variances
std   = np.sqrt(diag)                              # per-band std deviations
corr  = cov / np.outer(std, std)                   # (n_bands, n_bands)
corr[~valid] = NODATA

# find highly correlated band pairs (r > 0.95) — candidates for redundancy removal
hi  = np.argwhere((corr > 0.95) & (corr != 1.0) & (corr != NODATA))
for i, j in hi:
    print(f"  Band {i} ↔ Band {j}  r = {corr[i, j]:.3f}")
```

### Full PCA sharpening pipeline

```python
import numpy as np
from sentinel_processor.analysis._band_covariance_bridge import band_covariance, NODATA

# stack: (n_bands, rows, cols) float64, NODATA already set
cov  = band_covariance(stack)                       # (n_bands, n_bands)

cov_clean = np.where(cov == NODATA, 0.0, cov)
eigenvalues, eigenvectors = np.linalg.eigh(cov_clean)

# sort descending by explained variance
order       = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[order]
eigenvectors = eigenvectors[:, order]

explained = eigenvalues / eigenvalues.sum()
n_components = int(np.searchsorted(np.cumsum(explained), 0.99)) + 1
print(f"Components for 99 % variance: {n_components} of {stack.shape[0]}")

# forward transform
flat    = stack.reshape(stack.shape[0], -1)         # (n_bands, pixels)
reduced = eigenvectors[:, :n_components].T @ flat   # (n_components, pixels)
```

---

## Method selection guide

```
Band-to-band relationships      →  band_covariance  (full covariance matrix)
PCA / feature reduction         →  band_covariance  →  np.linalg.eigh
Mahalanobis anomaly detection   →  band_covariance  →  np.linalg.inv  →  distance
Spectral redundancy check       →  band_covariance  →  correlation matrix
Temporal pixel correlation      →  pearson_map  (from _sentinel_stats_bridge)
```

---

## Numerical notes

Kahan compensated summation reduces the accumulated floating-point error of a naive
sum from O(n · ε) to O(ε), where ε is machine epsilon (~2.2 × 10⁻¹⁶ for `float64`).
This matters for large tiles (10 000 × 10 000 pixels) where naive summation of
reflectance values can accumulate errors of the same order as the variance itself.

The two-pass structure (means first, cross-products second) is the standard
numerically stable algorithm for covariance. The single-pass Youngs-Cramer
online algorithm is not used here because the Fortran kernel fits entirely in
L2 cache at typical tile sizes, making the two-pass approach both simpler and
faster in practice.