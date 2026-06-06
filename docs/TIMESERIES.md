# Time Series

Build and quality-filter temporal stacks of Sentinel-2 scenes.

## Install

```bash
pip install sentinel-processor
pip install "sentinel-processor[netcdf]"   # for .nc output
```

## Import

```python
from sentinel_processor import stack_timeseries, TimeSeriesConfig, StackResult
```

Or explicitly:

```python
from sentinel_processor.input.timeseries import stack_timeseries, TimeSeriesConfig, StackResult
```

## Module layout

```
sentinel_processor/
└── input/
    └── timeseries.py
```

Internally reuses — no logic duplicated:

| Dependency | Used for |
|---|---|
| `validation._fortran_bridge` | `call_check_dimensions`, `call_validate_scl`, `call_check_radiometry` |
| `processing._raster_ops_bridge` | `reproject_nearest` — Fortran grid alignment |
| `indices.compute._BAND_ALIASES` | canonical band name resolution for `require_bands` |

## Output

`stack_timeseries` returns a `StackResult`:

```
StackResult
├── .stack          xr.DataArray  (time, band, y, x)  float32
└── .scenes         list[SceneInfo]  — all candidates, accepted + rejected
```

When `save_dir` is set in `TimeSeriesConfig`, files are also written to disk:

```
<save_dir>/
├── <name>_stack.nc          ← fmt="nc" (default) — single NetCDF-4 file
└── <name>_stack_<ts>.tif    ← fmt="tif" — one GeoTIFF per time step
```

## API reference

### TimeSeriesConfig

| Parameter | Type | Default | Description |
|---|---|---|---|
| `align` | `bool` | `True` | Reproject all scenes to the reference grid. Same CRS → Fortran `reproject_nearest`; different CRS → `rioxarray.reproject_match` fallback |
| `reference` | `str \| Path \| None` | `None` | Explicit reference scene for alignment. `None` → first source used |
| `max_cloud_fraction` | `float` | `0.30` | Max allowed cloud fraction (0–1). Mirrors `DownloadConfig.max_cloud_threshold` |
| `min_confidence` | `float` | `0.01` | Min SCL confidence score (0–1). Same scoring table as Fortran `validate_scl`. Mirrors `DownloadConfig.min_confidence` |
| `max_snow_fraction` | `float` | `1.0` | Max allowed snow/ice fraction (0–1). `1.0` = disabled |
| `require_bands` | `list[str] \| None` | `None` | Reject scenes missing any of these bands. Accepts any alias: `"nir"`, `"B08"`, `"swir16"`, … |
| `use_sidecar_report` | `bool` | `True` | Read `*_report.json` sidecars written by the downloader before re-reading SCL files |
| `run_radiometry_check` | `bool` | `True` | Apply Fortran `check_radiometry` (< 1 % pixels above 15 000 DN) |
| `nodata` | `float` | `nan` | Fill value written for rejected scenes when `fill_rejected=True` |
| `fill_rejected` | `bool` | `False` | Keep rejected scenes as nodata planes to preserve a contiguous time axis |
| `save_dir` | `str \| Path \| None` | `None` | If set, the finished stack is saved here automatically |
| `save_format` | `str` | `"nc"` | Output format: `"nc"` = single NetCDF-4 file; `"tif"` = one GeoTIFF per time step |
| `save_name` | `str \| None` | `None` | Base filename without extension. `None` → derived from first scene, e.g. `"budapest_stack"` |

### stack_timeseries

```python
def stack_timeseries(
    sources: Sequence[str | Path],
    scl_dir: str | Path | None = None,
    cfg: TimeSeriesConfig | None = None,
    progress: bool = True,
) -> StackResult:
```

Searches for `*_report.json` sidecars first (fast path — no SCL re-read).
Falls back to reading the SCL file and calling the three Fortran validation
routines in the same order as the downloader:
`check_dimensions` → `validate_scl` → `check_radiometry`.

Sources are sorted chronologically by the `YYYYMMDDTHHMMSS` timestamp
embedded in each filename regardless of the order they are passed in.
Falls back to file mtime when no timestamp is found in the filename.

Returns `StackResult`. Raises `ValueError` when sources is empty, no scene
passes the quality filter, or accepted scenes have inconsistent shapes after
alignment.

### StackResult

```python
@dataclass
class StackResult:
    stack:  xr.DataArray     # (time, band, y, x), float32
    scenes: list[SceneInfo]

    n_accepted: int          # property
    n_rejected: int          # property

    def summary(self) -> str: ...
    def to_json(self, path) -> None: ...
    def save(self, output_dir, name=None, fmt="nc") -> list[str]: ...
```

#### StackResult.save

```python
def save(
    self,
    output_dir: str | Path,
    name: str | None = None,
    fmt: str = "nc",
) -> list[str]:
```

Saves the stack to `output_dir`. Returns a list of absolute paths written.

| `fmt` | Behaviour |
|---|---|
| `"nc"` | Single NetCDF-4 file `<name>.nc`. Time encoded as `int64` seconds since epoch — no precision loss |
| `"tif"` | One GeoTIFF per time step `<name>_<YYYYMMDDTHHMMSS>.tif`. Each plane saved with full CRS via rioxarray |

### SceneInfo

```python
@dataclass
class SceneInfo:
    path:             Path
    timestamp:        datetime
    accepted:         bool
    cloud_fraction:   float   # 0–1, water-excluded
    snow_fraction:    float   # 0–1
    confidence_score: float   # 0 | 0.5 | 0.75 | 1.0
    dimension_pass:   bool
    radiometry_pass:  bool
    reject_reason:    str     # "" when accepted
    report_source:    str     # "sidecar_json" | "scl_computed" | "none"
```

## Quality pipeline

Each scene goes through the same pipeline as the downloader, in the same order:

```
scene file
    │
    ▼
Priority 1: *_report.json sidecar   ← read if use_sidecar_report=True
    │                                  (free — already computed at download)
    │   not found ↓
    ▼
Priority 2: SCL file  →  Fortran validation
    │
    ├── check_dimensions()     min 32 px, max 10 980 px, aspect ≤ 4:1
    ├── validate_scl()         cloud/snow ratios → confidence score
    └── check_radiometry()     < 1 % pixels above 15 000 DN
    │
    ▼
apply TimeSeriesConfig thresholds
    ├── dimension_pass=False       → rejected
    ├── radiometry_pass=False      → rejected
    ├── cloud_fraction > max       → rejected
    ├── snow_fraction > max        → rejected
    └── confidence_score < min     → rejected
    │
    └── pass  →  open scene, check require_bands, align, add to stack
```

**Confidence score table** (identical to `validate_scl` Fortran routine):

| Condition | Score |
|---|---|
| No valid pixels after water exclusion | 0.0 |
| snow_ratio > 0.50 | 0.0 |
| cloud_ratio > max_cloud_fraction | 0.0 |
| cloud_ratio < 0.10 | 1.0 |
| cloud_ratio < 0.30 | 0.75 |
| cloud_ratio < 0.40 | 0.50 |
| otherwise | 0.0 |

## Examples

### Minimal

```python
from pathlib import Path
from sentinel_processor import stack_timeseries

result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
)
print(result.summary())
da = result.stack   # xr.DataArray (time, band, y, x)
```

### Strict quality filter

```python
from sentinel_processor import stack_timeseries, TimeSeriesConfig

result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg = TimeSeriesConfig(
        max_cloud_fraction = 0.05,
        min_confidence     = 0.75,
        require_bands      = ["red", "nir", "swir16"],
    ),
)
```

### Save stack automatically

```python
from sentinel_processor import stack_timeseries, TimeSeriesConfig

result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg = TimeSeriesConfig(
        max_cloud_fraction = 0.10,
        save_dir           = "data/stacks",
        save_format        = "nc",
        save_name          = "budapest_may",
    ),
)
# writes: data/stacks/budapest_may.nc
```

### Save as per-scene GeoTIFFs

```python
cfg = TimeSeriesConfig(
    save_dir    = "data/stacks",
    save_format = "tif",
)
result = stack_timeseries(sources, scl_dir="data/technical", cfg=cfg)
# writes: data/stacks/budapest_stack_20260520T094746.tif  …
```

### Save explicitly after inspection

```python
result = stack_timeseries(sources, scl_dir="data/technical")
print(result.summary())

if result.n_rejected == 0:
    paths = result.save("data/stacks", fmt="nc")
    print("Saved →", paths[0])
```

### Keep rejected scenes as nodata planes

```python
cfg = TimeSeriesConfig(
    max_cloud_fraction = 0.10,
    fill_rejected      = True,   # keeps time axis contiguous
    nodata             = float("nan"),
)
result = stack_timeseries(sources, scl_dir="data/technical", cfg=cfg)
# result.stack has shape (n_total, band, y, x) including rejected scenes as NaN
```

### Per-scene metadata

```python
for scene in result.scenes:
    status = "✓" if scene.accepted else f"✗ {scene.reject_reason}"
    print(f"{scene.timestamp:%Y-%m-%d}  cloud={scene.cloud_fraction:.1%}  {status}")

# export to JSON
result.to_json("data/stacks/budapest_may_log.json")
```

### Integrate with downloader

```python
import sentinel_processor as sp
from sentinel_processor import stack_timeseries, TimeSeriesConfig

cfg_dl = sp.DownloadConfig(
    bands               = sp.SpectralBands.ALL,
    tech_bands          = sp.TechnicalLayers.SCL,
    keep_items          = 20,
    max_cloud_threshold = 0.15,
    min_confidence      = 0.75,
    save_report         = True,   # ← enables fast sidecar path in stack_timeseries
)
sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=cfg_dl,
)

# stack_timeseries reads the *_report.json sidecars — no SCL re-read needed
result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg = TimeSeriesConfig(
        max_cloud_fraction = 0.15,
        min_confidence     = 0.75,
        save_dir           = "data/stacks",
    ),
)
```

### Stack from .tif files

```python
result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.tif")),
    scl_dir = "data/technical",
)
```

### No alignment (files already co-registered)

```python
cfg = TimeSeriesConfig(align=False)
result = stack_timeseries(sources, cfg=cfg)
```

### Recompute quality from SCL (ignore sidecars)

```python
cfg = TimeSeriesConfig(
    use_sidecar_report   = False,
    run_radiometry_check = True,
    max_cloud_fraction   = 0.10,
)
result = stack_timeseries(sources, scl_dir="data/technical", cfg=cfg)
```

## Notes

- **CRS** — CRS is carried from the first accepted scene to the output stack. If no CRS is present, the stack is written without spatial reference; use `rioxarray.write_crs()` to assign one manually.
- **Memory** — the entire stack is assembled in RAM. For very large stacks or many scenes, process a shorter time window or subset spatially with `rioxarray.clip_box` before stacking.
- **Band names** — `require_bands` accepts any alias from the indices module: `"nir"`, `"B08"`, `"nir_broad"` all resolve to the same band. The alias table is the single source of truth shared with `compute_indices`.
- **Fortran libraries** — if `libsentinel_validation.so/.dll` or `libsentinel_raster_ops.so/.dll` are not compiled, the module falls back to Python/NumPy equivalents automatically. Build instructions: see [VALIDATION.md](VALIDATION.md) and [RASTER_OPS.md](RASTER_OPS.md).

---

# Gap Filling

Fill cloud-masked holes in a (n_times, rows, cols) time stack using Fortran-accelerated routines.

## Module layout

```
sentinel_processor/
└── processing/
    ├── _timeseries_bridge.py    ← Python bridge (public API)
    └── fortran/
        ├── timeseries_mod.f90
        └── libsentinel_timeseries.dll / .so
```

## Build

```bash
# Linux / macOS
gfortran -O2 -shared -fPIC \
    -o sentinel_processor/processing/fortran/libsentinel_timeseries.so \
    sentinel_processor/processing/fortran/timeseries_mod.f90

# Windows (MSYS2 / MinGW)
gfortran -O2 -shared \
    -o sentinel_processor/processing/fortran/libsentinel_timeseries.dll \
    sentinel_processor/processing/fortran/timeseries_mod.f90
```

## Import

```python
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps
```

## API reference

### interpolate_gaps

```python
def interpolate_gaps(
    arr:    np.ndarray,          # (n_times, rows, cols)  float64
    mask:   np.ndarray,          # (n_times, rows, cols)  int or bool
    method: str  = "linear",     # see method table below
    window: int  = 5,            # S-G / Gaussian window width
) -> np.ndarray:                 # (n_times, rows, cols)  float64
```

Returns a **copy** of `arr` with gap positions filled. The original array is never mutated.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `arr` | `np.ndarray` | — | 3-D array `(n_times, rows, cols)`, any numeric dtype. Internally cast to `float64`. |
| `mask` | `np.ndarray` | — | Same shape as `arr`. Non-zero = valid observation, `0` = gap to fill. |
| `method` | `str` | `"linear"` | Gap-filling algorithm. See method table below. |
| `window` | `int` | `5` | Window width for `savgol` and `gauss`. Silently forced odd and ≥ 3. Ignored for `linear`, `pchip`, `ets`. |

**Method table**

| `method` | Fortran flag | Algorithm | Best for |
|---|---|---|---|
| `"linear"` | 0 | Piecewise linear between nearest valid neighbours | Fast baseline, short gaps |
| `"savgol"` | 1 | Savitzky-Golay quadratic smoother | Noisy series; preserves peaks and troughs |
| `"pchip"` | 2 | Monotone cubic Hermite spline (Fritsch-Carlson) | NDVI / EVI / LAI — no overshoot, physically plausible |
| `"ets"` | 3 | Holt double-exponential smoothing (auto α, β) | Series with persistent seasonal trend |
| `"gauss"` | 4 | Gaussian-weighted moving average | Gentle smoothing; symmetric, renormalised at edges |

**Edge cases (all methods)**

| Situation | Result |
|---|---|
| Pixel has no valid observations | All time steps set to `-9999.0` (nodata) |
| Gap at the start of the series | Filled with the first valid value (nearest-neighbour) |
| Gap at the end of the series | Filled with the last valid value (nearest-neighbour) |
| Single valid observation | Constant fill across the entire time axis |

**Raises**

- `ValueError` — `arr` is not 3-D, `mask` shape does not match `arr`, or `method` is not recognised.
- `FileNotFoundError` — shared library not compiled (see Build above).

## Method notes

**`linear`** — The safest default. Splits each gap into its bounding valid observations and draws a straight line between them. No parameters, no smoothing, exact at the boundary points.

**`savgol`** — Fits a local quadratic polynomial over a sliding `window` of time steps (Savitzky-Golay). Gaps are pre-filled linearly before smoothing, so valid observations are always the anchor; only gap positions receive the smoothed estimate. Prefer `window=5` for 10-day composites, `window=7..11` for daily data with many cloudy days.

**`pchip`** — Piecewise Cubic Hermite Interpolating Polynomial with Fritsch-Carlson monotone tangents. Unlike a natural cubic spline it cannot overshoot between two observations — if NDVI drops from 0.8 to 0.3, no interpolated value exceeds 0.8 or goes below 0.3. Use this when physical value bounds must be respected.

**`ets`** — Holt's double-exponential (level + linear trend). The smoothing parameters α and β are chosen automatically by a 5×5 grid search that minimises MSE on the valid observations. Because this is a forward-running filter rather than a true interpolant, it works best for long gaps in series with a clear upward or downward trend (post-fire recovery, spring green-up). For short symmetric gaps, `pchip` is more accurate.

**`gauss`** — Convolution with a Gaussian kernel of σ = `window / 4` time steps. The kernel is renormalised at the series edges so there is no amplitude loss near the boundaries. Gaps are pre-filled linearly before the convolution. Produces softer transitions than S-G with the same window.

## Examples

### Basic usage — one method

```python
import numpy as np
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps

# arr: (n_times, rows, cols) float64 — your time stack
# valid: True where pixel is not cloud-masked
filled = interpolate_gaps(arr, valid_mask, method="pchip")
```

### Build mask from a stack with NaN nodata

```python
import numpy as np

# stack from stack_timeseries() with fill_rejected=True, nodata=float("nan")
arr  = result.stack.values.astype(np.float64)   # (time, band, y, x)
mask = np.isfinite(arr).astype(np.int32)

# fill each band independently
filled_bands = np.stack([
    interpolate_gaps(arr[:, b, :, :], mask[:, b, :, :], method="pchip")
    for b in range(arr.shape[1])
], axis=1)
```

### Build mask from a separate cloud mask array

```python
# cloud_mask: (n_times, rows, cols) — True where cloudy
valid_mask = (~cloud_mask).astype(np.int32)

filled = interpolate_gaps(
    arr.astype(np.float64),
    valid_mask,
    method="savgol",
    window=7,
)
```

### Compare methods on a single band

```python
import numpy as np
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps

ndvi = arr[:, nir_idx, :, :].astype(np.float64)   # (n_times, rows, cols)
mask = valid_mask[:, nir_idx, :, :]

results = {
    m: interpolate_gaps(ndvi, mask, method=m)
    for m in ("linear", "savgol", "pchip", "ets", "gauss")
}
# results["pchip"] has no overshoot; results["ets"] tracks the trend best
```

### Full end-to-end pipeline

```python
import numpy as np
from pathlib import Path
from sentinel_processor import stack_timeseries, TimeSeriesConfig
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps

# 1. Build stack — keep rejected scenes as nodata planes
cfg = TimeSeriesConfig(
    max_cloud_fraction = 0.30,
    fill_rejected      = True,
    nodata             = float("nan"),
)
result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg     = cfg,
)
print(result.summary())

# 2. Prepare inputs — stack is (time, band, y, x)
arr  = result.stack.values.astype(np.float64)
mask = np.isfinite(arr).astype(np.int32)

# 3. Fill each band
n_bands = arr.shape[1]
filled  = np.empty_like(arr)
for b in range(n_bands):
    filled[:, b, :, :] = interpolate_gaps(
        arr [:, b, :, :],
        mask[:, b, :, :],
        method = "pchip",
    )

# 4. Wrap back into an xarray DataArray (preserves coords and CRS)
import xarray as xr
filled_da = xr.DataArray(
    filled.astype(np.float32),
    dims   = result.stack.dims,
    coords = result.stack.coords,
    attrs  = result.stack.attrs,
)
```

### Savitzky-Golay — choose window by temporal resolution

```python
# Sentinel-2 revisit ~5 days → 10-day composite → window=5 (25 days)
filled_10d = interpolate_gaps(ndvi, mask, method="savgol", window=5)

# Daily data with high cloud frequency → wider window
filled_1d  = interpolate_gaps(ndvi, mask, method="savgol", window=11)
```

### Gaussian — wider window for smoother phenology curves

```python
# sigma = window/4 = 15/4 = 3.75 time steps
filled_smooth = interpolate_gaps(ndvi, mask, method="gauss", window=15)
```

### Holt ETS — post-fire recovery trend

```python
# Series with strong upward trend after disturbance
filled_recovery = interpolate_gaps(
    nbr_postfire.astype(np.float64),
    valid_mask,
    method = "ets",
)
```

### Handle nodata pixels explicitly

```python
NODATA = -9999.0

filled = interpolate_gaps(arr, mask, method="pchip")

# Pixels with no valid observations at all are marked -9999.0
no_data_pixels = filled == NODATA
print(f"Pixels with zero valid obs: {no_data_pixels.any(axis=0).sum()}")
```

## Method selection guide

```
Short gaps (1–3 steps), any series     →  linear
Noisy series, need smoothing           →  savgol   (window 5–11)
NDVI / EVI / LAI, physical bounds      →  pchip
Recovery / phenology with clear trend  →  ets
Need very smooth phenology curve       →  gauss    (window 9–15)
```