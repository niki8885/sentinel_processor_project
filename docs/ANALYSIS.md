# Analysis

Fortran-accelerated per-pixel temporal statistics for `(n_times, rows, cols)` raster stacks.

## Install

```bash
pip install sentinel-processor
```

Build the Fortran library once (required before any function can be called):

```bash
# Linux / macOS
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/analysis/fortran/libsentinel_stats.so \
  sentinel_processor/analysis/fortran/sentinel_stats.f90

# Windows (MSYS2 / MinGW ucrt64)
gfortran -O2 -shared \
  -static-libgfortran -static-libgcc \
  -o sentinel_processor\analysis\fortran\libsentinel_stats.dll \
  sentinel_processor\analysis\fortran\sentinel_stats.f90
```

Or use the build script:

```bash
python build_sentinel_stats.py          # Release
python build_sentinel_stats.py --debug  # -g -O0 -fcheck=all
python build_sentinel_stats.py --check  # check dependencies only, no compile
```

## Import

```python
from sentinel_processor.analysis._sentinel_stats_bridge import (
    # quality
    valid_obs_count,
    temporal_gap_stats,
    # distribution
    pixel_quantiles,
    pixel_iqr,
    # core statistics
    time_window_stats,
    anomaly_zscore,
    # trend
    trend_theil_sen,
    mann_kendall,
    # change detection
    bfast_breakpoint,
    # regression
    pixel_regression,
    # phenology
    phenology_doy,
    phenology_metrics,
    save_phenology,
    # correlation
    pearson_map,
    # sentinel value
    NODATA,
)
```

## Module layout

```
sentinel_processor/
└── analysis/
    ├── __init__.py
    ├── _sentinel_stats_bridge.py
    └── fortran/
        ├── sentinel_stats.f90
        ├── libsentinel_stats.dll   ← Windows
        └── libsentinel_stats.so    ← Linux / macOS
```

## Array layout

All functions share the same conventions:

| Convention | Detail |
|---|---|
| Input stack | `(n_times, rows, cols)` `float64`, C-contiguous |
| `dates` | `list[datetime]` — length must equal `arr.shape[0]`; timezone-aware datetimes accepted |
| 2-D outputs | `(rows, cols)` `float64` |
| 3-D outputs | `(n_times, rows, cols)` `float64` |
| Missing value | `NODATA = -9999.0` — used in both input and output |

Arrays of any numeric dtype are accepted and cast to `float64` internally.
The input array is never mutated — all functions return new buffers.

## NODATA handling

All functions skip `NODATA` values (`-9999.0`) in computation and write `NODATA`
to output pixels where the result is undefined (too few valid observations,
zero variance, etc.). Check for `NODATA` explicitly — do not rely on `np.nan`:

```python
valid = result["slope"] != NODATA
mean_slope = result["slope"][valid].mean()
```

---

## API reference

### valid_obs_count

```python
def valid_obs_count(arr: np.ndarray) -> ObsCountResult:
```

Count and fraction of non-NODATA observations per pixel.

**Returns** `ObsCountResult`:

| Key | Shape | Description |
|---|---|---|
| `count` | `(rows, cols)` | Number of valid scenes, cast to `float64` |
| `fraction` | `(rows, cols)` | `count / n_times` in `[0, 1]` |

Use `fraction < 0.3` to mask poorly-covered pixels before analysis.

---

### temporal_gap_stats

```python
def temporal_gap_stats(
    arr:   np.ndarray,
    dates: list[datetime],
) -> GapStatsResult:
```

Maximum and mean gap in days between consecutive valid observations.

**Returns** `GapStatsResult`:

| Key | Shape | Description |
|---|---|---|
| `max_gap` | `(rows, cols)` | Longest gap between consecutive valid obs |
| `mean_gap` | `(rows, cols)` | Mean gap between consecutive valid obs |

`NODATA` for pixels with fewer than 2 valid observations.
Use `max_gap > 60` to flag pixels with dangerous data holes before gap-filling.

---

### pixel_quantiles

```python
def pixel_quantiles(arr: np.ndarray) -> QuantileResult:
```

Empirical quantiles of the time series per pixel (linear interpolation).

**Returns** `QuantileResult` — keys `p10`, `p25`, `p50`, `p75`, `p90`, each `(rows, cols)`.

`NODATA` for pixels with no valid observations.

---

### pixel_iqr

```python
def pixel_iqr(arr: np.ndarray) -> IQRResult:
```

Interquartile range and Tukey outlier mask (`k = 1.5`) per pixel.

**Returns** `IQRResult`:

| Key | Shape | Description |
|---|---|---|
| `iqr` | `(rows, cols)` | Q75 − Q25 |
| `outlier` | `(n_times, rows, cols)` | `1.0` = outlier · `0.0` = within fence · `NODATA` = missing |

Feed `1 - outlier` as a validity mask into `interpolate_gaps` to replace cloud
residuals that bypassed SCL filtering.

---

### time_window_stats

```python
def time_window_stats(
    arr:         np.ndarray,
    dates:       list[datetime],
    window_days: int = 30,
) -> WindowStatsResult:
```

Rolling-window mean, sample std, and OLS linear slope per pixel.

The window is centred on each pixel's own valid-date mid-point; all observations
within `± window_days / 2` days of that centre are included. The window is
therefore adaptive per pixel rather than anchored to a single calendar date.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `arr` | `np.ndarray` | — | `(n_times, rows, cols)` stack |
| `dates` | `list[datetime]` | — | Acquisition timestamps |
| `window_days` | `int` | `30` | Total window width in days (must be ≥ 1) |

**Returns** `WindowStatsResult`:

| Key | Shape | Description |
|---|---|---|
| `mean` | `(rows, cols)` | Mean of valid observations inside the window |
| `std` | `(rows, cols)` | Sample standard deviation (N-1); `0.0` when only one valid obs |
| `slope` | `(rows, cols)` | OLS linear slope in `[units / day]` |

**Raises** `ValueError` when `window_days < 1` or `len(dates) != arr.shape[0]`.

---

### anomaly_zscore

```python
def anomaly_zscore(
    arr:         np.ndarray,
    dates:       list[datetime],
    window_days: int = 30,
    *,
    background:  WindowStatsResult | None = None,
) -> np.ndarray:
```

Per-pixel, per-scene z-score relative to the temporal background:
`z(t, r, c) = (arr[t, r, c] − mean[r, c]) / std[r, c]`.

When `background` is `None` it is computed internally with `time_window_stats`
using `window_days`. Pass a pre-computed `WindowStatsResult` to reuse the same
background across multiple stacks (e.g. comparing two years).

**Returns** `(n_times, rows, cols)` `float64`.
`NODATA` where `std == 0`, the input value is `NODATA`, or `mean` / `std` is `NODATA`.

---

### trend_theil_sen

```python
def trend_theil_sen(
    arr:   np.ndarray,
    dates: list[datetime],
) -> TheilSenResult:
```

Robust Theil-Sen median slope and intercept per pixel.

Slope is the median of all pairwise slopes between valid observations.
Intercept is `median(y) − slope × median(x)` (Conover formula).
Insensitive to outliers and cloud residuals.

**Returns** `TheilSenResult`:

| Key | Shape | Description |
|---|---|---|
| `slope` | `(rows, cols)` | Robust slope in `[units / day]` |
| `intercept` | `(rows, cols)` | Intercept at the temporal mid-point |

Time complexity is O(n²) per pixel. For stacks larger than ~150 scenes the OLS
slope from `time_window_stats` is faster and typically sufficient.

---

### mann_kendall

```python
def mann_kendall(
    arr:   np.ndarray,
    dates: list[datetime],
) -> MannKendallResult:
```

Non-parametric Mann-Kendall monotonic trend test per pixel (p ≈ 0.05 threshold).

Variance is tie-corrected; the Z statistic uses a continuity correction.

**Returns** `MannKendallResult`:

| Key | Shape | Description |
|---|---|---|
| `S` | `(rows, cols)` | Raw S statistic |
| `varS` | `(rows, cols)` | Variance of S (tie-corrected) |
| `Z` | `(rows, cols)` | Standardised test statistic |
| `trend` | `(rows, cols)` | `+1.0` increasing · `-1.0` decreasing · `0.0` no significant trend |

`NODATA` for pixels with fewer than 4 valid observations.

Prefer Mann-Kendall over OLS slope when the series violates normality assumptions
or the trend is non-linear (stepped change). Combine with `trend_theil_sen` for the
magnitude of a detected trend.

---

### bfast_breakpoint

```python
def bfast_breakpoint(
    arr:   np.ndarray,
    dates: list[datetime],
) -> BreakpointResult:
```

Single structural break detection — simplified BFAST exhaustive scan.

Finds the split of the valid-observation sequence that minimises the combined OLS
residual sum of squares of two line segments.

**Returns** `BreakpointResult`:

| Key | Shape | Description |
|---|---|---|
| `break_day` | `(rows, cols)` | Day number of the detected break |
| `magnitude` | `(rows, cols)` | Signed jump in fitted value at the break (same units as `arr`) |
| `rss_ratio` | `(rows, cols)` | `RSS_two_segments / RSS_one_segment` |

`rss_ratio` interpretation:

| Value | Meaning |
|---|---|
| Near 0 | Sharp structural break — two segments fit much better than one |
| Near 1 | No meaningful break — single line is just as good |

`NODATA` for pixels with fewer than 6 valid observations.
Rule of thumb: `rss_ratio < 0.7` indicates a meaningful structural change.

---

### pixel_regression

```python
def pixel_regression(
    y:      np.ndarray,
    x:      np.ndarray,
    *,
    nodata: float = NODATA,
) -> RegressionResult:
```

Per-pixel OLS linear regression against an explicit predictor vector.

Unlike `time_window_stats` (which always uses day offsets as the predictor),
here `x` is caller-supplied and fully decoupled from dates — pass day-of-year,
cumulative temperature, precipitation index, or any other 1-D predictor.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `y` | `(n_times, rows, cols)` | — | Dependent variable stack |
| `x` | `(n_times,)` | — | Predictor — same for every pixel; must not contain `NODATA` |
| `nodata` | `float` | `-9999.0` | Sentinel value in `y`; mapped to `NODATA` internally |

**Returns** `RegressionResult`:

| Key | Shape | Description |
|---|---|---|
| `slope` | `(rows, cols)` | OLS slope in `[units / x-unit]` |
| `intercept` | `(rows, cols)` | Regression intercept |
| `r2` | `(rows, cols)` | Coefficient of determination R² in `[0, 1]` |

**Edge cases**

| Condition | Result |
|---|---|
| `var(x over valid steps) < 1e-12` | `slope=0`, `intercept=mean(y)`, `r2=0` |
| Fewer than 2 valid `y` observations | All outputs `= NODATA` |

Use `r2 > 0.5` to mask pixels where the linear model is meaningful before
visualising slope maps.

**Raises** `ValueError` when `x.shape != (n_times,)` or `y.ndim != 3`.

---

### phenology_doy

```python
def phenology_doy(
    arr:         np.ndarray,
    dates:       list[datetime],
    rising_pct:  int = 20,
    falling_pct: int = 20,
) -> PhenologyResult:
```

Low-level phenology extractor — thin wrapper around the Fortran `phenology_doy`
subroutine. For most use cases prefer `phenology_metrics` (see below), which
adds built-in Savitzky-Golay smoothing and returns four metrics including
`peak_val`.

**Returns** `PhenologyResult` — keys `sos`, `pos`, `eos`, each `(rows, cols)`.

`NODATA` when `n_valid < 4`, the signal is flat (`amplitude < 1e-6`), or no
threshold crossing is found.

---

### phenology_metrics

```python
def phenology_metrics(
    ndvi_stack:   np.ndarray | xr.DataArray,
    dates:        list[datetime],
    *,
    smooth:       bool = True,
    savgol_window: int = 5,
    rising_pct:   int = 20,
    falling_pct:  int = 20,
    min_valid:    int = 5,
) -> PhenologyMetricsResult:
```

Extract per-pixel phenological season boundaries from a smoothed vegetation index
stack. Implements the TIMESAT amplitude-threshold convention.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ndvi_stack` | `np.ndarray` or `xr.DataArray` | — | `(n_times, rows, cols)` vegetation index (NDVI, EVI, LAI …) |
| `dates` | `list[datetime]` | — | Acquisition timestamps |
| `smooth` | `bool` | `True` | Apply Fortran Savitzky-Golay before threshold detection |
| `savgol_window` | `int` | `5` | S-G window width (odd, ≥ 3; forced odd if even is passed). Ignored when `smooth=False` |
| `rising_pct` | `int` | `20` | SOS threshold — % of seasonal amplitude above the minimum (TIMESAT default) |
| `falling_pct` | `int` | `20` | EOS threshold — same convention; may differ from `rising_pct` |
| `min_valid` | `int` | `5` | Minimum valid (non-NODATA) observations; pixels below this threshold return `NODATA` in all outputs |

**Returns** `PhenologyMetricsResult`:

| Key | Shape | Description |
|---|---|---|
| `sos_doy` | `(rows, cols)` | Start of Season — fractional days from `dates[0]` |
| `eos_doy` | `(rows, cols)` | End of Season — same units |
| `peak_doy` | `(rows, cols)` | Day of peak |
| `peak_val` | `(rows, cols)` | Peak vegetation index value |

Convert back to calendar dates:

```python
from datetime import timedelta
sos_date = dates[0] + timedelta(days=float(metrics["sos_doy"][r, c]))
```

**Season boundary definitions**

| Metric | Definition |
|---|---|
| SOS | First upward crossing of `min + rising_pct% × amplitude` before the peak |
| EOS | Last downward crossing of `min + falling_pct% × amplitude` after the peak |
| POS | Day of the global maximum in the smoothed series |

`NODATA` for: fewer than `min_valid` valid obs, flat signal (`amplitude < 1e-6`),
or no threshold crossing found on a given side.

Pass a **gap-filled** stack (from `interpolate_gaps`) for reliable results.
Cloud gaps produce spurious threshold crossings even after smoothing.

---

### save_phenology

```python
def save_phenology(
    metrics:    PhenologyMetricsResult,
    output_dir: str | Path,
    *,
    crs:        str | None = None,
    transform:  Any | None = None,
    nodata:     float = NODATA,
    compress:   str = "deflate",
) -> dict[str, Path]:
```

Write phenology metrics to four single-band GeoTIFFs. Requires `rasterio`.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `metrics` | `PhenologyMetricsResult` | — | Dict returned by `phenology_metrics` |
| `output_dir` | `str \| Path` | — | Output directory; created if it does not exist |
| `crs` | `str \| None` | `None` | EPSG code or WKT string, e.g. `"EPSG:32633"`. `None` omits spatial reference |
| `transform` | `affine.Affine \| None` | `None` | Geo-transform. `None` → pixel-coordinate identity |
| `nodata` | `float` | `-9999.0` | No-data value written to GeoTIFF metadata |
| `compress` | `str` | `"deflate"` | GDAL codec: `"deflate"` / `"lzw"` / `"none"` |

**Returns** `dict[str, Path]` — `{"sos_doy": Path, "eos_doy": Path, "peak_doy": Path, "peak_val": Path}`.

**Raises** `ImportError` if `rasterio` is not installed.

---

### pearson_map

```python
def pearson_map(
    arr_a: np.ndarray,
    arr_b: np.ndarray,
    dates: list[datetime],
) -> np.ndarray:
```

Pixel-wise Pearson r between two co-registered time stacks.

A time step contributes only when both stacks have a non-NODATA value at that
pixel — no imputation is performed.

**Returns** `(rows, cols)` `float64` — Pearson r in `[−1, 1]`.
`NODATA` when fewer than 3 joint valid observations.

**Raises** `ValueError` when `arr_b.shape != arr_a.shape`.

---

## Quality pipelines

### Coverage and gap check

Run these before any statistical analysis to decide which pixels are usable:

```
valid_obs_count    ──  fraction < 0.3  →  mask pixel
temporal_gap_stats ──  max_gap > 60 d  →  flag pixel (interpolation risky)
pixel_iqr          ──  outlier == 1.0  →  rebuild mask → re-run interpolate_gaps
```

### Trend analysis

```
time_window_stats   ──  mean / std / OLS slope   (fast, O(n) per pixel)
mann_kendall        ──  significance test, p ≈ 0.05
trend_theil_sen     ──  robust slope + intercept  (O(n²), use for < 150 scenes)
```

Use `mann_kendall` to decide where a trend is statistically significant, then
`trend_theil_sen` (or `pixel_regression`) for its magnitude.

### Anomaly detection

```
time_window_stats   ──  compute background mean / std
anomaly_zscore      ──  z = (val − mean) / std per scene
                         |z| > 2.0  ≈ p < 0.05  per pixel
```

### Change detection

```
bfast_breakpoint    ──  rss_ratio < 0.7  →  structural break detected
                        magnitude > 0    →  abrupt increase (e.g. deforestation recovery)
                        magnitude < 0    →  abrupt decrease (e.g. fire, deforestation)
```

---

## Fortran subroutines

All 14 C-exported subroutines live in `sentinel_stats.f90` in a single module `sentinel_stats`.

| Subroutine | Python function | Group |
|---|---|---|
| `valid_obs_count` | `valid_obs_count` | Quality |
| `temporal_gap_stats` | `temporal_gap_stats` | Quality |
| `pixel_quantiles` | `pixel_quantiles` | Distribution |
| `pixel_iqr` | `pixel_iqr` | Distribution |
| `time_window_stats` | `time_window_stats` | Core statistics |
| `anomaly_zscore` | `anomaly_zscore` | Core statistics |
| `trend_theil_sen` | `trend_theil_sen` | Trend |
| `mann_kendall` | `mann_kendall` | Trend test |
| `bfast_breakpoint` | `bfast_breakpoint` | Change detection |
| `pixel_regression` | `pixel_regression` | Regression |
| `phenology_doy` | `phenology_doy` | Phenology |
| `phenology_doy_v2` | `phenology_metrics` (internal) | Phenology |
| `savgol_smooth_stack` | `phenology_metrics` (internal) | Phenology |
| `pearson_map` | `pearson_map` | Correlation |

The internal flat C/row-major array layout:

```
arr[t, r, c]  →  flat index  =  t * rows * cols  +  r * cols  +  c   (1-based in Fortran)
```

All `dates_days` vectors use fractional day offsets from an arbitrary origin — only
differences between values are used, so the origin does not matter.

---

## Examples

### Coverage and quality check before analysis

```python
import numpy as np
from sentinel_processor.analysis._sentinel_stats_bridge import (
    valid_obs_count, temporal_gap_stats, pixel_iqr, NODATA,
)

# arr: (n_times, rows, cols) float64 gap-filled stack
oc = valid_obs_count(arr)
gs = temporal_gap_stats(arr, dates)
iq = pixel_iqr(arr)

# pixels safe to analyse
usable = (oc["fraction"] >= 0.3) & (gs["max_gap"] != NODATA) & (gs["max_gap"] <= 60)

# rebuild mask to re-fill IQR outliers
mask_clean = (iq["outlier"] == 0.0).astype(np.int32)
```

### Rolling statistics and anomaly map

```python
from sentinel_processor.analysis._sentinel_stats_bridge import (
    time_window_stats, anomaly_zscore, NODATA,
)

bg = time_window_stats(ndvi_stack, dates, window_days=90)
z  = anomaly_zscore(ndvi_stack, dates, background=bg)

# map of scenes with significant negative anomalies (drought / browning)
drought_mask = (z < -2.0) & (z != NODATA)   # shape (n_times, rows, cols)
print(f"Anomaly pixels per scene: {drought_mask.sum(axis=(1, 2))}")
```

### Trend significance map

```python
from sentinel_processor.analysis._sentinel_stats_bridge import (
    mann_kendall, trend_theil_sen,
)

mk  = mann_kendall(ndvi_stack, dates)
ts  = trend_theil_sen(ndvi_stack, dates)

# significant increasing pixels only
sig_increase = mk["trend"] == 1.0
print(f"Greening pixels: {sig_increase.sum():,}")
print(f"Mean slope on those pixels: {ts['slope'][sig_increase].mean():.5f} NDVI/day")
```

### OLS regression against explicit predictor

```python
import numpy as np
from sentinel_processor.analysis._sentinel_stats_bridge import pixel_regression

# Regress NIR reflectance against cumulative precipitation
precip_cumsum = np.array([...])   # (n_times,) mm

reg = pixel_regression(nir_stack, precip_cumsum)

# meaningful fit: R² > 0.5
good = reg["r2"] > 0.5
print(f"Pixels with R²>0.5: {good.sum():,} ({100*good.mean():.1f} %)")
```

### Structural break — forest disturbance

```python
from sentinel_processor.analysis._sentinel_stats_bridge import bfast_breakpoint

bp = bfast_breakpoint(ndvi_stack, dates)

# detect significant drops (deforestation, fire)
disturbed = (bp["rss_ratio"] < 0.7) & (bp["magnitude"] < -0.1)
print(f"Disturbance pixels: {disturbed.sum():,}")
```

### Phenology extraction

```python
from datetime import timedelta
from sentinel_processor.analysis._sentinel_stats_bridge import (
    phenology_metrics, save_phenology, NODATA,
)

# ndvi_stack: (n_times, rows, cols) — gap-filled, normalised [0, 1]
metrics = phenology_metrics(
    ndvi_stack,
    dates,
    smooth       = True,
    savgol_window = 7,
    rising_pct   = 20,
    falling_pct  = 20,
)

# convert peak day to calendar date
valid = metrics["peak_doy"] != NODATA
peak_dates = [
    dates[0] + timedelta(days=float(d))
    for d in metrics["peak_doy"][valid].ravel()
]
print(f"Median peak date: {sorted(peak_dates)[len(peak_dates)//2].strftime('%d %b')}")

# save all four maps
saved = save_phenology(metrics, "analysis_output/phenology", crs="EPSG:32633")
for name, path in saved.items():
    print(f"  {name}: {path}")
```

### Full end-to-end pipeline

```python
import numpy as np
from pathlib import Path
from sentinel_processor import stack_timeseries, TimeSeriesConfig
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps
from sentinel_processor.analysis._sentinel_stats_bridge import (
    valid_obs_count, temporal_gap_stats,
    time_window_stats, anomaly_zscore,
    mann_kendall, trend_theil_sen,
    pixel_regression, phenology_metrics, save_phenology,
    NODATA,
)

# 1. Stack scenes
result = stack_timeseries(
    sources = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir = "data/technical",
    cfg     = TimeSeriesConfig(fill_rejected=True, nodata=float("nan")),
)
dates = [s.timestamp for s in result.scenes]

# 2. Gap-fill each band
arr  = result.stack.values.astype(np.float64)
arr[~np.isfinite(arr)] = NODATA
for b in range(arr.shape[1]):
    arr[:, b] = interpolate_gaps(arr[:, b], (arr[:, b] != NODATA).astype(np.int32), method="pchip")
np.clip(arr, 0, None, out=arr)   # clip negative DN artefacts

# 3. Quality check
nir_idx = list(result.stack.coords["band"].values).index("nir")
nir     = arr[:, nir_idx]

oc = valid_obs_count(nir)
gs = temporal_gap_stats(nir, dates)
print(f"Mean coverage   : {oc['fraction'].mean():.1%}")
print(f"Max gap (median): {np.median(gs['max_gap'][gs['max_gap'] != NODATA]):.1f} days")

# 4. Trend analysis
mk  = mann_kendall(nir, dates)
ts  = trend_theil_sen(nir, dates)
sig = mk["trend"] != 0.0
print(f"Pixels with significant trend: {sig.sum():,} ({100*sig.mean():.1f} %)")

# 5. Anomaly z-scores
bg = time_window_stats(nir, dates, window_days=90)
z  = anomaly_zscore(nir, dates, background=bg)

# 6. Phenology from NDVI
ri  = list(result.stack.coords["band"].values).index("red")
ndvi = np.where(
    (nir + arr[:, ri]) > 0,
    (nir - arr[:, ri]) / (nir + arr[:, ri]),
    NODATA,
)
metrics = phenology_metrics(ndvi / 10_000.0, dates)   # normalise DN → [0,1]
save_phenology(metrics, "analysis_output/phenology", crs="EPSG:32633")
```

---

## Method selection guide

```
Coverage / quality           →  valid_obs_count · temporal_gap_stats
Outlier detection            →  pixel_iqr  (feed mask back to interpolate_gaps)
Distribution overview        →  pixel_quantiles  (p10 / p50 / p90, amplitude)
Fast rolling statistics      →  time_window_stats  (mean / std / slope, O(n))
Anomaly map per scene        →  anomaly_zscore  (|z| > 2.0 ≈ p < 0.05)
Trend significance test      →  mann_kendall
Robust trend magnitude       →  trend_theil_sen  (< ~150 scenes)
Fast trend magnitude         →  time_window_stats slope  or  pixel_regression
Arbitrary predictor          →  pixel_regression  (x = precipitation, DOY, …)
Structural break detection   →  bfast_breakpoint  (rss_ratio < 0.7)
Phenology (smoothed NDVI)    →  phenology_metrics + save_phenology
Band-to-band correlation     →  pearson_map
```