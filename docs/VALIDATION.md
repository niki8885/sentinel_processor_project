# Validation

## Overview

Validation runs automatically before saving, controlled by `DownloadConfig.validate`.
It uses three Fortran routines compiled into `libsentinel_validation.dll` (Windows) or `.so` (Linux/Mac).

```
SCL raster
    │
    ▼
check_dimensions()  ──  min/max side, aspect ratio
    │                   rejects degenerate clips like 15×1152
    ├── fail  ──►  rejected, nothing saved
    └── pass
            │
            ▼
        validate_scl()  ──  confidence_score, cloud_ratio, snow_ratio, issues[]
            │
            ├── confidence_score < min_confidence  ──►  rejected
            └── pass
                    │
                    ▼
                check_radiometry()  ──  bool  (< 1% pixels > 15 000 DN)
                    │
                    ├── False  ──►  rejected
                    └── True   ──►  files saved
```

## Install

The Fortran library must be compiled before validation can run.
Without it, set `validate=False` — everything else works normally.

**Linux / Mac**
```bash
cd sentinel_processor/validation/fortran
gfortran -O2 -shared -fPIC -o libsentinel_validation.so validation.f90
```

**Windows** (MSYS2 UCRT64 terminal)
```bat
cd sentinel_processor\validation\fortran
gfortran -O2 -shared -o libsentinel_validation.dll validation.f90
```

> Do not use `-static-libgfortran -static-libgcc` on Windows — it causes linker errors with `strndup` on GCC 16+.

## Fortran routines

### check_dimensions

Validates the pixel dimensions of the clipped SCL raster.
Runs first; a failure skips all further checks.

**Thresholds (compile-time constants in `validation.f90`)**

| Constant | Value | Meaning |
|---|---|---|
| `MIN_SIDE` | 32 px | Minimum pixels on each side |
| `MAX_SIDE` | 10 980 px | Maximum pixels on each side (full Sentinel-2 tile) |
| `MAX_ASPECT_RATIO` | 4.0 | Maximum ratio of longer side to shorter side |

A clip of 15×1152 has an aspect ratio of ~77 and is rejected immediately.

**Issues emitted**

| Message | Cause |
|---|---|
| `"Scene too small: side below minimum threshold"` | rows < 32 or cols < 32 |
| `"Scene too large: side exceeds maximum threshold"` | rows > 10 980 or cols > 10 980 |
| `"Degenerate shape: aspect ratio exceeds limit"` | max(rows,cols) / min(rows,cols) > 4.0 |

### validate_scl

Analyses SCL pixel values and returns a confidence score.

**SCL value mapping**

| Value | Meaning | Treatment |
|---|---|---|
| 0, 6 | No data / water | Excluded from ratio calculation |
| 3, 8, 9, 10 | Cloud shadow / cloud (med/high/cirrus) | Counted as bad |
| 11 | Snow / ice | Counted as snow |
| 1, 2, 4, 5, 7 | Defective, dark area, vegetation, bare soil, unclassified | Counted as valid |

**Confidence score**

| Condition | Score |
|---|---|
| No valid pixels after filtering | 0.0 |
| snow_ratio > 0.50 | 0.0 |
| cloud_ratio > max_cloud_thr | 0.0 |
| cloud_ratio < 0.10 | 1.0 |
| cloud_ratio < 0.30 | 0.75 |
| cloud_ratio < 0.40 | 0.50 |
| otherwise | 0.0 |

**Issues emitted**

| Message | Cause |
|---|---|
| `"No valid pixels after filtering"` | filtered_total == 0 |
| `"High cloud cover"` | cloud_ratio > 0.30 |
| `"Excessive snow/ice"` | snow_ratio > 0.50 |

### check_radiometry

Returns `True` when fewer than 1% of pixels exceed **15 000 DN**.
Applied to the SCL band values of the clipped raster.

## Python API

### Validate during download

Controlled via `DownloadConfig`:

```python
import sentinel_processor as sp

cfg = sp.DownloadConfig(
    validate            = True,
    max_cloud_threshold = 0.20,   # stricter than default 0.30
    min_confidence      = 0.75,
    save_report         = True,
)
```

### Validate an existing file

```python
from sentinel_processor import validate_file

report = validate_file(
    "downloads/technical/scl_budapest_20260501T103000.tif",
    max_cloud_threshold = 0.30,
    min_confidence      = 0.50,
)

print(report["passed"])            # True / False
print(report["confidence_score"])  # 1.0
print(report["cloud_ratio"])       # 0.04
print(report["issues"])            # []
print(report["radiometry_pass"])   # True
```

Accepts any SCL GeoTIFF or NetCDF file. Does not require a STAC connection.
Note: `validate_file` does not run `check_dimensions` — it operates on a file you have already chosen.

### Low-level calls

```python
import rioxarray
from sentinel_processor import call_validate_scl, call_check_radiometry, call_check_dimensions

da  = rioxarray.open_rasterio("scl.tif")
arr = da.values

# dimension check
dim = call_check_dimensions(arr.shape[-2], arr.shape[-1])
# {"passed": True, "issues": []}

# SCL quality check
scl_result = call_validate_scl(arr.flatten().astype(int).tolist(), max_cloud_threshold=0.30)
# {"confidence_score": 1.0, "cloud_ratio": 0.04, "snow_ratio": 0.0,
#  "water_excluded": True, "issues": []}

# radiometry check
radio = call_check_radiometry(arr.flatten().astype(float))
# True
```

## Validation report JSON

Written to `spectral/<name>_report.json` when `save_report=True`.

**Passed:**
```json
{
  "item_id": "S2B_34TCT_20260501_1_L2A",
  "passed": true,
  "rows": 1100,
  "cols": 1100,
  "dimension_pass": true,
  "confidence_score": 1.0,
  "cloud_ratio": 0.04,
  "snow_ratio": 0.0,
  "water_excluded": true,
  "issues": [],
  "radiometry_pass": true
}
```

**Rejected — degenerate shape:**
```json
{
  "item_id": "S2A_34TCT_20260510_1_L2A",
  "passed": false,
  "rows": 15,
  "cols": 1152,
  "dimension_pass": false,
  "issues": ["Degenerate shape: aspect ratio exceeds limit"]
}
```

**Rejected — clouds:**
```json
{
  "item_id": "S2A_34TCT_20260515_1_L2A",
  "passed": false,
  "rows": 1100,
  "cols": 1100,
  "dimension_pass": true,
  "confidence_score": 0.0,
  "cloud_ratio": 0.61,
  "snow_ratio": 0.0,
  "water_excluded": true,
  "issues": ["High cloud cover"],
  "radiometry_pass": true
}
```

**No SCL asset** (validation skipped, scene saved anyway):
```json
{
  "item_id": "S2X_...",
  "passed": true,
  "warning": "SCL asset missing; validation skipped"
}
```