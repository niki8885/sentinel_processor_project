# Validation

## Overview

Validation runs automatically before saving, controlled by `DownloadConfig.validate`.  
It uses two Fortran routines compiled into `libsentinel_validation.dll` (Windows) or `.so` (Linux/Mac).

```
SCL pixels
    │
    ▼
validate_scl()  ──►  confidence_score, cloud_ratio, snow_ratio, issues[]
    │
    ├── confidence_score < min_confidence  ──►  rejected, nothing saved
    └── passes
            │
            ▼
        check_radiometry()  ──►  bool  (< 1% pixels > 15 000 DN)
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

**Windows** (MSYS2 UCRT64 with `mingw-w64-ucrt-x86_64-gcc-fortran`)
```bat
cd sentinel_processor\validation\fortran
gfortran -O2 -shared -o libsentinel_validation.dll validation.f90
```

## Fortran routines

### validate_scl

Analyses SCL pixel values.

**SCL value mapping**

| Value | Meaning | Treatment |
|---|---|---|
| 0, 6 | No data / water | Excluded |
| 3, 8, 9, 10 | Cloud shadow / cloud (med/high/cirrus) | Bad |
| 11 | Snow / ice | Snow |
| others | Vegetation, bare soil, etc. | Valid |

**Confidence score**

| Condition | Score |
|---|---|
| No valid pixels | 0.0 |
| snow_ratio > 0.50 | 0.0 |
| cloud_ratio > max_cloud_thr | 0.0 |
| cloud_ratio < 0.10 | 1.0 |
| cloud_ratio < 0.30 | 0.75 |
| cloud_ratio < 0.40 | 0.50 |
| otherwise | 0.0 |

Issues: `"No valid pixels after filtering"`, `"High cloud cover"`, `"Excessive snow/ice"`

### check_radiometry

Returns `True` when fewer than 1% of pixels exceed **15 000 DN**.

## Python API

### Validate during download

Controlled by `DownloadConfig`:

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

print(report["passed"])          # True / False
print(report["confidence_score"]) # 1.0
print(report["cloud_ratio"])     # 0.04
print(report["issues"])          # []
print(report["radiometry_pass"]) # True
```

Accepts any SCL GeoTIFF or NetCDF file. Does not require a STAC connection.

### Low-level calls

```python
import rioxarray
from sentinel_processor import call_validate_scl, call_check_radiometry

da  = rioxarray.open_rasterio("scl.tif")
scl = da.values.flatten().astype(int).tolist()

result = call_validate_scl(scl, max_cloud_threshold=0.30)
# {"confidence_score": 1.0, "cloud_ratio": 0.04, "snow_ratio": 0.0,
#  "water_excluded": True, "issues": []}

radio = call_check_radiometry(da.values.flatten().astype(float))
# True
```

## Validation report JSON

Written to `spectral/<name>_report.json` when `save_report=True`.

**Passed:**
```json
{
  "item_id": "S2B_34TCT_20260501_1_L2A",
  "passed": true,
  "confidence_score": 1.0,
  "cloud_ratio": 0.04,
  "snow_ratio": 0.0,
  "water_excluded": true,
  "issues": [],
  "radiometry_pass": true
}
```

**Rejected:**
```json
{
  "item_id": "S2A_34TCT_20260510_1_L2A",
  "passed": false,
  "confidence_score": 0.0,
  "cloud_ratio": 0.61,
  "snow_ratio": 0.0,
  "water_excluded": true,
  "issues": ["High cloud cover"],
  "radiometry_pass": true
}
```

**No SCL asset** (validation skipped, item saved anyway):
```json
{
  "item_id": "S2X_...",
  "passed": true,
  "warning": "SCL asset missing; validation skipped"
}
```