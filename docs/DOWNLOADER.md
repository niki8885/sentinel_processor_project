# Downloader

## Install

```bash
pip install sentinel-processor
pip install "sentinel-processor[netcdf]"   # for .nc output
```

## Import

All public symbols are available at the top-level package:

```python
import sentinel_processor as sp

sp.download_sentinel2(...)
sp.DownloadConfig(...)
sp.LocationSpec(...)
sp.SpectralBands.ALL
sp.TechnicalLayers.ALL
sp.validate_file(...)
```

## Module layout

```
sentinel_processor/
├── __init__.py
├── config.py
├── input/
│   └── downloader.py
├── utils/
│   └── data_utils.py
├── processing/
│   ├── _fortran_bridge.py       ← pansharpening bridge
│   ├── _raster_ops_bridge.py    ← raster ops bridge
│   └── fortran/
│       ├── pansharpening.f90
│       ├── raster_ops.f90
│       ├── libsentinel_processing.dll / .so
│       └── libsentinel_raster_ops.dll / .so
└── validation/
    ├── _fortran_bridge.py
    └── fortran/
        ├── validation.f90
        └── libsentinel_validation.dll / .so
```

## Output structure

```
<output_dir>/
├── spectral/
│   ├── <name>_<timestamp>.tif
│   ├── <name>_<timestamp>.nc          ← requires [netcdf] extra
│   └── <name>_<timestamp>_report.json ← if save_report=True
├── technical/
│   ├── scl_<name>_<timestamp>.tif/.nc
│   ├── aot_<name>_<timestamp>.tif/.nc
│   └── wvp_<name>_<timestamp>.tif/.nc
└── visual/
    └── vis_<name>_<timestamp>.tif/.nc
```

File naming: `<name>_<YYYYMMDDTHHMMSS>` when `LocationSpec.name` is set,
`<YYYY-MM-DD>_<lat>_<lon>` otherwise. Files are always overwritten.

## API reference

### LocationSpec

```python
@dataclass
class LocationSpec:
    lon:  float
    lat:  float
    name: str | None = None
```

### DownloadConfig

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bands` | `_BandGroup \| list[str]` | `SpectralBands.ALL` | Spectral bands to download |
| `tech_bands` | same \| `None` | `TechnicalLayers.ALL` | Quality layers; `None` disables |
| `visual` | `bool` | `True` | Download RGB visual overview |
| `output_dir` | `str` | `"data"` | Root output directory |
| `bbox_half_deg` | `float` | `0.05` | Half-size of bounding box in degrees |
| `max_items` | `int` | `50` | Max STAC items per location search |
| `keep_items` | `int` | `10` | Most-recent items to process per location |
| `start_date` | `datetime \| None` | now − lookback_days | Search window start |
| `end_date` | `datetime \| None` | now | Search window end |
| `lookback_days` | `int` | `30` | Used when `start_date` is `None` |
| `validate` | `bool` | `True` | Run Fortran quality checks before saving |
| `band_workers` | `int` | `4` | Parallel threads per scene (bands, tech layers, visuals) |
| `scene_workers` | `int` | `2` | Parallel threads across scenes |
| `pansharpen_algorithm` | `"gram_schmidt" \| "ihs" \| "wavelet" \| None` | `None` | Pansharpening algorithm; `None` disables |
| `pansharpen_pan_key` | `str` | `"visual"` | STAC asset key used as the PAN source |
| `max_cloud_threshold` | `float` | `0.30` | Max allowed cloud fraction (0–1) |
| `min_confidence` | `float` | `0.01` | Min confidence score to accept scene |
| `save_report` | `bool` | `True` | Write JSON validation sidecar per scene |

### download_sentinel2

```python
def download_sentinel2(
    locations: Sequence[LocationSpec],
    cfg: DownloadConfig | None = None,
    progress: bool = True,
) -> dict[str, list[str]]:
```

Searches the Microsoft Planetary Computer STAC API for Sentinel-2 L2A scenes,
validates them, and writes spectral, technical, and visual rasters to disk.

Returns `dict[base_name → list_of_written_paths]`.
Rejected scenes return an empty list and are not saved.
`progress=False` suppresses the terminal progress bar.

## Processing pipeline

```
STAC search
    │
    ▼
per scene  (scene_workers threads)
    │
    ├── fetch spectral bands   ──  band_workers threads, parallel
    ├── fetch tech layers      ──  band_workers threads, parallel
    └── fetch visual assets    ──  band_workers threads, parallel
    │
    ▼
align bands to reference grid
    ├── same CRS   →  reproject_nearest()    Fortran  [raster_ops]
    └── other CRS  →  reproject_match()      rioxarray fallback
    │
    ▼
pansharpening  ──  optional, Fortran  [pansharpening]
    │   ├── rgb_to_luminance()   collapse visual RGB → PAN   [raster_ops]
    │   └── gram_schmidt / ihs / wavelet
    │
    ▼
validation  ──  optional, Fortran  [validation]
    │   check_dimensions → validate_scl → check_radiometry
    │
    └── pass  →  save .tif / .nc + report.json
```

See [PANSHARPENING.md](PANSHARPENING.md) for pansharpening details.
See [RASTER_OPS.md](RASTER_OPS.md) for individual Fortran raster routine details.
See [VALIDATION.md](VALIDATION.md) for validation thresholds and Fortran routines.

## Fortran acceleration summary

| Operation | Before | After |
|---|---|---|
| Align 20 m / 60 m bands to 10 m reference | `reproject_match` (GDAL warp) | `reproject_nearest` (Fortran) |
| Align inside `_clip` to reference grid | `reproject_match` (GDAL warp) | `reproject_nearest` (Fortran) |
| Collapse RGB visual → PAN luminance | `np.tensordot` (Python) | `rgb_to_luminance` (Fortran) |
| Pansharpening algorithms | — | `gram_schmidt` / `ihs` / `wavelet` (Fortran) |
| SCL cloud/snow validation | — | `validate_scl` (Fortran) |
| Dimension check | — | `check_dimensions` (Fortran) |
| Radiometry check | — | `check_radiometry` (Fortran) |

All Fortran paths have NumPy/rioxarray fallbacks — if a library is not compiled
the downloader continues to work correctly, only slower.

## Validation pipeline

When `validate=True` each scene passes three checks before being saved.
A failure at any stage rejects the scene immediately.

```
SCL raster
    │
    ▼
check_dimensions()   ── min side ≥ 32 px, max side ≤ 10 980 px,
    │                    aspect ratio ≤ 4:1  (rejects e.g. 15×1152)
    ▼
validate_scl()       ── cloud / snow analysis → confidence_score
    │
    ▼
check_radiometry()   ── < 1% pixels above 15 000 DN
    │
    └── all pass → files saved
```

See [VALIDATION.md](VALIDATION.md) for full details.

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
  "item_id": "S2A_34TCT_20260510_1_L2A",
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

## Band presets

```python
SpectralBands.RGB           # blue, green, red
SpectralBands.RGB_NIR       # blue, green, red, nir
SpectralBands.VEGETATION    # red, nir, rededge1, rededge2, rededge3
SpectralBands.AGRICULTURE   # blue, green, red, nir, rededge1, swir16, swir22
SpectralBands.ALL_10M       # blue, green, red, nir
SpectralBands.ALL_20M       # rededge1-3, nir08, swir16, swir22
SpectralBands.ALL           # ALL_10M + ALL_20M

TechnicalLayers.SCL         # Scene Classification Layer only
TechnicalLayers.ALL         # scl + aot + wvp
```

Custom list: `DownloadConfig(bands=["red", "nir", "swir16"])`

## Examples

### Minimal

```python
import sentinel_processor as sp

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17)],
)
```

### Custom window and bands

```python
import datetime
import sentinel_processor as sp

cfg = sp.DownloadConfig(
    bands               = sp.SpectralBands.VEGETATION,
    tech_bands          = None,
    visual              = False,
    output_dir          = "/mnt/sentinel",
    start_date          = datetime.datetime(2025, 1, 1, tzinfo=datetime.UTC),
    end_date            = datetime.datetime(2025, 6, 1, tzinfo=datetime.UTC),
    validate            = True,
    max_cloud_threshold = 0.15,
    min_confidence      = 0.75,
)

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="my_field")],
    cfg=cfg,
)
```

### Multiple locations

```python
import sentinel_processor as sp

results = sp.download_sentinel2(
    [
        sp.LocationSpec(lat=47.56, lon=19.17, name="budapest"),
        sp.LocationSpec(lat=52.52, lon=13.40, name="berlin"),
    ],
    cfg=sp.DownloadConfig(bands=sp.SpectralBands.RGB_NIR),
)
```

### With pansharpening

```python
import sentinel_processor as sp

# Gram-Schmidt — any number of bands
cfg = sp.DownloadConfig(
    bands                = sp.SpectralBands.ALL,
    pansharpen_algorithm = "gram_schmidt",
)

# IHS — exactly 3 bands
cfg = sp.DownloadConfig(
    bands                = sp.SpectralBands.RGB,
    pansharpen_algorithm = "ihs",
)

# Wavelet — any number of bands
cfg = sp.DownloadConfig(
    bands                = sp.SpectralBands.AGRICULTURE,
    pansharpen_algorithm = "wavelet",
)

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=cfg,
)
```

See [PANSHARPENING.md](PANSHARPENING.md) for algorithm details and build instructions.

### Tune parallelism

```python
cfg = sp.DownloadConfig(
    band_workers  = 8,   # more threads per scene (fast connection)
    scene_workers = 1,   # one scene at a time (memory-constrained)
)
```

`band_workers` controls parallel downloads within a single scene (spectral bands,
tech layers, and visual assets each get their own pool of this size).
`scene_workers` controls how many scenes are processed concurrently.
Total thread ceiling is roughly `scene_workers × band_workers × 3`.

### Disable validation

```python
import sentinel_processor as sp

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=sp.DownloadConfig(validate=False),
)
```