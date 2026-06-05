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