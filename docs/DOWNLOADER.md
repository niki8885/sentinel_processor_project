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
├── __init__.py          ← public API
├── config.py            ← constants & defaults
├── input/
│   └── downloader.py
├── utils/
│   └── data_utils.py
└── validation/
    ├── _fortran_bridge.py
    └── fortran/
        └── validation.f90
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

File naming: `<name>_<YYYYMMDDTHHMMSS>` when name is set, `<YYYY-MM-DD>_<lat>_<lon>` otherwise.  
Files are always overwritten.

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
| `bands` | `_BandGroup \| list[str]` | `SpectralBands.ALL` | Spectral bands |
| `tech_bands` | same \| `None` | `TechnicalLayers.ALL` | Quality layers; `None` disables |
| `visual` | `bool` | `True` | RGB visual overview |
| `output_dir` | `str` | `"data"` | Root output directory |
| `bbox_half_deg` | `float` | `0.05` | Half-size of bounding box in degrees |
| `max_items` | `int` | `50` | Max STAC items per search |
| `keep_items` | `int` | `10` | Most-recent items to process |
| `start_date` | `datetime \| None` | now − lookback_days | Search window start |
| `end_date` | `datetime \| None` | now | Search window end |
| `lookback_days` | `int` | `30` | Used when `start_date` is None |
| `validate` | `bool` | `True` | Run Fortran quality checks |
| `max_cloud_threshold` | `float` | `0.30` | Max allowed cloud fraction |
| `min_confidence` | `float` | `0.01` | Min confidence score to accept |
| `save_report` | `bool` | `True` | Write JSON sidecar per scene |

### download_sentinel2

```python
def download_sentinel2(
    locations: Sequence[LocationSpec],
    cfg: DownloadConfig | None = None,
    progress: bool = True,
) -> dict[str, list[str]]:
```

Returns `dict[base_name → list_of_written_paths]`.  
`progress=False` suppresses the terminal progress bar.

## Band presets

```python
SpectralBands.RGB           # blue, green, red
SpectralBands.RGB_NIR       # blue, green, red, NIR
SpectralBands.VEGETATION    # red, NIR, rededge 1-2-3
SpectralBands.AGRICULTURE   # blue, green, red, NIR, rededge1, SWIR1, SWIR2
SpectralBands.ALL_10M       # blue, green, red, NIR
SpectralBands.ALL_20M       # rededge 1-3, NIR narrow, SWIR 1-2
SpectralBands.ALL           # ALL_10M + ALL_20M

TechnicalLayers.SCL         # Scene Classification Layer only
TechnicalLayers.ALL         # SCL + AOT + WVP
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
import datetime, sentinel_processor as sp

cfg = sp.DownloadConfig(
    bands      = sp.SpectralBands.VEGETATION,
    tech_bands = None,
    visual     = False,
    output_dir = "/mnt/sentinel",
    start_date = datetime.datetime(2025, 1, 1, tzinfo=datetime.UTC),
    end_date   = datetime.datetime(2025, 6, 1, tzinfo=datetime.UTC),
    validate   = True,
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