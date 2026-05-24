# sentinel-processor

Sentinel-2 L2A downloader built on the [Element84 STAC API](https://earth-search.aws.element84.com/v1).  
Downloads spectral bands, quality layers, and visual overviews for any coordinate, with cloud/radiometry validation backed by a compiled Fortran library.

---

## Installation

```bash
pip install sentinel-processor
```

For NetCDF output (`.nc` files in addition to `.tif`):

```bash
pip install "sentinel-processor[netcdf]"
```

Full optional extras:

```bash
pip install "sentinel-processor[all]"   # netcdf + shapely + geoalchemy2
```

### Fortran validation library

The quality validation step requires a compiled shared library.  
Pre-built binaries are not yet distributed — build from source:

**Linux / Mac**
```bash
cd sentinel_processor/validation/fortran
gfortran -O2 -shared -fPIC -o libsentinel_validation.so validation.f90
```

**Windows** (requires [MSYS2](https://www.msys2.org/) with `mingw-w64-ucrt-x86_64-gcc-fortran`)
```bat
cd sentinel_processor\validation\fortran
gfortran -O2 -shared -o libsentinel_validation.dll validation.f90
```

If the library is not built, `validate=False` bypasses validation and everything else works normally.

---

## Quick start

```python
import sentinel_processor as sp

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.562938, lon=19.169805, name="budapest")],
)
# results["budapest_20260501T103000"] → ["downloads/spectral/budapest_20260501T103000.tif", ...]
```

---

## Usage

### Download with defaults

```python
import sentinel_processor as sp

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="my_field")],
)
```

Downloads the 5 most recent cloud-free Sentinel-2 scenes within 5 km of the point, saves to `data/spectral/`, `data/technical/`, `data/visual/`.

### Custom config

```python
import datetime
import sentinel_processor as sp

cfg = sp.DownloadConfig(
    bands               = sp.SpectralBands.VEGETATION,  # red, NIR, rededge 1-2-3
    tech_bands          = sp.TechnicalLayers.SCL,       # SCL only
    visual              = False,
    output_dir          = "/mnt/sentinel",
    bbox_half_deg       = 0.045,                        # ~5 km
    start_date          = datetime.datetime(2025, 3, 1, tzinfo=datetime.UTC),
    end_date            = datetime.datetime(2025, 6, 1, tzinfo=datetime.UTC),
    keep_items          = 3,
    validate            = True,
    max_cloud_threshold = 0.15,
    min_confidence      = 0.75,
    save_report         = True,
)

results = sp.download_sentinel2(
    [
        sp.LocationSpec(lat=47.56, lon=19.17, name="budapest"),
        sp.LocationSpec(lat=52.52, lon=13.40, name="berlin"),
    ],
    cfg=cfg,
)
```

### Validate an existing SCL file

```python
from sentinel_processor import validate_file

report = validate_file(
    "downloads/technical/scl_budapest_20260501T103000.tif",
    max_cloud_threshold = 0.30,
    min_confidence      = 0.50,
)
print(report["passed"])        # True
print(report["cloud_ratio"])   # 0.04
print(report["issues"])        # []
```

### Low-level validation calls

```python
import rioxarray
from sentinel_processor import call_validate_scl, call_check_radiometry

da = rioxarray.open_rasterio("scl.tif")
pixels = da.values.flatten().astype(int).tolist()

result = call_validate_scl(pixels, max_cloud_threshold=0.30)
radio  = call_check_radiometry(da.values.flatten().astype(float))
```

### Suppress progress bar

```python
results = sp.download_sentinel2(locations, cfg=cfg, progress=False)
```

---

## Output structure

```
<output_dir>/
├── spectral/
│   ├── <name>_<timestamp>.tif
│   ├── <name>_<timestamp>.nc       ← requires sentinel-processor[netcdf]
│   └── <name>_<timestamp>_report.json
├── technical/
│   ├── scl_<name>_<timestamp>.tif/.nc
│   ├── aot_<name>_<timestamp>.tif/.nc
│   └── wvp_<name>_<timestamp>.tif/.nc
└── visual/
    └── vis_<name>_<timestamp>.tif/.nc
```

File naming: `<name>_<YYYYMMDDTHHMMSS>` when `name` is set, otherwise `<YYYY-MM-DD>_<lat>_<lon>`.

---

## Band presets

| Preset | Bands |
|---|---|
| `SpectralBands.RGB` | blue, green, red |
| `SpectralBands.RGB_NIR` | blue, green, red, NIR |
| `SpectralBands.VEGETATION` | red, NIR, rededge 1-2-3 |
| `SpectralBands.AGRICULTURE` | blue, green, red, NIR, rededge1, SWIR1, SWIR2 |
| `SpectralBands.ALL_10M` | blue, green, red, NIR |
| `SpectralBands.ALL_20M` | rededge 1-3, NIR narrow, SWIR 1-2 |
| `SpectralBands.ALL` | all 10 m + 20 m bands |
| `TechnicalLayers.SCL` | Scene Classification Layer |
| `TechnicalLayers.ALL` | SCL + AOT + WVP |

Custom band lists are also accepted:
```python
cfg = sp.DownloadConfig(bands=["red", "nir", "swir16"])
```

---

## DownloadConfig reference

| Parameter | Default | Description |
|---|---|---|
| `bands` | `SpectralBands.ALL` | Spectral bands |
| `tech_bands` | `TechnicalLayers.ALL` | Quality layers; `None` to disable |
| `visual` | `True` | RGB visual overview |
| `output_dir` | `"data"` | Root output directory |
| `bbox_half_deg` | `0.05` | Bounding-box half-size in degrees |
| `max_items` | `50` | Max STAC items per search |
| `keep_items` | `10` | Most-recent items to process |
| `lookback_days` | `30` | Default search window when `start_date` is None |
| `validate` | `True` | Run Fortran quality checks |
| `max_cloud_threshold` | `0.30` | Max allowed cloud fraction |
| `min_confidence` | `0.01` | Min confidence score to save |
| `save_report` | `True` | Write JSON report per scene |

---

## Requirements

- Python ≥ 3.10
- `pystac-client`, `rioxarray`, `xarray`, `numpy`, `rasterio`
- `netCDF4` or `h5netcdf` for `.nc` output (optional)
- `gfortran` ≥ 9 to build the Fortran validation library

---

## Development

```bash
git clone https://github.com/your-org/sentinel-processor
cd sentinel-processor
pip install -e ".[dev]"
python run_test.py
```