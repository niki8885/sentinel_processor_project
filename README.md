# sentinel-processor

Sentinel-2 L2A downloader built on the [Element84 STAC API](https://earth-search.aws.element84.com/v1).  
Downloads spectral bands, quality layers, and visual overviews for any coordinate, with cloud/radiometry validation backed by a compiled Fortran library, spectral index computation, and interactive Plotly visualisation.

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

### Fortran libraries

Two Fortran kernels must be compiled before validation and index computation are available.  
Without them, set `validate=False` and avoid calling `compute_indices` — everything else works normally.

**Linux / Mac**
```bash
# validation
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/validation/fortran/libsentinel_validation.so \
  sentinel_processor/validation/fortran/validation.f90

# indices
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/indices/fortran/libsentinel_indices.so \
  sentinel_processor/indices/fortran/indices_mod.f90
```

**Windows** (MSYS2 UCRT64 terminal, requires `mingw-w64-ucrt-x86_64-gcc-fortran`)
```bash
# validation
gfortran -O2 -shared -static-libgfortran -static-libgcc \
  -o sentinel_processor/validation/fortran/libsentinel_validation.dll \
  sentinel_processor/validation/fortran/validation.f90

# indices
gfortran -O2 -shared -fPIC -static-libgfortran -static-libgcc \
  -o sentinel_processor/indices/fortran/libsentinel_indices.dll \
  sentinel_processor/indices/fortran/indices_mod.f90
```

---

## Quick start

```python
import sentinel_processor as sp
from sentinel_processor.indices.compute import compute_indices
from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask

# 1. Download
results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
)

scene = "data/spectral/budapest_20260526T095725.nc"
vis   = "data/visual/vis_budapest_20260526T095725.nc"
scl   = "data/technical/scl_budapest_20260526T095725.nc"

# 2. Compute indices
idx = compute_indices(scene, ["ndvi", "ndwi", "ndbi"])

# 3. Visualise
plot_rgb(vis).show()
plot_band(idx["ndvi"], colorscale="RdYlGn").show()
plot_mask(scl).show()
plot_grid([
    {"file": scene,        "band": "nir",  "label": "NIR"},
    {"file": idx["ndvi"],                  "label": "NDVI", "colorscale": "RdYlGn"},
    {"file": idx["ndwi"],                  "label": "NDWI", "colorscale": "Blues"},
], ncols=3).show()
```

---

## Modules

| Module | Description | Docs |
|---|---|---|
| `sentinel_processor` | Download, STAC search, validation | [DOWNLOADER.md](docs/DOWNLOADER.md) |
| `sentinel_processor.indices` | Spectral index computation (Fortran kernel) | [INDICES.md](docs/INDICES.md) |
| `sentinel_processor.validation` | SCL + radiometry quality checks (Fortran) | [VALIDATION.md](docs/VALIDATION.md) |
| `sentinel_processor.visualisation` | Interactive Plotly figures | [VISUALISATION.md](docs/VISUALISATION.md) |

---

## Usage

### Download with defaults

```python
import sentinel_processor as sp

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="my_field")],
)
```

Downloads the 10 most recent cloud-free Sentinel-2 scenes within 5 km of the point.  
Saves to `data/spectral/`, `data/technical/`, `data/visual/`.

### Custom config

```python
import datetime
import sentinel_processor as sp

cfg = sp.DownloadConfig(
    bands               = sp.SpectralBands.VEGETATION,
    tech_bands          = sp.TechnicalLayers.SCL,
    visual              = False,
    output_dir          = "/mnt/sentinel",
    start_date          = datetime.datetime(2025, 3, 1, tzinfo=datetime.UTC),
    end_date            = datetime.datetime(2025, 6, 1, tzinfo=datetime.UTC),
    keep_items          = 3,
    max_cloud_threshold = 0.15,
    min_confidence      = 0.75,
)

results = sp.download_sentinel2(
    [
        sp.LocationSpec(lat=47.56, lon=19.17, name="budapest"),
        sp.LocationSpec(lat=52.52, lon=13.40, name="berlin"),
    ],
    cfg=cfg,
)
```

### Compute spectral indices

```python
from sentinel_processor.indices.compute import compute_indices

results = compute_indices(
    source="data/spectral/budapest_20260526T095725.nc",
    indices=["ndvi", "evi", "ndwi", "ndbi", "nbr"],
    output_dir="data/indices",
)
# {"ndvi": "data/indices/indices_budapest_..._ndvi.tif", ...}
```

Supported indices: `ndvi`, `evi`, `savi`, `ndwi`, `mndwi`, `ndbi`, `nbr`, `ndsi`, `cig`, `arvi`.

### Visualise

```python
from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask

plot_band("data/spectral/scene.nc", band="nir", colorscale="Plasma").show()
plot_rgb("data/visual/vis_scene.nc").show()
plot_rgb("data/spectral/scene.nc", "nir", "red", "green").show()   # false colour
plot_mask("data/technical/scl_scene.nc", bad_classes=[8, 9, 10]).show()
plot_grid([
    {"file": "data/spectral/scene.nc", "band": "red",   "label": "Red"},
    {"file": "data/spectral/scene.nc", "band": "nir",   "label": "NIR"},
    {"file": "data/indices/indices_scene_ndvi.tif",     "label": "NDVI", "colorscale": "RdYlGn"},
], ncols=3).show()
```

### Validate an existing SCL file

```python
from sentinel_processor import validate_file

report = validate_file(
    "data/technical/scl_budapest_20260526T095725.tif",
    max_cloud_threshold=0.30,
    min_confidence=0.50,
)
print(report["passed"])           # True
print(report["cloud_ratio"])      # 0.04
print(report["confidence_score"]) # 1.0
```

---

## Output structure

```
<output_dir>/
├── spectral/
│   ├── <name>_<timestamp>.tif
│   ├── <name>_<timestamp>.nc
│   ├── <name>_<timestamp>_report.json
│   └── indices/
│       └── indices_<name>_<timestamp>_<index>.tif
├── technical/
│   ├── scl_<name>_<timestamp>.tif/.nc
│   ├── aot_<name>_<timestamp>.tif/.nc
│   └── wvp_<name>_<timestamp>.tif/.nc
└── visual/
    └── vis_<name>_<timestamp>.tif/.nc
```

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

Custom band list: `DownloadConfig(bands=["red", "nir", "swir16"])`

---

## Requirements

- Python ≥ 3.10
- `pystac-client`, `rioxarray`, `xarray`, `numpy`, `rasterio`, `plotly`
- `netCDF4` or `h5netcdf` for `.nc` output (optional)
- `gfortran` ≥ 9 to build the Fortran validation and indices libraries

---

## Development

```bash
git clone https://github.com/niki8885/sentinel-processor
cd sentinel-processor
pip install -e ".[dev]"
```