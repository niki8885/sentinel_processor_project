# Indices

Spectral index computation for Sentinel-2 scenes.  
Pixel-level math is executed by a compiled Fortran kernel (`indices_mod.f90`) via ctypes; Python handles I/O, band routing, and output formatting.

## Build the Fortran library

The shared library must be compiled once before the module can be used.

**Linux / macOS**
```bash
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/indices/fortran/libsentinel_indices.so \
  sentinel_processor/indices/fortran/indices_mod.f90
```

**Windows (MSYS2 UCRT64 terminal)**
```bash
gfortran -O2 -shared -fPIC -static-libgfortran -static-libgcc \
  -o sentinel_processor/indices/fortran/libsentinel_indices.dll \
  sentinel_processor/indices/fortran/indices_mod.f90
```

> Run from the project root. The library is loaded lazily on first use; a `FileNotFoundError` with build instructions is raised if it is missing.

## Import

```python
from sentinel_processor.indices.compute import compute_indices, list_indices
```

## Module layout

```
sentinel_processor/
└── indices/
    ├── __init__.py
    ├── compute.py           ← public API
    ├── _indices_bridge.py   ← ctypes bridge (internal)
    └── fortran/
        ├── indices_mod.f90
        └── libsentinel_indices.dll / .so
```

## Output structure

Results are written next to the source file by default:

```
<output_dir>/                          ← defaults to <source_dir>/indices/
└── indices_<stem>_<index_name>.tif    ← one file per index
```

For example, processing `downloads/spectral/budapest_20260520T094746.nc` with
`["ndvi", "ndwi"]` produces:

```
downloads/spectral/indices/
├── indices_budapest_20260520T094746_ndvi.tif
└── indices_budapest_20260520T094746_ndwi.tif
```

Output files are always `float32`. Pixels with a zero-sum denominator are written as `0.0`; no-data fill value defaults to `-9999.0`.

## API reference

### compute_indices

```python
def compute_indices(
    source: str | Path,
    indices: Sequence[str],
    output_dir: str | Path | None = None,
    scale_factor: float | None = None,
    nodata: float = -9999.0,
    overwrite: bool = True,
    output_format: str = "tif",
) -> dict[str, str]:
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source` | `str \| Path` | — | Path to a `.tif`, `.tiff`, or `.nc` scene file |
| `indices` | `Sequence[str]` | — | Index names to compute, e.g. `["ndvi", "evi"]` |
| `output_dir` | `str \| Path \| None` | `<source_dir>/indices/` | Directory for output files |
| `scale_factor` | `float \| None` | `None` | Multiplied against raw DN values before computation (e.g. `1e-4` for Sentinel-2 L2A integer reflectance) |
| `nodata` | `float` | `-9999.0` | Fill value written to output rasters |
| `overwrite` | `bool` | `True` | If `False`, existing files are skipped |
| `output_format` | `str` | `"tif"` | `"tif"` or `"nc"` |

Returns `dict[index_name → absolute_output_path]` for every index successfully saved.  
Indices whose required bands are absent in the source are skipped with a `WARNING` log; no exception is raised.

### list_indices

```python
def list_indices() -> dict[str, dict]:
```

Returns a summary of all registered indices and the bands they require.

```python
>>> from sentinel_processor.indices.compute import list_indices
>>> list_indices()
{
  "arvi":  {"bands_required": ["B08", "B04", "B02"], "long_name": "Atmospherically Resistant Vegetation Index"},
  "cig":   {"bands_required": ["B08", "B03"], "long_name": "Chlorophyll Index Green"},
  "evi":   {"bands_required": ["B08", "B04", "B02"], "long_name": "Enhanced Vegetation Index"},
  ...
}
```

## Supported indices

| Name | Long name | Formula | Required bands |
|------|-----------|---------|---------------|
| `ndvi` | Normalised Difference Vegetation Index | (NIR − Red) / (NIR + Red) | B08, B04 |
| `evi` | Enhanced Vegetation Index | 2.5 · (NIR − Red) / (NIR + 6·Red − 7.5·Blue + 1) | B08, B04, B02 |
| `savi` | Soil Adjusted Vegetation Index | (NIR − Red) · 1.5 / (NIR + Red + 0.5) | B08, B04 |
| `ndwi` | Normalised Difference Water Index | (Green − NIR) / (Green + NIR) | B03, B08 |
| `mndwi` | Modified NDWI (Xu 2006) | (Green − SWIR1) / (Green + SWIR1) | B03, B11 |
| `ndbi` | Normalised Difference Built-up Index | (SWIR1 − NIR) / (SWIR1 + NIR) | B11, B08 |
| `nbr` | Normalised Burn Ratio | (NIR − SWIR2) / (NIR + SWIR2) | B08, B12 |
| `ndsi` | Normalised Difference Snow Index | (Green − SWIR1) / (Green + SWIR1) | B03, B11 |
| `cig` | Chlorophyll Index Green | (NIR / Green) − 1 | B08, B03 |
| `arvi` | Atmospherically Resistant Vegetation Index | (NIR − (2·Red − Blue)) / (NIR + (2·Red − Blue)) | B08, B04, B02 |

## Band name aliases

The module accepts any of the following naming conventions for the `band` coordinate — no manual remapping needed:

| Canonical | Accepted aliases |
|-----------|-----------------|
| B02 | `b02`, `blue` |
| B03 | `b03`, `green` |
| B04 | `b04`, `red` |
| B05 | `b05`, `rededge1`, `red_edge1` |
| B06 | `b06`, `rededge2`, `red_edge2` |
| B07 | `b07`, `rededge3`, `red_edge3` |
| B08 | `b08`, `nir`, `nir_broad` |
| B8A | `b8a`, `nir08`, `nir_narrow` |
| B11 | `b11`, `swir16`, `swir1` |
| B12 | `b12`, `swir22`, `swir2` |

Files produced by `downloader.py` (which uses human-readable names such as `blue`, `nir`, `swir16`) are read correctly without any extra configuration.

## Examples

### Minimal

```python
from sentinel_processor.indices.compute import compute_indices

results = compute_indices(
    source="downloads/spectral/budapest_20260520T094746.nc",
    indices=["ndvi", "ndwi"],
)
# {"ndvi": "downloads/spectral/indices/indices_budapest_20260520T094746_ndvi.tif",
#  "ndwi": "downloads/spectral/indices/indices_budapest_20260520T094746_ndwi.tif"}
```

### Custom output directory and scale factor

```python
from sentinel_processor.indices.compute import compute_indices

results = compute_indices(
    source="downloads/spectral/budapest_20260520T094746.nc",
    indices=["ndvi", "evi", "savi", "ndwi", "mndwi", "ndbi", "nbr"],
    output_dir="outputs/indices/budapest",
    scale_factor=1e-4,   # Sentinel-2 L2A integer DN → reflectance
)
```

### Save as NetCDF

```python
results = compute_indices(
    source="scene.nc",
    indices=["ndvi", "nbr"],
    output_format="nc",
)
```

### Skip existing files

```python
results = compute_indices(
    source="scene.nc",
    indices=["ndvi", "evi", "ndwi"],
    overwrite=False,   # re-runs only missing indices
)
```

### Integrate with downloader

```python
import sentinel_processor as sp
from sentinel_processor.indices.compute import compute_indices

# 1. download
download_results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=sp.DownloadConfig(bands=sp.SpectralBands.ALL),
)

# 2. compute indices for every downloaded scene
for base_name, paths in download_results.items():
    nc_files = [p for p in paths if p.endswith(".nc")]
    for nc in nc_files:
        compute_indices(
            source=nc,
            indices=["ndvi", "evi", "ndwi", "ndbi", "nbr"],
        )
```

### List available indices

```python
from sentinel_processor.indices.compute import list_indices

for name, info in list_indices().items():
    print(f"{name:8s}  needs {info['bands_required']}")
```

## Notes

- **CRS**: if the source file has no CRS (e.g. `.nc` files saved by older versions of the downloader), output `.tif` files are still written using the x/y coordinate extent as a spatial transform. Upgrade to the current downloader to preserve CRS end-to-end.
- **Memory**: all band arrays for the requested indices are loaded into RAM at once. For very large scenes consider processing a spatial subset first with `rioxarray.clip_box`.
- **Fortran kernel**: division-by-zero is handled inside Fortran — pixels where `|a + b| < 1e-9` are written as `0.0` rather than `NaN`, keeping output arrays fully finite.