# sentinel_processor.dl — DL Preprocessing Reference

Pre-processing utilities for feeding Sentinel-2 imagery into deep-learning
models.  The module covers three concerns:

1. **Normalisation** — per-band scaling to the range a pre-trained encoder expects.
2. **Dataset statistics** — computing and serialising global mean/std across a
   collection of scenes for z-score normalisation.
3. **Tiling** — cutting full scenes into overlapping patches and stitching
   model outputs back into a seamless mosaic.

---

## Quick start — 4-band RGBN tensor for a torchvision model

```python
import numpy as np
import torch
from sentinel_processor.dl import normalize_for_dl, extract_tiles, stitch_tiles

# ── 1.  Load your scene (C=4, H, W) in DN units [0, 10000] ──────────────────
# e.g. from rioxarray, gdal, or your own loader:
scene = np.random.randint(0, 10000, (4, 1024, 1024)).astype(np.float64)
# band order must be: B04 (red), B03 (green), B02 (blue), B08 (NIR)
bands = ["B04", "B03", "B02", "B08"]

# ── 2.  Normalise with the SSL4EO-S12 RGBN preset ───────────────────────────
norm = normalize_for_dl(scene, method="sentinel2_rgbn", bands=bands)
# norm.shape == (4, 1024, 1024), float64

# ── 3.  Convert to torch tensor (C, H, W) ───────────────────────────────────
tensor = torch.from_numpy(norm).float()   # (4, 1024, 1024)

# ── 4.  (Optional) Cut into 256×256 tiles with 32 px overlap ────────────────
tiles, meta = extract_tiles(norm, tile_size=256, overlap=32)
# tiles.shape == (N, 4, 256, 256)

# ── 5.  Run model on each tile ───────────────────────────────────────────────
predictions = np.stack([
    my_model(torch.from_numpy(t[None]).float())[0].numpy()
    for t in tiles
])  # (N, num_classes, 256, 256)

# ── 6.  Stitch predictions back ──────────────────────────────────────────────
output_map = stitch_tiles(predictions, meta)   # (num_classes, 1024, 1024)
```

---

## Normalisation

### `normalize_for_dl(arr, method, bands=None, stats=None, nodata=-9999.0)`

| Argument | Type | Description |
|---|---|---|
| `arr` | `np.ndarray` | `(C, H, W)`, `(C, N)`, or `(C,)` float64 array |
| `method` | `str` | See table below |
| `bands` | `list[str]` | Band IDs matching axis-0 order; required for preset/keyed-zscore methods |
| `stats` | `dict` | Band-keyed stats dict; required for `zscore` when not computing from scene |
| `nodata` | `float` | Sentinel value (default −9999); preserved in output |

Returns `np.ndarray` of the same shape, float64.

#### Methods

| Method | Description | `bands` required | `stats` required |
|---|---|---|---|
| `minmax` | Per-band linear stretch to **[0, 1]** using scene min/max | No | No |
| `zscore` | `(x − mean) / std`; uses provided `stats` or computes scene-level | Optional | Optional |
| `sentinel2_rgb` | SeCo mean/std for **B04, B03, B02** | Yes | No |
| `sentinel2_rgbn` | SSL4EO-S12 mean/std for **B04, B03, B02, B08** | Yes | No |
| `sentinel2_all` | SSL4EO-S12 mean/std for all **10 bands** (B02–B12) | Yes | No |

#### Preset sources

| Preset | Paper | Bands | DN scale |
|---|---|---|---|
| `sentinel2_rgb` | Mañas et al. (2021) *SeCo* | B04, B03, B02 | ×10 000 (L1C TOA) |
| `sentinel2_rgbn` | Wang et al. (2022) *SSL4EO-S12* | B04, B03, B02, B08 | ×10 000 (L1C TOA) |
| `sentinel2_all` | Wang et al. (2022) *SSL4EO-S12* | 10 MSI bands | ×10 000 (L1C TOA) |

> **Band order for `sentinel2_all`:** B02, B03, B04, B08, B05, B06, B07, B8A, B11, B12.
> Always pass a matching `bands` list so the correct statistics are applied.

---

## Dataset statistics

### `compute_dataset_stats(scene_paths, bands) → dict`

Computes global mean and std across a list of scene files.  Per-scene
statistics are computed in Fortran (`band_stats`) using the Welford online
algorithm and then combined in Python using the parallel Welford formula, so
no scene needs to be held in memory simultaneously.

```python
from sentinel_processor.dl import compute_dataset_stats, save_stats, load_stats

stats = compute_dataset_stats(
    scene_paths=["scene_001.tif", "scene_002.tif", "scene_003.tif"],
    bands=["B04", "B03", "B02", "B08"],
)
# {"B04": {"mean": ..., "std": ..., "min": ..., "max": ..., "n_pixels": ...}, ...}

save_stats(stats, "my_dataset_stats.json")

# Later:
stats = load_stats("my_dataset_stats.json")
norm = normalize_for_dl(scene, method="zscore", bands=["B04","B03","B02","B08"],
                        stats=stats)
```

The returned dict (and JSON file) follow the schema:
```json
{
  "B04": { "mean": 1354.99, "std": 2173.26, "min": 0.0, "max": 10000.0, "n_pixels": 1234567 }
}
```

### `save_stats(stats, path)` / `load_stats(path) → dict`

JSON round-trip helpers.  `.json` extension is appended automatically if
absent.  Both functions log via the standard `logging` module.

---

## Tiling

### `extract_tiles(arr, tile_size=256, overlap=32, nodata=-9999.0) → (tiles, meta)`

Cuts a `(C, H, W)` array into `(N, C, tile_size, tile_size)` overlapping
patches.  The image is reflect-padded on the bottom/right if needed so the
grid covers the full extent.

| Argument | Default | Description |
|---|---|---|
| `tile_size` | 256 | Square tile edge in pixels |
| `overlap` | 32 | Overlap between adjacent tiles (must be < `tile_size / 2`) |

Returns `(tiles, meta)` where `meta` is a `TileMeta` NamedTuple containing
the full geometry (original shape, padded shape, grid dimensions, stride,
padding amounts).

### `stitch_tiles(tiles, meta, blend=True) → np.ndarray`

Reconstructs the original extent from model-output tiles.  With
`blend=True` (default), overlapping regions are averaged using a cosine
taper window that gives full weight at the tile centre and fades to near-zero
at each edge, eliminating block-boundary artefacts.

```python
tiles, meta = extract_tiles(norm, tile_size=256, overlap=32)
# ... run inference on tiles ...
output = stitch_tiles(predictions, meta, blend=True)
# output.shape == norm.shape  (padding is cropped automatically)
```

### `tile_grid_shape(image_h, image_w, tile_size, overlap) → (n_tiles, rows, cols)`

Computes the tiling geometry without allocating any arrays — useful for
pre-allocating output buffers before running inference.

---

## Fortran internals

The hot path lives in `fortran/normalize_mod.f90` (compiled to
`libsentinel_normalize.so` / `.dll`).

| Fortran symbol | Purpose |
|---|---|
| `minmax_band(arr, n, out_min, out_max, out)` | Linear stretch to [0,1] |
| `zscore_band(arr, n, mean, std, out)` | Standardise with fixed mean/std |
| `band_stats(arr, n, out_mean, out_std, out_min, out_max)` | One-pass Welford stats |

All routines skip NODATA pixels (sentinel −9999.0) and are safe with
all-NODATA inputs.

### Build commands

```bash
# Linux / macOS
gfortran -O2 -shared -fPIC \
  -o sentinel_processor/dl/fortran/libsentinel_normalize.so \
  sentinel_processor/dl/fortran/normalize_mod.f90

# Windows (MSYS2 / MinGW)
gfortran -O2 -shared \
  -o sentinel_processor/dl/fortran/libsentinel_normalize.dll \
  sentinel_processor/dl/fortran/normalize_mod.f90
```

---

## API surface — `from sentinel_processor.dl import …`

```
normalize_for_dl        # main normalisation entry point
compute_dataset_stats   # compute global statistics across scenes
save_stats              # write stats dict to JSON
load_stats              # read stats dict from JSON
AVAILABLE_METHODS       # list of valid method strings

extract_tiles           # cut image into overlapping patches
stitch_tiles            # reconstruct from patches with cosine blending
tile_grid_shape         # compute grid geometry without allocation
TileMeta                # NamedTuple describing tile geometry

get_preset              # look up a named statistics preset
preset_bands            # ordered band list for a preset
AVAILABLE_PRESETS       # list of preset names
SSL4EO_S12_ALL          # raw stats dict for all 10 bands
SECO_RGB                # raw stats dict for RGB (SeCo)
SSL4EO_RGBN             # raw stats dict for RGBN
```