from __future__ import annotations
import json
import logging
import math
from pathlib import Path
from typing import Sequence, Literal
import numpy as np
from sentinel_processor.dl._normalize_bridge import (
    minmax_band,
    zscore_band,
    band_stats as _band_stats_fortran,
)
from sentinel_processor.dl.presets import (
    AVAILABLE_PRESETS,
    get_preset,
)

logger = logging.getLogger(__name__)

NODATA: float = -9999.0

NormMethod = Literal[
    "minmax",
    "zscore",
    "sentinel2_rgb",
    "sentinel2_rgbn",
    "sentinel2_all",
]

AVAILABLE_METHODS: list[str] = ["minmax", "zscore"] + AVAILABLE_PRESETS


def normalize_for_dl(
        arr: np.ndarray,
        method: NormMethod,
        bands: Sequence[str] | None = None,
        stats: dict | None = None,
        nodata: float = NODATA,
) -> np.ndarray:
    """Normalise a (bands, H, W) or (bands,) array for DL inference.

    Parameters
    ----------
    arr : np.ndarray
        Input array.  Shape: ``(C, H, W)`` for multi-band spatial data,
        ``(C,)`` for a pixel vector, or ``(C, N)`` for a batch of pixels.
    method : str
        One of:

        ``'minmax'``
            Per-band linear stretch to [0, 1] using scene min/max.

        ``'zscore'``
            Per-band ``(x - mean) / std``.  Provide ``stats`` (output of
            :func:`compute_dataset_stats` or a custom dict) **or** set
            ``stats=None`` to compute scene-level statistics on the fly.

        ``'sentinel2_rgb'``
            Fixed SeCo mean/std for B04, B03, B02 (in that order).

        ``'sentinel2_rgbn'``
            SSL4EO-S12 mean/std for B04, B03, B02, B08.

        ``'sentinel2_all'``
            SSL4EO-S12 mean/std for all 10 bands: B02, B03, B04, B08,
            B05, B06, B07, B8A, B11, B12.

    bands : list[str] or None
        Band identifiers in the same order as axis-0 of ``arr``.
        Required for preset methods (``sentinel2_*``) and for ``zscore``
        when ``stats`` is a band-keyed dict.  Not needed for ``minmax``
        when ``stats=None``.

    stats : dict or None
        For ``zscore``: a dict returned by :func:`compute_dataset_stats`
        or manually constructed as ``{"B04": {"mean": ..., "std": ...}, ...}``.
        Pass ``None`` to compute scene statistics on the fly.

    nodata : float
        Pixels equal to this value are masked and preserved in the output.

    Returns
    -------
    np.ndarray  – float64, same shape as input.
    """
    method = method.strip().lower()
    if method not in AVAILABLE_METHODS:
        raise ValueError(
            f"Unknown normalisation method '{method}'. "
            f"Available: {AVAILABLE_METHODS}"
        )

    arr = np.asarray(arr, dtype=np.float64)
    original_shape = arr.shape

    # The Fortran kernels recognise only NODATA (-9999); remap a custom
    # sentinel onto it before the call and restore it in the output.
    nodata_mask: np.ndarray | None = None
    if nodata != NODATA:
        nodata_mask = np.isclose(arr, nodata, atol=1e-3)
        arr = np.where(nodata_mask, NODATA, arr)

    if arr.ndim == 1:
        arr_2d = arr[:, np.newaxis]
    elif arr.ndim == 2:
        arr_2d = arr
    elif arr.ndim == 3:
        # (C, H, W) → (C, H*W)
        C, H, W = arr.shape
        arr_2d = arr.reshape(C, H * W)
    else:
        raise ValueError(f"arr must be 1-D, 2-D, or 3-D, got shape {arr.shape}")

    C, N = arr_2d.shape
    result = np.empty_like(arr_2d)

    if method == "minmax":
        for b in range(C):
            norm, _, _ = minmax_band(arr_2d[b])
            result[b] = norm

    elif method == "zscore":
        if stats is None:
            # Compute scene-level statistics on the fly
            for b in range(C):
                s = _band_stats_fortran(arr_2d[b])
                result[b] = zscore_band(arr_2d[b], s["mean"], s["std"])
        else:
            # stats keyed by band name or by integer index
            for b in range(C):
                key = bands[b] if bands is not None else str(b)
                if key not in stats:
                    raise KeyError(
                        f"Band '{key}' not found in stats dict. "
                        f"Available keys: {list(stats)}"
                    )
                bstats = stats[key]
                result[b] = zscore_band(
                    arr_2d[b],
                    float(bstats["mean"]),
                    float(bstats["std"]),
                )

    else:
        # Preset method: sentinel2_rgb / sentinel2_rgbn / sentinel2_all
        preset = get_preset(method)
        if bands is None:
            raise ValueError(
                f"'bands' must be provided for preset method '{method}'. "
                f"Expected bands: {list(preset.keys())}"
            )
        if len(bands) != C:
            raise ValueError(
                f"len(bands)={len(bands)} does not match arr.shape[0]={C}."
            )
        for b, band_key in enumerate(bands):
            if band_key not in preset:
                raise KeyError(
                    f"Band '{band_key}' not in preset '{method}'. "
                    f"Preset covers: {list(preset.keys())}"
                )
            bstats = preset[band_key]
            result[b] = zscore_band(
                arr_2d[b],
                float(bstats["mean"]),
                float(bstats["std"]),
            )

    result = result.reshape(original_shape)
    if nodata_mask is not None:
        result[nodata_mask] = nodata
    return result


# Dataset statistics

def compute_dataset_stats(
        scene_paths: Sequence[str | Path],
        bands: Sequence[str],
        nodata: float = NODATA,
) -> dict[str, dict[str, float]]:
    """Compute global mean and std across a list of scenes.

    Uses Welford's online algorithm to aggregate per-scene statistics
    without loading all scenes into memory simultaneously.  Per-scene
    per-band stats are computed in Fortran (``band_stats``); Python
    aggregates them using the parallel Welford combine formula.

    Parameters
    ----------
    scene_paths : list of str or Path
        Paths to GeoTIFF or NetCDF files.  Each file must expose the
        requested ``bands`` (compatible with the conventions used in
        ``sentinel_processor.indices.compute``).
    bands : list[str]
        Band identifiers to compute statistics for, e.g.
        ``['B04', 'B03', 'B02', 'B08']``.

    Returns
    -------
    dict  –  ``{"B04": {"mean": ..., "std": ..., "min": ..., "max": ...}, ...}``
    """
    try:
        import rioxarray  # noqa: F401  — availability check
        import xarray as xr  # noqa: F401  — availability check
    except ImportError as exc:
        raise ImportError(
            "rioxarray and xarray are required for compute_dataset_stats. "
            "Install them with: pip install rioxarray xarray"
        ) from exc

    bands = list(bands)

    agg: dict[str, dict] = {
        b: {"n": 0, "mean": 0.0, "M2": 0.0,
            "min": float("inf"), "max": float("-inf")}
        for b in bands
    }

    for scene_idx, scene_path in enumerate(scene_paths):
        scene_path = Path(scene_path)
        logger.info(
            f"[dl/stats] Scene {scene_idx + 1}/{len(scene_paths)}: {scene_path.name}"
        )
        try:
            band_arrays = _load_bands_for_stats(scene_path, bands, nodata)
        except Exception as exc:
            logger.warning(f"[dl/stats] Skipping {scene_path.name}: {exc}")
            continue

        for b in bands:
            if b not in band_arrays:
                logger.warning(
                    f"[dl/stats] Band {b} missing in {scene_path.name}, skipping band."
                )
                continue
            flat = band_arrays[b].ravel().astype(np.float64)
            if nodata != NODATA:
                # Fortran band_stats recognises only NODATA (-9999)
                flat = np.where(np.isclose(flat, nodata, atol=1e-3), NODATA, flat)
            s = _band_stats_fortran(flat)
            n_b = int(np.sum(
                ~np.isclose(flat, NODATA, atol=1e-3)
            ))
            if n_b == 0:
                continue
            _welford_combine(agg[b], n_b, s["mean"], s["std"], s["min"], s["max"])

    result: dict[str, dict[str, float]] = {}
    for b in bands:
        a = agg[b]
        if a["n"] == 0:
            logger.warning(f"[dl/stats] No valid pixels found for band {b}.")
            result[b] = {"mean": 0.0, "std": 1.0, "min": 0.0, "max": 0.0,
                         "n_pixels": 0}
        else:
            result[b] = {
                "mean": a["mean"],
                "std": math.sqrt(max(a["M2"] / a["n"], 0.0)),
                "min": a["min"],
                "max": a["max"],
                "n_pixels": a["n"],
            }
    return result


def _welford_combine(
        agg: dict,
        n_b: int,
        mean_b: float,
        std_b: float,
        min_b: float,
        max_b: float,
) -> None:
    """Merge a new batch (n_b, mean_b, var_b) into an existing Welford accumulator."""
    n_a = agg["n"]
    if n_a == 0:
        agg["n"] = n_b
        agg["mean"] = mean_b
        agg["M2"] = std_b * std_b * n_b
        agg["min"] = min_b
        agg["max"] = max_b
        return

    n_new = n_a + n_b
    delta = mean_b - agg["mean"]
    new_mean = agg["mean"] + delta * n_b / n_new
    new_M2 = (agg["M2"]
              + std_b * std_b * n_b
              + delta * delta * n_a * n_b / n_new)
    agg["n"] = n_new
    agg["mean"] = new_mean
    agg["M2"] = new_M2
    agg["min"] = min(agg["min"], min_b)
    agg["max"] = max(agg["max"], max_b)


def _load_bands_for_stats(
        path: Path,
        bands: list[str],
        nodata: float,
) -> dict[str, np.ndarray]:
    import rioxarray
    import xarray as xr
    from sentinel_processor.indices.compute import (
        _load_bands_from_multiband,
        _load_bands_from_dataset,
    )

    suffix = path.suffix.lower()
    needed = set(bands)

    if suffix in (".tif", ".tiff"):
        da = rioxarray.open_rasterio(path)
        result = _load_bands_from_multiband(da, needed)
        if result is None:
            raise ValueError(f"No recognised band coordinate in {path.name}")
        return result

    elif suffix == ".nc":
        ds = xr.open_dataset(path)
        result = _load_bands_from_dataset(ds, needed)
        if result is not None:
            ds.close()
            return result
        for var in ds.data_vars:
            candidate = ds[var]
            if "band" in (candidate.coords or {}):
                result = _load_bands_from_multiband(candidate, needed)
                if result is not None:
                    ds.close()
                    return result
        ds.close()
        raise ValueError(f"Could not extract bands from NetCDF {path.name}")

    else:
        raise ValueError(f"Unsupported file format '{suffix}'")


# Stats persistence

def save_stats(stats: dict, path: str | Path) -> None:
    """Serialise a stats dict to JSON.

    Parameters
    ----------
    stats : dict
        Output of :func:`compute_dataset_stats` or a custom compatible dict.
    path : str or Path
        Destination file path.  The ``.json`` extension is appended if absent.
    """
    path = Path(path)
    if path.suffix.lower() != ".json":
        path = path.with_suffix(".json")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2)
    logger.info(f"[dl/stats] Stats saved → {path}")


def load_stats(path: str | Path) -> dict:
    """Load a stats dict previously saved with :func:`save_stats`.

    Parameters
    ----------
    path : str or Path
        Path to the JSON stats file.

    Returns
    -------
    dict  –  ``{"B04": {"mean": ..., "std": ...}, ...}``
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Stats file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        stats = json.load(fh)
    logger.info(f"[dl/stats] Stats loaded ← {path}")
    return stats
