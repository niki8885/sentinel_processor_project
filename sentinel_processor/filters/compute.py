from __future__ import annotations
import logging
from typing import Any

import numpy as np

from sentinel_processor.filters._filters_bridge import (
    apply_to_bands,
    bilateral_filter,
    convolve2d,
    gaussian_blur,
    laplacian,
    median_filter,
    morpho_close,
    morpho_dilate,
    morpho_erode,
    morpho_open,
    sobel_direction,
    sobel_magnitude,
    top_hat_black,
    top_hat_white,
    unsharp_mask,
)

logger = logging.getLogger(__name__)

_FILTER_REGISTRY: dict[str, dict] = {
    # Blur / smoothing
    "gaussian": {
        "fn": gaussian_blur,
        "params": {"sigma": 1.0, "kradius": None},
        "description": "Separable Gaussian blur",
    },
    "bilateral": {
        "fn": bilateral_filter,
        "params": {"sigma_s": 2.0, "sigma_r": 0.1, "kradius": None},
        "description": "Edge-preserving bilateral filter",
    },
    "median": {
        "fn": median_filter,
        "params": {"radius": 1},
        "description": "Median filter (box structuring element)",
    },

    # Edge detection
    "sobel_magnitude": {
        "fn": sobel_magnitude,
        "params": {"norm": "l2"},
        "description": "Sobel edge magnitude (L1 or L2 norm)",
    },
    "sobel_direction": {
        "fn": sobel_direction,
        "params": {},
        "description": "Sobel edge direction in radians [-π, π]",
    },
    "laplacian": {
        "fn": laplacian,
        "params": {"connectivity": 4},
        "description": "Discrete Laplacian (4- or 8-connectivity)",
    },

    # Sharpening
    "unsharp_mask": {
        "fn": unsharp_mask,
        "params": {"sigma": 1.0, "amount": 1.0, "kradius": None, "threshold": 0.0},
        "description": "Unsharp masking — high-boost sharpening",
    },

    # Morphology
    "erode": {
        "fn": morpho_erode,
        "params": {"radius": 1},
        "description": "Morphological erosion (rectangular SE)",
    },
    "dilate": {
        "fn": morpho_dilate,
        "params": {"radius": 1},
        "description": "Morphological dilation (rectangular SE)",
    },
    "open": {
        "fn": morpho_open,
        "params": {"radius": 1},
        "description": "Morphological opening — removes small bright blobs",
    },
    "close": {
        "fn": morpho_close,
        "params": {"radius": 1},
        "description": "Morphological closing — fills small dark holes",
    },
    "top_hat_white": {
        "fn": top_hat_white,
        "params": {"radius": 3},
        "description": "White top-hat: highlights small bright features",
    },
    "top_hat_black": {
        "fn": top_hat_black,
        "params": {"radius": 3},
        "description": "Black top-hat: highlights small dark features",
    },

    # Custom kernel
    "convolve": {
        "fn": convolve2d,
        "params": {"kernel": None},  # kernel is required
        "description": "Arbitrary 2-D convolution with a user-supplied kernel",
    },
}

AVAILABLE_FILTERS: list[str] = sorted(_FILTER_REGISTRY)


# Core API

def apply_filter(
        arr: np.ndarray,
        filter_name: str,
        **kwargs: Any,
) -> np.ndarray:
    entry = _FILTER_REGISTRY.get(filter_name.lower())
    if entry is None:
        raise ValueError(
            f"Unknown filter '{filter_name}'. "
            f"Available: {AVAILABLE_FILTERS}"
        )

    fn = entry["fn"]

    if arr.ndim == 2:
        return fn(arr.astype(np.float64, copy=False), **kwargs)
    if arr.ndim == 3:
        return apply_to_bands(fn, arr.astype(np.float64, copy=False), **kwargs)

    raise ValueError(f"Expected 2-D or 3-D array, got shape {arr.shape}")


def apply_filter_da(
        da,
        filter_name: str,
        **kwargs: Any,
):
    import xarray as xr

    arr = np.asarray(da.values, dtype=np.float64)
    result_np = apply_filter(arr, filter_name, **kwargs)

    result_da = xr.DataArray(
        result_np,
        dims=da.dims,
        coords=da.coords,
        attrs={**da.attrs, "filter_applied": filter_name, **kwargs},
    )

    try:
        crs = da.rio.crs
        if crs is not None:
            result_da = result_da.rio.write_crs(crs)
    except Exception:
        pass

    return result_da


def list_filters() -> dict[str, dict]:
    return {
        name: {
            "description": entry["description"],
            "default_params": entry["params"],
        }
        for name, entry in _FILTER_REGISTRY.items()
    }