from __future__ import annotations
import ctypes
import sys
from pathlib import Path

import numpy as np

_LIB_NAME = (
    "libsentinel_normalize.dll" if sys.platform == "win32"
    else "libsentinel_normalize.so"
)
_LIB_PATH = Path(__file__).parent / "fortran" / _LIB_NAME

if sys.platform == "win32":
    import os as _os

    for _candidate in [
        r"C:\msys64\ucrt64\bin",
        r"C:\msys64\mingw64\bin",
        r"C:\mingw64\bin",
    ]:
        if _os.path.isdir(_candidate):
            _os.add_dll_directory(_candidate)
            break

_lib = None
_DBL_P = ctypes.POINTER(ctypes.c_double)


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    if not _LIB_PATH.exists():
        raise FileNotFoundError(
            f"Fortran normalise library not found: {_LIB_PATH}\n"
            "Build (Linux/Mac):\n"
            f"  gfortran -O2 -shared -fPIC -o {_LIB_PATH} normalize_mod.f90\n"
            "Build (Windows):\n"
            f"  gfortran -O2 -shared -o {_LIB_PATH} normalize_mod.f90"
        )
    lib = ctypes.CDLL(str(_LIB_PATH))
    c_int = ctypes.c_int
    c_dbl = ctypes.c_double

    # minmax_band(arr, n, out_min, out_max, out)
    lib.minmax_band.restype = None
    lib.minmax_band.argtypes = [
        _DBL_P,  # arr(n)  in
        c_int,  # n       in  (by value)
        ctypes.POINTER(c_dbl),  # out_min out
        ctypes.POINTER(c_dbl),  # out_max out
        _DBL_P,  # out(n)  out
    ]

    # zscore_band(arr, n, mean, std, out)
    lib.zscore_band.restype = None
    lib.zscore_band.argtypes = [
        _DBL_P,  # arr(n)  in
        c_int,  # n       in  (by value)
        c_dbl,  # mean    in  (by value)
        c_dbl,  # std     in  (by value)
        _DBL_P,  # out(n)  out
    ]

    # band_stats(arr, n, out_mean, out_std, out_min, out_max)
    lib.band_stats.restype = None
    lib.band_stats.argtypes = [
        _DBL_P,  # arr(n)   in
        c_int,  # n        in  (by value)
        ctypes.POINTER(c_dbl),  # out_mean out
        ctypes.POINTER(c_dbl),  # out_std  out
        ctypes.POINTER(c_dbl),  # out_min  out
        ctypes.POINTER(c_dbl),  # out_max  out
    ]

    _lib = lib
    return _lib


# Array helpers

def _f64c(arr) -> np.ndarray:
    """Return a flat, C-contiguous float64 view/copy."""
    a = np.asarray(arr, dtype=np.float64).ravel()
    if not a.flags["C_CONTIGUOUS"]:
        a = np.ascontiguousarray(a)
    return a


def _ptr(arr: np.ndarray) -> ctypes.POINTER:
    return arr.ctypes.data_as(_DBL_P)


def _alloc(n: int) -> tuple[np.ndarray, ctypes.POINTER]:
    buf = np.empty(n, dtype=np.float64)
    return buf, buf.ctypes.data_as(_DBL_P)


# Public API

def minmax_band(arr) -> tuple[np.ndarray, float, float]:
    """Min–max stretch to [0, 1].

    Returns
    -------
    (result, vmin, vmax)
    """
    a = _f64c(arr)
    n = a.size
    out, out_ptr = _alloc(n)
    vmin = ctypes.c_double(0.0)
    vmax = ctypes.c_double(0.0)
    _get_lib().minmax_band(
        _ptr(a), ctypes.c_int(n),
        ctypes.byref(vmin), ctypes.byref(vmax),
        out_ptr,
    )
    return out, vmin.value, vmax.value


def zscore_band(arr, mean: float, std: float) -> np.ndarray:
    """Z-score standardisation: (arr - mean) / std."""
    a = _f64c(arr)
    n = a.size
    out, out_ptr = _alloc(n)
    _get_lib().zscore_band(
        _ptr(a), ctypes.c_int(n),
        ctypes.c_double(mean), ctypes.c_double(std),
        out_ptr,
    )
    return out


def band_stats(arr) -> dict[str, float]:
    """Compute mean, std, min, max of a band (NODATA excluded)."""
    a = _f64c(arr)
    n = a.size
    mean = ctypes.c_double(0.0)
    std = ctypes.c_double(0.0)
    vmin = ctypes.c_double(0.0)
    vmax = ctypes.c_double(0.0)
    _get_lib().band_stats(
        _ptr(a), ctypes.c_int(n),
        ctypes.byref(mean), ctypes.byref(std),
        ctypes.byref(vmin), ctypes.byref(vmax),
    )
    return {
        "mean": mean.value,
        "std": std.value,
        "min": vmin.value,
        "max": vmax.value,
    }
