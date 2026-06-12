from __future__ import annotations
import ctypes
import sys
from pathlib import Path
from typing import Literal
import numpy as np

_LIB_NAME = (
    "libsentinel_timeseries.dll" if sys.platform == "win32"
    else "libsentinel_timeseries.so"
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
_INT_P = ctypes.POINTER(ctypes.c_int)

_METHOD_MAP: dict[str, int] = {
    "linear": 0,
    "savgol": 1,
    "pchip": 2,
    "ets": 3,
    "gauss": 4,
}

Method = Literal["linear", "savgol", "pchip", "ets", "gauss"]
NODATA: float = -9999.0


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    if not _LIB_PATH.exists():
        raise FileNotFoundError(
            f"Fortran shared library not found: {_LIB_PATH}\n"
            "Build with:\n"
            "  gfortran -O2 -shared -fPIC \\\n"
            "    -o sentinel_processor/processing/fortran/libsentinel_timeseries.so \\\n"
            "    sentinel_processor/processing/fortran/timeseries_mod.f90"
        )
    _lib = ctypes.CDLL(str(_LIB_PATH))
    c_int = ctypes.c_int
    _lib.interpolate_gaps.restype = None
    _lib.interpolate_gaps.argtypes = [
        _DBL_P,  # arr
        _INT_P,  # mask
        c_int,  # n_times
        c_int,  # rows
        c_int,  # cols
        c_int,  # method
        c_int,  # window
    ]
    return _lib


def _prepare_arr(arr: np.ndarray) -> tuple[ctypes.POINTER, np.ndarray]:
    """Return (pointer, C-contiguous float64 copy)."""
    a = np.array(arr, dtype=np.float64, order="C", copy=True)
    return a.ctypes.data_as(_DBL_P), a


def _prepare_mask(mask: np.ndarray, shape: tuple) -> tuple[ctypes.POINTER, np.ndarray]:
    m = np.asarray(mask, dtype=np.int32)
    if m.shape != shape:
        raise ValueError(f"mask shape {m.shape} != arr shape {shape}")
    m = np.ascontiguousarray(m)
    return m.ctypes.data_as(_INT_P), m


# Public API

def interpolate_gaps(
        arr: np.ndarray,
        mask: np.ndarray,
        method: Method = "linear",
        window: int = 5,
) -> np.ndarray:
    """Fill cloud-masked gaps in a (n_times, rows, cols) time stack.

    Parameters
    ----------
    arr    : (n_times, rows, cols) float64 array
    mask   : (n_times, rows, cols) – non-zero = valid, 0 = gap
    method : 'linear' | 'savgol' | 'pchip' | 'ets' | 'gauss'
    window : S-G / Gaussian window (odd, >= 3). Forced odd if even.

    Returns
    -------
    np.ndarray  (n_times, rows, cols) float64

    Method notes
    ------------
    linear  – Piecewise linear between nearest valid neighbours.
              Fast, no parameters, good baseline.
    savgol  – Savitzky-Golay quadratic smoother. Preserves peaks and troughs
              better than a simple moving average. Use window=5..11.
    pchip   – Monotone cubic Hermite spline (Fritsch-Carlson). Smooth, no
              overshoot, physically plausible. Best for NDVI/EVI/LAI.
    ets     – Holt double-exponential (level + trend). Alpha and beta chosen
              automatically. Good for series with persistent seasonal drift.
    gauss   – Gaussian-weighted moving average. Gentle smoothing, symmetric
              kernel, renormalised at edges. Use window=7..15.
    """
    if arr.ndim != 3:
        raise ValueError(f"arr must be 3-D (n_times, rows, cols), got {arr.shape}")

    if method not in _METHOD_MAP:
        raise ValueError(
            f"Unknown method {method!r}. Choose from: {list(_METHOD_MAP)}"
        )

    n_times, rows, cols = arr.shape
    method_flag = _METHOD_MAP[method]

    if method in ("savgol", "gauss"):
        if window < 3:
            window = 3
        if window % 2 == 0:
            window += 1

    arr_ptr, arr_buf = _prepare_arr(arr)
    mask_ptr, mask_buf = _prepare_mask(mask, arr.shape)

    _get_lib().interpolate_gaps(
        arr_ptr, mask_ptr,
        ctypes.c_int(n_times),
        ctypes.c_int(rows),
        ctypes.c_int(cols),
        ctypes.c_int(method_flag),
        ctypes.c_int(window),
    )

    return arr_buf
