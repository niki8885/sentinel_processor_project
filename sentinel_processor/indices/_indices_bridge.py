from __future__ import annotations
import ctypes
import sys
from pathlib import Path

import numpy as np

_LIB_NAME = (
    "libsentinel_indices.dll" if sys.platform == "win32"
    else "libsentinel_indices.so"
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
_PTR_DBL = ctypes.POINTER(ctypes.c_double)


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    if not _LIB_PATH.exists():
        raise FileNotFoundError(
            f"Fortran indices library not found: {_LIB_PATH}\n"
            f"Build (Linux/Mac):\n"
            f"  gfortran -O2 -shared -fPIC -o {_LIB_PATH} indices_mod.f90\n"
            f"Build (Windows):\n"
            f"  gfortran -O2 -shared -o {_LIB_PATH} indices_mod.f90"
        )
    lib = ctypes.CDLL(str(_LIB_PATH))

    _two_band = [
        "compute_ndvi", "compute_savi", "compute_ndwi", "compute_mndwi",
        "compute_ndbi", "compute_nbr", "compute_ndsi", "compute_cig",
    ]
    _three_band = ["compute_evi", "compute_arvi"]

    for name in _two_band:
        fn = getattr(lib, name)
        fn.restype = None
        fn.argtypes = [_PTR_DBL, _PTR_DBL, ctypes.c_int, _PTR_DBL]

    for name in _three_band:
        fn = getattr(lib, name)
        fn.restype = None
        fn.argtypes = [_PTR_DBL, _PTR_DBL, _PTR_DBL, ctypes.c_int, _PTR_DBL]

    _lib = lib
    return _lib


# Array helpers

def _ensure_f64c(arr) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    a = a.ravel()  # view if C-contiguous, else copy
    if not a.flags["C_CONTIGUOUS"]:
        a = np.ascontiguousarray(a)
    return a


def _prepare(*arrays) -> tuple[list[np.ndarray], int]:
    out = [_ensure_f64c(a) for a in arrays]
    sizes = {a.size for a in out}
    if len(sizes) != 1:
        raise ValueError(
            f"All band arrays must have the same number of elements; "
            f"got sizes {sizes}"
        )
    return out, out[0].size


def _c(arr: np.ndarray) -> ctypes.POINTER:
    return arr.ctypes.data_as(_PTR_DBL)


def _alloc(n: int) -> tuple[np.ndarray, ctypes.POINTER]:
    buf = np.empty(n, dtype=np.float64)
    return buf, buf.ctypes.data_as(_PTR_DBL)


# Public API

def compute_ndvi(nir, red) -> np.ndarray:
    (a, b), n = _prepare(nir, red)
    out, ptr = _alloc(n)
    _get_lib().compute_ndvi(_c(a), _c(b), ctypes.c_int(n), ptr)
    return out


def compute_evi(nir, red, blue) -> np.ndarray:
    (a, b, c), n = _prepare(nir, red, blue)
    out, ptr = _alloc(n)
    _get_lib().compute_evi(_c(a), _c(b), _c(c), ctypes.c_int(n), ptr)
    return out


def compute_savi(nir, red) -> np.ndarray:
    (a, b), n = _prepare(nir, red)
    out, ptr = _alloc(n)
    _get_lib().compute_savi(_c(a), _c(b), ctypes.c_int(n), ptr)
    return out


def compute_ndwi(green, nir) -> np.ndarray:
    (a, b), n = _prepare(green, nir)
    out, ptr = _alloc(n)
    _get_lib().compute_ndwi(_c(a), _c(b), ctypes.c_int(n), ptr)
    return out


def compute_mndwi(green, swir1) -> np.ndarray:
    (a, b), n = _prepare(green, swir1)
    out, ptr = _alloc(n)
    _get_lib().compute_mndwi(_c(a), _c(b), ctypes.c_int(n), ptr)
    return out


def compute_ndbi(swir1, nir) -> np.ndarray:
    (a, b), n = _prepare(swir1, nir)
    out, ptr = _alloc(n)
    _get_lib().compute_ndbi(_c(a), _c(b), ctypes.c_int(n), ptr)
    return out


def compute_nbr(nir, swir2) -> np.ndarray:
    (a, b), n = _prepare(nir, swir2)
    out, ptr = _alloc(n)
    _get_lib().compute_nbr(_c(a), _c(b), ctypes.c_int(n), ptr)
    return out


def compute_ndsi(green, swir1) -> np.ndarray:
    (a, b), n = _prepare(green, swir1)
    out, ptr = _alloc(n)
    _get_lib().compute_ndsi(_c(a), _c(b), ctypes.c_int(n), ptr)
    return out


def compute_cig(nir, green) -> np.ndarray:
    (a, b), n = _prepare(nir, green)
    out, ptr = _alloc(n)
    _get_lib().compute_cig(_c(a), _c(b), ctypes.c_int(n), ptr)
    return out


def compute_arvi(nir, red, blue) -> np.ndarray:
    (a, b, c), n = _prepare(nir, red, blue)
    out, ptr = _alloc(n)
    _get_lib().compute_arvi(_c(a), _c(b), _c(c), ctypes.c_int(n), ptr)
    return out
