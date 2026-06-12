from __future__ import annotations
import ctypes
import sys
from pathlib import Path
from typing import Literal

import numpy as np

_LIB_NAME = (
    "libsentinel_texture.dll" if sys.platform == "win32"
    else "libsentinel_texture.so"
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

NODATA: float = -9999.0
Angle = Literal[0, 45, 90, 135, -1]


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    if not _LIB_PATH.exists():
        raise FileNotFoundError(
            f"Fortran shared library not found: {_LIB_PATH}\n"
            "Build with:\n"
            "  gfortran -O2 -shared -fPIC \\\n"
            "    -o sentinel_processor/texture/fortran/libsentinel_texture.so \\\n"
            "    sentinel_processor/texture/fortran/texture_mod.f90"
        )
    _lib = ctypes.CDLL(str(_LIB_PATH))
    c_int = ctypes.c_int

    _lib.compute_glcm.restype = None
    _lib.compute_glcm.argtypes = [
        _DBL_P,     # arr          (rows, cols) Fortran order
        c_int,      # rows
        c_int,      # cols
        c_int,      # window
        c_int,      # distance
        c_int,      # angle_deg  (-1 = isotropic)
        _DBL_P,     # energy_out
        _DBL_P,     # contrast_out
        _DBL_P,     # homogeneity_out
    ]
    return _lib


def _f64_fortran(arr: np.ndarray) -> tuple[ctypes.POINTER, np.ndarray]:
    """Return (pointer, Fortran-contiguous float64 array)."""
    a = np.asarray(arr, dtype=np.float64)
    if not a.flags["F_CONTIGUOUS"]:
        a = np.asfortranarray(a)
    return a.ctypes.data_as(_DBL_P), a


def _empty_fortran(shape: tuple) -> tuple[ctypes.POINTER, np.ndarray]:
    a = np.empty(shape, dtype=np.float64, order="F")
    return a.ctypes.data_as(_DBL_P), a


# Public API


def compute_glcm(
        arr: np.ndarray,
        window: int = 7,
        distance: int = 1,
        angle: int = -1,
) -> dict[str, np.ndarray]:
    if arr.ndim != 2:
        raise ValueError(f"arr must be 2-D (rows, cols), got shape {arr.shape}")

    valid_angles = {-1, 0, 45, 90, 135}
    if angle not in valid_angles:
        raise ValueError(
            f"angle must be one of {valid_angles}, got {angle!r}"
        )

    if window < 3:
        window = 3
    if window % 2 == 0:
        window += 1

    rows, cols = arr.shape

    src_ptr, src_ref = _f64_fortran(arr)
    en_ptr,  en_out  = _empty_fortran((rows, cols))
    co_ptr,  co_out  = _empty_fortran((rows, cols))
    ho_ptr,  ho_out  = _empty_fortran((rows, cols))

    _get_lib().compute_glcm(
        src_ptr,
        ctypes.c_int(rows),
        ctypes.c_int(cols),
        ctypes.c_int(window),
        ctypes.c_int(distance),
        ctypes.c_int(angle),
        en_ptr,
        co_ptr,
        ho_ptr,
    )

    del src_ref

    return {
        "energy":       np.ascontiguousarray(en_out),
        "contrast":     np.ascontiguousarray(co_out),
        "homogeneity":  np.ascontiguousarray(ho_out),
    }