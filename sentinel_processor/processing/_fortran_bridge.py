from __future__ import annotations
import ctypes
import sys
from pathlib import Path
from typing import Literal
import numpy as np

_LIB_NAME = (
    "libsentinel_processing.dll" if sys.platform == "win32"
    else "libsentinel_processing.so"
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
            f"Fortran shared library not found: {_LIB_PATH}\n"
            "Build with the provided Makefile:\n"
            "  make -C sentinel_processor/processing/fortran\n"
            "or manually:\n"
            "  gfortran -O2 -shared -fPIC -o libsentinel_processing.so pansharpening.f90\n"
            "  (Windows) gfortran -O2 -shared -o libsentinel_processing.dll pansharpening.f90"
        )
    _lib = ctypes.CDLL(str(_LIB_PATH))

    c_int = ctypes.c_int

    _lib.gram_schmidt_sharpen.restype = None
    _lib.gram_schmidt_sharpen.argtypes = [
        _DBL_P, c_int, c_int,
        _DBL_P, c_int, c_int,
        c_int,
        _DBL_P,
    ]
    _lib.ihs_sharpen.restype = None
    _lib.ihs_sharpen.argtypes = [
        _DBL_P, c_int, c_int,
        _DBL_P, c_int, c_int,
        _DBL_P,
    ]
    _lib.wavelet_sharpen.restype = None
    _lib.wavelet_sharpen.argtypes = [
        _DBL_P, c_int, c_int,
        _DBL_P, c_int, c_int,
        c_int,
        _DBL_P,
    ]
    return _lib


# Array helpers

def _as_fortran_double(arr: np.ndarray) -> tuple[ctypes.POINTER, np.ndarray]:
    """Return (pointer, array) — copy only when dtype or order differ."""
    a = np.asarray(arr, dtype=np.float64)
    if not a.flags["F_CONTIGUOUS"]:
        a = np.asfortranarray(a)
    return a.ctypes.data_as(_DBL_P), a


def _ms_to_fortran(
        ms: np.ndarray,
) -> tuple[ctypes.POINTER, np.ndarray, int, int, int]:
    if ms.ndim != 3:
        raise ValueError(
            f"ms must be 3-D (n_bands, rows, cols), got shape {ms.shape}"
        )
    n_bands, ms_rows, ms_cols = ms.shape
    ms_f = np.asfortranarray(ms.transpose(1, 2, 0), dtype=np.float64)
    return ms_f.ctypes.data_as(_DBL_P), ms_f, n_bands, ms_rows, ms_cols


def _allocate_out(rows: int, cols: int, n_bands: int) -> tuple[ctypes.POINTER, np.ndarray]:
    out_f = np.zeros((rows, cols, n_bands), dtype=np.float64, order="F")
    return out_f.ctypes.data_as(_DBL_P), out_f


def _out_to_numpy(out_f: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(out_f.transpose(2, 0, 1))


# Public API

Algorithm = Literal["gram_schmidt", "ihs", "wavelet"]


def pansharpen(
        pan: np.ndarray,
        ms: np.ndarray,
        algorithm: Algorithm = "gram_schmidt",
) -> np.ndarray:
    if pan.ndim != 2:
        raise ValueError(
            f"pan must be 2-D (rows, cols), got shape {pan.shape}"
        )

    rows, cols = pan.shape
    pan_ptr, pan_ref = _as_fortran_double(pan)
    ms_ptr, ms_ref, n_bands, ms_rows, ms_cols = _ms_to_fortran(ms)

    if algorithm == "ihs" and n_bands != 3:
        raise ValueError(
            f"IHS pansharpening requires exactly 3 bands, got {n_bands}"
        )

    out_ptr, out_f = _allocate_out(rows, cols, n_bands)

    lib = _get_lib()

    if algorithm == "gram_schmidt":
        lib.gram_schmidt_sharpen(
            pan_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
            ms_ptr, ctypes.c_int(ms_rows), ctypes.c_int(ms_cols),
            ctypes.c_int(n_bands),
            out_ptr,
        )
    elif algorithm == "ihs":
        lib.ihs_sharpen(
            pan_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
            ms_ptr, ctypes.c_int(ms_rows), ctypes.c_int(ms_cols),
            out_ptr,
        )
    elif algorithm == "wavelet":
        lib.wavelet_sharpen(
            pan_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
            ms_ptr, ctypes.c_int(ms_rows), ctypes.c_int(ms_cols),
            ctypes.c_int(n_bands),
            out_ptr,
        )
    else:
        raise ValueError(
            f"Unknown algorithm '{algorithm}'. "
            "Choose from: gram_schmidt, ihs, wavelet"
        )

    del pan_ref, ms_ref

    return _out_to_numpy(out_f)


def pansharpening_gs(pan: np.ndarray, ms: np.ndarray) -> np.ndarray:
    return pansharpen(pan, ms, algorithm="gram_schmidt")


def pansharpening_ihs(pan: np.ndarray, ms: np.ndarray) -> np.ndarray:
    return pansharpen(pan, ms, algorithm="ihs")


def pansharpening_wavelet(pan: np.ndarray, ms: np.ndarray) -> np.ndarray:
    return pansharpen(pan, ms, algorithm="wavelet")
