from __future__ import annotations
import ctypes
import sys
from pathlib import Path
from typing import Literal

import numpy as np

_LIB_NAME = (
    "libsentinel_filters.dll" if sys.platform == "win32"
    else "libsentinel_filters.so"
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

_lib: ctypes.CDLL | None = None
_DBL_P = ctypes.POINTER(ctypes.c_double)


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    if not _LIB_PATH.exists():
        raise FileNotFoundError(
            f"Filters library not found: {_LIB_PATH}\n"
            "Build (Linux/Mac):\n"
            "  gfortran -O2 -march=native -shared -fPIC "
            "-o libsentinel_filters.so filters.f90\n"
            "Build (Windows):\n"
            "  gfortran -O2 -shared "
            "-o libsentinel_filters.dll filters.f90"
        )
    lib = ctypes.CDLL(str(_LIB_PATH))
    c_int = ctypes.c_int
    c_dbl = ctypes.c_double

    # convolve2d(src, rows, cols, kernel, krows, kcols, dst)
    lib.convolve2d.restype = None
    lib.convolve2d.argtypes = [
        _DBL_P, c_int, c_int, _DBL_P, c_int, c_int, _DBL_P,
    ]

    # gaussian_blur(src, rows, cols, sigma, kradius, dst)
    lib.gaussian_blur.restype = None
    lib.gaussian_blur.argtypes = [
        _DBL_P, c_int, c_int, c_dbl, c_int, _DBL_P,
    ]

    # sobel_magnitude(src, rows, cols, use_l2, dst)
    lib.sobel_magnitude.restype = None
    lib.sobel_magnitude.argtypes = [_DBL_P, c_int, c_int, c_int, _DBL_P]

    # sobel_direction(src, rows, cols, dst)
    lib.sobel_direction.restype = None
    lib.sobel_direction.argtypes = [_DBL_P, c_int, c_int, _DBL_P]

    # laplacian(src, rows, cols, connectivity, dst)
    lib.laplacian.restype = None
    lib.laplacian.argtypes = [_DBL_P, c_int, c_int, c_int, _DBL_P]

    # unsharp_mask(src, rows, cols, sigma, kradius, amount, threshold, dst)
    lib.unsharp_mask.restype = None
    lib.unsharp_mask.argtypes = [
        _DBL_P, c_int, c_int, c_dbl, c_int, c_dbl, c_dbl, _DBL_P,
    ]

    # median_filter(src, rows, cols, radius, dst)
    lib.median_filter.restype = None
    lib.median_filter.argtypes = [_DBL_P, c_int, c_int, c_int, _DBL_P]

    # bilateral_filter(src, rows, cols, sigma_s, sigma_r, kradius, dst)
    lib.bilateral_filter.restype = None
    lib.bilateral_filter.argtypes = [
        _DBL_P, c_int, c_int, c_dbl, c_dbl, c_int, _DBL_P,
    ]

    # morpho_erode / morpho_dilate (src, rows, cols, radius, dst)
    for name in ("morpho_erode", "morpho_dilate"):
        fn = getattr(lib, name)
        fn.restype = None
        fn.argtypes = [_DBL_P, c_int, c_int, c_int, _DBL_P]

    # top_hat_white / top_hat_black (src, rows, cols, radius, dst)
    for name in ("top_hat_white", "top_hat_black"):
        fn = getattr(lib, name)
        fn.restype = None
        fn.argtypes = [_DBL_P, c_int, c_int, c_int, _DBL_P]

    _lib = lib
    return _lib


# Internal helpers

def _to_f_f64(arr: np.ndarray) -> tuple[ctypes.POINTER, np.ndarray]:
    a = np.asarray(arr, dtype=np.float64)
    if not a.flags["F_CONTIGUOUS"]:
        a = np.asfortranarray(a)
    return a.ctypes.data_as(_DBL_P), a


def _empty_f(shape: tuple) -> tuple[ctypes.POINTER, np.ndarray]:
    a = np.empty(shape, dtype=np.float64, order="F")
    return a.ctypes.data_as(_DBL_P), a


def _prepare_2d(arr: np.ndarray) -> tuple[np.ndarray, int, int]:
    if arr.ndim != 2:
        raise ValueError(f"Expected 2-D band array, got shape {arr.shape}")
    rows, cols = arr.shape
    return arr, rows, cols


def _out_2d(out_f: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(out_f)


# Public API

def convolve2d(src: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))
    if kernel.ndim != 2:
        raise ValueError(f"kernel must be 2-D, got shape {kernel.shape}")
    krows, kcols = kernel.shape

    src_ptr, src_ref = _to_f_f64(src)
    ker_ptr, ker_ref = _to_f_f64(kernel)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().convolve2d(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ker_ptr, ctypes.c_int(krows), ctypes.c_int(kcols),
        dst_ptr,
    )
    del src_ref, ker_ref
    return _out_2d(dst)


def gaussian_blur(
        src: np.ndarray,
        sigma: float = 1.0,
        kradius: int | None = None,
) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))
    if kradius is None:
        import math
        kradius = max(1, math.ceil(3 * sigma))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().gaussian_blur(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_double(sigma), ctypes.c_int(kradius),
        dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def sobel_magnitude(
        src: np.ndarray,
        norm: Literal["l1", "l2"] = "l2",
) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))
    use_l2 = 1 if norm == "l2" else 0

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().sobel_magnitude(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(use_l2), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def sobel_direction(src: np.ndarray) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().sobel_direction(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def laplacian(
        src: np.ndarray,
        connectivity: Literal[4, 8] = 4,
) -> np.ndarray:
    if connectivity not in (4, 8):
        raise ValueError("connectivity must be 4 or 8")
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().laplacian(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(connectivity), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def unsharp_mask(
        src: np.ndarray,
        sigma: float = 1.0,
        amount: float = 1.0,
        kradius: int | None = None,
        threshold: float = 0.0,
) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))
    if kradius is None:
        import math
        kradius = max(1, math.ceil(3 * sigma))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().unsharp_mask(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_double(sigma), ctypes.c_int(kradius),
        ctypes.c_double(amount), ctypes.c_double(threshold),
        dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def median_filter(src: np.ndarray, radius: int = 1) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().median_filter(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(radius), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def bilateral_filter(
        src: np.ndarray,
        sigma_s: float = 2.0,
        sigma_r: float = 0.1,
        kradius: int | None = None,
) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))
    if kradius is None:
        import math
        kradius = max(1, math.ceil(2 * sigma_s))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().bilateral_filter(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_double(sigma_s), ctypes.c_double(sigma_r),
        ctypes.c_int(kradius), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def morpho_erode(src: np.ndarray, radius: int = 1) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().morpho_erode(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(radius), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def morpho_dilate(src: np.ndarray, radius: int = 1) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().morpho_dilate(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(radius), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def morpho_open(src: np.ndarray, radius: int = 1) -> np.ndarray:
    return morpho_dilate(morpho_erode(src, radius), radius)


def morpho_close(src: np.ndarray, radius: int = 1) -> np.ndarray:
    return morpho_erode(morpho_dilate(src, radius), radius)


def top_hat_white(src: np.ndarray, radius: int = 3) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().top_hat_white(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(radius), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def top_hat_black(src: np.ndarray, radius: int = 3) -> np.ndarray:
    src, rows, cols = _prepare_2d(np.asarray(src, dtype=np.float64))

    src_ptr, src_ref = _to_f_f64(src)
    dst_ptr, dst = _empty_f((rows, cols))

    _get_lib().top_hat_black(
        src_ptr, ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(radius), dst_ptr,
    )
    del src_ref
    return _out_2d(dst)


def apply_to_bands(
        fn,
        arr: np.ndarray,
        **kwargs,
) -> np.ndarray:
    if arr.ndim == 2:
        return fn(arr, **kwargs)
    if arr.ndim != 3:
        raise ValueError(f"Expected 2-D or 3-D array, got shape {arr.shape}")
    return np.stack([fn(arr[b], **kwargs) for b in range(arr.shape[0])], axis=0)
