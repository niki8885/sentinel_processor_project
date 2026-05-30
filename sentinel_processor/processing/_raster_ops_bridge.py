from __future__ import annotations
import ctypes
import sys
from pathlib import Path
import numpy as np


_LIB_NAME = (
    "libsentinel_raster_ops.dll" if sys.platform == "win32"
    else "libsentinel_raster_ops.so"
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


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    if not _LIB_PATH.exists():
        raise FileNotFoundError(
            f"Fortran shared library not found: {_LIB_PATH}\n"
            "Build:\n"
            "  gfortran -O2 -shared -fPIC -o libsentinel_raster_ops.so raster_ops.f90\n"
            "  (Windows) gfortran -O2 -shared -o libsentinel_raster_ops.dll raster_ops.f90"
        )
    _lib = ctypes.CDLL(str(_LIB_PATH))

    dbl_p = ctypes.POINTER(ctypes.c_double)
    int_p = ctypes.POINTER(ctypes.c_int)
    c_int = ctypes.c_int
    c_dbl = ctypes.c_double

    _lib.rgb_to_luminance.restype = None
    _lib.rgb_to_luminance.argtypes = [dbl_p, c_int, c_int, c_int, dbl_p]

    _lib.align_bands.restype = None
    _lib.align_bands.argtypes = [dbl_p, c_int, c_int, dbl_p, c_int, c_int]

    _lib.clip_box_indices.restype = None
    _lib.clip_box_indices.argtypes = [
        c_dbl, c_dbl, c_dbl, c_dbl,
        c_int, c_int,
        c_dbl, c_dbl, c_dbl, c_dbl,
        int_p, int_p, int_p, int_p,
    ]

    _lib.band_stats.restype = None
    _lib.band_stats.argtypes = [dbl_p, c_int, dbl_p, dbl_p, dbl_p, dbl_p]

    _lib.normalize_band.restype = None
    _lib.normalize_band.argtypes = [dbl_p, c_int, c_dbl, c_dbl, c_dbl, c_dbl]

    _lib.reproject_nearest.restype = None
    _lib.reproject_nearest.argtypes = [
        dbl_p, c_int, c_int, c_dbl, c_dbl, c_dbl, c_dbl,
        dbl_p, c_int, c_int, c_dbl, c_dbl, c_dbl, c_dbl,
    ]

    return _lib


def _f64_f(arr: np.ndarray) -> tuple[ctypes.POINTER, np.ndarray]:
    a = np.asfortranarray(arr, dtype=np.float64)
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), a


def _empty_f(shape: tuple) -> tuple[ctypes.POINTER, np.ndarray]:
    a = np.empty(shape, dtype=np.float64, order='F')
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), a


def rgb_to_luminance(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 2:
        return arr.astype(np.float64, copy=False)
    if arr.ndim != 3:
        raise ValueError(f"Expected 2-D or 3-D array, got shape {arr.shape}")

    n_bands, rows, cols = arr.shape

    if n_bands == 1:
        return arr[0].astype(np.float64, copy=False)

    try:
        lib = _get_lib()
    except FileNotFoundError:
        if n_bands == 3:
            w = np.array([0.299, 0.587, 0.114])
        else:
            w = np.full(n_bands, 1.0 / n_bands)
        return np.tensordot(w, arr.astype(np.float64), axes=([0], [0]))

    src_f = np.asfortranarray(arr.transpose(1, 2, 0), dtype=np.float64)
    src_ptr = src_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    dst_ptr, dst = _empty_f((rows, cols))
    lib.rgb_to_luminance(
        src_ptr,
        ctypes.c_int(rows), ctypes.c_int(cols), ctypes.c_int(n_bands),
        dst_ptr,
    )
    del src_f
    return np.ascontiguousarray(dst)


def align_bands(
    src: np.ndarray,
    dst_rows: int,
    dst_cols: int,
) -> np.ndarray:

    if src.ndim != 2:
        raise ValueError(f"src must be 2-D, got shape {src.shape}")

    src_rows, src_cols = src.shape

    try:
        lib = _get_lib()
    except FileNotFoundError:
        r_idx = np.minimum(
            (np.arange(dst_rows) * src_rows / dst_rows).astype(int), src_rows - 1
        )
        c_idx = np.minimum(
            (np.arange(dst_cols) * src_cols / dst_cols).astype(int), src_cols - 1
        )
        return src[np.ix_(r_idx, c_idx)].astype(np.float64)

    src_ptr, src_ref = _f64_f(src)
    dst_ptr, dst = _empty_f((dst_rows, dst_cols))
    lib.align_bands(
        src_ptr,
        ctypes.c_int(src_rows), ctypes.c_int(src_cols),
        dst_ptr,
        ctypes.c_int(dst_rows), ctypes.c_int(dst_cols),
    )
    del src_ref
    return np.ascontiguousarray(dst)


def clip_box_indices(
    origin_x: float, origin_y: float,
    pixel_w: float,  pixel_h: float,
    rows: int,       cols: int,
    min_lon: float,  max_lon: float,
    min_lat: float,  max_lat: float,
) -> tuple[int, int, int, int]:

    try:
        lib = _get_lib()
    except FileNotFoundError:
        c0 = max(int((min_lon - origin_x) / pixel_w) + 1, 1)
        c1 = min(int((max_lon - origin_x) / pixel_w) + 1, cols)
        if pixel_h < 0:
            r0 = max(int((max_lat - origin_y) / pixel_h) + 1, 1)
            r1 = min(int((min_lat - origin_y) / pixel_h) + 1, rows)
        else:
            r0 = max(int((min_lat - origin_y) / pixel_h) + 1, 1)
            r1 = min(int((max_lat - origin_y) / pixel_h) + 1, rows)
        return r0, r1, c0, c1

    r_min = ctypes.c_int(0); r_max = ctypes.c_int(0)
    c_min = ctypes.c_int(0); c_max = ctypes.c_int(0)
    lib.clip_box_indices(
        ctypes.c_double(origin_x), ctypes.c_double(origin_y),
        ctypes.c_double(pixel_w),  ctypes.c_double(pixel_h),
        ctypes.c_int(rows),        ctypes.c_int(cols),
        ctypes.c_double(min_lon),  ctypes.c_double(max_lon),
        ctypes.c_double(min_lat),  ctypes.c_double(max_lat),
        ctypes.byref(r_min), ctypes.byref(r_max),
        ctypes.byref(c_min), ctypes.byref(c_max),
    )
    return r_min.value, r_max.value, c_min.value, c_max.value


def band_stats(arr: np.ndarray) -> dict[str, float]:
    flat = arr.ravel().astype(np.float64)
    n = flat.size

    try:
        lib = _get_lib()
    except FileNotFoundError:
        return {
            "mean": float(flat.mean()),
            "std":  float(flat.std()),
            "min":  float(flat.min()),
            "max":  float(flat.max()),
        }

    ptr, ref = _f64_f(flat)
    mean = ctypes.c_double(0.0); std  = ctypes.c_double(0.0)
    mn   = ctypes.c_double(0.0); mx   = ctypes.c_double(0.0)
    lib.band_stats(
        ptr, ctypes.c_int(n),
        ctypes.byref(mean), ctypes.byref(std),
        ctypes.byref(mn),   ctypes.byref(mx),
    )
    del ref
    return {"mean": mean.value, "std": std.value,
            "min": mn.value,   "max": mx.value}


def normalize_band(
    arr: np.ndarray,
    src_min: float | None = None,
    src_max: float | None = None,
    out_min: float = 0.0,
    out_max: float = 1.0,
) -> np.ndarray:
    arr = np.asfortranarray(arr, dtype=np.float64)
    flat = arr.ravel()   # view
    if src_min is None or src_max is None:
        st = band_stats(flat)
        if src_min is None: src_min = st["min"]
        if src_max is None: src_max = st["max"]

    try:
        lib = _get_lib()
        ptr = flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lib.normalize_band(
            ptr, ctypes.c_int(flat.size),
            ctypes.c_double(src_min), ctypes.c_double(src_max),
            ctypes.c_double(out_min), ctypes.c_double(out_max),
        )
    except FileNotFoundError:
        rng = src_max - src_min
        if abs(rng) > 1e-12:
            flat[:] = out_min + (flat - src_min) * ((out_max - out_min) / rng)
        else:
            flat[:] = out_min

    return arr


def reproject_nearest(
    src: np.ndarray,
    src_affine: tuple[float, float, float, float],
    dst_rows: int,
    dst_cols: int,
    dst_affine: tuple[float, float, float, float],
) -> np.ndarray:

    if src.ndim != 2:
        raise ValueError(f"src must be 2-D, got shape {src.shape}")

    src_rows, src_cols = src.shape
    src_ox, src_oy, src_pw, src_ph = src_affine
    dst_ox, dst_oy, dst_pw, dst_ph = dst_affine

    try:
        lib = _get_lib()
    except FileNotFoundError:
        return align_bands(src, dst_rows, dst_cols)

    src_ptr, src_ref = _f64_f(src)
    dst_ptr, dst = _empty_f((dst_rows, dst_cols))

    lib.reproject_nearest(
        src_ptr,
        ctypes.c_int(src_rows), ctypes.c_int(src_cols),
        ctypes.c_double(src_ox), ctypes.c_double(src_oy),
        ctypes.c_double(src_pw), ctypes.c_double(src_ph),
        dst_ptr,
        ctypes.c_int(dst_rows), ctypes.c_int(dst_cols),
        ctypes.c_double(dst_ox), ctypes.c_double(dst_oy),
        ctypes.c_double(dst_pw), ctypes.c_double(dst_ph),
    )
    del src_ref
    return np.ascontiguousarray(dst)