from __future__ import annotations
import ctypes
import sys
from pathlib import Path
from typing import Dict, List, Literal
import numpy as np

_LIB_NAME = (
    "libsentinel_wavelet.dll" if sys.platform == "win32"
    else "libsentinel_wavelet.so"
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

_WAVELET_MAP: dict[str, int] = {
    "haar": 0,
    "db4": 1,
    "db6": 2,
    "coif1": 3,
    "sym4": 4,
    "sym6": 5,
}
Wavelet = Literal["haar", "db4", "db6", "coif1", "sym4", "sym6"]


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    if not _LIB_PATH.exists():
        raise FileNotFoundError(
            f"Wavelet library not found: {_LIB_PATH}\n"
            "Build (Linux/macOS):\n"
            "  gfortran -O2 -march=native -shared -fPIC \\\n"
            "    -o <pkg>/wavelet/fortran/libsentinel_wavelet.so \\\n"
            "       <pkg>/wavelet/fortran/wavelet_mod.f90\n"
            "Build (Windows):\n"
            "  gfortran -O2 -shared \\\n"
            "    -o <pkg>/wavelet/fortran/libsentinel_wavelet.dll \\\n"
            "       <pkg>/wavelet/fortran/wavelet_mod.f90"
        )
    lib = ctypes.CDLL(str(_LIB_PATH))
    c_int = ctypes.c_int
    c_dbl = ctypes.c_double

    # dwt2d / idwt2d
    for fn in (lib.dwt2d, lib.idwt2d):
        fn.restype = None
        fn.argtypes = [_DBL_P, c_int, c_int, c_int, c_int, _DBL_P]

    # estimate_sigma(coeffs, n, sigma_out)
    lib.estimate_sigma.restype = None
    lib.estimate_sigma.argtypes = [_DBL_P, c_int, _DBL_P]

    # bayes_threshold(coeffs, n, sigma_n, threshold_out)
    lib.bayes_threshold.restype = None
    lib.bayes_threshold.argtypes = [_DBL_P, c_int, c_dbl, _DBL_P]

    # band_energy(coeffs, n, energy_out)
    lib.band_energy.restype = None
    lib.band_energy.argtypes = [_DBL_P, c_int, _DBL_P]

    # band_stats(coeffs, n, mean, var, l1, linf)
    lib.band_stats.restype = None
    lib.band_stats.argtypes = [_DBL_P, c_int, _DBL_P, _DBL_P, _DBL_P, _DBL_P]

    # dwt2d_batch / idwt2d_batch  (add n_bands argument)
    for fn in (lib.dwt2d_batch, lib.idwt2d_batch):
        fn.restype = None
        fn.argtypes = [_DBL_P, c_int, c_int, c_int, c_int, c_int, _DBL_P]

    # dwt3d / idwt3d  (add n_times argument)
    for fn in (lib.dwt3d, lib.idwt3d):
        fn.restype = None
        fn.argtypes = [_DBL_P, c_int, c_int, c_int, c_int, c_int, _DBL_P]

    _lib = lib
    return _lib


# Internal helpers

def _to_f64_f(arr: np.ndarray) -> tuple[ctypes.POINTER, np.ndarray]:
    """Return (ctypes pointer, F-contiguous float64 copy)."""
    a = np.asfortranarray(arr, dtype=np.float64)
    return a.ctypes.data_as(_DBL_P), a


def _scalar_out() -> tuple[np.ndarray, ctypes.POINTER]:
    a = np.zeros(1, dtype=np.float64)
    return a, a.ctypes.data_as(_DBL_P)


def _validate_arr(arr: np.ndarray, levels: int) -> tuple[int, int]:
    if arr.ndim != 2:
        raise ValueError(
            f"arr must be 2-D (rows, cols), got shape {arr.shape}"
        )
    rows, cols = arr.shape
    factor = 2 ** levels
    if rows % factor != 0:
        raise ValueError(
            f"rows ({rows}) must be divisible by 2^levels = {factor}"
        )
    if cols % factor != 0:
        raise ValueError(
            f"cols ({cols}) must be divisible by 2^levels = {factor}"
        )
    return rows, cols


def _validate_wavelet(wavelet: str) -> int:
    if wavelet not in _WAVELET_MAP:
        raise ValueError(
            f"Unknown wavelet {wavelet!r}. "
            f"Choose from: {list(_WAVELET_MAP)}"
        )
    return _WAVELET_MAP[wavelet]


# Coefficient <-> flat-array layout conversion

def _flat_to_dict(
        flat: np.ndarray,
        rows: int,
        cols: int,
        levels: int,
) -> dict:
    """
    Unpack a flat Fortran-order coefficient array into a structured dict.

    Returns
    -------
    dict  { level (int, 1=finest) :
              { 'LL': ndarray,  (only at level == levels)
                'LH': ndarray, 'HL': ndarray, 'HH': ndarray } }
    """
    result: dict[int, dict[str, np.ndarray]] = {}

    for lv in range(1, levels + 1):
        rh = rows >> lv
        ch = cols >> lv
        rT = rh * 2
        cT = ch * 2

        def _sub(r_start: int, r_end: int, c_start: int, c_end: int) -> np.ndarray:
            block = np.empty((r_end - r_start, c_end - c_start),
                             dtype=np.float64, order="C")
            for j, c in enumerate(range(c_start, c_end)):
                block[:, j] = flat[c * rows + r_start: c * rows + r_end]
            return block

        entry: dict[str, np.ndarray] = {
            "LH": _sub(0, rh, ch, cT),
            "HL": _sub(rh, rT, 0, ch),
            "HH": _sub(rh, rT, ch, cT),
        }
        if lv == levels:
            entry["LL"] = _sub(0, rh, 0, ch)

        result[lv] = entry

    return result


def _dict_to_flat(
        coeffs: dict,
        rows: int,
        cols: int,
        levels: int,
) -> np.ndarray:
    """Pack a coefficient dict back into a flat F-order array."""
    flat = np.zeros(rows * cols, dtype=np.float64)

    def _put(arr: np.ndarray, r_start: int, c_start: int) -> None:
        r_sz, c_sz = arr.shape
        for j in range(c_sz):
            flat[(c_start + j) * rows + r_start:
                 (c_start + j) * rows + r_start + r_sz] = arr[:, j]

    for lv in range(1, levels + 1):
        rh = rows >> lv
        ch = cols >> lv
        d = coeffs[lv]
        _put(d["LH"], 0, ch)
        _put(d["HL"], rh, 0)
        _put(d["HH"], rh, ch)
        if lv == levels:
            _put(d["LL"], 0, 0)

    return flat


# Public API — core transform

def dwt2d(
        arr: np.ndarray,
        levels: int = 1,
        wavelet: Wavelet = "haar",
) -> dict:
    """
    Forward N-level 2-D Discrete Wavelet Transform.

    Parameters
    ----------
    arr     : (rows, cols) float64 array.
              rows and cols must each be divisible by 2^levels.
    levels  : Number of decomposition levels (>= 1).
    wavelet : Wavelet family; see module docstring for options.

    Returns
    -------
    dict  { level: {'LL': ndarray, 'LH': ndarray, 'HL': ndarray, 'HH': ndarray} }
    'LL' is present only at level == levels (coarsest approximation).
    All sub-band arrays are C-order float64.
    """
    rows, cols = _validate_arr(arr, levels)
    wav_flag = _validate_wavelet(wavelet)

    src_ptr, src_buf = _to_f64_f(arr)
    out_buf = np.empty((rows, cols), dtype=np.float64, order="F")
    out_ptr = out_buf.ctypes.data_as(_DBL_P)

    _get_lib().dwt2d(
        src_ptr,
        ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(levels), ctypes.c_int(wav_flag),
        out_ptr,
    )

    return _flat_to_dict(out_buf.ravel(order="F"), rows, cols, levels)


def idwt2d(
        coeffs: dict,
        wavelet: Wavelet = "haar",
) -> np.ndarray:
    """
    Inverse N-level 2-D Discrete Wavelet Transform.

    Parameters
    ----------
    coeffs  : dict returned by :func:`dwt2d` (may be thresholded).
    wavelet : Must match the forward transform.

    Returns
    -------
    np.ndarray  (rows, cols) float64, C-order.
    """
    wav_flag = _validate_wavelet(wavelet)
    levels = max(coeffs.keys())
    ll = coeffs[levels]["LL"]
    rows = ll.shape[0] * (2 ** levels)
    cols = ll.shape[1] * (2 ** levels)

    flat = _dict_to_flat(coeffs, rows, cols, levels)
    in_ptr = flat.ctypes.data_as(_DBL_P)
    out_buf = np.empty(rows * cols, dtype=np.float64)
    out_ptr = out_buf.ctypes.data_as(_DBL_P)

    _get_lib().idwt2d(
        in_ptr,
        ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(levels), ctypes.c_int(wav_flag),
        out_ptr,
    )

    return np.ascontiguousarray(
        out_buf.reshape((rows, cols), order="F")
    )


# Public API — thresholding

def threshold_coeffs(
        coeffs: dict,
        threshold: float,
        mode: Literal["soft", "hard"] = "soft",
) -> dict:
    """
    Apply soft or hard thresholding to all detail sub-bands.

    The LL approximation sub-band is never modified.

    Parameters
    ----------
    coeffs    : dict from :func:`dwt2d`.
    threshold : Non-negative threshold value.
    mode      : 'soft' — shrink toward zero by threshold amount.
                'hard' — zero out coefficients with |x| <= threshold.

    Returns
    -------
    dict  Same structure as input; arrays are new copies.
    """
    if threshold < 0:
        raise ValueError(f"threshold must be >= 0, got {threshold}")
    if mode not in ("soft", "hard"):
        raise ValueError(f"mode must be 'soft' or 'hard', got {mode!r}")

    result: dict = {}
    for lv, sub_bands in coeffs.items():
        result[lv] = {}
        for key, arr in sub_bands.items():
            if key == "LL":
                result[lv]["LL"] = arr.copy()
            else:
                a = arr.copy()
                if mode == "soft":
                    result[lv][key] = np.sign(a) * np.maximum(
                        np.abs(a) - threshold, 0.0
                    )
                else:
                    a[np.abs(a) <= threshold] = 0.0
                    result[lv][key] = a
    return result


# Public API — noise estimation & BayesShrink

def estimate_sigma(band: np.ndarray) -> float:
    """
    Robust noise-sigma estimate via the Donoho-Johnstone MAD estimator.

    Recommended usage: pass the finest-scale HH sub-band.
    sigma = median(|x|) / 0.6745

    Parameters
    ----------
    band : Any flat or 2-D float array (typically coeffs[1]['HH']).

    Returns
    -------
    float  Estimated noise standard deviation.
    """
    b = np.ascontiguousarray(band.ravel(), dtype=np.float64)
    b_ptr = b.ctypes.data_as(_DBL_P)
    sig, sig_ptr = _scalar_out()
    _get_lib().estimate_sigma(b_ptr, ctypes.c_int(len(b)), sig_ptr)
    return float(sig[0])


def bayes_threshold(band: np.ndarray, sigma_n: float) -> float:
    """
    Compute the BayesShrink optimal threshold for one sub-band.

    T = sigma_n^2 / sigma_s,  where sigma_s = sqrt(max(E[x^2] - sigma_n^2, 0))
    If sigma_s ≈ 0 (sub-band is pure noise) returns max(|x|).

    Parameters
    ----------
    band    : 2-D or 1-D sub-band array.
    sigma_n : Noise sigma (from :func:`estimate_sigma`).

    Returns
    -------
    float  Threshold value.
    """
    b = np.ascontiguousarray(band.ravel(), dtype=np.float64)
    b_ptr = b.ctypes.data_as(_DBL_P)
    thr, thr_ptr = _scalar_out()
    _get_lib().bayes_threshold(
        b_ptr, ctypes.c_int(len(b)), ctypes.c_double(sigma_n), thr_ptr
    )
    return float(thr[0])


def bayes_denoise(
        arr: np.ndarray,
        levels: int = 3,
        wavelet: Wavelet = "db4",
) -> np.ndarray:
    """
    Full BayesShrink denoising pipeline.

    Steps:
      1. Forward DWT
      2. Estimate noise sigma from finest HH sub-band
      3. Compute per-sub-band BayesShrink threshold
      4. Soft-threshold each detail sub-band
      5. Inverse DWT

    Parameters
    ----------
    arr     : (rows, cols) float64 image.
    levels  : Decomposition depth.
    wavelet : Wavelet family.

    Returns
    -------
    np.ndarray  Denoised image, same shape as input.
    """
    coeffs = dwt2d(arr, levels=levels, wavelet=wavelet)
    sigma_n = estimate_sigma(coeffs[1]["HH"])

    result: dict = {}
    for lv, sub_bands in coeffs.items():
        result[lv] = {}
        for key, band in sub_bands.items():
            if key == "LL":
                result[lv]["LL"] = band.copy()
            else:
                thr = bayes_threshold(band, sigma_n)
                result[lv][key] = np.sign(band) * np.maximum(
                    np.abs(band) - thr, 0.0
                )
    return idwt2d(result, wavelet=wavelet)


# Public API — sub-band statistics

def band_energy(band: np.ndarray) -> float:
    """
    Squared L2 norm (Parseval energy) of a sub-band array.

    Parameters
    ----------
    band : Any flat or 2-D sub-band array.

    Returns
    -------
    float  Sum of squared coefficients.
    """
    b = np.ascontiguousarray(band.ravel(), dtype=np.float64)
    b_ptr = b.ctypes.data_as(_DBL_P)
    eng, eng_ptr = _scalar_out()
    _get_lib().band_energy(b_ptr, ctypes.c_int(len(b)), eng_ptr)
    return float(eng[0])


def band_stats(band: np.ndarray) -> dict[str, float]:
    """
    Descriptive statistics for a sub-band array.

    Parameters
    ----------
    band : Any flat or 2-D sub-band array.

    Returns
    -------
    dict with keys:
      'mean'  — arithmetic mean
      'var'   — population variance
      'l1'    — mean absolute value (L1 norm / n)
      'linf'  — maximum absolute value (L∞ norm)
    """
    b = np.ascontiguousarray(band.ravel(), dtype=np.float64)
    b_ptr = b.ctypes.data_as(_DBL_P)
    mn, mn_ptr = _scalar_out()
    vr, vr_ptr = _scalar_out()
    l1, l1_ptr = _scalar_out()
    li, li_ptr = _scalar_out()
    _get_lib().band_stats(
        b_ptr, ctypes.c_int(len(b)),
        mn_ptr, vr_ptr, l1_ptr, li_ptr,
    )
    return {
        "mean": float(mn[0]),
        "var": float(vr[0]),
        "l1": float(l1[0]),
        "linf": float(li[0]),
    }


# Public API — batch processing

def dwt2d_batch(
        stack: np.ndarray,
        levels: int = 1,
        wavelet: Wavelet = "haar",
) -> list[dict]:
    """
    Apply 2-D DWT to each band in a multi-spectral stack.

    Parameters
    ----------
    stack   : (n_bands, rows, cols) or (rows, cols, n_bands) float64 array.
              Assumed to be (n_bands, rows, cols) unless the first axis
              size equals the last, in which case band-last is assumed
              only if ndim==3 and shape[0] == shape[1] (ambiguous case
              defaults to band-first).
    levels  : Decomposition depth.
    wavelet : Wavelet family.

    Returns
    -------
    list of dicts, one per band (band-first order).
    """
    if stack.ndim != 3:
        raise ValueError(
            f"stack must be 3-D, got shape {stack.shape}"
        )
    n_bands, rows, cols = stack.shape
    _ = _validate_arr(stack[0], levels)
    wav_flag = _validate_wavelet(wavelet)
    npix = rows * cols

    flat = np.empty(n_bands * npix, dtype=np.float64)
    for b in range(n_bands):
        flat[b * npix: (b + 1) * npix] = (
            np.asfortranarray(stack[b], dtype=np.float64).ravel(order="F")
        )
    out = np.empty_like(flat)
    flat_ptr = flat.ctypes.data_as(_DBL_P)
    out_ptr = out.ctypes.data_as(_DBL_P)

    _get_lib().dwt2d_batch(
        flat_ptr,
        ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(n_bands), ctypes.c_int(levels),
        ctypes.c_int(wav_flag),
        out_ptr,
    )

    result = []
    for b in range(n_bands):
        result.append(_flat_to_dict(out[b * npix: (b + 1) * npix], rows, cols, levels))
    return result


def idwt2d_batch(
        coeffs_list: list[dict],
        wavelet: Wavelet = "haar",
) -> np.ndarray:
    """
    Inverse 2-D DWT applied to a list of coefficient dicts.

    Parameters
    ----------
    coeffs_list : list of dicts from :func:`dwt2d_batch`.
    wavelet     : Must match the forward transform.

    Returns
    -------
    np.ndarray  (n_bands, rows, cols) float64, C-order.
    """
    wav_flag = _validate_wavelet(wavelet)
    n_bands = len(coeffs_list)
    levels = max(coeffs_list[0].keys())
    ll = coeffs_list[0][levels]["LL"]
    rows = ll.shape[0] * (2 ** levels)
    cols = ll.shape[1] * (2 ** levels)
    npix = rows * cols

    flat = np.zeros(n_bands * npix, dtype=np.float64)
    for b, cd in enumerate(coeffs_list):
        flat[b * npix: (b + 1) * npix] = _dict_to_flat(cd, rows, cols, levels)

    out = np.empty_like(flat)
    flat_ptr = flat.ctypes.data_as(_DBL_P)
    out_ptr = out.ctypes.data_as(_DBL_P)

    _get_lib().idwt2d_batch(
        flat_ptr,
        ctypes.c_int(rows), ctypes.c_int(cols),
        ctypes.c_int(n_bands), ctypes.c_int(levels),
        ctypes.c_int(wav_flag),
        out_ptr,
    )

    # Each band b occupies out[b*npix:(b+1)*npix] in F column-major order
    result = np.empty((n_bands, rows, cols), dtype=np.float64)
    for b in range(n_bands):
        result[b] = out[b * npix: (b + 1) * npix].reshape((rows, cols), order="F")
    return np.ascontiguousarray(result)


# Public API — 3-D DWT (spatial + temporal)

def dwt3d(
        arr: np.ndarray,
        levels: int = 1,
        wavelet: Wavelet = "haar",
) -> np.ndarray:
    """
    Forward separable 3-D DWT: 2-D spatial per slice + 1-D temporal.

    Parameters
    ----------
    arr     : (n_times, rows, cols) float64 array.
              rows, cols must each be divisible by 2^levels.
              n_times must also be divisible by 2^levels.
    levels  : Decomposition depth (spatial and temporal).
    wavelet : Wavelet family.

    Returns
    -------
    np.ndarray  (n_times, rows, cols) coefficient array.
    Layout: spatial quad-tree per time slice; temporal lo sub-bands in
    slices 0..n_times//2-1, hi in n_times//2..n_times-1 (per level).
    """
    if arr.ndim != 3:
        raise ValueError(
            f"arr must be 3-D (n_times, rows, cols), got {arr.shape}"
        )
    n_times, rows, cols = arr.shape
    _ = _validate_arr(arr[0], levels)
    wav_flag = _validate_wavelet(wavelet)

    factor = 2 ** levels
    if n_times % factor != 0:
        raise ValueError(
            f"n_times ({n_times}) must be divisible by 2^levels = {factor}"
        )

    # Fortran layout: time slice t (1-based) at flat[(t-1)*npix .. t*npix],
    # each slice in column-major (F) order.
    npix = rows * cols
    flat = np.empty(n_times * npix, dtype=np.float64)
    for t in range(n_times):
        flat[t * npix: (t + 1) * npix] = (
            np.asfortranarray(arr[t], dtype=np.float64).ravel(order="F")
        )
    out = np.empty_like(flat)
    flat_ptr = flat.ctypes.data_as(_DBL_P)
    out_ptr = out.ctypes.data_as(_DBL_P)

    _get_lib().dwt3d(
        flat_ptr,
        ctypes.c_int(rows), ctypes.c_int(cols), ctypes.c_int(n_times),
        ctypes.c_int(levels), ctypes.c_int(wav_flag),
        out_ptr,
    )

    # Unpack: each time slice at out[t*npix:(t+1)*npix], F column-major
    result = np.empty((n_times, rows, cols), dtype=np.float64)
    for t in range(n_times):
        result[t] = out[t * npix: (t + 1) * npix].reshape((rows, cols), order="F")
    return np.ascontiguousarray(result)


def idwt3d(
        coeffs: np.ndarray,
        wavelet: Wavelet = "haar",
        levels: int | None = None,
) -> np.ndarray:
    """
    Inverse separable 3-D DWT.

    Parameters
    ----------
    coeffs  : (n_times, rows, cols) coefficient array from :func:`dwt3d`.
    wavelet : Must match the forward transform.
    levels  : Decomposition depth used in the forward :func:`dwt3d` call.
              Defaults to ``int(log2(n_times))``; pass explicitly when
              ``dwt3d`` was called with a different value.

    Returns
    -------
    np.ndarray  (n_times, rows, cols) reconstructed array.
    """
    if coeffs.ndim != 3:
        raise ValueError(
            f"coeffs must be 3-D (n_times, rows, cols), got {coeffs.shape}"
        )
    n_times, rows, cols = coeffs.shape
    wav_flag = _validate_wavelet(wavelet)

    npix = rows * cols
    if levels is None:
        levels = int(np.log2(n_times))
    flat = np.empty(n_times * npix, dtype=np.float64)
    for t in range(n_times):
        flat[t * npix: (t + 1) * npix] = (
            np.asfortranarray(coeffs[t], dtype=np.float64).ravel(order="F")
        )
    out = np.empty_like(flat)
    flat_ptr = flat.ctypes.data_as(_DBL_P)
    out_ptr = out.ctypes.data_as(_DBL_P)

    _get_lib().idwt3d(
        flat_ptr,
        ctypes.c_int(rows), ctypes.c_int(cols), ctypes.c_int(n_times),
        ctypes.c_int(levels), ctypes.c_int(wav_flag),
        out_ptr,
    )

    result = np.empty((n_times, rows, cols), dtype=np.float64)
    for t in range(n_times):
        result[t] = out[t * npix: (t + 1) * npix].reshape((rows, cols), order="F")
    return np.ascontiguousarray(result)
