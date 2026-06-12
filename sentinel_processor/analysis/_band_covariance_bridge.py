from __future__ import annotations
import ctypes
import sys
from pathlib import Path
import numpy as np

_LIB_NAME = (
    "libband_covariance.dll" if sys.platform == "win32"
    else "libband_covariance.so"
)
_LIB_PATH = Path(__file__).parent / "fortran" / _LIB_NAME


def _register_dll_directories() -> None:
    import os
    added: list[str] = []

    def _add(path: str) -> None:
        if os.path.isdir(path) and path not in added:
            os.add_dll_directory(path)
            added.append(path)

    _add(str(_LIB_PATH.parent))

    _msys2_roots = [
        r"C:\msys64", r"C:\msys2",
        r"D:\msys64", r"D:\msys2",
    ]
    _mingw_suffixes = [
        r"\ucrt64\bin", r"\mingw64\bin", r"\mingw32\bin", r"\clang64\bin",
    ]
    _standalone = [
        r"C:\mingw64\bin", r"C:\mingw32\bin",
        r"C:\Program Files\mingw-w64\bin",
        r"C:\Program Files (x86)\mingw-w64\bin",
    ]
    for root in _msys2_roots:
        for suf in _mingw_suffixes:
            _add(root + suf)
    for d in _standalone:
        _add(d)

    import shutil
    gfc = shutil.which("gfortran")
    if gfc:
        _add(os.path.dirname(gfc))


if sys.platform == "win32":
    _register_dll_directories()

NODATA: float = -9999.0

_lib: ctypes.CDLL | None = None
_DBL = ctypes.POINTER(ctypes.c_double)
_INT = ctypes.c_int


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is not None:
        return _lib
    if not _LIB_PATH.exists():
        _out = (
            "sentinel_processor/analysis/fortran/libband_covariance.dll"
            if sys.platform == "win32"
            else "sentinel_processor/analysis/fortran/libband_covariance.so"
        )
        raise FileNotFoundError(
            f"Band covariance library not found: {_LIB_PATH}\n"
            "Build with:\n"
            f"  gfortran -O2 -shared -fPIC -o {_out} "
            "sentinel_processor/analysis/fortran/band_covariance.f90\n"
            "Or use the build script:\n"
            "  python build_band_covariance.py"
        )
    try:
        lib = ctypes.CDLL(str(_LIB_PATH))
    except OSError as _e:
        raise OSError(
            f"Failed to load {_LIB_PATH}\n"
            f"Underlying error: {_e}\n"
            "On Windows this usually means a MinGW runtime DLL is missing.\n"
            "Fix: ensure gfortran is in PATH, or install MSYS2 ucrt64:\n"
            "  pacman -S mingw-w64-ucrt-x86_64-gcc-fortran\n"
            "Then recompile with:\n"
            "  python build_band_covariance.py"
        ) from _e

    fn = lib.band_covariance
    fn.restype = None
    fn.argtypes = [
        _DBL,  # arr (n_bands * rows * cols)
        _INT,  # rows
        _INT,  # cols
        _INT,  # n_bands
        _DBL,  # cov_out (n_bands * n_bands)
    ]

    _lib = lib
    return _lib


def _c64(a: np.ndarray) -> np.ndarray:
    """Return a C-contiguous float64 copy."""
    return np.array(a, dtype=np.float64, order="C", copy=True)


def _ptr(a: np.ndarray):
    return a.ctypes.data_as(_DBL)


def _validate_cube(arr: np.ndarray, name: str = "arr") -> tuple[int, int, int]:
    if arr.ndim != 3:
        raise ValueError(
            f"{name} must be 3-D (n_bands, rows, cols), got shape {arr.shape}"
        )
    n_bands, rows, cols = arr.shape
    return n_bands, rows, cols


def band_covariance(arr: np.ndarray) -> np.ndarray:
    """Per-band sample covariance matrix over all spatial pixels.

    Computes the (n_bands x n_bands) covariance matrix of band reflectances
    across the full spatial extent of the image in a single pass through
    the data.  Kahan compensated summation is applied in both the mean and
    the cross-product accumulation steps to minimise floating-point
    cancellation errors.

    Parameters
    ----------
    arr : np.ndarray, shape (n_bands, rows, cols)
        Multi-band raster cube.  Values equal to NODATA (-9999) are excluded
        from all statistics.  For an off-diagonal entry [i, j], a pixel
        contributes only when *both* band i and band j are valid there.

    Returns
    -------
    np.ndarray, shape (n_bands, n_bands), dtype float64
        Symmetric sample covariance matrix (divided by N - 1).
        An entry is NODATA (-9999) when fewer than 2 pixels are jointly
        valid in the corresponding band pair.

    Notes
    -----
    Band means are computed over each band's own valid pixels, while
    cross-products use only jointly valid pixels (pairwise deletion).
    When validity masks differ strongly between bands the result may not
    be positive semi-definite — check eigenvalues before inverting.

    Typical downstream use — PCA via eigendecomposition::

        cov = band_covariance(stack)               # (n_bands, n_bands)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        # eigenvectors[:, -1] is the first principal component direction

    For anomaly detection the inverse covariance (precision matrix) feeds
    the Mahalanobis distance calculation::

        prec = np.linalg.inv(cov)

    Raises
    ------
    ValueError
        If ``arr`` is not exactly 3-D.
    FileNotFoundError
        If the compiled Fortran shared library is not present.
    """
    n_bands, rows, cols = _validate_cube(arr)

    a = _c64(arr)  # band-major, C-contiguous
    cov_buf = np.empty(n_bands * n_bands, dtype=np.float64)

    _get_lib().band_covariance(
        _ptr(a),
        _INT(rows),
        _INT(cols),
        _INT(n_bands),
        _ptr(cov_buf),
    )

    return cov_buf.reshape(n_bands, n_bands)
