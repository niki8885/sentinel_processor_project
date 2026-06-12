from __future__ import annotations
import ctypes
import sys
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypedDict
import logging
import numpy as np

if TYPE_CHECKING:
    import xarray as xr

_LIB_NAME = (
    "libsentinel_stats.dll" if sys.platform == "win32"
    else "libsentinel_stats.so"
)
_LIB_PATH = Path(__file__).parent / "fortran" / _LIB_NAME


def _register_dll_directories() -> None:
    """Add MinGW/MSYS2 runtime directories so Windows can resolve DLL dependencies.

    Searches in order:
      1. The fortran/ directory itself (the .dll lives there)
      2. Known MSYS2 / MinGW install locations
      3. Every directory in PATH that contains gfortran.exe
    All found directories are registered; execution continues even if none exist.
    """
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
            "sentinel_processor/processing/fortran/libsentinel_stats.dll"
            if sys.platform == "win32"
            else "sentinel_processor/processing/fortran/libsentinel_stats.so"
        )
        raise FileNotFoundError(
            f"Sentinel stats library not found: {_LIB_PATH}\n"
            "Build with:\n"
            f"  gfortran -O2 -shared -fPIC -o {_out} sentinel_processor/analysis/fortran/sentinel_stats.f90\n"
            "Or use the build script:\n"
            "  python build_sentinel_stats.py"
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
            "  python build_sentinel_stats.py"
        ) from _e

    def _reg(name: str, *argtypes):
        fn = getattr(lib, name)
        fn.restype = None
        fn.argtypes = list(argtypes)

    _reg("time_window_stats",
         _DBL, _DBL, _INT, _INT, _INT, _INT,  # arr, dates, n,r,c, window_days
         _DBL, _DBL, _DBL)  # mean_out, std_out, slope_out

    _reg("anomaly_zscore",
         _DBL, _DBL, _DBL, _INT, _INT, _INT,  # arr, mean_in, std_in, n,r,c
         _DBL)  # zscore_out

    _reg("trend_theil_sen",
         _DBL, _DBL, _INT, _INT, _INT,  # arr, dates, n,r,c
         _DBL, _DBL)  # slope_out, intercept_out

    _reg("valid_obs_count",
         _DBL, _INT, _INT, _INT,  # arr, n,r,c
         _DBL, _DBL)  # count_out, fraction_out

    _reg("temporal_gap_stats",
         _DBL, _DBL, _INT, _INT, _INT,  # arr, dates, n,r,c
         _DBL, _DBL)  # max_gap_out, mean_gap_out

    _reg("pixel_quantiles",
         _DBL, _INT, _INT, _INT,  # arr, n,r,c
         _DBL, _DBL, _DBL, _DBL, _DBL)  # p10,p25,p50,p75,p90

    _reg("pixel_iqr",
         _DBL, _INT, _INT, _INT,  # arr, n,r,c
         _DBL, _DBL)  # iqr_out, outlier_out

    _reg("mann_kendall",
         _DBL, _INT, _INT, _INT,  # arr, n,r,c
         _DBL, _DBL, _DBL, _DBL)  # S, varS, Z, trend

    _reg("bfast_breakpoint",
         _DBL, _DBL, _INT, _INT, _INT,  # arr, dates, n,r,c
         _DBL, _DBL, _DBL)  # break_day, magnitude, rss_ratio

    _reg("phenology_doy",
         _DBL, _DBL, _INT, _INT, _INT,  # arr, dates, n,r,c
         _INT, _INT,  # rising_pct, falling_pct
         _DBL, _DBL, _DBL)  # sos, pos, eos

    _reg("pearson_map",
         _DBL, _DBL, _INT, _INT, _INT,  # arr_a, arr_b, n,r,c
         _DBL)  # r_out

    _reg("pixel_regression",
         _DBL, _DBL, _INT, _INT, _INT,  # y_arr, x_arr, n,r,c
         _DBL, _DBL, _DBL)  # slope_out, intercept_out, r2_out

    _reg("phenology_doy_v2",
         _DBL, _DBL, _INT, _INT, _INT,  # arr, dates, n,r,c
         _INT, _INT,  # rising_pct, falling_pct
         _DBL, _DBL, _DBL, _DBL)  # sos, eos, peak_doy, peak_val

    _reg("savgol_smooth_stack",
         _DBL, _INT, _INT, _INT, _INT)  # arr (inout), n,r,c, window

    _lib = lib
    return _lib


# Internal helpers

def _c64(a) -> np.ndarray:
    """Return a C-contiguous float64 copy."""
    return np.array(a, dtype=np.float64, order="C", copy=True)


def _ptr(a: np.ndarray):
    return a.ctypes.data_as(_DBL)


def _days(dates: list[datetime]) -> np.ndarray:
    t0 = dates[0].replace(tzinfo=None) if dates[0].tzinfo else dates[0]
    return np.array(
        [((d.replace(tzinfo=None) if d.tzinfo else d) - t0).total_seconds() / 86_400.0
         for d in dates],
        dtype=np.float64,
    )


def _validate(arr: np.ndarray, dates: list[datetime] | None = None,
              name: str = "arr") -> tuple[int, int, int]:
    if arr.ndim != 3:
        raise ValueError(f"{name} must be 3-D (n_times, rows, cols), got {arr.shape}")
    n, r, c = arr.shape
    if dates is not None and len(dates) != n:
        raise ValueError(f"len(dates)={len(dates)} != {name}.shape[0]={n}")
    return n, r, c


def _buf2(r, c) -> np.ndarray:
    return np.empty(r * c, np.float64)


def _buf3(n, r, c) -> np.ndarray:
    return np.empty(n * r * c, np.float64)


# TypedDict result types


class WindowStatsResult(TypedDict):
    mean: np.ndarray  # (rows, cols)
    std: np.ndarray  # (rows, cols)
    slope: np.ndarray  # (rows, cols)  units/day


class TheilSenResult(TypedDict):
    slope: np.ndarray  # (rows, cols)  units/day
    intercept: np.ndarray  # (rows, cols)


class ObsCountResult(TypedDict):
    count: np.ndarray  # (rows, cols)
    fraction: np.ndarray  # (rows, cols)  in [0, 1]


class GapStatsResult(TypedDict):
    max_gap: np.ndarray  # (rows, cols)  days
    mean_gap: np.ndarray  # (rows, cols)  days


class QuantileResult(TypedDict):
    p10: np.ndarray  # (rows, cols)
    p25: np.ndarray  # (rows, cols)
    p50: np.ndarray  # (rows, cols)
    p75: np.ndarray  # (rows, cols)
    p90: np.ndarray  # (rows, cols)


class IQRResult(TypedDict):
    iqr: np.ndarray  # (rows, cols)
    outlier: np.ndarray  # (n_times, rows, cols)  1=outlier / 0=ok / NODATA


class MannKendallResult(TypedDict):
    S: np.ndarray  # (rows, cols)
    varS: np.ndarray  # (rows, cols)
    Z: np.ndarray  # (rows, cols)
    trend: np.ndarray  # (rows, cols)  +1 / -1 / 0 / NODATA


class BreakpointResult(TypedDict):
    break_day: np.ndarray  # (rows, cols)
    magnitude: np.ndarray  # (rows, cols)
    rss_ratio: np.ndarray  # (rows, cols)  < 1 = improvement


class PhenologyResult(TypedDict):
    sos: np.ndarray  # (rows, cols)
    pos: np.ndarray  # (rows, cols)
    eos: np.ndarray  # (rows, cols)


# Public API

def time_window_stats(
        arr: np.ndarray,
        dates: list[datetime],
        window_days: int = 30,
) -> WindowStatsResult:
    """Rolling-window mean, std, and OLS slope.

    The window is centred on the pixel's own valid-date mid-point; all
    observations within ±(window_days/2) days of that centre are included.

    Returns
    -------
    WindowStatsResult  keys: ``mean``, ``std``, ``slope`` — each (rows, cols).
    NODATA pixels have no valid observations in the window.
    """
    if window_days < 1:
        raise ValueError(f"window_days must be >= 1, got {window_days}")
    n, r, c = _validate(arr, dates)
    a = _c64(arr)
    d = np.ascontiguousarray(_days(dates))
    mo, so, slo = _buf2(r, c), _buf2(r, c), _buf2(r, c)
    _get_lib().time_window_stats(
        _ptr(a), _ptr(d), _INT(n), _INT(r), _INT(c), _INT(window_days),
        _ptr(mo), _ptr(so), _ptr(slo),
    )
    sh = (r, c)
    return WindowStatsResult(mean=mo.reshape(sh), std=so.reshape(sh),
                             slope=slo.reshape(sh))


def anomaly_zscore(
        arr: np.ndarray,
        dates: list[datetime],
        window_days: int = 30,
        *,
        background: WindowStatsResult | None = None,
) -> np.ndarray:
    """Per-pixel, per-scene z-score.

    z(t, r, c) = (arr[t, r, c] - mean[r, c]) / std[r, c]

    Parameters
    ----------
    background
        Pre-computed result from :func:`time_window_stats`.  When *None*,
        the background is computed internally with *window_days*.

    Returns
    -------
    np.ndarray  (n_times, rows, cols) — NODATA where std==0 or any value is NODATA.
    """
    n, r, c = _validate(arr, dates)
    if background is None:
        background = time_window_stats(arr, dates, window_days)
    a = _c64(arr)
    mu = np.ascontiguousarray(background["mean"].ravel(), np.float64)
    sg = np.ascontiguousarray(background["std"].ravel(), np.float64)
    zo = _buf3(n, r, c)
    _get_lib().anomaly_zscore(
        _ptr(a), _ptr(mu), _ptr(sg), _INT(n), _INT(r), _INT(c), _ptr(zo),
    )
    return zo.reshape(n, r, c)


def trend_theil_sen(
        arr: np.ndarray,
        dates: list[datetime],
) -> TheilSenResult:
    """Robust Theil-Sen slope and intercept per pixel.

    Slope = median of all pairwise slopes; insensitive to outliers.
    O(n²) per pixel — for stacks > ~150 scenes use ``time_window_stats``
    OLS slope instead.

    Returns
    -------
    TheilSenResult  keys: ``slope`` [units/day], ``intercept`` — each (rows, cols).
    """
    n, r, c = _validate(arr, dates)
    a = _c64(arr)
    d = np.ascontiguousarray(_days(dates))
    sl, ic = _buf2(r, c), _buf2(r, c)
    _get_lib().trend_theil_sen(
        _ptr(a), _ptr(d), _INT(n), _INT(r), _INT(c), _ptr(sl), _ptr(ic),
    )
    sh = (r, c)
    return TheilSenResult(slope=sl.reshape(sh), intercept=ic.reshape(sh))


def valid_obs_count(
        arr: np.ndarray,
) -> ObsCountResult:
    """Count and fraction of valid (non-NODATA) observations per pixel.

    Returns
    -------
    ObsCountResult  keys: ``count`` (int as float64), ``fraction`` in [0,1].
    Use ``fraction < 0.3`` to mask poorly covered pixels before analysis.
    """
    n, r, c = _validate(arr)
    a = _c64(arr)
    co, fo = _buf2(r, c), _buf2(r, c)
    _get_lib().valid_obs_count(_ptr(a), _INT(n), _INT(r), _INT(c), _ptr(co), _ptr(fo))
    sh = (r, c)
    return ObsCountResult(count=co.reshape(sh), fraction=fo.reshape(sh))


def temporal_gap_stats(
        arr: np.ndarray,
        dates: list[datetime],
) -> GapStatsResult:
    """Maximum and mean gap between consecutive valid observations.

    Returns
    -------
    GapStatsResult  keys: ``max_gap``, ``mean_gap`` in days — each (rows, cols).
    NODATA for pixels with fewer than 2 valid observations.
    Use ``max_gap > 60`` to flag risky pixels before interpolation.
    """
    n, r, c = _validate(arr, dates)
    a = _c64(arr)
    d = np.ascontiguousarray(_days(dates))
    mg, mng = _buf2(r, c), _buf2(r, c)
    _get_lib().temporal_gap_stats(
        _ptr(a), _ptr(d), _INT(n), _INT(r), _INT(c), _ptr(mg), _ptr(mng),
    )
    sh = (r, c)
    return GapStatsResult(max_gap=mg.reshape(sh), mean_gap=mng.reshape(sh))


def pixel_quantiles(
        arr: np.ndarray,
) -> QuantileResult:
    """Empirical quantiles p10/p25/p50/p75/p90 per pixel.

    Returns
    -------
    QuantileResult  keys: ``p10``, ``p25``, ``p50``, ``p75``, ``p90``.
    NODATA for pixels with no valid observations.
    """
    n, r, c = _validate(arr)
    a = _c64(arr)
    q10, q25, q50, q75, q90 = (_buf2(r, c) for _ in range(5))
    _get_lib().pixel_quantiles(
        _ptr(a), _INT(n), _INT(r), _INT(c),
        _ptr(q10), _ptr(q25), _ptr(q50), _ptr(q75), _ptr(q90),
    )
    sh = (r, c)
    return QuantileResult(p10=q10.reshape(sh), p25=q25.reshape(sh),
                          p50=q50.reshape(sh), p75=q75.reshape(sh),
                          p90=q90.reshape(sh))


def pixel_iqr(
        arr: np.ndarray,
) -> IQRResult:
    """Interquartile range and Tukey outlier mask (k=1.5) per pixel.

    Returns
    -------
    IQRResult  keys:
        ``iqr``     (rows, cols)            Q75 − Q25
        ``outlier`` (n_times, rows, cols)   1.0 = outlier / 0.0 = ok / NODATA

    Feed ``1 - outlier`` as a validity mask into ``interpolate_gaps`` to
    replace cloud residuals that bypassed SCL filtering.
    """
    n, r, c = _validate(arr)
    a = _c64(arr)
    iq = _buf2(r, c)
    om = _buf3(n, r, c)
    _get_lib().pixel_iqr(_ptr(a), _INT(n), _INT(r), _INT(c), _ptr(iq), _ptr(om))
    return IQRResult(iqr=iq.reshape(r, c), outlier=om.reshape(n, r, c))


def mann_kendall(
        arr: np.ndarray,
        dates: list[datetime],
) -> MannKendallResult:
    """Non-parametric Mann-Kendall monotonic trend test (p ≈ 0.05).

    Returns
    -------
    MannKendallResult  keys:
        ``S``     raw S statistic
        ``varS``  variance of S (tie-corrected)
        ``Z``     standardised Z (continuity-corrected)
        ``trend`` +1.0 increasing / -1.0 decreasing / 0.0 none / NODATA

    NODATA for pixels with fewer than 4 valid observations.
    Pair with ``trend_theil_sen`` to quantify detected trends.
    """
    n, r, c = _validate(arr, dates)
    a = _c64(arr)
    so, vo, zo, to = (_buf2(r, c) for _ in range(4))
    _get_lib().mann_kendall(
        _ptr(a), _INT(n), _INT(r), _INT(c),
        _ptr(so), _ptr(vo), _ptr(zo), _ptr(to),
    )
    sh = (r, c)
    return MannKendallResult(S=so.reshape(sh), varS=vo.reshape(sh),
                             Z=zo.reshape(sh), trend=to.reshape(sh))


def bfast_breakpoint(
        arr: np.ndarray,
        dates: list[datetime],
) -> BreakpointResult:
    """Single structural break detection (simplified BFAST scan).

    Finds the split of the valid-obs sequence that minimises the combined
    OLS RSS of two line segments.

    Returns
    -------
    BreakpointResult  keys:
        ``break_day``  day number of the detected break
        ``magnitude``  signed jump in fitted value (same units as arr)
        ``rss_ratio``  RSS_two / RSS_one  — 0 = sharp break, 1 = no improvement

    NODATA for pixels with fewer than 6 valid observations.
    Rule of thumb: ``rss_ratio < 0.7`` indicates a meaningful structural change.
    """
    n, r, c = _validate(arr, dates)
    a = _c64(arr)
    d = np.ascontiguousarray(_days(dates))
    bdo, mago, rrso = (_buf2(r, c) for _ in range(3))
    _get_lib().bfast_breakpoint(
        _ptr(a), _ptr(d), _INT(n), _INT(r), _INT(c),
        _ptr(bdo), _ptr(mago), _ptr(rrso),
    )
    sh = (r, c)
    return BreakpointResult(break_day=bdo.reshape(sh),
                            magnitude=mago.reshape(sh),
                            rss_ratio=rrso.reshape(sh))


def phenology_doy(
        arr: np.ndarray,
        dates: list[datetime],
        rising_pct: int = 20,
        falling_pct: int = 20,
) -> PhenologyResult:
    """Extract SOS / POS / EOS from a smoothed vegetation index stack.

    Parameters
    ----------
    arr         : smoothed NDVI / EVI / LAI stack (use ``interpolate_gaps``
                  from ``timeseries_mod`` before calling this function)
    rising_pct  : SOS threshold — % of seasonal amplitude above the minimum.
                  20 % is the TIMESAT standard.
    falling_pct : EOS threshold — same convention, may differ from rising_pct.

    Returns
    -------
    PhenologyResult  keys: ``sos``, ``pos``, ``eos`` — each (rows, cols).
    Values are fractional day offsets from ``dates[0]``.

    Convert back to calendar dates::

        from datetime import timedelta
        sos_date = dates[0] + timedelta(days=float(result["sos"][r, c]))
    """
    n, r, c = _validate(arr, dates)
    a = _c64(arr)
    d = np.ascontiguousarray(_days(dates))
    sos, pos, eos = (_buf2(r, c) for _ in range(3))
    _get_lib().phenology_doy(
        _ptr(a), _ptr(d), _INT(n), _INT(r), _INT(c),
        _INT(rising_pct), _INT(falling_pct),
        _ptr(sos), _ptr(pos), _ptr(eos),
    )
    sh = (r, c)
    return PhenologyResult(sos=sos.reshape(sh), pos=pos.reshape(sh),
                           eos=eos.reshape(sh))


def pearson_map(
        arr_a: np.ndarray,
        arr_b: np.ndarray,
        dates: list[datetime],
) -> np.ndarray:
    """Pixel-wise Pearson r between two co-registered time stacks.

    Only time steps where *both* stacks have non-NODATA values at the pixel
    are used.

    Returns
    -------
    np.ndarray  (rows, cols) — Pearson r in [-1, 1].
    NODATA when fewer than 3 joint valid observations.
    """
    n, r, c = _validate(arr_a, dates, "arr_a")
    if arr_b.shape != arr_a.shape:
        raise ValueError(f"arr_b.shape {arr_b.shape} != arr_a.shape {arr_a.shape}")
    a = _c64(arr_a)
    b = _c64(arr_b)
    ro = _buf2(r, c)
    _get_lib().pearson_map(
        _ptr(a), _ptr(b), _INT(n), _INT(r), _INT(c), _ptr(ro),
    )
    return ro.reshape(r, c)


class RegressionResult(TypedDict):
    slope: np.ndarray  # (rows, cols)  units / x-unit
    intercept: np.ndarray  # (rows, cols)  units
    r2: np.ndarray  # (rows, cols)  [0, 1]


class PhenologyMetricsResult(TypedDict):
    sos_doy: np.ndarray  # (rows, cols)  start-of-season day number
    eos_doy: np.ndarray  # (rows, cols)  end-of-season day number
    peak_doy: np.ndarray  # (rows, cols)  peak day number
    peak_val: np.ndarray  # (rows, cols)  peak NDVI (or index) value


# pixel_regression

def pixel_regression(
        y: np.ndarray,
        x: np.ndarray,
        *,
        nodata: float = NODATA,
) -> RegressionResult:
    """Per-pixel OLS linear regression against an explicit predictor.

    Parameters
    ----------
    y : (n_times, rows, cols) float64
        Dependent variable (e.g. NDVI stack). Values equal to *nodata*
        are excluded from the regression for that pixel.
    x : (n_times,) float64
        Shared predictor vector — same for every pixel (e.g. day-of-year,
        days-since-start, cumulative temperature).  Must not contain NODATA.

    Returns
    -------
    RegressionResult
        ``slope``     – (rows, cols) regression slope  [units / x-unit]
        ``intercept`` – (rows, cols) regression intercept
        ``r2``        – (rows, cols) coefficient of determination in [0, 1]

    Edge cases
    ----------
    * Pixels with var(x) < 1e-12 (all valid x identical) → slope=0,
      intercept=mean(y), r2=0.
    * Pixels with fewer than 2 valid y observations → all outputs = NODATA.

    Use ``r2 > 0.5`` to mask pixels where the linear model is meaningful
    before visualising slope maps.
    """
    n, r, c = _validate(y, name="y")
    x_arr = np.ascontiguousarray(x, dtype=np.float64)
    if x_arr.shape != (n,):
        raise ValueError(f"x must have shape ({n},), got {x_arr.shape}")

    y_arr = _c64(y)
    if nodata != NODATA:
        y_arr[np.abs(y_arr - nodata) <= 1e-4] = NODATA

    sl = _buf2(r, c)
    ic = _buf2(r, c)
    r2 = _buf2(r, c)

    _get_lib().pixel_regression(
        _ptr(y_arr), _ptr(x_arr),
        _INT(n), _INT(r), _INT(c),
        _ptr(sl), _ptr(ic), _ptr(r2),
    )
    sh = (r, c)
    return RegressionResult(
        slope=sl.reshape(sh),
        intercept=ic.reshape(sh),
        r2=r2.reshape(sh),
    )


# phenology_metrics

_log = logging.getLogger(__name__)


def phenology_metrics(
        ndvi_stack: "np.ndarray | xr.DataArray",
        dates: list[datetime],
        *,
        smooth: bool = True,
        savgol_window: int = 5,
        rising_pct: int = 20,
        falling_pct: int = 20,
        min_valid: int = 5,
) -> PhenologyMetricsResult:
    """Extract per-pixel phenological metrics from a smoothed NDVI stack.

    Implements the TIMESAT convention:
    * **SOS** — first date the smoothed NDVI crosses *rising_pct* % of the
      seasonal amplitude above the season minimum.
    * **EOS** — last date the smoothed NDVI crosses *falling_pct* % of the
      amplitude above the season minimum (after the peak).
    * **peak_doy** — day of the global NDVI maximum.
    * **peak_val** — NDVI value at the peak.

    Parameters
    ----------
    ndvi_stack : (n_times, rows, cols) array or xr.DataArray
        Vegetation index time stack (NDVI, EVI, LAI …).  Values equal to
        NODATA (-9999) are treated as cloud/missing.
    dates : list[datetime]
        Acquisition timestamps, length must equal ``ndvi_stack.shape[0]``.
    smooth : bool, default True
        Apply a quadratic Savitzky-Golay smoother (SP-3 kernel from
        ``sentinel_stats.f90``) before detecting thresholds.
        Set ``False`` if the stack is already smoothed.
    savgol_window : int, default 5
        S-G window width (odd, ≥ 3). Forced odd if even is passed.
        Ignored when ``smooth=False``.
    rising_pct : int, default 20
        SOS threshold as a percentage of seasonal amplitude (TIMESAT default).
    falling_pct : int, default 20
        EOS threshold as a percentage of seasonal amplitude.
    min_valid : int, default 5
        Pixels with fewer valid (non-NODATA) observations return NODATA in
        all outputs.  Must be ≥ 5 per SCRUM-207 acceptance criteria.

    Returns
    -------
    PhenologyMetricsResult
        ``sos_doy``  – (rows, cols) start-of-season, fractional days from
                       ``dates[0]``
        ``eos_doy``  – (rows, cols) end-of-season, same units
        ``peak_doy`` – (rows, cols) day of peak
        ``peak_val`` – (rows, cols) peak index value

        Convert to calendar dates::

            from datetime import timedelta
            sos_date = dates[0] + timedelta(days=float(metrics["sos_doy"][r, c]))

    Notes
    -----
    Pass a gap-filled (not raw) stack for best results.  Cloud gaps produce
    spurious crossings even after smoothing.  Use ``interpolate_gaps`` from
    ``timeseries_mod`` or equivalent before calling this function.
    """
    # ── coerce xarray → numpy ──────────────────────────────────────────────
    try:
        import xarray as xr
        if isinstance(ndvi_stack, xr.DataArray):
            ndvi_stack = ndvi_stack.values
    except ImportError:
        pass

    n, r, c = _validate(ndvi_stack, dates)

    if min_valid < 5:
        _log.warning("min_valid=%d < 5; overriding to 5 per SCRUM-207 AC", min_valid)
        min_valid = 5

    # S-G window
    if savgol_window < 3:
        savgol_window = 3
    if savgol_window % 2 == 0:
        savgol_window += 1

    arr = _c64(ndvi_stack)
    d = np.ascontiguousarray(_days(dates))

    # optional Savitzky-Golay smoothing

    if smooth:
        _get_lib().savgol_smooth_stack(
            _ptr(arr), _INT(n), _INT(r), _INT(c), _INT(savgol_window),
        )

    # phenology extraction

    sos_buf = _buf2(r, c)
    eos_buf = _buf2(r, c)
    pdoy_buf = _buf2(r, c)
    pval_buf = _buf2(r, c)

    _get_lib().phenology_doy_v2(
        _ptr(arr), _ptr(d),
        _INT(n), _INT(r), _INT(c),
        _INT(rising_pct), _INT(falling_pct),
        _ptr(sos_buf), _ptr(eos_buf), _ptr(pdoy_buf), _ptr(pval_buf),
    )

    sh = (r, c)
    return PhenologyMetricsResult(
        sos_doy=sos_buf.reshape(sh),
        eos_doy=eos_buf.reshape(sh),
        peak_doy=pdoy_buf.reshape(sh),
        peak_val=pval_buf.reshape(sh),
    )


# save_phenology

def save_phenology(
        metrics: PhenologyMetricsResult,
        output_dir: "str | Path",
        *,
        crs: "str | None" = None,
        transform: "Any | None" = None,
        nodata: float = NODATA,
        compress: str = "deflate",
) -> dict[str, Path]:
    """Write phenology metrics to four GeoTIFF files.

    Parameters
    ----------
    metrics : PhenologyMetricsResult
        Dict returned by :func:`phenology_metrics`.
    output_dir : str or Path
        Directory where the four GeoTIFFs are written.  Created if missing.
    crs : str or None
        EPSG code or WKT string (e.g. ``"EPSG:32632"``).  Pass ``None`` to
        omit spatial reference — useful when rasterio is unavailable.
    transform : affine.Affine or None
        Geo-transform.  ``None`` → pixel-coordinate identity transform.
    nodata : float
        No-data value written to GeoTIFF metadata (default -9999).
    compress : str
        GDAL compression codec (``"deflate"`` / ``"lzw"`` / ``"none"``).

    Returns
    -------
    dict[str, Path]
        Mapping ``{"sos_doy": Path, "eos_doy": Path, "peak_doy": Path,
        "peak_val": Path}``.

    Raises
    ------
    ImportError
        If ``rasterio`` is not installed.
    """
    try:
        import rasterio
    except ImportError as exc:
        raise ImportError(
            "rasterio is required for save_phenology. "
            "Install with: pip install rasterio"
        ) from exc

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved: dict[str, Path] = {}
    band_map = {
        "sos_doy": metrics["sos_doy"],
        "eos_doy": metrics["eos_doy"],
        "peak_doy": metrics["peak_doy"],
        "peak_val": metrics["peak_val"],
    }

    rows, cols = next(iter(band_map.values())).shape

    default_transform = rasterio.transform.from_bounds(0, 0, cols, rows, cols, rows)
    geo_transform = transform if transform is not None else default_transform

    profile = {
        "driver": "GTiff",
        "dtype": "float64",
        "width": cols,
        "height": rows,
        "count": 1,
        "nodata": nodata,
        "transform": geo_transform,
        "compress": compress,
    }
    if crs is not None:
        profile["crs"] = crs

    for name, data in band_map.items():
        path = output_dir / f"{name}.tif"
        with rasterio.open(path, "w", **profile) as dst:
            dst.write(data.astype(np.float64), 1)
        saved[name] = path
        _log.info("Saved %s → %s", name, path)

    return saved
