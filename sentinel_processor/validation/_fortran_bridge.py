from __future__ import annotations
import ctypes
import sys
from pathlib import Path

_MAX_ISSUE_LEN = 128
_MAX_ISSUES = 4

_LIB_NAME = "libsentinel_validation.dll" if sys.platform == "win32" else "libsentinel_validation.so"
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
            f"Build (Linux/Mac): gfortran -O2 -shared -fPIC -o libsentinel_validation.so validation.f90\n"
            f"Build (Windows):   gfortran -O2 -shared -o libsentinel_validation.dll validation.f90"
        )
    _lib = ctypes.CDLL(str(_LIB_PATH))
    _lib.validate_scl.restype = None
    _lib.validate_scl.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    _lib.check_radiometry.restype = None
    _lib.check_radiometry.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    _lib.check_dimensions.restype = None
    _lib.check_dimensions.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    return _lib


def _to_c_int_array(arr):
    buf = (ctypes.c_int * len(arr))(*arr)
    return buf, len(arr)


def _to_c_double_array(arr):
    try:
        import numpy as np
        if isinstance(arr, np.ndarray):
            flat = arr.flatten().astype(np.float64)
            return flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), flat.size
    except ImportError:
        pass
    buf = (ctypes.c_double * len(arr))(*arr)
    return buf, len(arr)


def call_validate_scl(scl_values, max_cloud_threshold: float = 0.30) -> dict:
    lib = _get_lib()
    c_arr, n = _to_c_int_array(scl_values)
    confidence = ctypes.c_double(0.0)
    cloud_r = ctypes.c_double(0.0)
    snow_r = ctypes.c_double(0.0)
    water_excl = ctypes.c_int(0)
    n_issues = ctypes.c_int(0)
    issues_buf = ctypes.create_string_buffer(_MAX_ISSUE_LEN * _MAX_ISSUES)
    lib.validate_scl(
        c_arr, ctypes.c_int(n), ctypes.c_double(max_cloud_threshold),
        ctypes.byref(confidence), ctypes.byref(cloud_r), ctypes.byref(snow_r),
        ctypes.byref(water_excl), issues_buf, ctypes.c_int(_MAX_ISSUE_LEN),
        ctypes.byref(n_issues),
    )
    issues = []
    for i in range(n_issues.value):
        chunk = issues_buf.raw[i * _MAX_ISSUE_LEN: (i + 1) * _MAX_ISSUE_LEN]
        text = chunk.split(b"\x00", 1)[0].decode("utf-8", errors="replace")
        if text:
            issues.append(text)
    return {
        "confidence_score": confidence.value,
        "cloud_ratio": cloud_r.value,
        "snow_ratio": snow_r.value,
        "water_excluded": bool(water_excl.value),
        "issues": issues,
    }


def call_check_radiometry(pixels) -> bool:
    lib = _get_lib()
    c_arr, n = _to_c_double_array(pixels)
    result = ctypes.c_int(0)
    lib.check_radiometry(c_arr, ctypes.c_int(n), ctypes.byref(result))
    return bool(result.value)


def call_check_dimensions(rows: int, cols: int) -> dict:
    lib = _get_lib()
    result = ctypes.c_int(0)
    n_issues = ctypes.c_int(0)
    issues_buf = ctypes.create_string_buffer(_MAX_ISSUE_LEN * _MAX_ISSUES)
    lib.check_dimensions(
        ctypes.c_int(rows),
        ctypes.c_int(cols),
        ctypes.byref(result),
        issues_buf,
        ctypes.c_int(_MAX_ISSUE_LEN),
        ctypes.byref(n_issues),
    )
    issues = []
    for i in range(n_issues.value):
        chunk = issues_buf.raw[i * _MAX_ISSUE_LEN: (i + 1) * _MAX_ISSUE_LEN]
        text = chunk.split(b"\x00", 1)[0].decode("utf-8", errors="replace")
        if text:
            issues.append(text)
    return {
        "passed": bool(result.value),
        "issues": issues,
    }


def validate_file(
    scl_path: str,
    max_cloud_threshold: float = 0.30,
    min_confidence: float = 0.50,
) -> dict:

    import rioxarray
    try:
        da = rioxarray.open_rasterio(scl_path)
    except Exception as exc:
        return {"file": scl_path, "passed": False, "error": str(exc)}

    scl_arr = da.values.flatten().astype(int).tolist()
    scl_result = call_validate_scl(scl_arr, max_cloud_threshold)

    radio_pass = True
    try:
        radio_pass = call_check_radiometry(da.values.flatten().astype(float))
    except Exception as exc:
        scl_result["radiometry_error"] = str(exc)

    passed = scl_result["confidence_score"] >= min_confidence and radio_pass
    return {
        "file": scl_path,
        **scl_result,
        "radiometry_pass": radio_pass,
        "passed": passed,
    }