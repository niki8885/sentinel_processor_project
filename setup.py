"""
setup.py — compiles all Fortran shared libraries before the wheel is packaged.

cibuildwheel calls  `pip wheel .`  inside each platform container, which
triggers  build_py → BuildPyWithFortran.run() → _compile_fortran().

For development installs (`pip install -e .`) the same hook fires, but
failures are non-fatal — the package works without Fortran (validate=False).
"""

import os
import platform
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py

ROOT = Path(__file__).parent

def _is_win():
    return sys.platform == "win32"

def _is_mac():
    return sys.platform == "darwin"

def _shared_ext():
    return ".dll" if _is_win() else ".so"

def _base_flags():
    flags = ["-O2", "-shared"]
    if not _is_win():
        flags.append("-fPIC")
    if _is_mac() and os.environ.get("ARCHFLAGS"):
        for arch in os.environ["ARCHFLAGS"].split():
            if arch.startswith("-arch"):
                flags += [arch]
    return flags


_TARGETS = [
    ("sentinel_processor/validation/fortran",  "validation.f90",    "libsentinel_validation"),
    ("sentinel_processor/indices/fortran",      "indices_mod.f90",   "libsentinel_indices"),
    ("sentinel_processor/processing/fortran",   "raster_ops.f90",    "libsentinel_raster_ops"),
    ("sentinel_processor/processing/fortran",   "pansharpening.f90", "libsentinel_processing"),
    ("sentinel_processor/filters/fortran",      "filters.f90",       "libsentinel_filters"),
]


_RUNTIME_DLLS = [
    "libgfortran-5.dll",
    "libgcc_s_seh-1.dll",
    "libwinpthread-1.dll",
]

def _copy_runtime_dlls(dest_dir: Path):
    import shutil

    candidates = [
        Path(r"C:\msys64\ucrt64\bin"),
        Path(r"C:\msys64\mingw64\bin"),
        Path(r"C:\mingw64\bin"),
    ]
    msys_root = os.environ.get("MSYS2_UCRT64_PATH", "")
    if msys_root:
        candidates.insert(0, Path(msys_root) / "bin")

    for dll in _RUNTIME_DLLS:
        for candidate in candidates:
            src = candidate / dll
            if src.exists():
                dst = dest_dir / dll
                if not dst.exists():
                    shutil.copy2(str(src), str(dst))
                    print(f"  Copied runtime DLL: {dll}")
                break
        else:
            print(f"  WARNING: runtime DLL not found: {dll}")

def _compile_fortran():
    ext   = _shared_ext()
    flags = _base_flags()

    print(f"\n{'='*60}")
    print(f"Compiling Fortran libraries  (platform={sys.platform})")
    print(f"{'='*60}")

    for subdir, src_name, lib_stem in _TARGETS:
        src_path = ROOT / subdir / src_name
        out_path = ROOT / subdir / f"{lib_stem}{ext}"

        if not src_path.exists():
            print(f"  SKIP  {src_name}  (source not found)")
            continue

        cmd = ["gfortran"] + flags + ["-o", str(out_path), str(src_path)]
        print(f"  {src_name}  →  {lib_stem}{ext}")

        try:
            subprocess.check_call(cmd, stderr=subprocess.STDOUT)
        except FileNotFoundError:
            print(
                "  WARNING: gfortran not found. "
                "Install gfortran and rebuild, or set validate=False at runtime."
            )
            return
        except subprocess.CalledProcessError as exc:
            print(f"  WARNING: compilation failed ({exc}). Skipping.")
            continue

        if _is_win():
            _copy_runtime_dlls(ROOT / subdir)

    print("="*60 + "\n")

class BuildPyWithFortran(build_py):
    def run(self):
        _compile_fortran()
        super().run()


setup(cmdclass={"build_py": BuildPyWithFortran})