import os
import struct
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

ROOT = Path(__file__).parent

_TARGETS = [
    ("sentinel_processor/validation/fortran",  "validation.f90",    "libsentinel_validation"),
    ("sentinel_processor/indices/fortran",      "indices_mod.f90",   "libsentinel_indices"),
    ("sentinel_processor/processing/fortran",   "raster_ops.f90",    "libsentinel_raster_ops"),
    ("sentinel_processor/processing/fortran",   "pansharpening.f90", "libsentinel_processing"),
    ("sentinel_processor/processing/fortran",   "timeseries_mod.f90","libsentinel_timeseries"),
    ("sentinel_processor/filters/fortran",      "filters.f90",       "libsentinel_filters"),
]

_RUNTIME_DLLS = [
    "libgfortran-5.dll",
    "libgcc_s_seh-1.dll",
    "libwinpthread-1.dll",
]


def _find_gfortran() -> str | None:
    # On macOS, Homebrew installs gfortran-14 (or -13, -12) without a plain symlink
    candidates = [
        "gfortran",
        "gfortran-14", "gfortran-13", "gfortran-12", "gfortran-11",
        "/usr/local/bin/gfortran-14", "/usr/local/bin/gfortran-13",
        "/opt/homebrew/bin/gfortran-14", "/opt/homebrew/bin/gfortran-13",
        "/opt/homebrew/bin/gfortran-12",
    ]
    for name in candidates:
        try:
            r = subprocess.run(
                [name, "--version"],
                capture_output=True, text=True
            )
            if r.returncode == 0:
                return name
        except FileNotFoundError:
            continue
    return None


def _copy_runtime_dlls(dest: Path) -> None:
    import shutil
    search = [
        Path(os.environ.get("MSYS2_UCRT64_PATH", r"C:\msys64\ucrt64")) / "bin",
        Path(r"C:\msys64\ucrt64\bin"),
        Path(r"C:\msys64\mingw64\bin"),
        Path(r"C:\mingw64\bin"),
    ]
    for dll in _RUNTIME_DLLS:
        for d in search:
            src = d / dll
            if src.exists():
                dst = dest / dll
                if not dst.exists():
                    shutil.copy2(str(src), str(dst))
                break


def _compile_fortran(root: Path) -> None:
    gfortran = _find_gfortran()
    if gfortran is None:
        print(
            "\n  [sentinel-processor] gfortran not found — Fortran libraries "
            "will not be compiled.\n"
            "  Install gfortran and run: sentinel-processor-compile\n"
        )
        return

    is_win  = sys.platform == "win32"
    is_mac  = sys.platform == "darwin"
    ext     = ".dll" if is_win else ".so"
    flags   = ["-O2", "-shared"] + ([] if is_win else ["-fPIC"])

    if is_mac:
        # Pass ARCHFLAGS so cross-compilation (arm64/x86_64) works
        for token in os.environ.get("ARCHFLAGS", "").split():
            if token.startswith("-arch"):
                flags.append(token)

    for subdir, src_name, lib_stem in _TARGETS:
        src_path = root / subdir / src_name
        out_path = root / subdir / f"{lib_stem}{ext}"

        if not src_path.exists():
            print(f"  [sentinel-processor] SKIP {src_name} (not found)")
            continue

        cmd = [gfortran] + flags + ["-o", str(out_path), str(src_path)]
        print(f"  [sentinel-processor] Compiling {src_name} ...", end=" ", flush=True)
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0:
            print("OK")
            if is_win:
                _copy_runtime_dlls(root / subdir)
        else:
            print(f"FAILED\n  {r.stderr.strip()}")


class FortranBuildExt(build_ext):
    def build_extension(self, ext):
        _compile_fortran(ROOT)

        import shutil
        build_lib = Path(self.build_lib)
        is_win    = sys.platform == "win32"
        lib_ext   = ".dll" if is_win else ".so"

        for subdir, _src, lib_stem in _TARGETS:
            src_lib = ROOT / subdir / f"{lib_stem}{lib_ext}"
            if not src_lib.exists():
                continue
            dst_dir = build_lib / subdir
            dst_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(str(src_lib), str(dst_dir / src_lib.name))
            if is_win:
                for dll_name in _RUNTIME_DLLS:
                    src_dll = ROOT / subdir / dll_name
                    if src_dll.exists():
                        shutil.copy2(str(src_dll), str(dst_dir / dll_name))

        ext_path = Path(self.get_ext_fullpath(ext.name))
        ext_path.parent.mkdir(parents=True, exist_ok=True)

        if sys.platform == "darwin":
            # Compile a minimal real dylib so delocate --require-archs passes
            import tempfile
            c_src = """
#include <Python.h>
PyMODINIT_FUNC PyInit__sentinel_fortran(void) { return NULL; }
"""
            with tempfile.NamedTemporaryFile(suffix=".c", delete=False, mode="w") as tf:
                tf.write(c_src)
                c_path = tf.name
            try:
                import sysconfig
                inc = sysconfig.get_path("include")
                archflags = os.environ.get("ARCHFLAGS", "").split()
                cmd = ["cc", f"-I{inc}", "-shared", "-fPIC",
                       "-undefined", "dynamic_lookup",
                       "-o", str(ext_path), c_path] + archflags
                r = subprocess.run(cmd, capture_output=True, text=True)
                if r.returncode != 0:
                    ext_path.write_bytes(b"\xca\xfe\xba\xbe" + b"\x00" * 28)
            finally:
                os.unlink(c_path)
        elif sys.platform == "win32":
            ext_path.write_bytes(b"MZ" + b"\x00" * 62)
        else:
            ext_path.write_bytes(
                b"\x7fELF\x02\x01\x01\x00" + b"\x00" * 8 +
                b"\x03\x00" + b"\x3e\x00" + b"\x01\x00\x00\x00" + b"\x00" * 24
            )

_dummy_ext = Extension(
    name="sentinel_processor._sentinel_fortran",
    sources=[],
)

setup(
    ext_modules=[_dummy_ext],
    cmdclass={"build_ext": FortranBuildExt},
)