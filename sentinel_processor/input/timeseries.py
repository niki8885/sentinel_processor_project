from __future__ import annotations
import json
import logging
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence
import numpy as np
import rioxarray  # noqa: F401
import xarray as xr

from sentinel_processor.validation._fortran_bridge import (
    call_check_dimensions,
    call_check_radiometry,
    call_validate_scl,
)
from sentinel_processor.processing._raster_ops_bridge import (
    reproject_nearest as _ft_reproject_nearest,
)

from sentinel_processor.indices.compute import _BAND_ALIASES  # noqa: WPS450

logger = logging.getLogger(__name__)

_TS_RE = re.compile(r"(\d{8}T\d{6})")


def _ts_from_path(p: Path) -> datetime:
    m = _TS_RE.search(p.stem)
    if m:
        return datetime.strptime(m.group(1), "%Y%m%dT%H%M%S").replace(
            tzinfo=timezone.utc
        )
    logger.warning(
        f"[timeseries] No timestamp in filename '{p.name}' — using file mtime"
    )
    return datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc)


def _extract_affine(da: xr.DataArray) -> tuple[float, float, float, float] | None:
    try:
        xs = da.coords["x"].values
        ys = da.coords["y"].values
        if xs.size < 2 or ys.size < 2:
            return None
        pw = float(xs[1] - xs[0])
        ph = float(ys[1] - ys[0])
        if abs(pw) < 1e-12 or abs(ph) < 1e-12:
            return None
        return float(xs[0]), float(ys[0]), pw, ph
    except Exception:
        return None


def _align_to_reference(da: xr.DataArray, ref: xr.DataArray) -> xr.DataArray:
    if da.shape[-2:] == ref.shape[-2:]:
        return da

    same_crs = (
            da.rio.crs is not None
            and ref.rio.crs is not None
            and da.rio.crs == ref.rio.crs
    )
    dst_rows = ref.shape[-2]
    dst_cols = ref.shape[-1]
    dst_affine = _extract_affine(ref)

    if same_crs and dst_affine is not None:
        src_affine = _extract_affine(da)
        if src_affine is not None:
            band_coords = da.coords["band"].values if "band" in da.coords else None
            src_np = np.asarray(da.values, dtype=np.float64)
            if src_np.ndim == 2:
                src_np = src_np[np.newaxis]

            out_np = np.stack(
                [
                    _ft_reproject_nearest(
                        src_np[b], src_affine, dst_rows, dst_cols, dst_affine
                    )
                    for b in range(src_np.shape[0])
                ],
                axis=0,
            )

            coords: dict = {"y": ref.coords["y"], "x": ref.coords["x"]}
            if band_coords is not None:
                coords["band"] = band_coords

            result = xr.DataArray(
                out_np, dims=["band", "y", "x"], coords=coords, attrs=da.attrs,
            )
            if da.rio.crs is not None:
                result = result.rio.write_crs(da.rio.crs)
            return result

    logger.debug(f"[timeseries] reproject_match fallback (same_crs={same_crs})")
    return da.rio.reproject_match(ref)


def _scl_metrics(scl_path: Path, max_cloud_threshold: float) -> dict:
    try:
        da = rioxarray.open_rasterio(str(scl_path), lock=False)
        arr = da.values
    except Exception as exc:
        logger.warning(f"[timeseries] Cannot read SCL {scl_path}: {exc}")
        return {
            "cloud_fraction": 0.0, "snow_fraction": 0.0,
            "confidence_score": 1.0, "dimension_pass": True,
            "radiometry_pass": True, "issues": [],
            "warning": f"SCL unreadable ({exc}); quality check skipped",
        }

    rows = int(arr.shape[-2])
    cols = int(arr.shape[-1])

    dim_result = call_check_dimensions(rows, cols)
    if not dim_result["passed"]:
        return {
            "cloud_fraction": 0.0, "snow_fraction": 0.0,
            "confidence_score": 0.0, "dimension_pass": False,
            "radiometry_pass": False, "issues": dim_result["issues"],
        }

    flat_int = np.asarray(arr, dtype=np.int32).ravel()
    flat_dbl = flat_int.astype(np.float64)

    scl_result = call_validate_scl(flat_int, max_cloud_threshold)

    radio_pass = True
    try:
        radio_pass = call_check_radiometry(flat_dbl)
    except Exception as exc:
        logger.debug(f"[timeseries] Radiometry check error {scl_path.name}: {exc}")

    return {
        "cloud_fraction": scl_result["cloud_ratio"],
        "snow_fraction": scl_result["snow_ratio"],
        "confidence_score": scl_result["confidence_score"],
        "dimension_pass": True,
        "radiometry_pass": radio_pass,
        "issues": scl_result["issues"],
    }


def _read_report(report_path: Path) -> dict | None:
    try:
        with open(report_path) as fh:
            return json.load(fh)
    except Exception:
        return None


def _report_for_scene(spectral_path: Path, scl_dir: str | Path | None) -> dict | None:
    candidates = [spectral_path.parent / f"{spectral_path.stem}_report.json"]
    if scl_dir:
        candidates.append(Path(scl_dir) / f"{spectral_path.stem}_report.json")
    for p in candidates:
        r = _read_report(p)
        if r is not None:
            return r
    return None


def _scl_path_for_scene(spectral_path: Path, scl_dir: str | Path | None) -> Path | None:
    stem = spectral_path.stem
    for ext in (".tif", ".nc"):
        candidates: list[Path] = []
        if scl_dir:
            candidates.append(Path(scl_dir) / f"scl_{stem}{ext}")
        candidates.append(spectral_path.parent / f"scl_{stem}{ext}")
        for p in candidates:
            if p.exists():
                return p
    return None


def _resolve_bands(da: xr.DataArray) -> list[str]:
    if "band" not in da.coords:
        return []
    return [_BAND_ALIASES.get(str(raw), str(raw)) for raw in da.coords["band"].values]


@dataclass
class SceneInfo:
    path: Path
    timestamp: datetime
    accepted: bool
    cloud_fraction: float = 0.0
    snow_fraction: float = 0.0
    confidence_score: float = 1.0
    dimension_pass: bool = True
    radiometry_pass: bool = True
    reject_reason: str = ""
    report_source: str = ""

    def as_dict(self) -> dict:
        return {
            "file": str(self.path),
            "timestamp": self.timestamp.isoformat(),
            "accepted": self.accepted,
            "cloud_fraction": self.cloud_fraction,
            "snow_fraction": self.snow_fraction,
            "confidence_score": self.confidence_score,
            "dimension_pass": self.dimension_pass,
            "radiometry_pass": self.radiometry_pass,
            "reject_reason": self.reject_reason,
            "report_source": self.report_source,
        }


# Configuration


@dataclass
class TimeSeriesConfig:
    """All knobs for stack_timeseries().

    Thresholds mirror DownloadConfig so you can reuse the same values
    in both places without surprises.

    Parameters
    ----------
    align : bool
        Reproject every scene onto the reference grid.
        Same CRS  → Fortran reproject_nearest.
        Diff CRS  → rioxarray.reproject_match fallback.
    reference : str | Path | None
        Explicit reference scene for alignment. None = first source.
    max_cloud_fraction : float
        Max allowed cloud fraction (0–1). Mirrors
        DownloadConfig.max_cloud_threshold.
    min_confidence : float
        Min SCL confidence score (0–1). Same scoring table as Fortran
        validate_scl. Mirrors DownloadConfig.min_confidence.
    max_snow_fraction : float
        Max allowed snow/ice fraction (0–1). 1.0 = disabled.
    require_bands : list[str] | None
        Reject scenes missing any of these bands. Accepts any alias
        understood by the indices module ("nir", "B08", "swir16", …).
    use_sidecar_report : bool
        Read *_report.json sidecars before re-reading SCL files.
        False = always recompute from SCL (slower but always fresh).
    run_radiometry_check : bool
        Apply Fortran check_radiometry (< 1 % pixels > 15 000 DN).
        Mirrors DownloadConfig.validate.
    nodata : float
        Fill value for rejected scenes when fill_rejected=True.
    fill_rejected : bool
        Include rejected scenes as nodata planes to keep the time axis
        contiguous. False (default) = drop rejected scenes.
    save_dir : str | Path | None
        If set, the finished stack is saved here automatically after assembly.
        None (default) = do not save.
    save_format : str
        Output format when save_dir is set.
        ``"nc"`` (default) — single NetCDF-4 file ``<save_name>.nc``.
        ``"tif"`` — one GeoTIFF per time step ``<save_name>_<timestamp>.tif``.
    save_name : str | None
        Base filename without extension. None → derived from the first
        source file, e.g. ``"budapest_stack"``.
    """

    align: bool = True
    reference: str | Path | None = None
    max_cloud_fraction: float = 0.30
    min_confidence: float = 0.01
    max_snow_fraction: float = 1.0
    require_bands: list[str] | None = None
    use_sidecar_report: bool = True
    run_radiometry_check: bool = True
    nodata: float = float("nan")
    fill_rejected: bool = False
    save_dir: str | Path | None = None
    save_format: str = "nc"
    save_name: str | None = None


# StackResult


@dataclass
class StackResult:
    """Return value of stack_timeseries().
    stack : xr.DataArray
        Shape (time, band, y, x), dtype float32.
        time coordinate is numpy datetime64[ns] UTC.
    scenes : list[SceneInfo]
        Per-scene metadata for *all* candidates (accepted and rejected).
    """

    stack: xr.DataArray
    scenes: list[SceneInfo]

    @property
    def n_accepted(self) -> int:
        return sum(1 for s in self.scenes if s.accepted)

    @property
    def n_rejected(self) -> int:
        return sum(1 for s in self.scenes if not s.accepted)

    def summary(self) -> str:
        hdr = (
            f"{'File':<45} {'Timestamp':>20} {'Cloud':>7} "
            f"{'Snow':>6} {'Conf':>6} {'Dim':>4} {'Radio':>6} {'':>3}"
        )
        sep = "─" * len(hdr)
        rows = [hdr, sep]
        for s in self.scenes:
            ok = "✓" if s.accepted else "✗"
            note = f"  ← {s.reject_reason}" if not s.accepted else ""
            dim_ok = "✓" if s.dimension_pass else "✗"
            radio_ok = "✓" if s.radiometry_pass else "✗"
            rows.append(
                f"{s.path.name:<45} "
                f"{s.timestamp.strftime('%Y-%m-%dT%H:%M:%S'):>20} "
                f"{s.cloud_fraction:>7.2%} {s.snow_fraction:>6.2%} "
                f"{s.confidence_score:>6.2f} {dim_ok:>4} {radio_ok:>6} "
                f"{ok:>3}{note}"
            )
        rows += [sep,
                 f"Accepted {self.n_accepted} / {len(self.scenes)}  "
                 f"({self.n_rejected} rejected)"]
        return "\n".join(rows)

    def to_json(self, path: str | Path) -> None:
        with open(path, "w") as fh:
            json.dump({"scenes": [s.as_dict() for s in self.scenes]},
                      fh, indent=2, default=str)

    def save(
            self,
            output_dir: str | Path,
            name: str | None = None,
            fmt: str = "nc",
    ) -> list[str]:
        """Save the stack to disk.

        Parameters
        ----------
        output_dir : str | Path
            Directory to write into (created if absent).
        name : str | None
            Base filename without extension.
            None → ``"<first_scene_stem_prefix>_stack"``, e.g. ``"budapest_stack"``.
        fmt : str
            ``"nc"``  — single NetCDF-4 file  ``<name>.nc``  (default).
            ``"tif"`` — one GeoTIFF per time step ``<name>_<YYYYMMDDTHHMMSS>.tif``.

        Returns
        -------
        list[str]
            Absolute paths of written files.
        """
        import os

        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        fmt = fmt.lower().lstrip(".")
        if fmt not in ("nc", "netcdf", "tif", "tiff"):
            raise ValueError(f"[timeseries] save: unsupported fmt '{fmt}'. Use 'nc' or 'tif'.")

        if name is None:
            first = next((s for s in self.scenes if s.accepted), self.scenes[0])
            stem = first.path.stem
            m = re.match(r"^([A-Za-z0-9_]+?)_\d{8}T\d{6}", stem)
            prefix = m.group(1) if m else stem
            name = f"{prefix}_stack"

        written: list[str] = []

        if fmt in ("nc", "netcdf"):
            out_path = out_dir / f"{name}.nc"
            encoding = {
                "time": {
                    "units": "seconds since 1970-01-01",
                    "calendar": "proleptic_gregorian",
                    "dtype": "int64",
                }
            }
            self.stack.to_netcdf(str(out_path), encoding=encoding)
            written.append(str(out_path.resolve()))
            logger.info(f"[timeseries] Stack saved → {out_path}")

        else:
            times = self.stack.coords["time"].values
            for i, t in enumerate(times):
                ts_str = str(t)[:19].replace("-", "").replace("T", "T").replace(":", "")
                out_path = out_dir / f"{name}_{ts_str}.tif"
                plane = self.stack.isel(time=i)
                try:
                    plane.rio.to_raster(str(out_path), dtype="float32")
                except Exception as exc:
                    logger.warning(f"[timeseries] Could not save time step {ts_str}: {exc}")
                    continue
                written.append(str(out_path.resolve()))
            logger.info(f"[timeseries] Stack saved — {len(written)} GeoTIFF(s) → {out_dir}")

        return written


# Scene I/O


def _open_scene(path: Path) -> xr.DataArray:
    suffix = path.suffix.lower()
    if suffix in (".tif", ".tiff"):
        return rioxarray.open_rasterio(str(path), lock=False)
    ds = xr.open_dataset(str(path))
    for var in ds.data_vars:
        da = ds[var]
        if da.ndim >= 3:
            if "band" not in da.dims and da.ndim == 3:
                da = da.rename({da.dims[0]: "band"})
            return da
    raise ValueError(f"[timeseries] No 3-D variable found in {path}")


def _clean(da: xr.DataArray) -> xr.DataArray:
    return da.drop_vars(
        [v for v in da.coords if v == "spatial_ref"], errors="ignore"
    )


# Public entry point


def stack_timeseries(
        sources: Sequence[str | Path],
        scl_dir: str | Path | None = None,
        cfg: TimeSeriesConfig | None = None,
        progress: bool = True,
) -> StackResult:
    """Build a quality-filtered temporal stack from a list of scene files.

    Parameters
    ----------
    sources : Sequence[str | Path]
        Paths to .tif or .nc spectral scenes produced by download_sentinel2().
        Unsorted globs are fine — files are sorted chronologically by the
        timestamp embedded in their filename.
    scl_dir : str | Path | None
        Directory holding SCL files (typically <output_dir>/technical/).
        Used as a fallback when no JSON sidecar is found next to the scene.
        Expected naming: scl_<scene_stem>.tif or .nc
    cfg : TimeSeriesConfig | None
        Quality thresholds. None → defaults (match DownloadConfig defaults).
    progress : bool
        Print per-scene status lines to stdout.

    Returns
    -------
    StackResult
        .stack  xr.DataArray  (time, band, y, x)  float32
        .scenes list[SceneInfo] — all candidates, accepted and rejected

    Raises
    ------
    ValueError
        sources is empty, no scene passed quality filter, or accepted scenes
        have inconsistent shapes after alignment.
    """
    if cfg is None:
        cfg = TimeSeriesConfig()

    paths = [Path(p) for p in sources]
    if not paths:
        raise ValueError("[timeseries] sources list is empty")

    paths.sort(key=_ts_from_path)

    ref_path = Path(cfg.reference) if cfg.reference else paths[0]
    ref_da = _clean(_open_scene(ref_path))

    canonical_required: list[str] | None = None
    if cfg.require_bands:
        canonical_required = [_BAND_ALIASES.get(b, b) for b in cfg.require_bands]

    n_total = len(paths)
    scenes_meta: list[SceneInfo] = []
    accepted_arrays: list[xr.DataArray | None] = []
    accepted_times: list[np.datetime64] = []

    if progress:
        print(f"\n  Stacking {n_total} scene(s)\n")

    for idx, path in enumerate(paths):
        ts = _ts_from_path(path)
        time_np = np.datetime64(ts.replace(tzinfo=None), "ns")

        # quality gate
        quality: dict = {}
        report_source: str = "none"

        if cfg.use_sidecar_report:
            report = _report_for_scene(path, scl_dir)
            if report is not None:
                quality = {
                    "cloud_fraction": report.get("cloud_ratio", 0.0),
                    "snow_fraction": report.get("snow_ratio", 0.0),
                    "confidence_score": report.get("confidence_score", 1.0),
                    "dimension_pass": report.get("dimension_pass", True),
                    "radiometry_pass": report.get("radiometry_pass", True),
                    "issues": report.get("issues", []),
                }
                report_source = "sidecar_json"

        if not quality:
            scl_path = _scl_path_for_scene(path, scl_dir)
            if scl_path:
                quality = _scl_metrics(scl_path, cfg.max_cloud_fraction)
                if not cfg.run_radiometry_check:
                    quality["radiometry_pass"] = True
                report_source = "scl_computed"
            else:
                quality = {
                    "cloud_fraction": 0.0, "snow_fraction": 0.0,
                    "confidence_score": 1.0, "dimension_pass": True,
                    "radiometry_pass": True, "issues": [],
                }
                report_source = "none"
                logger.warning(
                    f"[timeseries] No SCL / report for {path.name} — "
                    "accepted without quality check"
                )

        cloud_frac = float(quality.get("cloud_fraction", 0.0))
        snow_frac = float(quality.get("snow_fraction", 0.0))
        conf_score = float(quality.get("confidence_score", 1.0))
        dim_pass = bool(quality.get("dimension_pass", True))
        radio_pass = bool(quality.get("radiometry_pass", True))

        reject_reason = ""
        if not dim_pass:
            reject_reason = "dimension check failed"
        elif not radio_pass:
            reject_reason = "radiometry check failed"
        elif cloud_frac > cfg.max_cloud_fraction:
            reject_reason = f"cloud {cloud_frac:.1%} > {cfg.max_cloud_fraction:.1%}"
        elif snow_frac > cfg.max_snow_fraction:
            reject_reason = f"snow {snow_frac:.1%} > {cfg.max_snow_fraction:.1%}"
        elif conf_score < cfg.min_confidence:
            reject_reason = f"confidence {conf_score:.2f} < {cfg.min_confidence:.2f}"

        accepted = not bool(reject_reason)

        info = SceneInfo(
            path=path, timestamp=ts, accepted=accepted,
            cloud_fraction=cloud_frac, snow_fraction=snow_frac,
            confidence_score=conf_score, dimension_pass=dim_pass,
            radiometry_pass=radio_pass, reject_reason=reject_reason,
            report_source=report_source,
        )
        scenes_meta.append(info)

        if progress:
            status = "✓" if accepted else f"✗  {reject_reason}"
            label = path.name if len(path.name) <= 40 else path.name[:37] + "..."
            pad = len(str(n_total))
            print(f"  [{idx + 1:>{pad}}/{n_total}]  {label:<40}  {status}")

        if not accepted:
            if cfg.fill_rejected:
                accepted_arrays.append(None)
                accepted_times.append(time_np)
            continue

        # open and validate bands
        try:
            da = _clean(_open_scene(path))
        except Exception as exc:
            logger.warning(f"[timeseries] Cannot open {path}: {exc}")
            info.accepted = False
            info.reject_reason = f"read error: {exc}"
            continue

        if canonical_required:
            scene_canonical = _resolve_bands(da)
            missing = [b for b in canonical_required if b not in scene_canonical]
            if missing:
                info.accepted = False
                info.reject_reason = f"missing bands: {missing}"
                logger.warning(f"[timeseries] {path.name} — missing bands {missing}")
                continue

        # align to reference gridd
        if cfg.align and path != ref_path:
            try:
                da = _align_to_reference(da, ref_da)
            except Exception as exc:
                logger.warning(f"[timeseries] Alignment failed {path.name}: {exc}")
                info.accepted = False
                info.reject_reason = f"alignment error: {exc}"
                continue

        accepted_arrays.append(da.astype(np.float32))
        accepted_times.append(time_np)

    # assemble stack
    if not any(a is not None for a in accepted_arrays):
        raise ValueError(
            "[timeseries] No scenes passed quality filter. "
            "Relax TimeSeriesConfig thresholds or check SCL / report files."
        )

    shapes = {a.shape for a in accepted_arrays if a is not None}
    if len(shapes) > 1:
        raise ValueError(
            f"[timeseries] Accepted scenes have inconsistent shapes: {shapes}. "
            "Set align=True or provide a reference scene."
        )

    template_da = next(a for a in accepted_arrays if a is not None)
    nodata_plane = np.full(template_da.shape, cfg.nodata, dtype=np.float32)

    planes: list[xr.DataArray] = [
        xr.DataArray(
            nodata_plane.copy(),
            dims=template_da.dims,
            coords={k: v for k, v in template_da.coords.items() if k != "time"},
            attrs=template_da.attrs,
        ) if plane is None else plane
        for plane in accepted_arrays
    ]

    stack = xr.concat(planes, dim="time").assign_coords(
        time=xr.DataArray(
            np.array(accepted_times, dtype="datetime64[ns]"), dims=["time"]
        )
    )

    if stack.dims != ("time", "band", "y", "x"):
        try:
            stack = stack.transpose("time", "band", "y", "x")
        except ValueError:
            pass

    try:
        if template_da.rio.crs is not None:
            stack = stack.rio.write_crs(template_da.rio.crs)
    except Exception:
        pass

    result = StackResult(stack=stack, scenes=scenes_meta)

    # save
    if cfg.save_dir is not None:
        fmt = (cfg.save_format or "nc").lower().lstrip(".")
        saved = result.save(
            output_dir=cfg.save_dir,
            name=cfg.save_name,
            fmt=fmt,
        )
        if progress and saved:
            label = saved[0] if len(saved) == 1 else f"{len(saved)} file(s) in {cfg.save_dir}"
            print(f"  Saved → {label}")

    if progress:
        n_ok = sum(1 for s in scenes_meta if s.accepted)
        print(f"\n  Stack ready — {n_ok}/{n_total} accepted,  shape {tuple(stack.shape)}\n")

    return result
