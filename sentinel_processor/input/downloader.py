from __future__ import annotations
import datetime
import json
import logging
import os
import warnings
from dataclasses import dataclass, field
from typing import Sequence, Union
from concurrent.futures import ThreadPoolExecutor, as_completed, Future
import numpy as np
import rioxarray
import xarray as xr
from pystac_client import Client

_GDAL_ENV = {
    "GDAL_HTTP_MULTIPLEX": "YES",
    "GDAL_HTTP_VERSION": "2",
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_CACHE_SIZE": "200000000",
    "GDAL_CACHEMAX": "512",
    "GDAL_HTTP_MERGE_CONSECUTIVE_RANGES": "YES",
    "GDAL_MAX_CONNECTIONS": "8",
}
for _k, _v in _GDAL_ENV.items():
    os.environ.setdefault(_k, _v)

from sentinel_processor.config import (  # noqa: E402
    DEFAULT_BBOX_HALF_DEG,
    DEFAULT_KEEP_ITEMS,
    DEFAULT_LOOKBACK_DAYS,
    DEFAULT_MAX_ITEMS,
    DEFAULT_OUTPUT_DIR,
    MAX_CLOUD_THRESHOLD,
    MIN_CONFIDENCE_TO_SAVE,
    SAVE_VALIDATION_REPORT,
    STAC_API_URL,
)
from sentinel_processor.validation._fortran_bridge import (  # noqa: E402
    call_check_radiometry,
    call_check_dimensions,
    call_validate_scl,
)
from sentinel_processor.processing._fortran_bridge import (  # noqa: E402
    pansharpen,
    Algorithm as PanAlgorithm,
)
from sentinel_processor.processing._raster_ops_bridge import (  # noqa: E402
    rgb_to_luminance as _ft_rgb_to_luminance,
    reproject_nearest as _ft_reproject_nearest,
)
from sentinel_processor.utils.data_utils import (  # noqa: E402
    LocationSpec, SpectralBands, TechnicalLayers, VisualAssets, _BandGroup,
)

logger = logging.getLogger(__name__)
logging.getLogger("rasterio.session").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=FutureWarning, module="xarray")

try:
    import rasterio as _r

    _proj = os.path.join(os.path.dirname(_r.__file__), "proj_data")
    if os.path.isdir(_proj):
        os.environ.setdefault("PROJ_DATA", _proj)
        os.environ.setdefault("PROJ_LIB", _proj)
except Exception:
    pass

BandInput = Union[_BandGroup, list[str]]


# Configuration

@dataclass
class DownloadConfig:
    bands: BandInput = field(default_factory=lambda: SpectralBands.ALL)
    tech_bands: BandInput | None = field(default_factory=lambda: TechnicalLayers.ALL)
    visual: bool = True
    output_dir: str = DEFAULT_OUTPUT_DIR
    bbox_half_deg: float = DEFAULT_BBOX_HALF_DEG
    max_items: int = DEFAULT_MAX_ITEMS
    keep_items: int = DEFAULT_KEEP_ITEMS
    start_date: datetime.datetime | None = None
    end_date: datetime.datetime | None = None
    lookback_days: int = DEFAULT_LOOKBACK_DAYS
    validate: bool = True
    band_workers: int = 8
    scene_workers: int = 4
    pansharpen_algorithm: PanAlgorithm | None = None
    pansharpen_pan_key: str = "visual"
    max_cloud_threshold: float = MAX_CLOUD_THRESHOLD
    min_confidence: float = MIN_CONFIDENCE_TO_SAVE
    save_report: bool = SAVE_VALIDATION_REPORT

    def resolved_end(self) -> datetime.datetime:
        return self.end_date or datetime.datetime.now(datetime.UTC)

    def resolved_start(self) -> datetime.datetime:
        return (
                self.start_date
                or (self.resolved_end() - datetime.timedelta(days=self.lookback_days))
        )

    def date_range_str(self) -> str:
        fmt = "%Y-%m-%dT%H:%M:%SZ"
        return (
            f"{self.resolved_start().strftime(fmt)}/"
            f"{self.resolved_end().strftime(fmt)}"
        )

    def band_keys(self) -> list[str]:
        return list(self.bands) if isinstance(self.bands, _BandGroup) else self.bands

    def tech_keys(self) -> list[str]:
        if self.tech_bands is None:
            return []
        return (
            list(self.tech_bands)
            if isinstance(self.tech_bands, _BandGroup)
            else self.tech_bands
        )

    def subdir(self, kind: str) -> str:
        path = os.path.join(self.output_dir, kind)
        os.makedirs(path, exist_ok=True)
        return path


# Geometry helpers

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


def _clip(
        da: xr.DataArray,
        lon: float,
        lat: float,
        half: float,
        reference_da: xr.DataArray | None = None,
) -> xr.DataArray:
    clipped = da.rio.clip_box(
        minx=lon - half, miny=lat - half,
        maxx=lon + half, maxy=lat + half,
        crs="EPSG:4326",
        allow_one_dimensional_raster=True,
    )
    clipped = clipped.squeeze().drop_vars(["band"], errors="ignore")

    if reference_da is None:
        return clipped
    if clipped.rio.crs is None or reference_da.rio.crs is None:
        return clipped

    if clipped.rio.crs == reference_da.rio.crs:
        src_affine = _extract_affine(clipped)
        dst_affine = _extract_affine(reference_da)
        if src_affine is not None and dst_affine is not None:
            dst_rows = reference_da.shape[-2]
            dst_cols = reference_da.shape[-1]
            src_np = np.asarray(clipped.values, dtype=np.float64)
            if src_np.ndim == 3 and src_np.shape[0] == 1:
                src_np = src_np[0]
            if src_np.ndim == 2:
                out_np = _ft_reproject_nearest(
                    src_np, src_affine, dst_rows, dst_cols, dst_affine
                )
                return (
                    xr.DataArray(
                        out_np,
                        dims=["y", "x"],
                        coords={
                            "y": reference_da.coords["y"],
                            "x": reference_da.coords["x"],
                        },
                        attrs=clipped.attrs,
                    ).rio.write_crs(clipped.rio.crs)
                )

    return clipped.rio.reproject_match(reference_da)


# I/O helpers

def _save_both(
        da_or_ds: xr.DataArray | xr.Dataset,
        directory: str,
        base: str,
) -> list[str]:
    written: list[str] = []
    tif_path = os.path.join(directory, base + ".tif")
    raster = (
        da_or_ds.to_array(dim="band")
        if isinstance(da_or_ds, xr.Dataset)
        else da_or_ds
    )
    raster.rio.to_raster(tif_path)
    written.append(tif_path)

    nc_path = os.path.join(directory, base + ".nc")
    try:
        nc_obj = da_or_ds
        try:
            crs = (
                da_or_ds.rio.crs
                if isinstance(da_or_ds, xr.DataArray)
                else next(iter(da_or_ds.data_vars.values())).rio.crs
            )
            if crs is not None:
                if isinstance(nc_obj, xr.DataArray):
                    nc_obj = nc_obj.rio.write_crs(crs, grid_mapping_name="spatial_ref")
                else:
                    first_var = next(iter(nc_obj.data_vars))
                    nc_obj = nc_obj.copy()
                    nc_obj[first_var] = nc_obj[first_var].rio.write_crs(
                        crs, grid_mapping_name="spatial_ref"
                    )
        except Exception as crs_exc:
            logger.debug(f"[sentinel] CRS encoding skipped for {nc_path}: {crs_exc}")
        nc_obj.to_netcdf(nc_path)
        written.append(nc_path)
    except ValueError as exc:
        if "backend" in str(exc).lower() or "netcdf" in str(exc).lower():
            logger.debug(
                f"[sentinel] NetCDF backend not available — skipping {nc_path}"
            )
        else:
            raise
    return written


def _write_report(report: dict, directory: str, base_name: str) -> None:
    path = os.path.join(directory, f"{base_name}_report.json")
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    logger.debug(f"[validation] Report → {path}")


# Validation

def _run_validation(
        item,
        lon: float,
        lat: float,
        half: float,
        reference_da: xr.DataArray | None,
        cfg: DownloadConfig,
) -> tuple[bool, dict, xr.DataArray | None]:
    report: dict = {"item_id": item.id, "passed": False}

    scl_asset = item.assets.get("scl")
    if not scl_asset:
        logger.warning(
            f"[validation] No SCL asset for {item.id} — skipping validation"
        )
        report["warning"] = "SCL asset missing; validation skipped"
        report["passed"] = True
        return True, report, None

    stac_shape = (
            scl_asset.extra_fields.get("proj:shape")
            or item.properties.get("proj:shape")
    )
    if stac_shape and len(stac_shape) >= 2:
        pre_rows, pre_cols = int(stac_shape[0]), int(stac_shape[1])
        dim_result = call_check_dimensions(pre_rows, pre_cols)
        if not dim_result["passed"]:
            report.update({
                "rows": pre_rows, "cols": pre_cols,
                "dimension_pass": False,
                "issues": dim_result["issues"],
            })
            logger.info(
                f"[validation] {item.id} rejected by dimension pre-check "
                f"({pre_rows}×{pre_cols}): {dim_result['issues']}"
            )
            return False, report, None

    try:
        da = rioxarray.open_rasterio(scl_asset.href, lock=False)
        clipped = _clip(da, lon, lat, half, reference_da)
    except Exception as exc:
        logger.warning(f"[validation] SCL read failed for {item.id}: {exc}")
        report["error"] = str(exc)
        report["passed"] = True
        return True, report, None

    clipped_shape = clipped.values.shape
    rows = clipped_shape[-2] if clipped.values.ndim >= 2 else 1
    cols = clipped_shape[-1] if clipped.values.ndim >= 1 else 1

    dim_result = call_check_dimensions(rows, cols)
    report.update({
        "rows": rows, "cols": cols,
        "dimension_pass": dim_result["passed"],
    })
    if not dim_result["passed"]:
        report["issues"] = dim_result["issues"]
        report["passed"] = False
        logger.info(
            f"[validation] {item.id} rejected by dimension check "
            f"({rows}×{cols}): {dim_result['issues']}"
        )
        return False, report, None

    scl_np = np.asarray(clipped.values, dtype=np.int32).ravel()
    scl_result = call_validate_scl(scl_np, cfg.max_cloud_threshold)
    report.update(scl_result)

    radio_pass = True
    try:
        radio_pass = call_check_radiometry(scl_np.astype(np.float64, copy=False))
    except Exception as exc:
        logger.warning(
            f"[validation] Radiometry check failed for {item.id}: {exc}"
        )
        report["radiometry_error"] = str(exc)

    report["radiometry_pass"] = radio_pass
    passes = scl_result["confidence_score"] >= cfg.min_confidence and radio_pass
    report["passed"] = passes
    return passes, report, clipped


# Pansharpening

def _apply_pansharpening(
        ds: xr.DataArray,
        item,
        lon: float,
        lat: float,
        half: float,
        cfg: DownloadConfig,
) -> xr.DataArray:
    pan_asset = item.assets.get(cfg.pansharpen_pan_key)
    if not pan_asset:
        logger.warning(
            f"[pansharpen] PAN asset '{cfg.pansharpen_pan_key}' not found "
            f"in {item.id} — skipping pansharpening"
        )
        return ds

    try:
        pan_da = rioxarray.open_rasterio(pan_asset.href, lock=False)
        pan_clipped = _clip(pan_da, lon, lat, half)
    except Exception as exc:
        logger.warning(f"[pansharpen] Failed to load PAN for {item.id}: {exc}")
        return ds

    try:
        pan_np = np.asarray(pan_clipped.values, dtype=np.float64)
        if pan_np.ndim == 3:
            pan_np = (
                pan_np[0]
                if pan_np.shape[0] == 1
                else _ft_rgb_to_luminance(pan_np)
            )

        ms_np = np.asarray(ds.values, dtype=np.float64)
        if ms_np.ndim == 2:
            ms_np = ms_np[np.newaxis]

        sharpened_np = pansharpen(pan_np, ms_np, algorithm=cfg.pansharpen_algorithm)
    except Exception as exc:
        logger.warning(f"[pansharpen] Fortran call failed for {item.id}: {exc}")
        return ds

    pan_rows, pan_cols = pan_np.shape
    sharpened_da = xr.DataArray(
        sharpened_np,
        dims=["band", "y", "x"],
        coords={
            "band": ds.coords["band"],
            "y": (
                pan_clipped.coords["y"]
                if "y" in pan_clipped.coords
                else np.arange(pan_rows)
            ),
            "x": (
                pan_clipped.coords["x"]
                if "x" in pan_clipped.coords
                else np.arange(pan_cols)
            ),
        },
    )
    if pan_clipped.rio.crs is not None:
        sharpened_da = (
            sharpened_da
            .rio.set_spatial_dims(x_dim="x", y_dim="y")
            .rio.write_crs(pan_clipped.rio.crs)
        )
    logger.debug(
        f"[pansharpen] {item.id} sharpened "
        f"{ms_np.shape[1]}×{ms_np.shape[2]} → {pan_rows}×{pan_cols} "
        f"via {cfg.pansharpen_algorithm}"
    )
    return sharpened_da


# Per-band fetch

def _fetch_band(
        key: str,
        href: str,
        lon: float,
        lat: float,
        half: float,
) -> tuple[str, xr.DataArray]:
    da = rioxarray.open_rasterio(href, lock=False)
    return key, _clip(da, lon, lat, half)


# Core scene download

def _download_item(
        item,
        lon: float,
        lat: float,
        base_name: str,
        cfg: DownloadConfig,
) -> list[str]:
    half = cfg.bbox_half_deg
    written: list[str] = []

    band_tasks: dict[str, str] = {}
    for key in cfg.band_keys():
        asset = item.assets.get(key)
        if asset:
            band_tasks[key] = asset.href
        else:
            logger.warning(f"[sentinel] Missing band '{key}' in {item.id}")

    with ThreadPoolExecutor(max_workers=cfg.band_workers) as pool:

        val_future: Future | None = None
        if cfg.validate:
            val_future = pool.submit(
                _run_validation, item, lon, lat, half, None, cfg
            )

        band_futures: dict[Future, str] = {
            pool.submit(_fetch_band, key, href, lon, lat, half): key
            for key, href in band_tasks.items()
        }

        passes = True
        report: dict = {"item_id": item.id, "passed": True}
        cached_scl: xr.DataArray | None = None

        if val_future is not None:
            try:
                passes, report, cached_scl = val_future.result()
            except Exception as exc:
                logger.warning(
                    f"[validation] Unexpected error for {item.id}: {exc}"
                )
                passes = True

            if not passes:
                for f in band_futures:
                    f.cancel()

                _conf = report.get("confidence_score")
                _cloud = report.get("cloud_ratio")
                _conf_s = f"{_conf:.2f}" if isinstance(_conf, float) else "?"
                _cloud_s = f"{_cloud:.2f}" if isinstance(_cloud, float) else "?"
                logger.info(
                    f"[validation] {base_name} rejected "
                    f"(confidence={_conf_s}, "
                    f"cloud={_cloud_s}, "
                    f"rows={report.get('rows', '?')}, "
                    f"cols={report.get('cols', '?')}, "
                    f"issues={report.get('issues', [])})"
                )
                if cfg.save_report:
                    _write_report(report, cfg.subdir("spectral"), base_name)
                return []

        band_results: dict[str, xr.DataArray] = {}
        for future in as_completed(band_futures):
            key = band_futures[future]
            try:
                _, clipped = future.result()
                band_results[key] = clipped
            except Exception as exc:
                logger.warning(f"[sentinel] Band '{key}' failed: {exc}")

        ok_keys = [k for k in cfg.band_keys() if k in band_results]
        arrays = [band_results[k] for k in ok_keys]

        if not arrays:
            logger.warning(
                f"[sentinel] No spectral bands for {base_name} — item skipped."
            )
            return []

        reference_da = arrays[0]

        if reference_da.rio.crs is not None:
            dst_rows = reference_da.shape[-2]
            dst_cols = reference_da.shape[-1]
            dst_affine = _extract_affine(reference_da)
            ref_crs = reference_da.rio.crs
            aligned: list[xr.DataArray] = []

            for da in arrays:
                if da is reference_da:
                    aligned.append(da)
                    continue

                if (
                        dst_affine is not None
                        and da.rio.crs is not None
                        and da.rio.crs == ref_crs
                ):
                    src_affine = _extract_affine(da)
                    src_np = np.asarray(da.values, dtype=np.float64)
                    if src_np.ndim == 3 and src_np.shape[0] == 1:
                        src_np = src_np[0]
                    if src_affine is not None and src_np.ndim == 2:
                        out_np = _ft_reproject_nearest(
                            src_np, src_affine, dst_rows, dst_cols, dst_affine
                        )
                        aligned.append(
                            xr.DataArray(
                                out_np,
                                dims=["y", "x"],
                                coords={
                                    "y": reference_da.coords["y"],
                                    "x": reference_da.coords["x"],
                                },
                                attrs=da.attrs,
                            ).rio.write_crs(ref_crs)
                        )
                        continue

                if da.rio.crs is not None:
                    aligned.append(da.rio.reproject_match(reference_da))
                else:
                    aligned.append(
                        da.interp(
                            x=reference_da.coords["x"],
                            y=reference_da.coords["y"],
                            method="nearest",
                        )
                    )
            arrays = aligned

        ds = xr.concat(arrays, dim="band").assign_coords(band=ok_keys)
        if reference_da.rio.crs:
            ds = (
                ds.rio.set_spatial_dims(x_dim="x", y_dim="y")
                .rio.write_crs(reference_da.rio.crs)
            )

        if cfg.pansharpen_algorithm is not None:
            ds = _apply_pansharpening(ds, item, lon, lat, half, cfg)

        written.extend(_save_both(ds, cfg.subdir("spectral"), base_name))
        if cfg.save_report:
            _write_report(report, cfg.subdir("spectral"), base_name)

        # Technical layers
        def _fetch_tech(key: str, href: str) -> tuple[str, list[str]]:
            if key == "scl" and cached_scl is not None:
                paths = _save_both(
                    cached_scl, cfg.subdir("technical"), f"{key}_{base_name}"
                )
                return key, paths
            da = rioxarray.open_rasterio(href, lock=False)
            layer = _clip(da, lon, lat, half, reference_da)
            return key, _save_both(layer, cfg.subdir("technical"), f"{key}_{base_name}")

        tech_tasks = {
            key: item.assets[key].href
            for key in cfg.tech_keys()
            if key in item.assets
        }
        tech_futures = {
            pool.submit(_fetch_tech, k, h): k
            for k, h in tech_tasks.items()
        }
        for future in as_completed(tech_futures):
            key = tech_futures[future]
            try:
                _, paths = future.result()
                written.extend(paths)
            except Exception as exc:
                logger.warning(f"[sentinel] Layer '{key}' failed: {exc}")

        # Visual
        if cfg.visual:
            def _fetch_visual(vis_key: str, href: str) -> list[str]:
                da = rioxarray.open_rasterio(href, lock=False)
                vis = _clip(da, lon, lat, half)
                return _save_both(vis, cfg.subdir("visual"), f"vis_{base_name}")

            vis_tasks = {
                vis_key: item.assets[vis_key].href
                for vis_key in list(VisualAssets.VISUAL)
                if vis_key in item.assets
            }
            vis_futures = {
                pool.submit(_fetch_visual, k, h): k
                for k, h in vis_tasks.items()
            }
            for future in as_completed(vis_futures):
                try:
                    written.extend(future.result())
                except Exception as exc:
                    logger.warning(f"[sentinel] Visual failed: {exc}")

    return written


# Public entry point

def download_sentinel2(
        locations: Sequence[LocationSpec],
        cfg: DownloadConfig | None = None,
        progress: bool = True,
) -> dict[str, list[str]]:
    if cfg is None:
        cfg = DownloadConfig()

    client = Client.open(STAC_API_URL)
    date_range = cfg.date_range_str()
    results: dict[str, list[str]] = {}

    all_items: list[tuple[object, LocationSpec]] = []
    for loc in locations:
        half = cfg.bbox_half_deg
        bbox = [loc.lon - half, loc.lat - half, loc.lon + half, loc.lat + half]
        search = client.search(
            collections=["sentinel-2-l2a"],
            bbox=bbox,
            datetime=date_range,
            max_items=cfg.max_items,
            sortby=[{"field": "properties.datetime", "direction": "desc"}],
        )
        items = list(search.items())
        if not items:
            logger.warning(f"[sentinel] No items for ({loc.lat}, {loc.lon})")
            continue
        items = sorted(
            items,
            key=lambda x: (
                    x.datetime or datetime.datetime.min.replace(tzinfo=datetime.UTC)
            ),
            reverse=True,
        )[: cfg.keep_items]
        for item in items:
            all_items.append((item, loc))

    total = len(all_items)
    if total == 0:
        return results

    _progress_start(total, progress)

    scene_tasks: list[tuple[object, LocationSpec, str]] = []
    for item, loc in all_items:
        timestamp = (
                item.datetime
                or datetime.datetime.min.replace(tzinfo=datetime.UTC)
        )
        base_name = (
            f"{loc.name}_{timestamp.strftime('%Y%m%dT%H%M%S')}"
            if loc.name
            else loc.base_name(timestamp)
        )
        scene_tasks.append((item, loc, base_name))

    completed_count = 0
    import threading
    lock = threading.Lock()

    def _process_scene(args: tuple) -> tuple[str, list[str]]:
        nonlocal completed_count
        item, loc, base_name = args
        written = _download_item(item, loc.lon, loc.lat, base_name, cfg)
        with lock:
            completed_count += 1
            _progress_update(completed_count, total, base_name, progress)
        return base_name, written

    with ThreadPoolExecutor(max_workers=cfg.scene_workers) as pool:
        for base_name, written in pool.map(_process_scene, scene_tasks):
            results.setdefault(base_name, []).extend(written)

    _progress_done(total, results, progress)
    return results


def _progress_start(total: int, enabled: bool) -> None:
    if enabled:
        print(f"\n  Fetching {total} scene(s)\n")


def _progress_update(idx: int, total: int, name: str, enabled: bool) -> None:
    if not enabled:
        return
    bar_width = 30
    filled = int(bar_width * idx / total)
    bar = "█" * filled + "░" * (bar_width - filled)
    pct = int(100 * idx / total)
    label = name if len(name) <= 35 else name[:32] + "..."
    print(f"\r  [{bar}] {pct:3d}%  {idx}/{total}  {label:<35}", end="", flush=True)


def _progress_done(total: int, results: dict, enabled: bool) -> None:
    if not enabled:
        return
    saved = sum(len(v) for v in results.values())
    scenes = len(results)
    print(f"\r  {'█' * 30}  100%  {total}/{total}{'':40}")
    print(f"\n  Done — {scenes} scene(s), {saved} file(s) written\n")