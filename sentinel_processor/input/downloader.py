from __future__ import annotations
import datetime
import json
import logging
import os
import warnings
from dataclasses import dataclass, field
from typing import Sequence, Union
import numpy as np
import rioxarray
import xarray as xr
from pystac_client import Client

from sentinel_processor.config import (
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
from sentinel_processor.validation._fortran_bridge import call_check_radiometry, call_validate_scl
from sentinel_processor.utils.data_utils import LocationSpec, SpectralBands, TechnicalLayers, VisualAssets, _BandGroup

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
    max_cloud_threshold: float = MAX_CLOUD_THRESHOLD
    min_confidence: float = MIN_CONFIDENCE_TO_SAVE
    save_report: bool = SAVE_VALIDATION_REPORT

    def resolved_end(self) -> datetime.datetime:
        return self.end_date or datetime.datetime.now(datetime.UTC)

    def resolved_start(self) -> datetime.datetime:
        return self.start_date or (self.resolved_end() - datetime.timedelta(days=self.lookback_days))

    def date_range_str(self) -> str:
        fmt = "%Y-%m-%dT%H:%M:%SZ"
        return f"{self.resolved_start().strftime(fmt)}/{self.resolved_end().strftime(fmt)}"

    def band_keys(self) -> list[str]:
        return list(self.bands) if isinstance(self.bands, _BandGroup) else self.bands

    def tech_keys(self) -> list[str]:
        if self.tech_bands is None:
            return []
        return list(self.tech_bands) if isinstance(self.tech_bands, _BandGroup) else self.tech_bands

    def subdir(self, kind: str) -> str:
        path = os.path.join(self.output_dir, kind)
        os.makedirs(path, exist_ok=True)
        return path


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
    if (
        reference_da is not None
        and clipped.rio.crs is not None
        and reference_da.rio.crs is not None
    ):
        clipped = clipped.rio.reproject_match(reference_da)
    return clipped.squeeze().drop_vars(["band", "spatial_ref"], errors="ignore")


def _save_both(da_or_ds: xr.DataArray | xr.Dataset, directory: str, base: str) -> list[str]:
    written: list[str] = []
    tif_path = os.path.join(directory, base + ".tif")
    raster = da_or_ds.to_array(dim="band") if isinstance(da_or_ds, xr.Dataset) else da_or_ds
    raster.rio.to_raster(tif_path)
    written.append(tif_path)
    nc_path = os.path.join(directory, base + ".nc")
    try:
        da_or_ds.to_netcdf(nc_path)
        written.append(nc_path)
    except ValueError as exc:
        if "backend" in str(exc).lower() or "netcdf" in str(exc).lower():
            logger.debug(f"[sentinel] NetCDF backend not available — skipping {nc_path}")
        else:
            raise
    return written


def _run_validation(
    item,
    lon: float,
    lat: float,
    half: float,
    reference_da: xr.DataArray | None,
    cfg: DownloadConfig,
) -> tuple[bool, dict]:
    report: dict = {"item_id": item.id, "passed": False}
    scl_asset = item.assets.get("scl")
    if not scl_asset:
        logger.warning(f"[validation] No SCL asset for {item.id} — skipping validation")
        report["warning"] = "SCL asset missing; validation skipped"
        report["passed"] = True
        return True, report
    try:
        da = rioxarray.open_rasterio(scl_asset.href)
        clipped = _clip(da, lon, lat, half, reference_da)
        scl_arr = clipped.values.flatten().astype(int).tolist()
    except Exception as exc:
        logger.warning(f"[validation] SCL read failed for {item.id}: {exc}")
        report["error"] = str(exc)
        report["passed"] = True
        return True, report
    scl_result = call_validate_scl(scl_arr, cfg.max_cloud_threshold)
    report.update(scl_result)
    radio_pass = True
    try:
        pixel_vals = clipped.values.flatten().astype(float)
        radio_pass = call_check_radiometry(pixel_vals)
    except Exception as exc:
        logger.warning(f"[validation] Radiometry check failed for {item.id}: {exc}")
        report["radiometry_error"] = str(exc)
    report["radiometry_pass"] = radio_pass
    passes = scl_result["confidence_score"] >= cfg.min_confidence and radio_pass
    report["passed"] = passes
    return passes, report


def _download_item(
    item,
    lon: float,
    lat: float,
    base_name: str,
    cfg: DownloadConfig,
) -> list[str]:
    half = cfg.bbox_half_deg
    written: list[str] = []
    arrays: list[xr.DataArray] = []
    ok_keys: list[str] = []
    reference_da = None

    for key in cfg.band_keys():
        asset = item.assets.get(key)
        if not asset:
            logger.warning(f"[sentinel] Missing band '{key}' in {item.id}")
            continue
        try:
            da = rioxarray.open_rasterio(asset.href)
            clipped = _clip(da, lon, lat, half)
            if reference_da is None:
                reference_da = clipped
            arrays.append(clipped)
            ok_keys.append(key)
        except Exception as exc:
            logger.warning(f"[sentinel] Band '{key}' failed: {exc}")

    if not arrays:
        logger.warning(f"[sentinel] No spectral bands for {base_name} — item skipped.")
        return []

    if reference_da is not None and reference_da.rio.crs is not None:
        aligned = []
        for da in arrays:
            if da is reference_da:
                aligned.append(da)
            elif da.rio.crs is not None:
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
    if reference_da is not None and reference_da.rio.crs:
        ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs(reference_da.rio.crs)

    report: dict = {"item_id": item.id, "passed": True}
    if cfg.validate:
        passes, report = _run_validation(item, lon, lat, half, reference_da, cfg)
        if not passes:
            logger.info(
                f"[validation] {base_name} rejected "
                f"(confidence={report.get('confidence_score', '?'):.2f}, "
                f"cloud={report.get('cloud_ratio', '?'):.2f}, "
                f"issues={report.get('issues', [])})"
            )
            if cfg.save_report:
                _write_report(report, cfg.subdir("spectral"), base_name)
            return []

    written.extend(_save_both(ds, cfg.subdir("spectral"), base_name))
    if cfg.save_report:
        _write_report(report, cfg.subdir("spectral"), base_name)

    for key in cfg.tech_keys():
        asset = item.assets.get(key)
        if not asset:
            continue
        try:
            da = rioxarray.open_rasterio(asset.href)
            layer = _clip(da, lon, lat, half, reference_da)
            written.extend(_save_both(layer, cfg.subdir("technical"), f"{key}_{base_name}"))
        except Exception as exc:
            logger.warning(f"[sentinel] Layer '{key}' failed: {exc}")

    if cfg.visual:
        for vis_key in list(VisualAssets.VISUAL):
            asset = item.assets.get(vis_key)
            if not asset:
                continue
            try:
                da = rioxarray.open_rasterio(asset.href)
                vis = _clip(da, lon, lat, half)
                written.extend(_save_both(vis, cfg.subdir("visual"), f"vis_{base_name}"))
            except Exception as exc:
                logger.warning(f"[sentinel] Visual failed: {exc}")

    return written


def _write_report(report: dict, directory: str, base_name: str) -> None:
    path = os.path.join(directory, f"{base_name}_report.json")
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    logger.debug(f"[validation] Report → {path}")


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
            key=lambda x: x.datetime or datetime.datetime.min.replace(tzinfo=datetime.UTC),
            reverse=True,
        )[: cfg.keep_items]
        for item in items:
            all_items.append((item, loc))

    total = len(all_items)
    if total == 0:
        return results

    _progress_start(total, progress)

    for idx, (item, loc) in enumerate(all_items, 1):
        timestamp = item.datetime or datetime.datetime.min.replace(tzinfo=datetime.UTC)
        if loc.name:
            base_name = f"{loc.name}_{timestamp.strftime('%Y%m%dT%H%M%S')}"
        else:
            base_name = loc.base_name(timestamp)

        _progress_update(idx, total, base_name, progress)

        written = _download_item(item, loc.lon, loc.lat, base_name, cfg)
        results.setdefault(base_name, []).extend(written)

    _progress_done(total, results, progress)
    return results


def _progress_start(total: int, enabled: bool) -> None:
    if not enabled:
        return
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