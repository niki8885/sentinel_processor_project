from __future__ import annotations
import logging
import os
from pathlib import Path
from typing import Sequence
import numpy as np
import rioxarray
import xarray as xr

from sentinel_processor.indices._indices_bridge import (
    compute_arvi,
    compute_cig,
    compute_evi,
    compute_mndwi,
    compute_nbr,
    compute_ndbi,
    compute_ndsi,
    compute_ndvi,
    compute_ndwi,
    compute_savi,
)

logger = logging.getLogger(__name__)

IndexName = str
BandKey = str


# Band name aliases

_BAND_ALIASES: dict[str, str] = {
    # canonical (pass-through)
    "B02": "B02", "B03": "B03", "B04": "B04", "B05": "B05",
    "B06": "B06", "B07": "B07", "B08": "B08", "B8A": "B8A",
    "B11": "B11", "B12": "B12",
    # lowercase canonical
    "b02": "B02", "b03": "B03", "b04": "B04", "b05": "B05",
    "b06": "B06", "b07": "B07", "b08": "B08", "b8a": "B8A",
    "b11": "B11", "b12": "B12",
    # human-readable names
    "blue":      "B02",
    "green":     "B03",
    "red":       "B04",
    "rededge1":  "B05",
    "rededge2":  "B06",
    "rededge3":  "B07",
    "nir":       "B08",
    "nir08":     "B8A",
    "swir16":    "B11",
    "swir22":    "B12",
    # alternative spellings
    "nir_broad": "B08",
    "nir_narrow":"B8A",
    "swir1":     "B11",
    "swir2":     "B12",
    "red_edge1": "B05",
    "red_edge2": "B06",
    "red_edge3": "B07",
}

_INDEX_REGISTRY: dict[IndexName, dict] = {
    "ndvi": {
        "bands": ["B08", "B04"],
        "fn": lambda bands: compute_ndvi(bands["B08"], bands["B04"]),
        "long_name": "Normalised Difference Vegetation Index",
        "valid_range": (-1.0, 1.0),
    },
    "evi": {
        "bands": ["B08", "B04", "B02"],
        "fn": lambda bands: compute_evi(bands["B08"], bands["B04"], bands["B02"]),
        "long_name": "Enhanced Vegetation Index",
        "valid_range": (-1.0, 1.0),
    },
    "savi": {
        "bands": ["B08", "B04"],
        "fn": lambda bands: compute_savi(bands["B08"], bands["B04"]),
        "long_name": "Soil Adjusted Vegetation Index",
        "valid_range": (-1.5, 1.5),
    },
    "ndwi": {
        "bands": ["B03", "B08"],
        "fn": lambda bands: compute_ndwi(bands["B03"], bands["B08"]),
        "long_name": "Normalised Difference Water Index",
        "valid_range": (-1.0, 1.0),
    },
    "mndwi": {
        "bands": ["B03", "B11"],
        "fn": lambda bands: compute_mndwi(bands["B03"], bands["B11"]),
        "long_name": "Modified NDWI (Xu 2006)",
        "valid_range": (-1.0, 1.0),
    },
    "ndbi": {
        "bands": ["B11", "B08"],
        "fn": lambda bands: compute_ndbi(bands["B11"], bands["B08"]),
        "long_name": "Normalised Difference Built-up Index",
        "valid_range": (-1.0, 1.0),
    },
    "nbr": {
        "bands": ["B08", "B12"],
        "fn": lambda bands: compute_nbr(bands["B08"], bands["B12"]),
        "long_name": "Normalised Burn Ratio",
        "valid_range": (-1.0, 1.0),
    },
    "ndsi": {
        "bands": ["B03", "B11"],
        "fn": lambda bands: compute_ndsi(bands["B03"], bands["B11"]),
        "long_name": "Normalised Difference Snow Index",
        "valid_range": (-1.0, 1.0),
    },
    "cig": {
        "bands": ["B08", "B03"],
        "fn": lambda bands: compute_cig(bands["B08"], bands["B03"]),
        "long_name": "Chlorophyll Index Green",
        "valid_range": (0.0, 30.0),
    },
    "arvi": {
        "bands": ["B08", "B04", "B02"],
        "fn": lambda bands: compute_arvi(bands["B08"], bands["B04"], bands["B02"]),
        "long_name": "Atmospherically Resistant Vegetation Index",
        "valid_range": (-1.0, 1.0),
    },
}

AVAILABLE_INDICES: list[str] = sorted(_INDEX_REGISTRY)


def _resolve_band_coords(da: xr.DataArray) -> dict[str, str]:
    if "band" not in da.coords:
        return {}
    result = {}
    for raw in da.coords["band"].values:
        canonical = _BAND_ALIASES.get(str(raw))
        if canonical:
            result[canonical] = str(raw)
    return result


def _load_bands_from_multiband(
    da: xr.DataArray,
    needed: set[BandKey],
) -> dict[BandKey, np.ndarray] | None:
    if "band" not in da.coords:
        return None
    coord_map = _resolve_band_coords(da)
    if needed and not needed.issubset(coord_map.keys()):
        return None
    keys = needed if needed else set(coord_map.keys())
    return {
        key: da.sel(band=coord_map[key]).values.astype(np.float64)
        for key in keys
        if key in coord_map
    }


def _load_bands_from_dataset(
    ds: xr.Dataset,
    needed: set[BandKey],
) -> dict[BandKey, np.ndarray] | None:

    var_map = {v.upper(): v for v in ds.data_vars}
    result = {}
    for key in needed:
        if key in ds.data_vars:
            result[key] = ds[key].values.astype(np.float64)
        elif key.upper() in var_map:
            result[key] = ds[var_map[key.upper()]].values.astype(np.float64)
        else:
            return None
    return result


def _check_missing_bands(
    available_bands: set[BandKey],
    requested_indices: list[str],
) -> dict[str, list[BandKey]]:

    missing: dict[str, list[BandKey]] = {}
    for idx_name in requested_indices:
        entry = _INDEX_REGISTRY.get(idx_name)
        if entry is None:
            missing[idx_name] = ["<unknown index>"]
            continue
        m = [b for b in entry["bands"] if b not in available_bands]
        if m:
            missing[idx_name] = m
    return missing

def compute_indices(
    source: str | Path,
    indices: Sequence[str],
    output_dir: str | Path | None = None,
    scale_factor: float | None = None,
    nodata: float = -9999.0,
    overwrite: bool = True,
    output_format: str = "tif",
) -> dict[str, str]:

    source = Path(source)
    if not source.exists():
        raise FileNotFoundError(f"Source file not found: {source}")
    if not indices:
        raise ValueError("'indices' list must not be empty.")

    indices = [i.lower().strip() for i in indices]
    unknown = [i for i in indices if i not in _INDEX_REGISTRY]
    if unknown:
        raise ValueError(
            f"Unknown index name(s): {unknown}. "
            f"Available: {AVAILABLE_INDICES}"
        )

    output_format = output_format.lower().lstrip(".")
    if output_format not in ("tif", "tiff", "nc", "netcdf"):
        raise ValueError(f"Unsupported output_format '{output_format}'. Use 'tif' or 'nc'.")
    ext = "tif" if output_format in ("tif", "tiff") else "nc"

    if output_dir is None:
        output_dir = source.parent / "indices"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    suffix = source.suffix.lower()
    logger.info(f"[indices] Loading {source.name} …")

    bands_dict: dict[BandKey, np.ndarray] | None = None
    template_da: xr.DataArray | None = None

    if suffix in (".tif", ".tiff"):
        da = rioxarray.open_rasterio(source)
        bands_dict = _load_bands_from_multiband(da, set())
        needed_all: set[BandKey] = set()
        for idx_name in indices:
            needed_all |= set(_INDEX_REGISTRY[idx_name]["bands"])
        bands_dict = _load_bands_from_multiband(da, needed_all)
        if bands_dict is None:
            logger.warning(
                "[indices] No named 'band' coordinate found in the DataArray. "
                "Cannot auto-match bands. Ensure the file has a 'band' coordinate "
                "with Sentinel-2 band names (B02, B03, …)."
            )
            return {}
        template_da = da.isel(band=0).squeeze(drop=True)
        if "band" in template_da.dims:
            template_da = template_da.isel(band=0, drop=True)

    elif suffix in (".nc",):
        ds = xr.open_dataset(source)
        needed_all = set()
        for idx_name in indices:
            needed_all |= set(_INDEX_REGISTRY[idx_name]["bands"])

        bands_dict = _load_bands_from_dataset(ds, needed_all)
        if bands_dict is not None:
            first_key = next(iter(bands_dict))
            template_da = ds[first_key] if first_key in ds.data_vars else None
        else:
            for var in ds.data_vars:
                candidate = ds[var]
                if "band" in (candidate.coords or {}):
                    bands_dict = _load_bands_from_multiband(candidate, needed_all)
                    if bands_dict is not None:
                        tmpl = candidate.isel(band=0, drop=True).squeeze()
                        template_da = tmpl
                        try:
                            import rioxarray as _rio
                            crs = candidate.rio.crs
                            if crs is None and "spatial_ref" in ds:
                                crs = ds["spatial_ref"].attrs.get("crs_wkt") or ds["spatial_ref"].attrs.get("grid_mapping_name")
                            if crs is not None:
                                template_da = template_da.rio.write_crs(crs)
                        except Exception:
                            pass
                        break

        if bands_dict is None:
            logger.warning(
                "[indices] Could not extract named Sentinel-2 bands from NetCDF. "
                "Expected variable names or 'band' coordinate with aliases: "
                "B02/blue, B03/green, B04/red, B08/nir, B11/swir16, B12/swir22, …"
            )
            return {}
        ds.close()

    else:
        raise ValueError(f"Unsupported file format '{suffix}'. Use .tif, .tiff, or .nc.")

    if scale_factor is not None:
        logger.info(f"[indices] Applying scale_factor={scale_factor}")
        bands_dict = {k: v * scale_factor for k, v in bands_dict.items()}

    available_bands = set(bands_dict.keys())
    missing_map = _check_missing_bands(available_bands, indices)

    skipped = {k: v for k, v in missing_map.items() if v}
    computable = [i for i in indices if not missing_map.get(i)]

    if skipped:
        for idx_name, missing_bands in skipped.items():
            logger.warning(
                f"[indices] Skipping '{idx_name}': "
                f"missing band(s) {missing_bands} in source file."
            )

    if not computable:
        logger.warning("[indices] No indices could be computed — all required bands missing.")
        return {}

    base_stem = source.stem
    results: dict[str, str] = {}

    for idx_name in computable:
        entry = _INDEX_REGISTRY[idx_name]
        out_filename = f"indices_{base_stem}_{idx_name}.{ext}"
        out_path = output_dir / out_filename

        if not overwrite and out_path.exists():
            logger.info(f"[indices] {out_filename} exists, skipping (overwrite=False).")
            results[idx_name] = str(out_path)
            continue

        logger.info(f"[indices] Computing {idx_name.upper()} …")
        flat_result: np.ndarray = entry["fn"](bands_dict)

        if template_da is not None:
            spatial_shape = template_da.squeeze().shape
            if len(spatial_shape) > 2:
                spatial_shape = spatial_shape[-2:]
            try:
                result_2d = flat_result.reshape(spatial_shape)
            except ValueError:
                logger.error(
                    f"[indices] Cannot reshape {idx_name} result "
                    f"({flat_result.size} elements) to {spatial_shape}."
                )
                continue

            tmpl = template_da.squeeze()
            spatial_dim_names = [
                d for d in tmpl.dims
                if d in ("x", "y", "lat", "lon", "latitude", "longitude")
            ]
            if len(spatial_dim_names) < 2:
                spatial_dim_names = list(tmpl.dims[-2:])
            out_da = xr.DataArray(
                result_2d,
                coords={
                    k: tmpl.coords[k]
                    for k in tmpl.coords
                    if k in ("x", "y", "latitude", "longitude", "lat", "lon")
                },
                dims=spatial_dim_names,
                attrs={
                    "long_name": entry["long_name"],
                    "valid_range": list(entry["valid_range"]),
                    "_FillValue": nodata,
                    "source_file": str(source),
                },
            )
            try:
                crs = template_da.rio.crs
                if crs is not None:
                    x_dim = "x" if "x" in out_da.dims else spatial_dim_names[-1]
                    y_dim = "y" if "y" in out_da.dims else spatial_dim_names[-2]
                    out_da = out_da.rio.set_spatial_dims(
                        x_dim=x_dim, y_dim=y_dim, inplace=False
                    ).rio.write_crs(crs)
            except Exception as crs_err:
                logger.debug(f"[indices] CRS not written for {idx_name}: {crs_err}")
        else:
            out_da = xr.DataArray(flat_result, attrs={"long_name": entry["long_name"]})

        try:
            if ext == "tif":
                try:
                    out_da.rio.to_raster(str(out_path), dtype="float32")
                except Exception as raster_err:
                    import rasterio
                    from rasterio.transform import from_bounds
                    arr = out_da.values.astype("float32")
                    h, w = arr.shape
                    xs = template_da.coords.get("x")
                    ys = template_da.coords.get("y")
                    if xs is not None and ys is not None:
                        transform = from_bounds(
                            float(xs.min()), float(ys.min()),
                            float(xs.max()), float(ys.max()),
                            w, h,
                        )
                    else:
                        from rasterio.transform import Affine
                        transform = Affine.identity()
                    with rasterio.open(
                        str(out_path), "w",
                        driver="GTiff", height=h, width=w,
                        count=1, dtype="float32",
                        nodata=nodata,
                        transform=transform,
                    ) as dst:
                        dst.write(arr, 1)
                    logger.debug(f"[indices] Saved without CRS via rasterio fallback → {out_path}")
            else:
                out_da.to_netcdf(str(out_path))
            results[idx_name] = str(out_path)
            logger.info(f"[indices] Saved → {out_path}")
        except Exception as exc:
            logger.error(f"[indices] Failed to save {idx_name}: {exc}")

    return results


def list_indices() -> dict[str, dict]:
    return {
        name: {
            "bands_required": entry["bands"],
            "long_name": entry["long_name"],
        }
        for name, entry in _INDEX_REGISTRY.items()
    }