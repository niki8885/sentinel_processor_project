from __future__ import annotations
import argparse
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore", category=RuntimeWarning)

try:
    import sentinel_processor as sp
    from sentinel_processor import stack_timeseries, TimeSeriesConfig
    from sentinel_processor.processing._timeseries_bridge import interpolate_gaps
except ImportError as e:
    print(f"[ERROR] sentinel_processor not found: {e}")
    print("  Run: pip install -e .  from the project root")
    sys.exit(1)

try:
    from sentinel_processor.analysis._sentinel_stats_bridge import (
        NODATA,
        valid_obs_count,
        temporal_gap_stats,
        pixel_quantiles,
        pixel_iqr,
        time_window_stats,
        anomaly_zscore,
        mann_kendall,
        trend_theil_sen,
        bfast_breakpoint,
        pixel_regression,
        phenology_metrics,
        save_phenology,
        pearson_map,
    )
except ImportError as e:
    print(f"[ERROR] _sentinel_stats_bridge not found: {e}")
    print("  Ensure the file is in sentinel_processor/analysis/")
    print("  and libsentinel_stats is compiled: python build_sentinel_stats.py")
    sys.exit(1)

try:
    import rasterio
    from rasterio.transform import from_bounds

    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False
    print("[WARN] rasterio not installed — GeoTIFF writing skipped")

OUT_DIR = Path("analysis_output")
OUT_DIR.mkdir(exist_ok=True)


def _sep(title: str) -> None:
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print('=' * 60)


def _stats_line(arr: np.ndarray, name: str, unit: str = "") -> None:
    valid = arr[arr != NODATA]
    valid = valid[np.isfinite(valid)]
    if valid.size == 0:
        print(f"  {name:<28}  (all NODATA)")
        return
    print(f"  {name:<28}  min={valid.min():.4f}  "
          f"mean={valid.mean():.4f}  max={valid.max():.4f}  {unit}")


def _save_tif(arr: np.ndarray, path: Path,
              transform=None, crs: str = "EPSG:32633") -> None:
    """Write a (rows, cols) float64 array to a single-band GeoTIFF."""
    if not HAS_RASTERIO:
        return
    rows, cols = arr.shape
    t = transform or from_bounds(0, 0, cols, rows, cols, rows)
    with rasterio.open(
            path, "w",
            driver="GTiff", dtype="float64",
            width=cols, height=rows, count=1,
            nodata=NODATA, transform=t, crs=crs,
            compress="deflate",
    ) as dst:
        dst.write(arr.astype(np.float64), 1)


def _ndvi(arr64: np.ndarray, band_names: list[str]) -> np.ndarray | None:
    """Compute NDVI = (NIR-RED)/(NIR+RED). Returns None if bands are missing."""
    try:
        ri = band_names.index("red")
        ni = band_names.index("nir")
    except ValueError:
        return None
    red = arr64[:, ri]
    nir = arr64[:, ni]
    with np.errstate(invalid="ignore", divide="ignore"):
        ndvi = np.where(
            (nir + red) > 0,
            (nir - red) / (nir + red),
            NODATA,
        )
    ndvi[~np.isfinite(ndvi)] = NODATA
    return ndvi


def _quick_crop(arr: np.ndarray, size: int = 64) -> np.ndarray:
    """Return a centre crop of size x size pixels for quick testing."""
    _, _, rows, cols = arr.shape
    r0 = max(0, (rows - size) // 2)
    c0 = max(0, (cols - size) // 2)
    return arr[:, :, r0:r0 + size, c0:c0 + size].copy()


# Step 1 — Download and stack

def step_download_and_stack(skip_dl: bool) -> tuple[xr.DataArray, list[datetime]]:
    _sep("Step 1 — Download and stack")

    if not skip_dl:
        cfg = sp.DownloadConfig(
            bands=sp.SpectralBands.ALL,
            tech_bands=sp.TechnicalLayers.SCL,
            visual=True,
            keep_items=20,
            max_cloud_threshold=0.80,
            min_confidence=0.0,
            pansharpen_algorithm="gram_schmidt",
            save_report=True,
        )
        sp.download_sentinel2(
            [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
            cfg=cfg,
        )
        print("  Download complete")
    else:
        print("  --skip-dl: skipping download, using cached data")

    result = stack_timeseries(
        sources=sorted(Path("data/spectral").glob("budapest_*.nc")),
        scl_dir="data/technical",
        cfg=TimeSeriesConfig(
            max_cloud_fraction=0.80,
            min_confidence=0.0,
            require_bands=["red", "nir"],
            fill_rejected=True,
            nodata=float("nan"),
            save_dir="data/stacks",
            save_format="nc",
            save_name="budapest_cloudy",
        ),
        progress=True,
    )

    print(result.summary())
    result.to_json("stack_log.json")
    da = result.stack

    # Extract acquisition datetimes from xarray time coordinate
    dates: list[datetime] = [
        datetime.utcfromtimestamp(int(t) / 1e9)
        for t in da.coords["time"].values.astype("int64")
    ]

    print(f"\n  Stack shape : {tuple(da.shape)}  (time, band, y, x)")
    print(f"  Scenes      : {len(dates)}")
    print(f"  Date range  : {dates[0].date()} → {dates[-1].date()}")
    return da, dates


# Step 2 — Gap-fill PCHIP

def step_gapfill(da: xr.DataArray) -> tuple[xr.DataArray, np.ndarray, list[str]]:
    _sep("Step 2 — Gap-fill (PCHIP)")
    arr = da.values.astype(np.float64)
    arr[~np.isfinite(arr)] = NODATA
    mask = (arr != NODATA).astype(np.int32)

    filled = arr.copy()
    for b in range(arr.shape[1]):
        filled[:, b] = interpolate_gaps(arr[:, b], mask[:, b], method="pchip")

    da_filled = xr.DataArray(
        filled.astype(np.float32),
        dims=da.dims, coords=da.coords, attrs=da.attrs,
    )
    da_filled = da_filled.where(da_filled != NODATA)

    band_names = (
        [str(b) for b in da.coords["band"].values]
        if "band" in da.coords
        else [f"b{i}" for i in range(da.shape[1])]
    )

    gaps_before = (arr == NODATA).sum()
    gaps_after = np.isnan(da_filled.values).sum()
    print(f"  Gaps before : {gaps_before:,}")
    print(f"  Gaps after  : {gaps_after:,}  (pixels with no valid observation)")
    print(f"  Bands       : {band_names}")

    return da_filled, filled, band_names


# Step 3 — Analysis via sentinel_stats

def step_analysis(
        filled64: np.ndarray,
        dates: list[datetime],
        band_names: list[str],
        quick: bool,
) -> None:
    """Run all 14 statistics module functions against the gap-filled stack."""

    if quick:
        filled64 = _quick_crop(filled64)
        print(f"  [--quick] Crop: {filled64.shape}")

    n_times, n_bands, rows, cols = filled64.shape
    days = np.array(
        [(d - dates[0]).total_seconds() / 86_400.0 for d in dates],
        dtype=np.float64,
    )

    # Replace NaN with NODATA sentinel before passing arrays to Fortran
    f64 = filled64.copy()
    f64[~np.isfinite(f64)] = NODATA

    # 3.1  valid_obs_count
    _sep("3.1  valid_obs_count")
    for b, bname in enumerate(band_names):
        oc = valid_obs_count(f64[:, b])
        _stats_line(oc["fraction"], f"{bname} fraction_valid")
        _save_tif(oc["count"].astype(np.float64), OUT_DIR / "coverage.tif")
    # coverage.tif from the first band is sufficient — geometry is shared

    # 3.2  temporal_gap_stats
    _sep("3.2  temporal_gap_stats")
    gs = temporal_gap_stats(f64[:, 0], dates)
    _stats_line(gs["max_gap"], "max_gap  [days]")
    _stats_line(gs["mean_gap"], "mean_gap [days]")
    _save_tif(gs["max_gap"], OUT_DIR / "max_gap.tif")

    pct_risky = (gs["max_gap"][gs["max_gap"] != NODATA] > 60).mean() * 100
    print(f"  Pixels with gap > 60 days : {pct_risky:.1f} %")

    # 3.3  pixel_quantiles
    _sep("3.3  pixel_quantiles")
    for b, bname in enumerate(band_names[:3]):  # first three bands
        q = pixel_quantiles(f64[:, b])
        amplitude = q["p90"] - q["p10"]
        _stats_line(q["p50"], f"{bname} median")
        _stats_line(amplitude, f"{bname} p90-p10 amplitude")
        _save_tif(q["p50"], OUT_DIR / f"p50_{bname}.tif")

    # 3.4  pixel_iqr
    _sep("3.4  pixel_iqr  (outlier mask)")
    for b, bname in enumerate(band_names[:3]):
        iq = pixel_iqr(f64[:, b])
        n_outliers = (iq["outlier"] == 1.0).sum()
        _stats_line(iq["iqr"], f"{bname} IQR")
        print(f"  {bname}: {n_outliers:,} outlier pixels across all scenes")
        _save_tif(iq["iqr"], OUT_DIR / f"iqr_{bname}.tif")

    # 3.5  time_window_stats + anomaly_zscore
    _sep("3.5  time_window_stats  +  anomaly_zscore")
    for b, bname in enumerate(band_names[:3]):
        ws = time_window_stats(f64[:, b], dates, window_days=90)
        _stats_line(ws["slope"], f"{bname} OLS slope [/day]")
        _stats_line(ws["std"], f"{bname} std")

        # z-score against the last scene as a sample
        z = anomaly_zscore(f64[:, b], dates, background=ws)
        last_z = z[-1]
        n_sig = (np.abs(last_z[last_z != NODATA]) > 2.0).sum()
        print(f"  {bname}: last-scene anomalies |z|>2.0 : {n_sig:,} pixels")

    # 3.6  mann_kendall
    _sep("3.6  mann_kendall")
    try:
        ri = band_names.index("nir")
    except ValueError:
        ri = 0
    bname = band_names[ri]
    mk = mann_kendall(f64[:, ri], dates)
    n_inc = (mk["trend"] == 1.0).sum()
    n_dec = (mk["trend"] == -1.0).sum()
    n_ns = (mk["trend"] == 0.0).sum()
    print(f"  {bname}  Increasing: {n_inc:,}  Decreasing: {n_dec:,}  No trend: {n_ns:,}")
    _stats_line(mk["Z"], f"{bname} Z-statistic")
    _save_tif(mk["trend"], OUT_DIR / f"mk_trend_{bname}.tif")

    # 3.7  trend_theil_sen
    _sep("3.7  trend_theil_sen")
    ts = trend_theil_sen(f64[:, ri], dates)
    _stats_line(ts["slope"], f"{bname} Theil-Sen slope [/day]")
    _stats_line(ts["intercept"], f"{bname} intercept")
    _save_tif(ts["slope"], OUT_DIR / f"ts_slope_{bname}.tif")

    # 3.8  bfast_breakpoint
    _sep("3.8  bfast_breakpoint")
    bp = bfast_breakpoint(f64[:, ri], dates)
    valid_bp = bp["rss_ratio"][bp["rss_ratio"] != NODATA]
    n_breaks = (valid_bp < 0.7).sum()
    print(f"  {bname}: pixels with structural break (rss_ratio<0.7): {n_breaks:,}")
    _stats_line(bp["magnitude"], f"{bname} break magnitude")

    # 3.9  pixel_regression
    _sep("3.9  pixel_regression  (slope / intercept / R²)")
    for b, bname in enumerate(band_names[:3]):
        reg = pixel_regression(f64[:, b], days)
        _stats_line(reg["slope"], f"{bname} OLS slope [/day]")
        _stats_line(reg["r2"], f"{bname} R²")
        # fraction of pixels where the linear model is meaningful
        good = (reg["r2"][reg["r2"] != NODATA] > 0.5).mean() * 100
        print(f"  {bname}: R2>0.5 : {good:.1f} % of pixels")
        _save_tif(reg["slope"], OUT_DIR / f"slope_{bname}.tif")
        _save_tif(reg["r2"], OUT_DIR / f"r2_{bname}.tif")

    # 3.10  phenology_metrics + save_phenology
    _sep("3.10  phenology_metrics  (NDVI → SOS / EOS / peak)")
    ndvi = _ndvi(f64, band_names)
    if ndvi is not None:
        pm = phenology_metrics(ndvi, dates, smooth=True, savgol_window=5)
        _stats_line(pm["sos_doy"], "SOS [days from start]")
        _stats_line(pm["eos_doy"], "EOS [days from start]")
        _stats_line(pm["peak_doy"], "peak_doy")
        _stats_line(pm["peak_val"], "peak NDVI")

        # Convert day offsets to calendar dates for interpretation
        t0 = dates[0]
        sos_valid = pm["sos_doy"][pm["sos_doy"] != NODATA]
        if sos_valid.size > 0:
            med_sos = float(np.median(sos_valid))
            sos_date = t0 + timedelta(days=med_sos)
            print(f"\n  Median SOS : ~{sos_date.strftime('%d %b')}")
        eos_valid = pm["eos_doy"][pm["eos_doy"] != NODATA]
        if eos_valid.size > 0:
            med_eos = float(np.median(eos_valid))
            eos_date = t0 + timedelta(days=med_eos)
            print(f"  Median EOS : ~{eos_date.strftime('%d %b')}")

        saved_ph = save_phenology(
            pm, OUT_DIR / "phenology",
            crs="EPSG:32633",
        )
        print("\n  GeoTIFFs written:")
        for name, p in saved_ph.items():
            print(f"    {p}")
    else:
        print("  NDVI not computed — no red/nir bands in stack")

    # 3.11  pearson_map
    _sep("3.11  pearson_map  (RED × NIR)")
    try:
        ri_r = band_names.index("red")
        ri_n = band_names.index("nir")
        r_map = pearson_map(f64[:, ri_r], f64[:, ri_n], dates)
        _stats_line(r_map, "Pearson r  RED-NIR")
        neg_frac = (r_map[r_map != NODATA] < 0).mean() * 100
        print(f"  Pixels with r < 0 : {neg_frac:.1f} %  (expected ~0 for healthy vegetation)")
        _save_tif(r_map, OUT_DIR / "pearson_red_nir.tif")
    except ValueError:
        print("  No red/nir bands — pearson_map skipped")

    _sep("Done")
    print(f"  All outputs written to: {OUT_DIR.resolve()}\n")


# Entry point

def main() -> None:
    parser = argparse.ArgumentParser(description="Sentinel-2 analysis via sentinel_stats")
    parser.add_argument("--skip-dl", action="store_true",
                        help="Skip download step (stack already cached)")
    parser.add_argument("--quick", action="store_true",
                        help="64x64 centre crop for a quick smoke-test")
    args = parser.parse_args()

    da, dates = step_download_and_stack(args.skip_dl)
    da_filled, filled64, band_names = step_gapfill(da)
    step_analysis(filled64, dates, band_names, args.quick)


if __name__ == "__main__":
    main()
