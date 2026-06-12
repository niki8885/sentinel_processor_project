import datetime
from pathlib import Path
import numpy as np
import xarray as xr
import sentinel_processor as sp
from sentinel_processor import stack_timeseries, TimeSeriesConfig
from sentinel_processor.processing._timeseries_bridge import interpolate_gaps

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

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name="budapest")],
    cfg=cfg,
)

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
)
print(result.summary())
result.to_json("stack_log.json")
da = result.stack

arr = da.values.astype(np.float64)
mask = np.isfinite(arr).astype(np.int32)

filled = np.empty_like(arr)
for b in range(arr.shape[1]):
    filled[:, b] = interpolate_gaps(arr[:, b], mask[:, b], method="pchip")

da_filled = xr.DataArray(
    filled.astype(np.float32),
    dims=da.dims,
    coords=da.coords,
    attrs=da.attrs,
)
da_filled = da_filled.where(da_filled != -9999.0)

gaps_before = (~np.isfinite(da.values)).sum()
gaps_after = np.isnan(da_filled.values).sum()
band_names = list(da.coords["band"].values) if "band" in da.coords else [f"b{i}" for i in range(da.shape[1])]

print(f"\nGaps before : {gaps_before:,}")
print(f"Gaps after  : {gaps_after:,}  (pixels with no valid observations at all)")
print(f"Has NaN     : {np.isnan(da_filled.values).any()}")
print(f"Has Inf     : {np.isinf(da_filled.values).any()}")
print(f"Value range : {np.nanmin(da_filled.values):.4f} ... {np.nanmax(da_filled.values):.4f}")

print("\nPer-band fill rate:")
for i, name in enumerate(band_names):
    was_gap = ~np.isfinite(da.values[:, i])
    now_filled = was_gap & np.isfinite(da_filled.values[:, i])
    total = was_gap.sum()
    print(f"  {name:>8}:  {now_filled.sum():,} / {total:,}  ({100 * now_filled.sum() / max(total, 1):.1f} %)")

cy, cx = da.shape[2] // 2, da.shape[3] // 2
print(f"\nCentre pixel [{cy}, {cx}], band 0:")
for t in range(da.shape[0]):
    orig = da.values[t, 0, cy, cx]
    fill = da_filled.values[t, 0, cy, cx]
    flag = "<-- filled" if not np.isfinite(orig) else ""
    print(f"  t={t:02d}  before={orig:8.4f}   after={fill:8.4f}  {flag}")
