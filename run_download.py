from __future__ import annotations
import datetime
import os
import sys
from pathlib import Path

ROOT = Path(__file__).parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    import rasterio as _r
    _proj = os.path.join(os.path.dirname(_r.__file__), "proj_data")
    if os.path.isdir(_proj):
        os.environ.setdefault("PROJ_DATA", _proj)
        os.environ.setdefault("PROJ_LIB", _proj)
except Exception:
    pass

import logging
import sentinel_processor as sp

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(name)s │ %(message)s")
logging.getLogger("sentinel_processor").setLevel(logging.WARNING)

LOCATION = sp.LocationSpec(lat=47.562938, lon=19.169805, name="budapest_target")

cfg = sp.DownloadConfig(
    bands               = sp.SpectralBands.ALL,
    tech_bands          = sp.TechnicalLayers.ALL,
    visual              = True,
    output_dir          = str(ROOT / "downloads"),
    bbox_half_deg       = 0.045,
    start_date          = datetime.datetime.now(datetime.UTC) - datetime.timedelta(days=90),
    max_items           = 20,
    keep_items          = 5,
    validate            = True,
    max_cloud_threshold = 0.30,
    min_confidence      = 0.50,
    save_report         = True,
)


def main() -> None:
    print(f"\n  Target   {LOCATION.lat}, {LOCATION.lon}  →  '{LOCATION.name}'")
    print(f"  BBox     ±{cfg.bbox_half_deg}°  (~{cfg.bbox_half_deg * 111:.1f} km)")
    print(f"  Window   {cfg.resolved_start().date()} → {cfg.resolved_end().date()}")
    print(f"  Bands    {cfg.band_keys()}")
    print(f"  Validate cloud ≤ {cfg.max_cloud_threshold:.0%}  confidence ≥ {cfg.min_confidence:.0%}")
    print(f"  Output   {cfg.output_dir}\n")

    results = sp.download_sentinel2([LOCATION], cfg=cfg, progress=True)

    if not results:
        print("  No items passed validation.")
        return

    print("  Files written:")
    for base_name, paths in results.items():
        print(f"\n  [{base_name}]")
        for p in paths:
            rel = Path(p).relative_to(ROOT) if Path(p).is_relative_to(ROOT) else Path(p)
            print(f"    {rel}")


if __name__ == "__main__":
    main()