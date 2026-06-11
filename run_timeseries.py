from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np

STACK = Path("data/stacks/budapest_cloudy.nc")
SCL_DIR = Path("data/technical")
OUT_DIR = Path("data/vis")

DEFAULT_LON = 19.040
DEFAULT_LAT = 47.498


def _make_test_stack():
    import xarray as xr

    n_times = 24
    times = np.array(
        [f"202{3 + i // 12}-{(i % 12) + 1:02d}-01" for i in range(n_times)],
        dtype="datetime64[ns]",
    )

    lons = np.linspace(19.00, 19.10, 32)
    lats = np.linspace(47.44, 47.54, 32)

    rng = np.random.default_rng(0)

    month_idx = np.arange(n_times) % 12
    base = 0.3 + 0.45 * np.sin(np.pi * month_idx / 11) ** 2  # (n_times,)

    ndvi = (base[:, None, None]
            + rng.normal(0, 0.03, (n_times, 32, 32))).clip(0.05, 0.9).astype(np.float32)
    evi = (ndvi * 0.80
           + rng.normal(0, 0.02, (n_times, 32, 32))).clip(0.02, 0.7).astype(np.float32)

    for cloudy_t in (4, 10, 16):
        ndvi[cloudy_t] = np.nan
        evi[cloudy_t] = np.nan

    data = np.stack([ndvi, evi], axis=1)  # (time, band=2, y, x)

    da = xr.DataArray(
        data,
        dims=["time", "band", "y", "x"],
        coords={
            "time": times,
            "band": ["ndvi", "evi"],
            "y": lats,
            "x": lons,
        },
        attrs={"description": "Synthetic test stack — budapest_cloudy"},
    )
    return da


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot time series from budapest_cloudy.nc",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Modes:\n"
            "  live  Read from data/stacks/budapest_cloudy.nc  (default)\n"
            "  test  Build a 24-step synthetic stack in memory — no files needed\n"
        ),
    )
    p.add_argument("--mode", choices=["live", "test"], default="live",
                   help="'live' reads the real .nc file; 'test' uses synthetic data (default: live)")
    p.add_argument("--lon", type=float, default=DEFAULT_LON,
                   help="WGS-84 longitude  (default: %(default)s)")
    p.add_argument("--lat", type=float, default=DEFAULT_LAT,
                   help="WGS-84 latitude   (default: %(default)s)")
    p.add_argument("--bands", nargs="+", default=None,
                   help="Band names to plot, e.g. --bands ndvi evi  (default: all)")
    p.add_argument("--agg", type=float, default=None,
                   help="Spatial averaging half-width in degrees, e.g. 0.002")
    p.add_argument("--scl", type=Path, default=SCL_DIR,
                   help="SCL file or directory for cloud markers  (default: %(default)s)")
    p.add_argument("--out", type=Path, default=None,
                   help="HTML output path (default: data/vis/budapest_cloudy_ts.html or …_test.html)")
    p.add_argument("--no-show", action="store_true",
                   help="Do not open the browser after rendering")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    from sentinel_processor.visualisation.plot import plot_timeseries

    if args.mode == "test":
        print("  [test mode] Building synthetic 24-step stack in memory …")
        stack = _make_test_stack()
        scl_path = None
        out = args.out or (OUT_DIR / "budapest_cloudy_ts_test.html")
    else:
        if not STACK.exists():
            raise FileNotFoundError(
                f"Stack not found: {STACK}\n"
                "  Run with --mode test to use synthetic data instead."
            )
        stack = STACK
        scl_path = args.scl if args.scl.exists() else None
        if scl_path is None:
            print(f"  [info] SCL path not found ({args.scl}) — cloud markers disabled")
        out = args.out or (OUT_DIR / "budapest_cloudy_ts.html")

    out.parent.mkdir(parents=True, exist_ok=True)

    source_label = "synthetic (test)" if args.mode == "test" else str(STACK)
    print(f"  Stack   : {source_label}")
    print(f"  Point   : lon={args.lon}  lat={args.lat}")
    print(f"  Bands   : {args.bands or 'all'}")
    print(f"  Mode    : {'region ±' + str(args.agg) + '°' if args.agg else 'single pixel'}")
    print(f"  SCL     : {scl_path or 'none'}")
    print(f"  Output  : {out}")
    print()

    fig = plot_timeseries(
        stack=stack,
        lon=args.lon,
        lat=args.lat,
        bands=args.bands,
        scl_path=scl_path,
        agg_bbox=args.agg,
        save_html=str(out),
    )

    print(f"  Saved → {out}")

    if not args.no_show:
        fig.show()


if __name__ == "__main__":
    main()
