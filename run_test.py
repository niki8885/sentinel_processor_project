import datetime
from pathlib import Path
import sentinel_processor as sp

# cfg = sp.DownloadConfig(
#     bands               = sp.SpectralBands.ALL,
#     tech_bands          = sp.TechnicalLayers.SCL,
#     visual              = True,
#     keep_items          = 20,
#     max_cloud_threshold = 0.15,
#     min_confidence      = 0.75,
#     pansharpen_algorithm = "gram_schmidt",
# )
#
# results = sp.download_sentinel2(
#     [
#         sp.LocationSpec(lat=47.56, lon=19.17, name="budapest"),
#     ],
#     cfg=cfg,)


from sentinel_processor import stack_timeseries, TimeSeriesConfig

result = stack_timeseries(
    sources  = sorted(Path("data/spectral").glob("budapest_*.nc")),
    scl_dir  = "data/technical",
    cfg      = TimeSeriesConfig(
        max_cloud_fraction = 0.10,
        min_confidence     = 0.75,
        require_bands      = ["red", "nir"],
        fill_rejected      = False,
        save_dir="data/stacks",
        save_format="nc",
        save_name="budapest_may",
    ),
)
print(result.summary())
result.to_json("stack_log.json")
da = result.stack