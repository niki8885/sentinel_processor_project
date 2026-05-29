import datetime
import sentinel_processor as sp

cfg = sp.DownloadConfig(
    bands               = sp.SpectralBands.ALL,
    tech_bands          = sp.TechnicalLayers.SCL,
    visual              = True,
    keep_items          = 5,
    max_cloud_threshold = 0.15,
    min_confidence      = 0.75,
)

results = sp.download_sentinel2(
    [
        sp.LocationSpec(lat=47.56, lon=19.17, name="budapest"),
    ],
    cfg=cfg,)