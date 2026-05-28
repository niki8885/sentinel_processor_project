
import sentinel_processor as sp

results = sp.download_sentinel2(
    [sp.LocationSpec(lat=47.56, lon=19.17, name='budapest_test')],
    cfg=sp.DownloadConfig(
        bands=sp.SpectralBands.RGB_NIR,
        tech_bands=None,
        visual=False,
        keep_items=1,
        lookback_days=14,
    ),
)
print(results)
