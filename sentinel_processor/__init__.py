from __future__ import annotations

__version__ = "0.1.1"

from sentinel_processor.input.downloader import DownloadConfig, download_sentinel2
from sentinel_processor.utils.data_utils import (
    LocationSpec,
    SpectralBands,
    TechnicalLayers,
    VisualAssets,
)
from sentinel_processor.validation._fortran_bridge import (
    call_check_radiometry,
    call_validate_scl,
    validate_file,
)

__all__ = [
    # download
    "download_sentinel2",
    "DownloadConfig",
    # location & bands
    "LocationSpec",
    "SpectralBands",
    "TechnicalLayers",
    "VisualAssets",
    # validation
    "call_validate_scl",
    "call_check_radiometry",
    "validate_file",
    "stack_timeseries",
    "TimeSeriesConfig",
    "StackResult",
]


def __getattr__(name: str):

    _timeseries_symbols = {"stack_timeseries", "TimeSeriesConfig", "StackResult"}
    if name in _timeseries_symbols:
        from sentinel_processor.input import timeseries as _ts
        return getattr(_ts, name)
    raise AttributeError(f"module 'sentinel_processor' has no attribute {name!r}")