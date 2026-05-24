from __future__ import annotations

__version__ = "0.1.0"

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
]