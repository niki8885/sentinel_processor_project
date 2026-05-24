from __future__ import annotations
from dataclasses import dataclass
from typing import List

class _BandGroup:

    def __init__(self, keys: List[str]) -> None:
        self.keys = keys
    def __iter__(self):
        return iter(self.keys)
    def __add__(self, other: _BandGroup) -> _BandGroup:
        return _BandGroup(self.keys + other.keys)
    def __repr__(self) -> str:
        return f"_BandGroup(keys={self.keys!r})"


class SpectralBands:
    # 10 m
    BLUE = _BandGroup(["blue"]) # B02 490 nm
    GREEN = _BandGroup(["green"]) # B03 560 nm
    RED = _BandGroup(["red"]) # B04 665 nm
    NIR = _BandGroup(["nir"]) # B08 842 nm

    # 20 m
    RED_EDGE_1 = _BandGroup(["rededge1"]) # B05 705 nm
    RED_EDGE_2 = _BandGroup(["rededge2"]) # B06 740 nm
    RED_EDGE_3 = _BandGroup(["rededge3"]) # B07 783 nm
    NIR_NARROW = _BandGroup(["nir08"]) # B8A 865 nm
    SWIR_1 = _BandGroup(["swir16"]) # B11 1610 nm
    SWIR_2 = _BandGroup(["swir22"]) # B12 2190 nm

    # 60 m
    COASTAL_AEROSOL = _BandGroup(["coastal"]) # B01 443 nm
    WATER_VAPOUR = _BandGroup(["nir09"]) # B09 945 nm

    # Presets
    RGB = BLUE + GREEN + RED
    RGB_NIR = RGB + NIR
    VEGETATION = RED + NIR + RED_EDGE_1 + RED_EDGE_2 + RED_EDGE_3
    AGRICULTURE = BLUE + GREEN + RED + NIR + RED_EDGE_1 + SWIR_1 + SWIR_2
    ALL_10M = BLUE + GREEN + RED + NIR
    ALL_20M = RED_EDGE_1 + RED_EDGE_2 + RED_EDGE_3 + NIR_NARROW + SWIR_1 + SWIR_2
    ALL = ALL_10M + ALL_20M


class TechnicalLayers:
    SCL = _BandGroup(["scl"]) # Scene Classification Layer
    AOT = _BandGroup(["aot"]) # Aerosol Optical Thickness
    WVP = _BandGroup(["wvp"]) # Water Vapour Map
    ALL = SCL + AOT + WVP


class VisualAssets:
    VISUAL = _BandGroup(["visual"]) # True-colour 10 m
    THUMBNAIL = _BandGroup(["thumbnail"]) # Low-res preview

@dataclass
class LocationSpec:
    lon: float
    lat: float
    name: str | None = None

    def base_name(self, timestamp) -> str:
        if self.name:
            return self.name
        return f"{timestamp.strftime('%Y-%m-%d')}_{self.lat:.4f}_{self.lon:.4f}"