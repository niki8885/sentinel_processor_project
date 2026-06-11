from __future__ import annotations
from typing import Sequence

SSL4EO_S12_ALL: dict[str, dict[str, float]] = {
    "B02": {"mean": 1605.57, "std": 1390.58},
    "B03": {"mean": 1390.42, "std": 1244.92},
    "B04": {"mean": 1314.56, "std": 1347.54},
    "B05": {"mean": 1484.54, "std": 1228.58},
    "B06": {"mean": 2183.85, "std": 1316.56},
    "B07": {"mean": 2398.07, "std": 1355.18},
    "B08": {"mean": 2298.36, "std": 1383.44},
    "B8A": {"mean": 2525.60, "std": 1368.69},
    "B11": {"mean": 1796.35, "std": 1198.10},
    "B12": {"mean": 1241.28, "std": 990.22},
}

SECO_RGB: dict[str, dict[str, float]] = {
    "B04": {"mean": 1354.99, "std": 2173.26},
    "B03": {"mean": 1117.36, "std": 2065.40},
    "B02": {"mean": 990.22, "std": 2057.84},
}

SSL4EO_RGBN: dict[str, dict[str, float]] = {
    k: SSL4EO_S12_ALL[k] for k in ("B04", "B03", "B02", "B08")
}

# Registry: method name → band-keyed stats dict

PRESET_REGISTRY: dict[str, dict[str, dict[str, float]]] = {
    "sentinel2_all": SSL4EO_S12_ALL,
    "sentinel2_rgb": SECO_RGB,
    "sentinel2_rgbn": SSL4EO_RGBN,
}

AVAILABLE_PRESETS: list[str] = sorted(PRESET_REGISTRY)


def get_preset(name: str) -> dict[str, dict[str, float]]:
    """Return the stats dict for a named preset.

    Raises
    ------
    ValueError  if the preset name is unknown.
    """
    if name not in PRESET_REGISTRY:
        raise ValueError(
            f"Unknown preset '{name}'. "
            f"Available presets: {AVAILABLE_PRESETS}"
        )
    return PRESET_REGISTRY[name]


def preset_bands(name: str) -> list[str]:
    """Return the ordered list of bands covered by a preset."""
    return list(get_preset(name).keys())
