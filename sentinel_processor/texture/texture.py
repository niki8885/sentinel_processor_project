from __future__ import annotations
from typing import Literal
import numpy as np
from ._texture_bridge import NODATA, Angle

__all__ = ["compute_glcm", "NODATA"]

_N_LEVELS = 64


def _quantise(arr: np.ndarray) -> np.ndarray:
    valid = arr > NODATA + 1.0
    out = np.zeros(arr.shape, dtype=np.int32)
    if not valid.any():
        return out
    bmin = float(arr[valid].min())
    bmax = float(arr[valid].max())
    rng = bmax - bmin
    if rng < 1e-12:
        rng = 1e-12
    lev = ((arr - bmin) / rng * (_N_LEVELS - 1)).astype(np.int32) + 1
    lev = np.clip(lev, 1, _N_LEVELS)
    out[valid] = lev[valid]
    return out


def _glcm_numpy(
        qlev: np.ndarray,
        r0: int, c0: int,
        half_w: int,
        dr: int, dc: int,
) -> tuple[float, float, float]:
    rows, cols = qlev.shape
    r_lo, r_hi = max(r0 - half_w, 0), min(r0 + half_w + 1, rows)
    c_lo, c_hi = max(c0 - half_w, 0), min(c0 + half_w + 1, cols)

    patch_i = qlev[r_lo:r_hi, c_lo:c_hi]
    r1_lo = r_lo + dr;
    r1_hi = r_hi + dr
    c1_lo = c_lo + dc;
    c1_hi = c_hi + dc

    # Clip neighbour patch to valid raster extent
    pr_lo = max(r1_lo, 0);
    pr_hi = min(r1_hi, rows)
    pc_lo = max(c1_lo, 0);
    pc_hi = min(c1_hi, cols)

    # Matching offsets in patch_i
    pi_r_lo = pr_lo - r1_lo;
    pi_r_hi = pr_hi - r1_lo
    pi_c_lo = pc_lo - c1_lo;
    pi_c_hi = pc_hi - c1_lo

    patch_j = qlev[pr_lo:pr_hi, pc_lo:pc_hi]
    patch_i = patch_i[pi_r_lo:pi_r_hi, pi_c_lo:pi_c_hi]

    valid = (patch_i > 0) & (patch_j > 0)
    if not valid.any():
        return 0.0, 0.0, 0.0

    pi_v = patch_i[valid].ravel()
    pj_v = patch_j[valid].ravel()

    glcm = np.zeros((_N_LEVELS + 1, _N_LEVELS + 1), dtype=np.float64)
    np.add.at(glcm, (pi_v, pj_v), 1.0)
    np.add.at(glcm, (pj_v, pi_v), 1.0)  # symmetrise
    glcm = glcm[1:, 1:]  # drop row/col 0

    total = glcm.sum()
    if total < 1e-12:
        return 0.0, 0.0, 0.0
    glcm /= total

    ii, jj = np.meshgrid(np.arange(1, _N_LEVELS + 1),
                         np.arange(1, _N_LEVELS + 1), indexing="ij")
    diff2 = (ii - jj) ** 2

    energy = float((glcm ** 2).sum())
    contrast = float((diff2 * glcm).sum())
    homogeneity = float((glcm / (1.0 + diff2)).sum())
    return energy, contrast, homogeneity


_ANGLE_OFFSETS: dict[int, tuple[int, int]] = {
    0: (0, 1),
    45: (-1, 1),
    90: (-1, 0),
    135: (-1, -1),
}


def _compute_glcm_numpy(
        arr: np.ndarray,
        window: int,
        distance: int,
        angle: int,
) -> dict[str, np.ndarray]:
    """NumPy fallback"""
    rows, cols = arr.shape
    half_w = window // 2
    border = half_w + distance

    qlev = _quantise(arr)

    energy_out = np.full((rows, cols), NODATA, dtype=np.float64)
    contrast_out = np.full((rows, cols), NODATA, dtype=np.float64)
    homogeneity_out = np.full((rows, cols), NODATA, dtype=np.float64)

    angles = list(_ANGLE_OFFSETS.keys()) if angle == -1 else [angle]

    for r in range(border, rows - border):
        for c in range(border, cols - border):
            en_sum = co_sum = ho_sum = 0.0
            for a in angles:
                dr, dc = _ANGLE_OFFSETS[a]
                en, co, ho = _glcm_numpy(qlev, r, c, half_w,
                                         dr * distance, dc * distance)
                en_sum += en
                co_sum += co
                ho_sum += ho
            n = len(angles)
            energy_out[r, c] = en_sum / n
            contrast_out[r, c] = co_sum / n
            homogeneity_out[r, c] = ho_sum / n

    return {
        "energy": energy_out,
        "contrast": contrast_out,
        "homogeneity": homogeneity_out,
    }


# Public entry point

def compute_glcm(
        arr: np.ndarray,
        window: int = 7,
        distance: int = 1,
        angle: int = -1,
) -> dict[str, np.ndarray]:
    arr = np.asarray(arr, dtype=np.float64)

    if arr.ndim != 2:
        raise ValueError(f"arr must be 2-D (rows, cols), got shape {arr.shape}")

    valid_angles = {-1, 0, 45, 90, 135}
    if angle not in valid_angles:
        raise ValueError(
            f"angle must be one of {valid_angles}, got {angle!r}"
        )

    if window < 3:
        window = 3
    if window % 2 == 0:
        window += 1

    try:
        from ._texture_bridge import compute_glcm as _fortran_compute_glcm
        return _fortran_compute_glcm(arr, window=window,
                                     distance=distance, angle=angle)
    except FileNotFoundError as exc:
        import warnings
        warnings.warn(
            f"{exc}\n"
            "Falling back to pure-NumPy GLCM implementation.  "
            "This is significantly slower on large tiles.\n"
            "See the module docstring for build instructions.",
            RuntimeWarning,
            stacklevel=2,
        )
        return _compute_glcm_numpy(arr, window=window,
                                   distance=distance, angle=angle)
