from __future__ import annotations
import math
from typing import NamedTuple
import numpy as np


class TileMeta(NamedTuple):
    """Geometry descriptor returned by :func:`extract_tiles`."""
    original_shape: tuple[int, int, int]  # (C, H, W) before padding
    padded_shape: tuple[int, int, int]  # (C, H_pad, W_pad) after padding
    tile_size: int
    overlap: int
    stride: int
    grid_rows: int  # number of tile rows
    grid_cols: int  # number of tile columns
    pad_h: int  # rows added at bottom
    pad_w: int  # cols added at right
    nodata: float


# Core functions

def extract_tiles(
        arr: np.ndarray,
        tile_size: int = 256,
        overlap: int = 32,
        nodata: float = -9999.0,
) -> tuple[np.ndarray, TileMeta]:
    """Cut a (C, H, W) array into overlapping square tiles.

    The image is padded (reflect mode) on the bottom and right edges so
    that the tile grid covers the full image.  Tiles are returned as a
    contiguous float64 array.

    Parameters
    ----------
    arr : np.ndarray
        Shape ``(C, H, W)``.  A 2-D ``(H, W)`` array is treated as
        single-band and returned as single-band tiles.
    tile_size : int
        Spatial size of each tile (pixels).  Must be > ``2 * overlap``.
    overlap : int
        Number of pixels of overlap between adjacent tiles.  Must be >= 0.
        The effective stride is ``tile_size - overlap``.
    nodata : float
        Sentinel value used when padding with ``reflect`` is unsuitable;
        carried through to :class:`TileMeta` for use in
        :func:`stitch_tiles`.

    Returns
    -------
    tiles : np.ndarray
        Shape ``(N_tiles, C, tile_size, tile_size)``, float64.
    meta : TileMeta
        Everything needed to reconstruct the original extent via
        :func:`stitch_tiles`.
    """
    if tile_size <= 0:
        raise ValueError(f"tile_size must be > 0, got {tile_size}")
    if overlap < 0:
        raise ValueError(f"overlap must be >= 0, got {overlap}")
    if tile_size <= 2 * overlap:
        raise ValueError(
            f"tile_size ({tile_size}) must be > 2 * overlap ({2 * overlap})"
        )

    # Normalise to (C, H, W)
    squeeze = False
    if arr.ndim == 2:
        arr = arr[np.newaxis]
        squeeze = True
    if arr.ndim != 3:
        raise ValueError(f"arr must be 2-D or 3-D, got shape {arr.shape}")

    arr = np.asarray(arr, dtype=np.float64)
    C, H, W = arr.shape
    stride = tile_size - overlap

    def _pad_amount(size: int) -> int:
        if size <= tile_size:
            return tile_size - size
        n_steps = math.ceil((size - tile_size) / stride)
        return n_steps * stride + tile_size - size

    pad_h = _pad_amount(H)
    pad_w = _pad_amount(W)
    H_pad = H + pad_h
    W_pad = W + pad_w

    # Reflect-pad
    arr_pad = np.pad(arr, ((0, 0), (0, pad_h), (0, pad_w)), mode="reflect")

    grid_rows = (H_pad - tile_size) // stride + 1
    grid_cols = (W_pad - tile_size) // stride + 1
    N = grid_rows * grid_cols

    tiles = np.empty((N, C, tile_size, tile_size), dtype=np.float64)
    idx = 0
    for r in range(grid_rows):
        for c in range(grid_cols):
            y0 = r * stride
            x0 = c * stride
            tiles[idx] = arr_pad[:, y0: y0 + tile_size, x0: x0 + tile_size]
            idx += 1

    meta = TileMeta(
        original_shape=(C, H, W),
        padded_shape=(C, H_pad, W_pad),
        tile_size=tile_size,
        overlap=overlap,
        stride=stride,
        grid_rows=grid_rows,
        grid_cols=grid_cols,
        pad_h=pad_h,
        pad_w=pad_w,
        nodata=nodata,
    )

    if squeeze:
        tiles = tiles[:, 0]

    return tiles, meta


def stitch_tiles(
        tiles: np.ndarray,
        meta: TileMeta,
        blend: bool = True,
) -> np.ndarray:
    """Reconstruct an image from overlapping tiles.

    Parameters
    ----------
    tiles : np.ndarray
        Shape ``(N_tiles, C, tile_size, tile_size)`` or
        ``(N_tiles, tile_size, tile_size)`` for single-band output.
        This is typically the output of your model applied to tiles
        produced by :func:`extract_tiles`.
    meta : TileMeta
        The :class:`TileMeta` object returned alongside the tiles.
    blend : bool
        If ``True`` (default), overlapping regions are averaged using
        a smooth cosine blending weight.  This removes seam artefacts.
        If ``False``, later tiles overwrite earlier ones in overlap zones.

    Returns
    -------
    np.ndarray
        Reconstructed array cropped to the original ``(C, H, W)`` shape
        (padding introduced by :func:`extract_tiles` is removed).

    """
    C_orig, H_orig, W_orig = meta.original_shape
    _, H_pad, W_pad = meta.padded_shape
    tile_size = meta.tile_size
    stride = meta.stride
    overlap = meta.overlap

    # Handle single-band tiles (N, T, T) → (N, 1, T, T)
    squeeze = False
    if tiles.ndim == 3:
        tiles = tiles[:, np.newaxis]
        squeeze = True

    N, C_tile, T, _ = tiles.shape
    assert T == tile_size, f"Tile spatial size mismatch: {T} vs {tile_size}"

    canvas = np.zeros((C_tile, H_pad, W_pad), dtype=np.float64)
    weights = np.zeros((H_pad, W_pad), dtype=np.float64)

    _WEIGHT_EPS = 1e-6
    if blend and overlap > 0:
        weight_1d = _cosine_blend_1d(tile_size, overlap)
        weight_2d = np.outer(weight_1d, weight_1d)
        weight_2d = np.clip(weight_2d, _WEIGHT_EPS, None)
    else:
        weight_2d = np.ones((tile_size, tile_size), dtype=np.float64)

    idx = 0
    for r in range(meta.grid_rows):
        for c in range(meta.grid_cols):
            y0 = r * stride
            x0 = c * stride
            canvas[:, y0: y0 + tile_size, x0: x0 + tile_size] += (
                    tiles[idx] * weight_2d
            )
            weights[y0: y0 + tile_size, x0: x0 + tile_size] += weight_2d
            idx += 1

    mask = weights > 1e-12
    canvas[:, mask] /= weights[mask]

    result = canvas[:, :H_orig, :W_orig]

    if squeeze:
        result = result[0]  # (H, W)

    return result


# Helpers

def _cosine_blend_1d(tile_size: int, overlap: int) -> np.ndarray:
    """1-D cosine blend: tapers from 0→1 over *overlap* pixels at each edge."""
    w = np.ones(tile_size, dtype=np.float64)
    for i in range(overlap):
        fade = 0.5 * (1.0 - math.cos(math.pi * i / overlap))
        w[i] = fade
        w[tile_size - 1 - i] = fade
    return w


# Geometry utility

def tile_grid_shape(
        image_h: int,
        image_w: int,
        tile_size: int = 256,
        overlap: int = 32,
) -> tuple[int, int, int]:
    """Return (n_tiles, grid_rows, grid_cols) without allocating memory."""
    if tile_size <= 0 or overlap < 0 or tile_size <= 2 * overlap:
        raise ValueError("Invalid tile_size / overlap combination.")
    stride = tile_size - overlap

    def _grid(size: int) -> int:
        if size <= tile_size:
            return 1
        pad = math.ceil((size - tile_size) / stride) * stride + tile_size - size
        return (size + pad - tile_size) // stride + 1

    rows = _grid(image_h)
    cols = _grid(image_w)
    return rows * cols, rows, cols
