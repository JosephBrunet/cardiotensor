"""Utilities module."""

from typing import Generator

import numpy as np
import numpy.typing as npt

try:
    import cupy as cp
except Exception as ex:
    cp = None


def _calculate_kernel_radius(window, truncate):
    return int(window * truncate + 1)


def get_block_count(data: npt.NDArray, block_size: int = 512) -> int:
    """Gets the number of blocks that will be created for the given input."""
    return np.prod(np.ceil(np.array(data.shape) / block_size).astype(int)).item()


def get_block(
    i: int,
    data: npt.NDArray,
    window: float,
    block_size: int = 512,
    truncate: float = 4.0,
    copy: bool = False,
) -> tuple[npt.NDArray, np.ndarray, np.ndarray]:
    """Gets the ith block."""

    kernel_radius = _calculate_kernel_radius(window, truncate)
    count = 0

    for x0 in range(0, data.shape[0], block_size):
        for y0 in range(0, data.shape[1], block_size):
            for z0 in range(0, data.shape[2], block_size):
                if count == i:
                    x1 = x0 + block_size
                    y1 = y0 + block_size
                    z1 = z0 + block_size

                    block = data[
                        max(0, x0 - kernel_radius) : x1 + kernel_radius,
                        max(0, y0 - kernel_radius) : y1 + kernel_radius,
                        max(0, z0 - kernel_radius) : z1 + kernel_radius,
                    ]

                    cx0 = kernel_radius + min(0, x0 - kernel_radius)
                    cy0 = kernel_radius + min(0, y0 - kernel_radius)
                    cz0 = kernel_radius + min(0, z0 - kernel_radius)
                    cx1 = max(0, min(kernel_radius, data.shape[0] - x1))
                    cy1 = max(0, min(kernel_radius, data.shape[1] - y1))
                    cz1 = max(0, min(kernel_radius, data.shape[2] - z1))

                    if copy:
                        block = np.array(block)

                    return (
                        block,
                        np.array(((x0, x1), (y0, y1), (z0, z1))),
                        np.array(((cx0, cx1), (cy0, cy1), (cz0, cz1))),
                    )

                count += 1

    raise IndexError(f"Index {i} is out of bounds for {count} blocks.")

def get_block_generator(
    data: npt.NDArray,
    sigma: float,
    block_size: int = 512,
    truncate: float = 4.0,
    copy: bool = False,
) -> Generator[tuple[npt.NDArray, np.ndarray, np.ndarray], None, None]:
    """
    Generator yielding (block, position, padding) for vector field data in (3, z, y, x).

    Parameters
    ----------
    data : np.ndarray
        Input data of shape (3, z, y, x).
    sigma : float
        Smoothing kernel sigma.
    block_size : int
        Spatial size of each block (applied to z, y, x).
    truncate : float
        Kernel truncation factor.
    copy : bool
        Whether to copy the block data.

    Yields
    ------
    block : np.ndarray
        The extracted block including kernel padding, shape (3, bz, by, bx)
    pos : np.ndarray
        Block spatial indices in the original volume: ((z0, z1), (y0, y1), (x0, x1))
    pad : np.ndarray
        Padding applied on each side: ((cz0, cz1), (cy0, cy1), (cx0, cx1))
    """
    kernel_radius = _calculate_kernel_radius(sigma, truncate)
    _, Z, Y, X = data.shape

    for z0 in range(0, Z, block_size):
        z1 = min(z0 + block_size, Z)
        for y0 in range(0, Y, block_size):
            y1 = min(y0 + block_size, Y)
            for x0 in range(0, X, block_size):
                x1 = min(x0 + block_size, X)

                # Slice including kernel padding
                block = data[
                    :,
                    max(0, z0 - kernel_radius) : z1 + kernel_radius,
                    max(0, y0 - kernel_radius) : y1 + kernel_radius,
                    max(0, x0 - kernel_radius) : x1 + kernel_radius,
                ]

                # Compute actual padding applied (if near edges)
                cz0 = kernel_radius - max(0, kernel_radius - z0)
                cy0 = kernel_radius - max(0, kernel_radius - y0)
                cx0 = kernel_radius - max(0, kernel_radius - x0)
                cz1 = kernel_radius - max(0, kernel_radius - (Z - z1))
                cy1 = kernel_radius - max(0, kernel_radius - (Y - y1))
                cx1 = kernel_radius - max(0, kernel_radius - (X - x1))

                if copy:
                    block = np.array(block)

                yield (
                    block,
                    np.array(((z0, z1), (y0, y1), (x0, x1))),
                    np.array(((cz0, cz1), (cy0, cy1), (cx0, cx1))),
                )

def get_blocks(
    data: npt.NDArray,
    window: float,
    block_size: int = 512,
    truncate: float = 4.0,
    copy: bool = False,
) -> tuple[list[npt.NDArray], np.ndarray, np.ndarray]:
    """Gets a tuple of blocks, block positions and block paddings."""

    blocks = []
    block_positions = []
    block_paddings = []

    for block, pos, pad in get_block_generator(data, window, block_size=block_size, truncate=truncate, copy=copy):
        blocks.append(block)
        block_positions.append(pos)
        block_paddings.append(pad)

    return blocks, np.array(block_positions), np.array(block_paddings)


def remove_padding(block: npt.NDArray, pad: npt.NDArray[np.integer]) -> npt.NDArray:
    """Slices away the block padding."""

    block = block[
        ...,
        pad[0, 0] : block.shape[-3] - pad[0, 1],
        pad[1, 0] : block.shape[-2] - pad[1, 1],
        pad[2, 0] : block.shape[-1] - pad[2, 1],
    ]
    return block


def insert_block(
    volume: npt.NDArray,
    block: npt.NDArray,
    pos: npt.NDArray[np.integer],
    pad: npt.NDArray[np.integer] | None = None,
    mask: npt.NDArray[np.bool_] | None = None,
):
    """
    Inserts a block into a volume at a specific position, removing padding if needed.

    Parameters
    ----------
    volume : np.ndarray
        The output volume, shape (z, y, x).
    block : np.ndarray
        Block to insert. Shape can be padded: (z_pad, y_pad, x_pad).
    pos : np.ndarray
        Start/end positions of block in volume: shape ((x0,x1),(y0,y1),(z0,z1)).
    pad : np.ndarray | None
        Optional padding applied to block: shape ((cx0,cx1),(cy0,cy1),(cz0,cz1)).
    mask : np.ndarray | None
        Optional mask for selective insertion.

    Notes
    -----
    Handles channel-first blocks (e.g., (3, z, y, x)) automatically by inserting
    only the spatial dimensions.
    """

    # Remove padding if necessary
    if pad is not None:
        block = remove_padding(block, pad)

    # Determine slicing for spatial dimensions
    z0, z1 = pos[2]
    y0, y1 = pos[1]
    x0, x1 = pos[0]

    if block.ndim == volume.ndim + 1:
        # Input block has extra channel dimension (e.g., (3, z, y, x)), take spatial
        block_spatial = block[..., :, :, :]
        volume_view = volume[x0:x1, y0:y1, z0:z1]  # volume is (x,y,z) order
    else:
        block_spatial = block
        volume_view = volume[x0:x1, y0:y1, z0:z1]

    # Handle potential GPU arrays
    if cp is not None and isinstance(block_spatial, cp.ndarray):
        block_spatial = cp.asnumpy(block_spatial.astype(volume.dtype))

    # Insert block
    if mask is None:
        volume_view[:] = block_spatial
    else:
        volume_view[..., mask] = block_spatial