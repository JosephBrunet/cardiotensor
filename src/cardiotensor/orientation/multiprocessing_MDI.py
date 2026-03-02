import logging
import os
from dataclasses import dataclass
from multiprocessing import Pool, RawArray, SimpleQueue, cpu_count
from multiprocessing.pool import ThreadPool
from typing import Callable, Literal, Sequence
import numpy as np
import numpy.typing as npt
from scipy.ndimage import uniform_filter

from . import util

logger = logging.getLogger(__name__)

DEFAULT_POOL_TYPE = "thread" if os.name == "nt" else "process"


# -----------------------------------------------------------------------------
# Shared memory helpers (unchanged pattern)
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class _ArrayArgs:
    array: np.ndarray
    def get_array(self) -> np.ndarray:
        return self.array


@dataclass(frozen=True)
class _RawArrayArgs:
    array: any
    shape: tuple[int, ...]
    dtype: npt.DTypeLike
    def get_array(self) -> np.ndarray:
        return np.frombuffer(self.array, dtype=self.dtype).reshape(self.shape)


def _create_raw_array(shape: tuple[int, ...], dtype: npt.DTypeLike):
    raw = RawArray(
        "b",
        np.prod(np.asarray(shape), dtype=np.int64).item()
        * np.dtype(dtype).itemsize,
    )
    a = np.frombuffer(raw, dtype=dtype).reshape(shape)
    return raw, a


# -----------------------------------------------------------------------------
# Init args
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class _InitArgs:
    data_args: _RawArrayArgs | _ArrayArgs
    output_args: _RawArrayArgs | _ArrayArgs
    window_size: int
    block_size: int
    devices: SimpleQueue

    def get_data_sources(self):
        return (
            self.data_args.get_array(),
            self.output_args.get_array(),
            self.devices.get(),
        )


_worker_args = None
_data = None
_output = None
_device = None


# -----------------------------------------------------------------------------
# Worker init
# -----------------------------------------------------------------------------

def _init_worker(init_args: _InitArgs):
    global _worker_args, _data, _output, _device
    _worker_args = init_args
    _data, _output, _device = init_args.get_data_sources()


# -----------------------------------------------------------------------------
# Core MDI computation (block-level)
# -----------------------------------------------------------------------------

def _compute_mdi_block(block: np.ndarray, window_size: int):

    vx, vy, vz = block

    norm = np.sqrt(vx**2 + vy**2 + vz**2)
    mask = norm > 0

    vx = np.where(mask, vx / norm, 0)
    vy = np.where(mask, vy / norm, 0)
    vz = np.where(mask, vz / norm, 0)

    Sxx = uniform_filter(vx * vx, size=window_size)
    Syy = uniform_filter(vy * vy, size=window_size)
    Szz = uniform_filter(vz * vz, size=window_size)

    Sxy = uniform_filter(vx * vy, size=window_size)
    Sxz = uniform_filter(vx * vz, size=window_size)
    Syz = uniform_filter(vy * vz, size=window_size)

    S = np.stack([
        np.stack([Sxx, Sxy, Sxz], axis=-1),
        np.stack([Sxy, Syy, Syz], axis=-1),
        np.stack([Sxz, Syz, Szz], axis=-1),
    ], axis=-2)

    eigvals = np.linalg.eigvalsh(S)[..., ::-1]

    l1 = eigvals[..., 0]
    l2 = eigvals[..., 1]
    l3 = eigvals[..., 2]

    l0 = (l1 + l2 + l3) / 3.0

    mdi = np.zeros_like(l0)
    valid = l0 > 0

    mdi[valid] = np.sqrt(
        (l1[valid] - l0[valid]) ** 2 +
        (l2[valid] - l0[valid]) ** 2 +
        (l3[valid] - l0[valid]) ** 2
    ) / (np.sqrt(6.0) * l0[valid])

    return mdi.astype(np.float32)


# -----------------------------------------------------------------------------
# Worker
# -----------------------------------------------------------------------------

def _do_work(block_info: tuple[npt.NDArray, np.ndarray, np.ndarray]):
    """
    Worker function for parallel MDI calculation.
    
    Parameters
    ----------
    block_info : tuple
        (block, pos, pad) from get_block_generator
    
    Returns
    -------
    tuple
        (pos, mdi_block) to insert into the output array
    """
    if _worker_args is None:
        raise RuntimeError("Worker not initialized")

    block, pos, pad = block_info

    # Ensure block is float64 for calculations
    block = np.asarray(block, dtype=np.float64)

    # Compute MDI for this block
    mdi_block = _compute_mdi_block(block, _worker_args.window_size)

    # Insert the MDI block into the shared output array
    util.insert_block(_output, mdi_block, pos, pad)

    return pos  # returning pos can help track progress




def parallel_mdi_analysis(
    v3: np.ndarray,
    window_size: int = 9,
    block_size: int = 128,
    output: np.memmap | npt.DTypeLike | None = np.float32,
    devices: Sequence[str] | None = None,
    progress_callback_fn: Callable[[int, int], None] | None = None,
    pool_type: Literal["process", "thread"] = DEFAULT_POOL_TYPE,
):
    """
    Parallel MDI calculation using block-wise processing.

    Parameters
    ----------
    v3 : np.ndarray
        Vector field with shape (3, z, y, x).
    window_size : int
        Window size for MDI computation.
    block_size : int
        Block size for splitting the volume.
    output : np.memmap | dtype | None
        Output array or dtype. If None, defaults to float32.
    devices : Sequence[str] | None
        List of devices to use, e.g., ["cpu", "cpu", "cuda:0"].
    progress_callback_fn : callable | None
        Called with (count, total_blocks) after each block is done.
    pool_type : {"process", "thread"}
        Pool type for parallel execution.

    Returns
    -------
    np.ndarray
        The computed MDI volume with shape (z, y, x).
    """

    if v3.shape[0] != 3:
        raise ValueError("Input must have shape (3, z, y, x)")

    devices = devices or ["cpu"] * cpu_count()
    use_process_pool = pool_type == "process"

    # Output array
    output_shape = v3.shape[1:]  # (z, y, x)
    if isinstance(output, np.memmap):
        output_array = output
    else:
        _, output_array = _create_raw_array(output_shape, output)

    # Shared globals for workers
    global _data, _output, _worker_args
    _data = v3
    _output = output_array
    _worker_args = type("WorkerArgs", (), {})()
    _worker_args.window_size = window_size
    _worker_args.block_size = block_size

    # Generator over blocks
    block_gen = util.get_block_generator(
        v3,
        sigma=window_size,
        block_size=block_size,
    )

    pool_ctor = Pool if use_process_pool else ThreadPool

    total_blocks = util.get_block_count(v3[0], block_size)
    count = 0

    with pool_ctor(processes=len(devices)) as pool:
        for _ in pool.imap_unordered(_do_work, block_gen, chunksize=1):
            count += 1
            if progress_callback_fn:
                progress_callback_fn(count, total_blocks)

    return output_array