from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from cardiotensor.utils.downsampling import process_vector_block

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_npy_slices(directory: Path, n: int, shape=(3, 8, 8)) -> list[Path]:
    """Write n random (3, H, W) .npy slices and return sorted paths."""
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for i in range(n):
        arr = np.random.rand(*shape).astype(np.float32)
        # normalise so vectors are unit length
        norms = np.linalg.norm(arr, axis=0, keepdims=True)
        arr = arr / np.maximum(norms, 1e-12)
        p = directory / f"slice_{i:04d}.npy"
        np.save(p, arr)
        paths.append(p)
    return sorted(paths)


# ---------------------------------------------------------------------------
# process_vector_block
# ---------------------------------------------------------------------------


def test_process_vector_block_creates_file(tmp_path: Path):
    H, W = 8, 8
    bin_factor = 2
    paths = _write_npy_slices(tmp_path / "input", n=bin_factor)
    out_dir = tmp_path / "output"
    process_vector_block(
        paths, bin_factor=bin_factor, h=H, w=W, output_dir=out_dir, idx=0
    )
    out_file = out_dir / "eigen_vec" / "eigen_vec_000000.npy"
    assert out_file.exists()
    result = np.load(out_file)
    assert result.shape == (3, math.ceil(H / bin_factor), math.ceil(W / bin_factor))


def test_process_vector_block_skips_existing(tmp_path: Path):
    H, W = 8, 8
    bin_factor = 2
    paths = _write_npy_slices(tmp_path / "input", n=bin_factor)
    out_dir = tmp_path / "output"
    out_file = out_dir / "eigen_vec" / "eigen_vec_000000.npy"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    # Write a sentinel file
    sentinel = np.zeros((3, H // bin_factor, W // bin_factor), dtype=np.float32)
    np.save(out_file, sentinel)
    process_vector_block(
        paths, bin_factor=bin_factor, h=H, w=W, output_dir=out_dir, idx=0
    )
    # File should still be the sentinel (not overwritten)
    result = np.load(out_file)
    np.testing.assert_array_equal(result, sentinel)


def test_process_vector_block_output_dtype(tmp_path: Path):
    H, W = 8, 8
    bin_factor = 2
    paths = _write_npy_slices(tmp_path / "input", n=bin_factor, shape=(3, H, W))
    out_dir = tmp_path / "output"
    process_vector_block(
        paths, bin_factor=bin_factor, h=H, w=W, output_dir=out_dir, idx=3
    )
    out_file = out_dir / "eigen_vec" / "eigen_vec_000003.npy"
    result = np.load(out_file)
    assert result.dtype == np.float32

