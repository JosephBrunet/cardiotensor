from __future__ import annotations

import math
import importlib
from pathlib import Path

import numpy as np
import pytest

from cardiotensor.utils.downsampling import (
    _downsample_vector_axes,
    _process_vector_range,
    downsample_vector_volume,
    downsample_volume,
    process_vector_block,
)
from cardiotensor.utils.image_io import initialize_zarr_vector_field

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


def test_downsample_vector_axes_does_not_cancel_opposite_signs():
    vectors = np.zeros((3, 2, 2, 2), dtype=np.float32)
    vectors[0, 0] = 1.0
    vectors[0, 1] = -1.0

    result = _downsample_vector_axes(vectors, bin_factor=2)

    assert result.shape == (3, 1, 1)
    np.testing.assert_allclose(result[:, 0, 0], [1.0, 0.0, 0.0], atol=1e-6)
    np.testing.assert_allclose(np.linalg.norm(result[:, 0, 0]), 1.0, atol=1e-6)


def test_downsample_vector_volume_reads_zarr(tmp_path: Path):
    store = initialize_zarr_vector_field(tmp_path / "source", (2, 4, 4))
    vector_slice = np.zeros((3, 4, 4), dtype=np.float32)
    vector_slice[0] = 1.0
    store.write_slice(vector_slice, 0)
    store.write_slice(-vector_slice, 1)

    output_path = downsample_vector_volume(
        store.path, bin_factor=2, output_dir=tmp_path / "downsampled"
    )

    result = np.load(output_path / "eigen_vec_000000.npy")
    assert result.shape == (3, 2, 2)
    np.testing.assert_allclose(result[0], 1.0, atol=1e-6)
    np.testing.assert_allclose(result[1:], 0.0, atol=1e-6)


def test_downsample_vector_volume_applies_mask_before_binning(tmp_path: Path):
    store = initialize_zarr_vector_field(tmp_path / "source", (2, 4, 4))
    vector_slice = np.zeros((3, 4, 4), dtype=np.float32)
    vector_slice[1] = 1.0
    vector_slice[0, :2, :2] = 1.0
    vector_slice[1, :2, :2] = 0.0
    store.write_slice(vector_slice, 0)
    store.write_slice(vector_slice, 1)

    mask_dir = tmp_path / "myocardium_mask"
    mask_dir.mkdir()
    mask_slice = np.zeros((2, 2), dtype=np.uint8)
    mask_slice[0, 0] = 1
    np.save(mask_dir / "mask_000000.npy", mask_slice)

    output_dir = tmp_path / "downsampled"
    output_path = downsample_vector_volume(store.path, 2, output_dir)
    unmasked = np.load(output_path / "eigen_vec_000000.npy")
    np.testing.assert_allclose(unmasked[:, 1, 1], [0.0, 1.0, 0.0], atol=1e-6)

    output_path = downsample_vector_volume(
        store.path, 2, output_dir, mask_path=mask_dir
    )
    masked = np.load(output_path / "eigen_vec_000000.npy")
    np.testing.assert_allclose(masked[:, 0, 0], [1.0, 0.0, 0.0], atol=1e-6)
    np.testing.assert_allclose(masked[:, 1, 1], 0.0, atol=1e-6)


def test_empty_mask_block_skips_vector_read(tmp_path: Path):
    class VectorReader:
        shape = (3, 2, 4, 4)

        def load_volume(self, *args, **kwargs):
            raise AssertionError("empty mask block should not read vectors")

    class MaskReader:
        shape = (1, 2, 2)

        def load_volume(self, *args, **kwargs):
            return np.zeros((1, 2, 2), dtype=np.uint8)

    _process_vector_range(
        VectorReader(),
        start_index=0,
        end_index=2,
        bin_factor=2,
        output_dir=tmp_path,
        output_index=0,
        mask_reader=MaskReader(),
        mask_y_indices=np.array([0, 0, 1, 1]),
        mask_x_indices=np.array([0, 0, 1, 1]),
    )

    result = np.load(tmp_path / "eigen_vec_000000.npy")
    np.testing.assert_array_equal(result, np.zeros((3, 2, 2), dtype=np.float32))


def test_streamline_seed_selection_uses_mask(tmp_path: Path, monkeypatch):
    tractography = importlib.import_module(
        "cardiotensor.tractography.generate_streamlines"
    )
    vector_path = tmp_path / "vectors"
    fa_path = tmp_path / "FA"
    angle_path = tmp_path / "HA"
    mask_path = tmp_path / "mask"
    for path in (vector_path, fa_path, angle_path, mask_path):
        path.mkdir()

    myocardium = np.zeros((2, 2, 2), dtype=np.uint8)
    myocardium[:, 0, 0] = 1

    class FakeReader:
        def __init__(self, path):
            self.path = Path(path)
            self.shape = (3, 2, 2, 2) if self.path == vector_path else (2, 2, 2)

        def load_region(self, **kwargs):
            return np.ones((3, 2, 2, 2), dtype=np.float32)

        def load_volume(self, **kwargs):
            if self.path == mask_path:
                return myocardium.copy()
            return np.ones((2, 2, 2), dtype=np.float32)

    captured = {}

    def capture_seed_mask(seed_mask, num_seeds, random_seed):
        captured["seed_mask"] = seed_mask.copy()
        raise RuntimeError("seed mask captured")

    monkeypatch.setattr(tractography, "DataReader", FakeReader)
    monkeypatch.setattr(tractography, "_select_seed_points", capture_seed_mask)

    with pytest.raises(RuntimeError, match="seed mask captured"):
        tractography.generate_streamlines_from_params(
            vector_field_dir=vector_path,
            output_dir=tmp_path / "output",
            fa_dir=fa_path,
            angle_dir=angle_path,
            mask_path=mask_path,
            fa_seed_min=0.5,
        )

    np.testing.assert_array_equal(captured["seed_mask"], myocardium.astype(bool))


def test_streamline_generation_rejects_incomplete_scalar_outputs(
    tmp_path: Path, monkeypatch
):
    tractography = importlib.import_module(
        "cardiotensor.tractography.generate_streamlines"
    )
    vector_path = tmp_path / "vectors"
    fa_path = tmp_path / "FA"
    angle_path = tmp_path / "HA"
    for path in (vector_path, fa_path, angle_path):
        path.mkdir()

    class FakeReader:
        def __init__(self, path):
            self.path = Path(path)
            self.volume_info = {"type": "npy"}
            self.shape = (3, 8, 4, 4) if self.path == vector_path else (7, 4, 4)

    monkeypatch.setattr(tractography, "DataReader", FakeReader)

    with pytest.raises(RuntimeError, match="Orientation outputs are incomplete"):
        tractography.generate_streamlines_from_params(
            vector_field_dir=vector_path,
            output_dir=tmp_path / "output",
            fa_dir=fa_path,
            angle_dir=angle_path,
        )


def test_scalar_downsampling_rebuilds_unmarked_cache(tmp_path: Path):
    source = tmp_path / "source"
    source.mkdir()
    for index in range(4):
        np.save(source / f"FA_{index:06d}.npy", np.full((4, 4), 255, np.uint8))

    stale_dir = tmp_path / "output" / "bin2" / "FA"
    stale_dir.mkdir(parents=True)
    stale_path = stale_dir / "FA_000000.tif"
    stale_path.write_bytes(b"stale")

    downsample_volume(source, 2, tmp_path / "output", subfolder="FA", out_ext="tif")

    assert stale_path.stat().st_size > len(b"stale")
    assert (stale_dir / "FA_000001.tif").exists()
    assert (stale_dir / ".source_complete").exists()
