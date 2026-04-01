from pathlib import Path

import cv2
import numpy as np
import pytest
import tifffile

from cardiotensor.utils.DataReader import DataReader


@pytest.fixture
def test_stack(tmp_path: Path):
    # Create a fake image stack
    for i in range(5):
        img = np.full((10, 10), i * 50, dtype=np.uint8)
        cv2.imwrite(str(tmp_path / f"img_{i:03d}.tif"), img)
    return tmp_path


def test_read_volume_shape(test_stack):
    reader = DataReader(test_stack)
    vol = reader.load_volume()
    assert vol.shape[0] == 5
    assert vol.ndim == 3


def test_mask_optional(test_stack):
    reader = DataReader(test_stack)
    # Should not fail without mask
    vol = reader.load_volume()
    assert isinstance(vol, np.ndarray)


def test_read_single_tiff_stack(tmp_path: Path):
    stack = np.stack(
        [np.full((10, 10), i * 50, dtype=np.uint8) for i in range(5)], axis=0
    )
    stack_path = tmp_path / "stack.tif"
    tifffile.imwrite(stack_path, stack)

    reader = DataReader(stack_path)
    vol = reader.load_volume()

    assert vol.shape == (5, 10, 10)
    assert np.array_equal(vol[2], stack[2])
