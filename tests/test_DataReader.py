import tempfile
from pathlib import Path
import numpy as np
import cv2
import pytest

from cardiotensor.utils.DataReader import DataReader

@pytest.fixture
def test_stack(tmp_path: Path):
    # Create a fake image stack
    for i in range(5):
        img = np.full((10, 10), i*50, dtype=np.uint8)
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
