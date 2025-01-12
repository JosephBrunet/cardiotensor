import pytest
import numpy as np
import cv2
from pathlib import Path
from unittest.mock import patch, MagicMock
from cardiotensor.utils.downsampling import (
    downsample_vector_volume,
    downsample_volume,
    process_vector_block,
    process_image_block,
)


@pytest.fixture
def mock_vector_files(tmp_path):
    """
    Create mock .npy files for testing vector downsampling.
    """
    vector_dir = tmp_path / "vectors"
    vector_dir.mkdir()
    for i in range(10):  # Create 10 mock files
        np.save(vector_dir / f"eigen_vec_{i:06d}.npy", np.random.rand(3, 100, 100))
    return vector_dir


@pytest.fixture
def mock_image_files(tmp_path):
    """
    Create mock image files for testing image downsampling.
    """
    image_dir = tmp_path / "images"
    image_dir.mkdir()
    for i in range(10):  # Create 10 mock files
        img = (np.random.rand(100, 100) * 255).astype(np.uint8)
        cv2.imwrite(str(image_dir / f"HA_{i:06d}.tif"), img)
    return image_dir



def test_process_vector_block(tmp_path, mock_vector_files):
    """
    Test the process_vector_block function.
    """
    bin_factor = 2
    output_dir = tmp_path / f"output/bin{bin_factor}"
    output_dir.mkdir(parents=True, exist_ok=True)
    block = sorted(mock_vector_files.glob("*.npy"))[0:bin_factor] 
    
    process_vector_block(
        block=block,
        bin_factor=bin_factor,
        h=100,
        w=100,
        output_dir=output_dir,
        idx=0,
    )

    # Check that the downsampled file exists
    output_file = output_dir / "eigen_vec/eigen_vec_000000.npy"
    assert output_file.exists()

    # Verify the shape of the output
    data = np.load(output_file)
    assert data.shape == (3, 50, 50)  # Original (100x100) downsampled by a factor of 2


def test_downsample_vector_volume(tmp_path, mock_vector_files):
    """
    Test the downsample_vector_volume function.
    """
    output_dir = tmp_path / "output"
    downsample_vector_volume(mock_vector_files, bin_factor=2, output_dir=output_dir)

    # Check that output directory exists
    assert (output_dir / "bin2/eigen_vec").exists()

    # Verify that all blocks are processed
    files = list((output_dir / "bin2/eigen_vec").glob("*.npy"))
    assert len(files) == 5  # 10 original files processed in blocks of 2


def test_process_image_block(tmp_path, mock_image_files):
    """
    Test the process_image_block function.
    """
    bin_factor = 2
    output_dir = tmp_path / f"output/bin{bin_factor}"
    output_dir.mkdir(parents=True, exist_ok=True)
    block = sorted(mock_image_files.glob("*.tif"))[0:bin_factor] 
    process_image_block(
        block=block,
        bin_factor=2,
        h=100,
        w=100,
        output_dir=output_dir,
        idx=0,
    )

    # Check that the downsampled file exists
    output_file = output_dir / "HA/HA_000000.tif"
    assert output_file.exists()

    # Verify the shape of the output
    img = cv2.imread(str(output_file), cv2.IMREAD_UNCHANGED)
    assert img.shape == (50, 50)  # Original (100x100) downsampled by a factor of 2


def test_downsample_volume(tmp_path, mock_image_files):
    """
    Test the downsample_volume function.
    """
    output_dir = tmp_path / "output"
    downsample_volume(
        input_path=mock_image_files,
        bin_factor=2,
        output_dir=output_dir,
        file_format="tif",
    )

    # Check that output directory exists
    assert (output_dir / "bin2/HA").exists()

    # Verify that all blocks are processed
    files = list((output_dir / "bin2/HA").glob("*.tif"))
    assert len(files) == 5  # 10 original files processed in blocks of 2
