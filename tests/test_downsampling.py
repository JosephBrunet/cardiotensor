import cv2
import numpy as np
import pytest

from cardiotensor.utils.downsampling import (
    downsample_vector_volume,
    downsample_volume,
    process_image_block,
    process_vector_block,
)


@pytest.fixture
def mock_vector_files(tmp_path):
    """
    Create mock .npy files for testing vector downsampling.
    """
    vector_dir = tmp_path / "vectors"
    vector_dir.mkdir()
    for i in range(10):
        np.save(vector_dir / f"eigen_vec_{i:06d}.npy", np.random.rand(3, 100, 100))
    return vector_dir


@pytest.fixture
def mock_image_files(tmp_path):
    """
    Create mock image files for testing image downsampling.
    """
    image_dir = tmp_path / "images"
    image_dir.mkdir()
    for i in range(10):
        img = (np.random.rand(100, 100) * 255).astype(np.uint8)
        cv2.imwrite(str(image_dir / f"HA_{i:06d}.tif"), img)
    return image_dir


def test_process_vector_block(tmp_path, mock_vector_files):
    """
    Test the process_vector_block function.
    """
    bin_factor = 2
    output_dir = tmp_path / "output"

    block = sorted(mock_vector_files.glob("*.npy"))[:bin_factor]

    process_vector_block(
        block=block,
        bin_factor=bin_factor,
        h=100,
        w=100,
        output_dir=output_dir,
        idx=0,
    )

    # Output file should be created inside output_dir/bin2/eigen_vec/
    output_file = output_dir / "eigen_vec/eigen_vec_000000.npy"
    assert output_file.exists()

    data = np.load(output_file)
    assert data.shape == (3, 50, 50)


def test_downsample_vector_volume(tmp_path, mock_vector_files):
    """
    Test the downsample_vector_volume function.
    """
    output_dir = tmp_path / "output"
    downsample_vector_volume(mock_vector_files, bin_factor=2, output_dir=output_dir)

    files = list((output_dir / "bin2/eigen_vec").glob("*.npy"))
    assert len(files) == 5  # 10 original files processed in blocks of 2


def test_process_image_block(tmp_path, mock_image_files):
    """
    Test the process_image_block function.
    """
    bin_factor = 2
    output_dir = tmp_path / "output"

    block = sorted(mock_image_files.glob("*.tif"))[:bin_factor]

    process_image_block(
        block=block,
        bin_factor=bin_factor,
        h=100,
        w=100,
        output_dir=output_dir,
        idx=0,
    )

    output_file = output_dir / "HA/HA_000000.tif"
    assert output_file.exists()

    img = cv2.imread(str(output_file), cv2.IMREAD_UNCHANGED)
    assert img.shape == (50, 50)


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

    files = list((output_dir / "bin2/HA").glob("*.tif"))
    assert len(files) == 5
