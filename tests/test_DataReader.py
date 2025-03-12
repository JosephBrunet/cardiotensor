import tempfile
import warnings
from pathlib import Path

import cv2
import numpy as np

warnings.simplefilter("ignore", category=DeprecationWarning)

from cardiotensor.utils.DataReader import DataReader


def create_test_stack(directory: Path, num_images: int, shape: tuple[int, int]):
    """
    Create a test image stack in the given directory.

    Args:
        directory (Path): Directory to create the images in.
        num_images (int): Number of images to create.
        shape (tuple[int, int]): Shape of each image (height, width).
    """
    directory.mkdir(exist_ok=True)
    for i in range(num_images):
        img = (np.random.rand(*shape) * 255).astype(np.uint8)
        img_path = directory / f"image_{i:03d}.tif"
        cv2.imwrite(str(img_path), img)


def create_test_mhd(
    directory: Path, shape: tuple[int, int, int], element_type="MET_UCHAR"
):
    """
    Create a test .mhd file and its corresponding raw data file.

    Args:
        directory (Path): Directory to create the .mhd and .raw files in.
        shape (tuple[int, int, int]): Shape of the volume (z, y, x).
        element_type (str): The data type for the .mhd file.
    """
    directory.mkdir(exist_ok=True)
    volume = (np.random.rand(*shape) * 255).astype(np.uint8)
    raw_file = directory / "test.raw"
    volume.tofile(raw_file)

    mhd_content = f"""ObjectType = Image
NDims = 3
DimSize = {" ".join(map(str, shape))}
ElementType = {element_type}
ElementSpacing = 1 1 1
ElementDataFile = {raw_file.name}
"""
    mhd_file = directory / "test.mhd"
    mhd_file.write_text(mhd_content)


def test_datareader_with_stack():
    """Test DataReader with an image stack."""
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        stack_dir = temp_dir / "stack"
        create_test_stack(stack_dir, num_images=10, shape=(128, 128))

        # Initialize DataReader and validate volume info
        reader = DataReader(stack_dir)
        assert reader.volume_info["type"] == "tif"
        assert len(reader.volume_info["file_list"]) == 10
        assert reader.shape == (10, 128, 128)

        # Load the volume
        volume = reader.load_volume()
        assert volume.shape == (10, 128, 128)
        print("Stack test passed.")


def test_datareader_with_mhd():
    """Test DataReader with an MHD file."""
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        mhd_dir = temp_dir / "mhd"
        create_test_mhd(mhd_dir, shape=(10, 128, 128))

        # Initialize DataReader and validate volume info
        reader = DataReader(mhd_dir / "test.mhd")
        assert reader.volume_info["type"] == "mhd"
        assert reader.shape == (10, 128, 128)

        # Load the volume
        volume = reader.load_volume()
        assert volume.shape == (10, 128, 128)
        print("MHD test passed.")


if __name__ == "__main__":
    test_datareader_with_stack()
    test_datareader_with_mhd()
