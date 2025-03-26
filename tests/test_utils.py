import os
import tempfile

import numpy as np

from cardiotensor.utils.utils import convert_to_8bit, read_conf_file


def test_read_conf_file():
    """Test the read_conf_file function."""
    conf_content = """
    [DATASET]
    IMAGES_PATH = /path/to/images
    MASK_PATH = /path/to/masks
    VOXEL_SIZE = 0.5

    [OUTPUT]
    OUTPUT_PATH = /path/to/output
    OUTPUT_FORMAT = jp2
    OUTPUT_TYPE = volume

    [STRUCTURE TENSOR CALCULATION]
    SIGMA = 1.5
    RHO = 2.0
    N_CHUNK = 50
    USE_GPU = True
    WRITE_VECTORS = True
    REVERSE = True

    [ANGLE CALCULATION]
    WRITE_ANGLES = True
    AXIS_POINTS = (10,20,30), (40,50,60)

    [TEST]
    TEST = True
    N_SLICE_TEST = 10
    """
    with tempfile.NamedTemporaryFile(suffix=".conf", delete=False) as temp_file:
        temp_file.write(conf_content.encode())
        temp_file_path = temp_file.name

    try:
        config = read_conf_file(temp_file_path)

        assert config["IMAGES_PATH"] == "/path/to/images"
        assert config["MASK_PATH"] == "/path/to/masks"
        assert config["VOXEL_SIZE"] == 0.5

        assert config["OUTPUT_PATH"] == "/path/to/output"
        assert config["OUTPUT_FORMAT"] == "jp2"
        assert config["OUTPUT_TYPE"] == "volume"

        assert config["SIGMA"] == 1.5
        assert config["RHO"] == 2.0
        assert config["N_CHUNK"] == 50
        assert config["USE_GPU"] is True
        assert config["WRITE_VECTORS"] is True
        assert config["REVERSE"] is True

        assert config["WRITE_ANGLES"] is True
        assert config["AXIS_POINTS"] == [(10, 20, 30), (40, 50, 60)]

        assert config["TEST"] is True
        assert config["N_SLICE_TEST"] == 10

        print("✅ read_conf_file test passed.")
    finally:
        os.remove(temp_file_path)


def test_convert_to_8bit():
    """Test the convert_to_8bit function."""
    img = np.array([[0, 50, 100], [150, 200, 250]], dtype=np.float32)

    # Test default behavior
    img_8bit = convert_to_8bit(img)
    assert img_8bit.min() == 0
    assert img_8bit.max() == 255
    print("✅ convert_to_8bit default test passed.")

    # Test with percentiles
    img_8bit = convert_to_8bit(img, perc_min=0, perc_max=100)
    assert img_8bit.min() == 0
    assert img_8bit.max() == 255
    print("✅ convert_to_8bit percentile test passed.")

    # Test with explicit min/max
    img_8bit = convert_to_8bit(img, min_value=50, max_value=200)
    assert img_8bit.min() == 0
    assert img_8bit.max() == 255
    print("✅ convert_to_8bit explicit range test passed.")


if __name__ == "__main__":
    test_read_conf_file()
    test_convert_to_8bit()
