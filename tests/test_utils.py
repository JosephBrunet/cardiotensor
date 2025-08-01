import numpy as np
import pytest
from pathlib import Path

from cardiotensor.utils.utils import convert_to_8bit, read_conf_file

def test_convert_to_8bit():
    arr = np.array([0, 128, 255], dtype=np.uint16)
    out = convert_to_8bit(arr)
    assert out.dtype == np.uint8
    assert out.min() >= 0 and out.max() <= 255

def test_read_conf_file(tmp_path):
    # Create dummy directories for IMAGES_PATH and MASK_PATH
    images_dir = tmp_path / "images"
    images_dir.mkdir()
    mask_file = tmp_path / "mask.tif"
    mask_file.write_text("dummy")

    # Create a dummy .conf file
    conf_file = tmp_path / "test.conf"
    conf_file.write_text(
        "[DATASET]\n"
        f"IMAGES_PATH = {images_dir}\n"
        f"MASK_PATH = {mask_file}\n"
        "VOXEL_SIZE = 0.5\n\n"
        "[STRUCTURE TENSOR CALCULATION]\n"
        "SIGMA = 2.0\n"
        "RHO = 1.5\n"
        "TRUNCATE = 4.0\n"
        "VERTICAL_PADDING = 10.0\n"
        "N_CHUNK = 50\n"
        "USE_GPU = True\n"
        "WRITE_VECTORS = True\n"
        "REVERSE = False\n\n"
        "[ANGLE CALCULATION]\n"
        "WRITE_ANGLES = True\n"
        "AXIS_POINTS = (0,0,0), (1,1,1)\n\n"
        "[TEST]\n"
        "TEST = True\n"
        "N_SLICE_TEST = 5\n\n"
        "[OUTPUT]\n"
        "OUTPUT_PATH = ./output\n"
        "OUTPUT_FORMAT = tif\n"
        "OUTPUT_TYPE = rgb\n"
    )

    config = read_conf_file(conf_file)

    # --- Assertions ---
    # Dataset paths
    assert config["IMAGES_PATH"] == str(images_dir)
    assert config["MASK_PATH"] == str(mask_file)
    assert pytest.approx(config["VOXEL_SIZE"]) == 0.5

    # Structure tensor
    assert pytest.approx(config["SIGMA"]) == 2.0
    assert pytest.approx(config["RHO"]) == 1.5
    assert config["N_CHUNK"] == 50
    assert config["USE_GPU"] is True
    assert config["WRITE_VECTORS"] is True
    assert config["REVERSE"] is False

    # Angles
    assert config["WRITE_ANGLES"] is True
    assert isinstance(config["AXIS_POINTS"], list)
    assert config["AXIS_POINTS"] == [(0,0,0), (1,1,1)]

    # Test
    assert config["TEST"] is True
    assert config["N_SLICE_TEST"] == 5

    # Output
    assert config["OUTPUT_PATH"] == "./output"
    assert config["OUTPUT_FORMAT"] == "tif"
    assert config["OUTPUT_TYPE"] == "rgb"