import warnings
from pathlib import Path

import numpy as np
import pytest
import tifffile

warnings.filterwarnings("ignore", category=DeprecationWarning)

from cardiotensor.export.amira_writer import amira_writer, write_am_file


@pytest.fixture
def mock_data(tmp_path):
    """
    Create a mock vector field (3D numpy array) for testing.
    """
    data_dir = tmp_path / "data/volume"
    data_dir.mkdir(parents=True, exist_ok=True)
    for i in range(5):
        tifffile.imwrite(
            str(data_dir / f"slice_{i:06d}.tif"),
            np.random.rand(50, 50).astype(np.uint8),
        )
    return data_dir


@pytest.fixture
def mock_mask(tmp_path):
    """
    Create a mock vector field (3D numpy array) for testing.
    """
    mask_dir = tmp_path / "data/mask"
    mask_dir.mkdir(parents=True, exist_ok=True)
    for i in range(5):
        tifffile.imwrite(
            str(mask_dir / f"mask_{i:06d}.tif"),
            (np.random.rand(50, 50) * 255).astype(np.uint8),
        )
    return mask_dir


@pytest.fixture
def mock_vector_field(tmp_path):
    """
    Create a mock vector field (3D numpy array) for testing.
    """
    vector_dir = tmp_path / "output" / "eigen_vec"
    vector_dir.mkdir(parents=True, exist_ok=True)
    for i in range(5):
        np.save(vector_dir / f"vec_{i:06d}.npy", np.random.rand(3, 50, 50))
    return vector_dir


@pytest.fixture
def mock_helix_angles(tmp_path):
    """Create mock helix angle data as .tif files."""
    helix_dir = tmp_path / "output/HA"
    helix_dir.mkdir(parents=True, exist_ok=True)

    for i in range(5):  # Create 5 mock helix angle slices
        mock_data = np.random.rand(128, 128) * 180 - 90  # Range (-90, 90)
        tifffile.imwrite(str(helix_dir / f"HA_{i:06d}.tif"), mock_data.astype(np.uint8))

    return helix_dir


@pytest.fixture
def mock_configuration_file(tmp_path, mock_data, mock_mask, mock_vector_field):
    """
    Create a mock configuration file for testing.
    """
    conf_file = tmp_path / "parameters.conf"
    with open(conf_file, "w") as f:
        f.write(
            f"""
            [DATASET]
            IMAGES_PATH = {mock_data}
            VOXEL_SIZE = 1.0
            MASK_PATH = {mock_mask}
            FLIP = False

            [OUTPUT]
            OUTPUT_PATH = {mock_vector_field.parent}
            OUTPUT_FORMAT = tif
            OUTPUT_TYPE = 8bit
            VECTORS = True

            [RUN]
            TEST = False
            N_SLICE_TEST = 0
            """
        )
    return conf_file


def test_amira_writer_valid_input(
    mock_configuration_file,
    mock_data,
    mock_mask,
    mock_vector_field,
    mock_helix_angles,
    tmp_path,
):
    """
    Test amira_writer with valid input.
    """
    output_dir = tmp_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    amira_writer(
        conf_file_path=str(mock_configuration_file),
        start_index=0,
        end_index=5,
        bin_factor=None,
        num_ini_points=1000,
        num_steps=10,
        segment_length=5.0,
        angle_threshold=30.0,
        segment_min_length_threshold=3,
    )

    output_am_file = output_dir / "output.am"
    assert output_am_file.exists(), "Output AM file was not created."


def test_amira_writer_write_am_file():
    """
    Test the write_am_file function directly.
    """
    consecutive_points_list = [
        [(0, 0, 0), (10, 0, 0)],
        [(20, 0, 0), (30, 0, 0)],
    ]
    HA_angle = [45.0, 30.0]
    z_angle = [10.0, 20.0]
    output_file = "test_output.am"

    write_am_file(consecutive_points_list, HA_angle, z_angle, output_file)

    with open(output_file) as f:
        content = f.read()
        assert "AmiraMesh 3D ASCII 3.0" in content
        assert "define VERTEX 4" in content
        assert "define EDGE 2" in content

    Path(output_file).unlink()  # Clean up


def test_amira_writer_invalid_config_file(tmp_path):
    """
    Test amira_writer with an invalid configuration file.
    """
    invalid_conf = tmp_path / "invalid.conf"
    invalid_conf.write_text("[INVALID_SECTION]\nINVALID_KEY=VALUE\n")

    with pytest.raises(SystemExit) as exc_info:
        amira_writer(conf_file_path=str(invalid_conf))
    assert "Error reading parameter file" in str(exc_info.value)
