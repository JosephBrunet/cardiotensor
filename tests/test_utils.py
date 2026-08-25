from pathlib import Path
from types import SimpleNamespace
from unittest import TestCase

import numpy as np
import pytest

from cardiotensor.utils import utils
from cardiotensor.utils.utils import (
    _memory_available_from_cgroup,
    convert_to_8bit,
    get_available_memory_bytes,
    read_conf_file,
)


def test_memory_available_from_cgroup_v2(tmp_path: Path):
    proc_cgroup = tmp_path / "cgroup"
    proc_cgroup.write_text("0::/slurm/job_123\n")
    job_cgroup = tmp_path / "sys_fs_cgroup" / "slurm" / "job_123"
    job_cgroup.mkdir(parents=True)
    (job_cgroup / "memory.max").write_text("1000")
    (job_cgroup / "memory.current").write_text("250")

    assert _memory_available_from_cgroup(proc_cgroup, tmp_path / "sys_fs_cgroup") == 750

    (job_cgroup / "memory.current").write_text("1000")
    assert _memory_available_from_cgroup(proc_cgroup, tmp_path / "sys_fs_cgroup") == 0


def test_available_memory_uses_cgroup_limit(monkeypatch):
    monkeypatch.setattr(
        utils.psutil, "virtual_memory", lambda: SimpleNamespace(available=1000)
    )
    monkeypatch.setattr(utils, "_memory_available_from_cgroup", lambda: 400)

    assert get_available_memory_bytes() == 400


def test_convert_to_8bit():
    arr = np.array([0, 128, 255], dtype=np.uint16)
    out = convert_to_8bit(arr)
    assert out.dtype == np.uint8
    assert out.min() >= 0 and out.max() <= 255


def test_convert_to_8bit_skips_percentiles_for_explicit_range(monkeypatch):
    def fail_if_called(*args, **kwargs):
        raise AssertionError("nanpercentile should not be called")

    monkeypatch.setattr(np, "nanpercentile", fail_if_called)
    out = convert_to_8bit(
        np.array([0.0, 0.5, 1.0], dtype=np.float32),
        min_value=0.0,
        max_value=1.0,
    )
    np.testing.assert_array_equal(out, [0, 127, 255])


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
        "IMAGES_PATH = images\n"
        "MASK_PATH = mask.tif\n"
        "VOXEL_SIZE = 0.5\n\n"
        "[STRUCTURE TENSOR CALCULATION]\n"
        "SIGMA = 2.0\n"
        "RHO = 1.5\n"
        "TRUNCATE = 4.0\n"
        "VERTICAL_PADDING = 10.0\n"
        "N_CHUNK = 50\n"
        "USE_GPU = True\n"
        "GPU_WORKERS_PER_DEVICE = 2\n"
        "LOW_MEMORY = True\n"
        "LOW_MEMORY_DIR = scratch\n"
        "WRITE_VECTORS = True\n"
        "REVERSE = False\n\n"
        "[ANGLE CALCULATION]\n"
        "WRITE_ANGLES = True\n"
        "PROJECTED_ANGLES = True\n"
        "AXIS_POINTS = (0,0,0), (1,1,1)\n\n"
        "[TEST]\n"
        "TEST = True\n"
        "N_SLICE_TEST = 5\n"
        "SHOW_QUIVER = True\n\n"
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
    assert config["GPU_WORKERS_PER_DEVICE"] == 2
    assert config["LOW_MEMORY"] is True
    assert config["LOW_MEMORY_DIR"] == str(tmp_path / "scratch")
    assert config["WRITE_VECTORS"] is True
    assert config["VECTOR_FORMAT"] == "zarr"
    assert config["REVERSE"] is False

    # Angles
    assert config["WRITE_ANGLES"] is True
    TestCase().assertIs(config["PROJECTED_ANGLES"], True)
    assert isinstance(config["AXIS_POINTS"], list)
    assert config["AXIS_POINTS"] == [(0, 0, 0), (1, 1, 1)]

    # Test
    assert config["TEST"] is True
    assert config["N_SLICE_TEST"] == 5
    assert config["SHOW_QUIVER"] is True

    # Output
    assert config["OUTPUT_PATH"] == str(tmp_path / "output")
    assert config["OUTPUT_FORMAT"] == "tif"
    assert config["OUTPUT_TYPE"] == "rgb"


def test_read_conf_file_rejects_no_requested_outputs(tmp_path):
    images_dir = tmp_path / "images"
    images_dir.mkdir()
    conf_file = tmp_path / "invalid.conf"
    conf_file.write_text(
        "[DATASET]\n"
        "IMAGES_PATH = images\n\n"
        "[STRUCTURE TENSOR CALCULATION]\n"
        "WRITE_VECTORS = False\n\n"
        "[ANGLE CALCULATION]\n"
        "WRITE_ANGLES = False\n"
    )

    with pytest.raises(ValueError, match="At least one"):
        read_conf_file(conf_file)
