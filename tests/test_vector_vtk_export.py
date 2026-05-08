from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cardiotensor.utils.vector_vtk_export import (
    _writeFieldInVtk,
    export_vector_field_to_vtk,
    writeStructuredVTK,
)

# ---------------------------------------------------------------------------
# writeStructuredVTK
# ---------------------------------------------------------------------------


def test_writeStructuredVTK_scalar(tmp_path: Path):
    out = str(tmp_path / "out.vtk")
    data = np.ones((2, 3, 4), dtype=np.float32)
    writeStructuredVTK(cellData={"myScalar": data}, fileName=out)
    content = Path(out).read_text()
    assert "vtk DataFile" in content
    assert "SCALARS myScalar" in content
    assert "CELL_DATA" in content


def test_writeStructuredVTK_vector(tmp_path: Path):
    out = str(tmp_path / "out.vtk")
    data = np.ones((2, 3, 4, 3), dtype=np.float32)
    writeStructuredVTK(cellData={"vecs": data}, fileName=out)
    content = Path(out).read_text()
    assert "VECTORS vecs" in content


def test_writeStructuredVTK_no_data(tmp_path: Path, capsys):
    out = str(tmp_path / "out.vtk")
    writeStructuredVTK(fileName=out)
    captured = capsys.readouterr()
    assert "No data" in captured.out


def test_writeStructuredVTK_inconsistent_shapes(tmp_path: Path):
    out = str(tmp_path / "out.vtk")
    a = np.ones((2, 3, 4), dtype=np.float32)
    b = np.ones((2, 3, 5), dtype=np.float32)
    with pytest.raises(ValueError, match="Inconsistent shape"):
        writeStructuredVTK(cellData={"a": a, "b": b}, fileName=out)


def test_writeStructuredVTK_dimensions(tmp_path: Path):
    out = str(tmp_path / "out.vtk")
    data = np.zeros((2, 3, 4), dtype=np.float32)
    writeStructuredVTK(cellData={"x": data}, fileName=out)
    content = Path(out).read_text()
    # dimensions = shape + 1 in each axis
    assert "DIMENSIONS" in content


# ---------------------------------------------------------------------------
# _writeFieldInVtk
# ---------------------------------------------------------------------------


def test_writeFieldInVtk_scalar(tmp_path: Path):
    out = tmp_path / "field.txt"
    with open(out, "w") as f:
        _writeFieldInVtk({"s": np.ones((2, 2, 2))}, f)
    content = out.read_text()
    assert "SCALARS s float" in content
    assert "LOOKUP_TABLE default" in content


def test_writeFieldInVtk_vector(tmp_path: Path):
    out = tmp_path / "field.txt"
    data = np.ones((2, 2, 2, 3))
    with open(out, "w") as f:
        _writeFieldInVtk({"v": data}, f)
    content = out.read_text()
    assert "VECTORS v float" in content


def test_writeFieldInVtk_unsupported_shape(tmp_path: Path):
    out = tmp_path / "field.txt"
    with open(out, "w") as f:
        with pytest.raises(ValueError, match="Unsupported shape"):
            _writeFieldInVtk({"bad": np.ones((2, 3))}, f)


def test_writeFieldInVtk_flat_raises(tmp_path: Path):
    out = tmp_path / "field.txt"
    with open(out, "w") as f:
        with pytest.raises(NotImplementedError):
            _writeFieldInVtk({"x": np.ones((2, 2, 2))}, f, flat=True)


# ---------------------------------------------------------------------------
# export_vector_field_to_vtk
# ---------------------------------------------------------------------------


def test_export_vector_field_to_vtk_basic(tmp_path: Path):
    vf = np.random.rand(4, 4, 4, 3).astype(np.float32)
    color = np.random.rand(4, 4, 4).astype(np.float32)
    out = tmp_path / "paraview.vtk"
    result = export_vector_field_to_vtk(
        vf, color, voxel_size=1.0, downsample=2, save_path=out
    )
    assert result == out
    assert out.exists()
    assert out.stat().st_size > 0


def test_export_vector_field_default_path(tmp_path: Path):
    import os

    old_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        vf = np.ones((4, 4, 4, 3), dtype=np.float32)
        color = np.ones((4, 4, 4), dtype=np.float32)
        result = export_vector_field_to_vtk(vf, color, voxel_size=1.0, downsample=2)
        assert result.name == "paraview.vtk"
    finally:
        os.chdir(old_cwd)


def test_export_wrong_shape_raises():
    bad = np.ones((4, 4, 4), dtype=np.float32)  # missing last dim
    color = np.ones((4, 4, 4), dtype=np.float32)
    with pytest.raises(ValueError, match="shape"):
        export_vector_field_to_vtk(bad, color, voxel_size=1.0)


def test_export_negative_z_flipped(tmp_path: Path):
    vf = np.zeros((4, 4, 4, 3), dtype=np.float32)
    vf[..., 2] = -1.0  # all Z negative
    color = np.zeros((4, 4, 4), dtype=np.float32)
    out = tmp_path / "out.vtk"
    export_vector_field_to_vtk(vf, color, voxel_size=1.0, downsample=1, save_path=out)
    # After export the file should exist (flip happened in-place but no crash)
    assert out.exists()
