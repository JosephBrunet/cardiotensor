from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cardiotensor.utils.am_utils import (
    _normalize_edge_attribute,
    _normalize_point_attribute,
    _sanitize_field_name,
    write_spatialgraph_am,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _two_streamlines():
    sl0 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=float)
    sl1 = np.array([[0, 1, 0], [0, 2, 0]], dtype=float)
    return [sl0, sl1]


# ---------------------------------------------------------------------------
# _sanitize_field_name
# ---------------------------------------------------------------------------

def test_sanitize_field_name_simple():
    assert _sanitize_field_name("HA", "p") == "HA"


def test_sanitize_field_name_replaces_spaces():
    result = _sanitize_field_name("my field", "p")
    assert " " not in result
    assert "_" in result


def test_sanitize_field_name_empty_raises():
    with pytest.raises(ValueError):
        _sanitize_field_name("", "p")


def test_sanitize_field_name_whitespace_only_raises():
    with pytest.raises(ValueError):
        _sanitize_field_name("   ", "p")


def test_sanitize_field_name_non_string_raises():
    with pytest.raises(ValueError):
        _sanitize_field_name(123, "p")  # type: ignore


# ---------------------------------------------------------------------------
# _normalize_edge_attribute
# ---------------------------------------------------------------------------

def test_normalize_edge_attribute_list():
    arr = _normalize_edge_attribute([1.0, 2.0, 3.0], n_edges=3, name="x")
    assert arr.shape == (3,)
    np.testing.assert_allclose(arr, [1.0, 2.0, 3.0])


def test_normalize_edge_attribute_ndarray():
    arr = _normalize_edge_attribute(np.array([0.5, 1.5]), n_edges=2, name="x")
    assert arr.dtype == float


def test_normalize_edge_attribute_wrong_length():
    with pytest.raises(ValueError, match="n_edges"):
        _normalize_edge_attribute([1.0, 2.0], n_edges=3, name="x")


# ---------------------------------------------------------------------------
# _normalize_point_attribute
# ---------------------------------------------------------------------------

def test_normalize_point_attribute_flat_array():
    num_pts = [3, 2]
    attr = np.ones(5, dtype=float)
    result = _normalize_point_attribute(attr, num_pts, n_points_total=5, name="t")
    assert result.shape == (5,)


def test_normalize_point_attribute_list_of_arrays():
    num_pts = [3, 2]
    attr = [np.ones(3), np.ones(2)]
    result = _normalize_point_attribute(attr, num_pts, n_points_total=5, name="t")
    assert result.shape == (5,)


def test_normalize_point_attribute_flat_wrong_length():
    with pytest.raises(ValueError, match="total points"):
        _normalize_point_attribute(np.ones(4), [3, 2], n_points_total=5, name="t")


def test_normalize_point_attribute_list_wrong_count():
    with pytest.raises(ValueError, match="one array per streamline"):
        _normalize_point_attribute([np.ones(3)], [3, 2], n_points_total=5, name="t")


def test_normalize_point_attribute_list_wrong_lengths():
    with pytest.raises(ValueError, match="streamline length"):
        _normalize_point_attribute([np.ones(2), np.ones(2)], [3, 2], n_points_total=5, name="t")


# ---------------------------------------------------------------------------
# write_spatialgraph_am
# ---------------------------------------------------------------------------

def test_write_basic(tmp_path: Path):
    out = tmp_path / "out.am"
    write_spatialgraph_am(out, _two_streamlines())
    assert out.exists()
    content = out.read_text()
    assert "AmiraMesh" in content
    assert "@1" in content
    assert "VERTEX" in content
    assert "EDGE" in content


def test_write_creates_parent_dirs(tmp_path: Path):
    out = tmp_path / "deep" / "nested" / "out.am"
    write_spatialgraph_am(out, _two_streamlines())
    assert out.exists()


def test_write_empty_raises():
    with pytest.raises(ValueError, match="No streamlines"):
        write_spatialgraph_am(Path("/dev/null"), [])


def test_write_short_streamline_raises(tmp_path: Path):
    sl = np.array([[0, 0, 0]], dtype=float)  # only 1 point
    with pytest.raises(ValueError, match="at least 2 points"):
        write_spatialgraph_am(tmp_path / "out.am", [sl])


def test_write_with_thickness_flat(tmp_path: Path):
    sls = _two_streamlines()
    thickness = np.ones(5, dtype=float)  # 3 + 2 = 5 points total
    out = tmp_path / "out.am"
    write_spatialgraph_am(out, sls, point_thickness=thickness)
    assert "@5" in out.read_text()


def test_write_with_thickness_list(tmp_path: Path):
    sls = _two_streamlines()
    thickness = [np.ones(3), np.ones(2)]
    out = tmp_path / "out.am"
    write_spatialgraph_am(out, sls, point_thickness=thickness)
    assert "@5" in out.read_text()


def test_write_with_single_edge_scalar(tmp_path: Path):
    sls = _two_streamlines()
    out = tmp_path / "out.am"
    write_spatialgraph_am(
        out, sls, edge_scalar=[0.1, 0.2], edge_scalar_name="HA"
    )
    content = out.read_text()
    assert "HA" in content
    assert "@6" in content


def test_write_single_edge_scalar_without_name_raises(tmp_path: Path):
    sls = _two_streamlines()
    with pytest.raises(ValueError, match="edge_scalar_name"):
        write_spatialgraph_am(tmp_path / "out.am", sls, edge_scalar=[0.1, 0.2])


def test_write_with_dict_edge_scalars(tmp_path: Path):
    sls = _two_streamlines()
    out = tmp_path / "out.am"
    write_spatialgraph_am(
        out,
        sls,
        edge_scalar={"HA": [10.0, 20.0], "IA": [-5.0, 5.0]},
    )
    content = out.read_text()
    assert "HA" in content
    assert "IA" in content
    assert "@6" in content
    assert "@7" in content


def test_written_vertex_count(tmp_path: Path):
    sls = _two_streamlines()
    out = tmp_path / "out.am"
    write_spatialgraph_am(out, sls)
    content = out.read_text()
    # 2 streamlines → 4 vertices (start+end for each)
    assert "define VERTEX 4" in content


def test_written_edge_count(tmp_path: Path):
    sls = _two_streamlines()
    out = tmp_path / "out.am"
    write_spatialgraph_am(out, sls)
    assert "define EDGE 2" in out.read_text()
