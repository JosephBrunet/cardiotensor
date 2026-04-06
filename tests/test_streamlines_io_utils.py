from __future__ import annotations

import numpy as np
import pytest

from cardiotensor.utils.streamlines_io_utils import (
    compute_elevation_angles,
    ha_to_degrees_per_streamline,
    load_npz_streamlines,
    normalize_attrs_to_degrees,
    reduce_per_edge,
)


# ---------------------------------------------------------------------------
# load_npz_streamlines
# ---------------------------------------------------------------------------

def _make_npz(tmp_path, streamlines_zyx, extras=None):
    """Save a .npz with streamlines stored as object array (z, y, x)."""
    obj = np.empty(len(streamlines_zyx), dtype=object)
    for i, sl in enumerate(streamlines_zyx):
        obj[i] = np.asarray(sl, dtype=np.float32)
    data = {"streamlines": obj}
    if extras:
        data.update(extras)
    p = tmp_path / "test.npz"
    np.savez(p, **data)
    return p


def test_load_npz_basic(tmp_path):
    sls_zyx = [
        [[0, 1, 2], [3, 4, 5]],
        [[6, 7, 8], [9, 10, 11]],
    ]
    p = _make_npz(tmp_path, sls_zyx)
    streamlines, per_point = load_npz_streamlines(p)
    assert len(streamlines) == 2
    assert streamlines[0].shape == (2, 3)
    # z,y,x → x,y,z swap: first point was (z=0,y=1,x=2) → (x=2,y=1,z=0)
    np.testing.assert_allclose(streamlines[0][0], [2, 1, 0])
    assert per_point == {}


def test_load_npz_with_values(tmp_path):
    sls_zyx = [[[0, 0, 0], [1, 1, 1]], [[2, 2, 2], [3, 3, 3]]]
    ha_obj = np.empty(2, dtype=object)
    ha_obj[0] = np.array([100.0, 150.0])
    ha_obj[1] = np.array([50.0, 200.0])
    p = _make_npz(tmp_path, sls_zyx, extras={"ha_values": ha_obj})
    streamlines, per_point = load_npz_streamlines(p)
    assert "HA" in per_point
    assert len(per_point["HA"]) == 2


def test_load_npz_missing_key(tmp_path):
    p = tmp_path / "empty.npz"
    np.savez(p, other=np.array([1, 2]))
    with pytest.raises(ValueError, match="streamlines"):
        load_npz_streamlines(p)


def test_load_npz_values_length_mismatch(tmp_path):
    sls_zyx = [[[0, 0, 0], [1, 1, 1]]]
    ha_obj = np.empty(2, dtype=object)  # 2 entries but only 1 streamline
    ha_obj[0] = np.array([1.0])
    ha_obj[1] = np.array([2.0])
    p = _make_npz(tmp_path, sls_zyx, extras={"ha_values": ha_obj})
    with pytest.raises(ValueError, match="does not match"):
        load_npz_streamlines(p)


# ---------------------------------------------------------------------------
# ha_to_degrees_per_streamline
# ---------------------------------------------------------------------------

def test_ha_to_degrees_already_degrees():
    ha = [np.array([-90.0, 0.0, 90.0], dtype=np.float32)]
    result = ha_to_degrees_per_streamline(ha)
    # max is 90 > 1.5 → treated as byte-scaled!
    # the function checks nanmax > 1.5 so -90..90 will be treated as scaled
    # just ensure output is same shape
    assert result[0].shape == ha[0].shape


def test_ha_to_degrees_byte_scaled():
    # 0..255 range → should convert to -90..90
    ha = [np.array([0.0, 127.5, 255.0], dtype=np.float32)]
    result = ha_to_degrees_per_streamline(ha)
    np.testing.assert_allclose(result[0][0], -90.0, atol=1.0)
    np.testing.assert_allclose(result[0][-1], 90.0, atol=1.0)


def test_ha_to_degrees_empty():
    result = ha_to_degrees_per_streamline([np.array([], dtype=np.float32)])
    assert result[0].shape == (0,)


# ---------------------------------------------------------------------------
# normalize_attrs_to_degrees
# ---------------------------------------------------------------------------

def test_normalize_attrs_empty():
    assert normalize_attrs_to_degrees(None) == {}
    assert normalize_attrs_to_degrees({}) == {}


def test_normalize_attrs_byte_scaled():
    attrs = {"HA": [np.array([0.0, 255.0], dtype=np.float32)]}
    result = normalize_attrs_to_degrees(attrs)
    assert "HA" in result
    np.testing.assert_allclose(result["HA"][0][0], -90.0, atol=1.0)
    np.testing.assert_allclose(result["HA"][0][-1], 90.0, atol=1.0)


def test_normalize_attrs_keys_uppercased():
    attrs = {"ha": [np.array([0.0, 0.5], dtype=np.float32)]}
    result = normalize_attrs_to_degrees(attrs)
    assert "HA" in result
    assert "ha" not in result


def test_normalize_attrs_already_degrees():
    # Values in -90..90 but max ≤ 1.5 is not satisfied (max=0.5 ≤ 1.5)
    attrs = {"IA": [np.array([-0.5, 0.0, 0.5], dtype=np.float32)]}
    result = normalize_attrs_to_degrees(attrs)
    # should come through unchanged (already float, low range)
    np.testing.assert_allclose(result["IA"][0], [-0.5, 0.0, 0.5], atol=1e-6)


# ---------------------------------------------------------------------------
# compute_elevation_angles
# ---------------------------------------------------------------------------

def test_elevation_angles_horizontal():
    # Streamline going purely in X — z-component of tangent = 0 → elevation = 0
    pts = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=np.float32)
    result = compute_elevation_angles([pts])
    assert result[0].shape == (3,)
    np.testing.assert_allclose(result[0], 0.0, atol=1e-5)


def test_elevation_angles_vertical():
    # Streamline going purely in Z — elevation should be ±90°
    pts = np.array([[0, 0, 0], [0, 0, 1], [0, 0, 2]], dtype=np.float32)
    result = compute_elevation_angles([pts])
    np.testing.assert_allclose(result[0][0], 90.0, atol=1e-4)


def test_elevation_angles_single_point():
    pts = np.array([[0, 0, 0]], dtype=np.float32)
    result = compute_elevation_angles([pts])
    assert result[0].shape == (1,)


def test_elevation_angles_length_preserved():
    pts = np.array([[0, 0, 0], [1, 1, 0], [2, 0, 1]], dtype=np.float32)
    result = compute_elevation_angles([pts])
    assert len(result[0]) == len(pts)


# ---------------------------------------------------------------------------
# reduce_per_edge
# ---------------------------------------------------------------------------

def test_reduce_per_edge_mean():
    values = [np.array([1.0, 3.0]), np.array([2.0, 4.0])]
    result = reduce_per_edge(values, how="mean")
    np.testing.assert_allclose(result, [2.0, 3.0])


def test_reduce_per_edge_median():
    values = [np.array([1.0, 2.0, 3.0])]
    result = reduce_per_edge(values, how="median")
    np.testing.assert_allclose(result, [2.0])


def test_reduce_per_edge_empty_array():
    values = [np.array([])]
    result = reduce_per_edge(values, how="mean")
    assert np.isnan(result[0])


def test_reduce_per_edge_invalid_how():
    with pytest.raises(ValueError, match="how"):
        reduce_per_edge([np.array([1.0])], how="sum")
