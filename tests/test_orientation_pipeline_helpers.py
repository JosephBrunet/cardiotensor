import numpy as np

from cardiotensor.orientation.orientation_computation_pipeline import (
    _centerline_neighborhood,
    _normalize_vectors_in_place,
)


def test_centerline_neighborhood_uses_global_slice_index():
    center_line = np.column_stack(
        (
            np.arange(20, dtype=np.float64),
            np.zeros(20),
            np.arange(20, dtype=np.float64),
        )
    )

    result = _centerline_neighborhood(center_line, global_slice_idx=10, buffer=2)

    np.testing.assert_array_equal(result, center_line[8:13])


def test_normalize_vectors_keeps_masked_zero_vectors_finite():
    vectors = np.zeros((3, 1, 1, 2), dtype=np.float32)
    vectors[:, 0, 0, 1] = [3.0, 4.0, 0.0]

    returned = _normalize_vectors_in_place(vectors)

    assert returned is vectors
    np.testing.assert_array_equal(vectors[:, 0, 0, 0], [0.0, 0.0, 0.0])
    np.testing.assert_allclose(vectors[:, 0, 0, 1], [0.6, 0.8, 0.0])
    assert np.all(np.isfinite(vectors))


def test_normalize_vectors_allocates_norms_one_slice_at_a_time(monkeypatch):
    vectors = np.ones((3, 3, 2, 4), dtype=np.float32)
    original_norm = np.linalg.norm
    shapes = []

    def recording_norm(array, *args, **kwargs):
        shapes.append(array.shape)
        return original_norm(array, *args, **kwargs)

    monkeypatch.setattr(np.linalg, "norm", recording_norm)
    _normalize_vectors_in_place(vectors)

    assert shapes == [(3, 2, 4)] * 3
