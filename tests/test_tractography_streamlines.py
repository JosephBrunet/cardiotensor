import numpy as np

from cardiotensor.tractography.generate_streamlines import (
    _fa_to_unit_interval,
    generate_streamlines_from_vector_field,
    trace_streamline,
    trilinear_interpolate_vector,
)


def test_trilinear_interpolate_vector_aligns_reversed_corners_to_reference():
    vector_field = np.zeros((3, 1, 1, 2), dtype=np.float32)
    vector_field[:, 0, 0, 0] = [1.0, 0.0, 0.0]
    vector_field[:, 0, 0, 1] = [-1.0, 0.0, 0.0]

    raw = trilinear_interpolate_vector(vector_field, (0.0, 0.0, 0.5))
    aligned = trilinear_interpolate_vector(
        vector_field,
        (0.0, 0.0, 0.5),
        reference_vector=np.array([1.0, 0.0, 0.0]),
    )

    assert np.allclose(raw, [0.0, 0.0, 0.0])
    assert np.allclose(aligned, [1.0, 0.0, 0.0])


def test_trace_streamline_continues_across_reversed_axis_region():
    vector_field = np.zeros((3, 3, 3, 8), dtype=np.float32)
    vector_field[0, :, :, :4] = 1.0
    vector_field[0, :, :, 4:] = -1.0

    streamline = trace_streamline(
        start_pt=(1.0, 1.0, 1.0),
        vector_field=vector_field,
        step_length=0.5,
        max_steps=10,
        angle_threshold=30.0,
    )

    assert len(streamline) > 6
    assert streamline[-1][2] > 4.0


def test_fa_to_unit_interval_decodes_uint8_and_preserves_float_data():
    encoded = np.array([0, 128, 255], dtype=np.uint8)
    decoded = _fa_to_unit_interval(encoded)
    np.testing.assert_allclose(decoded, [0.0, 128.0 / 255.0, 1.0])
    assert decoded.dtype == np.float32

    scientific = np.array([0.0, 0.25, 1.0], dtype=np.float32)
    np.testing.assert_array_equal(_fa_to_unit_interval(scientific), scientific)


def test_streamline_minimum_is_at_least_two_points():
    vector_field = np.zeros((3, 1, 1, 1), dtype=np.float32)
    streamlines = generate_streamlines_from_vector_field(
        vector_field,
        seed_points=np.array([[0, 0, 0]]),
        min_length_pts=0,
        max_steps=1,
        bidirectional=False,
    )

    assert streamlines == []
