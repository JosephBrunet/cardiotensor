import sys

import numpy as np

import cardiotensor.scripts.convert_trk_to_am as converter


def test_converter_skips_one_point_streamlines(tmp_path, monkeypatch):
    input_path = tmp_path / "streamlines.trk"
    input_path.touch()
    streamlines = [
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
    ]
    attrs = {
        "HA": [np.array([5.0]), np.array([10.0, 20.0])],
        "IA": [np.array([6.0]), np.array([30.0, 40.0])],
    }
    captured = {}

    monkeypatch.setattr(
        converter,
        "load_trk_streamlines",
        lambda path: (streamlines, attrs),
    )

    def capture_write(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(converter, "write_spatialgraph_am", capture_write)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cardio-trk_2_am",
            str(input_path),
            "--edge-scalar-sources",
            "ha,ia",
        ],
    )

    converter.script()

    assert len(captured["streamlines_xyz"]) == 1
    assert len(captured["streamlines_xyz"][0]) == 2
    np.testing.assert_allclose(captured["edge_scalar"]["HA"], [15.0])
    np.testing.assert_allclose(captured["edge_scalar"]["IA"], [35.0])
