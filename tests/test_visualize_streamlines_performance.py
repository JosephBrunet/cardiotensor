import json
from types import SimpleNamespace

import numpy as np
import pytest

import cardiotensor.scripts.visualize_streamlines as visualize_cli
import cardiotensor.visualization.streamlines as streamline_visualization
from cardiotensor.scripts.visualize_streamlines import (
    _apply_session_settings,
    _discover_trk_color_fields,
)
from cardiotensor.visualization import (
    fury_plotting_streamlines,
    pyvista_plotting_streamlines,
)


def test_color_field_discovery_uses_lazy_trk_loading(monkeypatch, tmp_path):
    load_calls = []
    tractogram = SimpleNamespace(data_per_point={"HA": [], "IA": []})

    def fake_load(path, **kwargs):
        load_calls.append((path, kwargs))
        return SimpleNamespace(tractogram=tractogram)

    monkeypatch.setattr("nibabel.streamlines.load", fake_load)

    fields = _discover_trk_color_fields(tmp_path / "streamlines.trk")

    assert fields == ["HA", "IA"]
    assert load_calls == [(str(tmp_path / "streamlines.trk"), {"lazy_load": True})]


def test_visualizer_accepts_custom_per_point_scalar(monkeypatch, tmp_path):
    trk_path = tmp_path / "centroids.trk"
    trk_path.write_bytes(b"trk")
    streamlines = [
        np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32),
        np.array([[0, 1, 0], [1, 1, 0]], dtype=np.float32),
    ]
    cluster_ids = [
        np.array([0, 0], dtype=np.float32),
        np.array([1, 1], dtype=np.float32),
    ]
    captured = {}

    monkeypatch.setattr(
        streamline_visualization,
        "load_trk_streamlines",
        lambda path, include_per_streamline=False: (
            streamlines,
            {"cluster_id": cluster_ids},
            {},
        ),
    )
    monkeypatch.setattr(
        fury_plotting_streamlines,
        "show_streamlines",
        lambda **kwargs: captured.update(kwargs),
    )

    streamline_visualization.visualize_streamlines(
        trk_path, color_by="cluster_id", interactive=False
    )

    np.testing.assert_array_equal(captured["color_values"][0], [0, 0])
    np.testing.assert_array_equal(captured["color_values"][1], [1, 1])
    assert captured["color_range"] == (0.0, 1.0)
    assert captured["color_label"] == "CLUSTER_ID"


def test_top_clusters_keeps_every_member_of_the_largest_clusters(monkeypatch, tmp_path):
    trk_path = tmp_path / "members.trk"
    trk_path.write_bytes(b"trk")
    streamlines = [
        np.array([[0, index, 0], [1, index, 0]], dtype=np.float32) for index in range(6)
    ]
    cluster_ids = [0, 0, 1, 1, 1, 2]
    cluster_values = [
        np.full(len(streamline), cluster_id, dtype=np.float32)
        for streamline, cluster_id in zip(streamlines, cluster_ids)
    ]
    captured = {}

    monkeypatch.setattr(
        streamline_visualization,
        "load_trk_streamlines",
        lambda path, include_per_streamline=False: (
            streamlines,
            {"cluster_id": cluster_values},
            {"cluster_size": np.array([[2], [2], [3], [3], [3], [1]])},
        ),
    )
    monkeypatch.setattr(
        fury_plotting_streamlines,
        "show_streamlines",
        lambda **kwargs: captured.update(kwargs),
    )

    streamline_visualization.visualize_streamlines(
        trk_path,
        color_by="cluster_id",
        top_clusters=2,
        interactive=False,
    )

    assert len(captured["streamlines_xyz"]) == 5
    assert {int(values[0]) for values in captured["color_values"]} == {0, 1}


def test_top_clusters_rejects_normal_tractogram(monkeypatch, tmp_path):
    trk_path = tmp_path / "normal.trk"
    trk_path.write_bytes(b"trk")
    streamlines = [np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)]
    monkeypatch.setattr(
        streamline_visualization,
        "load_trk_streamlines",
        lambda path, include_per_streamline=False: (
            streamlines,
            {"HA": [np.array([0, 0], dtype=np.float32)]},
            {},
        ),
    )

    with pytest.raises(ValueError, match="requires a clustered TRK"):
        streamline_visualization.visualize_streamlines(
            trk_path,
            color_by="HA",
            top_clusters=1,
            interactive=False,
        )


def test_cli_accepts_hyphenated_custom_color_field(monkeypatch, tmp_path):
    trk_path = tmp_path / "centroids.trk"
    trk_path.write_bytes(b"trk")
    captured = {}

    monkeypatch.setattr(
        visualize_cli, "_discover_trk_color_fields", lambda path: ["cluster_id"]
    )
    monkeypatch.setattr(
        visualize_cli, "visualize_streamlines", lambda **kwargs: captured.update(kwargs)
    )
    monkeypatch.setattr(
        visualize_cli.sys,
        "argv",
        [
            "cardio-visualize-streamlines",
            str(trk_path),
            "--color-by",
            "cluster-id",
            "--top-clusters",
            "3",
        ],
    )

    visualize_cli.script()

    assert captured["color_by"] == "cluster_id"
    assert captured["top_clusters"] == 3


def test_cli_uses_default_session_path_when_none_is_given(monkeypatch, tmp_path):
    trk_path = tmp_path / "streamlines.trk"
    trk_path.write_bytes(b"trk")
    captured = {}

    monkeypatch.setattr(visualize_cli, "_discover_trk_color_fields", lambda path: [])
    monkeypatch.setattr(
        visualize_cli, "visualize_streamlines", lambda **kwargs: captured.update(kwargs)
    )
    monkeypatch.setattr(
        visualize_cli.sys, "argv", ["cardio-visualize-streamlines", str(trk_path)]
    )

    visualize_cli.script()

    assert captured["session_path"] == tmp_path / "streamlines_session.json"
    assert captured["restore_session"] is False


def test_cli_asks_before_restoring_default_session(monkeypatch, tmp_path):
    trk_path = tmp_path / "streamlines.trk"
    trk_path.write_bytes(b"trk")
    session_path = tmp_path / "streamlines_session.json"
    session_path.write_text(
        json.dumps(
            {
                "streamlines_file": str(trk_path),
                "streamlines_size": 3,
                "settings": {"line_width": 7.0},
            }
        )
    )
    monkeypatch.setattr(visualize_cli, "_discover_trk_color_fields", lambda path: [])
    monkeypatch.setattr(
        visualize_cli.sys, "stdin", SimpleNamespace(isatty=lambda: True)
    )
    monkeypatch.setattr(
        visualize_cli.sys, "argv", ["cardio-visualize-streamlines", str(trk_path)]
    )

    for answer, expected_restore, expected_width in (
        ("y", True, 7.0),
        ("n", False, 4.0),
    ):
        captured = {}
        monkeypatch.setattr("builtins.input", lambda prompt, value=answer: value)
        monkeypatch.setattr(
            visualize_cli,
            "visualize_streamlines",
            lambda **kwargs: captured.update(kwargs),
        )

        visualize_cli.script()

        assert captured["session_path"] == session_path
        assert captured["restore_session"] is expected_restore
        assert captured["line_width"] == expected_width


def test_pyvista_connectivity_is_built_correctly(monkeypatch):
    class FakePolyData:
        def __init__(self, points):
            self.points = points
            self.lines = None
            self.point_data = {}
            self.active_scalars = None

        def set_active_scalars(self, name):
            self.active_scalars = name

    fake_pyvista = SimpleNamespace(PolyData=FakePolyData)
    monkeypatch.setattr(
        pyvista_plotting_streamlines, "_load_pyvista", lambda: fake_pyvista
    )

    streamlines = [
        np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32),
        np.array([[2, 0, 0], [2, 1, 0], [2, 2, 0]], dtype=np.float32),
    ]
    colors = [
        np.array([10, 20], dtype=np.float32),
        np.array([30, 40, 50], dtype=np.float32),
    ]

    poly = pyvista_plotting_streamlines._streamlines_to_polydata(
        streamlines, colors, "HA"
    )

    np.testing.assert_array_equal(poly.points, np.concatenate(streamlines))
    np.testing.assert_array_equal(poly.lines, [2, 0, 1, 3, 2, 3, 4])
    np.testing.assert_array_equal(poly.point_data["HA"], [10, 20, 30, 40, 50])
    assert poly.active_scalars == "HA"


def test_pyvista_bounds_ignore_empty_streamlines():
    streamlines = [
        np.empty((0, 3), dtype=np.float32),
        np.array([[1, 5, -2], [4, 3, 8]], dtype=np.float32),
        np.array([[-1, 7, 0], [2, 6, 3]], dtype=np.float32),
    ]

    mins, maxs = pyvista_plotting_streamlines._compute_streamline_bounds(streamlines)

    np.testing.assert_array_equal(mins, [-1, 3, -2])
    np.testing.assert_array_equal(maxs, [4, 7, 8])


def test_saved_settings_are_defaults_but_cli_options_win():
    args = SimpleNamespace(
        line_width=2.0,
        subsample=1,
        crop_x=None,
        crop_y=None,
        crop_z=None,
        width=800,
        height=800,
        hide_axes=False,
        show_bounds=False,
        hide_bounds=False,
        shadows=False,
        no_shadows=False,
    )
    settings = {
        "line_width": 6.0,
        "tube_thickness": 0.25,
        "subsample_factor": 4,
        "crop_bounds": [[1, 2], [3, 4], [5, 6]],
        "window_size": [1600, 1000],
        "show_axes": False,
        "show_bounds": True,
        "shadows": True,
    }

    _apply_session_settings(args, settings, {"--line-width"})

    assert args.line_width == 2.0
    assert args.subsample == 4
    assert (args.crop_x, args.crop_y, args.crop_z) == ([1, 2], [3, 4], [5, 6])
    assert (args.width, args.height) == (1600, 1000)
    assert args.hide_axes is True
    assert args.show_bounds is True
    assert args.shadows is True

    _apply_session_settings(args, settings, set())
    assert args.line_width == 0.25


def test_pyvista_subsampling_is_reproducible():
    streamlines = [
        np.array([[index, 0, 0], [index, 1, 0]], dtype=np.float32)
        for index in range(20)
    ]
    colors = [np.array([0, 1], dtype=np.float32) for _ in streamlines]

    first, _, _ = pyvista_plotting_streamlines._prepare_streamlines(
        streamlines,
        colors,
        downsample_factor=1,
        max_streamlines=None,
        filter_min_len=None,
        subsample_factor=4,
        crop_bounds=None,
        random_seed=1234,
    )
    second, _, _ = pyvista_plotting_streamlines._prepare_streamlines(
        streamlines,
        colors,
        downsample_factor=1,
        max_streamlines=None,
        filter_min_len=None,
        subsample_factor=4,
        crop_bounds=None,
        random_seed=1234,
    )

    assert [line[0, 0] for line in first] == [line[0, 0] for line in second]


def make_fury_viewer(tmp_path=None, **kwargs):
    points = np.array([[10, 0, 0], [11, 1, 0], [12, 0, 1]], dtype=np.float32)
    return fury_plotting_streamlines.StreamlineViewer(
        [points],
        [np.array([0.0, 1.0, 2.0])],
        "line",
        2,
        (320, 240),
        fury_plotting_streamlines.ScalarColors("turbo", (0, 2)),
        session_path=tmp_path / "session.json" if tmp_path else None,
        **kwargs,
    )


def test_fury_orbit_frames_use_video_output_scratch(monkeypatch, tmp_path):
    viewer = make_fury_viewer()
    initial = viewer.camera.local.position.copy()
    monkeypatch.setattr(
        viewer, "_save_screenshot", lambda path: path.write_bytes(b"png")
    )
    video_path = tmp_path / "analysis" / "orbit.mp4"

    def writer(paths, output, fps):
        assert all(
            path.parent.parent == video_path.parent / ".cardiotensor_scratch"
            for path in paths
        )
        assert len(paths) == 2
        assert fps == 30
        output.write_bytes(b"video")

    monkeypatch.setattr(
        fury_plotting_streamlines, "_write_frame_sequence_to_video", writer
    )
    viewer._record_orbit(video_path, 2, 30)
    assert video_path.exists()
    np.testing.assert_allclose(viewer.camera.local.position, initial)
    assert not list(
        (video_path.parent / ".cardiotensor_scratch").glob("cardiotensor_fury_orbit_*")
    )


def test_fury_window_close_saves_once_and_terminates():
    viewer = make_fury_viewer()
    calls = []
    viewer.showm = SimpleNamespace(close=lambda: calls.append("closed"))
    viewer._autosave_session = lambda: calls.append("saved")
    viewer._close_window()
    viewer._close_window()
    assert calls == ["saved", "closed"]


def test_fury_screenshot_uses_current_viewer_window(monkeypatch, tmp_path):
    from cardiotensor.visualization.fury_support import capture

    calls = []
    renderer = SimpleNamespace(pixel_ratio=1)
    showm = SimpleNamespace(
        renderer=renderer,
        render=lambda: calls.append("render"),
        window=SimpleNamespace(draw=lambda: calls.append("draw")),
        snapshot=lambda fname: np.zeros(
            (int(180 * renderer.pixel_ratio), int(320 * renderer.pixel_ratio), 4),
            dtype=np.uint8,
        ),
    )
    assert capture(showm, tmp_path / "image.png", scale=2) == (640, 360)
    assert renderer.pixel_ratio == 1
    assert calls == ["render", "draw"]


def test_fury_spline_subdivision_and_crop_box_rotation():
    import pylinalg as la

    from cardiotensor.visualization.fury_support import materials

    viewer = make_fury_viewer(spline_subdiv=2)
    assert viewer.render_lines[0].shape == (3, 3)
    np.testing.assert_allclose(viewer.actor0.local.position, [0, 0, 0])
    viewer.box_object.local.rotation = la.quat_from_axis_angle((0, 0, 1), np.pi / 4)
    viewer.box_clipping_active = True
    viewer._apply_clipping_planes()
    assert len(next(materials(viewer.actor0)).clipping_planes) == 6
    equations = np.array(viewer._box_equations())
    assert np.all(equations[:, :3] @ viewer.center >= equations[:, 3] - 1e-6)
    outside = viewer.center + [100, 0, 0]
    assert np.any(equations[:, :3] @ outside < equations[:, 3])
    viewer.clipping_active = True
    viewer._apply_clipping_planes()
    assert len(next(materials(viewer.actor0)).clipping_planes) == 7


def test_ctrl_s_saves_complete_fury_session(tmp_path):
    viewer = make_fury_viewer(tmp_path)
    viewer.clipping_active = True
    viewer.box_clipping_active = True
    viewer.box_object.local.scale = (0.5, 0.7, 0.9)
    viewer._on_keypress(SimpleNamespace(key="s", modifiers=["Control"]))
    saved = json.loads(viewer.session_path.read_text())
    assert saved["format"] == "cardiotensor-fury-session"
    assert saved["version"] == 2
    assert saved["view"]["clipping_plane"]["enabled"]
    assert saved["view"]["crop_box"]["enabled"]
    np.testing.assert_allclose(
        saved["view"]["crop_box"]["matrix"], viewer.box_object.local.matrix
    )
    assert saved["settings"]["line_width"] == 2


def test_fury_session_restores_camera_plane_and_controls(tmp_path):
    viewer = make_fury_viewer(tmp_path)
    viewer.camera.local.position = [10, 20, 30]
    viewer.plane_object.local.position = [4, 5, 6]
    viewer.clipping_active = viewer.box_clipping_active = True
    viewer.box_object.local.scale = [0.5, 0.6, 0.7]
    viewer.current_bg = (1, 1, 1)
    viewer.scale_bar_on = False
    viewer.controls_visible = True
    viewer._save_session()
    restored = make_fury_viewer(tmp_path)
    assert restored._restore_session()
    np.testing.assert_allclose(restored.camera.local.position, [10, 20, 30])
    np.testing.assert_allclose(restored.plane_object.local.position, [4, 5, 6])
    np.testing.assert_allclose(restored.box_object.local.scale, [0.5, 0.6, 0.7])
    assert restored.clipping_active and restored.box_clipping_active
    assert restored.current_bg == (1, 1, 1)
    assert not restored.scale_bar_on
    assert restored.controls_visible


def test_fury_legacy_session_migrates(tmp_path):
    viewer = make_fury_viewer(tmp_path)
    viewer.session_path.write_text(
        json.dumps(
            {
                "version": 1,
                "view": {
                    "camera": {
                        "position": [10, 20, 30],
                        "focal_point": [1, 2, 3],
                        "view_up": [0, 0, 1],
                        "parallel_projection": True,
                        "parallel_scale": 12,
                    },
                    "clipping_plane": {
                        "origin": [4, 5, 6],
                        "normal": [0, 1, 0],
                        "enabled": True,
                    },
                },
            }
        )
    )
    assert viewer._restore_session()
    np.testing.assert_allclose(viewer.camera.local.position, [10, 20, 30])
    assert viewer.camera.fov == 0
    assert viewer.camera.height == 24
    np.testing.assert_allclose(viewer._plane_equation(), [0, 1, 0, 5], atol=1e-6)


def test_fury_smoothing_preserves_scalar_endpoints():
    viewer = make_fury_viewer(spline_subdiv=20)
    assert len(viewer.render_lines[0]) == 21
    np.testing.assert_allclose(viewer.render_colors[0][[0, -1]], viewer.lut.map([0, 2]))
