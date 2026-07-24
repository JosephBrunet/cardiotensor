import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import cardiotensor.scripts.visualize_streamlines as visualize_cli
import cardiotensor.visualization.streamlines as streamline_visualization
from cardiotensor.scripts.visualize_streamlines import (
    _apply_session_settings,
    _discover_trk_color_fields,
)
from cardiotensor.visualization import fury_plotting_streamlines
from cardiotensor.visualization import pyvista_plotting_streamlines


def test_fury_orbit_frames_use_video_output_scratch(monkeypatch, tmp_path):
    viewer = object.__new__(fury_plotting_streamlines.StreamlineViewer)
    viewer.actor_fast = SimpleNamespace(SetVisibility=lambda visible: None)
    viewer.actor0 = SimpleNamespace(SetVisibility=lambda visible: None)
    viewer.scene = SimpleNamespace(azimuth=lambda degrees: None)
    viewer.window_size = (100, 100)
    viewer._set_default_camera = lambda: None

    monkeypatch.setattr(
        fury_plotting_streamlines.fury.window,
        "record",
        lambda scene, out_path, size, reset_camera: Path(out_path).write_bytes(b"png"),
    )

    video_path = tmp_path / "analysis" / "orbit.mp4"

    def fake_video_writer(frame_paths, out_path, fps):
        scratch_root = video_path.parent / ".cardiotensor_scratch"
        assert all(path.parent.parent == scratch_root for path in frame_paths)
        Path(out_path).write_bytes(b"video")

    monkeypatch.setattr(
        fury_plotting_streamlines,
        "_write_frame_sequence_to_video",
        fake_video_writer,
    )

    viewer._record_orbit(video_path, video_frames=2, video_fps=30)

    assert video_path.exists()
    assert not list(
        (video_path.parent / ".cardiotensor_scratch").glob(
            "cardiotensor_fury_orbit_*"
        )
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


def test_top_clusters_keeps_every_member_of_the_largest_clusters(
    monkeypatch, tmp_path
):
    trk_path = tmp_path / "members.trk"
    trk_path.write_bytes(b"trk")
    streamlines = [
        np.array([[0, index, 0], [1, index, 0]], dtype=np.float32)
        for index in range(6)
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
    assert {
        int(values[0]) for values in captured["color_values"]
    } == {0, 1}


def test_top_clusters_rejects_normal_tractogram(monkeypatch, tmp_path):
    trk_path = tmp_path / "normal.trk"
    trk_path.write_bytes(b"trk")
    streamlines = [
        np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)
    ]
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


def test_fury_window_close_saves_once_and_terminates():
    calls = {"saved": 0, "exited": 0, "widgets_disabled": 0}

    class FakeWidget:
        def EnabledOff(self):
            calls["widgets_disabled"] += 1

    viewer = fury_plotting_streamlines.StreamlineViewer.__new__(
        fury_plotting_streamlines.StreamlineViewer
    )
    viewer._closing = False
    viewer.plane_widget = FakeWidget()
    viewer.box_widget = FakeWidget()
    viewer.showm = SimpleNamespace(
        exit=lambda: calls.__setitem__("exited", calls["exited"] + 1)
    )
    viewer._autosave_session = lambda: calls.__setitem__("saved", calls["saved"] + 1)

    viewer._close_window()
    viewer._close_window()

    assert calls == {"saved": 1, "exited": 1, "widgets_disabled": 2}


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


def test_fury_screenshot_uses_current_viewer_window(monkeypatch, tmp_path):
    calls = {}

    class FakeRenderWindow:
        def Render(self):
            calls["rendered"] = True

    class FakeOutput:
        def __init__(self, scale):
            self.scale = scale

        def GetDimensions(self):
            return (3200 * self.scale, 1800 * self.scale, 1)

    class FakeCapture:
        scale = 1

        def SetInput(self, render_window):
            calls["capture_input"] = render_window

        def SetInputBufferTypeToRGB(self):
            calls["rgb"] = True

        def Update(self):
            calls["updated"] = True

        def SetScale(self, scale):
            self.scale = scale
            calls["scale"] = scale

        def GetOutputPort(self):
            return self

        def GetOutput(self):
            return FakeOutput(self.scale)

    class FakeWriter:
        def SetFileName(self, filename):
            calls["filename"] = filename

        def SetInputConnection(self, connection):
            calls["writer_input"] = connection

        def Write(self):
            calls["written"] = True

    capture = FakeCapture()
    monkeypatch.setattr(
        fury_plotting_streamlines.vtk,
        "vtkWindowToImageFilter",
        lambda: capture,
    )
    monkeypatch.setattr(fury_plotting_streamlines.vtk, "vtkPNGWriter", FakeWriter)

    render_window = FakeRenderWindow()
    viewer = fury_plotting_streamlines.StreamlineViewer.__new__(
        fury_plotting_streamlines.StreamlineViewer
    )
    viewer.showm = SimpleNamespace(window=render_window)
    out_path = tmp_path / "view.png"

    size = viewer._save_screenshot(out_path)

    assert size == (3200, 1800)
    assert calls["capture_input"] is render_window
    assert calls["filename"] == str(out_path)
    assert calls["rendered"] and calls["updated"] and calls["written"]

    size = viewer._save_screenshot(out_path, scale=2)
    assert size == (6400, 3600)
    assert calls["scale"] == 2


def test_fury_spline_subdivision_and_crop_box_rotation():
    streamlines = [
        np.array([[index, index % 2, 0] for index in range(10)], dtype=np.float32)
    ]
    colors = [np.arange(10, dtype=np.float32)]
    lut = fury_plotting_streamlines.actor.colormap_lookup_table(scale_range=(0, 9))

    viewer = fury_plotting_streamlines.StreamlineViewer(
        streamlines,
        colors,
        "line",
        2,
        (800, 600),
        lut,
        spline_subdiv=2,
    )
    assert viewer.actor_spline_subdiv == 2
    viewer.actor0.GetMapper().Update()
    assert viewer.actor0.GetMapper().GetInput().GetNumberOfPoints() == 3

    viewer.showm = fury_plotting_streamlines.window.ShowManager(
        scene=viewer.scene, size=(800, 600), reset_camera=False
    )
    viewer._setup_box_widget()
    assert viewer.box_widget.GetRotationEnabled() == 1

    viewer.box_clipping_active = True
    viewer._apply_clipping_planes()
    assert viewer.actor0.GetMapper().GetNumberOfClippingPlanes() == 6

    viewer.clipping_active = True
    viewer._apply_clipping_planes()
    assert viewer.actor0.GetMapper().GetNumberOfClippingPlanes() == 7


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


def test_ctrl_s_saves_complete_fury_session(tmp_path):
    class FakeCamera:
        def GetPosition(self):
            return (10, 20, 30)

        def GetFocalPoint(self):
            return (1, 2, 3)

        def GetViewUp(self):
            return (0, 0, 1)

        def GetClippingRange(self):
            return (0.1, 5000)

        def GetViewAngle(self):
            return 30

        def GetParallelProjection(self):
            return False

        def GetParallelScale(self):
            return 12

    class FakePlane:
        def GetOrigin(self, output):
            output[:] = [4, 5, 6]

        def GetNormal(self, output):
            output[:] = [0, 1, 0]

    trk_path = tmp_path / "streamlines.trk"
    trk_path.write_bytes(b"trk")
    session_path = tmp_path / "my_view.json"
    viewer = fury_plotting_streamlines.StreamlineViewer.__new__(
        fury_plotting_streamlines.StreamlineViewer
    )
    viewer.session_path = session_path
    viewer.session_settings = {"subsample_factor": 4, "random_seed": 1234}
    viewer.streamlines_file = trk_path
    viewer.scene = SimpleNamespace(GetActiveCamera=lambda: FakeCamera())
    viewer.plane_rep = FakePlane()
    viewer.plane_widget = SimpleNamespace(GetEnabled=lambda: True)
    viewer.showm = SimpleNamespace(window=SimpleNamespace(GetSize=lambda: (1800, 1200)))
    viewer.window_size = (800, 800)
    viewer.linewidth = 2.5
    viewer.current_bg = (0.0, 0.0, 0.0)
    viewer.clipping_active = True
    viewer.scale_bar_on = False
    viewer.material = {
        "ambient": 0.4,
        "diffuse": 0.7,
        "specular": 0.2,
        "opacity": 0.9,
    }
    viewer.quality = "publication"
    viewer.box_rep = None
    viewer.box_widget = None
    viewer.box_clipping_active = False
    viewer.controls_visible = False

    class ControlS:
        def GetKeySym(self):
            return "s"

        def GetControlKey(self):
            return 1

    viewer._on_keypress(ControlS(), None)
    saved = json.loads(session_path.read_text())

    assert saved["format"] == "cardiotensor-fury-session"
    assert saved["streamlines_size"] == 3
    assert saved["settings"]["random_seed"] == 1234
    assert saved["settings"]["line_width"] == 2.5
    assert saved["settings"]["tube_thickness"] == 2.5
    assert saved["view"]["camera"]["position"] == [10, 20, 30]
    assert saved["view"]["clipping_plane"] == {
        "origin": [4, 5, 6],
        "normal": [0, 1, 0],
        "enabled": True,
        "gizmo_visible": True,
    }
    assert saved["view"]["window_size"] == [1800, 1200]


def test_fury_session_restores_camera_plane_and_controls(tmp_path):
    calls = {}
    session_path = tmp_path / "my_view.json"
    session_path.write_text(
        json.dumps(
            {
                "view": {
                    "camera": {
                        "position": [10, 20, 30],
                        "focal_point": [1, 2, 3],
                        "view_up": [0, 0, 1],
                        "clipping_range": [0.1, 5000],
                        "view_angle": 25,
                        "parallel_projection": True,
                        "parallel_scale": 12,
                    },
                    "clipping_plane": {
                        "origin": [4, 5, 6],
                        "normal": [0, 1, 0],
                        "enabled": True,
                        "gizmo_visible": True,
                    },
                    "background_color": [1, 1, 1],
                    "scale_bar_visible": False,
                }
            }
        )
    )

    class FakeCamera:
        def __getattr__(self, name):
            if name.startswith("Set"):
                return lambda *values: calls.__setitem__(name, values)
            raise AttributeError(name)

    class FakePlaneRepresentation:
        origin = [0, 0, 0]
        normal = [1, 0, 0]

        def SetOrigin(self, *values):
            self.origin = list(values)

        def SetNormal(self, *values):
            self.normal = list(values)

        def GetOrigin(self, output):
            output[:] = self.origin

        def GetNormal(self, output):
            output[:] = self.normal

        def UpdatePlacement(self):
            calls["plane_updated"] = True

    class FakePlane:
        def SetOrigin(self, *values):
            calls["plane_origin"] = values

        def SetNormal(self, *values):
            calls["plane_normal"] = values

    class FakeMapper:
        def RemoveAllClippingPlanes(self):
            calls["planes_removed"] = calls.get("planes_removed", 0) + 1

        def AddClippingPlane(self, plane):
            calls["planes_added"] = calls.get("planes_added", 0) + 1

    class FakeWidget:
        def EnabledOn(self):
            calls["gizmo"] = True

        def EnabledOff(self):
            calls["gizmo"] = False

    mapper_a = FakeMapper()
    mapper_b = FakeMapper()
    camera = FakeCamera()
    scene = SimpleNamespace(
        GetActiveCamera=lambda: camera,
        SetBackground=lambda *values: calls.__setitem__("background", values),
        add=lambda actor: calls.__setitem__("scale_bar_added", actor),
        rm=lambda actor: calls.__setitem__("scale_bar_removed", actor),
    )
    viewer = fury_plotting_streamlines.StreamlineViewer.__new__(
        fury_plotting_streamlines.StreamlineViewer
    )
    viewer.session_path = session_path
    viewer.scene = scene
    viewer.plane_rep = FakePlaneRepresentation()
    viewer.plane_fn = FakePlane()
    viewer.actor0 = SimpleNamespace(GetMapper=lambda: mapper_a)
    viewer.actor_fast = SimpleNamespace(GetMapper=lambda: mapper_b)
    viewer.plane_widget = FakeWidget()
    viewer.box_rep = None
    viewer.box_widget = None
    viewer.box_clipping_active = False
    viewer.clipping_active = False
    viewer.current_bg = (0, 0, 0)
    viewer.scale_bar = object()
    viewer.scale_bar_on = True
    viewer.showm = SimpleNamespace(render=lambda: calls.__setitem__("rendered", True))

    assert viewer._restore_session() is True
    assert calls["SetPosition"] == (10, 20, 30)
    assert calls["SetParallelProjection"] == (True,)
    assert calls["plane_origin"] == ([4, 5, 6],)
    assert calls["plane_normal"] == ([0, 1, 0],)
    assert calls["planes_added"] == 2
    assert calls["gizmo"] is True
    assert calls["background"] == (1.0, 1.0, 1.0)
    assert calls["scale_bar_removed"] is viewer.scale_bar
    assert calls["rendered"] is True
