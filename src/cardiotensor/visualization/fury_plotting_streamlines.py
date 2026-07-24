#!/usr/bin/env python3
from __future__ import annotations

import datetime
import inspect
import json
import random
import tempfile
from pathlib import Path

import fury
import matplotlib.pyplot as plt
import numpy as np
import vtk
from fury import actor, ui, window

# ---------------------------
# small utilities
# ---------------------------


def parse_background_color(color) -> tuple[float, float, float]:
    NAMED = {
        "white": (1.0, 1.0, 1.0),
        "black": (0.0, 0.0, 0.0),
        "gray": (0.5, 0.5, 0.5),
        "lightgray": (0.9, 0.9, 0.9),
        "red": (1.0, 0.0, 0.0),
        "green": (0.0, 1.0, 0.0),
        "blue": (0.0, 0.0, 1.0),
    }
    if isinstance(color, str):
        key = color.lower()
        if key not in NAMED:
            raise ValueError(
                f"Unknown background color '{color}'. Available: {list(NAMED.keys())}"
            )
        return NAMED[key]
    elif isinstance(color, tuple | list) and len(color) == 3:
        return tuple(float(c) for c in color)
    else:
        raise TypeError("background_color must be str or 3-tuple")


def downsample_streamline(streamline: np.ndarray, factor: int = 2) -> np.ndarray:
    return streamline if len(streamline) < 3 or factor <= 1 else streamline[::factor]


def _supported_actor_kwargs(func, **kwargs):
    try:
        params = inspect.signature(func).parameters
    except (TypeError, ValueError):
        return kwargs

    if any(param.kind == inspect.Parameter.VAR_KEYWORD for param in params.values()):
        return kwargs
    return {name: value for name, value in kwargs.items() if name in params}


def _write_frame_sequence_to_video(
    frame_paths: list[Path], output_path: Path, fps: int
) -> None:
    if not frame_paths:
        raise ValueError("No FURY frames were recorded for the video.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = output_path.suffix.lower()
    if suffix == ".gif":
        try:
            import imageio.v2 as imageio
        except ImportError as err:
            raise RuntimeError(
                "GIF export requires imageio: pip install imageio"
            ) from err

        with imageio.get_writer(
            str(output_path), mode="I", duration=1.0 / fps
        ) as writer:
            for frame_path in frame_paths:
                writer.append_data(imageio.imread(frame_path))
        return

    try:
        import cv2
    except ImportError as err:
        raise RuntimeError(
            "MP4 export requires OpenCV: pip install opencv-python"
        ) from err

    first = cv2.imread(str(frame_paths[0]), cv2.IMREAD_COLOR)
    if first is None:
        raise RuntimeError(f"Could not read recorded frame: {frame_paths[0]}")

    height, width = first.shape[:2]
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    writer = cv2.VideoWriter(
        str(output_path), fourcc, float(fps), (width, height), isColor=True
    )
    if not writer.isOpened():
        raise RuntimeError(f"Could not open video writer for {output_path}")

    try:
        for frame_path in frame_paths:
            frame = cv2.imread(str(frame_path), cv2.IMREAD_COLOR)
            if frame is None:
                raise RuntimeError(f"Could not read recorded frame: {frame_path}")
            if frame.shape[:2] != (height, width):
                frame = cv2.resize(frame, (width, height), interpolation=cv2.INTER_AREA)
            writer.write(frame)
    finally:
        writer.release()


def _compute_streamline_bounds(
    streamlines_xyz: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    if not streamlines_xyz:
        raise ValueError("No streamlines available to compute bounds.")

    mins = np.full(3, np.inf)
    maxs = np.full(3, -np.inf)
    for streamline in streamlines_xyz:
        if len(streamline) == 0:
            continue
        mins = np.minimum(mins, np.min(streamline, axis=0))
        maxs = np.maximum(maxs, np.max(streamline, axis=0))
    if not np.all(np.isfinite(mins)):
        raise ValueError("No non-empty streamlines available to compute bounds.")
    return mins, maxs


def _print_box_shape(label: str, mins: np.ndarray, maxs: np.ndarray) -> None:
    size = maxs - mins
    print(
        f"{label}: \n"
        f"x=[{mins[0]:.2f}, {maxs[0]:.2f}], \n"
        f"y=[{mins[1]:.2f}, {maxs[1]:.2f}], \n"
        f"z=[{mins[2]:.2f}, {maxs[2]:.2f}] \n"
        f"-> shape=({size[0]:.2f}, {size[1]:.2f}, {size[2]:.2f})\n"
    )


def matplotlib_cmap_to_fury_lut(
    cmap, value_range=(-1, 1), n_colors=256
) -> vtk.vtkLookupTable:
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    colors = cmap(np.linspace(0, 1, n_colors))
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(n_colors)
    lut.SetRange(*value_range)
    lut.Build()
    for i in range(n_colors):
        r, g, b, a = colors[i]
        lut.SetTableValue(i, r, g, b, a)
    return lut


def _split_streamline_by_bounds(
    sl: np.ndarray, cl: np.ndarray, x_min, x_max, y_min, y_max, z_min, z_max
):
    within = (
        (sl[:, 0] >= x_min)
        & (sl[:, 0] <= x_max)
        & (sl[:, 1] >= y_min)
        & (sl[:, 1] <= y_max)
        & (sl[:, 2] >= z_min)
        & (sl[:, 2] <= z_max)
    )
    if not np.any(within):
        return [], []

    w = within.astype(np.int8)
    trans = np.diff(np.pad(w, (1, 1), constant_values=0))
    starts = np.where(trans == +1)[0]
    ends = np.where(trans == -1)[0]

    segs, cols = [], []
    for s, e in zip(starts, ends):
        seg = sl[s:e]
        col = cl[s:e]
        if len(seg) > 0:
            segs.append(seg)
            cols.append(col)
    return segs, cols


# ===========================
# Class-based viewer
# ===========================
class StreamlineViewer:
    def __init__(
        self,
        streamlines_xyz,
        color_values,
        mode,
        line_width,
        window_size,
        lut,
        background_color="black",
        spline_subdiv=16,
        tube_sides: int = 9,
        color_range: tuple[float, float] | None = None,
        color_label: str = "Angle",
        flat_values: np.ndarray | None = None,
        session_path: str | Path | None = None,
        restore_session: bool = True,
        session_settings: dict | None = None,
        streamlines_file: str | Path | None = None,
        opacity: float = 1.0,
        quality: str = "interactive",
    ):
        self.streamlines_xyz = streamlines_xyz
        self.color_values = color_values
        self.mode = mode
        self.window_size = window_size
        self.lut = lut
        self.color_label = color_label
        self.session_path = (
            Path(session_path).expanduser().resolve() if session_path else None
        )
        self.restore_session = bool(restore_session)
        self.session_settings = dict(session_settings or {})
        self.streamlines_file = (
            Path(streamlines_file).resolve() if streamlines_file else None
        )
        self.quality = quality.lower().strip()
        if self.quality not in {"interactive", "publication"}:
            raise ValueError("quality must be 'interactive' or 'publication'")

        # Scene
        self.scene = fury.window.Scene()
        self.current_bg = parse_background_color(background_color)
        self.scene.SetBackground(*self.current_bg)
        self._setup_lighting()

        # thickness state
        self.linewidth = max(0.001, float(line_width))  # used for both line and tube
        self.spline_subdiv = max(1, int(spline_subdiv))
        self.actor_spline_subdiv = (
            self.spline_subdiv if self.spline_subdiv > 1 else None
        )
        self.tube_sides = max(3, int(tube_sides))
        if self.quality == "publication":
            self.material = {
                "ambient": 0.30,
                "diffuse": 0.80,
                "specular": 0.25,
                "opacity": float(np.clip(opacity, 0.05, 1.0)),
            }
        else:
            self.material = {
                "ambient": 0.45,
                "diffuse": 0.70,
                "specular": 0.12,
                "opacity": float(np.clip(opacity, 0.05, 1.0)),
            }
        self.controls_panel = None
        self.controls_visible = False
        self.control_sliders = {}

        self.scale_bar = None
        self.scale_bar_on = False
        self._add_scale_bar()

        self.clipping_active = False  # Track clipping state

        # VTK/FURY objects
        self.showm: window.ShowManager | None = None
        self._closing = False

        # clipped branch objects
        self.plane_rep = None
        self.plane_fn = None
        self.plane_widget = None
        self.box_widget = None
        self.box_rep = None
        self.box_planes = vtk.vtkPlanes()
        self.box_clipping_active = False
        self.mapper0 = None
        self.actor0 = None

        # precompute flat scalars for LUT mapping
        if flat_values is None:
            flat_values = np.concatenate(
                [np.asarray(c).ravel() for c in self.color_values]
            )
        self.flat_vals = np.asarray(flat_values, dtype=np.float32).reshape(-1)
        if color_range is None:
            self.vmin = float(np.nanmin(self.flat_vals))
            self.vmax = float(np.nanmax(self.flat_vals))
        else:
            self.vmin, self.vmax = (float(color_range[0]), float(color_range[1]))
        self.lut.SetRange(self.vmin, self.vmax)

        # bounds and center from NumPy
        mins, maxs = _compute_streamline_bounds(self.streamlines_xyz)
        self.center = (mins + maxs) / 2.0
        self.bounds = [mins[0], maxs[0], mins[1], maxs[1], mins[2], maxs[2]]

        self._build_pipeline()

        # self._add_origin_marker()

        self._add_scalar_bar()

    def _build_pipeline(self):
        # create actor, pass scalars and LUT directly
        if self.mode == "tube":
            self.actor0 = actor.streamtube(
                self.streamlines_xyz,
                colors=self.flat_vals,
                linewidth=self.linewidth,
                spline_subdiv=self.actor_spline_subdiv,
                lookup_colormap=self.lut,
                tube_sides=self.tube_sides,
                **_supported_actor_kwargs(actor.streamtube, lod=False),
            )
        else:
            self.actor0 = actor.line(
                self.streamlines_xyz,
                colors=self.flat_vals,  # scalars
                linewidth=self.linewidth,
                spline_subdiv=self.actor_spline_subdiv,
                lookup_colormap=self.lut,
                **_supported_actor_kwargs(actor.line, lod=False),
            )

        self.scene.add(self.actor0)
        self._style_streamline_actor()

        # fast actor for interaction (cheap line rendering)
        self.actor_fast = actor.line(
            self.streamlines_xyz,
            colors=self.flat_vals,  # pass scalars
            linewidth=1.0,  # very light
            lookup_colormap=self.lut,
            fake_tube=True,  # a bit of shading to hint tubes
            **_supported_actor_kwargs(actor.line, lod=False),
        )
        self.actor_fast.SetVisibility(False)
        self.actor_fast.GetProperty().SetOpacity(self.material["opacity"])
        self.scene.add(self.actor_fast)

        # clipping plane setup, start disabled
        self.plane_rep = vtk.vtkImplicitPlaneRepresentation()
        self.plane_rep.SetPlaceFactor(1.25)
        self.plane_rep.PlaceWidget(self.bounds)
        self.plane_rep.SetOrigin(*self.center)
        self.plane_rep.SetNormal(1, 0, 0)

        self.plane_fn = vtk.vtkPlane()
        origin = [0.0, 0.0, 0.0]
        normal = [1.0, 0.0, 0.0]
        self.plane_rep.GetOrigin(origin)
        self.plane_rep.GetNormal(normal)
        self.plane_fn.SetOrigin(origin)
        self.plane_fn.SetNormal(normal)

        self.mapper0 = self.actor0.GetMapper()
        self.mapper0.RemoveAllClippingPlanes()  # start with clipping off

    def _style_streamline_actor(self):
        prop = self.actor0.GetProperty()
        prop.SetInterpolationToPhong()
        prop.SetAmbient(self.material["ambient"])
        prop.SetDiffuse(self.material["diffuse"])
        prop.SetSpecular(self.material["specular"])
        prop.SetSpecularPower(20 if self.quality == "publication" else 10)
        prop.SetOpacity(self.material["opacity"])

    def _setup_lighting(self):
        """Use stable camera-relative key, fill and rim lights."""
        try:
            self.scene.RemoveAllLights()
            headlight = vtk.vtkLight()
            headlight.SetLightTypeToHeadlight()
            headlight.SetIntensity(0.80)

            fill = vtk.vtkLight()
            fill.SetLightTypeToCameraLight()
            fill.SetPosition(-1.0, 1.0, 0.5)
            fill.SetFocalPoint(0.0, 0.0, 0.0)
            fill.SetIntensity(0.35)

            rim = vtk.vtkLight()
            rim.SetLightTypeToCameraLight()
            rim.SetPosition(1.0, -1.0, 0.75)
            rim.SetFocalPoint(0.0, 0.0, 0.0)
            rim.SetIntensity(0.25)

            self.lights = [headlight, fill, rim]
            for light in self.lights:
                self.scene.AddLight(light)
        except Exception:
            self.lights = []

    def _apply_material(self):
        self._style_streamline_actor()
        if self.actor_fast is not None:
            self.actor_fast.GetProperty().SetOpacity(self.material["opacity"])
        self._render_now()

    def _set_material_value(self, name: str, value: float):
        self.material[name] = float(value)
        self._apply_material()

    def _build_controls_panel(self):
        panel = ui.Panel2D(
            size=(280, 255),
            position=(12, 12),
            color=(0.08, 0.09, 0.12),
            opacity=0.88,
            has_border=True,
            border_color=(0.55, 0.58, 0.65),
            border_width=1,
        )
        title = ui.TextBlock2D(
            text="FURY material controls (L)",
            font_size=16,
            bold=True,
            color=(1, 1, 1),
        )
        panel.add_element(title, (12, 225))

        specs = (
            ("ambient", "Ambient", 0.0, 1.0),
            ("diffuse", "Diffuse", 0.0, 1.0),
            ("specular", "Specular", 0.0, 1.0),
            ("opacity", "Opacity", 0.05, 1.0),
        )
        for row, (name, label_text, minimum, maximum) in enumerate(specs):
            y = 185 - row * 52
            label = ui.TextBlock2D(
                text=label_text, font_size=13, color=(0.92, 0.94, 1.0)
            )
            slider = ui.LineSlider2D(
                initial_value=self.material[name],
                min_value=minimum,
                max_value=maximum,
                length=145,
                line_width=3,
                outer_radius=7,
                font_size=12,
                text_template="{value:.2f}",
            )
            slider.on_change = (
                lambda current_slider, material_name=name: self._set_material_value(
                    material_name, current_slider.value
                )
            )
            panel.add_element(label, (12, y))
            panel.add_element(slider, (180, y + 5), anchor="center")
            self.control_sliders[name] = slider

        self.controls_panel = panel
        self.scene.add(panel)
        panel.set_visibility(False)

    def _toggle_controls_panel(self):
        if self.controls_panel is None:
            return
        self.controls_visible = not self.controls_visible
        self.controls_panel.set_visibility(self.controls_visible)
        self._render_now()
        print(f"Material controls {'ON' if self.controls_visible else 'OFF'}")

    def _update_overlay_colors(self):
        luminance = (
            0.2126 * self.current_bg[0]
            + 0.7152 * self.current_bg[1]
            + 0.0722 * self.current_bg[2]
        )
        color = (1.0, 1.0, 1.0) if luminance < 0.5 else (0.0, 0.0, 0.0)
        if getattr(self, "scalar_bar", None) is not None:
            self.scalar_bar.GetTitleTextProperty().SetColor(*color)
            self.scalar_bar.GetLabelTextProperty().SetColor(*color)
            try:
                self.scalar_bar.GetAnnotationTextProperty().SetColor(*color)
            except Exception:
                pass
        if self.scale_bar is not None:
            for axis_name in ("Bottom", "Left", "Right", "Top"):
                try:
                    axis = getattr(self.scale_bar, f"Get{axis_name}Axis")()
                    axis.GetLabelTextProperty().SetColor(*color)
                    axis.GetTitleTextProperty().SetColor(*color)
                except Exception:
                    pass

    def _add_origin_marker(self):
        bounds_size = np.array(
            [
                self.bounds[1] - self.bounds[0],
                self.bounds[3] - self.bounds[2],
                self.bounds[5] - self.bounds[4],
            ],
            dtype=float,
        )
        diagonal = float(np.linalg.norm(bounds_size))
        marker_radius = max(1.0, 0.01 * diagonal)

        self.origin_actor = actor.sphere(
            centers=np.array([[0.0, 0.0, 0.0]], dtype=np.float32),
            colors=np.array([[1.0, 0.0, 0.0]], dtype=np.float32),
            radii=np.array([marker_radius], dtype=np.float32),
        )
        try:
            self.origin_actor.GetProperty().SetAmbient(1.0)
            self.origin_actor.GetProperty().SetDiffuse(0.0)
        except Exception:
            pass
        self.scene.add(self.origin_actor)
        origin_inside = (
            self.bounds[0] <= 0.0 <= self.bounds[1]
            and self.bounds[2] <= 0.0 <= self.bounds[3]
            and self.bounds[4] <= 0.0 <= self.bounds[5]
        )
        print(
            f"Origin marker added at (0, 0, 0) with radius {marker_radius:.2f} "
            f"(inside current box: {origin_inside})"
        )

    def _add_scalar_bar(self):
        self.scalar_bar = fury.actor.scalar_bar(
            lookup_table=self.lut, title=self.color_label
        )
        self.scene.add(self.scalar_bar)
        self._update_overlay_colors()

    def _render_now(self):
        try:
            self.scene.ResetCameraClippingRange()
        except Exception:
            pass
        if self.showm is not None:
            self.showm.render()

    def _save_screenshot(self, out_path: Path, scale: int = 1) -> tuple[int, int]:
        """Capture the current viewer without creating another VTK window."""
        if self.showm is None or not hasattr(self.showm, "window"):
            raise RuntimeError("The interactive viewer is not initialized.")

        render_window = self.showm.window
        render_window.Render()

        capture = vtk.vtkWindowToImageFilter()
        capture.SetInput(render_window)
        capture.SetInputBufferTypeToRGB()
        if scale > 1:
            capture.SetScale(int(scale))
        capture.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName(str(out_path))
        writer.SetInputConnection(capture.GetOutputPort())
        writer.Write()

        width, height, _ = capture.GetOutput().GetDimensions()
        return int(width), int(height)

    def _save_session(self, announce: bool = True) -> Path:
        if self.session_path is None:
            self.session_path = (Path.cwd() / "streamlines_session.json").resolve()

        camera = self.scene.GetActiveCamera()
        origin = [0.0, 0.0, 0.0]
        normal = [1.0, 0.0, 0.0]
        self.plane_rep.GetOrigin(origin)
        self.plane_rep.GetNormal(normal)

        window_size = list(self.window_size)
        if self.showm is not None and hasattr(self.showm, "window"):
            window_size = [int(value) for value in self.showm.window.GetSize()]

        settings = dict(self.session_settings)
        settings.update(
            {
                "line_width": float(self.linewidth),
                "tube_thickness": float(self.linewidth),
                "window_size": window_size,
                "background_color": list(self.current_bg),
                "opacity": float(self.material["opacity"]),
                "quality": self.quality,
            }
        )
        box_transform = None
        if self.box_rep is not None:
            transform = vtk.vtkTransform()
            self.box_rep.GetTransform(transform)
            matrix = transform.GetMatrix()
            box_transform = [
                matrix.GetElement(row, column)
                for row in range(4)
                for column in range(4)
            ]
        payload = {
            "format": "cardiotensor-fury-session",
            "version": 1,
            "updated_at": datetime.datetime.now().isoformat(timespec="seconds"),
            "streamlines_file": (
                str(self.streamlines_file) if self.streamlines_file else None
            ),
            "streamlines_size": (
                self.streamlines_file.stat().st_size
                if self.streamlines_file and self.streamlines_file.exists()
                else None
            ),
            "settings": settings,
            "view": {
                "camera": {
                    "position": list(camera.GetPosition()),
                    "focal_point": list(camera.GetFocalPoint()),
                    "view_up": list(camera.GetViewUp()),
                    "clipping_range": list(camera.GetClippingRange()),
                    "view_angle": float(camera.GetViewAngle()),
                    "parallel_projection": bool(camera.GetParallelProjection()),
                    "parallel_scale": float(camera.GetParallelScale()),
                },
                "clipping_plane": {
                    "origin": origin,
                    "normal": normal,
                    "enabled": bool(self.clipping_active),
                    "gizmo_visible": bool(
                        self.plane_widget and self.plane_widget.GetEnabled()
                    ),
                },
                "crop_box": {
                    "transform": box_transform,
                    "enabled": bool(self.box_clipping_active),
                    "gizmo_visible": bool(
                        self.box_widget and self.box_widget.GetEnabled()
                    ),
                },
                "material": dict(self.material),
                "controls_visible": bool(self.controls_visible),
                "background_color": list(self.current_bg),
                "scale_bar_visible": bool(self.scale_bar_on),
                "window_size": window_size,
            },
        }

        self.session_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = self.session_path.with_suffix(
            self.session_path.suffix + ".tmp"
        )
        temporary_path.write_text(json.dumps(payload, indent=2) + "\n")
        temporary_path.replace(self.session_path)
        self.session_settings = settings
        if announce:
            print(f"Saved FURY session to {self.session_path}")
        return self.session_path

    def _restore_session(self) -> bool:
        if self.session_path is None or not self.session_path.exists():
            return False

        try:
            view = json.loads(self.session_path.read_text()).get("view", {})
            camera_state = view.get("camera", {})
            camera = self.scene.GetActiveCamera()
            if camera_state.get("position"):
                camera.SetPosition(*camera_state["position"])
            if camera_state.get("focal_point"):
                camera.SetFocalPoint(*camera_state["focal_point"])
            if camera_state.get("view_up"):
                camera.SetViewUp(*camera_state["view_up"])
            if camera_state.get("view_angle") is not None:
                camera.SetViewAngle(float(camera_state["view_angle"]))
            if camera_state.get("parallel_projection") is not None:
                camera.SetParallelProjection(bool(camera_state["parallel_projection"]))
            if camera_state.get("parallel_scale") is not None:
                camera.SetParallelScale(float(camera_state["parallel_scale"]))
            if camera_state.get("clipping_range"):
                camera.SetClippingRange(*camera_state["clipping_range"])

            plane = view.get("clipping_plane", {})
            if plane.get("origin"):
                self.plane_rep.SetOrigin(*plane["origin"])
            if plane.get("normal"):
                self.plane_rep.SetNormal(*plane["normal"])
            self.plane_rep.UpdatePlacement()
            self._sync_plane_from_widget()

            self.clipping_active = bool(plane.get("enabled", False))

            if self.plane_widget is not None:
                if plane.get("gizmo_visible", False):
                    self.plane_widget.EnabledOn()
                else:
                    self.plane_widget.EnabledOff()

            crop_box = view.get("crop_box", {})
            transform_values = crop_box.get("transform")
            if getattr(self, "box_rep", None) is not None and transform_values:
                matrix = vtk.vtkMatrix4x4()
                for row in range(4):
                    for column in range(4):
                        matrix.SetElement(
                            row, column, transform_values[row * 4 + column]
                        )
                transform = vtk.vtkTransform()
                transform.SetMatrix(matrix)
                self.box_rep.SetTransform(transform)
            self.box_clipping_active = bool(crop_box.get("enabled", False))
            if getattr(self, "box_widget", None) is not None:
                if crop_box.get("gizmo_visible", False):
                    self.box_widget.EnabledOn()
                else:
                    self.box_widget.EnabledOff()
            self._apply_clipping_planes()

            saved_material = view.get("material", {})
            if saved_material:
                for name in self.material:
                    if name in saved_material:
                        self.material[name] = float(saved_material[name])
                self._apply_material()
                for name, slider in self.control_sliders.items():
                    slider.value = self.material[name]
            self.controls_visible = bool(view.get("controls_visible", False))
            if getattr(self, "controls_panel", None) is not None:
                self.controls_panel.set_visibility(self.controls_visible)

            background = view.get("background_color")
            if background:
                self.current_bg = tuple(float(value) for value in background)
                self.scene.SetBackground(*self.current_bg)
                self._update_overlay_colors()

            show_scale_bar = bool(view.get("scale_bar_visible", True))
            if show_scale_bar and not self.scale_bar_on:
                self.scene.add(self.scale_bar)
                self.scale_bar_on = True
            elif not show_scale_bar and self.scale_bar_on:
                self.scene.rm(self.scale_bar)
                self.scale_bar_on = False

            if camera_state.get("clipping_range"):
                camera.SetClippingRange(*camera_state["clipping_range"])
            if self.showm is not None:
                self.showm.render()
            print(f"Restored FURY session from {self.session_path}")
            return True
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as err:
            print(f"Warning: could not restore session {self.session_path}: {err}")
            return False

    def _autosave_session(self, *_):
        if self.session_path is None:
            return
        try:
            self._save_session(announce=False)
            print(f"Updated FURY session: {self.session_path}")
        except (OSError, ValueError, TypeError) as err:
            print(f"Warning: could not update session {self.session_path}: {err}")

    def _close_window(self, *_):
        """Save once, then restore VTK's normal window-close behavior."""
        if self._closing:
            return
        self._closing = True
        self._autosave_session()
        for widget in (self.plane_widget, self.box_widget):
            if widget is not None:
                widget.EnabledOff()
        if self.showm is not None:
            self.showm.exit()

    def _sync_plane_from_widget(self, *_):
        origin = [0.0, 0.0, 0.0]
        normal = [1.0, 0.0, 0.0]
        self.plane_rep.GetOrigin(origin)
        self.plane_rep.GetNormal(normal)
        self.plane_fn.SetOrigin(origin)
        self.plane_fn.SetNormal(normal)
        if self.clipping_active:
            self._apply_clipping_planes()

    def _setup_box_widget(self):
        self.box_rep = vtk.vtkBoxRepresentation()
        self.box_rep.SetPlaceFactor(1.0)
        self.box_rep.PlaceWidget(self.bounds)
        try:
            self.box_rep.GetOutlineProperty().SetColor(0.2, 0.8, 1.0)
            self.box_rep.GetHandleProperty().SetColor(1.0, 0.75, 0.2)
        except Exception:
            pass

        self.box_widget = vtk.vtkBoxWidget2()
        self.box_widget.SetRepresentation(self.box_rep)
        self.box_widget.SetInteractor(self.showm.iren)
        self.box_widget.RotationEnabledOn()
        self.box_widget.EnabledOff()
        for event in (
            vtk.vtkCommand.StartInteractionEvent,
            vtk.vtkCommand.InteractionEvent,
            vtk.vtkCommand.EndInteractionEvent,
        ):
            self.box_widget.AddObserver(event, self._sync_box_from_widget)

    def _sync_box_from_widget(self, *_):
        if self.box_rep is not None:
            self.box_rep.GetPlanes(self.box_planes)
        if self.box_clipping_active:
            self._apply_clipping_planes()

    def _apply_clipping_planes(self):
        box_rep = getattr(self, "box_rep", None)
        box_active = bool(getattr(self, "box_clipping_active", False))
        if box_rep is not None:
            box_rep.GetPlanes(self.box_planes)
        for current_actor in (self.actor0, self.actor_fast):
            if current_actor is None:
                continue
            mapper = current_actor.GetMapper()
            mapper.RemoveAllClippingPlanes()
            if self.clipping_active:
                mapper.AddClippingPlane(self.plane_fn)
            if box_active:
                for plane_index in range(self.box_planes.GetNumberOfPlanes()):
                    mapper.AddClippingPlane(self.box_planes.GetPlane(plane_index))

    def _toggle_box_clipping(self):
        self.box_clipping_active = not self.box_clipping_active
        if self.box_widget is not None:
            if self.box_clipping_active:
                self.box_widget.EnabledOn()
            else:
                self.box_widget.EnabledOff()
        self._apply_clipping_planes()
        self._render_now()
        state = "ON" if self.box_clipping_active else "OFF"
        print(f"Rotatable crop box {state}")

    def _toggle_box_gizmo(self):
        if self.box_widget is None:
            return
        if self.box_widget.GetEnabled():
            self.box_widget.EnabledOff()
            print("Crop-box gizmo hidden; box clipping state unchanged")
        else:
            self.box_widget.EnabledOn()
            print("Crop-box gizmo shown; drag faces or handles and rotate it")
        self._render_now()

    def _reset_box(self):
        if self.box_rep is None:
            return
        self.box_rep.PlaceWidget(self.bounds)
        self._sync_box_from_widget()
        self._render_now()
        print("Crop box reset to the full streamline bounds")

    def _add_scale_bar(self):
        self.scale_bar = vtk.vtkLegendScaleActor()
        self.scale_bar.LeftAxisVisibilityOff()
        self.scale_bar.TopAxisVisibilityOff()
        self.scale_bar.RightAxisVisibilityOff()
        self.scale_bar.BottomAxisVisibilityOn()
        try:
            self.scale_bar.SetNumberOfLabels(5)
            self.scale_bar.SetCornerOffset(5)
        except Exception:
            pass
        self.scene.add(self.scale_bar)
        self.scale_bar_on = True

    def _toggle_clipping(self):
        """Toggle clipping state."""
        self.clipping_active = not self.clipping_active
        self._apply_clipping_planes()
        print(f"Clipping plane {'ON' if self.clipping_active else 'OFF'}")
        self._render_now()

    def _rebuild_unclipped_actor(self):
        if self.actor0 is not None:
            try:
                self.scene.rm(self.actor0)
            except Exception:
                pass

        if self.mode == "tube":
            self.actor0 = actor.streamtube(
                self.streamlines_xyz,
                colors=self.flat_vals,
                linewidth=self.linewidth,
                spline_subdiv=self.actor_spline_subdiv,
                tube_sides=self.tube_sides,
                lookup_colormap=self.lut,
                **_supported_actor_kwargs(actor.streamtube, lod=False),
            )
        else:
            self.actor0 = actor.line(
                self.streamlines_xyz,
                colors=self.flat_vals,
                linewidth=self.linewidth,
                spline_subdiv=self.actor_spline_subdiv,
                lookup_colormap=self.lut,
                **_supported_actor_kwargs(actor.line, lod=False),
            )

        self.scene.add(self.actor0)
        self._style_streamline_actor()

        self.mapper0 = self.actor0.GetMapper()
        self._apply_clipping_planes()

    # ---------------------------
    # key handling
    # ---------------------------
    def _on_keypress(self, obj, evt):
        key = obj.GetKeySym().lower()
        control_pressed = (
            bool(obj.GetControlKey()) if hasattr(obj, "GetControlKey") else False
        )
        shift_pressed = (
            bool(obj.GetShiftKey()) if hasattr(obj, "GetShiftKey") else False
        )

        if key == "s" and control_pressed:
            try:
                self._save_session()
            except (OSError, ValueError, TypeError) as err:
                print(f"Failed to save session: {err}")
            return

        if key == "o":
            self._toggle_clipping()

        elif key == "c":
            self._toggle_box_clipping()

        elif key == "g":
            self._toggle_box_gizmo()

        elif key == "x":
            self._reset_box()

        elif key == "l":
            self._toggle_controls_panel()

        elif key in ("1", "kp_1"):
            self._set_camera_preset("front")

        elif key in ("2", "kp_2"):
            self._set_camera_preset("side")

        elif key in ("3", "kp_3"):
            self._set_camera_preset("top")

        elif key in ("4", "kp_4"):
            self._set_camera_preset("isometric")

        elif key in ("5", "kp_5"):
            self._toggle_projection()

        elif key == "h":
            if self.plane_widget:
                currently_on = self.plane_widget.GetEnabled()
                if currently_on:
                    self.plane_widget.EnabledOff()
                    print("Plane gizmo hidden, clipping state unchanged")
                else:
                    self.plane_widget.EnabledOn()
                    print("Plane gizmo shown")
                self._render_now()

        elif key == "i":
            n = [0.0, 0.0, 0.0]
            self.plane_rep.GetNormal(n)
            self.plane_rep.SetNormal(-n[0], -n[1], -n[2])
            self._sync_plane_from_widget()

        elif key == "b":
            self.current_bg = (
                (0.0, 0.0, 0.0)
                if self.current_bg == (1.0, 1.0, 1.0)
                else (1.0, 1.0, 1.0)
            )
            self.scene.SetBackground(*self.current_bg)
            self._update_overlay_colors()
            self._render_now()
            print(f"Background set to {self.current_bg}")

        elif key == "s":
            if self.scale_bar is None:
                self._add_scale_bar()
            else:
                if self.scale_bar_on:
                    self.scene.rm(self.scale_bar)
                    self.scale_bar_on = False
                    print("Scale bar OFF")
                else:
                    self.scene.add(self.scale_bar)
                    self.scale_bar_on = True
                    print("Scale bar ON")
            self._render_now()

        elif key == "p":
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            scale = 2 if shift_pressed else 1
            suffix = "_2x" if scale == 2 else ""
            out_path = Path.cwd() / f"view_{ts}{suffix}.png"
            try:
                width, height = self._save_screenshot(out_path, scale=scale)
                print(
                    f"Saved screenshot to {out_path} "
                    f"({width}x{height}, {scale}x viewer resolution)"
                )
            except Exception as e:
                print(f"Failed to save screenshot: {e}")

        elif key == "v":
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            self._record_orbit(Path.cwd() / f"fury_orbit_{ts}.mp4", 120, 30)

        elif key == "r":
            self.plane_rep.SetOrigin(*self.center)
            self.plane_rep.SetNormal(1, 0, 0)
            self.plane_rep.UpdatePlacement()
            self._sync_plane_from_widget()
            print("Cropping plane reset to center, normal +X")

        elif key in ("plus", "equal", "kp_add"):
            self.linewidth = min(1000.0, self.linewidth * 1.25)
            if self.mode == "tube":
                self._rebuild_unclipped_actor()
            else:
                self.actor0.GetProperty().SetLineWidth(self.linewidth)
            self._render_now()
            print(f"Thickness up, lw={self.linewidth:.2f}")

        elif key in ("minus", "kp_subtract", "underscore"):
            self.linewidth = max(0.001, self.linewidth * 0.8)
            if self.mode == "tube":
                self._rebuild_unclipped_actor()
            else:
                self.actor0.GetProperty().SetLineWidth(self.linewidth)
            self._render_now()
            print(f"Thickness down, lw={self.linewidth:.2f}")

    def _lod_on(self, obj=None, evt=None):
        """Activate low-res actor for interaction, and apply clipping only if it was previously enabled."""
        try:
            # Hide the full-res actor and show the low-res actor
            self.actor0.SetVisibility(False)
            self.actor_fast.SetVisibility(True)

            self._apply_clipping_planes()
        except Exception:
            pass
        self._render_now()

    def _lod_off(self, obj=None, evt=None):
        """Switch back to full-res actor, applying clipping if it was previously enabled."""
        try:
            self.actor_fast.SetVisibility(False)  # Hide the fast actor
            self.actor0.SetVisibility(True)  # Show the full-res actor

            self._apply_clipping_planes()
        except Exception:
            pass
        self._render_now()

    def _set_default_camera(self) -> None:
        self.scene.reset_camera()
        self.scene.azimuth(15)
        self.scene.elevation(10)
        self.scene.zoom(1.1)

    def _set_camera_preset(self, preset: str) -> None:
        camera = self.scene.GetActiveCamera()
        center = np.asarray(self.center, dtype=float)
        diagonal = max(
            1.0,
            float(
                np.linalg.norm(
                    [
                        self.bounds[1] - self.bounds[0],
                        self.bounds[3] - self.bounds[2],
                        self.bounds[5] - self.bounds[4],
                    ]
                )
            ),
        )
        directions = {
            "front": (np.array([0.0, -1.0, 0.0]), (0.0, 0.0, 1.0)),
            "side": (np.array([1.0, 0.0, 0.0]), (0.0, 0.0, 1.0)),
            "top": (np.array([0.0, 0.0, 1.0]), (0.0, 1.0, 0.0)),
            "isometric": (np.array([1.0, -1.0, 1.0]), (0.0, 0.0, 1.0)),
        }
        direction, view_up = directions[preset]
        direction /= np.linalg.norm(direction)
        camera.SetFocalPoint(*center)
        camera.SetPosition(*(center + 1.8 * diagonal * direction))
        camera.SetViewUp(*view_up)
        camera.OrthogonalizeViewUp()
        self.scene.ResetCameraClippingRange()
        self._render_now()
        print(f"Camera preset: {preset}")

    def _toggle_projection(self) -> None:
        camera = self.scene.GetActiveCamera()
        use_parallel = not bool(camera.GetParallelProjection())
        camera.SetParallelProjection(use_parallel)
        if use_parallel:
            size = np.array(
                [
                    self.bounds[1] - self.bounds[0],
                    self.bounds[3] - self.bounds[2],
                    self.bounds[5] - self.bounds[4],
                ]
            )
            camera.SetParallelScale(0.55 * float(np.max(size)))
        self._render_now()
        print(f"Projection: {'parallel' if use_parallel else 'perspective'}")

    def _record_orbit(
        self, video_path: str | Path, video_frames: int, video_fps: int
    ) -> None:
        out_path = Path(video_path).resolve()
        n_frames = max(1, int(video_frames))
        fps = max(1, int(video_fps))
        step_degrees = 360.0 / float(n_frames)

        try:
            self.actor_fast.SetVisibility(False)
            self.actor0.SetVisibility(True)
        except Exception:
            pass

        self._set_default_camera()
        print(
            f"Recording FURY orbit video: {n_frames} frames at {fps} fps -> {out_path}"
        )

        scratch_root = out_path.parent / ".cardiotensor_scratch"
        scratch_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix="cardiotensor_fury_orbit_", dir=scratch_root
        ) as tmpdir:
            tmpdir_path = Path(tmpdir)
            frame_paths = []
            try:
                from tqdm import tqdm

                frame_iter = tqdm(
                    range(n_frames), total=n_frames, unit="frame", desc="Recording"
                )
            except ImportError:
                frame_iter = range(n_frames)

            for frame_idx in frame_iter:
                frame_path = tmpdir_path / f"frame_{frame_idx:05d}.png"
                fury.window.record(
                    scene=self.scene,
                    out_path=str(frame_path),
                    size=self.window_size,
                    reset_camera=False,
                )
                frame_paths.append(frame_path)
                self.scene.azimuth(step_degrees)

            _write_frame_sequence_to_video(frame_paths, out_path, fps)

        print(f"Saved FURY orbit video to {out_path}")

    # ---------------------------
    # main entry
    # ---------------------------
    def run(
        self,
        interactive: bool,
        screenshot_path: str | None,
        video_path: str | None = None,
        video_fps: int = 30,
        video_frames: int = 120,
    ):
        if interactive:
            self.showm = window.ShowManager(
                scene=self.scene, size=self.window_size, reset_camera=False
            )
            self.showm.initialize()
            self._build_controls_panel()

            iren = self.showm.iren
            iren.SetDesiredUpdateRate(60.0)
            iren.SetStillUpdateRate(30.0)  # higher so it returns to full-res

            # ensure anti-alias looks good when idle
            try:
                self.showm.window.SetMultiSamples(0)
                self.scene.enable_anti_aliasing("fxaa")
            except Exception:
                pass

            # switch LOD on drag only
            iren.AddObserver(vtk.vtkCommand.StartInteractionEvent, self._lod_on)
            iren.AddObserver(vtk.vtkCommand.EndInteractionEvent, self._lod_off)

            self.plane_widget = vtk.vtkImplicitPlaneWidget2()
            self.plane_widget.SetRepresentation(self.plane_rep)
            self.plane_widget.SetInteractor(self.showm.iren)
            self.plane_widget.EnabledOff()  # start with gizmo deactivated

            self.plane_widget.AddObserver(
                vtk.vtkCommand.StartInteractionEvent, self._sync_plane_from_widget
            )
            self.plane_widget.AddObserver(
                vtk.vtkCommand.InteractionEvent, self._sync_plane_from_widget
            )
            self.plane_widget.AddObserver(
                vtk.vtkCommand.EndInteractionEvent, self._sync_plane_from_widget
            )
            self._setup_box_widget()

            self.showm.iren.AddObserver("KeyPressEvent", self._on_keypress)
            # An ExitEvent observer replaces VTK's default TerminateApp callback,
            # so our handler must save the session and explicitly end interaction.
            self._closing = False
            self.showm.iren.AddObserver("ExitEvent", self._close_window)

            self._set_default_camera()
            if self.restore_session:
                self._restore_session()

            print(
                "Keys: O plane clipping, H plane gizmo, I flip plane, R reset plane; "
                "C rotatable crop box, G box gizmo, X reset box; L material controls; "
                "1 front, 2 side, 3 top, 4 isometric, 5 projection; +/- thickness; "
                "B background, S scale bar, Ctrl+S session, P PNG, Shift+P 2x PNG, "
                "V orbit"
            )
            self.showm.start()
        else:
            if not screenshot_path and not video_path:
                raise ValueError(
                    "Must specify screenshot_path or video_path when interactive=False."
                )
            self._set_default_camera()
            if screenshot_path:
                screenshot = Path(screenshot_path).resolve()
                screenshot.parent.mkdir(parents=True, exist_ok=True)
                fury.window.record(
                    scene=self.scene,
                    out_path=str(screenshot),
                    size=self.window_size,
                    reset_camera=False,
                )
            if video_path:
                self._record_orbit(video_path, video_frames, video_fps)


# ===========================
# Public API
# ===========================
def show_streamlines(
    streamlines_xyz: list[np.ndarray],
    color_values: list[np.ndarray],
    mode: str = "tube",
    line_width: float = 4,
    interactive: bool = True,
    screenshot_path: str | None = None,
    video_path: str | None = None,
    video_fps: int = 30,
    video_frames: int = 120,
    window_size: tuple[int, int] = (800, 800),
    downsample_factor: int = 2,
    max_streamlines: int | None = None,
    filter_min_len: int | None = None,
    subsample_factor: int = 1,
    crop_bounds: (
        tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None
    ) = None,
    colormap=None,
    background_color: str | tuple[float, float, float] = "black",
    spline_subdiv: int = 16,
    tube_sides: int = 9,
    color_range: tuple[float, float] | None = None,
    color_label: str = "Angle (deg)",
    random_seed: int | None = None,
    session_path: str | Path | None = None,
    restore_session: bool = True,
    session_settings: dict | None = None,
    streamlines_file: str | Path | None = None,
    opacity: float = 1.0,
    quality: str = "interactive",
):
    print(f"Initial number of streamlines: {len(streamlines_xyz)}")
    full_mins, full_maxs = _compute_streamline_bounds(streamlines_xyz)
    _print_box_shape("Full streamline box", full_mins, full_maxs)

    if crop_bounds is not None:
        (x_min, x_max), (y_min, y_max), (z_min, z_max) = crop_bounds
        print(f"Cropping streamlines within bounds: {crop_bounds}")
        _print_box_shape(
            "Requested crop box",
            np.array([x_min, y_min, z_min], dtype=float),
            np.array([x_max, y_max, z_max], dtype=float),
        )
        new_streamlines, new_colors = [], []
        for sl, cl in zip(streamlines_xyz, color_values):
            segs, cols = _split_streamline_by_bounds(
                sl, cl, x_min, x_max, y_min, y_max, z_min, z_max
            )
            if segs:
                new_streamlines.extend(segs)
                new_colors.extend(cols)
        streamlines_xyz, color_values = new_streamlines, new_colors
        if not streamlines_xyz:
            raise ValueError("No streamlines intersect the crop box.")
        cropped_mins, cropped_maxs = _compute_streamline_bounds(streamlines_xyz)
        _print_box_shape("Cropped streamline box", cropped_mins, cropped_maxs)
        print("Cropping applied.")
    else:
        print("No cropping applied.")

    print(f"Downsampling points by factor {downsample_factor}")
    if filter_min_len is not None:
        print(f"Filtering out streamlines shorter than {filter_min_len} points")

    ds_streamlines, ds_colors = [], []
    for sl, cl in zip(streamlines_xyz, color_values):
        ds_sl = downsample_streamline(sl, downsample_factor)
        ds_cl = downsample_streamline(cl, downsample_factor)
        if filter_min_len is None or len(ds_sl) >= filter_min_len:
            ds_streamlines.append(ds_sl)
            ds_colors.append(ds_cl)

    streamlines_xyz, color_values = ds_streamlines, ds_colors
    if not streamlines_xyz:
        raise ValueError("No streamlines left after downsampling or filtering.")

    rng = random.Random(random_seed)
    if subsample_factor > 1:
        print(f"Subsampling: keeping 1 in every {subsample_factor} streamlines")
        total = len(streamlines_xyz)
        keep_idx = sorted(rng.sample(range(total), max(1, total // subsample_factor)))
        streamlines_xyz = [streamlines_xyz[i] for i in keep_idx]
        color_values = [color_values[i] for i in keep_idx]

    if max_streamlines is not None and len(streamlines_xyz) > max_streamlines:
        print(f"Limiting to max {max_streamlines} streamlines")
        keep_idx = sorted(rng.sample(range(len(streamlines_xyz)), max_streamlines))
        streamlines_xyz = [streamlines_xyz[i] for i in keep_idx]
        color_values = [color_values[i] for i in keep_idx]

    print(f"Final number of streamlines to render: {len(streamlines_xyz)}")
    points_per_streamline = np.array(
        [len(streamline) for streamline in streamlines_xyz], dtype=float
    )
    avg_points_per_streamline = float(np.mean(points_per_streamline))
    std_points_per_streamline = float(np.std(points_per_streamline))
    print(
        f"Average points per streamline: {avg_points_per_streamline:.1f} +/- "
        f"{std_points_per_streamline:.1f} "
        "(more points usually means more detail, but heavier rendering)."
    )

    if not color_values:
        raise ValueError("No color arrays after filtering.")

    flat_colors = np.concatenate([np.asarray(c).ravel() for c in color_values]).astype(
        np.float32
    )
    min_val = float(np.nanmin(flat_colors))
    max_val = float(np.nanmax(flat_colors))
    print(f"Coloring range: min={min_val:.3f}, max={max_val:.3f}")
    if color_range is None:
        lut_range = (min_val, max_val)
    else:
        lut_range = (float(color_range[0]), float(color_range[1]))
        print(f"Colorbar range: min={lut_range[0]:.3f}, max={lut_range[1]:.3f}")
    print(f"Rendering mode: {mode}")
    if int(spline_subdiv) <= 1:
        print("Visual smoothing: off (original streamline points preserved)")
    else:
        print(
            "Visual smoothing: VTK cardinal spline with "
            f"{int(spline_subdiv)} segments per streamline"
        )
    print(f"FURY quality preset: {quality}")

    if colormap is None:
        lut = fury.actor.colormap_lookup_table(
            scale_range=lut_range,
            hue_range=(0.7, 0.0),
            saturation_range=(0.5, 1.0),
        )
    else:
        lut = matplotlib_cmap_to_fury_lut(
            cmap=colormap, value_range=lut_range, n_colors=256
        )
    lut.SetRange(*lut_range)

    viewer = StreamlineViewer(
        streamlines_xyz,
        color_values,
        mode,
        line_width,
        window_size,
        lut,
        background_color=background_color,
        spline_subdiv=spline_subdiv,
        tube_sides=tube_sides,
        color_range=lut_range,
        color_label=color_label,
        flat_values=flat_colors,
        session_path=session_path,
        restore_session=restore_session,
        session_settings=session_settings,
        streamlines_file=streamlines_file,
        opacity=opacity,
        quality=quality,
    )
    viewer.run(
        interactive=interactive,
        screenshot_path=screenshot_path,
        video_path=video_path,
        video_fps=video_fps,
        video_frames=video_frames,
    )


# ---------------------------
# Smoke test
# ---------------------------
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    sls, cols = [], []
    for k in range(50):
        t = np.linspace(0, 2 * np.pi, 150)
        r = 20 + 2 * rng.normal()
        x = r * np.cos(t) + 5 * rng.normal()
        y = r * np.sin(t) + 5 * rng.normal()
        z = np.linspace(-30, 30, t.size) + 3 * rng.normal()
        sl = np.c_[x, y, z].astype(np.float32)
        c = (np.degrees(np.arctan2(y, x))).astype(np.float32)
        sls.append(sl)
        cols.append(c)

    show_streamlines(
        streamlines_xyz=sls,
        color_values=cols,
        mode="tube",
        line_width=4,
        interactive=True,
        window_size=(900, 900),
        downsample_factor=1,
        subsample_factor=1,
        max_streamlines=None,
        filter_min_len=None,
        crop_bounds=None,
        colormap="turbo",
    )
