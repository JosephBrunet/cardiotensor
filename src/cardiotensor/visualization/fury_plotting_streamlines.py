#!/usr/bin/env python3
from __future__ import annotations

import datetime
import json
import random
import tempfile
from pathlib import Path

import numpy as np
import pygfx as gfx
import pylinalg as la
from fury import actor, ui, window

from cardiotensor.visualization.fury_support import (
    ScalarColors,
    camera_for_scene,
    capture,
    materials,
    scalar_bar,
)

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
    """Interactive FURY 2 viewer using native PyGfx actors and transform gizmos."""

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
        tube_sides=9,
        color_range=None,
        color_label="Angle",
        flat_values=None,
        session_path=None,
        restore_session=True,
        session_settings=None,
        streamlines_file=None,
        opacity=1.0,
        quality="interactive",
        highlight_flags=None,
        other_opacity=0.05,
    ):
        if mode not in {"line", "tube"}:
            raise ValueError("mode must be 'line' or 'tube'")
        if quality not in {"interactive", "publication"}:
            raise ValueError("quality must be 'interactive' or 'publication'")
        if highlight_flags is not None and len(highlight_flags) != len(streamlines_xyz):
            raise ValueError("highlight_flags must match streamlines")
        if len(streamlines_xyz) != len(color_values):
            raise ValueError("color_values must match streamlines")
        self.streamlines_xyz = streamlines_xyz
        self.color_values = color_values
        self.mode, self.window_size, self.lut = mode, tuple(window_size), lut
        self.color_label = color_label
        self.session_path = (
            Path(session_path).expanduser().resolve() if session_path else None
        )
        self.streamlines_file = (
            Path(streamlines_file).resolve() if streamlines_file else None
        )
        self.restore_session = restore_session
        self.session_settings = dict(session_settings or {})
        self.quality = quality
        self.linewidth = max(0.001, float(line_width))
        self.spline_subdiv = max(1, int(spline_subdiv))
        self.actor_spline_subdiv = (
            self.spline_subdiv if self.spline_subdiv > 1 else None
        )
        self.tube_sides = max(3, int(tube_sides))
        self.highlight_flags = highlight_flags
        self.other_opacity = float(other_opacity)
        self.material = dict(
            ambient=0.45, diffuse=0.70, specular=0.12, opacity=float(opacity)
        )
        if quality == "publication":
            self.material.update(ambient=0.30, diffuse=0.80, specular=0.25)
        self.current_bg = parse_background_color(background_color)
        self.scene = window.Scene(background=(*self.current_bg, 1))
        self.showm = None
        self._closing = False
        self.controls_panel = None
        self.controls_visible = False
        self.control_sliders = {}
        mins, maxs = _compute_streamline_bounds(streamlines_xyz)
        self.center = (mins + maxs) / 2
        self.extent = np.maximum(maxs - mins, 1e-3)
        self.bounds = np.column_stack((mins, maxs)).ravel().tolist()
        self.flat_vals = (
            np.concatenate(color_values) if flat_values is None else flat_values
        )
        self.vmin, self.vmax = color_range or (
            float(np.nanmin(self.flat_vals)),
            float(np.nanmax(self.flat_vals)),
        )
        self.lut.value_range = (self.vmin, self.vmax)
        self._prepare_geometry()
        self.actor0 = self.actor_highlight = None
        self.actor_fast = self._make_actor(range(len(self.render_lines)), fast=True)
        self.actor_fast.visible = False
        self.scene.add(self.actor_fast)
        self.clipping_active = self.box_clipping_active = False
        self.plane_widget = self.box_widget = None
        self._setup_clip_objects()
        self._rebuild_unclipped_actor()
        self.scalar_bar = scalar_bar(self.lut, self.color_label)
        self.scene.add(self.scalar_bar)
        self.scale_bar = gfx.Ruler(
            start_pos=(0, 0, 0),
            end_pos=(float(self.extent.max()) / 4, 0, 0),
            ticks=None,
            tick_format="0.3g",
            color="white",
        )
        self.scale_bar.local.position = mins
        self.scale_bar.color = "black" if sum(self.current_bg) > 1.5 else "white"
        self.scale_bar_on = True
        self.scene.add(self.scale_bar)
        self.camera = camera_for_scene(self.scene, self.window_size)
        self.ambient_light = gfx.AmbientLight(intensity=self.material["ambient"])
        self.scene.add(self.ambient_light)
        self._set_default_camera()

    def _prepare_geometry(self):
        """Resample geometry and scalar values together before RGB mapping."""
        from scipy.interpolate import CubicSpline

        self.render_lines, self.render_colors = [], []
        for line, values in zip(self.streamlines_xyz, self.color_values):
            line = np.asarray(line, dtype=np.float32)
            values = np.asarray(values).ravel()
            if len(line) != len(values):
                raise ValueError("Each streamline needs one scalar per point")
            if len(line) < 2:
                continue
            if self.spline_subdiv > 1:
                distances = np.r_[
                    0, np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))
                ]
                keep = np.r_[True, np.diff(distances) > 0]
                distances, points, scalars = distances[keep], line[keep], values[keep]
                if len(points) < 2:
                    continue
                samples = np.linspace(0, distances[-1], self.spline_subdiv + 1)
                line = CubicSpline(distances, points, bc_type="natural")(
                    samples
                ).astype(np.float32)
                values = np.interp(samples, distances, scalars)
            self.render_lines.append(line)
            self.render_colors.append(self.lut.map(values))
        if len(self.render_lines) != len(self.streamlines_xyz):
            raise ValueError("Streamlines must have at least two distinct points")

    def _make_actor(self, indices, fast=False):
        indices = list(indices)
        lines = [self.render_lines[i] for i in indices]
        colors = [self.render_colors[i] for i in indices]
        if self.mode == "tube" and not fast:
            obj = actor.streamtube(
                lines,
                colors=colors,
                radius=self.linewidth,
                segments=self.tube_sides,
                backend="gpu",
            )
        else:
            obj = actor.line(lines, colors=np.concatenate(colors))
            # FURY 2.0 line() adds the first position to already-world-space data.
            obj.local.position = (0, 0, 0)
            for material in materials(obj):
                material.thickness = 1.0 if fast else self.linewidth
        return obj

    def _rebuild_unclipped_actor(self):
        for obj in (self.actor0, self.actor_highlight):
            if obj is not None:
                self.scene.remove(obj)
        self.actor0 = self._make_actor(range(len(self.render_lines)))
        self.scene.add(self.actor0)
        self.actor_highlight = None
        if self.highlight_flags is not None:
            selected = np.flatnonzero(self.highlight_flags)
            if len(selected):
                self.actor_highlight = self._make_actor(selected)
                self.scene.add(self.actor_highlight)
        self._apply_material()
        self._apply_clipping_planes()

    def _apply_material(self):
        for obj in (self.actor0, self.actor_fast, self.actor_highlight):
            if obj is None:
                continue
            opacity = self.material["opacity"]
            if self.highlight_flags is not None:
                opacity = 1.0 if obj is self.actor_highlight else self.other_opacity
            for mat in materials(obj):
                mat.opacity = opacity
                if isinstance(mat, gfx.MeshPhongMaterial):
                    mat.shininess = 20 if self.quality == "publication" else 10
                    mat.specular = (self.material["specular"],) * 3
                    mat.color = (self.material["diffuse"],) * 3
        if hasattr(self, "ambient_light"):
            self.ambient_light.intensity = self.material["ambient"]
        self._render_now()

    def _set_material_value(self, name, value):
        self.material[name] = float(value)
        self._apply_material()

    def _build_controls_panel(self):
        panel = ui.Panel2D(size=(280, 255), position=(12, 95), opacity=0.9)
        panel.add_element(
            ui.TextBlock2D(text="Material controls (L)", size=(260, 25)), (10, 225)
        )
        for row, name in enumerate(self.material):
            slider = ui.LineSlider2D(
                initial_value=self.material[name],
                min_value=0.0,
                max_value=1.0,
                length=140,
                text_template="{value:.2f}",
            )
            slider.on_change = lambda current, key=name: self._set_material_value(
                key, current.value
            )
            panel.add_element(
                ui.TextBlock2D(text=name.title(), size=(95, 25), font_size=13),
                (10, 180 - row * 48),
            )
            panel.add_element(slider, (175, 190 - row * 48), anchor="center")
            self.control_sliders[name] = slider
        self.controls_panel = panel
        self.scene.add(panel)
        panel.set_visibility(self.controls_visible)

    def _setup_clip_objects(self):
        self.plane_object = gfx.Mesh(
            gfx.plane_geometry(float(self.extent.max()), float(self.extent.max())),
            gfx.MeshBasicMaterial(color=(1, 0.7, 0.2), opacity=0.2, side="both"),
        )
        self.box_object = gfx.Mesh(
            gfx.box_geometry(*self.extent),
            gfx.MeshBasicMaterial(color=(0.2, 0.8, 1), wireframe=True, opacity=0.5),
        )
        self.scene.add(self.plane_object, self.box_object)
        self._reset_plane()
        self._reset_box()
        self.plane_object.visible = self.box_object.visible = False

    def _setup_widgets(self):
        for obj, name in (
            (self.plane_object, "plane_widget"),
            (self.box_object, "box_widget"),
        ):
            widget = gfx.TransformGizmo(obj)
            widget.add_default_event_handlers(self.showm.renderer, self.camera)
            widget.visible = False
            self.scene.add(widget)
            setattr(self, name, widget)
        self.showm.renderer.add_event_handler(self._sync_widgets, "before_render")
        self.showm.renderer.add_event_handler(self._lod_on, "pointer_down")
        self.showm.renderer.add_event_handler(
            self._lod_off, "pointer_up", "pointer_leave"
        )

    def _lod_on(self, event=None):
        # Keep highlighted bundles visible with their correct opacity during drag.
        if self.quality == "interactive" and self.highlight_flags is None:
            self.actor0.visible = False
            self.actor_fast.visible = True
            self._render_now()

    def _lod_off(self, event=None):
        self.actor0.visible = True
        self.actor_fast.visible = False
        self._render_now()

    def _sync_widgets(self, event=None):
        self._apply_clipping_planes()
        if self.scale_bar_on:
            self.scale_bar.update(self.camera, self.showm.size)

    def _plane_equation(self):
        # PyGfx clips dot(normal, point) < distance.
        matrix = self.plane_object.world.matrix
        normal = np.linalg.inv(matrix[:3, :3]).T @ np.array([0.0, 0.0, 1.0])
        normal /= np.linalg.norm(normal)
        return (*normal, float(normal @ matrix[:3, 3]))

    def _box_equations(self):
        inverse = np.linalg.inv(self.box_object.world.matrix)
        planes = []
        for axis in range(3):
            for sign in (-1, 1):
                normal = sign * inverse[axis, :3]
                constant = sign * inverse[axis, 3] + self.extent[axis] / 2
                length = np.linalg.norm(normal)
                planes.append((*(normal / length), float(-constant / length)))
        return planes

    def _apply_clipping_planes(self):
        planes = [self._plane_equation()] if self.clipping_active else []
        if self.box_clipping_active:
            planes.extend(self._box_equations())
        for obj in (self.actor0, self.actor_fast, self.actor_highlight):
            if obj is not None:
                for material in materials(obj):
                    material.clipping_planes = planes

    def _set_gizmo_visible(self, name, visible):
        widget = getattr(self, name + "_widget")
        if widget is not None:
            widget.visible = visible
        getattr(self, name + "_object").visible = visible
        self._render_now()

    def _reset_plane(self):
        self.plane_object.local.position = self.center
        self.plane_object.local.rotation = la.quat_from_axis_angle((0, 1, 0), np.pi / 2)
        self.plane_object.local.scale = (1, 1, 1)

    def _reset_box(self):
        self.box_object.local.position = self.center
        self.box_object.local.rotation = (0, 0, 0, 1)
        self.box_object.local.scale = (1, 1, 1)

    def _set_default_camera(self):
        self._set_camera_preset("isometric")

    def _set_camera_preset(self, preset):
        directions = {
            "front": (0, -1, 0),
            "side": (1, 0, 0),
            "top": (0, 0, 1),
            "isometric": (1, -1, 1),
        }
        direction = np.array(directions[preset], dtype=float)
        direction /= np.linalg.norm(direction)
        self.camera.local.position = self.center + direction * 1.8 * np.linalg.norm(
            self.extent
        )
        self.camera.world.reference_up = (0, 1, 0) if preset == "top" else (0, 0, 1)
        self.camera.look_at(self.center)
        self.camera.width = self.camera.height = float(self.extent.max()) * 1.2
        self._render_now()

    def _toggle_projection(self):
        self.camera.fov = 45 if self.camera.fov == 0 else 0
        self._render_now()

    def _render_now(self):
        if self.showm is not None:
            self.showm.render()

    def _save_screenshot(self, out_path, scale=1):
        if self.showm is None:
            raise RuntimeError("The viewer is not initialized")
        self._lod_off()
        self._sync_widgets()
        return capture(self.showm, out_path, scale)

    def _save_session(self, announce=True):
        if self.session_path is None:
            self.session_path = Path.cwd() / "streamlines_session.json"
        size = list(self.showm.size if self.showm else self.window_size)
        settings = dict(self.session_settings)
        settings.update(
            line_width=self.linewidth,
            tube_thickness=self.linewidth,
            window_size=size,
            background_color=list(self.current_bg),
            opacity=self.material["opacity"],
            quality=self.quality,
        )
        payload = {
            "format": "cardiotensor-fury-session",
            "version": 2,
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
                "camera": self.camera.get_state(),
                "clipping_plane": {
                    "matrix": self.plane_object.local.matrix.tolist(),
                    "enabled": self.clipping_active,
                    "gizmo_visible": self.plane_object.visible,
                },
                "crop_box": {
                    "matrix": self.box_object.local.matrix.tolist(),
                    "enabled": self.box_clipping_active,
                    "gizmo_visible": self.box_object.visible,
                },
                "material": self.material,
                "controls_visible": self.controls_visible,
                "background_color": list(self.current_bg),
                "scale_bar_visible": self.scale_bar_on,
                "window_size": size,
            },
        }
        self.session_path.parent.mkdir(parents=True, exist_ok=True)
        temp = self.session_path.with_suffix(self.session_path.suffix + ".tmp")
        temp.write_text(
            json.dumps(
                payload, indent=2, default=lambda value: np.asarray(value).tolist()
            )
            + "\n"
        )
        temp.replace(self.session_path)
        self.session_settings = settings
        if announce:
            print(f"Saved FURY session to {self.session_path}")
        return self.session_path

    def _restore_session(self):
        if self.session_path is None or not self.session_path.exists():
            return False
        try:
            payload = json.loads(self.session_path.read_text())
            view = payload.get("view", {})
            state = view.get("camera", {})
            if payload.get("version", 1) >= 2:
                self.camera.set_state(state)
            else:
                # Migrate the camera and clip transforms from VTK sessions.
                self.camera.local.position = state.get(
                    "position", self.camera.local.position
                )
                self.camera.world.reference_up = state.get("view_up", (0, 0, 1))
                self.camera.look_at(state.get("focal_point", self.center))
                self.camera.fov = (
                    0
                    if state.get("parallel_projection")
                    else state.get("view_angle", 45)
                )
                self.camera.height = 2 * state.get(
                    "parallel_scale", self.extent.max() / 2
                )
            for name, key in (("plane", "clipping_plane"), ("box", "crop_box")):
                saved = view.get(key, {})
                obj = getattr(self, name + "_object")
                if "matrix" in saved:
                    obj.local.matrix = np.asarray(saved["matrix"], dtype=float)
                elif name == "plane":
                    obj.local.position = saved.get("origin", self.center)
                    normal = np.asarray(saved.get("normal", (1, 0, 0)), dtype=float)
                    obj.local.rotation = la.quat_from_vecs((0, 0, 1), normal)
                elif saved.get("transform"):
                    transform = np.asarray(saved["transform"]).reshape(4, 4)
                    placement = np.eye(4)
                    placement[:3, 3] = self.center
                    obj.local.matrix = transform @ placement
                setattr(
                    self,
                    "clipping_active" if name == "plane" else "box_clipping_active",
                    bool(saved.get("enabled", False)),
                )
                self._set_gizmo_visible(name, bool(saved.get("gizmo_visible", False)))
            self.material.update(view.get("material", {}))
            self._apply_material()
            self._apply_clipping_planes()
            self.current_bg = tuple(view.get("background_color", self.current_bg))
            self.scene.background = (*self.current_bg, 1)
            self.scale_bar.color = "black" if sum(self.current_bg) > 1.5 else "white"
            self.scale_bar_on = bool(view.get("scale_bar_visible", True))
            self.scale_bar.visible = self.scale_bar_on
            self.controls_visible = bool(view.get("controls_visible", False))
            if self.controls_panel is not None:
                self.controls_panel.set_visibility(self.controls_visible)
                for name, slider in self.control_sliders.items():
                    slider.value = self.material[name]
            self._render_now()
            print(f"Restored FURY session from {self.session_path}")
            return True
        except (OSError, ValueError, TypeError, KeyError) as err:
            print(f"Warning: could not restore session {self.session_path}: {err}")
            return False

    def _autosave_session(self, *_):
        if self.session_path is not None:
            try:
                self._save_session(announce=False)
            except (OSError, ValueError, TypeError) as err:
                print(f"Warning: could not save session: {err}")

    def _close_window(self, *_):
        if self._closing:
            return
        self._closing = True
        self._autosave_session()
        if self.showm is not None:
            self.showm.close()

    def _on_keypress(self, event):
        key = event.key.lower()
        modifiers = {str(mod).lower() for mod in event.modifiers}
        if key == "s" and ("control" in modifiers or "ctrl" in modifiers):
            self._save_session()
        elif key == "o":
            self.clipping_active = not self.clipping_active
        elif key == "c":
            self.box_clipping_active = not self.box_clipping_active
            self._set_gizmo_visible("box", self.box_clipping_active)
        elif key in {"h", "g"}:
            name = "plane" if key == "h" else "box"
            self._set_gizmo_visible(name, not getattr(self, name + "_object").visible)
        elif key == "i":
            self.plane_object.local.rotation = la.quat_mul(
                self.plane_object.local.rotation,
                la.quat_from_axis_angle((1, 0, 0), np.pi),
            )
        elif key == "r":
            self._reset_plane()
        elif key == "x":
            self._reset_box()
        elif key == "l":
            self.controls_visible = not self.controls_visible
            if self.controls_panel is not None:
                self.controls_panel.set_visibility(self.controls_visible)
        elif key in {"1", "2", "3", "4"}:
            self._set_camera_preset(
                {"1": "front", "2": "side", "3": "top", "4": "isometric"}[key]
            )
        elif key == "5":
            self._toggle_projection()
        elif key in {"+", "=", "-", "_"}:
            self.linewidth = max(
                0.001, self.linewidth * (1.25 if key in {"+", "="} else 0.8)
            )
            self._rebuild_unclipped_actor()
        elif key == "b":
            self.current_bg = (0, 0, 0) if self.current_bg == (1, 1, 1) else (1, 1, 1)
            self.scene.background = (*self.current_bg, 1)
            self.scale_bar.color = "black" if sum(self.current_bg) > 1.5 else "white"
        elif key == "s":
            self.scale_bar_on = not self.scale_bar_on
            self.scale_bar.visible = self.scale_bar_on
        elif key == "p":
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            self._save_screenshot(
                Path.cwd() / f"view_{timestamp}.png", 2 if "shift" in modifiers else 1
            )
        elif key == "v":
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            self._record_orbit(Path.cwd() / f"fury_orbit_{timestamp}.mp4", 120, 30)
        self._apply_clipping_planes()
        self._render_now()

    def _record_orbit(self, video_path, video_frames, video_fps):
        out_path = Path(video_path).resolve()
        count, fps = max(1, int(video_frames)), max(1, int(video_fps))
        scratch = out_path.parent / ".cardiotensor_scratch"
        scratch.mkdir(parents=True, exist_ok=True)
        state = self.camera.get_state()
        offset = np.asarray(self.camera.local.position) - self.center
        try:
            with tempfile.TemporaryDirectory(
                prefix="cardiotensor_fury_orbit_", dir=scratch
            ) as tmpdir:
                paths = []
                for index in range(count):
                    angle = 2 * np.pi * index / count
                    rotation = np.array(
                        [
                            [np.cos(angle), -np.sin(angle), 0],
                            [np.sin(angle), np.cos(angle), 0],
                            [0, 0, 1],
                        ]
                    )
                    self.camera.local.position = self.center + rotation @ offset
                    self.camera.look_at(self.center)
                    path = Path(tmpdir) / f"frame_{index:05d}.png"
                    self._save_screenshot(path)
                    paths.append(path)
                _write_frame_sequence_to_video(paths, out_path, fps)
        finally:
            self.camera.set_state(state)
            self._render_now()
        print(f"Saved FURY orbit video to {out_path}")

    def run(
        self,
        interactive,
        screenshot_path,
        video_path=None,
        video_fps=30,
        video_frames=120,
    ):
        if not interactive and not screenshot_path and not video_path:
            raise ValueError(
                "Specify screenshot_path or video_path when interactive=False"
            )
        self.showm = window.ShowManager(
            scene=self.scene,
            camera=self.camera,
            size=self.window_size,
            window_type="default" if interactive else "offscreen",
        )
        self._set_default_camera()
        if interactive:
            self._build_controls_panel()
            self._setup_widgets()
            self.showm.renderer.add_event_handler(self._on_keypress, "key_down")
            self.showm.window.add_event_handler(self._close_window, "close")
        if self.restore_session:
            self._restore_session()
        self._sync_widgets()
        try:
            if screenshot_path:
                self._save_screenshot(screenshot_path)
            if video_path:
                self._record_orbit(video_path, video_frames, video_fps)
            if interactive:
                print(
                    "Keys: O clip plane, H plane gizmo, I flip, R reset plane; C crop box, G box gizmo, X reset box; L material controls; 1–4 camera, 5 projection; +/- thickness; B background, S scale, Ctrl+S save, P PNG, Shift+P 2x PNG, V orbit"
                )
                self.showm.start()
        finally:
            self._close_window()


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
    highlight_flags: list[bool] | None = None,
    other_opacity: float = 0.05,
):
    print(f"Initial number of streamlines: {len(streamlines_xyz)}")
    highlighting = highlight_flags is not None
    if highlight_flags is None:
        highlight_flags = [False] * len(streamlines_xyz)
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
        new_streamlines, new_colors, new_flags = [], [], []
        for sl, cl, flag in zip(streamlines_xyz, color_values, highlight_flags):
            segs, cols = _split_streamline_by_bounds(
                sl, cl, x_min, x_max, y_min, y_max, z_min, z_max
            )
            if segs:
                new_streamlines.extend(segs)
                new_colors.extend(cols)
                new_flags.extend([flag] * len(segs))
        streamlines_xyz, color_values = new_streamlines, new_colors
        highlight_flags = new_flags
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

    ds_streamlines, ds_colors, ds_flags = [], [], []
    for sl, cl, flag in zip(streamlines_xyz, color_values, highlight_flags):
        ds_sl = downsample_streamline(sl, downsample_factor)
        ds_cl = downsample_streamline(cl, downsample_factor)
        if filter_min_len is None or len(ds_sl) >= filter_min_len:
            ds_streamlines.append(ds_sl)
            ds_colors.append(ds_cl)
            ds_flags.append(flag)

    streamlines_xyz, color_values = ds_streamlines, ds_colors
    highlight_flags = ds_flags
    if not streamlines_xyz:
        raise ValueError("No streamlines left after downsampling or filtering.")

    rng = random.Random(random_seed)
    if subsample_factor > 1:
        print(f"Subsampling: keeping 1 in every {subsample_factor} streamlines")
        total = len(streamlines_xyz)
        keep_idx = sorted(rng.sample(range(total), max(1, total // subsample_factor)))
        streamlines_xyz = [streamlines_xyz[i] for i in keep_idx]
        color_values = [color_values[i] for i in keep_idx]
        highlight_flags = [highlight_flags[i] for i in keep_idx]

    if max_streamlines is not None and len(streamlines_xyz) > max_streamlines:
        print(f"Limiting to max {max_streamlines} streamlines")
        keep_idx = sorted(rng.sample(range(len(streamlines_xyz)), max_streamlines))
        streamlines_xyz = [streamlines_xyz[i] for i in keep_idx]
        color_values = [color_values[i] for i in keep_idx]
        highlight_flags = [highlight_flags[i] for i in keep_idx]

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
            "Visual smoothing: natural cubic spline with "
            f"{int(spline_subdiv)} segments per streamline"
        )
    print(f"FURY quality preset: {quality}")

    lut = ScalarColors(colormap if colormap is not None else "turbo", lut_range)

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
        highlight_flags=highlight_flags if highlighting else None,
        other_opacity=other_opacity,
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
