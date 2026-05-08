#!/usr/bin/env python3
from __future__ import annotations

import random

import fury
import matplotlib.pyplot as plt
import numpy as np
import vtk
from fury import actor, window

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
    elif isinstance(color, (tuple, list)) and len(color) == 3:
        return tuple(float(c) for c in color)
    else:
        raise TypeError("background_color must be str or 3-tuple")


def downsample_streamline(streamline: np.ndarray, factor: int = 2) -> np.ndarray:
    return streamline if len(streamline) < 3 or factor <= 1 else streamline[::factor]


def _compute_streamline_bounds(
    streamlines_xyz: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    if not streamlines_xyz:
        raise ValueError("No streamlines available to compute bounds.")

    mins = np.min([sl.min(axis=0) for sl in streamlines_xyz], axis=0)
    maxs = np.max([sl.max(axis=0) for sl in streamlines_xyz], axis=0)
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
        color_range: tuple[float, float] | None = None,
        color_label: str = "Angle",
    ):
        self.streamlines_xyz = streamlines_xyz
        self.color_values = color_values
        self.mode = mode
        self.window_size = window_size
        self.lut = lut
        self.color_label = color_label

        # Scene
        self.scene = fury.window.Scene()
        self.current_bg = parse_background_color(background_color)
        self.scene.SetBackground(*self.current_bg)

        # thickness state
        self.linewidth = max(1.0, float(line_width))  # used for both line and tube
        self.spline_subdiv = spline_subdiv

        self.scale_bar = None
        self.scale_bar_on = False
        self._add_scale_bar()

        self.clipping_active = False  # Track clipping state

        # VTK/FURY objects
        self.showm: window.ShowManager | None = None

        # clipped branch objects
        self.plane_rep = None
        self.plane_fn = None
        self.plane_widget = None
        self.mapper0 = None
        self.actor0 = None

        # precompute flat scalars for LUT mapping
        self.flat_vals = np.concatenate(
            [np.asarray(c).ravel() for c in self.color_values]
        ).astype(np.float32)
        if color_range is None:
            self.vmin = float(np.nanmin(self.flat_vals))
            self.vmax = float(np.nanmax(self.flat_vals))
        else:
            self.vmin, self.vmax = (float(color_range[0]), float(color_range[1]))
        self.lut.SetRange(self.vmin, self.vmax)

        # bounds and center from NumPy
        mins = np.min([sl.min(axis=0) for sl in self.streamlines_xyz], axis=0)
        maxs = np.max([sl.max(axis=0) for sl in self.streamlines_xyz], axis=0)
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
                spline_subdiv=self.spline_subdiv,
                lookup_colormap=self.lut,
                tube_sides=20,
                # lod=False,  # <—
                # lod_points=20000,  # optional, tune
                # lod_points_size=2,  # optional, tune
            )
        else:
            self.actor0 = actor.line(
                self.streamlines_xyz,
                colors=self.flat_vals,  # scalars
                linewidth=self.linewidth,
                lookup_colormap=self.lut,
            )

        self.scene.add(self.actor0)
        self._style_streamline_actor()

        # fast actor for interaction (cheap line rendering)
        self.actor_fast = actor.line(
            self.streamlines_xyz,
            colors=self.flat_vals,  # pass scalars
            linewidth=1.0,  # very light
            lookup_colormap=self.lut,
            # lod=False,  # make it deterministic
            fake_tube=True,  # a bit of shading to hint tubes
        )
        self.actor_fast.SetVisibility(False)
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
        prop.SetAmbient(0.1)
        prop.SetDiffuse(0.95)
        prop.SetSpecular(0.25)
        prop.SetSpecularPower(12)

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
        self.scene.add(
            fury.actor.scalar_bar(lookup_table=self.lut, title=self.color_label)
        )

    def _render_now(self):
        try:
            self.scene.ResetCameraClippingRange()
        except Exception:
            pass
        if self.showm is not None:
            self.showm.render()

    def _sync_plane_from_widget(self, *_):
        origin = [0.0, 0.0, 0.0]
        normal = [1.0, 0.0, 0.0]
        self.plane_rep.GetOrigin(origin)
        self.plane_rep.GetNormal(normal)
        self.plane_fn.SetOrigin(origin)
        self.plane_fn.SetNormal(normal)

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
        if self.clipping_active:
            # Deactivate clipping
            self.actor0.GetMapper().RemoveAllClippingPlanes()
            self.actor_fast.GetMapper().RemoveAllClippingPlanes()
            self.clipping_active = False
            print("Clipping OFF")
        else:
            # Activate clipping
            self.actor0.GetMapper().AddClippingPlane(self.plane_fn)
            self.actor_fast.GetMapper().AddClippingPlane(self.plane_fn)
            self.clipping_active = True
            print("Clipping ON")

        self._render_now()

    def _rebuild_unclipped_actor(self):
        clipping_on = (
            self.mapper0.GetNumberOfClippingPlanes() > 0
            if self.mapper0 is not None
            else False
        )

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
                spline_subdiv=self.spline_subdiv,
                tube_sides=20,
                lookup_colormap=self.lut,
                # lod=False,  # <—
                # lod_points=20000,  # optional, tune
                # lod_points_size=2,  # optional, tune
            )
        else:
            self.actor0 = actor.line(
                self.streamlines_xyz,
                colors=self.flat_vals,
                linewidth=self.linewidth,
                lookup_colormap=self.lut,
            )

        self.scene.add(self.actor0)
        self._style_streamline_actor()

        self.mapper0 = self.actor0.GetMapper()
        if clipping_on:
            self.mapper0.RemoveAllClippingPlanes()
            self.mapper0.AddClippingPlane(self.plane_fn)

    # ---------------------------
    # key handling
    # ---------------------------
    def _on_keypress(self, obj, evt):
        key = obj.GetKeySym().lower()

        if key == "o":
            self._toggle_clipping()

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
            import datetime
            import os

            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            out_path = os.path.abspath(f"view_{ts}_hires.png")
            try:
                highres_size = (2000, 2000)
                fury.window.record(
                    scene=self.scene,
                    out_path=out_path,
                    size=highres_size,
                    reset_camera=False,
                )
                print(
                    f"Saved high-resolution screenshot to {out_path} ({highres_size[0]}x{highres_size[1]})"
                )
            except Exception as e:
                print(f"Failed to save screenshot: {e}")

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

            # Apply clipping only if it was previously enabled
            if self.clipping_active:
                self.actor_fast.GetMapper().AddClippingPlane(
                    self.plane_fn
                )  # Apply clipping to the fast actor
        except Exception:
            pass
        self._render_now()

    def _lod_off(self, obj=None, evt=None):
        """Switch back to full-res actor, applying clipping if it was previously enabled."""
        # Rebuild actor to ensure it shows at full resolution with clipping applied
        self._rebuild_unclipped_actor()

        try:
            self.actor_fast.SetVisibility(False)  # Hide the fast actor
            self.actor0.SetVisibility(True)  # Show the full-res actor

            # Apply clipping to the full-res actor if it was previously enabled
            if self.clipping_active:
                self.actor0.GetMapper().AddClippingPlane(
                    self.plane_fn
                )  # Reapply clipping
        except Exception:
            pass
        self._render_now()

    # ---------------------------
    # main entry
    # ---------------------------
    def run(self, interactive: bool, screenshot_path: str | None):
        if interactive:
            self.showm = window.ShowManager(
                scene=self.scene, size=self.window_size, reset_camera=False
            )
            self.showm.initialize()

            iren = self.showm.iren
            iren.SetDesiredUpdateRate(60.0)
            iren.SetStillUpdateRate(30.0)  # higher so it returns to full-res

            # ensure anti-alias looks good when idle
            try:
                self.showm.renwin.SetMultiSamples(0)
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

            self.showm.iren.AddObserver("KeyPressEvent", self._on_keypress)

            self.scene.reset_camera()

            self.scene.azimuth(15)
            self.scene.elevation(10)
            self.scene.zoom(1.1)

            print(
                "Keys: O toggle plane, H hide gizmo, I flip side, R reset plane, +/- thickness, B background, S scale bar, P save PNG"
            )
            self.showm.start()
        else:
            if not screenshot_path:
                raise ValueError("Must specify screenshot_path when interactive=False.")
            self.scene.reset_camera()
            fury.window.record(
                scene=self.scene, out_path=screenshot_path, size=self.window_size
            )


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
    window_size: tuple[int, int] = (800, 800),
    downsample_factor: int = 2,
    max_streamlines: int | None = None,
    filter_min_len: int | None = None,
    subsample_factor: int = 1,
    crop_bounds: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    | None = None,
    colormap=None,
    background_color: str | tuple[float, float, float] = "black",
    spline_subdiv: int = 16,
    color_range: tuple[float, float] | None = None,
    color_label: str = "Angle (deg)",
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

    if subsample_factor > 1:
        print(f"Subsampling: keeping 1 in every {subsample_factor} streamlines")
        total = len(streamlines_xyz)
        keep_idx = sorted(
            random.sample(range(total), max(1, total // subsample_factor))
        )
        streamlines_xyz = [streamlines_xyz[i] for i in keep_idx]
        color_values = [color_values[i] for i in keep_idx]

    if max_streamlines is not None and len(streamlines_xyz) > max_streamlines:
        print(f"Limiting to max {max_streamlines} streamlines")
        keep_idx = sorted(random.sample(range(len(streamlines_xyz)), max_streamlines))
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
        color_range=lut_range,
        color_label=color_label,
    )
    viewer.run(interactive=interactive, screenshot_path=screenshot_path)


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
