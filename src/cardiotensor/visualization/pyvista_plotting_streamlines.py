#!/usr/bin/env python3
from __future__ import annotations

import datetime
import random
from pathlib import Path

import numpy as np

_NAMED_COLORS = {
    "white": (1.0, 1.0, 1.0),
    "black": (0.0, 0.0, 0.0),
    "gray": (0.5, 0.5, 0.5),
    "grey": (0.5, 0.5, 0.5),
    "lightgray": (0.9, 0.9, 0.9),
    "lightgrey": (0.9, 0.9, 0.9),
    "charcoal": (0.05, 0.06, 0.08),
    "midnight": (0.02, 0.03, 0.07),
}


def _load_pyvista():
    try:
        import pyvista as pv
    except ImportError as err:
        raise ImportError(
            "PyVista backend requested but pyvista is not installed. "
            "Install the project with the updated dependencies, or run `pip install pyvista`."
        ) from err
    return pv


def _parse_background_color(color) -> tuple[float, float, float] | str:
    if color is None:
        return _NAMED_COLORS["white"]
    if isinstance(color, str):
        return _NAMED_COLORS.get(color.lower(), color)
    if isinstance(color, tuple | list) and len(color) == 3:
        return tuple(float(c) for c in color)
    raise TypeError("background_color must be a string or 3-tuple")


def _is_dark_color(color: tuple[float, float, float] | str) -> bool:
    if isinstance(color, str):
        color = _NAMED_COLORS.get(color.lower(), (0.0, 0.0, 0.0))
    r, g, b = color
    return 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.45


def _downsample_streamline(streamline: np.ndarray, factor: int = 2) -> np.ndarray:
    return streamline if len(streamline) < 3 or factor <= 1 else streamline[::factor]


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
        if len(seg) > 1:
            segs.append(seg)
            cols.append(col)
    return segs, cols


def _prepare_streamlines(
    streamlines_xyz: list[np.ndarray],
    color_values: list[np.ndarray],
    downsample_factor: int,
    max_streamlines: int | None,
    filter_min_len: int | None,
    subsample_factor: int,
    crop_bounds: (
        tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None
    ),
    random_seed: int | None = None,
) -> tuple[list[np.ndarray], list[np.ndarray], tuple[float, float]]:
    if len(streamlines_xyz) != len(color_values):
        raise ValueError("color_values must contain one array per streamline")

    print(f"Initial number of streamlines: {len(streamlines_xyz)}")
    full_mins, full_maxs = _compute_streamline_bounds(streamlines_xyz)
    _print_box_shape("Full streamline box", full_mins, full_maxs)

    if crop_bounds is not None:
        (x_min, x_max), (y_min, y_max), (z_min, z_max) = crop_bounds
        print(f"Cropping streamlines within bounds: {crop_bounds}")
        new_streamlines, new_colors = [], []
        for sl, cl in zip(streamlines_xyz, color_values):
            segs, cols = _split_streamline_by_bounds(
                np.asarray(sl), np.asarray(cl), x_min, x_max, y_min, y_max, z_min, z_max
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
        sl_arr = np.asarray(sl, dtype=np.float32)
        cl_arr = np.asarray(cl, dtype=np.float32).ravel()
        if len(sl_arr) != len(cl_arr):
            raise ValueError("Each color array must match its streamline length")
        ds_sl = _downsample_streamline(sl_arr, downsample_factor)
        ds_cl = _downsample_streamline(cl_arr, downsample_factor)
        if len(ds_sl) > 1 and (filter_min_len is None or len(ds_sl) >= filter_min_len):
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

    flat_colors = np.concatenate([np.asarray(c).ravel() for c in color_values]).astype(
        np.float32
    )
    min_val = float(np.nanmin(flat_colors))
    max_val = float(np.nanmax(flat_colors))
    print(f"Final number of streamlines to render: {len(streamlines_xyz)}")
    print(f"Coloring range: min={min_val:.3f}, max={max_val:.3f}")
    return streamlines_xyz, color_values, (min_val, max_val)


def _crop_streamlines_by_bounds(
    streamlines_xyz: list[np.ndarray],
    color_values: list[np.ndarray],
    bounds: tuple[float, float, float, float, float, float],
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    cropped_streamlines, cropped_colors = [], []
    for sl, cl in zip(streamlines_xyz, color_values):
        segs, cols = _split_streamline_by_bounds(
            np.asarray(sl), np.asarray(cl), x_min, x_max, y_min, y_max, z_min, z_max
        )
        cropped_streamlines.extend(segs)
        cropped_colors.extend(cols)
    return cropped_streamlines, cropped_colors


def _format_crop_cli_bounds(
    bounds: tuple[float, float, float, float, float, float],
) -> str:
    return (
        f"--crop-x {bounds[0]:.2f} {bounds[1]:.2f} "
        f"--crop-y {bounds[2]:.2f} {bounds[3]:.2f} "
        f"--crop-z {bounds[4]:.2f} {bounds[5]:.2f}"
    )


def _smooth_polydata_with_vtk_spline(
    poly,
    pv,
    subdivisions: int,
):
    if subdivisions <= 1:
        return poly
    try:
        from vtkmodules.vtkFiltersGeneral import vtkSplineFilter
    except ImportError:
        try:
            from vtk import vtkSplineFilter
        except ImportError:
            print("Visual smoothing unavailable: VTK spline filter is not installed.")
            return poly

    spline_filter = vtkSplineFilter()
    spline_filter.SetInputData(poly)
    spline_filter.SetSubdivideToSpecified()
    spline_filter.SetNumberOfSubdivisions(int(subdivisions))
    spline_filter.Update()
    return pv.wrap(spline_filter.GetOutput()).copy(deep=True)


def _print_render_summary(
    streamlines_xyz: list[np.ndarray],
    color_values: list[np.ndarray],
    data_range: tuple[float, float],
    color_range: tuple[float, float] | None,
    mode: str,
    line_width: float,
    tube_sides: int,
    downsample_factor: int,
    subsample_factor: int,
    max_streamlines: int | None,
    filter_min_len: int | None,
    crop_bounds,
    spline_subdiv: int,
) -> None:
    lengths = np.array([len(sl) for sl in streamlines_xyz], dtype=np.int64)
    total_points = int(lengths.sum())
    total_segments = int(np.maximum(lengths - 1, 0).sum())
    mins, maxs = _compute_streamline_bounds(streamlines_xyz)
    box_size = maxs - mins
    diagonal = float(np.linalg.norm(box_size))
    scalar_points = int(sum(len(np.asarray(c).ravel()) for c in color_values))
    effective_range = color_range or data_range
    approx_tube_polys = total_segments * max(3, int(tube_sides)) * 2

    print("\nPyVista render summary")
    print("-" * 24)
    print(f"Backend/mode: pyvista / {mode}")
    print(f"Streamlines: {len(streamlines_xyz):,}")
    print(
        "Points per streamline: "
        f"mean={float(lengths.mean()):.1f}, "
        f"median={float(np.median(lengths)):.1f}, "
        f"min={int(lengths.min())}, max={int(lengths.max())}"
    )
    print(f"Total points: {total_points:,}")
    print(f"Total segments: {total_segments:,}")
    print(f"Scalar samples: {scalar_points:,}")
    print(
        "Bounds: "
        f"x=[{mins[0]:.2f}, {maxs[0]:.2f}], "
        f"y=[{mins[1]:.2f}, {maxs[1]:.2f}], "
        f"z=[{mins[2]:.2f}, {maxs[2]:.2f}]"
    )
    print(
        "Box size: "
        f"dx={box_size[0]:.2f}, dy={box_size[1]:.2f}, dz={box_size[2]:.2f}, "
        f"diagonal={diagonal:.2f}"
    )
    print(f"Data color range: min={data_range[0]:.3f}, max={data_range[1]:.3f}")
    print(
        f"Colorbar range: min={effective_range[0]:.3f}, "
        f"max={effective_range[1]:.3f}"
    )
    print(
        f"Requested crop bounds: {crop_bounds if crop_bounds is not None else 'none'}"
    )
    print(
        "Filtering: "
        f"downsample_factor={downsample_factor}, "
        f"subsample_factor={subsample_factor}, "
        f"max_streamlines={max_streamlines}, "
        f"min_length={filter_min_len}, "
        f"spline_subdiv={spline_subdiv}"
    )
    print(
        "Geometry: "
        f"line_width/radius={line_width}, tube_sides={tube_sides}, "
        f"approx_tube_triangles={approx_tube_polys:,}"
    )
    if mode == "tube" and approx_tube_polys > 1_000_000:
        print(
            "Performance hint: many streamlines. Press M to switch to line mode for faster navigation."
        )
    elif mode == "line":
        print(
            "Performance hint: line mode is best for navigation; tube mode is best for final figures."
        )
    print("-" * 24)


def _streamlines_to_polydata(
    streamlines_xyz: list[np.ndarray],
    color_values: list[np.ndarray],
    scalar_name: str,
):
    pv = _load_pyvista()
    lengths = np.fromiter(
        (len(streamline) for streamline in streamlines_xyz),
        dtype=np.int64,
        count=len(streamlines_xyz),
    )
    n_points = int(lengths.sum())
    points = np.concatenate(streamlines_xyz, axis=0).astype(np.float32, copy=False)
    scalars = np.concatenate(color_values).astype(np.float32, copy=False)

    # VTK stores lines as [length, point ids..., length, point ids...]. Build
    # that array in NumPy instead of appending millions of Python integers.
    line_starts = np.empty(len(lengths), dtype=np.int64)
    line_starts[0] = 0
    if len(lengths) > 1:
        np.cumsum(lengths[:-1] + 1, out=line_starts[1:])
    lines = np.empty(n_points + len(lengths), dtype=np.int64)
    lines[line_starts] = lengths
    point_slots = np.ones(len(lines), dtype=bool)
    point_slots[line_starts] = False
    lines[point_slots] = np.arange(n_points, dtype=np.int64)

    poly = pv.PolyData(points)
    poly.lines = lines
    poly.point_data[scalar_name] = scalars
    poly.set_active_scalars(scalar_name)
    return poly


def show_streamlines_pyvista(
    streamlines_xyz: list[np.ndarray],
    color_values: list[np.ndarray],
    mode: str = "tube",
    line_width: float = 2.0,
    interactive: bool = True,
    screenshot_path: str | None = None,
    video_path: str | None = None,
    video_fps: int = 30,
    video_frames: int = 120,
    window_size: tuple[int, int] = (1000, 900),
    downsample_factor: int = 2,
    max_streamlines: int | None = None,
    filter_min_len: int | None = None,
    subsample_factor: int = 1,
    crop_bounds: (
        tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None
    ) = None,
    colormap=None,
    background_color: str | tuple[float, float, float] | None = None,
    tube_sides: int = 12,
    color_range: tuple[float, float] | None = None,
    color_label: str = "Angle (deg)",
    opacity: float = 1.0,
    show_axes: bool = True,
    show_bounds: bool = False,
    shadows: bool = False,
    spline_subdiv: int = 1,
    random_seed: int | None = None,
) -> None:
    pv = _load_pyvista()
    mode = mode.lower().strip()
    if mode not in {"tube", "line"}:
        raise ValueError("mode must be 'tube' or 'line'")

    streamlines_xyz, color_values, data_range = _prepare_streamlines(
        streamlines_xyz=streamlines_xyz,
        color_values=color_values,
        downsample_factor=downsample_factor,
        max_streamlines=max_streamlines,
        filter_min_len=filter_min_len,
        subsample_factor=subsample_factor,
        crop_bounds=crop_bounds,
        random_seed=random_seed,
    )
    spline_subdiv = max(1, int(spline_subdiv))
    if spline_subdiv <= 1:
        print("Visual smoothing: off")
    clim = color_range or data_range
    print(f"Rendering backend: pyvista, mode: {mode}")
    print(f"Colorbar range: min={clim[0]:.3f}, max={clim[1]:.3f}")
    _print_render_summary(
        streamlines_xyz=streamlines_xyz,
        color_values=color_values,
        data_range=data_range,
        color_range=color_range,
        mode=mode,
        line_width=line_width,
        tube_sides=tube_sides,
        downsample_factor=downsample_factor,
        subsample_factor=subsample_factor,
        max_streamlines=max_streamlines,
        filter_min_len=filter_min_len,
        crop_bounds=crop_bounds,
        spline_subdiv=spline_subdiv,
    )

    base_streamlines_xyz = [np.asarray(sl, dtype=np.float32) for sl in streamlines_xyz]
    base_color_values = [
        np.asarray(cl, dtype=np.float32).ravel() for cl in color_values
    ]

    scalar_name = color_label or "Scalar"

    def _make_display_polydata(
        display_streamlines: list[np.ndarray], display_colors: list[np.ndarray]
    ):
        display_poly = _streamlines_to_polydata(
            display_streamlines, display_colors, scalar_name
        )
        if spline_subdiv <= 1:
            return display_poly

        before_points = int(display_poly.n_points)
        display_poly = _smooth_polydata_with_vtk_spline(
            display_poly,
            pv,
            subdivisions=spline_subdiv,
        )
        display_poly.set_active_scalars(scalar_name)
        after_points = int(display_poly.n_points)
        display_segments = max(0, after_points - len(display_streamlines))
        approx_tube_triangles = display_segments * max(3, int(tube_sides)) * 2
        print(
            "Visual smoothing: "
            "VTK vtkSplineFilter/cardinal spline, "
            f"{spline_subdiv} segments per streamline "
            f"({before_points:,} -> {after_points:,} points)"
        )
        if mode == "tube":
            print(
                "Display geometry after smoothing: "
                f"approx_tube_triangles={approx_tube_triangles:,}"
            )
        return display_poly

    poly = _make_display_polydata(streamlines_xyz, color_values)

    bg = _parse_background_color(background_color)
    is_dark = _is_dark_color(bg)
    text_color = "white" if is_dark else "black"
    bounds_color = (0.75, 0.78, 0.82) if is_dark else (0.25, 0.28, 0.32)

    plotter = pv.Plotter(window_size=window_size, off_screen=not interactive)
    if is_dark:
        plotter.set_background(bg, top=(0.12, 0.14, 0.18))
    else:
        plotter.set_background(bg, top=(0.96, 0.97, 0.99))

    # Headlight: attached to the camera, always illuminates from the viewer's
    # direction — identical behaviour to FURY's default lighting.  Even brightness
    # across the whole sample regardless of viewing angle.
    try:
        plotter.clear_lights()
        plotter.add_light(pv.Light(light_type="headlight", intensity=0.85))
    except Exception:
        pass

    scalar_bar_args = {
        "title": color_label,
        "vertical": True,
        "title_font_size": 16,
        "label_font_size": 12,
        "color": text_color,
        "position_x": 0.88,
        "position_y": 0.22,
        "width": 0.08,
        "height": 0.55,
    }

    state = {
        "axes": False,
        "bounds": False,
        "box": False,
        "box_actor": None,
        "background_dark": is_dark,
        "shadows": False,
        "streamline_actor": None,
        "actors": {"tube": None, "line": None},
        "render_mode": mode,
        "tube_thickness": max(float(line_width), 0.5),
        "line_thickness": max(float(line_width) * 0.4, 0.5),  # thinner LOD lines
        "crop_widget": None,
        "crop_widget_active": False,
        "crop_bounds": None,
        "auto_lod": interactive and mode == "tube",
        "lod_active": False,
        # lighting params — kept in state so sliders can read/write them
        "ambient": 0.60,
        "diffuse": 0.45,
        "specular": 0.08,
        "light_intensity": 0.85,
        "light_sliders_visible": False,
    }
    tube_sides = max(3, int(tube_sides))

    bounds = poly.bounds
    full_bounds = tuple(float(v) for v in bounds)

    def _enable_shadows() -> None:
        try:
            plotter.enable_shadows()
            state["shadows"] = True
            print("Shadows ON")
        except Exception as err:
            print(f"Shadows unavailable in this PyVista/VTK setup: {err}")

    def _disable_shadows() -> None:
        try:
            plotter.disable_shadows()
        except Exception:
            pass
        state["shadows"] = False
        print("Shadows OFF")

    def _toggle_shadows() -> None:
        if state["shadows"]:
            _disable_shadows()
        else:
            _enable_shadows()
        _render()

    if shadows:
        _enable_shadows()

    def _add_streamline_actor(render_mode: str, show_scalar_bar: bool) -> None:
        # Both modes use the same line PolyData with render_lines_as_tubes — no
        # poly.tube() call, so startup and thickness changes are instant.
        # Tube mode enables lighting so the GPU shader computes per-pixel tube normals:
        # silhouette edges are dark (→ natural black outline), lit faces are bright.
        lw = (
            state["tube_thickness"]
            if render_mode == "tube"
            else state["line_thickness"]
        )
        mesh_kwargs = {
            "scalars": scalar_name,
            "cmap": colormap,
            "clim": clim,
            "opacity": opacity,
            "name": f"streamlines_{render_mode}",
            "line_width": lw,
            "render_lines_as_tubes": True,
        }
        if render_mode == "tube":
            # Read current lighting values from state so interactive sliders take effect.
            mesh_kwargs.update(
                {
                    "lighting": True,
                    "ambient": state["ambient"],
                    "diffuse": state["diffuse"],
                    "specular": state["specular"],
                    "specular_power": 10,
                }
            )
        else:
            # Line mode: flat colours, no shading, used for LOD during camera motion.
            mesh_kwargs["lighting"] = False

        if show_scalar_bar:
            mesh_kwargs["scalar_bar_args"] = scalar_bar_args
        else:
            mesh_kwargs["show_scalar_bar"] = False

        actor = plotter.add_mesh(poly, reset_camera=False, **mesh_kwargs)
        try:
            actor.GetProperty().SetLineWidth(lw)
        except Exception:
            pass
        actor.SetVisibility(render_mode == state["render_mode"])
        state["actors"][render_mode] = actor
        if render_mode == state["render_mode"]:
            state["streamline_actor"] = actor

    def _ensure_streamline_actor(
        render_mode: str, show_scalar_bar: bool = False
    ) -> None:
        if state["actors"][render_mode] is None:
            _add_streamline_actor(render_mode, show_scalar_bar=show_scalar_bar)
        else:
            lw = (
                state["tube_thickness"]
                if render_mode == "tube"
                else state["line_thickness"]
            )
            try:
                state["actors"][render_mode].GetProperty().SetLineWidth(lw)
            except Exception:
                pass

    def _show_streamline_mode(render_mode: str) -> None:
        _ensure_streamline_actor(render_mode)
        for mode_name, actor in state["actors"].items():
            if actor is not None:
                actor.SetVisibility(mode_name == render_mode)
        state["streamline_actor"] = state["actors"][render_mode]
        state["lod_active"] = False

    def _show_interaction_lod(*args) -> None:
        if not state["auto_lod"] or state["render_mode"] != "tube":
            return
        _ensure_streamline_actor("line")
        tube_actor = state["actors"].get("tube")
        line_actor = state["actors"].get("line")
        if tube_actor is None or line_actor is None:
            return
        tube_actor.SetVisibility(False)
        line_actor.SetVisibility(True)
        state["streamline_actor"] = line_actor
        state["lod_active"] = True
        _render()

    def _hide_interaction_lod(*args) -> None:
        if not state["lod_active"]:
            return
        tube_actor = state["actors"].get("tube")
        line_actor = state["actors"].get("line")
        if tube_actor is not None:
            tube_actor.SetVisibility(True)
            state["streamline_actor"] = tube_actor
        if line_actor is not None:
            line_actor.SetVisibility(False)
        state["lod_active"] = False
        _render()

    def _add_interaction_observer(event_name: str, callback) -> bool:
        iren = getattr(plotter, "iren", None)
        candidates = [iren, getattr(iren, "interactor", None)]
        for candidate in candidates:
            if candidate is None:
                continue
            add_observer = getattr(candidate, "add_observer", None)
            if add_observer is not None:
                add_observer(event_name, callback)
                return True
            add_observer = getattr(candidate, "AddObserver", None)
            if add_observer is not None:
                add_observer(event_name, callback)
                return True
        return False

    _add_streamline_actor(state["render_mode"], show_scalar_bar=True)
    if state["render_mode"] == "tube":
        _add_streamline_actor("line", show_scalar_bar=False)

    def _make_outline_mesh():
        xmin, xmax, ymin, ymax, zmin, zmax = poly.bounds
        points = np.array(
            [
                [xmin, ymin, zmin],
                [xmax, ymin, zmin],
                [xmax, ymax, zmin],
                [xmin, ymax, zmin],
                [xmin, ymin, zmax],
                [xmax, ymin, zmax],
                [xmax, ymax, zmax],
                [xmin, ymax, zmax],
            ],
            dtype=np.float32,
        )
        edges = np.array(
            [
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 0],
                [4, 5],
                [5, 6],
                [6, 7],
                [7, 4],
                [0, 4],
                [1, 5],
                [2, 6],
                [3, 7],
            ],
            dtype=np.int64,
        )
        lines = np.column_stack([np.full(len(edges), 2), edges]).ravel()
        outline = pv.PolyData(points)
        outline.lines = lines
        return outline

    outline_mesh = _make_outline_mesh()

    def _show_bounds() -> None:
        try:
            plotter.show_bounds(
                grid="front",
                location="outer",
                all_edges=False,
                color=bounds_color,
                font_size=10,
                xtitle="X",
                ytitle="Y",
                ztitle="Z",
            )
        except TypeError:
            plotter.show_bounds(
                grid="front",
                location="outer",
                all_edges=False,
                color=bounds_color,
                font_size=10,
                xlabel="X",
                ylabel="Y",
                zlabel="Z",
            )
        state["bounds"] = True

    def _hide_bounds() -> None:
        try:
            plotter.remove_bounds_axes()
        except Exception:
            try:
                plotter.remove_bounds_axis()
            except Exception:
                pass
        state["bounds"] = False

    def _show_axes() -> None:
        try:
            plotter.add_axes(line_width=2, color=text_color)
        except TypeError:
            plotter.add_axes()
        state["axes"] = True

    def _hide_axes() -> None:
        try:
            plotter.hide_axes()
        except Exception:
            pass
        state["axes"] = False

    def _render() -> None:
        try:
            plotter.render()
        except Exception:
            pass

    def _toggle_background() -> None:
        if state["background_dark"]:
            plotter.set_background("white", top=(0.96, 0.97, 0.99))
            state["background_dark"] = False
            print("Background set to white")
        else:
            plotter.set_background("black", top=(0.12, 0.14, 0.18))
            state["background_dark"] = True
            print("Background set to black")
        _render()

    def _toggle_bounds() -> None:
        if state["bounds"]:
            _hide_bounds()
            print("Bounds grid OFF")
        else:
            _show_bounds()
            print("Bounds grid ON")
        _render()

    def _toggle_axes() -> None:
        if state["axes"]:
            _hide_axes()
            print("Axes OFF")
        else:
            _show_axes()
            print("Axes ON")
        _render()

    def _add_box_actor():
        return plotter.add_mesh(
            outline_mesh,
            color=bounds_color,
            line_width=3,
            render_lines_as_tubes=True,
            name="streamline_bounding_box",
            reset_camera=False,
        )

    def _toggle_box() -> None:
        if state["box"]:
            try:
                plotter.remove_actor(state["box_actor"])
            except Exception:
                pass
            state["box_actor"] = None
            state["box"] = False
            print("Bounding box OFF")
        else:
            state["box_actor"] = _add_box_actor()
            state["box"] = True
            print("Bounding box ON")
        _render()

    def _capture_camera_state():
        camera = plotter.camera
        try:
            return {
                "position": camera.GetPosition(),
                "focal_point": camera.GetFocalPoint(),
                "view_up": camera.GetViewUp(),
                "parallel_scale": camera.GetParallelScale(),
                "view_angle": camera.GetViewAngle(),
                "parallel_projection": camera.GetParallelProjection(),
            }
        except Exception:
            return None

    def _restore_camera_state(camera_state) -> None:
        if camera_state is None:
            return
        camera = plotter.camera
        try:
            camera.SetPosition(camera_state["position"])
            camera.SetFocalPoint(camera_state["focal_point"])
            camera.SetViewUp(camera_state["view_up"])
            camera.SetParallelScale(camera_state["parallel_scale"])
            camera.SetViewAngle(camera_state["view_angle"])
            if camera_state["parallel_projection"]:
                camera.ParallelProjectionOn()
            else:
                camera.ParallelProjectionOff()
            plotter.reset_camera_clipping_range()
        except Exception:
            pass

    def _set_displayed_streamlines(
        new_streamlines: list[np.ndarray], new_colors: list[np.ndarray]
    ) -> None:
        nonlocal poly, outline_mesh
        camera_state = _capture_camera_state()

        for actor in state["actors"].values():
            if actor is not None:
                try:
                    plotter.remove_actor(actor)
                except Exception:
                    pass
        state["actors"] = {"tube": None, "line": None}
        state["streamline_actor"] = None

        poly = _make_display_polydata(new_streamlines, new_colors)
        outline_mesh = _make_outline_mesh()

        if state["box"]:
            try:
                plotter.remove_actor(state["box_actor"])
            except Exception:
                pass
            state["box_actor"] = _add_box_actor()

        if state["bounds"]:
            _hide_bounds()
            _show_bounds()

        _ensure_streamline_actor(state["render_mode"], show_scalar_bar=True)
        if state["auto_lod"]:
            _ensure_streamline_actor("line")
        _show_streamline_mode(state["render_mode"])
        _restore_camera_state(camera_state)
        _render()

    def _print_interactive_crop_summary(
        bounds: tuple[float, float, float, float, float, float],
        cropped_streamlines: list[np.ndarray],
    ) -> None:
        lengths = np.array([len(sl) for sl in cropped_streamlines], dtype=np.int64)
        mins, maxs = _compute_streamline_bounds(cropped_streamlines)
        box_size = maxs - mins
        print("\nInteractive crop applied")
        print("-" * 24)
        print(f"Streamlines: {len(cropped_streamlines):,}")
        print(
            "Points per streamline: "
            f"mean={float(lengths.mean()):.1f}, "
            f"median={float(np.median(lengths)):.1f}, "
            f"min={int(lengths.min())}, max={int(lengths.max())}"
        )
        print(
            "Crop box: "
            f"x=[{bounds[0]:.2f}, {bounds[1]:.2f}], "
            f"y=[{bounds[2]:.2f}, {bounds[3]:.2f}], "
            f"z=[{bounds[4]:.2f}, {bounds[5]:.2f}]"
        )
        print(
            "Displayed bounds: "
            f"x=[{mins[0]:.2f}, {maxs[0]:.2f}], "
            f"y=[{mins[1]:.2f}, {maxs[1]:.2f}], "
            f"z=[{mins[2]:.2f}, {maxs[2]:.2f}]"
        )
        print(
            "Displayed size: "
            f"dx={box_size[0]:.2f}, dy={box_size[1]:.2f}, dz={box_size[2]:.2f}"
        )
        print(f"Reproduce from CLI: {_format_crop_cli_bounds(bounds)}")
        print("-" * 24)

    def _apply_interactive_crop(box) -> None:
        try:
            selected_bounds = tuple(float(v) for v in box.bounds)
        except Exception as err:
            print(f"Interactive crop failed to read box bounds: {err}")
            return

        cropped_streamlines, cropped_colors = _crop_streamlines_by_bounds(
            base_streamlines_xyz, base_color_values, selected_bounds
        )
        if not cropped_streamlines:
            print("Interactive crop ignored: no streamlines inside the selected box.")
            return

        state["crop_bounds"] = selected_bounds
        _set_displayed_streamlines(cropped_streamlines, cropped_colors)
        _print_interactive_crop_summary(selected_bounds, cropped_streamlines)

    def _disable_crop_widget(quiet: bool = False) -> None:
        try:
            plotter.clear_box_widgets()
        except Exception:
            widget = state.get("crop_widget")
            if widget is not None:
                try:
                    widget.Off()
                except Exception:
                    pass
        state["crop_widget"] = None
        state["crop_widget_active"] = False
        if not quiet:
            print(
                "Interactive crop box OFF; current crop is kept. Press X to clear it."
            )
        _render()

    def _enable_crop_widget() -> None:
        if state["crop_widget_active"]:
            return
        widget_bounds = state["crop_bounds"] or full_bounds
        try:
            state["crop_widget"] = plotter.add_box_widget(
                callback=_apply_interactive_crop,
                bounds=widget_bounds,
                factor=1.0,
                rotation_enabled=False,
                color=bounds_color,
                use_planes=False,
                outline_translation=True,
                interaction_event="end",
            )
        except TypeError:
            state["crop_widget"] = plotter.add_box_widget(
                callback=_apply_interactive_crop,
                bounds=widget_bounds,
                factor=1.0,
                rotation_enabled=False,
                color=bounds_color,
                use_planes=False,
                outline_translation=True,
            )
        except Exception as err:
            print(f"Interactive crop box unavailable: {err}")
            return
        state["crop_widget_active"] = True
        print(
            "Interactive crop box ON: drag the box faces, crop updates when you release. "
            "Press X to clear crop."
        )
        _render()

    def _toggle_interactive_crop() -> None:
        if state["crop_widget_active"]:
            _disable_crop_widget()
        else:
            _enable_crop_widget()

    def _clear_interactive_crop() -> None:
        was_active = state["crop_widget_active"]
        if was_active:
            _disable_crop_widget(quiet=True)
        state["crop_bounds"] = None
        _set_displayed_streamlines(base_streamlines_xyz, base_color_values)
        print(
            "Interactive crop cleared; showing all streamlines after command-line filters."
        )
        if was_active:
            _enable_crop_widget()

    def _redraw_current_streamline_actor() -> None:
        render_mode = state["render_mode"]
        camera_state = _capture_camera_state()
        actor = state["actors"].get(render_mode)
        if actor is not None:
            try:
                plotter.remove_actor(actor)
            except Exception:
                pass
            state["actors"][render_mode] = None
        _ensure_streamline_actor(render_mode)
        _show_streamline_mode(render_mode)
        _restore_camera_state(camera_state)

    def _change_thickness(factor: float) -> None:
        # Fast path: just update line_width on the existing actor — no geometry rebuild.
        if state["render_mode"] == "tube":
            lower, upper = 0.5, 200.0
            state["tube_thickness"] = min(
                upper, max(lower, state["tube_thickness"] * factor)
            )
            lw = state["tube_thickness"]
            print(f"Thickness: {lw:.1f} px")
        else:
            lower, upper = 0.5, 100.0
            state["line_thickness"] = min(
                upper, max(lower, state["line_thickness"] * factor)
            )
            lw = state["line_thickness"]
            print(f"Thickness: {lw:.1f} px")

        # Update all relevant actors in-place (no rebuild).
        for mode_name, act in state["actors"].items():
            if act is None:
                continue
            expected_lw = (
                state["tube_thickness"]
                if mode_name == "tube"
                else state["line_thickness"]
            )
            try:
                act.GetProperty().SetLineWidth(expected_lw)
            except Exception:
                pass
        _render()

    def _increase_thickness() -> None:
        _change_thickness(1.25)

    def _decrease_thickness() -> None:
        _change_thickness(0.8)

    def _toggle_fast_mode() -> None:
        # Toggle between tube mode (thick, Phong-shaded, dark outlines) and
        # line mode (thin, flat colours, fastest for navigation).
        camera_state = _capture_camera_state()
        state["render_mode"] = "line" if state["render_mode"] == "tube" else "tube"
        state["auto_lod"] = interactive and state["render_mode"] == "tube"
        _show_streamline_mode(state["render_mode"])
        _restore_camera_state(camera_state)
        _render()
        if state["render_mode"] == "line":
            print(
                f"Line mode ON  (flat colours, width={state['line_thickness']:.1f} px)"
            )
        else:
            print(
                f"Tube mode ON  (Phong shading, width={state['tube_thickness']:.1f} px)"
            )

    def _reset_view() -> None:
        plotter.view_isometric()
        plotter.camera.zoom(1.15)
        _render()
        print("View reset")

    def _save_screenshot() -> None:
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        out_path = Path.cwd() / f"pyvista_view_{ts}.png"
        plotter.screenshot(str(out_path))
        print(f"Saved screenshot to {out_path}")

    def _record_orbit(out_path: str | Path, n_frames: int, fps: int) -> None:
        """Orbit the camera once and write every frame.

        Uses plotter.open_movie() / open_gif() + manual write_frame() loop,
        following the PyVista MP4 Movie example pattern:
          open_movie → (show if off-screen) → write_frame × N → close_movie
        This avoids the orbit_on_path(write_frames=True) internal bug.
        .mp4 requires imageio-ffmpeg:  pip install imageio-ffmpeg
        """
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        camera_state = _capture_camera_state()

        viewup = [0, 0, 1]
        focal_point = (
            0.5 * (full_bounds[0] + full_bounds[1]),
            0.5 * (full_bounds[2] + full_bounds[3]),
            0.5 * (full_bounds[4] + full_bounds[5]),
        )

        try:
            orbit_positions = plotter.generate_orbital_path(
                n_points=n_frames, viewup=viewup, factor=2.0
            ).points
        except Exception as err:
            print(f"Could not generate orbital path: {err}")
            return

        suffix = out_path.suffix.lower()
        try:
            if suffix == ".gif":
                plotter.open_gif(str(out_path))
            else:
                plotter.open_movie(str(out_path), framerate=fps)
        except Exception as err:
            print(
                f"Could not open writer for {out_path}: {err}\n"
                "Tip: for .mp4 install imageio-ffmpeg:  pip install imageio-ffmpeg"
            )
            _restore_camera_state(camera_state)
            return

        print(f"Recording {n_frames} frames at {fps} fps → {out_path} …")
        try:
            from tqdm import tqdm

            frame_iter = tqdm(
                orbit_positions, total=n_frames, unit="frame", desc="Encoding"
            )
        except ImportError:
            frame_iter = orbit_positions  # plain iterator, no progress bar

        try:
            for i, pos in enumerate(frame_iter):
                plotter.camera.position = tuple(float(v) for v in pos)
                plotter.camera.focal_point = focal_point
                plotter.camera.up = viewup
                plotter.reset_camera_clipping_range()
                plotter.render()
                plotter.write_frame()
        except Exception as err:
            print(f"Error while recording frame {i}: {err}")
        finally:
            try:
                plotter.close_movie()
            except Exception:
                pass

        _restore_camera_state(camera_state)
        _render()
        print(f"Saved orbit video to {out_path}")

    def _record_orbit_interactive() -> None:
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        out_path = Path.cwd() / f"pyvista_orbit_{ts}.mp4"
        _record_orbit(out_path, n_frames=video_frames, fps=video_fps)

    # ------------------------------------------------------------------ lighting
    def _apply_lighting_to_actors() -> None:
        """Push current state ambient/diffuse/specular onto every tube actor."""
        for act in state["actors"].values():
            if act is None:
                continue
            try:
                prop = act.GetProperty()
                prop.SetAmbient(state["ambient"])
                prop.SetDiffuse(state["diffuse"])
                prop.SetSpecular(state["specular"])
            except Exception:
                pass

    def _apply_light_intensity(value: float) -> None:
        """Set intensity on all lights in the renderer."""
        try:
            lights = plotter.renderer.GetLights()
            lights.InitTraversal()
            light = lights.GetNextItem()
            while light:
                light.SetIntensity(value)
                light = lights.GetNextItem()
        except Exception:
            pass

    def _toggle_light_sliders() -> None:
        if state["light_sliders_visible"]:
            # Hide: remove all slider widgets (crop widget is not a slider, safe to clear)
            try:
                plotter.clear_slider_widgets()
            except Exception:
                pass
            state["light_sliders_visible"] = False
            print("Lighting sliders OFF")
            _render()
            return

        # Show three sliders on the left edge of the window.
        slider_style = {
            "slider_width": 0.03,
            "tube_width": 0.005,
        }

        def _cb_ambient(value: float) -> None:
            state["ambient"] = round(value, 3)
            _apply_lighting_to_actors()
            _render()

        def _cb_diffuse(value: float) -> None:
            state["diffuse"] = round(value, 3)
            _apply_lighting_to_actors()
            _render()

        def _cb_light(value: float) -> None:
            state["light_intensity"] = round(value, 3)
            _apply_light_intensity(value)
            _render()

        try:
            plotter.add_slider_widget(
                _cb_ambient,
                rng=[0.0, 1.0],
                value=state["ambient"],
                title="Ambient",
                pointa=(0.02, 0.74),
                pointb=(0.22, 0.74),
                color=text_color,
                **slider_style,
            )
            plotter.add_slider_widget(
                _cb_diffuse,
                rng=[0.0, 1.0],
                value=state["diffuse"],
                title="Diffuse",
                pointa=(0.02, 0.55),
                pointb=(0.22, 0.55),
                color=text_color,
                **slider_style,
            )
            plotter.add_slider_widget(
                _cb_light,
                rng=[0.0, 1.5],
                value=state["light_intensity"],
                title="Light",
                pointa=(0.02, 0.36),
                pointb=(0.22, 0.36),
                color=text_color,
                **slider_style,
            )
            state["light_sliders_visible"] = True
            print(
                "Lighting sliders ON  (Ambient / Diffuse / Light intensity) — press I to hide"
            )
        except Exception as err:
            print(f"Could not create lighting sliders: {err}")
        _render()

    # ------------------------------------------------------------------

    def _print_shortcuts() -> None:
        print(
            "PyVista keys: B background, G bounds grid, S bounding box, "
            "C crop box, X clear crop, M tube/line toggle, L shadows, I lighting sliders, "
            "+/- thickness, A axes, R reset view, P save PNG, V record orbit, H show keys"
        )

    if show_bounds:
        _show_bounds()

    if show_axes:
        _show_axes()

    try:
        plotter.enable_anti_aliasing("fxaa")
    except Exception:
        pass

    plotter.add_key_event("b", _toggle_background)
    plotter.add_key_event("g", _toggle_bounds)
    plotter.add_key_event("s", _toggle_box)
    plotter.add_key_event("c", _toggle_interactive_crop)
    plotter.add_key_event("x", _clear_interactive_crop)
    plotter.add_key_event("l", _toggle_shadows)
    plotter.add_key_event("i", _toggle_light_sliders)
    plotter.add_key_event("m", _toggle_fast_mode)
    plotter.add_key_event("plus", _increase_thickness)
    plotter.add_key_event("equal", _increase_thickness)
    plotter.add_key_event("minus", _decrease_thickness)
    plotter.add_key_event("a", _toggle_axes)
    plotter.add_key_event("r", _reset_view)
    plotter.add_key_event("p", _save_screenshot)
    plotter.add_key_event("v", _record_orbit_interactive)
    plotter.add_key_event("h", _print_shortcuts)

    if state["auto_lod"]:
        start_ok = _add_interaction_observer(
            "StartInteractionEvent", _show_interaction_lod
        )
        end_ok = _add_interaction_observer("EndInteractionEvent", _hide_interaction_lod)
        if start_ok and end_ok:
            print(
                "PyVista auto-fast interaction ON: using lines while moving, tubes when idle."
            )
        else:
            print(
                "PyVista auto-fast interaction unavailable in this PyVista/VTK setup."
            )

    plotter.view_isometric()
    plotter.camera.zoom(1.15)

    if interactive:
        _print_shortcuts()

    if not interactive and not screenshot_path and not video_path:
        raise ValueError(
            "Must specify screenshot_path or video_path when interactive=False."
        )

    if screenshot_path:
        Path(screenshot_path).parent.mkdir(parents=True, exist_ok=True)

    if video_path and not interactive:
        # PyVista MP4 pattern: open_movie → show(auto_close=False) → write_frame × N → close
        # show() must be called first to initialise the renderer before any write_frame().
        plotter.show(auto_close=False)
        _record_orbit(video_path, n_frames=video_frames, fps=video_fps)
        plotter.close()
    else:
        plotter.show(screenshot=screenshot_path, auto_close=True)
