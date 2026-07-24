from pathlib import Path

import matplotlib.cm as cm
import numpy as np

from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.utils.streamlines_io_utils import load_trk_streamlines

ANGLE_RANGES = {
    "HA": (-90.0, 90.0),
    "IA": (-90.0, 90.0),
    "AZ": (-180.0, 180.0),
    "EL": (0.0, 90.0),
}


def _compute_az_el_from_streamlines(
    streamlines_xyz: list[np.ndarray],
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """
    Derive per-vertex azimuth and elevation from streamline tangents.
    - elevation = arcsin(ẑ) in degrees
    - azimuth   = atan2(ŷ, x̂) in degrees
    The last vertex repeats the previous angle to keep lengths aligned.
    """
    az_list: list[np.ndarray] = []
    el_list: list[np.ndarray] = []
    for pts in streamlines_xyz:
        pts = np.asarray(pts, dtype=np.float32)
        n = len(pts)
        if n < 2:
            az_list.append(np.zeros((n,), dtype=np.float32))
            el_list.append(np.zeros((n,), dtype=np.float32))
            continue

        vecs = np.diff(pts, axis=0)
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        unit = np.divide(vecs, norms, out=np.zeros_like(vecs), where=norms != 0)

        el = np.degrees(np.arcsin(np.clip(unit[:, 2], -1.0, 1.0)))
        az = np.degrees(np.arctan2(unit[:, 1], unit[:, 0]))

        # repeat last value to match vertex count
        el = np.concatenate([el, [el[-1]]]).astype(np.float32)
        az = np.concatenate([az, [az[-1]]]).astype(np.float32)

        el_list.append(el)
        az_list.append(az)
    return az_list, el_list


def visualize_streamlines(
    streamlines_file: str | Path,
    color_by: str = "ha",  # {"ha","ia","az","el","elevation","azimuth"}
    line_width: float = 4.0,
    subsample_factor: int = 1,
    filter_min_len: int | None = None,
    downsample_factor: int = 1,
    max_streamlines: int | None = None,
    top_clusters: int | None = None,
    crop_bounds: tuple | None = None,  # ((xmin,xmax),(ymin,ymax),(zmin,zmax))
    interactive: bool = True,
    screenshot_path: str | None = None,
    video_path: str | None = None,
    video_fps: int = 30,
    video_frames: int = 120,
    window_size: tuple[int, int] = (800, 800),
    colormap=None,
    spline_subdiv: int = 2,
    backend: str = "fury",
    mode: str = "tube",
    background_color: str | tuple[float, float, float] | None = None,
    tube_sides: int = 9,
    pyvista_opacity: float = 1.0,
    pyvista_show_axes: bool = True,
    pyvista_show_bounds: bool = False,
    pyvista_shadows: bool = False,
    random_seed: int | None = None,
    session_path: str | Path | None = None,
    restore_session: bool = True,
    session_settings: dict | None = None,
    fury_quality: str = "interactive",
):
    """
    Visualize .trk streamlines with per-point angle-based coloring.

    Parameters
    ----------
    backend
        Rendering backend: "fury" for tractography-oriented interaction or
        "pyvista" for polished scientific plots and screenshots.
    mode
        "tube" for explicit tube geometry or "line" for faster line rendering.
    """
    p = Path(streamlines_file)
    if not p.exists():
        raise FileNotFoundError(f"Streamlines file not found: {p}")
    if p.suffix.lower() != ".trk":
        raise ValueError("Only .trk input is supported here")

    print(f"Loading .trk streamlines: {p}")
    streamlines_xyz, attrs, streamline_attrs = load_trk_streamlines(
        p, include_per_streamline=True
    )
    attrs = {name.upper(): seq for name, seq in attrs.items()}
    streamline_attrs = {
        name.upper(): values for name, values in streamline_attrs.items()
    }

    if top_clusters is not None:
        if top_clusters < 1:
            raise ValueError("top_clusters must be at least 1")
        if "CLUSTER_ID" not in attrs:
            raise ValueError(
                "--top-clusters requires a clustered TRK with a per-point "
                "cluster_id field"
            )

        cluster_ids = []
        for values in attrs["CLUSTER_ID"]:
            values = np.asarray(values, dtype=np.float64).reshape(-1)
            if (
                not values.size
                or not np.all(np.isfinite(values))
                or not np.allclose(values, values[0])
                or not np.isclose(values[0], round(float(values[0])))
            ):
                raise ValueError(
                    "Each streamline must have one constant integer cluster_id"
                )
            cluster_ids.append(int(round(float(values[0]))))

        cluster_ids = np.asarray(cluster_ids, dtype=np.int64)
        unique_ids, counts = np.unique(cluster_ids, return_counts=True)
        cluster_sizes = {
            int(cluster_id): int(count)
            for cluster_id, count in zip(unique_ids, counts)
        }
        if "CLUSTER_SIZE" in streamline_attrs:
            stored_sizes = np.asarray(
                streamline_attrs["CLUSTER_SIZE"], dtype=np.float64
            ).reshape(len(streamlines_xyz), -1)[:, 0]
            for cluster_id, size in zip(cluster_ids, stored_sizes):
                if np.isfinite(size) and size > 0:
                    cluster_sizes[int(cluster_id)] = max(
                        cluster_sizes[int(cluster_id)], int(round(float(size)))
                    )

        ranked_ids = sorted(
            (int(cluster_id) for cluster_id in unique_ids),
            key=lambda cluster_id: (-cluster_sizes[cluster_id], cluster_id),
        )
        selected_ids = set(ranked_ids[:top_clusters])
        keep = [
            index
            for index, cluster_id in enumerate(cluster_ids)
            if int(cluster_id) in selected_ids
        ]
        print(
            f"Top-cluster filter: keeping {len(selected_ids):,}/"
            f"{len(unique_ids):,} largest clusters and {len(keep):,}/"
            f"{len(streamlines_xyz):,} streamlines"
        )
        streamlines_xyz = [streamlines_xyz[index] for index in keep]
        attrs = {
            name: [values[index] for index in keep]
            for name, values in attrs.items()
        }

    # ---- Inform the user of available stored scalar fields ----
    available = list(attrs.keys())

    print(
        "\n🎨  Available per-point fields in this .trk:",
        available if available else "None stored",
    )
    print(
        "💡 Note: 'az' and 'el' can still be computed on-the-fly from streamline geometry."
    )
    options = ["elevation", "azimuth", "az", "el", *available]
    print(f"🧭 You can use: color_by = {', '.join(options)}\n")

    # Decide the color scalar
    color_mode = color_by.lower().strip()
    stored_fields = {name.lower(): name for name in available}
    color_values: list[np.ndarray] | None = None
    color_range: tuple[float, float] | None = None
    color_label = "Angle (deg)"

    if color_mode in {"ha", "ia", "az", "el"}:
        key = color_mode.upper()
        color_range = ANGLE_RANGES[key]
        color_label = f"{key} (deg)"
        if key in attrs:
            color_values = [
                np.asarray(arr, dtype=np.float32).reshape(-1) for arr in attrs[key]
            ]
        else:
            if key in {"AZ", "EL"}:
                az_list, el_list = _compute_az_el_from_streamlines(streamlines_xyz)
                color_values = (
                    az_list
                    if key == "AZ"
                    else [np.abs(el).astype(np.float32) for el in el_list]
                )
            else:
                raise ValueError(f"No per-point attribute '{key}' found in .trk")
    elif color_mode in {"elevation", "azimuth"}:
        az_list, el_list = _compute_az_el_from_streamlines(streamlines_xyz)
        color_values = el_list if color_mode == "elevation" else az_list
        color_range = (-90.0, 90.0) if color_mode == "elevation" else ANGLE_RANGES["AZ"]
        color_label = (
            "Elevation (deg)" if color_mode == "elevation" else "Azimuth (deg)"
        )
    elif color_mode in stored_fields:
        key = stored_fields[color_mode]
        color_values = []
        field_min = np.inf
        field_max = -np.inf
        for values in attrs[key]:
            values = np.asarray(values, dtype=np.float32)
            if values.ndim != 1:
                raise ValueError(f"Per-point field '{key}' must contain one scalar")
            color_values.append(values)
            finite = values[np.isfinite(values)]
            if finite.size:
                field_min = min(field_min, float(finite.min()))
                field_max = max(field_max, float(finite.max()))
        if not np.isfinite(field_min):
            raise ValueError(f"Per-point field '{key}' contains no finite values")
        if field_min == field_max:
            field_max = field_min + 1.0
        color_range = (field_min, field_max)
        color_label = key
    else:
        raise ValueError(
            f"Unknown color_by '{color_by}'. Available: {', '.join(options)}"
        )

    # Default colormap selection
    if colormap is None:
        if color_mode == "el":
            colormap = cm.viridis
        elif color_mode in {"ha", "ia", "elevation"}:
            colormap = helix_angle_cmap
        elif color_mode in {"az", "azimuth"}:
            colormap = cm.hsv
        else:
            colormap = cm.viridis

    backend = backend.lower().strip()
    render_mode = mode.lower().strip()
    if render_mode not in {"tube", "line"}:
        raise ValueError("mode must be one of: tube, line")

    if backend == "fury":
        from cardiotensor.visualization.fury_plotting_streamlines import (
            show_streamlines,
        )

        show_streamlines(
            streamlines_xyz=streamlines_xyz,
            color_values=color_values,
            mode=render_mode,
            line_width=line_width,
            interactive=interactive,
            screenshot_path=screenshot_path,
            video_path=video_path,
            video_fps=video_fps,
            video_frames=video_frames,
            window_size=window_size,
            downsample_factor=downsample_factor,
            max_streamlines=max_streamlines,
            filter_min_len=filter_min_len,
            subsample_factor=subsample_factor,
            crop_bounds=crop_bounds,
            colormap=colormap,
            background_color=background_color or "black",
            spline_subdiv=spline_subdiv,
            tube_sides=tube_sides,
            color_range=color_range,
            color_label=color_label,
            random_seed=random_seed,
            session_path=session_path,
            restore_session=restore_session,
            session_settings=session_settings,
            streamlines_file=p,
            opacity=pyvista_opacity,
            quality=fury_quality,
        )
    elif backend == "pyvista":
        from cardiotensor.visualization.pyvista_plotting_streamlines import (
            show_streamlines_pyvista,
        )

        show_streamlines_pyvista(
            streamlines_xyz=streamlines_xyz,
            color_values=color_values,
            mode=render_mode,
            line_width=line_width,
            interactive=interactive,
            screenshot_path=screenshot_path,
            video_path=video_path,
            video_fps=video_fps,
            video_frames=video_frames,
            window_size=window_size,
            downsample_factor=downsample_factor,
            max_streamlines=max_streamlines,
            filter_min_len=filter_min_len,
            subsample_factor=subsample_factor,
            crop_bounds=crop_bounds,
            colormap=colormap,
            background_color=background_color or "white",
            tube_sides=tube_sides,
            color_range=color_range,
            color_label=color_label,
            opacity=pyvista_opacity,
            show_axes=pyvista_show_axes,
            show_bounds=pyvista_show_bounds,
            shadows=pyvista_shadows,
            spline_subdiv=spline_subdiv,
            random_seed=random_seed,
        )
    else:
        raise ValueError("backend must be one of: fury, pyvista")
