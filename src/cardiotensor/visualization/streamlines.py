from pathlib import Path

import matplotlib.cm as cm
import numpy as np

from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.utils.streamlines_io_utils import load_trk_streamlines
from cardiotensor.visualization.fury_plotting_streamlines import show_streamlines


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
    crop_bounds: tuple | None = None,  # ((xmin,xmax),(ymin,ymax),(zmin,zmax))
    interactive: bool = True,
    screenshot_path: str | None = None,
    window_size: tuple[int, int] = (800, 800),
    colormap=None,
    spline_subdiv: int = 16,
):
    """
    Visualize .trk streamlines with per-point angle-based coloring.
    Always renders in tube mode.
    """
    p = Path(streamlines_file)
    if not p.exists():
        raise FileNotFoundError(f"Streamlines file not found: {p}")
    if p.suffix.lower() != ".trk":
        raise ValueError("Only .trk input is supported here")

    print(f"Loading .trk streamlines: {p}")
    streamlines_xyz, attrs = load_trk_streamlines(
        p
    )  # attrs is dict[str, List[np.ndarray]]
    attrs = {
        name.upper(): [np.asarray(arr, dtype=np.float32) for arr in seq]
        for name, seq in attrs.items()
    }

    # ---- Inform the user of available angle fields ----
    available = list(attrs.keys())

    # Also say that AZ and EL can be computed even if missing
    print(
        "\n🎨  Available angle fields in this .trk:",
        available if available else "None stored",
    )
    print(
        "💡 Note: 'az' and 'el' can still be computed on-the-fly from streamline geometry."
    )
    print("🧭 You can use: color_by = ha, ia, az, el, elevation, azimuth\n")

    # Decide the color scalar
    mode = color_by.lower().strip()
    color_values: list[np.ndarray] | None = None
    color_range: tuple[float, float] | None = None
    color_label = "Angle (deg)"

    if mode in {"ha", "ia", "az", "el"}:
        key = mode.upper()
        color_range = ANGLE_RANGES[key]
        color_label = f"{key} (deg)"
        if key in attrs:
            color_values = attrs[key]
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
    elif mode in {"elevation", "azimuth"}:
        az_list, el_list = _compute_az_el_from_streamlines(streamlines_xyz)
        color_values = el_list if mode == "elevation" else az_list
        color_range = (-90.0, 90.0) if mode == "elevation" else ANGLE_RANGES["AZ"]
        color_label = "Elevation (deg)" if mode == "elevation" else "Azimuth (deg)"
    else:
        raise ValueError("color_by must be one of: ha, ia, az, el, elevation, azimuth")

    # Default colormap selection
    if colormap is None:
        if mode == "el":
            colormap = cm.viridis
        elif mode in {"ha", "ia", "elevation"}:
            colormap = helix_angle_cmap
        else:
            colormap = cm.hsv

    # Always use tube mode
    show_streamlines(
        streamlines_xyz=streamlines_xyz,
        color_values=color_values,
        mode="tube",
        line_width=line_width,
        interactive=interactive,
        screenshot_path=screenshot_path,
        window_size=window_size,
        downsample_factor=downsample_factor,
        max_streamlines=max_streamlines,
        filter_min_len=filter_min_len,
        subsample_factor=subsample_factor,
        crop_bounds=crop_bounds,
        colormap=colormap,
        spline_subdiv=spline_subdiv,
        color_range=color_range,
        color_label=color_label,
    )
