from pathlib import Path

import matplotlib.cm as cm
import numpy as np

from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.utils.streamlines_io_utils import load_trk_streamlines
from cardiotensor.visualization.fury_plotting_streamlines import show_streamlines


def _normalize_attrs_to_degrees(
    attrs: dict[str, list[np.ndarray]],
) -> dict[str, list[np.ndarray]]:
    """
    Normalize known angle fields to degrees.

    Rules
    - HA, IA, EL: if values look like 0..255, map to -90..90 via (v/255)*180 - 90
                  else keep as float32
    - AZ:  if values look like 0..255, map to 0..360 via (v/255)*360
           else keep as float32
    Other fields are passed through unchanged.
    """
    out: dict[str, list[np.ndarray]] = {}
    for name, seq in attrs.items():
        key = name.upper()
        norm_list: list[np.ndarray] = []
        for arr in seq:
            a = np.asarray(arr)
            if a.size == 0:
                norm_list.append(a.astype(np.float32))
                continue

            mx = float(np.nanmax(a))
            if key in {"HA", "IA", "EL"}:
                if mx > 1.5:  # looks like byte scale
                    a = (a.astype(np.float32) / 255.0) * 180.0 - 90.0
                else:
                    a = a.astype(np.float32)
            elif key == "AZ":
                if mx > 6.0:  # likely 0..255
                    a = (a.astype(np.float32) / 255.0) * 360.0
                else:
                    a = a.astype(np.float32)
            else:
                a = a.astype(np.float32)

            norm_list.append(a)
        out[key] = norm_list
    return out


def _compute_az_el_from_streamlines(
    streamlines_xyz: list[np.ndarray],
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """
    Derive per-vertex azimuth and elevation from streamline tangents.
    - elevation = arcsin(ẑ) in degrees
    - azimuth   = atan2(ŷ, x̂) in degrees mapped to [0, 360)
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
        az = np.mod(az, 360.0)

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

    # Normalize known attributes to degrees
    attrs_deg: dict[str, list[np.ndarray]] = _normalize_attrs_to_degrees(attrs)

    # ---- Inform the user of available angle fields ----
    available = list(attrs_deg.keys())

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

    if mode in {"ha", "ia", "az", "el"}:
        key = mode.upper()
        if key in attrs_deg:
            color_values = attrs_deg[key]
        else:
            if key in {"AZ", "EL"}:
                az_list, el_list = _compute_az_el_from_streamlines(streamlines_xyz)
                color_values = az_list if key == "AZ" else el_list
            else:
                raise ValueError(f"No per-point attribute '{key}' found in .trk")
    elif mode in {"elevation", "azimuth"}:
        az_list, el_list = _compute_az_el_from_streamlines(streamlines_xyz)
        color_values = el_list if mode == "elevation" else az_list
    else:
        raise ValueError("color_by must be one of: ha, ia, az, el, elevation, azimuth")

    # Default colormap selection
    if colormap is None:
        if mode in {"ha", "ia", "el", "elevation"}:
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
    )
