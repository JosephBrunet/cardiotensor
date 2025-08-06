from pathlib import Path

import numpy as np

from cardiotensor.visualization.fury_plotting_streamlines import (
    show_streamlines,
)


def visualize_streamlines(
    streamlines_file: str | Path,
    color_by: str = "ha",
    mode: str = "tube",
    line_width: float = 4.0,
    subsample_factor: int = 1,
    filter_min_len: int | None = None,
    downsample_factor: int = 1,
    max_streamlines: int | None = None,
    crop_bounds: tuple | None = None,
    interactive: bool = True,
    screenshot_path: str | None = None,
    window_size: tuple[int, int] = (800, 800),
):
    """
    Visualize precomputed streamlines with optional coloring.

    Parameters
    ----------
    streamlines_file : str or Path
        Path to a .npz file containing streamlines and optional ha_values.
    color_by : {"ha", "elevation"}
        Color streamlines by Helix Angle (HA) or elevation angle.
    mode : {"tube", "fake_tube", "line"}
        Rendering mode for streamlines.
    line_width : float
        Tube/line thickness.
    subsample_factor : int
        Factor to subsample streamline points for rendering speed.
    filter_min_len : int, optional
        Minimum streamline length to display.
    downsample_factor : int
        Downsampling factor for streamlines (visual only).
    max_streamlines : int, optional
        Maximum number of streamlines to visualize.
    crop_bounds : tuple, optional
        ((x_min, x_max), (y_min, y_max), (z_min, z_max)) for cropping.
    interactive : bool
        If False, window will close immediately after rendering or screenshot.
    screenshot_path : str or Path, optional
        Save a screenshot instead of opening interactive window.
    window_size : (int, int)
        Window width and height in pixels.
    """
    streamlines_file = Path(streamlines_file)
    if not streamlines_file.exists():
        raise FileNotFoundError(f"❌ Streamlines file not found: {streamlines_file}")

    data = np.load(streamlines_file, allow_pickle=True)
    raw_streamlines = data.get("streamlines")
    if raw_streamlines is None:
        raise ValueError("'streamlines' array missing in .npz.")

    # Convert (z, y, x) to (x, y, z)
    streamlines_xyz = [
        np.array([(pt[2], pt[1], pt[0]) for pt in sl], dtype=np.float32)
        for sl in raw_streamlines.tolist()
    ]

    # Compute coloring
    if color_by == "elevation":
        color_values = compute_elevation_angles(streamlines_xyz)
    else:  # color_by == "ha"
        ha_values = data.get("ha_values")
        if ha_values is None:
            raise ValueError("'ha_values' array missing in .npz.")
        color_values = ha_values.astype(np.float32)

    # Render streamlines
    show_streamlines(
        streamlines_xyz=streamlines_xyz,
        color_values=color_values,
        mode=mode,
        line_width=line_width,
        interactive=interactive,
        screenshot_path=screenshot_path,
        window_size=window_size,
        downsample_factor=downsample_factor,
        max_streamlines=max_streamlines,
        filter_min_len=filter_min_len,
        subsample_factor=subsample_factor,
        crop_bounds=crop_bounds,
    )


def compute_elevation_angles(streamlines):
    """Compute per-vertex elevation angle as a flat 1D array."""
    all_angles = []
    for pts in streamlines:
        if len(pts) < 2:
            all_angles.append(np.zeros((len(pts),), dtype=np.float32))
            continue
        vecs = np.diff(pts, axis=0)
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        normalized = np.divide(vecs, norms, where=norms != 0)
        z_components = normalized[:, 2]
        elev = np.arcsin(z_components) * 180.0 / np.pi
        elev = np.concatenate([elev, [elev[-1]]])  # match #points
        all_angles.append(elev.astype(np.float32))

    # flatten into a single array
    return np.hstack(all_angles).astype(np.float32)
