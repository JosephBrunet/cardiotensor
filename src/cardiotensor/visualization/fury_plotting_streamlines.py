import random

import fury
import matplotlib.pyplot as plt
import numpy as np
import vtk


def downsample_streamline(streamline, factor=2):
    return streamline if len(streamline) < 3 else streamline[::factor]


def matplotlib_cmap_to_fury_lut(cmap, value_range=(-1, 1), n_colors=256):
    """
    Convert a Matplotlib colormap to a VTK LookupTable for FURY.

    Parameters
    ----------
    cmap : str or matplotlib colormap
        Name or object of Matplotlib colormap.
    value_range : (float, float)
        Tuple specifying the (min, max) of scalar values.
    n_colors : int
        Number of LUT entries.

    Returns
    -------
    lut : vtk.vtkLookupTable
        A VTK-compatible lookup table.
    """
    import numpy as np

    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)

    colors = cmap(np.linspace(0, 1, n_colors))  # RGBA in [0,1]
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(n_colors)
    lut.SetRange(*value_range)
    lut.Build()

    for i in range(n_colors):
        r, g, b, a = colors[i]
        lut.SetTableValue(i, r, g, b, a)

    return lut


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
):
    """
    Render 3D streamlines in FURY with per-vertex coloring.

    Parameters
    ----------
    streamlines_xyz : list of np.ndarray
        Each streamline as an array of (x, y, z) points.
    color_values : list of np.ndarray
        Flattened per-vertex scalar values for coloring.
    mode : {"tube", "fake_tube", "line"}
        Rendering mode for streamlines.
    line_width : float
        Thickness of tubes/lines.
    interactive : bool
        Whether to open an interactive FURY window.
    screenshot_path : str, optional
        Save screenshot if not interactive.
    window_size : (int, int)
        Window size for display or screenshot.
    downsample_factor : int
        Step size to reduce points in each streamline.
    max_streamlines : int, optional
        Maximum number of streamlines to display.
    filter_min_len : int, optional
        Minimum streamline length (after downsampling) to keep.
    subsample_factor : int
        Randomly keep 1 of every N streamlines.
    crop_bounds : tuple of ((x_min,x_max),...), optional
        Crop streamlines within bounds.
    colormap : callable or str, optional
        Function or name of Matplotlib colormap to map scalar values to RGB.
        If None, defaults to FURY HSV mapping.
    """
    print(f"Initial number of streamlines: {len(streamlines_xyz)}")
    if filter_min_len:
        print(f"Filtering out streamlines shorter than {filter_min_len} points")
    if downsample_factor > 1:
        print(f"Downsampling each streamline by factor {downsample_factor}")
    if subsample_factor > 1:
        print(f"Subsampling: keeping 1 in every {subsample_factor} streamlines")
    if max_streamlines:
        print(f"Limiting to max {max_streamlines} streamlines")

    # --- Cropping
    print(
        f"Cropping streamlines within bounds: {crop_bounds}"
        if crop_bounds
        else "No cropping applied."
    )
    if crop_bounds is not None:
        z_min, z_max = crop_bounds[2]
        y_min, y_max = crop_bounds[1]
        x_min, x_max = crop_bounds[0]

        new_streamlines = []
        new_color_values = []
        color_idx = 0

        for sl in streamlines_xyz:
            n_pts = len(sl)
            cl = color_values[color_idx : color_idx + n_pts]

            sl = np.asarray(sl)
            cl = np.asarray(cl)

            within = (
                (sl[:, 0] >= x_min)
                & (sl[:, 0] <= x_max)
                & (sl[:, 1] >= y_min)
                & (sl[:, 1] <= y_max)
                & (sl[:, 2] >= z_min)
                & (sl[:, 2] <= z_max)
            )

            if np.any(within):  # Keep only remaining points
                new_sl = sl[within]
                new_cl = cl[within]
                if len(new_sl) > 0:
                    new_streamlines.append(new_sl)
                    new_color_values.append(new_cl)

            color_idx += n_pts

        streamlines_xyz = new_streamlines
        color_values = (
            np.concatenate(new_color_values) if new_color_values else np.array([])
        )

    # --- Downsample and filter
    downsampled_streamlines = []
    downsampled_colors = []
    idx = 0
    for sl in streamlines_xyz:
        color_slice = color_values[idx : idx + len(sl)]
        ds_sl = downsample_streamline(sl, downsample_factor)
        ds_cl = downsample_streamline(color_slice, downsample_factor)

        if filter_min_len is None or len(ds_sl) >= filter_min_len:
            downsampled_streamlines.append(ds_sl)
            downsampled_colors.append(ds_cl)

        idx += len(sl)

    streamlines_xyz = downsampled_streamlines
    color_values = downsampled_colors

    if not streamlines_xyz:
        raise ValueError("❌ No streamlines left after downsampling and filtering.")

    # --- Subsample
    if subsample_factor > 1:
        total = len(streamlines_xyz)
        selected_idx = sorted(random.sample(range(total), total // subsample_factor))
        streamlines_xyz = [streamlines_xyz[i] for i in selected_idx]
        color_values = [color_values[i] for i in selected_idx]

    # --- Cap max
    if max_streamlines is not None and len(streamlines_xyz) > max_streamlines:
        selected_idx = sorted(
            random.sample(range(len(streamlines_xyz)), max_streamlines)
        )
        streamlines_xyz = [streamlines_xyz[i] for i in selected_idx]
        color_values = [color_values[i] for i in selected_idx]

    print(f"Final number of streamlines to render: {len(streamlines_xyz)}")

    if not color_values:
        raise ValueError(
            "❌ No streamlines left after filtering and cropping. Adjust parameters like --crop or --min-length."
        )

    flat_colors = np.concatenate(color_values)
    print(f"Coloring mode: min={flat_colors.min():.2f}, max={flat_colors.max():.2f}")
    print(f"Rendering mode: {mode}")

    min_val = float(flat_colors.min())
    max_val = float(flat_colors.max())

    # Convert your 0–255 values into real [-90, 90] degrees
    flat_colors = flat_colors / 255 * 180 - 90

    if colormap is None:
        # Default HSV colormap from FURY
        lut = fury.actor.colormap_lookup_table(
            scale_range=(min_val, max_val),
            hue_range=(0.7, 0.0),
            saturation_range=(0.5, 1.0),
        )
    else:
        # Convert Matplotlib cmap to VTK LUT
        lut = matplotlib_cmap_to_fury_lut(
            cmap=colormap,
            value_range=(-90, 90),
            n_colors=256,
        )
    # Map scalar values to LUT indices
    lookup_colors = flat_colors

    scene = fury.window.Scene()

    # --- Render according to mode
    if mode == "tube":
        actor = fury.actor.streamtube(
            streamlines_xyz,
            lookup_colors,
            linewidth=line_width,
            opacity=1,
            lookup_colormap=lut,
        )
    elif mode == "fake_tube":
        actor = fury.actor.line(
            streamlines_xyz,
            lookup_colors,
            linewidth=line_width,
            fake_tube=True,
            depth_cue=True,
        )
    elif mode == "line":
        actor = fury.actor.line(
            streamlines_xyz,
            lookup_colors,
            linewidth=line_width,
            fake_tube=False,
            depth_cue=False,
            lookup_colormap=lut,
        )
    else:
        raise ValueError(f"Unknown mode: {mode}")

    scene.add(actor)
    if lut is not None:
        scene.add(fury.actor.scalar_bar(lut))
    scene.reset_camera()

    # radii = 7
    # # Example coordinates of the point (x, y, z)
    # highlight_point = np.array([[6548-6000,8001-7500,17296-17041]])/2  # 3D coords in same space as streamlines
    # # Create a sphere actor for the point
    # sphere_actor = fury.actor.sphere(centers=highlight_point, colors=(1, 0, 0), radii=radii)
    # # Add the point to the scene
    # scene.add(sphere_actor)

    # # Example coordinates of the point (x, y, z)
    # highlight_point = np.array([[0,0,0]])  # 3D coords in same space as streamlines
    # # Create a sphere actor for the point
    # sphere_actor = fury.actor.sphere(centers=highlight_point, colors=(1, 0, 0), radii=radii)
    # # Add the point to the scene
    # scene.add(sphere_actor)

    # # Example coordinates of the point (x, y, z)
    # highlight_point = np.array([[6570-6000,7970-7500,17142-17041]])/2  # 3D coords in same space as streamlines
    # # Create a sphere actor for the point
    # sphere_actor = fury.actor.sphere(centers=highlight_point, colors=(1, 0, 0), radii=radii)
    # # Add the point to the scene
    # scene.add(sphere_actor)

    # # Example coordinates of the point (x, y, z)
    # highlight_point = np.array([[6548-6000,7997-7500,17394-17041]])/2  # 3D coords in same space as streamlines
    # # Create a sphere actor for the point
    # sphere_actor = fury.actor.sphere(centers=highlight_point, colors=(1, 0, 0), radii=radii)
    # # Add the point to the scene
    # scene.add(sphere_actor)

    # --- Display or screenshot
    if interactive:
        print("🕹️ Opening interactive window...")
        fury.window.show(scene, size=window_size, reset_camera=False)
    else:
        if not screenshot_path:
            raise ValueError("Must specify screenshot_path when interactive=False.")
        print(f"📸 Saving screenshot to: {screenshot_path}")
        fury.window.record(scene=scene, out_path=screenshot_path, size=window_size)
