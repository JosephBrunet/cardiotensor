import numpy as np
from matplotlib.colors import to_rgb
import vtk
from fury import actor, window

from cardiotensor.colormaps.helix_angle import helix_angle_cmap


def matplotlib_cmap_to_fury_lut(
    cmap, value_range: tuple[float, float], n_colors: int = 256
) -> vtk.vtkLookupTable:
    """Convert a Matplotlib colormap to a VTK lookup table for FURY scalar bars."""
    colors = cmap(np.linspace(0, 1, n_colors))
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(n_colors)
    lut.SetRange(*value_range)
    lut.Build()
    for i, (r, g, b, a) in enumerate(colors):
        lut.SetTableValue(i, float(r), float(g), float(b), float(a))
    return lut


def plot_vector_field_fury(
    vector_field: np.ndarray,
    size: float = 1.0,
    radius: float = 0.5,
    color_volume: np.ndarray = None,
    downsample: int = 10,
    voxel_size: float = 1.0,
    mode: str = "arrow",  # "arrow" or "cylinder"
    save_path: str = None,
    colormap=None,
    opacity_mask: np.ndarray | None = None,
    outside_opacity: float = 0.05,
    mask_color: str | None = None,
    outside_color: str | None = None,
):
    """
    Visualize a 3D vector field using FURY as arrows or cylinders.

    Parameters
    ----------
    vector_field : np.ndarray
        4D array (Z, Y, X, 3) of vectors.
    size : float
        Scaling factor for arrow/cylinder lengths.
    radius : float
        Radius of the cylinders (ignored in arrow mode).
    color_volume : np.ndarray, optional
        3D array (Z, Y, X) of scalar values for coloring.
    downsample : int
        Display one vector every `downsample` voxels in Z, Y, and X.
    voxel_size : float
        Physical voxel size for proper scaling.
    mode : str
        Visualization mode: "arrow" or "cylinder".
    save_path : Path or str, optional
        If provided, save the screenshot to this path.
    colormap : matplotlib colormap, optional
        Colormap to use for coloring the vectors.
        Default is helix_angle_cmap.
    opacity_mask : np.ndarray, optional
        Binary (Z, Y, X) mask. Inside glyphs are opaque.
    outside_opacity : float, optional
        Opacity of glyphs outside opacity_mask.
    mask_color : str, optional
        Fixed color for glyphs inside opacity_mask.
    outside_color : str, optional
        Fixed color for glyphs outside opacity_mask.
    """
    print("Starting FURY vector field visualization...")

    # Default colormap
    if colormap is None:
        colormap = helix_angle_cmap

    downsample = max(1, int(downsample))
    Z, Y, X, _ = vector_field.shape

    # Downsample grid
    zz, yy, xx = np.mgrid[0:Z:downsample, 0:Y:downsample, 0:X:downsample]
    coords = np.stack((xx, yy, zz), axis=-1)
    vector_field = vector_field[0:Z:downsample, 0:Y:downsample, 0:X:downsample]

    # Flatten
    coords_flat = coords.reshape(-1, 3)
    vectors_flat = vector_field.reshape(-1, 3)
    del vector_field

    # Filter valid vectors
    norms = np.linalg.norm(vectors_flat, axis=1)
    valid_mask = np.isfinite(norms) & (norms > 0)
    centers = coords_flat[valid_mask] * voxel_size
    directions = vectors_flat[valid_mask]
    norms = norms[valid_mask]
    directions /= norms[:, None]  # normalize

    print(f"Number of vectors to display: {centers.shape[0]}")

    if centers.shape[0] == 0:
        print(
            "No finite, nonzero vectors remain in the sampled selection. "
            "Check the --start/--end slice range and mask coverage, or try "
            "a smaller --downsample value."
        )
        return

    # Colors
    scalar_lut = None
    if color_volume is not None:
        color_sub = color_volume[0:Z:downsample, 0:Y:downsample, 0:X:downsample]
        color_flat = color_sub.reshape(-1)
        color_values = color_flat[valid_mask]

        # Normalize to [0, 1]
        cmin, cmax = np.nanmin(color_values), np.nanmax(color_values)
        scalar_lut = matplotlib_cmap_to_fury_lut(colormap, (float(cmin), float(cmax)))
        color_values = (color_values - cmin) / (cmax - cmin + 1e-8)

        # Map to RGB using the chosen colormap
        color_array = colormap(color_values)[:, :3]  # drop alpha
    else:
        color_array = np.tile([1.0, 0.0, 0.0], (centers.shape[0], 1))

    highlighted = None
    if opacity_mask is not None:
        if opacity_mask.shape != (Z, Y, X):
            raise ValueError(
                f"opacity_mask shape {opacity_mask.shape} does not match {(Z, Y, X)}"
            )
        mask_sub = opacity_mask[0:Z:downsample, 0:Y:downsample, 0:X:downsample]
        highlighted = mask_sub.reshape(-1)[valid_mask].astype(bool)
        print(
            f"Opaque vectors inside mask: {np.count_nonzero(highlighted):,}/"
            f"{len(highlighted):,}"
        )

    # Create scene
    scene = window.Scene()
    current_size = float(size)
    current_radius = float(radius)
    vector_actor = None
    highlight_actor = None

    def rebuild_actor() -> None:
        nonlocal vector_actor, highlight_actor
        for current_actor in (vector_actor, highlight_actor):
            if current_actor is not None:
                scene.rm(current_actor)
        highlight_actor = None

        def make_actor(selected, color_override=None):
            selected_centers = centers[selected]
            selected_directions = directions[selected]
            selected_colors = color_array[selected]
            if color_override is not None:
                selected_colors = np.tile(
                    to_rgb(color_override), (np.count_nonzero(selected), 1)
                )
            selected_norms = norms[selected]
            if mode == "arrow":
                scales = selected_norms * voxel_size * current_size
                scales = np.repeat(scales[:, None], 3, axis=1)
                return actor.arrow(
                    selected_centers, selected_directions,
                    colors=selected_colors, scales=10.0 * scales,
                )
            if mode == "cylinder":
                heights = selected_norms * voxel_size * current_size
                return actor.cylinder(
                    centers=selected_centers, directions=selected_directions,
                    colors=selected_colors, heights=heights * 10,
                    radius=current_radius * voxel_size, capped=True,
                )
            raise ValueError("Mode must be 'arrow' or 'cylinder'")

        if highlighted is None:
            vector_actor = make_actor(
                np.ones(len(centers), dtype=bool), outside_color
            )
            scene.add(vector_actor)
            return

        outside = ~highlighted
        if np.any(outside):
            vector_actor = make_actor(outside, outside_color)
            vector_actor.GetProperty().SetOpacity(float(outside_opacity))
            scene.add(vector_actor)
        else:
            vector_actor = None

        if np.any(highlighted):
            highlight_actor = make_actor(highlighted, mask_color)
            highlight_actor.GetProperty().SetOpacity(1.0)
            scene.add(highlight_actor)

    print(f"Rendering as {mode}s...")
    rebuild_actor()
    if scalar_lut is not None and not (mask_color and outside_color):
        scene.add(actor.scalar_bar(lookup_table=scalar_lut, title="Value"))

    # Show or save
    if save_path:
        print(f"Saving FURY vector plot to: {save_path}")
        window.record(scene, out_path=str(save_path), size=(800, 800))
    else:
        print("Displaying interactive scene...")
        showm = window.ShowManager(scene=scene, size=(800, 800), reset_camera=True)

        def on_keypress(obj, _evt):
            nonlocal current_size, current_radius
            key = obj.GetKeySym().lower()
            if key in ("plus", "equal", "kp_add"):
                current_size *= 1.25
            elif key in ("minus", "underscore", "kp_subtract"):
                current_size = max(0.001, current_size * 0.8)
            elif mode == "cylinder" and key in ("bracketright", "rightbracket"):
                current_radius *= 1.25
            elif mode == "cylinder" and key in ("bracketleft", "leftbracket"):
                current_radius = max(0.001, current_radius * 0.8)
            else:
                return

            rebuild_actor()
            scene.ResetCameraClippingRange()
            showm.render()
            if mode == "cylinder":
                print(
                    f"Cylinder length scale: {current_size:.3f}; "
                    f"radius: {current_radius:.3f}; "
                    f"diameter: {2 * current_radius:.3f}"
                )
            else:
                print(f"Arrow size: {current_size:.3f}")

        showm.initialize()
        showm.iren.AddObserver("KeyPressEvent", on_keypress)
        if mode == "cylinder":
            print("Keys: +/- cylinder length; [/] cylinder diameter")
        else:
            print("Keys: +/- arrow size")
        showm.start()
