import numpy as np
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
    valid_mask = norms > 0
    centers = coords_flat[valid_mask] * voxel_size
    directions = vectors_flat[valid_mask]
    norms = norms[valid_mask]
    directions /= norms[:, None]  # normalize

    print(f"Number of vectors to display: {centers.shape[0]}")

    # Colors
    scalar_lut = None
    if color_volume is not None:
        color_sub = color_volume[
            0:Z:downsample, 0:Y:downsample, 0:X:downsample
        ]
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

    # Create scene
    scene = window.Scene()
    current_size = float(size)
    vector_actor = None

    def rebuild_actor() -> None:
        nonlocal vector_actor
        if vector_actor is not None:
            scene.rm(vector_actor)

        if mode == "arrow":
            scales = norms * voxel_size * current_size
            scales = np.repeat(scales[:, None], 3, axis=1)  # Nx3 for arrows
            vector_actor = actor.arrow(
                centers,
                directions,
                colors=color_array,
                scales=10.0 * scales,
            )
        elif mode == "cylinder":
            heights = norms * voxel_size * current_size
            vector_actor = actor.cylinder(
                centers=centers,
                directions=directions,
                colors=color_array,
                heights=heights * 10,
                radius=radius * voxel_size,
                capped=True,
            )
        else:
            raise ValueError("Mode must be 'arrow' or 'cylinder'")

        scene.add(vector_actor)

    print(f"Rendering as {mode}s...")
    rebuild_actor()
    if scalar_lut is not None:
        scene.add(actor.scalar_bar(lookup_table=scalar_lut, title="Value"))

    # Show or save
    if save_path:
        print(f"Saving FURY vector plot to: {save_path}")
        window.record(scene, out_path=str(save_path), size=(800, 800))
    else:
        print("Displaying interactive scene...")
        showm = window.ShowManager(scene=scene, size=(800, 800), reset_camera=True)

        def on_keypress(obj, _evt):
            nonlocal current_size
            key = obj.GetKeySym().lower()
            if key in ("plus", "equal", "kp_add"):
                current_size *= 1.25
            elif key in ("minus", "underscore", "kp_subtract"):
                current_size = max(0.001, current_size * 0.8)
            else:
                return

            rebuild_actor()
            scene.ResetCameraClippingRange()
            showm.render()
            print(f"Vector size: {current_size:.3f}")

        showm.initialize()
        showm.iren.AddObserver("KeyPressEvent", on_keypress)
        print("Keys: +/- vector size")
        showm.start()
