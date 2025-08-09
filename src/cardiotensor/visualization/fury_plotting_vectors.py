import numpy as np
from fury import actor, window

from cardiotensor.colormaps.helix_angle import helix_angle_cmap


def plot_vector_field_fury(
    vector_field: np.ndarray,
    size: float = 1.0,
    radius: float = 0.5,
    color_volume: np.ndarray = None,
    stride: int = 10,
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
    stride : int
        Downsampling stride to reduce number of vectors.
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

    Z, Y, X, _ = vector_field.shape

    # Downsample grid
    zz, yy, xx = np.mgrid[0:Z:stride, 0:Y:stride, 0:X:stride]
    coords = np.stack((zz, yy, xx), axis=-1)
    vector_field = vector_field[0:Z:stride, 0:Y:stride, 0:X:stride]

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
    if color_volume is not None:
        color_sub = color_volume[0:Z:stride, 0:Y:stride, 0:X:stride]
        color_flat = color_sub.reshape(-1)
        color_values = color_flat[valid_mask]

        # Normalize to [0, 1]
        cmin, cmax = np.nanmin(color_values), np.nanmax(color_values)
        color_values = (color_values - cmin) / (cmax - cmin + 1e-8)

        # Map to RGB using the chosen colormap
        color_array = colormap(color_values)[:, :3]  # drop alpha
    else:
        color_array = np.tile([1.0, 0.0, 0.0], (centers.shape[0], 1))

    # Create scene
    scene = window.Scene()

    if mode == "arrow":
        print("Rendering as arrows...")
        scales = norms * voxel_size * size
        scales = np.repeat(scales[:, None], 3, axis=1)  # Nx3 for arrows
        arrow_actor = actor.arrow(
            centers,
            directions,
            colors=color_array,
            scales=10 * scales,
        )
        scene.add(arrow_actor)

    elif mode == "cylinder":
        print("Rendering as cylinders...")
        heights = norms * voxel_size * size
        cylinder_actor = actor.cylinder(
            centers=centers,
            directions=directions,
            colors=color_array,
            heights=heights * 10,
            radius=radius * voxel_size * 0.0002,
            capped=True,
        )
        scene.add(cylinder_actor)

    else:
        raise ValueError("Mode must be 'arrow' or 'cylinder'")

    # Show or save
    if save_path:
        print(f"Saving FURY vector plot to: {save_path}")
        window.record(scene, out_path=str(save_path), size=(800, 800))
    else:
        print("Displaying interactive scene...")
        window.show(scene, size=(800, 800))
