import numpy as np
from fury import actor, colormap, window


def plot_vector_field_fury(
    vector_field,
    size=1.0,
    radius=0.5,
    color_volume=None,
    stride=10,
    voxel_size=1.0,
    mode="arrow",  # "arrow" or "cylinder"
    save_path=None,
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
    """
    print("Starting FURY vector field visualization...")
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
        color_array = colormap.create_colormap(color_values, name="hsv", auto=True)
    else:
        color_array = np.tile([1.0, 0.0, 0.0], (centers.shape[0], 1))

    # Create scene
    scene = window.Scene()

    if mode == "arrow":
        print("Rendering as arrows...")
        arrow_actor = actor.arrow(
            centers,
            directions,
            colors=color_array,
            scales=10000 * size,
        )
        scene.add(arrow_actor)

    elif mode == "cylinder":
        print("Rendering as cylinders...")
        lengths = norms * size
        cylinder_actor = actor.cylinder(
            centers=centers,
            directions=directions,
            colors=color_array,
            heights=size*5000,
            radius=radius*0.1,
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
