import numpy as np
from fury import actor, colormap, window


def plot_vector_field_with_fury(
    vector_field,
    size_arrow=1,
    ha_volume=None,
    stride=10,
    voxel_size=1.0,
    save_path=None,
):
    """
    Visualize a 3D vector field using FURY with optional HA-based coloring.

    Parameters:
        vector_field (np.ndarray): 4D array (Z, Y, X, 3) of vectors.
        ha_volume (np.ndarray, optional): 3D array (Z, Y, X) of helix angles in degrees.
        stride (int): Downsampling stride to reduce number of vectors.
        voxel_size (float): Physical voxel size for proper arrow scaling.
        save_path (Path, optional): If provided, save the image to this path.
    """
    print("Starting FURY vector field visualization...")
    Z, Y, X, _ = vector_field.shape

    print("Creating coordinate grid...")
    zz, yy, xx = np.mgrid[0:Z:stride, 0:Y:stride, 0:X:stride]
    coords = np.stack((zz, yy, xx), axis=-1)
    vector_field = vector_field[0:Z:stride, 0:Y:stride, 0:X:stride]

    # Flatten coordinates and vectors
    coords_flat = coords.reshape(-1, 3)
    vectors_flat = vector_field.reshape(-1, 3)
    del vector_field

    print("Extracting and filtering vectors...")
    norms = np.linalg.norm(vectors_flat, axis=1)
    valid_mask = norms > 0

    centers = coords_flat[valid_mask] * voxel_size
    directions = vectors_flat[valid_mask] / norms[valid_mask, None]
    del coords_flat, vectors_flat, norms

    print("Generating colors...")
    if ha_volume is not None:
        ha_sub = ha_volume[0:Z:stride, 0:Y:stride, 0:X:stride]
        ha_flat = ha_sub.reshape(-1)
        ha_values = ha_flat[valid_mask]
        color_array = colormap.create_colormap(ha_values, name="hsv", auto=True)
    else:
        color_array = np.tile([1.0, 0.0, 0.0], (centers.shape[0], 1))

    if centers.shape[0] == 0:
        print("⚠ No arrows to display.")
        return

    print("Creating arrow actor...")
    arrow_actor = actor.arrow(
        centers,
        directions,
        colors=color_array,
        scales=10 * size_arrow,
    )

    scene = window.Scene()
    scene.add(arrow_actor)

    if save_path:
        print(f"Saving FURY vector plot to: {save_path}")
        window.record(scene, out_path=str(save_path), size=(800, 800))
    else:
        print("Displaying interactive scene...")
        window.show(scene, size=(800, 800))
