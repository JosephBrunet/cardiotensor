import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from structure_tensor.multiprocessing import parallel_structure_tensor_analysis
from tqdm import tqdm

from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.utils.image_io import write_jp2, write_tif
from cardiotensor.utils.utils import (
    convert_to_8bit,
    get_available_cpu_count,
    get_gpu_count,
)


def interpolate_points(
    points: list[tuple[float, float, float]], N_img: int
) -> np.ndarray:
    """
    Generates interpolated points using cubic spline interpolation for a given set of 3D points.

    Args:
        points (list[tuple[float, float, float]]): A list of (x, y, z) points.
        N_img (int): The number of slices in the z-dimension.

    Returns:
        np.ndarray: Array of interpolated points.
    """
    if len(points) < 2:
        raise ValueError("At least two points are required for interpolation.")

    # Sort based on the third element (z-coordinate)
    points = sorted(points, key=lambda p: p[2])

    # Extract x, y, z coordinates separately
    points_array = np.array(points)
    x_vals, y_vals, z_vals = points_array[:, 0], points_array[:, 1], points_array[:, 2]

    # Define cubic splines for x and y based on given z values
    cs_x = CubicSpline(z_vals, x_vals, bc_type="natural")
    cs_y = CubicSpline(z_vals, y_vals, bc_type="natural")

    # Generate integer z-values from 1 to N_img
    z_interp = np.arange(0, N_img)

    # Compute interpolated x and y values at integer z positions
    x_interp = cs_x(z_interp)
    y_interp = cs_y(z_interp)

    # Stack into an Nx3 array
    interpolated_points = np.column_stack((x_interp, y_interp, z_interp))

    return interpolated_points


def calculate_center_vector(points: np.ndarray) -> np.ndarray:
    """Compute the linear regression vector for a given set of 3D points.

    Args:
        points (np.ndarray): An Nx3 array of (x, y, z) coordinates representing the curved line.

    Returns:
        np.ndarray: A single 3D unit vector representing the direction of the best-fit line.
    """
    if points.shape[1] != 3:
        raise ValueError("Input must be an Nx3 array of (x, y, z) coordinates.")

    # Compute the centroid (mean position of all points)
    centroid = np.mean(points, axis=0)

    # Center the points by subtracting the centroid
    centered_points = points - centroid

    # Perform Singular Value Decomposition (SVD)
    # This decomposes the data into principal components
    _, _, vh = np.linalg.svd(centered_points)

    center_vec = vh[0] / np.linalg.norm(vh[0])

    # Keep the fitted centerline direction in the same x/y/z frame as AXIS_POINTS
    # and the orientation vector components used by the angle helpers.
    for component in (2, 1, 0):
        if abs(center_vec[component]) > 1e-8:
            center_vec = center_vec if center_vec[component] > 0 else -center_vec
            break

    return center_vec


def adjust_start_end_index(
    start_index: int,
    end_index: int,
    N_img: int,
    padding_start: int = 0,
    padding_end: int = 0,
    is_test: bool = False,
    n_slice: int = 0,
) -> tuple[int, int]:
    """
    Adjusts start and end indices for image processing, considering padding and test mode.

    Args:
        start_index (int): The initial start index.
        end_index (int): The initial end index.
        N_img (int): Number of images in the volume data.
        padding_start (int): Padding to add at the start.
        padding_end (int): Padding to add at the end.
        is_test (bool): Flag indicating whether in test mode.
        n_slice (int): Test slice index.

    Returns:
        Tuple[int, int]: Adjusted start and end indices.
    """

    # Adjust indices for test condition
    if is_test:
        test_index = n_slice
        # else:
        #     test_index = int(N_img / 1.68)
        #     # test_index = 1723

        start_index_padded = max(test_index - padding_start, 0)
        end_index_padded = min(test_index + 1 + padding_end, N_img)
    else:
        # Adjust start and end indices considering padding
        start_index_padded = max(start_index - padding_start, 0)
        end_index_padded = min(end_index + padding_end, N_img)

    return start_index_padded, end_index_padded


def calculate_structure_tensor(
    volume: np.ndarray,
    sigma: float,
    rho: float,
    truncate: float = 4.0,
    devices: list[str] | None = None,
    block_size: int = 200,
    use_gpu: bool = False,
    dtype: type = np.float32,  # Default to np.float64
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculates the structure tensor of a volume.

    Args:
        volume (np.ndarray): The 3D volume data.
        sigma (float): sigma value for Gaussian smoothing.
        rho (float): rho value for Gaussian smoothing.
        devices (Optional[list[str]]): List of devices for parallel processing (e.g., ['cpu', 'cuda:0']).
        block_size (int): Size of the blocks for processing. Default is 200.
        use_gpu (bool): If True, uses GPU for calculations. Default is False.

    Returns:
        tuple[np.ndarray, np.ndarray]: Eigenvalues and eigenvectors of the structure tensor.
    """
    # Filter or ignore specific warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    num_cpus = get_available_cpu_count(default=4)

    devices = devices or []
    num_gpus = 0

    if use_gpu:
        print("🔍 Checking for GPU support...")
        try:
            num_gpus = get_gpu_count()
            print(f"Detected {num_gpus} GPU(s)")
        except Exception as e:
            use_gpu = False
            print(
                f"⚠️ GPU not available or failed to initialize. Using CPU. Reason: {e}"
            )

    if not devices:
        if use_gpu and num_gpus > 0:
            print(f"Using {num_gpus} GPUs for computation")
            devices = []
            for i in range(num_gpus):
                devices.extend([f"cuda:{i}"] * 16)
        else:
            print(f"Using {num_cpus} CPU for computation")
            devices = ["cpu"] * num_cpus

    print("\nStarting structure tensor computation...")
    print(f"---  Volume shape: {volume.shape}")
    print(f"---  sigma: {sigma}, rho: {rho}, Block size: {block_size}")
    if use_gpu and num_gpus > 0:
        device_str = f"{num_gpus} GPU{'s' if num_gpus > 1 else ''}"
    else:
        device_str = f"{num_cpus} CPU{'s' if num_cpus > 1 else ''}"
    print(f"---  Devices: {device_str}")

    class TqdmTotal(tqdm):
        def update_with_total(self, n=1, total=None):
            if total is not None:
                self.total = total
            return self.update(1)

    def run_structure_tensor(selected_devices: list[str]):
        with TqdmTotal(desc="Computing structure tensors", unit="block") as t:
            return parallel_structure_tensor_analysis(
                volume,
                sigma,
                rho,
                devices=selected_devices,
                block_size=block_size,
                truncate=truncate,
                structure_tensor=None,
                eigenvectors=dtype,
                eigenvalues=dtype,
                progress_callback_fn=t.update_with_total,
            )

    try:
        S, val, vec = run_structure_tensor(devices)
    except Exception as e:
        if not any(str(device).startswith("cuda") for device in devices):
            raise

        print(
            f"⚠️ GPU structure tensor computation failed. Retrying on CPU. Reason: {e}"
        )
        devices = ["cpu"] * num_cpus
        S, val, vec = run_structure_tensor(devices)

    print("Structure tensor computation completed\n")

    # vec has shape =(3,z,y,x) in the order of (x,y,z)

    return val, vec


def remove_padding(
    volume: np.ndarray,
    val: np.ndarray,
    vec: np.ndarray,
    padding_start: int,
    padding_end: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Removes padding from the volume, eigenvalues, and eigenvectors.

    Args:
        volume (np.ndarray): The 3D volume data.
        val (np.ndarray): The eigenvalues.
        vec (np.ndarray): The eigenvectors.
        padding_start (int): Padding at the start to remove.
        padding_end (int): Padding at the end to remove.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Adjusted data without padding.
    """
    array_end = vec.shape[1] - padding_end
    volume = volume[padding_start:array_end, :, :]
    vec = vec[:, padding_start:array_end, :, :]
    val = val[:, padding_start:array_end, :, :]

    return volume, val, vec


def compute_fraction_anisotropy(eigenvalues_2d: np.ndarray) -> np.ndarray:
    """
    Computes Fractional Anisotropy (FA) from eigenvalues of a structure tensor.

    Args:
        eigenvalues_2d (np.ndarray): 2D array of eigenvalues (l1, l2, l3).

    Returns:
        np.ndarray: Fractional Anisotropy values.
    """
    l1 = eigenvalues_2d[0, :, :]
    l2 = eigenvalues_2d[1, :, :]
    l3 = eigenvalues_2d[2, :, :]
    mean_eigenvalue = (l1 + l2 + l3) / 3
    numerator = np.sqrt(
        (l1 - mean_eigenvalue) ** 2
        + (l2 - mean_eigenvalue) ** 2
        + (l3 - mean_eigenvalue) ** 2
    )
    denominator = np.sqrt(l1**2 + l2**2 + l3**2) + 1e-10
    img_FA = np.sqrt(3 / 2) * (numerator / denominator)

    return img_FA


def rotate_vectors_to_new_axis(
    vector_field_slice: np.ndarray, new_axis_vec: np.ndarray
) -> np.ndarray:
    """
    Rotate vectors into a local frame where ``new_axis_vec`` is aligned to +Z.

    Args:
        vector_field_slice (np.ndarray): Array of shape (3, Y, X) for a slice.
        new_axis_vec (np.ndarray): Local center axis direction in x/y/z order.

    Returns:
        np.ndarray: Rotated vectors with the same shape as input.
    """
    axis_norm = np.linalg.norm(new_axis_vec)
    if axis_norm == 0:
        raise ValueError("new_axis_vec must be non-zero")

    source_axis = np.asarray(new_axis_vec, dtype=np.float64) / axis_norm
    target_axis = np.array([0.0, 0.0, 1.0], dtype=np.float64)

    rotation_axis = np.cross(source_axis, target_axis)
    axis_length = np.linalg.norm(rotation_axis)
    cos_angle = float(np.clip(np.dot(source_axis, target_axis), -1.0, 1.0))

    if np.isclose(cos_angle, 1.0, atol=1e-10):
        rotation_matrix = np.eye(3, dtype=np.float64)
    elif np.isclose(cos_angle, -1.0, atol=1e-10):
        rotation_matrix = np.diag([1.0, -1.0, -1.0])
    else:
        kx, ky, kz = rotation_axis
        skew_matrix = np.array(
            [[0.0, -kz, ky], [kz, 0.0, -kx], [-ky, kx, 0.0]],
            dtype=np.float64,
        )
        rotation_matrix = (
            np.eye(3, dtype=np.float64)
            + skew_matrix
            + np.dot(skew_matrix, skew_matrix) * ((1.0 - cos_angle) / (axis_length**2))
        )

    # Flatten the vector field to shape (3, N)
    vec_reshaped = np.reshape(vector_field_slice, (3, -1))

    # Normalize safely
    norms = np.linalg.norm(vec_reshaped, axis=0)
    nonzero_mask = norms > 0
    vec_reshaped[:, nonzero_mask] /= norms[nonzero_mask]

    # Rotate
    rotated_vecs = np.dot(rotation_matrix, vec_reshaped)

    # Reshape back
    rotated_vecs = rotated_vecs.reshape(vector_field_slice.shape)

    return rotated_vecs


def orient_vectors_y_positive(vector_field_slice: np.ndarray) -> np.ndarray:
    """
    Flip unoriented vectors so their local Y component is non-negative.

    Structure-tensor eigenvectors are axes: ``v`` and ``-v`` represent the same
    local fiber direction. This helper chooses a deterministic polarity after
    vectors have been rotated into the analysis frame.
    """
    oriented = np.array(vector_field_slice, copy=True)
    flip_mask = oriented[1, :, :] < 0
    oriented[:, flip_mask] *= -1
    return oriented


def orient_vectors_z_positive(vector_field_slice: np.ndarray) -> np.ndarray:
    """
    Flip unoriented vectors so their local Z component is non-negative.

    This is useful for azimuth/elevation maps when the desired elevation is an
    unsigned angle above the local XY plane. Structure-tensor eigenvectors are
    axes, so ``v`` and ``-v`` represent the same local orientation.
    """
    oriented = np.array(vector_field_slice, copy=True)
    flip_mask = oriented[2, :, :] < 0
    oriented[:, flip_mask] *= -1
    return oriented


def compute_helix_and_transverse_angles(
    vector_field_2d: np.ndarray,
    center_point: tuple[int, int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes helix and transverse angles from a 2D vector field.

    Args:
        vector_field_2d (np.ndarray): 2D orientation vector field.
        center_point (Tuple[int, int, int]): Coordinates of the center point.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Helix and transverse angle arrays.
    """
    center = center_point[0:2]  # Replace with actual values
    rows, cols = vector_field_2d.shape[1:3]

    reshaped_vector_field = np.reshape(vector_field_2d, (3, -1))

    center_x, center_y = center[0], center[1]

    X, Y = np.meshgrid(np.arange(cols) - center_x, np.arange(rows) - center_y)

    theta = -np.arctan2(Y.flatten(), X.flatten())
    cos_angle = np.cos(theta)
    sin_angle = np.sin(theta)

    # Change coordinate system to cylindrical
    rotated_vector_field = np.copy(reshaped_vector_field)
    rotated_vector_field[0, :] = (
        cos_angle * reshaped_vector_field[0, :]
        - sin_angle * reshaped_vector_field[1, :]
    )
    rotated_vector_field[1, :] = (
        sin_angle * reshaped_vector_field[0, :]
        + cos_angle * reshaped_vector_field[1, :]
    )

    # Reshape rotated vector field to original image dimensions
    reshaped_rotated_vector_field = np.zeros((3, rows, cols))
    for i in range(3):
        reshaped_rotated_vector_field[i] = rotated_vector_field[i].reshape(rows, cols)

    reshaped_rotated_vector_field = orient_vectors_y_positive(
        reshaped_rotated_vector_field
    )

    # Calculate helix and transverse angles
    helix_angle = np.arctan2(
        reshaped_rotated_vector_field[2, :, :],
        reshaped_rotated_vector_field[1, :, :],
    )
    transverse_angle = np.arctan2(
        reshaped_rotated_vector_field[0, :, :],
        reshaped_rotated_vector_field[1, :, :],
    )
    helix_angle = np.rad2deg(helix_angle)
    transverse_angle = np.rad2deg(transverse_angle)

    return helix_angle, transverse_angle


def compute_azimuth_and_elevation(
    vector_field_2d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Azimuth = angle in XY plane from +X toward +Y, in degrees [-180, 180]
    Elevation = unsigned angle from XY plane toward +Z, in degrees [0, 90]
    """
    oriented = orient_vectors_z_positive(vector_field_2d)
    vx = oriented[0, :, :]
    vy = oriented[1, :, :]
    vz = oriented[2, :, :]

    az = np.rad2deg(np.arctan2(vy, vx))  # [-180, 180]
    el = np.rad2deg(np.arctan2(vz, np.hypot(vx, vy)))  # [0, 90]

    return az, el


def _plot_vector_overlay(
    ax,
    vector_field_slice: np.ndarray,
    quiver_step: int | None = None,
    quiver_length: float | None = None,
    quiver_color: str = "cyan",
    scalar_map: np.ndarray | None = None,
    colormap=None,
) -> None:
    """
    Overlay a subsampled in-plane vector field on an existing axes.

    The field is expected in shape ``(3, Y, X)`` where the first two channels
    encode the in-plane X and Y components.
    """
    if vector_field_slice.ndim != 3 or vector_field_slice.shape[0] != 3:
        raise ValueError("vector_field_slice must have shape (3, Y, X)")

    _, rows, cols = vector_field_slice.shape
    if scalar_map is not None and scalar_map.shape != (rows, cols):
        raise ValueError("scalar_map shape must match the vector field slice shape")

    if quiver_step is None:
        quiver_step = max(4, min(rows, cols) // 20)
    quiver_step = max(1, int(quiver_step))

    if quiver_length is None:
        quiver_length = 0.8 * quiver_step

    ys = np.arange(0, rows, quiver_step, dtype=int)
    xs = np.arange(0, cols, quiver_step, dtype=int)
    grid_x, grid_y = np.meshgrid(xs, ys)

    vx = vector_field_slice[0][np.ix_(ys, xs)]
    vy = vector_field_slice[1][np.ix_(ys, xs)]
    norms = np.hypot(vx, vy)
    valid = np.isfinite(vx) & np.isfinite(vy) & (norms > 1e-8)
    scalar_values = None
    if scalar_map is not None:
        scalar_values = scalar_map[np.ix_(ys, xs)]
        valid &= np.isfinite(scalar_values)
    if not np.any(valid):
        return

    vx_plot = np.zeros_like(vx, dtype=np.float32)
    vy_plot = np.zeros_like(vy, dtype=np.float32)
    vx_plot[valid] = quiver_length * vx[valid] / norms[valid]
    vy_plot[valid] = quiver_length * vy[valid] / norms[valid]

    if scalar_values is not None:
        vmin = float(np.nanmin(scalar_values[valid]))
        vmax = float(np.nanmax(scalar_values[valid]))
        denom = (vmax - vmin) if vmax != vmin else 1.0
        normalized = np.clip((scalar_values[valid] - vmin) / denom, 0.0, 1.0)
        quiver_colors = (
            colormap(normalized)[:, :3] if colormap is not None else quiver_color
        )
    else:
        quiver_colors = quiver_color

    ax.quiver(
        grid_x[valid],
        grid_y[valid],
        vx_plot[valid],
        vy_plot[valid],
        color=quiver_colors,
        angles="xy",
        scale_units="xy",
        scale=1.0,
        width=0.004,
        headwidth=3.5,
        headlength=5.0,
        headaxislength=4.5,
        pivot="mid",
    )


def plot_images(
    img: np.ndarray,
    img_angle1: np.ndarray,
    img_angle2: np.ndarray,
    img_FA: np.ndarray,
    center_point: tuple[int, int, int],
    vector_field_slice: np.ndarray | None = None,
    colormap_angle=None,
    colormap_angle2=None,
    colormap_FA=None,
    angle1_title: str = "Helix Angle",
    angle2_title: str = "Intrusion Angle",
    angle_ranges: tuple[tuple[float, float], tuple[float, float]] | None = None,
    save_path: str | None = None,
    show: bool = True,
    quiver_step: int | None = None,
    quiver_length: float | None = None,
    quiver_color: str = "cyan",
    overlay_scalar_map: np.ndarray | None = None,
) -> None:
    """
    Render a 2x2 figure of source, angle1, angle2, FA for a single slice.

    Parameters
    ----------
    img : np.ndarray
        2D grayscale slice of the anatomical image.
    img_angle1 : np.ndarray
        2D float map for the first angle, for example HA or AZ, in degrees.
    img_angle2 : np.ndarray
        2D float map for the second angle, for example IA or EL, in degrees.
    img_FA : np.ndarray
        2D float map of fractional anisotropy in [0, 1].
    center_point : tuple[int, int, int]
        Voxel coordinates ``(x, y, z)`` of the centerline point on this slice.
        Only ``(x, y)`` is used here for the marker.
    vector_field_slice : np.ndarray, optional
        Slice of the vector field with shape ``(3, Y, X)``. When provided, a
        subsampled arrow overlay is drawn on the source image.
    colormap_angle : matplotlib colormap, optional
        Colormap for the first angle. Defaults to helix_angle_cmap if None.
    colormap_angle2 : matplotlib colormap, optional
        Colormap for the second angle. Defaults to colormap_angle if None.
    colormap_FA : matplotlib colormap, optional
        Colormap for FA. Defaults to plt.cm.magma if None.
    angle1_title : str
        Title for the first angle panel.
    angle2_title : str
        Title for the second angle panel.
    angle_ranges : tuple[(float, float), (float, float)], optional
        Fixed display range for angle1 and angle2.
    save_path : str, optional
        If provided, save the figure to this location.
    show : bool
        If True, display the figure interactively. Useful to disable in tests.
    quiver_step : int, optional
        Spatial subsampling step for the arrow overlay.
    quiver_length : float, optional
        Arrow length in pixel units for the overlay.
    quiver_color : str
        Arrow color for the vector overlay.
    overlay_scalar_map : np.ndarray, optional
        Optional scalar map used to color the quiver arrows, typically the
        first angle image (HA or AZ).

    Notes
    -----
    This function shows the centerline marker on the source panel and can
    optionally overlay the in-plane vector direction.
    """
    if colormap_angle is None:
        colormap_angle = helix_angle_cmap
    if colormap_angle2 is None:
        colormap_angle2 = colormap_angle
    if colormap_FA is None:
        colormap_FA = plt.cm.magma
    if angle_ranges is None:
        angle_ranges = ((None, None), (None, None))

    if vector_field_slice is not None and vector_field_slice.shape[1:] != img.shape:
        raise ValueError("vector_field_slice shape must match the source image shape")
    if overlay_scalar_map is not None and overlay_scalar_map.shape != img.shape:
        raise ValueError("overlay_scalar_map shape must match the source image shape")

    x, y = center_point[0:2]

    fig, ax = plt.subplots(2, 2, figsize=(10, 9))

    valid = img[np.isfinite(img)]
    vmin, vmax = np.nanpercentile(valid, (2, 98))
    ax[0, 0].imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    ax[0, 0].plot(x, y, "r.", ms=6)
    if vector_field_slice is not None:
        _plot_vector_overlay(
            ax[0, 0],
            vector_field_slice,
            quiver_step=quiver_step,
            quiver_length=quiver_length,
            quiver_color=quiver_color,
            scalar_map=overlay_scalar_map,
            colormap=colormap_angle,
        )
    ax[0, 0].set_title(
        f"Source + vectors ({angle1_title})"
        if vector_field_slice is not None
        else "Source"
    )
    ax[0, 0].axis("off")

    im1 = ax[0, 1].imshow(
        img_angle1,
        cmap=colormap_angle,
        vmin=angle_ranges[0][0],
        vmax=angle_ranges[0][1],
    )
    ax[0, 1].set_title(angle1_title)
    ax[0, 1].axis("off")
    fig.colorbar(
        im1,
        ax=ax[0, 1],
        fraction=0.046,
        pad=0.04,
        ticks=np.linspace(angle_ranges[0][0], angle_ranges[0][1], 10),
    )

    im2 = ax[1, 0].imshow(
        img_angle2,
        cmap=colormap_angle2,
        vmin=angle_ranges[1][0],
        vmax=angle_ranges[1][1],
    )
    ax[1, 0].set_title(angle2_title)
    ax[1, 0].axis("off")
    fig.colorbar(
        im2,
        ax=ax[1, 0],
        fraction=0.046,
        pad=0.04,
        ticks=np.linspace(angle_ranges[1][0], angle_ranges[1][1], 10),
    )

    im3 = ax[1, 1].imshow(img_FA, cmap=colormap_FA, vmin=0.0, vmax=1.0)
    ax[1, 1].set_title("Fractional Anisotropy")
    ax[1, 1].axis("off")
    fig.colorbar(im3, ax=ax[1, 1], fraction=0.046, pad=0.04)

    fig.tight_layout()
    if save_path is not None:
        save_dir = os.path.dirname(save_path)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def write_img_rgb(
    img: np.ndarray,
    out_path: str,
    vmin: float,
    vmax: float,
    colormap: object | None = None,
    output_format: str = "jp2",
) -> None:
    """
    Write a single 2D float image as an RGB file using a matplotlib colormap.

    Parameters
    ----------
    img : np.ndarray
        2D float array to color map.
    out_path : str
        Full output path without extension handling.
    vmin : float
        Minimum for normalization.
    vmax : float
        Maximum for normalization.
    colormap : matplotlib colormap, optional
        Colormap to apply. If None, use helix_angle_cmap.
    output_format : {"jp2", "tif"}
        Output format. Uses glymur for jp2 and tifffile for tif.

    Notes
    -----
    The function normalizes to [0, 1], applies the colormap, then writes uint8 RGB.
    """
    if colormap is None:
        colormap = helix_angle_cmap

    denom = (vmax - vmin) if (vmax - vmin) != 0 else 1.0
    x = (img.astype(np.float32) - vmin) / denom
    x = np.clip(x, 0.0, 1.0)
    rgb = (colormap(x)[..., :3] * 255.0 + 0.5).astype(np.uint8)  # H, W, 3

    if output_format == "jp2":
        write_jp2(out_path, rgb)
    elif output_format == "tif":
        write_tif(out_path, rgb)
    else:
        raise ValueError(f"Unsupported output_format {output_format}")


def write_images(
    img_angle1: np.ndarray,
    img_angle2: np.ndarray,
    img_FA: np.ndarray,
    start_index: int,
    output_dir: str,
    output_format: str,
    output_type: str,
    z: int,
    colormap_angle=None,
    colormap_angle2=None,
    colormap_FA=None,
    angle_names: tuple[str, str] = ("HA", "IA"),
    angle_ranges: tuple[tuple[float, float], tuple[float, float]] = (
        (-90, 90),
        (-90, 90),
    ),
) -> None:
    """
    Write per-slice angle1, angle2, and FA images to disk with flexible naming and ranges.

    Parameters
    ----------
    img_angle1 : np.ndarray
        2D float array for the first angle, for example HA or AZ, in degrees.
    img_angle2 : np.ndarray
        2D float array for the second angle, for example IA or EL, in degrees.
    img_FA : np.ndarray
        2D float array for FA in [0, 1].
    start_index : int
        Global starting index for z numbering in filenames.
    output_dir : str
        Base output directory. Subfolders angle_names[0], angle_names[1], and FA are created.
    output_format : {"jp2", "tif"}
        Output file format. Uses glymur for jp2 and tifffile for tif.
    output_type : {"8bit", "rgb"}
        Write mode. "8bit" writes single channel uint8, "rgb" writes colormapped RGB.
    z : int
        Current z offset used to compute the running index in filenames.
    colormap_angle : matplotlib colormap, optional
        Colormap for the first angle in "rgb" mode. Defaults to helix_angle_cmap if None.
    colormap_angle2 : matplotlib colormap, optional
        Colormap for the second angle in "rgb" mode. Defaults to colormap_angle if None.
    colormap_FA : matplotlib colormap, optional
        Colormap for FA in "rgb" mode. Defaults to plt.cm.magma if None.
    angle_names : tuple[str, str]
        Names used for subfolders and file prefixes, for example ("HA", "IA") or ("AZ", "EL").
    angle_ranges : tuple[(float, float), (float, float)]
        Min and max for normalization of angle1 and angle2. For example
        HA and IA use (-90, 90), AZ uses (-180, 180) and EL uses (0, 90).

    Raises
    ------
    RuntimeError
        If required IO backends are missing for the chosen format.
    ValueError
        If output_format or output_type is unsupported.

    Notes
    -----
    The function creates subdirectories and writes files named like:
    {output_dir}/{name}/{name}_{index:06d}.{ext} and {output_dir}/FA/FA_{index:06d}.{ext}
    """
    name1, name2 = angle_names
    (a1_min, a1_max), (a2_min, a2_max) = angle_ranges
    idx = start_index + z

    os.makedirs(f"{output_dir}/{name1}", exist_ok=True)
    os.makedirs(f"{output_dir}/{name2}", exist_ok=True)
    os.makedirs(f"{output_dir}/FA", exist_ok=True)

    out1 = f"{output_dir}/{name1}/{name1}_{idx:06d}.{output_format}"
    out2 = f"{output_dir}/{name2}/{name2}_{idx:06d}.{output_format}"
    outF = f"{output_dir}/FA/FA_{idx:06d}.{output_format}"

    if output_type == "8bit":
        img1_8 = convert_to_8bit(img_angle1, min_value=a1_min, max_value=a1_max)
        img2_8 = convert_to_8bit(img_angle2, min_value=a2_min, max_value=a2_max)
        imgF_8 = convert_to_8bit(img_FA, min_value=0.0, max_value=1.0)

        if output_format == "jp2":
            write_jp2(out1, img1_8)
            write_jp2(out2, img2_8)
            write_jp2(outF, imgF_8)
        elif output_format == "tif":
            write_tif(out1, img1_8)
            write_tif(out2, img2_8)
            write_tif(outF, imgF_8)
        else:
            raise ValueError(f"I do not recognise output_format {output_format}")

    elif output_type == "rgb":
        if colormap_angle is None:
            colormap_angle = helix_angle_cmap
        if colormap_angle2 is None:
            colormap_angle2 = colormap_angle
        if colormap_FA is None:
            colormap_FA = plt.cm.magma

        write_img_rgb(img_angle1, out1, a1_min, a1_max, colormap_angle, output_format)
        write_img_rgb(img_angle2, out2, a2_min, a2_max, colormap_angle2, output_format)
        write_img_rgb(img_FA, outF, 0.0, 1.0, colormap_FA, output_format)

    else:
        raise ValueError(f"I do not recognise output_type {output_type}")
