import os
import sys
import warnings

import glymur
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from scipy.interpolate import CubicSpline

# Optional GPU support
try:
    import cupy as cp  # noqa: F401

    USE_GPU = True
    print("GPU support enabled.")
except ImportError:
    USE_GPU = False

print(f"USE_GPU: {USE_GPU}")

from structure_tensor.multiprocessing import parallel_structure_tensor_analysis

from cardiotensor.utils.utils import convert_to_8bit


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

    # Extract the Dominant Direction
    center_vec = -center_vec[[2, 1, 0]]

    return center_vec


def adjust_start_end_index(
    start_index: int,
    end_index: int,
    N_img: int,
    padding_start: int,
    padding_end: int,
    is_test: bool,
    n_slice: int,
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
    SIGMA: float,
    RHO: float,
    devices: list[str] | None = None,
    block_size: int = 200,
    use_gpu: bool = False,
    dtype: type = np.float32,  # Default to np.float64
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculates the structure tensor of a volume.

    Args:
        volume (np.ndarray): The 3D volume data.
        SIGMA (float): Sigma value for Gaussian smoothing.
        RHO (float): Rho value for Gaussian smoothing.
        devices (Optional[list[str]]): List of devices for parallel processing (e.g., ['cpu', 'cuda:0']).
        block_size (int): Size of the blocks for processing. Default is 200.
        use_gpu (bool): If True, uses GPU for calculations. Default is False.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: Structure tensor, eigenvalues, and eigenvectors.
    """
    # Filter or ignore specific warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    num_cpus = os.cpu_count() or 4  # Default to 4 if os.cpu_count() returns None
    num_cpus = max(num_cpus, 4)
    print(f"Number of CPUs used: {num_cpus}")

    if devices is None:  # Initialize devices if not provided
        devices = []

    if use_gpu:
        print("GPU activated")
        if not devices:  # Assign default GPU and CPU devices if the list is empty
            devices = 16 * ["cuda:0"] + 16 * ["cuda:1"] + num_cpus * ["cpu"]

        S, val, vec = parallel_structure_tensor_analysis(
            volume,
            SIGMA,
            RHO,
            devices=devices,
            block_size=block_size,
            truncate=4.0,
            structure_tensor=None,
            eigenvectors=dtype,
            eigenvalues=dtype,
        )
    else:
        print("GPU not activated")
        S, val, vec = parallel_structure_tensor_analysis(
            volume,
            SIGMA,
            RHO,
            devices=num_cpus * ["cpu"],
            block_size=block_size,
            truncate=4.0,
            structure_tensor=None,
            eigenvectors=dtype,
            eigenvalues=dtype,
        )  # vec has shape =(3,x,y,z) in the order of (z,y,x)

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
    denominator = np.sqrt(l1**2 + l2**2 + l3**2)
    img_FA = np.sqrt(3 / 2) * (numerator / denominator)

    return img_FA


def rotate_vectors_to_new_axis(
    vector_field_slice: np.ndarray, new_axis_vec: np.ndarray
) -> np.ndarray:
    """
    Rotates a vector field slice to align with a new axis.

    Args:
        vector_field_slice (np.ndarray): Array of vectors to rotate.
        new_axis_vec (np.ndarray): The new axis to align vectors with.

    Returns:
        np.ndarray: Rotated vectors aligned with the new axis.
    """
    # Ensure new_axis_vec is normalized
    new_axis_vec = new_axis_vec / np.linalg.norm(new_axis_vec)

    # Calculate the rotation matrix
    vec1 = np.array([1, 0, 0])  # Initial vertical axis

    vec1 = vec1 * np.sign(new_axis_vec[0])

    a = (vec1 / np.linalg.norm(vec1)).reshape(3)
    b = (new_axis_vec).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)

    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    if np.any(kmat):
        rotation_matrix = (
            np.eye(3) + kmat + np.dot(kmat, kmat) * ((1 - c) / (np.linalg.norm(v) ** 2))
        )
    else:
        rotation_matrix = np.eye(3)

    # Reshape vec_2D to (3, N) for matrix multiplication
    vec_2D_reshaped = np.reshape(vector_field_slice, (3, -1))

    vec_2D_reshaped = vec_2D_reshaped / np.linalg.norm(vec_2D_reshaped, axis=0)

    # Rotate the vectors
    rotated_vecs = np.dot(rotation_matrix, vec_2D_reshaped)

    # Reshape back to the original shape
    rotated_vecs = rotated_vecs.reshape(vector_field_slice.shape)

    # print(f"Rotation matrix:\n{rotation_matrix}")

    return rotated_vecs


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

    # Calculate helix and transverse angles
    helix_angle = np.arctan(
        reshaped_rotated_vector_field[2, :, :] / reshaped_rotated_vector_field[1, :, :]
    )
    transverse_angle = np.arctan(
        reshaped_rotated_vector_field[0, :, :] / reshaped_rotated_vector_field[1, :, :]
    )
    helix_angle = np.rad2deg(helix_angle)
    transverse_angle = np.rad2deg(transverse_angle)

    return helix_angle, transverse_angle


def plot_images(
    img: np.ndarray,
    img_helix: np.ndarray,
    img_intrusion: np.ndarray,
    img_FA: np.ndarray,
    center_point: tuple[int, int, int],
) -> None:
    """
    Plots images of the heart with helix, intrusion, and FA annotations.

    Args:
        img (np.ndarray): Grayscale image of the heart.
        img_helix (np.ndarray): Helix angle image.
        img_intrusion (np.ndarray): Intrusion angle image.
        img_FA (np.ndarray): Fractional Anisotropy (FA) image.
        center_point (Tuple[int, int, int]): Coordinates of the center point.

    Returns:
        None
    """

    img_vmin, img_vmax = np.nanpercentile(img, (5, 95))
    orig_map = plt.get_cmap("hsv")

    # Create a figure and axes
    fig, axes = plt.subplots(
        2, 2, figsize=(10, 8)
    )  # Adjust the size for better visibility
    ax = axes

    # Original Image with Red Point
    ax[0, 0].imshow(img, vmin=img_vmin, vmax=img_vmax, cmap=plt.cm.gray)
    x, y = center_point[0:2]
    ax[0, 0].scatter(x, y, c="red", s=50, marker="o", label="Axis Point")
    ax[0, 0].set_title("Original Image")
    ax[0, 0].legend(loc="upper right")

    # Helix Image
    tmp = ax[0, 1].imshow(img_helix, cmap=orig_map)
    ax[0, 1].set_title("Helix Angle")

    # Intrusion Image
    ax[1, 0].imshow(img_intrusion, cmap=orig_map)
    ax[1, 0].set_title("Intrusion Angle")

    # FA Image
    fa_plot = ax[1, 1].imshow(img_FA, cmap="inferno")
    ax[1, 1].set_title("Fractional Anisotropy")

    # Add colorbars for relevant subplots
    cbar1 = fig.colorbar(tmp, ax=ax[0, 1], orientation="vertical")
    cbar1.set_label("Helix Angle")
    cbar2 = fig.colorbar(fa_plot, ax=ax[1, 1], orientation="vertical")
    cbar2.set_label("Fractional Anisotropy")

    # Set common labels
    for axis in ax.flat:
        axis.axis("off")  # Turn off axes if desired

    # Adjust layout to prevent overlap
    fig.tight_layout()

    plt.show()


def write_images(
    img_helix: np.ndarray,
    img_intrusion: np.ndarray,
    img_FA: np.ndarray,
    start_index: int,
    OUTPUT_DIR: str,
    OUTPUT_FORMAT: str,
    OUTPUT_TYPE: str,
    z: int,
) -> None:
    """
    Writes processed images to the specified directory.

    Args:
        img_helix (np.ndarray): Image data for helix angles.
        img_intrusion (np.ndarray): Image data for intrusion angles.
        img_FA (np.ndarray): Image data for fractional anisotropy.
        start_index (int): Starting index for filenames.
        OUTPUT_DIR (str): Directory to save the images.
        OUTPUT_FORMAT (str): Format of the output files ('tif' or 'jp2').
        OUTPUT_TYPE (str): Type of output ('8bit' or 'rgb').
        z (int): Current slice index.

    Returns:
        None
    """

    try:
        os.makedirs(OUTPUT_DIR + "/HA", exist_ok=True)
        os.makedirs(OUTPUT_DIR + "/IA", exist_ok=True)
        os.makedirs(OUTPUT_DIR + "/FA", exist_ok=True)
    except PermissionError:
        print("⚠️ - Permission error during creation of output directories")

    # print(f"Saving image: {z}")

    if "8bit" in OUTPUT_TYPE:
        # Convert the float64 image to int8
        img_helix = convert_to_8bit(img_helix, min_value=-90, max_value=90)
        img_intrusion = convert_to_8bit(img_intrusion, min_value=-90, max_value=90)
        img_FA = convert_to_8bit(img_FA, min_value=0, max_value=1)

        if OUTPUT_FORMAT == "jp2":
            ratio_compression = 10

            # Define file paths
            ha_path = f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.jp2"
            ia_path = f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.jp2"
            fa_path = f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.jp2"

            # Remove existing files if they exist
            for file_path in [ha_path, ia_path, fa_path]:
                if os.path.exists(file_path):
                    os.remove(file_path)

            # Write new JP2 files
            glymur.Jp2k(
                ha_path,
                data=img_helix,
                cratios=[ratio_compression],
                numres=8,
                irreversible=True,
            )
            glymur.Jp2k(
                ia_path,
                data=img_intrusion,
                cratios=[ratio_compression],
                numres=8,
                irreversible=True,
            )
            glymur.Jp2k(
                fa_path,
                data=img_FA,
                cratios=[ratio_compression],
                numres=8,
                irreversible=True,
            )
        elif OUTPUT_FORMAT == "tif":
            tifffile.imwrite(
                f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.tif", img_helix
            )
            tifffile.imwrite(
                f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.tif", img_intrusion
            )
            tifffile.imwrite(f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.tif", img_FA)
        else:
            sys.exit(f"I don't recognise the OUTPUT_FORMAT ({OUTPUT_FORMAT})")

    elif "rgb" in OUTPUT_TYPE:

        def write_img_rgb(
            img: np.ndarray,
            output_path: str,
            cmap: plt.Colormap | None = plt.get_cmap("hsv"),
        ) -> None:
            """
            Writes an RGB image to the specified output path.

            Args:
                img (np.ndarray): The input image data to be converted and saved.
                output_path (str): The path where the output image will be saved.
                cmap (Optional[plt.Colormap]): The colormap to use for converting the image.
                                            Default is the 'hsv' colormap.

            Returns:
                None
            """
            minimum = np.nanmin(img)
            maximum = np.nanmax(img)
            img = (img + np.abs(minimum)) * (1 / (maximum - minimum))

            if cmap is not None:
                img = cmap(img)
            img = (img[:, :, :3] * 255).astype(np.uint8)

            print(f"Writing image to {output_path}")
            if OUTPUT_FORMAT == "jp2":
                glymur.Jp2k(
                    output_path,
                    data=img,
                    cratios=[ratio_compression],
                    numres=8,
                    irreversible=True,
                )
            elif OUTPUT_FORMAT == "tif":
                tifffile.imwrite(output_path, img)
            else:
                sys.exit(f"I don't recognise the OUTPUT_FORMAT ({OUTPUT_FORMAT})")

        if OUTPUT_FORMAT == "jp2":
            write_img_rgb(
                img_helix,
                f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.jp2",
                cmap=plt.get_cmap("hsv"),
            )
            write_img_rgb(
                img_intrusion,
                f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.jp2",
                cmap=plt.get_cmap("hsv"),
            )
            write_img_rgb(
                img_FA,
                f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.jp2",
                cmap=plt.get_cmap("inferno"),
            )
        elif OUTPUT_FORMAT == "tif":
            write_img_rgb(
                img_helix,
                f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.tif",
                cmap=plt.get_cmap("hsv"),
            )
            write_img_rgb(
                img_intrusion,
                f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.tif",
                cmap=plt.get_cmap("hsv"),
            )
            write_img_rgb(
                img_FA,
                f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.tif",
                cmap=plt.get_cmap("inferno"),
            )
        else:
            sys.exit(f"I don't recognise the OUTPUT_FORMAT ({OUTPUT_FORMAT})")


def write_vector_field(
    vector_field_slice: np.ndarray, start_index: int, output_dir: str, slice_idx: int
) -> None:
    """
    Saves a vector field slice to the specified directory in .npy format.

    Args:
        vector_field_slice (np.ndarray): Vector field data slice.
        start_index (int): Starting index for filenames.
        output_dir (str): Directory to save the vector field.
        slice_idx (int): Current slice index.

    Returns:
        None
    """
    os.makedirs(f"{output_dir}/eigen_vec", exist_ok=True)
    np.save(
        f"{output_dir}/eigen_vec/eigen_vec_{(start_index + slice_idx):06d}.npy",
        vector_field_slice,
    )
    # print(f"Vector field slice saved at index {slice_idx}")
