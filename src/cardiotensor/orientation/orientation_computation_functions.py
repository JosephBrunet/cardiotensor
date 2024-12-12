import os
import sys
import warnings

import cv2
import glymur
import matplotlib.pyplot as plt
import numpy as np
import tifffile

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


def interpolate_points(point1, point2, N_img):
    """
    Generates interpolated points between two points.

    Parameters:
    - point1: tuple
        The first point (x, y, z).
    - point2: tuple
        The second point (x, y, z).
    - N_img: int
        The number of images or points to interpolate.

    Returns:
    - np.array
        Array of interpolated points.
    """
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    z_values = list(range(N_img))

    points = []
    for z in z_values:
        t = (z - z1) / (z2 - z1)  # Calculate the interpolation parameter
        x = x1 + (x2 - x1) * t
        y = y1 + (y2 - y1) * t
        points.append((x, y, z))
    return np.array(points)


def calculate_center_vector(pt_mv, pt_apex, is_flip=True):
    """
    Calculates the center vector between two points.

    Parameters:
    - pt_mv: numpy.ndarray
        The mitral valve point.
    - pt_apex: numpy.ndarray
        The apex point.

    Returns:
    - center_vec: numpy.ndarray
        The normalized center vector.
    """
    # Calculate center vector

    center_vec = pt_mv - pt_apex
    center_vec = center_vec / np.linalg.norm(center_vec)

    if is_flip:
        center_vec[2] = -center_vec[
            2
        ]  # Invert z components (Don't know why but it works)

    # center_vec[0] = -center_vec[0]  # Invert z components (Don't know why but it works)
    # center_vec[1] = -center_vec[1]  # Invert z components (Don't know why but it works)

    # if center_vec[2] < 0:
    #     center_vec = -center_vec

    return center_vec


def adjust_start_end_index(
    start_index, end_index, N_img, padding_start, padding_end, is_test, n_slice
):
    """
    Adjusts the start and end indices for processing.

    Parameters:
    - start_index: int
        The initial start index.
    - end_index: int
        The initial end index.
    - N_img: int
        Number of images in the volume data.
    - padding_start: int
        Padding to add at the start.
    - padding_end: int
        Padding to add at the end.
    - is_test: bool
        Flag to indicate if it is a test run.

    Returns:
    - start_index_padded: int
        The adjusted start index.
    - end_index_padded: int
        The adjusted end index.
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
    volume, SIGMA, RHO, device="", block_size="", use_gpu=False
):
    """
    Calculates the structure tensor for the given volume data.

    Parameters:
    - volume: numpy.ndarray
        The volume data.
    - SIGMA: float
        Sigma value for structure tensor calculation.
    - RHO: float
        Rho value for structure tensor calculation.
    - use_gpu: bool
        Flag to use GPU if available.

    Returns:
    - tuple: (s, val, vec)
        The structure tensor (s), eigenvectors (vec), and eigenvalues (val).
    """
    # Filter or ignore specific warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    num_cpus = os.cpu_count()
    num_cpus = max(num_cpus, 4)
    print(f"Number of CPUs used: {num_cpus}")

    if use_gpu:
        print("GPU activated")
        if device == "":
            device = 16 * ["cuda:0"] + 16 * ["cuda:1"] + num_cpus * ["cpu"]
        if block_size == "":
            block_size = 400

        S, val, vec = parallel_structure_tensor_analysis(
            volume,
            SIGMA,
            RHO,
            devices=device,
            block_size=block_size,
            truncate=4.0,
            structure_tensor=np.float32,
        )
    else:
        print("GPU not activated")
        S, val, vec = parallel_structure_tensor_analysis(
            volume,
            SIGMA,
            RHO,
            devices=num_cpus * ["cpu"],
            truncate=4.0,
            structure_tensor=np.float32,
        )  # vec has shape =(3,x,y,z) in the order of (z,y,x)

    return S, val, vec


def remove_padding(volume, s, val, vec, padding_start, padding_end):
    """
    Removes padding from the processed data.

    Parameters:
    - volume: numpy.ndarray
        The volume data.
    - s: numpy.ndarray
        The structure tensor data.
    - vec: numpy.ndarray
        The eigenvectors.
    - val: numpy.ndarray
        The eigenvalues.
    - padding_start: int
        The start padding to remove.
    - padding_end: int
        The end padding to remove.

    Returns:
    - tuple: (volume, s, val, vec)
        The adjusted data without padding.
    """
    array_end = vec.shape[1] - padding_end
    volume = volume[padding_start:array_end, :, :]
    s = s[:, padding_start:array_end, :, :]
    vec = vec[:, padding_start:array_end, :, :]
    val = val[:, padding_start:array_end, :, :]

    return volume, s, val, vec


def compute_fraction_anisotropy(eigenvalues_2d):
    """
    Calculates Fractional Anisotropy from eigenvalues of a structure tensor.

    Parameters:
    - eigenvalues_2d (numpy.ndarray): Eigenvalues of a slice of the 3D volume (l1, l2, l3).

    Returns:
    - numpy.ndarray: An array of Fractional Anisotropy values.
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


def rotate_vectors_to_new_axis(vector_field_slice, new_axis_vec):
    """
    Rotates vectors to align with a new axis.

    Parameters:
    - vector_field_slice (numpy.ndarray): Array of vectors to be rotated.
    - new_axis_vec (numpy.ndarray): The new axis to align the vectors with.

    Returns:
    - numpy.ndarray: Rotated vectors aligned with the new axis.
    """
    # Ensure new_axis_vec is normalized
    new_axis_vec = new_axis_vec / np.linalg.norm(new_axis_vec)

    # Calculate the rotation matrix
    vec1 = np.array([0, 0, 1])  # Initial vertical axis

    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (new_axis_vec).reshape(3)
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


def compute_helix_and_transverse_angles(vector_field_2d, center_point):
    """
    Calculate the helix angles at each point of a 2D image using a vector field.

    Parameters:
    - image (numpy.ndarray): The 2D image array.
    - vector_field (numpy.ndarray): The 2D orientation vector field.
    - center (tuple): The center point (x, y) of the image.

    Returns:
    - tuple of numpy.ndarray: Arrays containing the helix and transverse angles.
    """
    center = center_point[0:2]  # Replace with actual values
    rows, cols = vector_field_2d.shape[1:3]

    reshaped_vector_field = np.reshape(vector_field_2d, (3, -1))

    center_x, center_y = center[0], center[1]

    X, Y = np.meshgrid(np.arange(cols) - center_x, np.arange(rows) - center_y)

    theta = -np.arctan2(Y.flatten(), X.flatten())
    cos_angle = np.cos(theta)
    sin_angle = np.sin(theta)

    # Rotate the vector field
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


def plot_images(img, img_helix, img_intrusion, img_FA, center_point, PT_MV, PT_APEX):
    """
    Plot images of the heart.

    Args:
        img (numpy.ndarray): The grayscale image of the heart.
        img_helix (numpy.ndarray): The helix image of the heart.
        img_intrusion (numpy.ndarray): The intrusion image of the heart.
        img_FA (numpy.ndarray): The FA (fractional anisotropy) image of the heart.
        center_point (tuple): The coordinates of the center point.
        PT_MV (tuple): The coordinates of the mitral valve point.
        PT_APEX (tuple): The coordinates of the apex point.

    Returns:
        None
    """

    print("\nPlotting images...")
    img_vmin, img_vmax = np.nanpercentile(img, (5, 95))
    orig_map = plt.get_cmap("hsv")
    # reversed_map = orig_map.reversed()
    fig, axes = plt.subplots(2, 2, figsize=(8, 4))
    ax = axes
    ax[0, 0].imshow(img, vmin=img_vmin, vmax=img_vmax, cmap=plt.cm.gray)
    # Draw a red point at the specified coordinate
    x, y = center_point[0:2]
    ax[0, 0].scatter(x, y, c="red", s=50, marker="o")
    # ax[0,0].scatter(PT_MV[0],PT_MV[1], c='green', s=50, marker='o')
    # ax[0,0].scatter(PT_APEX[0],PT_APEX[1], c='blue', s=50, marker='o')
    tmp = ax[0, 1].imshow(img_helix, cmap=orig_map)
    ax[1, 0].imshow(img_intrusion, cmap=orig_map)
    fig.colorbar(tmp)
    plt.show()


def write_images(
    img_helix,
    img_intrusion,
    img_FA,
    start_index,
    OUTPUT_DIR,
    OUTPUT_FORMAT,
    OUTPUT_TYPE,
    z,
):
    """
    Write images to the specified output directory.

    Args:
        img_helix (numpy.ndarray): The image data for helix.
        img_intrusion (numpy.ndarray): The image data for intrusion.
        img_FA (numpy.ndarray): The image data for FA (False Alarm).
        start_index (int): The starting index for the image filenames.
        OUTPUT_DIR (str): The output directory to save the images.

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
        img_helix = convert_to_8bit(img_helix, output_min=-90, output_max=90)
        img_intrusion = convert_to_8bit(img_intrusion, output_min=-90, output_max=90)
        img_FA = convert_to_8bit(img_FA, output_min=0, output_max=1)

        if OUTPUT_FORMAT == "jp2":
            ratio_compression = 10
            glymur.Jp2k(
                f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.jp2",
                data=img_helix,
                cratios=[ratio_compression],
                numres=8,
                irreversible=True,
            )
            glymur.Jp2k(
                f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.jp2",
                data=img_intrusion,
                cratios=[ratio_compression],
                numres=8,
                irreversible=True,
            )
            glymur.Jp2k(
                f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.jp2",
                data=img_FA,
                cratios=[ratio_compression],
                numres=8,
                irreversible=True,
            )
        elif OUTPUT_FORMAT == "tif":
            cv2.imwrite(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.tif", img_helix)
            cv2.imwrite(
                f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.tif", img_intrusion
            )
            cv2.imwrite(f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.tif", img_FA)
        else:
            sys.exit(f"I don't recognise the OUTPUT_FORMAT ({OUTPUT_FORMAT})")

    elif "rgb" in OUTPUT_TYPE:

        def write_img_rgb(img, output_path, cmap=plt.get_cmap("hsv")):
            minimum = np.nanmin(img)
            maximum = np.nanmax(img)
            img = (img + np.abs(minimum)) * (1 / (maximum - minimum))
            img = cmap(img)
            img = (img[:, :, :3] * 255).astype(np.uint8)

            # cv2.imwrite(output_path, img)
            print(f"Writing image to {output_path}")
            if OUTPUT_FORMAT == "jp2":
                cv2.imwrite(output_path, img)
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


def write_vector_field(vector_field_slice, start_index, output_dir, slice_idx):
    """
    Save the vector field slice to the specified directory in .npy format.

    Args:
        vector_field_slice (numpy.ndarray): The vector field data slice.
        start_index (int): Starting index for the slice filenames.
        output_dir (str): Output directory to save the vector field slice.
        slice_idx (int): Current slice index being saved.

    Returns:
        None
    """
    os.makedirs(f"{output_dir}/eigen_vec", exist_ok=True)
    np.save(
        f"{output_dir}/eigen_vec/eigen_vec_{(start_index + slice_idx):06d}.npy",
        vector_field_slice,
    )
    # print(f"Vector field slice saved at index {slice_idx}")
