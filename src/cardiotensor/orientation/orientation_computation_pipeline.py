import math
import multiprocessing as mp
import os
import sys
import time

import cv2
import numpy as np
from alive_progress import alive_bar

# from memory_profiler import profile
from cardiotensor.orientation.orientation_computation_functions import (
    adjust_start_end_index,
    calculate_center_vector,
    calculate_structure_tensor,
    compute_fraction_anisotropy,
    compute_helix_and_transverse_angles,
    interpolate_points,
    plot_images,
    remove_padding,
    rotate_vectors_to_new_axis,
    write_images,
    write_vector_field,
)
from cardiotensor.utils.utils import read_conf_file
from cardiotensor.utils.DataReader import DataReader


MULTIPROCESS = True


def is_tiff_image_valid(image_path: str) -> bool:
    """
    Check if the TIFF image is readable and valid.

    Args:
        image_path (str): Path to the TIFF image.

    Returns:
        bool: True if the image is valid, False otherwise.
    """

    try:
        image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
        if image is None:
            raise ValueError(
                "cv2.imread failed to load the image. Unsupported or corrupted TIFF file."
            )
        return True
    except Exception as e:
        print(f"Error reading TIFF image {image_path}: {e}")
        return False






# @profile
def compute_orientation(
    conf_file_path: str, start_index: int = 0, end_index: int | None = None, use_gpu: bool = False
) -> None:
    """
    Compute the orientation for a volume dataset based on the configuration.

    Args:
        conf_file_path (str): Path to the configuration file.
        start_index (int, optional): Start index for processing. Default is 0.
        end_index (int, optional): End index for processing. Default is None.
        use_gpu (bool, optional): Whether to use GPU for calculations. Default is False.

    Returns:
        None
    """

    print("\n---------------------------------")
    print(f"🤖 - Processing slices {start_index} to {end_index}")

    print("\n---------------------------------")
    print(f"READING PARAMETER FILE : {conf_file_path}\n")

    print(f"Start index, End index : {start_index}, {end_index}\n")

    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(e)
        sys.exit(f"⚠️  Error reading parameter file: {conf_file_path}")

    (
        VOLUME_PATH,
        MASK_PATH,
        IS_FLIP,
        OUTPUT_DIR,
        OUTPUT_FORMAT,
        OUTPUT_TYPE,
        IS_VECTORS,
        SIGMA,
        RHO,
        N_CHUNK,
        PT_MV,
        PT_APEX,
        IS_TEST,
        N_SLICE_TEST,
    ) = (
        params[key]
        for key in [
            "IMAGES_PATH",
            "MASK_PATH",
            "FLIP",
            "OUTPUT_PATH",
            "OUTPUT_FORMAT",
            "OUTPUT_TYPE",
            "VECTORS",
            "SIGMA",
            "RHO",
            "N_CHUNK",
            "POINT_MITRAL_VALVE",
            "POINT_APEX",
            "TEST",
            "N_SLICE_TEST",
        ]
    )

    #Check if alreayd done the processing
    if not IS_TEST:
        is_already_done = True
        if end_index is None:
            is_already_done = False
        for idx in range(start_index, end_index):
            if (
                not os.path.exists(f"{OUTPUT_DIR}/HA/HA_{idx:06d}.tif")
                or not os.path.exists(f"{OUTPUT_DIR}/IA/IA_{idx:06d}.tif")
                or not os.path.exists(f"{OUTPUT_DIR}/FA/FA_{idx:06d}.tif")
                or not os.path.exists(f"{OUTPUT_DIR}/eigen_vec/eigen_vec_{idx:06d}.npy")
            ):
                is_already_done = False
                break

        if is_already_done:
            print("All images are already done")
            return

    if MASK_PATH:
        is_mask = True
    else:
        print("Mask path not provided")
        is_mask = False

    print("\n---------------------------------")
    print("READING VOLUME INFORMATION\n")
    print(f"Volume path: {VOLUME_PATH}")
    
    data_reader = DataReader(VOLUME_PATH)
    
    print(f"Number of slices: {data_reader.shape[0]}")

    print("\n---------------------------------")
    print("CALCULATE CENTER LINE AND CENTER VECTOR\n")
    center_line = interpolate_points(PT_MV, PT_APEX, data_reader.shape[0])
    center_vec = calculate_center_vector(PT_MV, PT_APEX, IS_FLIP)
    print(f"Center vector: {center_vec}")

    print("\n---------------------------------")
    print("CALCULATE PADDING START AND ENDING INDEXES\n")
    padding_start = math.ceil(RHO)
    padding_end = math.ceil(RHO)
    if not IS_TEST:
        if padding_start > start_index:
            padding_start = start_index
        if padding_end > (data_reader.shape[0] - end_index):
            padding_end = data_reader.shape[0] - end_index
    if IS_TEST:
        if N_SLICE_TEST > data_reader.shape[0]:
            sys.exit("Error: N_SLICE_TEST > number of images")

    print(f"Padding start, Padding end : {padding_start}, {padding_end}")
    start_index_padded, end_index_padded = adjust_start_end_index(
        start_index, end_index, data_reader.shape[0], padding_start, padding_end, IS_TEST, N_SLICE_TEST
    )
    print(
        f"Start index padded, End index padded : {start_index_padded}, {end_index_padded}"
    )

    print("\n---------------------------------")
    print("LOAD DATASET\n")
    volume = data_reader.load_volume(start_index_padded, end_index_padded).astype(
        "float32"
    )
    print(f"Loaded volume shape {volume.shape}")

    if is_mask:
        mask_reader = DataReader(MASK_PATH)
        
        mask = mask_reader.load_volume(start_index_padded, end_index_padded, unbinned_shape=data_reader.shape).astype(
            "float32"
        )
        
        assert (
            mask.shape == volume.shape
        ), f"Mask shape {mask.shape} does not match volume shape {volume.shape}"
        
        volume[mask == 0] = 0

        
    print("\n---------------------------------")
    print("CALCULATING STRUCTURE TENSOR")
    t1 = time.perf_counter()  # start time
    val, vec = calculate_structure_tensor(volume, SIGMA, RHO, use_gpu=use_gpu)
    print(f"Vector shape: {vec.shape}")

    if is_mask:
        volume[mask == 0] = np.nan
        val[0, :, :, :][mask == 0] = np.nan
        val[1, :, :, :][mask == 0] = np.nan
        val[2, :, :, :][mask == 0] = np.nan
        vec[0, :, :, :][mask == 0] = np.nan
        vec[1, :, :, :][mask == 0] = np.nan
        vec[2, :, :, :][mask == 0] = np.nan

        print("Mask applied to image volume")

        del mask

    volume, val, vec = remove_padding(volume, val, vec, padding_start, padding_end)
    print(f"Vector shape after removing padding: {vec.shape}")

    center_line = center_line[start_index_padded:end_index_padded]

    # Putting all the vectors in positive direction
    # posdef = np.all(val >= 0, axis=0)  # Check if all elements are non-negative along the first axis
    vec = vec / np.linalg.norm(vec, axis=0)

    # Check for negative z component and flip if necessary
    # negative_z = vec[2, :] < 0
    # vec[:, negative_z] *= -1

    t2 = time.perf_counter()  # stop time
    print(f"finished calculating structure tensors in {t2 - t1} seconds")

    print("\nCalculating helix/intrusion angle and fractional anisotropy:")
    if MULTIPROCESS and not IS_TEST:
        print(f"Using {mp.cpu_count()} CPU cores")
        with mp.Pool(processes=mp.cpu_count()) as pool:
            # Initialize the alive-progress bar
            with alive_bar(
                vec.shape[1], title="Processing slices (Multiprocess)", bar="smooth"
            ) as bar:
                result = pool.starmap_async(
                    compute_slice_angles_and_anisotropy,
                    [
                        (
                            z,
                            vec[:, z, :, :],
                            volume[z, :, :],
                            np.around(center_line[z]),
                            val[:, z, :, :],
                            center_vec,
                            OUTPUT_DIR,
                            OUTPUT_FORMAT,
                            OUTPUT_TYPE,
                            start_index,
                            IS_VECTORS,
                            IS_TEST,
                        )
                        for z in range(vec.shape[1])
                    ],
                    callback=lambda _: bar(
                        vec.shape[1]
                    ),  # Update progress bar once all tasks are completed
                )
                result.wait()  # Wait for all tasks to complete

    else:
        # Add a progress bar for single-threaded processing
        with alive_bar(
            vec.shape[1], title="Processing slices (Single-thread)", bar="smooth"
        ) as bar:
            for z in range(vec.shape[1]):
                # Call the function directly
                compute_slice_angles_and_anisotropy(
                    z,
                    vec[:, z, :, :],
                    volume[z, :, :],
                    np.around(center_line[z]),
                    val[:, z, :, :],
                    center_vec,
                    OUTPUT_DIR,
                    OUTPUT_FORMAT,
                    OUTPUT_TYPE,
                    start_index,
                    IS_VECTORS,
                    IS_TEST,
                )
                bar()  # Update the progress bar for each slice

    # #Check images
    # for idx in range(start_index, end_index):
    #     ia_path = f"{OUTPUT_DIR}/IA/IA_{idx:06d}.tif"
    #     ha_path = f"{OUTPUT_DIR}/HA/HA_{idx:06d}.tif"
    #     fa_path = f"{OUTPUT_DIR}/FA/FA_{idx:06d}.tif"

    #     # Validate IA image
    #     if not is_tiff_image_valid(ia_path):
    #         print(f"Invalid or corrupt IA image: {ia_path}")
    #         if os.path.exists(ia_path):
    #             os.remove(ia_path)
    #         continue

    #     # Validate HA image
    #     if not is_tiff_image_valid(ha_path):
    #         print(f"Invalid or corrupt HA image: {ha_path}")
    #         if os.path.exists(ha_path):
    #             os.remove(ha_path)
    #         continue

    #     # Validate FA image
    #     if not is_tiff_image_valid(fa_path):
    #         print(f"Invalid or corrupt FA image: {fa_path}")
    #         if os.path.exists(fa_path):
    #             os.remove(fa_path)
    #         continue

    #     print(f"Valid image set: {idx:06d}")

    print(f"\n🤖 - Finished processing slices {start_index} - {end_index}")
    print("---------------------------------\n\n")

    return


def compute_slice_angles_and_anisotropy(
    z: int,
    vector_field_slice: np.ndarray,
    img_slice: np.ndarray,
    center_point: np.ndarray,
    eigen_val_slice: np.ndarray,
    center_vec: np.ndarray,
    OUTPUT_DIR: str,
    OUTPUT_FORMAT: str,
    OUTPUT_TYPE: str,
    start_index: int,
    IS_VECTORS: bool,
    IS_TEST: bool,
) -> None:
    """
    Compute helix angles, transverse angles, and fractional anisotropy for a slice.

    Args:
        z (int): Index of the slice.
        vector_field_slice (np.ndarray): Vector field for the slice.
        img_slice (np.ndarray): Image data for the slice.
        center_point (np.ndarray): Center point for alignment.
        eigen_val_slice (np.ndarray): Eigenvalues for the slice.
        center_vec (np.ndarray): Center vector for alignment.
        OUTPUT_DIR (str): Directory to save the output.
        OUTPUT_FORMAT (str): Format for the output files (e.g., "tif").
        OUTPUT_TYPE (str): Type of output (e.g., "HA", "IA").
        start_index (int): Start index of the slice.
        IS_VECTORS (bool): Whether to output vector fields.
        IS_TEST (bool): Whether in test mode.

    Returns:
        None
    """
    # print(f"Processing image: {start_index + z}")
    paths = [
        f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.tif",
        f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.tif",
        f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.tif",
    ]
    if not IS_TEST and all(os.path.exists(path) for path in paths):
        # print(f"File {(start_index + z):06d} already exists")
        return

    img_FA = compute_fraction_anisotropy(eigen_val_slice)
    vector_field_slice_rotated = rotate_vectors_to_new_axis(
        vector_field_slice, center_vec
    )
    img_helix, img_intrusion = compute_helix_and_transverse_angles(
        vector_field_slice_rotated, center_point
    )

    if IS_TEST:
        plot_images(img_slice, img_helix, img_intrusion, img_FA, center_point)
    else:
        write_images(
            img_helix,
            img_intrusion,
            img_FA,
            start_index,
            OUTPUT_DIR,
            OUTPUT_FORMAT,
            OUTPUT_TYPE,
            z,
        )
        if IS_VECTORS:
            write_vector_field(vector_field_slice, start_index, OUTPUT_DIR, z)
