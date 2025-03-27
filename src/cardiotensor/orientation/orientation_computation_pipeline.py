import math
import multiprocessing as mp
import os
import sys
import time

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
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import read_conf_file, remove_corrupted_files

MULTIPROCESS = True


# @profile
def compute_orientation(
    conf_file_path: str,
    start_index: int = 0,
    end_index: int | None = None,
    use_gpu: bool = False,
) -> None:
    """
    Compute the orientation for a volume dataset based on the configuration.

    Args:
        conf_file_path (str): Path to the configuration file.
        start_index (int, optional): Start index for processing. Default is 0.
        end_index (int | None): End index for processing. Default is None.
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
        print(f"⚠️  Error reading parameter file '{conf_file_path}': {e}")
        sys.exit(1)

    # Extracting parameters safely using .get() with defaults where necessary
    VOLUME_PATH = params.get("IMAGES_PATH", "")
    MASK_PATH = params.get("MASK_PATH", "")
    OUTPUT_DIR = params.get("OUTPUT_PATH", "./output")
    OUTPUT_FORMAT = params.get("OUTPUT_FORMAT", "jp2")
    OUTPUT_TYPE = params.get("OUTPUT_TYPE", "8bit")
    WRITE_VECTORS = params.get("WRITE_VECTORS", False)
    WRITE_ANGLES = params.get("WRITE_ANGLES", True)
    SIGMA = params.get("SIGMA", 3.0)
    RHO = params.get("RHO", 1.0)
    AXIS_PTS = params.get("AXIS_POINTS", None)
    IS_TEST = params.get("TEST", False)
    N_SLICE_TEST = params.get("N_SLICE_TEST", None)

    print("\n---------------------------------")
    print("READING VOLUME INFORMATION\n")
    print(f"Volume path: {VOLUME_PATH}")

    data_reader = DataReader(VOLUME_PATH)

    if end_index is None:
        end_index = data_reader.shape[0]

    print(f"Number of slices: {data_reader.shape[0]}")

    # Check if already done the processing
    if not IS_TEST:
        is_already_done = True

        if end_index is None:
            is_already_done = False

        for idx in range(start_index, end_index):
            existing_files = []

            if WRITE_ANGLES:
                existing_files.extend(
                    [
                        f"{OUTPUT_DIR}/HA/HA_{idx:06d}.{OUTPUT_FORMAT}",
                        f"{OUTPUT_DIR}/IA/IA_{idx:06d}.{OUTPUT_FORMAT}",
                        f"{OUTPUT_DIR}/FA/FA_{idx:06d}.{OUTPUT_FORMAT}",
                    ]
                )

            if WRITE_VECTORS:
                existing_files.append(f"{OUTPUT_DIR}/eigen_vec/eigen_vec_{idx:06d}.npy")

            # Remove small files before checking existence
            remove_corrupted_files(existing_files)

            # If any required file is missing, mark processing as needed
            if not all(os.path.exists(file) for file in existing_files):
                is_already_done = False

        if is_already_done:
            print("\nAll images are already processed. Skipping computation.\n")
            return

    if MASK_PATH:
        is_mask = True
    else:
        print("Mask path not provided")
        is_mask = False

    print("\n---------------------------------")
    print("CALCULATE CENTER LINE\n")
    center_line = interpolate_points(AXIS_PTS, data_reader.shape[0])

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
        start_index,
        end_index,
        data_reader.shape[0],
        padding_start,
        padding_end,
        IS_TEST,
        N_SLICE_TEST,
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
        print("\n---------------------------------")
        print("LOAD MASK\n")
        mask_reader = DataReader(MASK_PATH)

        mask = mask_reader.load_volume(
            start_index_padded, end_index_padded, unbinned_shape=data_reader.shape
        ).astype("float32")

        assert mask.shape == volume.shape, (
            f"Mask shape {mask.shape} does not match volume shape {volume.shape}"
        )

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
        num_slices = vec.shape[1]
        print(f"Using {mp.cpu_count()} CPU cores")

        def update_bar(_):
            """Callback function to update progress bar."""
            bar()

        with mp.Pool(processes=mp.cpu_count()) as pool:
            with alive_bar(
                num_slices, title="Processing slices (Multiprocess)", bar="smooth"
            ) as bar:
                results = []
                for z in range(num_slices):
                    result = pool.apply_async(
                        compute_slice_angles_and_anisotropy,
                        (
                            z,
                            vec[:, z, :, :],
                            volume[z, :, :],
                            np.around(center_line[z]),
                            val[:, z, :, :],
                            center_line,
                            OUTPUT_DIR,
                            OUTPUT_FORMAT,
                            OUTPUT_TYPE,
                            start_index,
                            WRITE_VECTORS,
                            WRITE_ANGLES,
                            IS_TEST,
                        ),
                        callback=update_bar,  # ✅ Update progress bar after each task
                    )
                    results.append(result)

                for result in results:
                    result.wait()  # Ensure all tasks are completed before exiting
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
                    center_line,
                    OUTPUT_DIR,
                    OUTPUT_FORMAT,
                    OUTPUT_TYPE,
                    start_index,
                    WRITE_VECTORS,
                    WRITE_ANGLES,
                    IS_TEST,
                )
                bar()  # Update the progress bar for each slice

    print(f"\n🤖 - Finished processing slices {start_index} - {end_index}")
    print("---------------------------------\n\n")

    return


def compute_slice_angles_and_anisotropy(
    z: int,
    vector_field_slice: np.ndarray,
    img_slice: np.ndarray,
    center_point: np.ndarray,
    eigen_val_slice: np.ndarray,
    center_line: np.ndarray,
    OUTPUT_DIR: str,
    OUTPUT_FORMAT: str,
    OUTPUT_TYPE: str,
    start_index: int,
    WRITE_VECTORS: bool,
    WRITE_ANGLES: bool,
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
        center_line (np.ndarray): Center line for alignment.
        OUTPUT_DIR (str): Directory to save the output.
        OUTPUT_FORMAT (str): Format for the output files (e.g., "tif").
        OUTPUT_TYPE (str): Type of output (e.g., "8bits", "rgb").
        start_index (int): Start index of the slice.
        WRITE_VECTORS (bool): Whether to output vector fields.
        WRITE_ANGLES (bool): Whether to output angles and fractional anisotropy.
        IS_TEST (bool): Whether in test mode.

    Returns:
        None
    """
    # print(f"Processing image: {start_index + z}")

    paths = []
    if WRITE_ANGLES:
        paths = [
            f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.{OUTPUT_FORMAT}",
            f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.{OUTPUT_FORMAT}",
            f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.{OUTPUT_FORMAT}",
        ]
    if WRITE_VECTORS:
        paths.append(f"{OUTPUT_DIR}/eigen_vec/eigen_vec_{(start_index + z):06d}.npy")

    # Skip processing if files already exist (unless in test mode)
    if not IS_TEST and all(os.path.exists(path) for path in paths):
        # print(f"File {(start_index + z):06d} already exists")
        return

    buffer = 5
    if z < buffer:
        VEC_PTS = center_line[: min(z + buffer, len(center_line))]
    elif z >= len(center_line) - buffer:
        VEC_PTS = center_line[max(z - buffer, 0) :]
    else:
        VEC_PTS = center_line[z - buffer : z + buffer]

    center_vec = calculate_center_vector(VEC_PTS)
    # print(f"(Center vector: {center_vec})")

    # Compute angles and FA if needed
    if WRITE_ANGLES or IS_TEST:
        img_FA = compute_fraction_anisotropy(eigen_val_slice)
        vector_field_slice_rotated = rotate_vectors_to_new_axis(
            vector_field_slice, center_vec
        )
        img_helix, img_intrusion = compute_helix_and_transverse_angles(
            vector_field_slice_rotated,
            center_point,
        )

    # Visualization in test mode
    if IS_TEST:
        plot_images(img_slice, img_helix, img_intrusion, img_FA, center_point)
        write_images(
            img_helix,
            img_intrusion,
            img_FA,
            start_index,
            OUTPUT_DIR + "/test_slice",
            OUTPUT_FORMAT,
            OUTPUT_TYPE,
            z,
        )
        return

    # Save results
    if WRITE_ANGLES:
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
    if WRITE_VECTORS:
        write_vector_field(vector_field_slice, start_index, OUTPUT_DIR, z)
