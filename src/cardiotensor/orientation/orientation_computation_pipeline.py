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
from cardiotensor.utils.utils import remove_corrupted_files



# --- small helpers ---
def check_already_processed(
    output_dir: str,
    start_index: int,
    end_index: int,
    write_vectors: bool,
    write_angles: bool,
    output_format: str
) -> bool:
    """
    Check whether all required output files already exist.

    Args:
        output_dir: Path to output directory.
        start_index: Start slice index.
        end_index: End slice index (exclusive).
        write_vectors: If True, eigenvector fields are expected.
        write_angles: If True, HA/IA/FA maps are expected.
        output_format: Output image format (e.g., "jp2", "tif").

    Returns:
        bool: True if all expected files exist and are valid, else False.
    """
    for idx in range(start_index, end_index):
        expected_files = []

        if write_angles:
            expected_files.extend([
                f"{output_dir}/HA/HA_{idx:06d}.{output_format}",
                f"{output_dir}/IA/IA_{idx:06d}.{output_format}",
                f"{output_dir}/FA/FA_{idx:06d}.{output_format}",
            ])

        if write_vectors:
            expected_files.append(f"{output_dir}/eigen_vec/eigen_vec_{idx:06d}.npy")

        # Remove small/corrupted files before checking
        remove_corrupted_files(expected_files)

        # If any file is missing, processing is required
        if not all(os.path.exists(file) for file in expected_files):
            return False

    print(f"Checking already processed files: All expected files exist in {output_dir}")
    return True



# --- main API ---
def compute_orientation(
    volume_path: str,
    mask_path: str | None = None,
    output_dir: str = "./output",
    output_format: str = "jp2",
    output_type: str = "8bit",
    sigma: float = 1.0,
    rho: float = 3.0,
    truncate: float = 4.0,
    axis_points: np.ndarray | None = None,
    vertical_padding: float | None = None,
    write_vectors: bool = False,
    write_angles: bool = True,
    use_gpu: bool = True,
    is_test: bool = False,
    n_slice_test: int | None = None,
    start_index: int = 0,
    end_index: int | None = None,
) -> None:
    """
    Compute the orientation for a volume dataset.

    Args:
        volume_path: Path to the 3D volume.
        mask_path: Optional binary mask path.
        output_dir: Output directory for results.
        output_format: Image format for results.
        output_type: Image type ("8bit" or "rgb").
        sigma: Noise scale for structure tensor.
        rho: Integration scale for structure tensor.
        truncate: Gaussian kernel truncation.
        axis_points: 3D points defining LV axis for cylindrical coordinates.
        vertical_padding: Padding slices for tensor computation.
        write_vectors: Whether to save eigenvectors.
        write_angles: Whether to save HA/IA/FA maps.
        use_gpu: Use GPU acceleration for tensor computation.
        is_test: If True, runs in test mode and outputs plots.
        n_slice_test: Number of slices to process in test mode.
        start_index: Start slice index.
        end_index: End slice index (None = last slice).
    """

    # --- Sanity checks ---
    if sigma > rho:
        raise ValueError("sigma must be <= rho")

    print(f"""
Parameters:
    - Volume path:    {volume_path}
    - Mask path:      {mask_path or "[None]"}
    - Output dir:     {output_dir}
    - Output format:  {output_format}
    - Output type:    {output_type}
    - sigma / rho:    {sigma} / {rho}
    - truncate:       {truncate}
    - Write angles:   {write_angles}
    - Write vectors:  {write_vectors}
    - Use GPU:        {use_gpu}
    - Test mode:      {is_test}
    """)


    print("\n" + "-" * 40)
    print("READING VOLUME INFORMATION")
    print("-" * 40 + "\n")

    print(f"Volume path: {volume_path}")

    data_reader = DataReader(volume_path)

    if end_index is None:
        end_index = data_reader.shape[0]

    print(f"Number of slices: {data_reader.shape[0]}")
    
        
    # --- Check if already processed ---
    print("Check if file is already processed...")
    if check_already_processed(
        output_dir,
        start_index,
        end_index,
        write_vectors,
        write_angles,
        output_format,
    ) and not is_test:
        print("\nAll images are already processed. Skipping computation.\n")
        return
    
    
    print("\n---------------------------------")
    print("CALCULATE CENTER LINE\n")
    center_line = interpolate_points(axis_points, data_reader.shape[0])

    print("\n---------------------------------")
    print("CALCULATE PADDING START AND ENDING INDEXES\n")

    if vertical_padding is None:
        vertical_padding = truncate * rho + 0.5

    padding_start = padding_end = math.ceil(vertical_padding)
    if not is_test:
        if padding_start > start_index:
            padding_start = start_index
        if padding_end > (data_reader.shape[0] - end_index):
            padding_end = data_reader.shape[0] - end_index
    if is_test:
        if n_slice_test > data_reader.shape[0]:
            sys.exit("Error: n_slice_test > number of images")

    print(f"Padding start, Padding end : {padding_start}, {padding_end}")
    start_index_padded, end_index_padded = adjust_start_end_index(
        start_index,
        end_index,
        data_reader.shape[0],
        padding_start,
        padding_end,
        is_test,
        n_slice_test,
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

    if mask_path is not None:
        print("\n---------------------------------")
        print("LOAD MASK\n")
        mask_reader = DataReader(mask_path)

        mask = mask_reader.load_volume(
            start_index_padded, end_index_padded, unbinned_shape=data_reader.shape
        )

        assert mask.shape == volume.shape, (
            f"Mask shape {mask.shape} does not match volume shape {volume.shape}"
        )

        volume[mask == 0] = 0

    print("\n" + "-" * 40)
    print("CALCULATING STRUCTURE TENSOR")
    print("-" * 40 + "\n")
    t1 = time.perf_counter()  # start time
    val, vec = calculate_structure_tensor(
        volume, sigma, rho, truncate=truncate, use_gpu=use_gpu
    )

    print(f"Vector field shape: {vec.shape}")

    if mask_path is not None:
        print("Applying mask to tensors and vectors...")

        volume[mask == 0] = np.nan
        val[0, :, :, :][mask == 0] = np.nan
        val[1, :, :, :][mask == 0] = np.nan
        val[2, :, :, :][mask == 0] = np.nan
        vec[0, :, :, :][mask == 0] = np.nan
        vec[1, :, :, :][mask == 0] = np.nan
        vec[2, :, :, :][mask == 0] = np.nan

        print("Masking complete")

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

    print("\n" + "-" * 40)
    print("ANGLE & ANISOTROPY CALCULATION")
    print("-" * 40 + "\n")

    if not is_test:
        num_slices = vec.shape[1]
        print(f"Using {mp.cpu_count()} CPU cores")

        def update_bar(_):
            """Callback function to update progress bar."""
            bar()

        # Limit the number of processors used to avoid exceeing the max number of handlers in Windows
        if sys.platform.startswith("win"):
            num_procs = min(mp.cpu_count(), 59)
        else:
            num_procs = mp.cpu_count()

        with mp.Pool(processes=num_procs) as pool:
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
                            output_dir,
                            output_format,
                            output_type,
                            start_index,
                            write_vectors,
                            write_angles,
                            is_test,
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
                    output_dir,
                    output_format,
                    output_type,
                    start_index,
                    write_vectors,
                    write_angles,
                    is_test,
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
    output_dir: str,
    output_format: str,
    output_type: str,
    start_index: int,
    write_vectors: bool,
    write_angles: bool,
    is_test: bool,
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
        output_dir (str): Directory to save the output.
        output_format (str): Format for the output files (e.g., "tif").
        output_type (str): Type of output (e.g., "8bits", "rgb").
        start_index (int): Start index of the slice.
        write_vectors (bool): Whether to output vector fields.
        write_angles (bool): Whether to output angles and fractional anisotropy.
        is_test (bool): Whether in test mode.

    Returns:
        None
    """
    # print(f"Processing image: {start_index + z}")

    paths = []
    if write_angles:
        paths = [
            f"{output_dir}/HA/HA_{(start_index + z):06d}.{output_format}",
            f"{output_dir}/IA/IA_{(start_index + z):06d}.{output_format}",
            f"{output_dir}/FA/FA_{(start_index + z):06d}.{output_format}",
        ]
    if write_vectors:
        paths.append(f"{output_dir}/eigen_vec/eigen_vec_{(start_index + z):06d}.npy")

    # Skip processing if files already exist (unless in test mode)
    if not is_test and all(os.path.exists(path) for path in paths):
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
    if write_angles or is_test:
        img_FA = compute_fraction_anisotropy(eigen_val_slice)
        vector_field_slice_rotated = rotate_vectors_to_new_axis(
            vector_field_slice, center_vec
        )
        img_helix, img_intrusion = compute_helix_and_transverse_angles(
            vector_field_slice_rotated,
            center_point,
        )

    # Visualization in test mode
    if is_test:
        plot_images(img_slice, img_helix, img_intrusion, img_FA, center_point)
        write_images(
            img_helix,
            img_intrusion,
            img_FA,
            start_index,
            output_dir + "/test_slice",
            output_format,
            output_type,
            z,
        )
        return

    # Save results
    if write_angles:
        write_images(
            img_helix,
            img_intrusion,
            img_FA,
            start_index,
            output_dir,
            output_format,
            output_type,
            z,
        )
    if write_vectors:
        write_vector_field(vector_field_slice, start_index, output_dir, z)
