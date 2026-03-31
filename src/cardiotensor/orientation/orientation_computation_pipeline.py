import math
import multiprocessing as mp
import os
import sys
import time
from collections.abc import Sequence

import numpy as np
from alive_progress import alive_bar
from matplotlib import pyplot as plt

# from memory_profiler import profile
from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.orientation.orientation_computation_functions import (
    adjust_start_end_index,
    calculate_center_vector,
    calculate_structure_tensor,
    compute_azimuth_and_elevation,
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
def _resolve_colormap(colormap: str | None):
    """Resolve colormap names, including the project-specific helix angle map."""
    if colormap is None:
        return None
    cmap_name = colormap.strip()
    if not cmap_name:
        return None
    if cmap_name.lower() in {"helix_angle", "helix_angle_cmap"}:
        return helix_angle_cmap
    try:
        return plt.get_cmap(cmap_name)
    except ValueError as err:
        raise ValueError(
            f"Unknown colormap '{colormap}'. Use a valid matplotlib colormap name "
            "or 'helix_angle'."
        ) from err


def check_already_processed(
    output_dir: str,
    start_index: int,
    end_index: int,
    write_vectors: bool,
    write_angles: bool,
    output_format: str,
    angle_names: tuple[str, str] = ("HA", "IA"),
    fa_name: str = "FA",
    extra_expected: Sequence[str] | None = None,
) -> bool:
    """
    Check whether all required output files already exist for every slice index.

    Parameters
    ----------
    output_dir : str
        Base output directory.
    start_index : int
        First global slice index to check (inclusive).
    end_index : int
        Last global slice index to check (exclusive).
    write_vectors : bool
        If True, expect eigenvector .npy files (e.g., eigen_vec_{idx:06d}.npy).
    write_angles : bool
        If True, expect angle images for angle_names[0], angle_names[1], and FA.
    output_format : str
        Image format/extension for angles, for example "jp2" or "tif".
    angle_names : tuple[str, str], optional
        Names of the two angle outputs, e.g. ("HA", "IA") or ("AZ", "EL").
    fa_name : str, optional
        Name of the FA subfolder, default "FA".
    extra_expected : sequence of str, optional
        Additional per-slice path templates to check. Each template must contain
        "{idx}" which will be formatted as a zero-padded integer (06d), and may
        also contain "{ext}" for the image extension.

    Returns
    -------
    bool
        True if all expected files for all indices exist (and pass the quick
        corruption filter), False otherwise.
    """
    if start_index >= end_index:
        # Nothing to check
        return True

    # Normalize extension
    ext = output_format.lstrip(".")
    if not ext:
        raise ValueError(
            "output_format must be a non-empty extension like 'jp2' or 'tif'."
        )

    # Prepare optional extras
    extra_expected = tuple(extra_expected or ())

    for idx in range(start_index, end_index):
        expected_files = []

        if write_angles:
            a1, a2 = angle_names
            expected_files += [
                os.path.join(output_dir, a1, f"{a1}_{idx:06d}.{ext}"),
                os.path.join(output_dir, a2, f"{a2}_{idx:06d}.{ext}"),
                os.path.join(output_dir, fa_name, f"{fa_name}_{idx:06d}.{ext}"),
            ]

        if write_vectors:
            expected_files.append(
                os.path.join(output_dir, "eigen_vec", f"eigen_vec_{idx:06d}.npy")
            )

        # User-specified extras, if any
        for tmpl in extra_expected:
            expected_files.append(tmpl.format(idx=f"{idx:06d}", ext=ext))

        # Remove small/corrupted files before checking (function defined elsewhere)
        remove_corrupted_files(expected_files)

        # If any file is missing, we need to process
        if not all(os.path.exists(p) for p in expected_files):
            return False

    print(f"Checking already processed files: all expected files exist in {output_dir}")
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
    angle_mode: str = "ha_ia",
    write_angles: bool = True,
    use_gpu: bool = True,
    is_test: bool = False,
    n_slice_test: int | None = None,
    show_quiver: bool = True,
    start_index: int = 0,
    end_index: int | None = None,
    colormap: str | None = None,
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
        show_quiver: If True, overlay the vector field on the test-slice figure.
        start_index: Start slice index.
        end_index: End slice index (None = last slice).
        colormap: Colormap name for RGB angle outputs. Defaults to helix-angle map.
    """

    # --- Sanity checks ---
    if sigma > rho:
        raise ValueError("sigma must be <= rho")

    if angle_mode.lower() == "ha_ia":
        angle_names = ("HA", "IA")
    elif angle_mode.lower() == "az_el":
        angle_names = ("AZ", "EL")
    else:
        raise ValueError("ANGLE_MODE must be 'ha_ia' or 'az_el'")

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
    - Angle mode:     {angle_mode}  -> {angle_names[0]}, {angle_names[1]}
    - Write vectors:  {write_vectors}
    - Use GPU:        {use_gpu}
    - Test mode:      {is_test}
    - Show quiver:    {show_quiver}
    - Colormap:       {colormap or "[default]"}
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
    if (
        check_already_processed(
            output_dir,
            start_index,
            end_index,
            write_vectors,
            write_angles,
            output_format,
            angle_names=angle_names,
        )
        and not is_test
    ):
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
            """Callback to tick the progress bar after each finished task."""
            bar()

        # Limit number of processes to avoid exceeding handler limits on Windows
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
                            show_quiver,
                            angle_mode,
                            colormap,
                        ),
                        callback=update_bar,
                    )
                    results.append(result)

                # Ensure all tasks complete before leaving the pool context
                for r in results:
                    r.wait()
    else:
        # Single threaded path with a progress bar
        with alive_bar(
            vec.shape[1], title="Processing slices (Single-thread)", bar="smooth"
        ) as bar:
            for z in range(vec.shape[1]):
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
                    show_quiver,
                    angle_mode,
                    colormap,
                )
                bar()

    end_index_local = start_index + vec.shape[1]
    print(f"\nFinished processing slices {start_index} to {end_index_local}")
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
    output_format: str = "jp2",
    output_type: str = "8bit",
    start_index: int = 0,
    write_vectors: bool = False,
    write_angles: bool = True,
    is_test: bool = False,
    show_quiver: bool = True,
    angle_mode: str = "ha_ia",
    colormap: str | None = None,
) -> None:
    """
    Compute either HA/IA or Azimuth/Elevation plus FA for a single slice,
    then plot and/or write outputs depending on flags.
    """
    # Decide angle labels and ranges based on mode
    mode = angle_mode.lower().strip()
    if mode == "ha_ia":
        angle_names = ("HA", "IA")
        angle_ranges = ((-90.0, 90.0), (-90.0, 90.0))
    elif mode == "az_el":
        angle_names = ("AZ", "EL")
        angle_ranges = ((-180.0, 180.0), (-90.0, 90.0))
    else:
        raise ValueError("ANGLE_MODE must be 'ha_ia' or 'az_el'")

    ext = output_format.lstrip(".")
    idx = start_index + z
    colormap_angle = _resolve_colormap(colormap)

    # Expected outputs for skip logic
    expected_paths = []
    if write_angles:
        a1, a2 = angle_names
        expected_paths = [
            os.path.join(output_dir, a1, f"{a1}_{idx:06d}.{ext}"),
            os.path.join(output_dir, a2, f"{a2}_{idx:06d}.{ext}"),
            os.path.join(output_dir, "FA", f"FA_{idx:06d}.{ext}"),
        ]
    if write_vectors:
        expected_paths.append(
            os.path.join(output_dir, "eigen_vec", f"eigen_vec_{idx:06d}.npy")
        )

    # Skip if all outputs are already present and we are not in test mode
    if (
        not is_test
        and expected_paths
        and all(os.path.exists(p) for p in expected_paths)
    ):
        return

    # Build a small window around the slice index to estimate the local axis direction
    buffer = 5
    if z < buffer:
        VEC_PTS = center_line[: min(z + buffer, len(center_line))]
    elif z >= len(center_line) - buffer:
        VEC_PTS = center_line[max(z - buffer, 0) :]
    else:
        VEC_PTS = center_line[z - buffer : z + buffer]

    center_vec = calculate_center_vector(VEC_PTS)

    # Compute FA and the chosen angle pair
    if write_angles or is_test:
        img_FA = compute_fraction_anisotropy(eigen_val_slice)
        vector_field_slice_rotated = rotate_vectors_to_new_axis(
            vector_field_slice, center_vec
        )

        if mode == "ha_ia":
            img_angle1, img_angle2 = compute_helix_and_transverse_angles(
                vector_field_slice_rotated, center_point
            )
        else:  # "az_el"
            img_angle1, img_angle2 = compute_azimuth_and_elevation(
                vector_field_slice_rotated
            )

    # Test mode: visualize a 2x2 figure and write to test subfolder
    if is_test:
        titles = (
            ("Helix Angle", "Intrusion Angle")
            if mode == "ha_ia"
            else ("Azimuth", "Elevation")
        )
        test_output_dir = os.path.join(output_dir, "test_slice")
        plot_images(
            img_slice,
            img_angle1,
            img_angle2,
            img_FA,
            center_point,
            vector_field_slice=vector_field_slice if show_quiver else None,
            colormap_angle=colormap_angle,
            angle1_title=titles[0],
            angle2_title=titles[1],
            save_path=os.path.join(test_output_dir, "result_test_slice.png"),
        )
        write_images(
            img_angle1,
            img_angle2,
            img_FA,
            start_index,
            test_output_dir,
            ext,
            output_type,
            z,
            colormap_angle=colormap_angle,
            angle_names=angle_names,
            angle_ranges=angle_ranges,
        )
        return

    # Persist outputs
    if write_angles:
        write_images(
            img_angle1,
            img_angle2,
            img_FA,
            start_index,
            output_dir,
            ext,
            output_type,
            z,
            colormap_angle=colormap_angle,
            angle_names=angle_names,
            angle_ranges=angle_ranges,
        )
    if write_vectors:
        write_vector_field(vector_field_slice, start_index, output_dir, z)
