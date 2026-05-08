import math
import os
import sys
from multiprocessing.pool import ThreadPool
from pathlib import Path

import numpy as np
from alive_progress import alive_bar
from skimage.measure import block_reduce

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.image_io import read_image_file, write_image
from cardiotensor.utils.utils import convert_to_8bit, get_available_cpu_count


def process_vector_block(
    block: list[Path],
    bin_factor: int,
    h: int,
    w: int,
    output_dir: Path,
    idx: int,
) -> None:
    """
    Processes a single block of numpy files and saves the downsampled output.

    Args:
        block (List[Path]): List of file paths to the numpy files in the block.
        bin_factor (int): Binning factor for downsampling.
        h (int): Height of the data block.
        w (int): Width of the data block.
        output_dir (Path): Path to the output directory.
        idx (int): Index of the current block.
    """

    try:
        output_file = output_dir / "eigen_vec" / f"eigen_vec_{idx:06d}.npy"
        output_file.parent.mkdir(parents=True, exist_ok=True)

        if output_file.exists():
            return

        array = np.empty((3, len(block), h, w), dtype=np.float32)
        bin_array = np.empty(
            (3, math.ceil(h / bin_factor), math.ceil(w / bin_factor)), dtype=np.float32
        )

        for i, p in enumerate(block):
            array[:, i, :, :] = np.load(p)

        array = array.mean(axis=1)

        block_size = (bin_factor, bin_factor)
        for comp in range(3):
            bin_array[comp] = block_reduce(
                array[comp], block_size=block_size, func=np.mean
            )

        np.save(output_file, bin_array.astype(np.float32))

    except Exception as e:
        print(f"Error processing block {idx}: {e}")
        # print(f"Failed to process files: {[str(p) for p in block]}")
        raise e


def downsample_vector_volume(
    input_npy: Path,
    bin_factor: int,
    output_dir: Path,
) -> None:
    """
    Downsamples a vector volume using a thread pool.

    Args:
        input_npy (Path): Path to the directory containing numpy files.
        bin_factor (int): Binning factor for downsampling.
        output_dir (Path): Path to the output directory.
    """
    bin_dir = output_dir / f"bin{bin_factor}"
    eig_out_dir = bin_dir / "eigen_vec"
    os.makedirs(eig_out_dir, exist_ok=True)

    npy_list = sorted(input_npy.glob("*.npy"))
    if len(npy_list) == 0:
        return

    # Determine block count
    blocks = [npy_list[i : i + bin_factor] for i in range(0, len(npy_list), bin_factor)]
    total_blocks = len(blocks)

    # Quick check: if all expected output files already exist, skip processing
    all_exist = True
    for idx in range(total_blocks):
        expected_file = eig_out_dir / f"eigen_vec_{idx:06d}.npy"
        if not expected_file.exists():
            all_exist = False
            break
    if all_exist:
        print("✅ Downsampled images for eigen_vec already exist. Skipping.")
        return

    # Load dimensions from the first npy in each block
    sample = np.load(npy_list[0])
    _, h, w = sample.shape

    tasks = [
        (block, bin_factor, h, w, bin_dir, idx) for idx, block in enumerate(blocks)
    ]

    with ThreadPool(processes=min(get_available_cpu_count(), 32)) as pool:
        with alive_bar(len(tasks), title="Downsampling vector volumes") as bar:
            results = [
                pool.apply_async(
                    process_vector_block, args=task, callback=lambda _: bar()
                )
                for task in tasks
            ]
            for result in results:
                result.wait()


def process_image_block(
    file_list, block_idx, bin_factor, h, w, out_file, min_value, max_value
):
    """
    Process a Z-block of images by averaging along the Z axis,
    downsampling in XY, converting to 8-bit, and writing to disk.

    Args:
        file_list (list): List of file paths (entire volume stack).
        bin_factor (int): Binning factor for XY downsampling.
        h (int): Input image height.
        w (int): Input image width.
        out_file (Path): Output file path for the downsampled image.
        min_value (float): Minimum intensity for 8-bit scaling.
        max_value (float): Maximum intensity for 8-bit scaling.
    """
    block = file_list[block_idx * bin_factor : (block_idx + 1) * bin_factor]
    array = np.full((len(block), h, w), np.nan, dtype=np.float32)

    for i, p in enumerate(block):
        img = read_image_file(p)
        if img is None or img.size == 0:
            print(f"⚠️ Warning: Failed to read image {p}, filling with NaNs")
            continue
        array[i] = img

    # Compute mean ignoring NaNs
    mean_z = np.nanmean(array, axis=0)

    # Optional: if all slices were NaN, handle gracefully
    if np.isnan(mean_z).all():
        raise ValueError("❌ All slices in this block are empty or unreadable")

    downsampled = block_reduce(
        mean_z, block_size=(bin_factor, bin_factor), func=np.mean
    )
    downsampled_8bit = convert_to_8bit(
        downsampled, min_value=min_value, max_value=max_value
    )
    write_image(out_file, downsampled_8bit)
    return True


def downsample_volume(
    input_path: Path,
    bin_factor: int,
    output_dir: Path,
    subfolder: str = "HA",
    out_ext: str = "tif",
    min_value: float = 0,
    max_value: float = 255,
) -> None:
    """
    Downsamples a 3D image volume along the Z and XY axes and saves as 8-bit images.

    This function reads a volumetric image dataset (e.g. TIFF stack) using DataReader,
    performs block averaging along the Z-axis and spatial downsampling in XY, then saves
    each resulting slice in a specified output directory as 8-bit images.

    Args:
        input_path (Path): Path to the directory containing the image stack.
        bin_factor (int): Factor to downsample in XY and the number of Z-slices to average per output slice.
        output_dir (Path): Path to the output root directory.
        subfolder (str): Subdirectory name under `binX/` to place results (default: "HA").
        out_ext (str): Output image format extension (e.g., 'tif', 'png').
        min_value (float): Minimum value for intensity normalization to 8-bit.
        max_value (float): Maximum value for intensity normalization to 8-bit.

    Returns:
        None
    """

    reader = DataReader(input_path)
    Z, H, W = reader.shape
    file_list = reader.volume_info["file_list"]

    num_blocks = math.ceil(Z / bin_factor)
    bin_dir = output_dir / f"bin{bin_factor}"
    out_dir = bin_dir / subfolder
    out_dir.mkdir(parents=True, exist_ok=True)

    # Early exit if all output files already exist
    expected_files = [
        out_dir / f"{subfolder}_{i:06d}.{out_ext}" for i in range(num_blocks)
    ]
    if all(f.exists() for f in expected_files):
        print(f"✅ Downsampled images for '{subfolder}' already exist. Skipping.")
        return

    sample = reader._custom_image_reader(file_list[0])
    if sample.ndim == 3 and sample.shape[2] == 1:
        sample = sample[:, :, 0]
    if sample.ndim != 2:
        raise ValueError(f"Expected 2D image slices for downsampling, got {sample.shape}")
    H, W = sample.shape

    tasks = []
    for block_idx in range(num_blocks):
        out_file = out_dir / f"{subfolder}_{block_idx:06d}.{out_ext}"
        if not out_file.exists():
            tasks.append(
                (file_list, block_idx, bin_factor, H, W, out_file, min_value, max_value)
            )

    if not tasks:
        print(f"✔️ All downsampled blocks already exist for '{subfolder}'. Skipping.")
        return

    cpu_count = min(get_available_cpu_count(), 32)

    with ThreadPool(processes=cpu_count) as pool:
        with alive_bar(len(tasks), title=f"Downsampling '{subfolder}' volume") as bar:
            results = [
                pool.apply_async(
                    process_image_block, args=task, callback=lambda _: bar()
                )
                for task in tasks
            ]
            for r in results:
                r.wait()
