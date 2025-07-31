import math
import multiprocessing as mp
import os
from pathlib import Path

import cv2
import numpy as np
from alive_progress import alive_bar
from skimage.measure import block_reduce

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import convert_to_8bit


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
    Downsamples a vector volume using multiprocessing.

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
        return

    # Load dimensions from the first npy in each block
    sample = np.load(npy_list[0])
    _, h, w = sample.shape

    tasks = [
        (block, bin_factor, h, w, bin_dir, idx) for idx, block in enumerate(blocks)
    ]

    with mp.Pool(processes=min(mp.cpu_count(), 16)) as pool:
        with alive_bar(len(tasks), title="Downsampling vector volumes") as bar:
            results = [
                pool.apply_async(
                    process_vector_block, args=task, callback=lambda _: bar()
                )
                for task in tasks
            ]
            for result in results:
                result.wait()


from pathlib import Path


def _process_chunk(
    slice_paths: list[Path],
    bin_factor: int,
    H: int,
    W: int,
    eig_out_dir: Path,
    block_idx: int,
) -> None:
    """
    Worker function for downsampling one block of fine‐scale slices.
    Loads each slice in `slice_paths` (each shape = (3, H, W)), averages them, downsamples in‐plane,
    renormalizes, and writes to eig_out_dir/eigen_vec_{block_idx:06d}.npy.
    """
    out_file = eig_out_dir / f"eigen_vec_{block_idx:06d}.npy"
    if out_file.exists():
        return

    # 1) Load exactly those `n_slices_in_block` into a small array of shape (3, n_slices, H, W)
    n_slices = len(slice_paths)
    arr_block = np.empty((3, n_slices, H, W), dtype=np.float32)
    for i, p in enumerate(slice_paths):
        arr_block[:, i, :, :] = np.load(p)

    # 2) Average across Z within this block: shape → (3, H, W)
    avg_block = arr_block.mean(axis=1)

    # 3) First renormalize each 3‐vector to unit length
    norms = np.linalg.norm(avg_block, axis=0, keepdims=True)  # (1, H, W)
    avg_unit = avg_block / np.maximum(norms, 1e-12)

    # 4) Downsample in‐plane (H, W) by block_reduce averaging each bin_factor×bin_factor patch
    Hc = math.ceil(H / bin_factor)
    Wc = math.ceil(W / bin_factor)
    downsampled = np.empty((3, Hc, Wc), dtype=np.float32)
    for comp in range(3):
        downsampled_comp = block_reduce(
            avg_unit[comp], block_size=(bin_factor, bin_factor), func=np.mean
        )
        downsampled[comp] = downsampled_comp

    # 5) Second renormalization
    norms2 = np.linalg.norm(downsampled, axis=0, keepdims=True)
    downsampled_unit = downsampled / np.maximum(norms2, 1e-12)

    # 6) Save coarse slice
    np.save(out_file, downsampled_unit)


def chunked_downsample_vector_volume_mp(
    input_npy_dir: Path,
    bin_factor: int,
    output_dir: Path,
) -> None:
    """
    Multi‐process + progress‐bar version of chunked_downsample_vector_volume.
    Reads a directory of per‐slice NumPy files (shape = (3, H, W) each),
    groups every `bin_factor` consecutive slices into blocks, averages, downsamples, renormalizes,
    and writes each block as one coarse slice in output_dir/bin{bin_factor}/eigen_vec/.

    Parameters
    ----------
    input_npy_dir : Path
        Directory containing fine‐scale “eigen_vec_*.npy” files, each shape (3, H, W).
    bin_factor : int
        Number of fine Z‐slices per block.
    output_dir : Path
        Base output directory. Will create output_dir/bin{bin_factor}/eigen_vec/.
    """
    all_files = sorted(input_npy_dir.glob("*.npy"))
    if not all_files:
        raise RuntimeError(f"No .npy files found in {input_npy_dir}")

    # Determine H, W from the first slice
    sample = np.load(all_files[0])
    if sample.ndim != 3 or sample.shape[0] != 3:
        raise RuntimeError(
            f"Expected each .npy to have shape (3, H, W), but got {sample.shape}"
        )
    _, H, W = sample.shape

    Z_full = len(all_files)
    num_blocks = math.ceil(Z_full / bin_factor)

    # Prepare output directory
    bin_dir = output_dir / f"bin{bin_factor}"
    eig_out_dir = bin_dir / "eigen_vec"
    eig_out_dir.mkdir(parents=True, exist_ok=True)

    # Quick skip: if all expected files exist, return immediately
    all_exist = True
    for block_idx in range(num_blocks):
        if not (eig_out_dir / f"eigen_vec_{block_idx:06d}.npy").exists():
            all_exist = False
            break
    if all_exist:
        return

    # Build task list: (slice_paths, bin_factor, H, W, eig_out_dir, block_idx)
    tasks = []
    for block_idx in range(num_blocks):
        z_start = block_idx * bin_factor
        z_end = min(z_start + bin_factor, Z_full)
        slice_paths = all_files[z_start:z_end]
        tasks.append((slice_paths, bin_factor, H, W, eig_out_dir, block_idx))

    # Launch multiprocessing pool with a progress bar
    cpu_count = min(mp.cpu_count(), len(tasks))
    with mp.Pool(processes=cpu_count) as pool:
        with alive_bar(len(tasks), title="Downsampling vector volumes") as bar:
            results = [
                pool.apply_async(_process_chunk, args=task, callback=lambda _: bar())
                for task in tasks
            ]
            for r in results:
                r.wait()


def process_image_block(
    file_list, block_idx, bin_factor, out_file, min_value, max_value
):
    """
    Process a Z-block of images by averaging along the Z axis,
    downsampling in XY, converting to 8-bit, and writing to disk.

    Args:
        file_list (list): List of file paths (entire volume stack).
        bin_factor (int): Binning factor for XY downsampling.
        out_file (Path): Output file path for the downsampled image.
        min_value (float): Minimum intensity for 8-bit scaling.
        max_value (float): Maximum intensity for 8-bit scaling.
    """
    h, w = cv2.imread(str(file_list[0]), cv2.IMREAD_UNCHANGED).shape
    block = file_list[block_idx * bin_factor : (block_idx + 1) * bin_factor]
    array = np.full((len(block), h, w), np.nan, dtype=np.float32)

    for i, p in enumerate(block):
        img = cv2.imread(str(p), cv2.IMREAD_UNCHANGED)
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
    cv2.imwrite(str(out_file), downsampled_8bit)
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

    tasks = []
    for block_idx in range(num_blocks):
        out_file = out_dir / f"{subfolder}_{block_idx:06d}.{out_ext}"
        if not out_file.exists():
            tasks.append(
                (file_list, block_idx, bin_factor, out_file, min_value, max_value)
            )

    if not tasks:
        print(f"✔️ All downsampled blocks already exist for '{subfolder}'. Skipping.")
        return

    with mp.Pool(processes=mp.cpu_count()) as pool:
        with alive_bar(len(tasks), title=f"Downsampling '{subfolder}' volume") as bar:
            results = [
                pool.apply_async(
                    process_image_block, args=task, callback=lambda _: bar()
                )
                for task in tasks
            ]
            for r in results:
                r.wait()
