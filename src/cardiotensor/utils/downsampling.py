import math
import os
from multiprocessing.pool import ThreadPool
from pathlib import Path

import numpy as np
from alive_progress import alive_bar
from skimage.measure import block_reduce

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.image_io import read_image_file, write_image
from cardiotensor.utils.utils import convert_to_8bit, get_available_cpu_count


def _downsample_vector_axes(vector_block: np.ndarray, bin_factor: int) -> np.ndarray:
    """Downsample unoriented vectors without cancelling equivalent v/-v axes."""
    if vector_block.ndim != 4 or vector_block.shape[0] != 3:
        raise ValueError("vector_block must have shape (3, Z, Y, X)")
    if bin_factor <= 0:
        raise ValueError("bin_factor must be > 0")

    vectors = np.asarray(vector_block, dtype=np.float32)
    vectors = np.nan_to_num(vectors, copy=False)
    _, depth, height, width = vectors.shape
    out_h = math.ceil(height / bin_factor)
    out_w = math.ceil(width / bin_factor)
    pad_h = out_h * bin_factor - height
    pad_w = out_w * bin_factor - width
    if pad_h or pad_w:
        vectors = np.pad(vectors, ((0, 0), (0, 0), (0, pad_h), (0, pad_w)))

    # Each output pixel is represented by a 3x3 orientation tensor:
    # sum(v @ v.T). Its principal eigenvector is invariant to vector sign.
    samples = vectors.reshape(3, depth, out_h, bin_factor, out_w, bin_factor).transpose(
        2, 4, 1, 3, 5, 0
    )
    samples = samples.reshape(out_h, out_w, -1, 3)
    valid = np.any(samples, axis=(2, 3))
    axes = np.zeros((out_h, out_w, 3), dtype=np.float32)
    if not np.any(valid):
        return np.moveaxis(axes, -1, 0)

    all_valid = np.all(valid)
    active_samples = samples if all_valid else samples[valid]
    tensors = np.einsum(
        "...ni,...nj->...ij", active_samples, active_samples, optimize=True
    )
    eigenvalues, eigenvectors = np.linalg.eigh(tensors)
    active_axes = eigenvectors[..., -1].astype(np.float32, copy=False)

    # Choose a deterministic polarity for storage. The physical orientation
    # remains an axis, so this does not change the result.
    eps = np.finfo(np.float32).eps
    z_zero = np.abs(active_axes[..., 2]) <= eps
    y_zero = np.abs(active_axes[..., 1]) <= eps
    flip = active_axes[..., 2] < -eps
    flip |= z_zero & (active_axes[..., 1] < -eps)
    flip |= z_zero & y_zero & (active_axes[..., 0] < 0)
    active_axes[flip] *= -1
    active_axes[eigenvalues[..., -1] <= eps] = 0
    if all_valid:
        axes = active_axes
    else:
        axes[valid] = active_axes

    return np.moveaxis(axes, -1, 0)


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
        for i, p in enumerate(block):
            array[:, i, :, :] = np.load(p)

        bin_array = _downsample_vector_axes(array, bin_factor)
        np.save(output_file, bin_array)

    except Exception as e:
        print(f"Error processing block {idx}: {e}")
        # print(f"Failed to process files: {[str(p) for p in block]}")
        raise


def _process_vector_range(
    vector_reader: DataReader,
    start_index: int,
    end_index: int,
    bin_factor: int,
    output_dir: Path,
    output_index: int,
    mask_reader: DataReader | None = None,
    mask_y_indices: np.ndarray | None = None,
    mask_x_indices: np.ndarray | None = None,
    overwrite: bool = False,
) -> None:
    """Load and downsample one Z range from any DataReader vector backend."""
    output_file = output_dir / f"eigen_vec_{output_index:06d}.npy"
    if output_file.exists() and not overwrite:
        return

    mask = None
    if mask_reader is not None:
        vector_depth = vector_reader.shape[1]
        mask_depth = mask_reader.shape[0]
        mask_z_indices = (
            np.arange(start_index, end_index, dtype=np.int64)
            * mask_depth
            // vector_depth
        )
        mask_start = int(mask_z_indices[0])
        mask = mask_reader.load_volume(
            mask_start,
            int(mask_z_indices[-1]) + 1,
            show_progress=False,
        )
        mask = np.take(mask, mask_y_indices, axis=1)
        mask = np.take(mask, mask_x_indices, axis=2)
        mask = mask[mask_z_indices - mask_start]
        expected_shape = (
            end_index - start_index,
            vector_reader.shape[2],
            vector_reader.shape[3],
        )
        if mask.shape != expected_shape:
            raise ValueError(
                f"Mask block shape {mask.shape} does not match vector block "
                f"shape {expected_shape}"
            )
        if not np.any(mask):
            output_shape = (
                3,
                math.ceil(vector_reader.shape[2] / bin_factor),
                math.ceil(vector_reader.shape[3] / bin_factor),
            )
            np.save(output_file, np.zeros(output_shape, dtype=np.float32))
            return

    vector_block = vector_reader.load_volume(
        start_index, end_index, show_progress=False
    )
    if mask is not None:
        vector_block[:, mask <= 0] = 0
    np.save(output_file, _downsample_vector_axes(vector_block, bin_factor))


def downsample_vector_volume(
    input_npy: Path,
    bin_factor: int,
    output_dir: Path,
    mask_path: str | Path | None = None,
) -> Path:
    """
    Downsamples a vector volume using a thread pool.

    Args:
        input_npy (Path): Path to the directory containing numpy files.
        bin_factor (int): Binning factor for downsampling.
        output_dir (Path): Path to the output directory.
        mask_path: Optional mask applied before vector downsampling.
    """
    bin_dir = output_dir / f"bin{bin_factor}"
    eig_out_dir = bin_dir / "eigen_vec"
    os.makedirs(eig_out_dir, exist_ok=True)

    input_npy = Path(input_npy)
    mask_path = None if mask_path is None else Path(mask_path)
    reader = DataReader(input_npy)
    if len(reader.shape) != 4 or reader.shape[0] != 3:
        raise ValueError(
            f"Expected vector field shape (3, Z, Y, X), got {reader.shape}"
        )
    total_slices = reader.shape[1]
    mask_reader = None
    mask_y_indices = None
    mask_x_indices = None
    if mask_path is not None:
        mask_reader = DataReader(mask_path)
        if len(mask_reader.shape) != 3:
            raise ValueError(f"Expected a 3D mask, got shape {mask_reader.shape}")
        _, vector_height, vector_width = reader.shape[1:]
        _, mask_height, mask_width = mask_reader.shape
        mask_y_indices = (
            np.arange(vector_height, dtype=np.int64) * mask_height // vector_height
        )
        mask_x_indices = (
            np.arange(vector_width, dtype=np.int64) * mask_width // vector_width
        )
        print(
            f"Applying mask shape {mask_reader.shape} to vector field "
            f"shape {reader.shape[1:]} before binning"
        )
    blocks = [
        (start, min(start + bin_factor, total_slices))
        for start in range(0, total_slices, bin_factor)
    ]
    total_blocks = len(blocks)

    # Do not reuse vectors produced with a different mask. A missing marker also
    # invalidates caches created by older Cardiotensor versions.
    mask_source = (
        "mask-index-v3:none"
        if mask_path is None
        else f"mask-index-v3:{mask_path.resolve()}"
    )
    cache_marker = eig_out_dir / ".mask_source"
    cache_matches = (
        cache_marker.exists() and cache_marker.read_text().strip() == mask_source
    )
    all_exist = cache_matches
    for idx in range(total_blocks):
        expected_file = eig_out_dir / f"eigen_vec_{idx:06d}.npy"
        if not expected_file.exists():
            all_exist = False
            break
    if all_exist:
        print("✅ Downsampled images for eigen_vec already exist. Skipping.")
        return eig_out_dir

    tasks = [
        (
            reader,
            start,
            end,
            bin_factor,
            eig_out_dir,
            idx,
            mask_reader,
            mask_y_indices,
            mask_x_indices,
            not cache_matches,
        )
        for idx, (start, end) in enumerate(blocks)
    ]

    with ThreadPool(processes=min(get_available_cpu_count(), 32)) as pool:
        with alive_bar(len(tasks), title="Downsampling vector volumes") as bar:
            results = [
                pool.apply_async(
                    _process_vector_range, args=task, callback=lambda _: bar()
                )
                for task in tasks
            ]
            for result in results:
                result.get()
    cache_marker.write_text(mask_source + "\n")
    return eig_out_dir


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

    cache_key = f"scalar-v2:{Path(input_path).resolve()}:{reader.shape}:bin{bin_factor}"
    cache_marker = out_dir / ".source_complete"
    cache_matches = (
        cache_marker.exists() and cache_marker.read_text().strip() == cache_key
    )
    expected_files = [
        out_dir / f"{subfolder}_{i:06d}.{out_ext}" for i in range(num_blocks)
    ]
    if cache_matches and all(f.exists() for f in expected_files):
        print(f"✅ Downsampled images for '{subfolder}' already exist. Skipping.")
        return

    sample = reader._custom_image_reader(file_list[0])
    if sample.ndim == 3 and sample.shape[2] == 1:
        sample = sample[:, :, 0]
    if sample.ndim != 2:
        raise ValueError(
            f"Expected 2D image slices for downsampling, got {sample.shape}"
        )
    H, W = sample.shape

    tasks = []
    for block_idx in range(num_blocks):
        out_file = out_dir / f"{subfolder}_{block_idx:06d}.{out_ext}"
        if not cache_matches or not out_file.exists():
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
                r.get()
    cache_marker.write_text(cache_key + "\n")
