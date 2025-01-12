import math
import multiprocessing as mp
import os
from pathlib import Path

import cv2
import numpy as np
from alive_progress import alive_bar
from skimage.measure import block_reduce

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
    # print(f"Processing block: {idx}")
    output_file = output_dir / f"eigen_vec/eigen_vec_{idx:06d}.npy"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if output_file.exists():
        # print(f"Skipping block {idx}, file already exists.")
        return

    array = np.empty((3, len(block), h, w))
    bin_array = np.empty((3, math.ceil(h / bin_factor), math.ceil(w / bin_factor)))

    for i, p in enumerate(block):
        print(f"Reading file: {p.name}")
        array[:, i, :, :] = np.load(p)

    array = array.mean(axis=1)

    block_size = (bin_factor, bin_factor)
    for i in range(3):
        bin_array[i, :, :] = block_reduce(
            array[i, :, :], block_size=block_size, func=np.mean
        )

    np.save(output_file, bin_array.astype(np.float32))
    print(f"Saved block {idx} to {output_file}")


def downsample_vector_volume(
    input_npy: Path, bin_factor: int, output_dir: Path
) -> None:
    """
    Downsamples a vector volume using multiprocessing.

    Args:
        input_npy (Path): Path to the directory containing numpy files.
        bin_factor (int): Binning factor for downsampling.
        output_dir (Path): Path to the output directory.
    """
    output_dir = output_dir / f"bin{bin_factor}"
    os.makedirs(output_dir / "eigen_vec", exist_ok=True)

    npy_list = sorted(input_npy.glob("*.npy"))
    _, h, w = np.load(npy_list[0]).shape

    blocks = [npy_list[i : i + bin_factor] for i in range(0, len(npy_list), bin_factor)]

    tasks = [
        (block, bin_factor, h, w, output_dir, idx) for idx, block in enumerate(blocks)
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
                result.wait()  # Wait for the task to complete


def process_image_block(
    block: list[Path],
    bin_factor: int,
    h: int,
    w: int,
    output_dir: Path,
    idx: int,
) -> None:
    """
    Processes a single block of image files and saves the downsampled output.

    Args:
        block (List[Path]): List of file paths to the image files in the block.
        bin_factor (int): Binning factor for downsampling.
        h (int): Height of the data block.
        w (int): Width of the data block.
        output_dir (Path): Path to the output directory.
        idx (int): Index of the current block.
    """
    # print(f"Processing block: {idx}")
    output_file = output_dir / f"HA/HA_{idx:06d}.tif"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if output_file.exists():
        # print(f"Skipping block {idx}, file already exists.")
        return

    array = np.empty((len(block), h, w))
    bin_array = np.empty((math.ceil(h / bin_factor), math.ceil(w / bin_factor)))

    for i, p in enumerate(block):
        print(f"Reading file: {p.name}")
        array[i, :, :] = cv2.imread(str(p), cv2.IMREAD_UNCHANGED)

    array = array.mean(axis=0)

    block_size = (bin_factor, bin_factor)
    bin_array[:, :] = block_reduce(array[:, :], block_size=block_size, func=np.mean)
    bin_array = convert_to_8bit(bin_array, min_value=-90, max_value=90)

    cv2.imwrite(str(output_file), bin_array)
    print(f"Saved block {idx} to {output_file}")


def downsample_volume(
    input_path: Path, bin_factor: int, output_dir: Path, file_format: str
) -> None:
    """
    Downsamples a volume of images using multiprocessing.

    Args:
        input_path (Path): Path to the directory containing image files.
        bin_factor (int): Binning factor for downsampling.
        output_dir (Path): Path to the output directory.
    """
    output_dir = output_dir / f"bin{bin_factor}"
    os.makedirs(output_dir / "HA", exist_ok=True)

    HA_list = sorted(input_path.glob(f"*.{file_format}"))
    h, w = cv2.imread(str(HA_list[0]), cv2.IMREAD_UNCHANGED).shape

    blocks = [HA_list[i : i + bin_factor] for i in range(0, len(HA_list), bin_factor)]
    tasks = [
        (block, bin_factor, h, w, output_dir, idx) for idx, block in enumerate(blocks)
    ]

    # print(f"Total blocks to process: {len(blocks)}")

    with mp.Pool(processes=min(mp.cpu_count(), 16)) as pool:
        with alive_bar(len(tasks), title="Downsampling image volumes") as bar:
            results = [
                pool.apply_async(
                    process_image_block, args=task, callback=lambda _: bar()
                )
                for task in tasks
            ]
            for result in results:
                result.wait()  # Wait for the task to complete
