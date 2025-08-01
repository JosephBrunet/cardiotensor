import sys
from os import PathLike
from pathlib import Path
from typing import Any

import cv2
import dask
import numpy as np
import psutil
import SimpleITK as sitk
from alive_progress import alive_bar
from scipy.ndimage import zoom


class DataReader:
    def __init__(self, path: str | Path):
        """
        Initializes the DataReader with a path to the volume.

        Args:
            path (str | Path): Path to the volume directory or file.
        """
        self.path = Path(path)
        self.supported_extensions = ["tif", "tiff", "jp2", "png", "npy"]
        self.volume_info = self._get_volume_info()

    # ---------------------------
    # Properties
    # ---------------------------
    @property
    def shape(self) -> tuple[int, ...]:
        """Returns the shape of the volume as (Z, Y, X) or (Z, Y, X, C)."""
        return self.volume_info["shape"]

    @property
    def dtype(self) -> np.dtype:
        """Returns the data type of the volume."""
        return self.volume_info["dtype"]

    @property
    def volume_size_gb(self) -> float:
        """Returns the total size of the volume in GB."""
        n_bytes = np.prod(self.shape) * np.dtype(self.dtype).itemsize
        return n_bytes / (1024**3)

    def _get_volume_info(self) -> dict:
        """
        Detects volume type, shape, and dtype.
        Returns a dict with keys: type, stack, file_list, shape, dtype
        """
        volume_info = {
            "type": "",
            "stack": False,
            "file_list": [],
            "shape": None,
            "dtype": None,
        }

        if not self.path.exists():
            raise ValueError(f"The path does not exist: {self.path}")

        # Case 1: Directory of images
        if self.path.is_dir():
            volume_info["stack"] = True
            image_files = {
                ext: sorted(self.path.glob(f"*.{ext}"))
                for ext in self.supported_extensions
            }
            volume_info["type"], volume_info["file_list"] = max(
                image_files.items(), key=lambda item: len(item[1])
            )
            if not volume_info["file_list"]:
                raise ValueError(
                    "No supported image files found in the specified directory."
                )

            # Inspect first file
            first_image = self._custom_image_reader(volume_info["file_list"][0])
            volume_info["dtype"] = first_image.dtype

            # Shape: handle scalar vs vector
            if (
                volume_info["type"] == "npy"
                and first_image.ndim == 3
                and first_image.shape[0] == 3
            ):
                # Vector field stored as (3, Y, X)
                volume_info["shape"] = (
                    3,
                    len(volume_info["file_list"]),
                    first_image.shape[1],
                    first_image.shape[2],
                )
            elif first_image.ndim == 3:
                # 4D scalar stack (Z from files)
                volume_info["shape"] = (
                    len(volume_info["file_list"]),
                    *first_image.shape,
                )
            else:
                # Standard 3D stack
                volume_info["shape"] = (
                    len(volume_info["file_list"]),
                    first_image.shape[0],
                    first_image.shape[1],
                )

        # Case 2: Single MHD file
        elif self.path.is_file() and self.path.suffix == ".mhd":
            volume_info["type"] = "mhd"
            img = sitk.ReadImage(str(self.path))
            arr = sitk.GetArrayFromImage(img)  # Z, Y, X
            volume_info["shape"] = arr.shape
            volume_info["dtype"] = arr.dtype

        else:
            raise ValueError(f"Unsupported volume type for path: {self.path}")

        return volume_info

    def load_volume(
        self,
        start_index: int = 0,
        end_index: int | None = None,
        unbinned_shape: tuple[int, int, int] | None = None,
    ) -> np.ndarray:
        """
        Loads the volume and resizes it to unbinned_shape if provided, using fast scipy.ndimage.zoom.

        Args:
            start_index (int): Start index for slicing (for stacks).
            end_index (int): End index for slicing (for stacks). If None, loads the entire stack.
            unbinned_shape (tuple): Desired shape (Z, Y, X). If None, no resizing is done.

        Returns:
            np.ndarray: Loaded volume.
        """
        if end_index is None:
            end_index = self.shape[0]

        self.check_memory_requirement(self.shape, self.dtype)

        # Decide if resize is needed
        need_resize = False
        if unbinned_shape is not None and self.shape != unbinned_shape:
            need_resize = True
            zoom_factors = tuple(u / s for u, s in zip(unbinned_shape, self.shape))
            print(f"Zoom factors: {zoom_factors}")
        else:
            zoom_factors = (1.0, 1.0, 1.0)

        # Optional padding if resizing
        if need_resize:
            start_index_ini, end_index_ini = start_index, end_index
            start_index = int(start_index_ini / zoom_factors[0]) - 1
            start_index = max(start_index, 0)
            end_index = int(end_index_ini / zoom_factors[0]) + 1
            end_index = min(end_index, self.shape[0])
            print(f"Volume start index padded: {start_index} - end: {end_index}")

        # Load volume from stack or mhd
        if not self.volume_info["stack"]:
            if self.volume_info["type"] == "mhd":
                volume, _ = _load_raw_data_with_mhd(self.path)
                volume = volume[start_index:end_index, :, :]
        else:
            volume = self._load_image_stack(
                self.volume_info["file_list"], start_index, end_index
            )

        if need_resize:
            print("Resizing with scipy.ndimage.zoom...")

            # Ensure float32 for better memory and speed
            volume = volume.astype(np.float32)
            volume = zoom(
                volume, zoom=zoom_factors, order=0
            )  # Nearest-neighbor for mask

            # Final crop to original range
            crop_start = int(abs(start_index * zoom_factors[0] - start_index_ini))
            crop_end = crop_start + (end_index_ini - start_index_ini)
            crop_start = max(crop_start, 0)
            crop_end = min(crop_end, volume.shape[0])

            volume = volume[crop_start:crop_end, :, :]

        return volume

    def _custom_image_reader(self, file_path: Path) -> np.ndarray:
        """
        Reads an image from the given file path.

        Args:
            file_path (Path): Path to the image file.

        Returns:
            np.ndarray: Image data as a NumPy array.
        """
        if file_path.suffix == ".npy":
            return np.load(file_path)
        return cv2.imread(str(file_path), cv2.IMREAD_UNCHANGED)

    def _load_image_stack(
        self, file_list: list[Path], start_index: int, end_index: int
    ) -> np.ndarray:
        """
        Loads a stack of images into a 3D NumPy array.

        Args:
            file_list (List[Path]): List of file paths to load.
            start_index (int): Start index for slicing.
            end_index (int): End index for slicing.

        Returns:
            np.ndarray: Loaded volume as a 3D array.
        """
        if end_index == 0:
            end_index = len(file_list)

        total_files = end_index - start_index
        # print(f"Loading {total_files} files...")

        if start_index < 0 or end_index > len(file_list):
            raise ValueError(
                f"Invalid indices: start_index={start_index}, end_index={end_index}, total_files={len(file_list)}"
            )

        with alive_bar(total_files, title="Loading Volume", length=40) as bar:

            def progress_bar_reader(file_path: Path) -> np.ndarray:
                bar()  # Update the progress bar
                return self._custom_image_reader(file_path)

            delayed_tasks = [
                dask.delayed(progress_bar_reader)(file_path)
                for file_path in sorted(file_list[start_index:end_index])
            ]

            # Compute the volume
            computed_data = dask.compute(*delayed_tasks)

            # Validate shape consistency
            first_shape = computed_data[0].shape
            for idx, data in enumerate(computed_data):
                if data.shape != first_shape:
                    raise ValueError(
                        f"Inconsistent file shape at index {idx}: Expected {first_shape}, got {data.shape}"
                    )

        # Combine into a NumPy array
        print("Stacking images into a 3D volume...")
        if self.volume_info["type"] == "npy":
            volume = np.stack(computed_data, axis=1)
        else:
            volume = np.stack(computed_data, axis=0)

        return volume

    def check_memory_requirement(self, shape, dtype, safety_factor=0.8):
        """
        Check if the dataset can fit in available memory.

        Args:
            shape (tuple[int]): Shape of the array.
            dtype (np.dtype): NumPy dtype of the array.
            safety_factor (float): Fraction of available memory allowed to be used.
        """
        # Compute dataset size in bytes
        n_bytes = np.prod(shape) * np.dtype(dtype).itemsize
        size_gb = n_bytes / (1024**3)

        # Check available memory
        available_gb = psutil.virtual_memory().available / (1024**3)

        print(
            f"Dataset size: {size_gb:.2f} GB | Available memory: {available_gb:.2f} GB"
        )

        if size_gb > available_gb * safety_factor:
            print("❌ Dataset is too large to safely load into memory.")
            sys.exit(1)


def _read_mhd(filename: PathLike[str]) -> dict[str, Any]:
    """
    Return a dictionary of meta data from an MHD meta header file.

    Args:
        filename (PathLike[str]): File path to the .mhd file.

    Returns:
        dict[str, Any]: A dictionary containing parsed metadata.
    """
    meta_dict: dict[str, Any] = {}
    tag_set = [
        "ObjectType",
        "NDims",
        "DimSize",
        "ElementType",
        "ElementDataFile",
        "ElementNumberOfChannels",
        "BinaryData",
        "BinaryDataByteOrderMSB",
        "CompressedData",
        "CompressedDataSize",
        "Offset",
        "CenterOfRotation",
        "AnatomicalOrientation",
        "ElementSpacing",
        "TransformMatrix",
        "Comment",
        "SeriesDescription",
        "AcquisitionDate",
        "AcquisitionTime",
        "StudyDate",
        "StudyTime",
    ]

    with open(filename) as fn:
        for line in fn:
            tags = line.split("=")
            if len(tags) < 2:
                continue
            key, content = tags[0].strip(), tags[1].strip()
            if key in tag_set:
                if key in [
                    "ElementSpacing",
                    "Offset",
                    "CenterOfRotation",
                    "TransformMatrix",
                ]:
                    # Parse as a list of floats
                    meta_dict[key] = [float(value) for value in content.split()]
                elif key in ["NDims", "ElementNumberOfChannels"]:
                    # Parse as an integer
                    meta_dict[key] = int(content)
                elif key == "DimSize":
                    # Parse as a list of integers
                    meta_dict[key] = [int(value) for value in content.split()]
                elif key in ["BinaryData", "BinaryDataByteOrderMSB", "CompressedData"]:
                    # Parse as a boolean
                    meta_dict[key] = content.lower() == "true"
                else:
                    # Parse as a string
                    meta_dict[key] = content
    return meta_dict


def _load_raw_data_with_mhd(
    filename: PathLike[str],
) -> tuple[np.ndarray, dict[str, Any]]:
    """
    Load a MHD file

    :param filename: file of type .mhd that should be loaded
    :returns: tuple with raw data and dictionary of meta data
    """
    meta_dict = _read_mhd(filename)
    dim = int(meta_dict["NDims"])
    if "ElementNumberOfChannels" in meta_dict:
        element_channels = int(meta_dict["ElementNumberOfChannels"])
    else:
        element_channels = 1

    if meta_dict["ElementType"] == "MET_FLOAT":
        np_type = np.float32
    elif meta_dict["ElementType"] == "MET_DOUBLE":
        np_type = np.float64
    elif meta_dict["ElementType"] == "MET_CHAR":
        np_type = np.byte
    elif meta_dict["ElementType"] == "MET_UCHAR":
        np_type = np.ubyte
    elif meta_dict["ElementType"] == "MET_SHORT":
        np_type = np.int16
    elif meta_dict["ElementType"] == "MET_USHORT":
        np_type = np.ushort
    elif meta_dict["ElementType"] == "MET_INT":
        np_type = np.int32
    elif meta_dict["ElementType"] == "MET_UINT":
        np_type = np.uint32
    else:
        raise NotImplementedError(
            "ElementType " + meta_dict["ElementType"] + " not understood."
        )
    arr = list(meta_dict["DimSize"])

    volume = np.prod(arr[0 : dim - 1])

    pwd = Path(filename).parents[0].resolve()
    data_file = Path(meta_dict["ElementDataFile"])
    if not data_file.is_absolute():
        data_file = pwd / data_file

    shape = (arr[dim - 1], volume, element_channels)
    with open(data_file, "rb") as f:
        data = np.fromfile(f, count=np.prod(shape), dtype=np_type)
    data.shape = shape

    # Adjust byte order in numpy array to match default system byte order
    if "BinaryDataByteOrderMSB" in meta_dict:
        sys_byteorder_msb = sys.byteorder == "big"
        file_byteorder_ms = meta_dict["BinaryDataByteOrderMSB"]
        if sys_byteorder_msb != file_byteorder_ms:
            data = data.byteswap()

    # Begin 3D fix
    # arr.reverse()
    if element_channels > 1:
        data = data.reshape(arr + [element_channels])
    else:
        data = data.reshape(arr)
    # End 3D fix

    return (data, meta_dict)
