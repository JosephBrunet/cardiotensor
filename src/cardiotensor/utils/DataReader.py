import sys
from os import PathLike
from pathlib import Path
from typing import Any

import cv2
import dask
import numpy as np
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
        self.shape = self._get_volume_shape()

    def _get_volume_info(self) -> dict:
        """
        Detects the type of volume based on the path and retrieves relevant information.

        Returns:
            dict: Volume information containing type, stack status, and file list (if applicable).
        """
        volume_info: dict[str, str | bool | list[Path]] = {
            "type": "",
            "stack": False,
            "file_list": [],
        }

        if not self.path.exists():
            raise ValueError(f"The path does not exist: {self.path}")

        if self.path.is_dir():
            volume_info["stack"] = True
            image_files = {
                ext: sorted(self.path.glob(f"*.{ext}"))
                for ext in self.supported_extensions
            }
            volume_info["type"], volume_info["file_list"] = max(
                image_files.items(), key=lambda item: len(item[1])
            )

            assert isinstance(
                volume_info["file_list"], list
            )  # Runtime check for type safety

            if not volume_info["file_list"]:
                raise ValueError(
                    "No supported image files found in the specified directory."
                )
            volume_info["file_list"] = sorted(volume_info["file_list"])

        elif self.path.is_file() and self.path.suffix == ".mhd":
            volume_info["type"] = "mhd"

        if volume_info["type"] == "":
            raise ValueError(f"Unsupported volume type for path: {self.path}")

        return volume_info

    def _get_volume_shape(self) -> tuple[int, int, int]:
        """
        Retrieves the size (dimensions) of the volume.

        Returns:
            Tuple[int, int, int]: Dimensions of the volume (z, y, x).
        """
        if not self.volume_info["stack"]:  # Single file (e.g., .mhd)
            if self.volume_info["type"] == "mhd":
                image = sitk.ReadImage(str(self.path))
                return tuple(image.GetSize())  # Return (z, y, x)

        elif self.volume_info["stack"]:  # Stack of images
            first_image = self._custom_image_reader(self.volume_info["file_list"][0])
            return (
                len(self.volume_info["file_list"]),
                first_image.shape[0],
                first_image.shape[1],
            )

        raise ValueError("Unable to determine volume dimensions.")

    def load_volume(
        self,
        start_index: int = 0,
        end_index: int | None = None,
        unbinned_shape: tuple[int, int, int] | None = None,
    ) -> np.ndarray:
        """
        Loads the volume data based on the detected volume type.

        Args:
            start_index (int): Start index for slicing (for stacks).
            end_index (int): End index for slicing (for stacks). If 0, loads the entire stack.
            unbinned_shape (tuple): Shape of the volume without downsampling. Default is None (no binning).

        Returns:
            np.ndarray: Loaded volume data.
        """

        if end_index is None:
            end_index = self.shape[0]

        binning_factor = 1.0
        if unbinned_shape is not None:
            binning_factor = unbinned_shape[0] / self.shape[0]
            print(f"Mask bining factor: {binning_factor}")

        if binning_factor != 1.0:
            start_index_ini = start_index
            end_index_ini = end_index
            start_index = int((start_index_ini / binning_factor) - 1)
            if start_index < 0:
                start_index = 0
            end_index = int((end_index_ini / binning_factor) + 1)
            if end_index > self.shape[0]:
                end_index = self.shape[0]

            print(
                f"Mask start index padded: {start_index} - Mask end index padded : {end_index}"
            )

        if self.volume_info["stack"] == False:
            if self.volume_info["type"] == "mhd":
                volume, _ = _load_raw_data_with_mhd(self.path)
                volume = volume[start_index:end_index, :, :]
        elif self.volume_info["stack"] == True:
            # if end_index is None:
            #     end_index = len(self.volume_info["file_list"])
            volume = self._load_image_stack(
                self.volume_info["file_list"], start_index, end_index
            )
        else:
            raise ValueError("Unsupported volume type.")

        if binning_factor != 1.0 and unbinned_shape is not None:
            print("Resizing mask")
            volume = zoom(
                volume,
                zoom=binning_factor,
                order=0,
            )

            start_index_bin_upscaled = int(
                np.abs(start_index * binning_factor - start_index_ini)
            )
            end_index_bin_upscaled = start_index_bin_upscaled + (
                end_index_ini - start_index_ini
            )

            if start_index_bin_upscaled < 0:
                start_index_bin_upscaled = 0
            if end_index_bin_upscaled > volume.shape[0]:
                end_index_bin_upscaled = volume.shape[0]

            volume = volume[start_index_bin_upscaled:end_index_bin_upscaled, :]

            volume_resized = np.empty_like(volume)
            for i in range(volume.shape[0]):
                # Resize the slice to match the corresponding slice of the volume
                volume_resized[i] = cv2.resize(
                    volume[i],
                    (unbinned_shape[2], unbinned_shape[1]),
                    interpolation=cv2.INTER_LINEAR,
                )

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
        print(f"Loading {total_files} files...")

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
            volume = np.stack(computed_data, axis=0)

        return volume


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
