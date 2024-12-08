import configparser
import os
import sys
from os import PathLike
from pathlib import Path
from typing import Any

import cv2
import dask
import numpy as np
import SimpleITK as sitk
from alive_progress import alive_bar


def read_conf_file(file_path: str) -> dict[str, Any]:
    """
    Reads and parses a configuration file into a dictionary.

    Args:
        file_path (str): Path to the configuration file.

    Returns:
        Dict[str, Any]: Parsed configuration parameters.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The configuration file {file_path} does not exist.")

    config = configparser.ConfigParser()
    config.read(file_path)

    return {
        "IMAGES_PATH": config.get("DATASET", "IMAGES_PATH", fallback="").strip(),
        "MASK_PATH": config.get("DATASET", "MASK_PATH", fallback="").strip(),
        "FLIP": config.getboolean("DATASET", "FLIP", fallback=True),
        "OUTPUT_PATH": config.get("OUTPUT", "OUTPUT_PATH", fallback="").strip(),
        "OUTPUT_FORMAT": config.get("OUTPUT", "OUTPUT_FORMAT", fallback="jp2").strip(),
        "OUTPUT_TYPE": config.get("OUTPUT", "OUTPUT_TYPE", fallback="").strip(),
        "VECTORS": config.getboolean("OUTPUT", "VECTORS", fallback=False),
        "SIGMA": config.getfloat("STRUCTURE TENSOR CALCULATION", "SIGMA", fallback=1.0),
        "RHO": config.getfloat("STRUCTURE TENSOR CALCULATION", "RHO", fallback=1.0),
        "N_CHUNK": config.getint(
            "STRUCTURE TENSOR CALCULATION", "N_CHUNK", fallback=100
        ),
        "POINT_MITRAL_VALVE": np.array(
            [
                float(value)
                for value in config.get(
                    "LV AXIS COORDINATES", "POINT_MITRAL_VALVE", fallback=""
                ).split(",")
            ]
        ),
        "POINT_APEX": np.array(
            [
                float(value)
                for value in config.get(
                    "LV AXIS COORDINATES", "POINT_APEX", fallback=""
                ).split(",")
            ]
        ),
        "REVERSE": config.getboolean("RUN", "REVERSE", fallback=False),
        "TEST": config.getboolean("TEST", "TEST", fallback=True),
        "N_SLICE_TEST": config.getint("TEST", "N_SLICE_TEST", fallback=None),
    }


def convert_to_8bit(
    img: np.ndarray,
    perc_min: int = 0,
    perc_max: int = 100,
    output_min: float | None = None,
    output_max: float | None = None,
) -> np.ndarray:
    """
    Converts a NumPy array to an 8-bit image.

    Args:
        img (np.ndarray): Input image array.
        perc_min (int): Minimum percentile for normalization. Default is 0.
        perc_max (int): Maximum percentile for normalization. Default is 100.
        output_min (Optional[float]): Optional explicit minimum value.
        output_max (Optional[float]): Optional explicit maximum value.

    Returns:
        np.ndarray: 8-bit converted image.
    """
    minimum, maximum = np.nanpercentile(img, (perc_min, perc_max))

    if output_min is not None and output_max is not None:
        minimum, maximum = output_min, output_max

    img_normalized = (img + abs(minimum)) * (255 / (maximum - minimum))
    return img_normalized.astype(np.uint8)


def get_image_list(directory: str | Path) -> tuple[list[Path], str]:
    """
    Identifies the predominant image type in a directory and returns file paths.

    Args:
        directory (Union[str, Path]): Path to the directory.

    Returns:
        Tuple[List[Path], str]: List of file paths and the determined image type.

    Raises:
        ValueError: If no supported image files are found.
    """
    directory = Path(directory)
    image_extensions = ["tif", "jp2", "png", "edf"]

    image_files = {ext: sorted(directory.glob(f"*.{ext}")) for ext in image_extensions}
    predominant_type, img_list = max(image_files.items(), key=lambda item: len(item[1]))

    if not img_list:
        raise ValueError("No supported image files found in the specified directory.")

    return img_list, predominant_type


def get_volume_shape(volume_path: str) -> tuple[int, int, int]:
    """
    Determines the shape of a 3D volume.

    Args:
        volume_path (str): Path to the volume (either a .mhd file or folder).

    Returns:
        Tuple[int, int, int]: Dimensions of the volume (x, y, z).

    Raises:
        ValueError: If no valid files are found in the folder or the path is invalid.
    """
    if os.path.isfile(volume_path) and volume_path.endswith(".mhd"):
        image = sitk.ReadImage(volume_path)
        return image.GetSize()

    elif os.path.isdir(volume_path):
        supported_extensions = (".tif", ".tiff", ".jp2")
        slice_files = [
            f
            for f in os.listdir(volume_path)
            if f.lower().endswith(supported_extensions)
        ]
        z_dim = len(slice_files)

        if z_dim == 0:
            raise ValueError(
                "No supported image slices found in the specified directory."
            )

        first_image = cv2.imread(
            os.path.join(volume_path, slice_files[0]), cv2.IMREAD_UNCHANGED
        )
        y_dim, x_dim = first_image.shape[:2]
        return x_dim, y_dim, z_dim

    else:
        raise ValueError(
            "Invalid volume path. Must be a .mhd file or folder containing slices."
        )


def load_volume(
    file_list: list[Path], start_index: int = 0, end_index: int | None = None
) -> np.ndarray:
    """
    Loads a 3D volume from a list of files.

    Args:
        file_list (List[Path]): List of file paths.
        start_index (int): Starting index for loading. Default is 0.
        end_index (Optional[int]): Ending index for loading. Default is None (loads all files).

    Returns:
        np.ndarray: The loaded 3D volume.
    """
    end_index = end_index or len(file_list)

    def read_image(file_path: Path) -> np.ndarray:
        if file_path.suffix == ".npy":
            return np.load(file_path)
        return cv2.imread(str(file_path), cv2.IMREAD_UNCHANGED)

    with alive_bar(end_index - start_index, title="Loading Volume") as bar:
        delayed_tasks = [
            dask.delayed(lambda x: np.array(read_image(x)))(file)
            for file in sorted(file_list[start_index:end_index])
        ]
        volume = np.array(dask.compute(*delayed_tasks))
        bar()

    return volume


def read_mhd(filename: PathLike) -> dict[str, Any]:
    """
    Return a dictionary of meta data from MHD meta header file

    :param filename: file of type .mhd that should be loaded
    :returns: dictionary of meta data
    """

    # Define tags
    meta_dict = {}
    tag_set = []
    tag_set.extend(
        [
            "ObjectType",
            "NDims",
            "DimSize",
            "ElementType",
            "ElementDataFile",
            "ElementNumberOfChannels",
        ]
    )
    tag_set.extend(
        ["BinaryData", "BinaryDataByteOrderMSB", "CompressedData", "CompressedDataSize"]
    )
    tag_set.extend(
        [
            "Offset",
            "CenterOfRotation",
            "AnatomicalOrientation",
            "ElementSpacing",
            "TransformMatrix",
        ]
    )
    tag_set.extend(
        [
            "Comment",
            "SeriesDescription",
            "AcquisitionDate",
            "AcquisitionTime",
            "StudyDate",
            "StudyTime",
        ]
    )

    tag_flag = [False] * len(tag_set)

    with open(filename) as fn:
        line = fn.readline()
        while line:
            tags = str.split(line, "=")
            # print(tags[0])
            for i in range(len(tag_set)):
                tag = tag_set[i]
                if (str.strip(tags[0]) == tag) and (not tag_flag[i]):
                    # print(tags[1])
                    content = str.strip(tags[1])
                    if tag in [
                        "ElementSpacing",
                        "Offset",
                        "CenterOfRotation",
                        "TransformMatrix",
                    ]:
                        meta_dict[tag] = [float(s) for s in content.split()]
                    elif tag in ["NDims", "ElementNumberOfChannels"]:
                        meta_dict[tag] = int(content)
                    elif tag in ["DimSize"]:
                        meta_dict[tag] = [int(s) for s in content.split()]
                    elif tag in [
                        "BinaryData",
                        "BinaryDataByteOrderMSB",
                        "CompressedData",
                    ]:
                        if content == "True":
                            meta_dict[tag] = True
                        else:
                            meta_dict[tag] = False
                    else:
                        meta_dict[tag] = content
                    tag_flag[i] = True
            line = fn.readline()
    return meta_dict


def load_raw_data_with_mhd(filename: PathLike) -> tuple[np.ndarray, dict[str, Any]]:
    """
    Load a MHD file

    :param filename: file of type .mhd that should be loaded
    :returns: tuple with raw data and dictionary of meta data
    """
    meta_dict = read_mhd(filename)
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
    arr.reverse()
    if element_channels > 1:
        data = data.reshape(arr + [element_channels])
    else:
        data = data.reshape(arr)
    # End 3D fix

    return (data, meta_dict)
