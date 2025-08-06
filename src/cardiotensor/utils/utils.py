import ast
import configparser
import os
from pathlib import Path
from typing import Any

import numpy as np


def read_conf_file(file_path: str) -> dict[str, Any]:
    """
    Reads and parses a configuration file into a dictionary.

    Args:
        file_path (str): Path to the configuration file.

    Returns:
        Dict[str, Any]: Parsed configuration parameters.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
        ValueError: If expected numerical or array values are incorrectly formatted.
    """

    file_path = str(Path(file_path).resolve())  # Ensure the file path is absolute
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The configuration file {file_path} does not exist.")

    if not file_path.endswith(".conf"):
        raise ValueError("The file is not a .conf file")

    config = configparser.ConfigParser()
    config.read(file_path)

    def parse_coordinates(section: str, option: str, fallback: str = ""):
        """
        Parses a coordinate string from a configuration file into a list of tuples.

        Parameters:
            config (ConfigParser): The configuration parser object.
            section (str): The section in the config.
            option (str): The option name.
            fallback (str): A fallback value if the option is missing.

        Returns:
            list[tuple[int, int, int]]: List of 3D coordinate tuples or an empty list.
        """
        value = config.get(section, option, fallback=fallback).strip().replace(" ", "")

        if not value:  # Return an empty list if the value is empty
            return []

        try:
            value = f"[{value}]"  # Ensure it's formatted as a list
            points_list = ast.literal_eval(value)  # Safely evaluate the string
            return [tuple(point) for point in points_list]  # Convert to list of tuples
        except (SyntaxError, ValueError) as e:
            raise ValueError(
                f"Invalid coordinate format for {option} in [{section}]: {value}"
            ) from e

    # Read the two paths
    images_path = config.get("DATASET", "IMAGES_PATH").strip()
    mask_path = config.get("DATASET", "MASK_PATH", fallback="").strip()

    # Existence check (file or directory)
    if not os.path.exists(images_path):
        raise FileNotFoundError(f"The IMAGES_PATH '{images_path}' does not exist.")
    if mask_path and not os.path.exists(mask_path):
        raise FileNotFoundError(f"The MASK_PATH '{mask_path}' does not exist.")

    # OUTPUT
    output_dir = config.get("OUTPUT", "OUTPUT_PATH", fallback="").strip()
    if not output_dir:
        output_dir = "./output"  # Default output directory if not specified


    return {
        # DATASET
        "IMAGES_PATH": images_path,
        "MASK_PATH": mask_path,
        "VOXEL_SIZE": config.getfloat("DATASET", "VOXEL_SIZE", fallback=1.0),
        # STRUCTURE TENSOR CALCULATION
        "SIGMA": config.getfloat("STRUCTURE TENSOR CALCULATION", "SIGMA", fallback=3.0),
        "RHO": config.getfloat("STRUCTURE TENSOR CALCULATION", "RHO", fallback=1.0),
        "TRUNCATE": config.getfloat(
            "STRUCTURE TENSOR CALCULATION", "TRUNCATE", fallback=4.0
        ),
        "VERTICAL_PADDING": config.getfloat(
            "STRUCTURE TENSOR CALCULATION", "VERTICAL_PADDING", fallback=None
        ),
        "N_CHUNK": config.getint(
            "STRUCTURE TENSOR CALCULATION", "N_CHUNK", fallback=100
        ),
        "USE_GPU": config.getboolean(
            "STRUCTURE TENSOR CALCULATION", "USE_GPU", fallback=False
        ),
        "WRITE_VECTORS": config.getboolean(
            "STRUCTURE TENSOR CALCULATION", "WRITE_VECTORS", fallback=False
        ),
        "REVERSE": config.getboolean(
            "STRUCTURE TENSOR CALCULATION", "REVERSE", fallback=False
        ),
        # ANGLE CALCULATION
        "WRITE_ANGLES": config.getboolean(
            "ANGLE CALCULATION", "WRITE_ANGLES", fallback=True
        ),
        "AXIS_POINTS": parse_coordinates("ANGLE CALCULATION", "AXIS_POINTS"),
        # TEST
        "TEST": config.getboolean("TEST", "TEST", fallback=False),
        "N_SLICE_TEST": config.getint("TEST", "N_SLICE_TEST", fallback=None),
        # OUTPUT
        "OUTPUT_PATH": output_dir,
        "OUTPUT_FORMAT": config.get("OUTPUT", "OUTPUT_FORMAT", fallback="jp2").strip(),
        "OUTPUT_TYPE": config.get("OUTPUT", "OUTPUT_TYPE", fallback="8bit").strip(),
    }


# Function to remove files smaller than 1KB
def remove_corrupted_files(file_paths, size_threshold=200):
    for file_path in file_paths:
        if os.path.exists(file_path) and os.path.getsize(file_path) < size_threshold:
            print("Corrupted file removed:", file_path)
            os.remove(file_path)


def convert_to_8bit(
    img: np.ndarray,
    perc_min: int = 0,
    perc_max: int = 100,
    min_value: float | None = None,
    max_value: float | None = None,
) -> np.ndarray:
    """
    Converts a NumPy array to an 8-bit image.

    Args:
        img (np.ndarray): Input image array.
        perc_min (int): Minimum percentile for normalization. Default is 0.
        perc_max (int): Maximum percentile for normalization. Default is 100.
        min_value (Optional[float]): Optional explicit minimum value.
        max_value (Optional[float]): Optional explicit maximum value.

    Returns:
        np.ndarray: 8-bit converted image.
    """
    # Compute percentiles
    minimum, maximum = np.nanpercentile(img, (perc_min, perc_max))

    # Override percentiles with explicit min/max if provided
    if min_value is not None and max_value is not None:
        minimum, maximum = min_value, max_value

    # Avoid division by zero
    if maximum == minimum:
        return np.zeros_like(img, dtype=np.uint8)

    # Normalize and scale to 8-bit range
    img_normalized = (img - minimum) / (maximum - minimum) * 255

    # Convert NaN and Inf values to numbers
    img_cleaned = np.nan_to_num(img_normalized, nan=0.0, posinf=255.0, neginf=0.0)

    # Clip values to ensure they are in 8-bit range
    img_clipped = np.clip(img_cleaned, 0, 255)

    return img_clipped.astype(np.uint8)
