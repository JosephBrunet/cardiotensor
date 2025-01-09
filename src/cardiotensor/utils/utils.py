import configparser
import os
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

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The configuration file {file_path} does not exist.")
    if not file_path[-5:] == ".conf":
        raise FileNotFoundError("The file is not a .conf file")

    config = configparser.ConfigParser()
    config.read(file_path)

    def parse_coordinates(section: str, option: str, fallback: str = "") -> np.ndarray:
        """Parse a comma-separated coordinate string into a NumPy array."""
        value = config.get(section, option, fallback=fallback).strip()
        if value:
            try:
                return np.array([float(coord) for coord in value.split(",")])
            except ValueError:
                raise ValueError(
                    f"Invalid coordinate format for {option} in [{section}]: {value}"
                )
        return np.array([])

    return {
        # DATASET
        "IMAGES_PATH": config.get("DATASET", "IMAGES_PATH", fallback="").strip(),
        "MASK_PATH": config.get("DATASET", "MASK_PATH", fallback="").strip(),
        "FLIP": config.getboolean("DATASET", "FLIP", fallback=True),
        "VOXEL_SIZE": config.getfloat("DATASET", "VOXEL_SIZE", fallback=1.0),
        # OUTPUT
        "OUTPUT_PATH": config.get("OUTPUT", "OUTPUT_PATH", fallback="").strip(),
        "OUTPUT_FORMAT": config.get("OUTPUT", "OUTPUT_FORMAT", fallback="jp2").strip(),
        "OUTPUT_TYPE": config.get("OUTPUT", "OUTPUT_TYPE", fallback="").strip(),
        "VECTORS": config.getboolean("OUTPUT", "VECTORS", fallback=False),
        # STRUCTURE TENSOR CALCULATION
        "SIGMA": config.getfloat("STRUCTURE TENSOR CALCULATION", "SIGMA", fallback=1.0),
        "RHO": config.getfloat("STRUCTURE TENSOR CALCULATION", "RHO", fallback=1.0),
        "N_CHUNK": config.getint(
            "STRUCTURE TENSOR CALCULATION", "N_CHUNK", fallback=100
        ),
        # LV AXIS COORDINATES
        "POINT_MITRAL_VALVE": parse_coordinates(
            "LV AXIS COORDINATES", "POINT_MITRAL_VALVE"
        ),
        "POINT_APEX": parse_coordinates("LV AXIS COORDINATES", "POINT_APEX"),
        # RUN
        "REVERSE": config.getboolean("RUN", "REVERSE", fallback=False),
        "MASK_REMOVAL": config.get("RUN", "MASK_REMOVAL", fallback="after").strip(),
        # TEST
        "TEST": config.getboolean("TEST", "TEST", fallback=False),
        "N_SLICE_TEST": config.getint("TEST", "N_SLICE_TEST", fallback=None),
    }


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

    # Clip values to ensure they are in 8-bit range
    img_clipped = np.clip(img_normalized, 0, 255)

    return img_clipped.astype(np.uint8)
